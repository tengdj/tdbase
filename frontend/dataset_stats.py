"""Read full-resolution dataset counts without executing a spatial query."""
from collections import OrderedDict
from concurrent.futures import Future
import hashlib
import json
import math
import mmap
from pathlib import Path
import struct
import tempfile
import threading


META_2D = struct.Struct("<IIQ4d")
COUNT = struct.Struct("<Q")
POINT_2D = struct.Struct("<2f")
POINT_3D = struct.Struct("<3f")
SCHEMA_VERSION = 1


def _require(condition, message="数据文件结构无效或与所选维度不符"):
    if not condition:
        raise ValueError(message)


def polygon_stats(data):
    """Count every boundary ring, excluding a duplicated closing vertex."""
    size = len(data)
    _require(size >= 8)
    objects = COUNT.unpack_from(data, size - 8)[0]
    _require(objects <= (size - 8) // META_2D.size)
    payload_end = size - 8 - objects * META_2D.size
    vertices = 0
    for index in range(objects):
        record_size, outer_count, offset, *_ = META_2D.unpack_from(
            data, payload_end + index * META_2D.size)
        end = offset + record_size
        _require(record_size >= 16 and end <= payload_end)
        holes = COUNT.unpack_from(data, offset)[0]
        _require(holes <= (record_size - 16) // 8)
        cursor = offset + 8
        for ring_index in range(holes + 1):
            _require(cursor + 8 <= end)
            count = COUNT.unpack_from(data, cursor)[0]
            cursor += 8
            _require(count >= 3 and count <= (end - cursor) // 8)
            _require(ring_index != 0 or count == outer_count)
            first = POINT_2D.unpack_from(data, cursor)
            last = POINT_2D.unpack_from(data, cursor + (count - 1) * 8)
            vertices += count - int(first == last)
            cursor += count * 8
        _require(cursor == end)
    return {"objects": objects, "vertices": vertices, "edges": vertices,
            "format": "farm2d", "counting": "all_rings_without_closing_duplicate"}


def mesh_stats(data, compressed_reader, path):
    """Count LOD 100 only; triangle soup shares vertices within each object."""
    _require(len(data) >= 1 and data[0] in (0, 1))
    compressed = data[0] == 0
    cursor, objects, vertices, faces = 1, 0, 0, 0
    while cursor < len(data):
        _require(cursor + 8 <= len(data))
        data_size = COUNT.unpack_from(data, cursor)[0]
        start = cursor + 8
        metadata = start + data_size
        _require(metadata + 8 <= len(data))
        voxel_count = COUNT.unpack_from(data, metadata)[0]
        voxel_start = metadata + 8 + (0 if compressed else 40)
        stride = 36 if compressed else 116
        _require(voxel_start <= len(data) and voxel_count <= (len(data) - voxel_start) // stride)
        if compressed:
            # The header contains only the base mesh, so native progressive
            # decoding is required for the complete vertex/triangle counts.
            _require(data_size >= 34)
        else:
            points = set()
            for voxel in range(voxel_count):
                offset, triangles = struct.unpack_from("<QQ", data, voxel_start + voxel * stride + 100)
                _require(offset <= data_size and triangles <= (data_size - offset) // 44)
                faces += triangles
                begin = start + offset
                for point in POINT_3D.iter_unpack(data[begin:begin + triangles * 36]):
                    _require(all(math.isfinite(value) for value in point), "三维文件含无效顶点坐标")
                    points.add(point)
            vertices += len(points)
        objects += 1
        cursor = voxel_start + voxel_count * stride
    _require(cursor == len(data))
    if compressed:
        _require(compressed_reader is not None, "压缩三维数据的统计需要包含 mesh_stats 命令的 TDBase 版本")
        result = compressed_reader(path)
        _require(result.get("objects") == objects and
                 all(type(result.get(key)) is int and result[key] >= 0 for key in ("vertices", "faces")),
                 "TDBase 返回的网格统计无效")
        vertices, faces = result["vertices"], result["faces"]
    return {"objects": objects, "vertices": vertices, "faces": faces,
            "format": "compressed" if compressed else "decoded", "lod": 100,
            "counting": "distinct_coordinates_per_object"}


class DatasetStats:
    """File-signature cache; concurrent requests for one file share one scan."""

    def __init__(self, root, compressed_reader=None, cache_dir=None):
        self.root = Path(root)
        self.compressed_reader = compressed_reader
        self.cache_dir = Path(cache_dir) if cache_dir else None
        self.cache = OrderedDict()
        self.pending = {}
        self.lock = threading.Lock()

    def get(self, dimension, value):
        _require(dimension in ("2d", "3d"), "dimension 必须是 2d 或 3d")
        _require(isinstance(value, str) and value.strip() and "\x00" not in value, "请先选择数据文件")
        path = Path(value.strip()).expanduser()
        path = (self.root / path).resolve() if not path.is_absolute() else path.resolve()
        if not path.is_file():
            raise FileNotFoundError(f"找不到数据文件：{value}")
        stat = path.stat()
        key = (SCHEMA_VERSION, dimension, str(path), stat.st_size, stat.st_mtime_ns)
        with self.lock:
            if key in self.cache:
                self.cache.move_to_end(key)
                return dict(self.cache[key])
            waiting = self.pending.get(key)
            if waiting is None:
                waiting = Future()
                self.pending[key] = waiting
                owner = True
            else:
                owner = False
        if not owner:
            return dict(waiting.result())
        try:
            result = self._read(path, dimension, key)
            current = path.stat()
            _require((current.st_size, current.st_mtime_ns) == key[-2:], "文件正在更新，请稍后重新选择")
            with self.lock:
                self.cache[key] = result
                while len(self.cache) > 64:
                    self.cache.popitem(last=False)
            waiting.set_result(result)
            return dict(result)
        except Exception as error:
            waiting.set_exception(error)
            raise
        finally:
            with self.lock:
                self.pending.pop(key, None)

    def _read(self, path, dimension, key):
        cache_path = None
        if self.cache_dir:
            digest = hashlib.sha256(json.dumps(key).encode()).hexdigest()
            cache_path = self.cache_dir / f"{digest}.json"
            try:
                cached = json.loads(cache_path.read_text(encoding="utf-8"))
                required = ("objects", "vertices", "edges" if dimension == "2d" else "faces")
                if cached.get("signature") == list(key) and all(
                        type(cached.get("result", {}).get(name)) is int and cached["result"][name] >= 0
                        for name in required):
                    return cached["result"]
            except (OSError, ValueError, TypeError, AttributeError):
                pass
        with path.open("rb") as source:
            _require(path.stat().st_size > 0)
            with mmap.mmap(source.fileno(), 0, access=mmap.ACCESS_READ) as data:
                counts = polygon_stats(data) if dimension == "2d" else mesh_stats(data, self.compressed_reader, path)
        result = {"dimension": dimension, "path": str(path), **counts}
        if cache_path:
            temporary = None
            try:
                self.cache_dir.mkdir(parents=True, exist_ok=True)
                with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=self.cache_dir, delete=False) as output:
                    temporary = Path(output.name)
                    json.dump({"signature": key, "result": result}, output)
                temporary.replace(cache_path)
            except OSError:
                pass  # A read-only cache directory does not prevent statistics.
            finally:
                if temporary:
                    temporary.unlink(missing_ok=True)
        return result
