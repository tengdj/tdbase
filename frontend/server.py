#!/usr/bin/env python3
"""Local web UI/API for the unified TDBase + FARM executable.

Run from WSL:
    python3 frontend/server.py --build build/gpu --port 8765
"""
from __future__ import annotations

import argparse
import atexit
from collections import OrderedDict, deque
import json
import math
import os
from pathlib import Path
import re
import struct
import subprocess
import sys
import tempfile
import threading
import time
import uuid
from http import HTTPStatus
from html import escape
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import shlex
from urllib.parse import urlparse


ROOT = Path(__file__).resolve().parents[1]
FRONTEND = ROOT / "frontend"
# Support both direct startup and the existing file-based API test loaders.
if str(FRONTEND) not in sys.path:
    sys.path.insert(0, str(FRONTEND))
from dataset_stats import DatasetStats
DEFAULT_GPU = ROOT / "build" / "gpu" / "tdbase"
DEFAULT_CPU = ROOT / "build" / "cpu" / "tdbase"
MAX_LOG_OUTPUT = 2_000_000
PREVIEW_PAIRS = 100
MAX_VISUAL_NEIGHBORS = 150
MAX_VISUAL_VERTICES = 240
RESULT_CACHE_LIMIT = 3
RESULT_CACHE_TTL = 1800
RESULT_CACHE = OrderedDict()
RESULT_CACHE_LOCK = threading.Lock()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="TDBase/FARM local web UI")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--build", type=Path, help="directory containing tdbase")
    parser.add_argument("--gpu-build", type=Path, help="GPU build directory (overrides --build for GPU)")
    parser.add_argument("--cpu-build", type=Path, help="CPU build directory (overrides --build for CPU)")
    parser.add_argument("--timeout", type=int, default=600, help="query timeout in seconds")
    return parser.parse_args()


ARGS = parse_args()
BUILD_OVERRIDE = ARGS.build.resolve() if ARGS.build else None
GPU_BUILD_OVERRIDE = ARGS.gpu_build.resolve() if ARGS.gpu_build else None
CPU_BUILD_OVERRIDE = ARGS.cpu_build.resolve() if ARGS.cpu_build else None


def executable(backend: str) -> Path:
    override = GPU_BUILD_OVERRIDE if backend == "gpu" else CPU_BUILD_OVERRIDE
    if override:
        candidate = override / "tdbase"
    elif BUILD_OVERRIDE and backend == "gpu":
        # A directory explicitly named cpu is a CPU-only selection.
        candidate = (BUILD_OVERRIDE / "tdbase.gpu-missing"
                     if "cpu" in BUILD_OVERRIDE.name.lower()
                     else BUILD_OVERRIDE / "tdbase")
    elif BUILD_OVERRIDE and backend == "cpu":
        # A GPU build contains the same CPU-capable unified executable, but
        # prefer the separately built CPU artifact when it exists.
        candidate = (ROOT / "build" / "cpu" / "tdbase") if (ROOT / "build" / "cpu" / "tdbase").is_file() else BUILD_OVERRIDE / "tdbase"
    else:
        candidate = DEFAULT_GPU if backend == "gpu" else DEFAULT_CPU
    return candidate.resolve()


def compressed_dataset_stats(path: Path) -> dict:
    binary = next((executable(backend) for backend in ("cpu", "gpu")
                   if executable(backend).is_file()), None)
    if binary is None:
        raise FileNotFoundError("找不到用于统计压缩三维数据的 TDBase 程序")
    try:
        process = subprocess.run([str(binary), "mesh_stats", str(path)], cwd=ROOT,
                                 capture_output=True, text=True, timeout=180)
    except subprocess.TimeoutExpired as error:
        raise TimeoutError("压缩网格统计超时，请稍后重试") from error
    if process.returncode:
        raise ValueError("无法读取完整精度网格：" + process.stderr.strip()[-300:])
    try:
        result = json.loads(process.stdout)
        if not isinstance(result, dict):
            raise ValueError()
        return result
    except ValueError as error:
        raise ValueError("当前 TDBase 程序不支持网格统计，请重新编译") from error


DATASET_STATS = DatasetStats(ROOT, compressed_dataset_stats, ROOT / "build" / "dataset-stats")


def clean_path(value: object, name: str, required: bool = True) -> str:
    text = str(value or "").strip()
    if required and not text:
        raise ValueError(f"{name} 不能为空")
    if "\x00" in text:
        raise ValueError(f"{name} 包含非法字符")
    return text


def absolute_dataset_path(value: object, name: str, required: bool = True) -> str:
    text = clean_path(value, name, required)
    if not text:
        return ""
    path = Path(text).expanduser()
    return str((path if path.is_absolute() else ROOT / path).resolve())


def integer(value: object, name: str, minimum: int, maximum: int | None = None, default: int | None = None) -> int:
    if value in (None, ""):
        if default is not None:
            return default
        raise ValueError(f"{name} 不能为空")
    try:
        number = int(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} 必须是整数") from error
    if number < minimum or (maximum is not None and number > maximum):
        limit = f"{minimum}..{maximum}" if maximum is not None else f">={minimum}"
        raise ValueError(f"{name} 必须在 {limit} 范围内")
    return number


def build_command(payload: dict) -> tuple[list[str], Path, dict]:
    dimension = payload.get("dimension", "2d")
    backend = payload.get("backend", "cpu")
    query = payload.get("query", "intersect")
    if dimension not in {"2d", "3d"}:
        raise ValueError("dimension 必须是 2d 或 3d")
    if backend not in {"cpu", "gpu"}:
        raise ValueError("backend 必须是 cpu 或 gpu")
    allowed = {"intersect", "intersection", "within"} if dimension == "2d" else {"intersect", "within", "nn"}
    if query not in allowed:
        raise ValueError(f"{dimension} 不支持 {query} 查询")

    binary = executable(backend)
    if not binary.is_file():
        raise FileNotFoundError(f"找不到 {backend.upper()} 版本: {binary}")

    source = absolute_dataset_path(payload.get("source"), "源数据")
    # Both engines support self joins. FARM's within query and all TDBase
    # queries accept an omitted target and reuse the source dataset.
    target_required = dimension == "2d" and query != "within"
    target = absolute_dataset_path(payload.get("target"), "目标数据", target_required)
    command = [str(binary)]

    if dimension == "2d":
        command += ["farm", query, "-s", source]
        if target:
            command += ["-t", target]
        target_object = payload.get("target_object")
        if target_object not in (None, ""):
            command += ["--target-object", str(integer(str(target_object), "目标物体 ID", 0, 2_147_483_647))]
        if backend == "gpu":
            command += ["-g", "-b", "100000"]
        if payload.get("hierarchy"):
            command.append("-h")
        if payload.get("approximation"):
            command.append("-a")
        command += ["-n", str(integer(payload.get("threads"), "线程数", 1, 256, 8))]
        command += ["-v", str(integer(payload.get("vpr"), "二维栅格粒度", 1, 10000, 10))]
        command += ["-w", str(integer(payload.get("bitwidth"), "位宽", 2, 8, 3))]
        if query == "within":
            command += ["-d", str(integer(payload.get("distance"), "距离", 0, 10_000_000, 10))]
    else:
        command += ["join", "--dim", "3", "-q", query, "--tile1", source]
        if target:
            command += ["--tile2", target]
        target_object = payload.get("target_object")
        if target_object not in (None, ""):
            command += ["--target", str(integer(target_object, "目标物体 ID", 0, 2_147_483_647))]
        if backend == "gpu":
            command.append("-g")
        command += ["--cn", str(integer(payload.get("compute_threads"), "几何计算线程数", 1, 256, 8))]
        # 多个 LOD 用于渐进式过滤；关闭时只在最高层级执行一次。
        if payload.get("progressive_filter", True):
            # 页面启用时使用完整层级；兼容未传开关的旧 API 请求。
            default_lods = "20,40,60,80,100" if payload.get("progressive_filter") else "100"
            lods = str(payload.get("lods", default_lods)).strip()
            if not re.fullmatch(r"(?:20|40|60|80|100)(?:\s*,\s*(?:20|40|60|80|100))*", lods):
                raise ValueError("LOD 只能由 20、40、60、80、100 组成，用逗号分隔")
            for lod in dict.fromkeys(part.strip() for part in lods.split(",")):
                command += ["-l", lod]
        else:
            command += ["-l", "100"]
        if query == "within":
            command += ["-w", str(integer(payload.get("distance"), "距离", 0, 10_000_000, 100))]
        elif query == "nn":
            command += ["-k", str(integer(payload.get("k"), "K", 1, 100_000, 1))]

    command.append("-p")
    return command, binary, {"dimension": dimension, "backend": backend, "query": query}


def timing_seconds(logs: str, stage: str) -> float | None:
    pattern = rf"^\s*{re.escape(stage)}:\s*([0-9]+(?:\.[0-9]+)?(?:e[+-]?[0-9]+)?)\s*(ms|s)\b"
    matches = re.findall(pattern, logs, re.MULTILINE | re.IGNORECASE)
    if not matches:
        return None
    value, unit = matches[-1]
    return float(value) / (1000 if unit.lower() == "ms" else 1)


def legacy_timing_seconds(logs: str, labels: tuple[str, ...]) -> float | None:
    """Read timings emitted by older binaries as ``... takes N ms`` lines."""
    label_pattern = "|".join(re.escape(label) for label in labels)
    pattern = rf"\b(?:{label_pattern})\b[^\n]*?\btakes\s+([0-9]+(?:\.[0-9]+)?(?:e[+-]?[0-9]+)?)\s*(ms|s)\b"
    matches = re.findall(pattern, logs, re.IGNORECASE)
    if not matches:
        return None
    return sum(float(value) / (1000 if unit.lower() == "ms" else 1)
               for value, unit in matches)


def parse_3d_statistics_timings(stdout: str) -> dict[str, float] | None:
    """Use TDBase's final 17-column statistics line when report logs are absent."""
    number = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
    for line in reversed(stdout.splitlines()):
        fields = line.strip().split()
        if len(fields) != 17 or not all(re.fullmatch(number, field) for field in fields):
            continue
        values = [float(field) for field in fields]
        # The final five fields are index, decode, prepare, compute, evaluate (ms).
        return {
            "index": values[12] / 1000,
            "decode": values[13] / 1000,
            "prepare": values[14] / 1000,
            "compute": values[15] / 1000,
        }
    return None


def sum_timings(*values: float | None) -> float | None:
    present = [value for value in values if value is not None]
    return sum(present) if present else None


def preview_results(lines, dimension: str, query: str) -> tuple[str, int, int]:
    number = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
    pattern = r"\d+\s+\d+"
    if dimension == "2d" and query == "intersection":
        pattern += rf"\s+{number}"
    pair_row = re.compile(rf"^\s*{pattern}\s*$")
    first, last, statistics = [], deque(maxlen=PREVIEW_PAIRS), []
    pair_count = 0
    for line in lines:
        if pair_row.fullmatch(line):
            pair_count += 1
            if len(first) < PREVIEW_PAIRS:
                first.append(line)
            else:
                last.append(line)
        elif line.strip():
            # Keep TDBase's 17-column LOD statistics separate from object pairs.
            statistics.append(line)
    omitted = max(0, pair_count - 2 * PREVIEW_PAIRS)
    if omitted:
        output = (f"【前{PREVIEW_PAIRS}对】\n" + "".join(first)
                  + f"\n…… 中间省略 {omitted:,} 对 ……\n\n【后{PREVIEW_PAIRS}对】\n"
                  + "".join(last))
    else:
        output = "".join(first) + "".join(last)
    if statistics:
        output += ("\n" if output else "") + "".join(statistics)
    return output, pair_count, omitted


def _remove_result_file(entry: dict) -> None:
    try:
        Path(entry["path"]).unlink(missing_ok=True)
    except OSError:
        pass


def _cleanup_result_cache() -> None:
    with RESULT_CACHE_LOCK:
        entries = list(RESULT_CACHE.values())
        RESULT_CACHE.clear()
    for entry in entries:
        _remove_result_file(entry)


atexit.register(_cleanup_result_cache)


def _prune_result_cache(now: float | None = None) -> None:
    current = time.time() if now is None else now
    expired = [key for key, entry in RESULT_CACHE.items()
               if current - entry["created"] > RESULT_CACHE_TTL]
    for key in expired:
        _remove_result_file(RESULT_CACHE.pop(key))
    while len(RESULT_CACHE) > RESULT_CACHE_LIMIT:
        _, entry = RESULT_CACHE.popitem(last=False)
        _remove_result_file(entry)


def _cache_result_output(output_file, payload: dict, meta: dict) -> str:
    output_file.flush()
    path = Path(output_file.name)
    token = uuid.uuid4().hex
    entry = {
        "path": str(path),
        "created": time.time(),
        "dimension": meta["dimension"],
        "query": meta["query"],
        "source": absolute_dataset_path(payload.get("source"), "源数据"),
        "target": absolute_dataset_path(payload.get("target"), "目标数据", required=False),
        "distance": payload.get("distance", 0),
    }
    with RESULT_CACHE_LOCK:
        _prune_result_cache()
        RESULT_CACHE[token] = entry
        _prune_result_cache()
    return token


def _get_result_snapshot(token: object) -> dict:
    key = str(token or "")
    with RESULT_CACHE_LOCK:
        _prune_result_cache()
        entry = RESULT_CACHE.get(key)
        if entry is None:
            raise ValueError("查询结果已过期，请重新运行查询后再可视化")
        RESULT_CACHE.move_to_end(key)
        return dict(entry)


def _result_neighbors(entry: dict, selected_id: int) -> tuple[list[int], int]:
    found = set()
    number = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
    with open(entry["path"], "r", encoding="utf-8", errors="replace") as stream:
        for line in stream:
            fields = line.split()
            if len(fields) != (3 if entry["dimension"] == "2d" and entry["query"] == "intersection" else 2):
                continue
            if not fields[0].isdigit() or not fields[1].isdigit():
                continue
            if int(fields[0]) != selected_id:
                continue
            if len(fields) == 3:
                if not re.fullmatch(number, fields[2]) or float(fields[2]) <= 0:
                    continue
            neighbor_id = int(fields[1])
            if neighbor_id != selected_id:
                found.add(neighbor_id)
    ordered = sorted(found)
    return ordered[:MAX_VISUAL_NEIGHBORS], len(ordered)


def _dataset_path(value: str) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = ROOT / path
    path = path.resolve()
    if not path.is_file():
        raise FileNotFoundError(f"找不到可视化数据文件: {path}")
    return path


def _ring_points(data: bytes, offset: int, max_vertices: int) -> tuple[list[list[float]], int, bool]:
    if offset + 8 > len(data):
        raise ValueError("FARM 多边形数据不完整")
    count = struct.unpack_from("<Q", data, offset)[0]
    offset += 8
    if count < 3 or count > (len(data) - offset) // 8:
        raise ValueError("FARM 多边形顶点数量无效")
    end = offset + count * 8
    points = [(float(x), float(y)) for x, y in struct.iter_unpack("<ff", data[offset:end])]
    offset = end
    if len(points) > 1 and points[0] == points[-1]:
        points.pop()
    simplified = len(points) > max_vertices
    if simplified:
        stride = math.ceil(len(points) / max_vertices)
        points = points[::stride]
        if len(points) < 3:
            raise ValueError("多边形顶点过少，无法可视化")
    return [[x, y] for x, y in points], offset, simplified


def _farm_polygon(path: Path, object_id: int) -> dict:
    size = path.stat().st_size
    if size < 8:
        raise ValueError("无效的 FARM 数据文件")
    with path.open("rb") as stream:
        stream.seek(-8, os.SEEK_END)
        count = struct.unpack("<Q", stream.read(8))[0]
        meta_size = struct.calcsize("<IIQ4d")
        footer_size = 8 + count * meta_size
        if count == 0 or footer_size > size or object_id >= count:
            raise ValueError(f"二维对象 ID {object_id} 超出数据范围")
        payload_size = size - footer_size
        stream.seek(payload_size + object_id * meta_size)
        record_size, expected_vertices, offset, *_ = struct.unpack("<IIQ4d", stream.read(meta_size))
        if record_size < 16 or offset + record_size > payload_size:
            raise ValueError("FARM 多边形元数据损坏")
        stream.seek(offset)
        data = stream.read(record_size)
    if len(data) != record_size or record_size < 16:
        raise ValueError("无法读取完整的 FARM 多边形")
    holes = struct.unpack_from("<Q", data, 0)[0]
    rings = []
    cursor = 8
    outer_count = struct.unpack_from("<Q", data, cursor)[0]
    ring, cursor, simplified = _ring_points(data, cursor, MAX_VISUAL_VERTICES)
    rings.append(ring)
    if outer_count != expected_vertices:
        raise ValueError("FARM 多边形顶点元数据不匹配")
    for _ in range(holes):
        ring, cursor, ring_simplified = _ring_points(data, cursor, MAX_VISUAL_VERTICES)
        rings.append(ring)
        simplified = simplified or ring_simplified
    xs = [point[0] for point in rings[0]]
    ys = [point[1] for point in rings[0]]
    center = [(min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2]
    return {"id": object_id, "rings": rings, "center": center, "simplified": simplified}


def _tile_boxes(path: Path, object_ids: set[int]) -> dict[int, dict]:
    if not object_ids:
        return {}
    found = {}
    with path.open("rb") as stream:
        dtype_byte = stream.read(1)
        if not dtype_byte or dtype_byte[0] not in (0, 1):
            raise ValueError("不支持该三维 tile 文件格式")
        dtype = dtype_byte[0]
        object_id = 0
        while object_id <= max(object_ids):
            raw_size = stream.read(8)
            if not raw_size:
                break
            if len(raw_size) != 8:
                raise ValueError("三维 tile 文件结构不完整")
            data_size = struct.unpack("<Q", raw_size)[0]
            stream.seek(data_size, os.SEEK_CUR)
            count_data = stream.read(8)
            if len(count_data) != 8:
                raise ValueError("三维 tile 对象元数据不完整")
            voxel_count = struct.unpack("<Q", count_data)[0]
            if voxel_count > 10_000_000:
                raise ValueError("三维 tile 的空间块数量异常")
            if dtype == 1:
                stream.seek(5 * 8, os.SEEK_CUR)
            if object_id not in object_ids:
                stream.seek(voxel_count * (36 + (80 if dtype == 1 else 0)), os.SEEK_CUR)
                object_id += 1
                continue
            lows = [math.inf, math.inf, math.inf]
            highs = [-math.inf, -math.inf, -math.inf]
            for _ in range(voxel_count):
                bounds = stream.read(36)
                if len(bounds) != 36:
                    raise ValueError("三维 tile 包围盒数据不完整")
                values = struct.unpack("<9f", bounds)
                for axis in range(3):
                    lows[axis] = min(lows[axis], values[axis])
                    highs[axis] = max(highs[axis], values[axis + 3])
                if dtype == 1:
                    stream.seek(5 * 16, os.SEEK_CUR)
            if object_id in object_ids:
                if voxel_count == 0 or not all(math.isfinite(value) for value in lows + highs):
                    raise ValueError(f"三维对象 {object_id} 没有有效空间范围")
                center = [(lows[i] + highs[i]) / 2 for i in range(3)]
                found[object_id] = {
                    "id": object_id,
                    "center": center,
                    "bounds": [lows, highs],
                }
            object_id += 1
    missing = sorted(object_ids.difference(found))
    if missing:
        raise ValueError(f"三维对象 ID 超出数据范围: {missing[0]}")
    return found


def visualize_result(payload: dict) -> dict:
    selected_id = integer(payload.get("id"), "对象 ID", 0, 2_147_483_647)
    entry = _get_result_snapshot(payload.get("result_id"))
    neighbor_ids, total_neighbors = _result_neighbors(entry, selected_id)
    source_path = _dataset_path(entry["source"])
    target_path = _dataset_path(entry["target"] or entry["source"])
    object_ids = set(neighbor_ids) | {selected_id}
    if entry["dimension"] == "2d":
        geometries = {object_id: _farm_polygon(path, object_id)
                      for object_id, path in [(selected_id, source_path)]}
        for object_id in neighbor_ids:
            geometries[object_id] = _farm_polygon(target_path, object_id)
        target = geometries.pop(selected_id)
        neighbors = [geometries[object_id] for object_id in neighbor_ids]
    else:
        if source_path == target_path:
            geometries = _tile_boxes(source_path, object_ids)
            target = geometries.pop(selected_id)
            neighbors = [geometries[object_id] for object_id in neighbor_ids]
        else:
            target = _tile_boxes(source_path, {selected_id})[selected_id]
            neighbors = list(_tile_boxes(target_path, set(neighbor_ids)).values()) if neighbor_ids else []
    distance = None
    if entry["query"] == "within":
        distance = integer(entry["distance"], "距离", 0, 10_000_000, 0)
    return {
        "dimension": entry["dimension"],
        "query": entry["query"],
        "target": target,
        "neighbors": neighbors,
        "total_neighbors": total_neighbors,
        "shown_neighbors": len(neighbors),
        "truncated": total_neighbors > len(neighbors),
        "distance": distance,
        "distance_unit": "km" if entry["dimension"] == "2d" else "三维单位",
        "simplified": bool(target.get("simplified")) or any(item.get("simplified") for item in neighbors),
    }


def run_query(payload: dict) -> dict:
    command, binary, meta = build_command(payload)
    started = time.monotonic()
    output_file_path = None
    try:
        # Retain only complete preview rows in memory, even for millions of pairs.
        with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8", errors="replace", delete=False,
                                         prefix="tdbase-visual-", suffix=".txt") as output_file, \
                tempfile.TemporaryFile(mode="w+", encoding="utf-8", errors="replace") as log_file:
            output_file_path = Path(output_file.name)
            process = subprocess.run(
                command,
                cwd=ROOT,
                stdout=output_file,
                stderr=log_file,
                timeout=ARGS.timeout,
                env=os.environ.copy(),
            )
            elapsed = round(time.monotonic() - started, 3)
            output_file.seek(0)
            stdout, output_pair_count, output_omitted_pairs = preview_results(
                output_file, meta["dimension"], meta["query"])
            result_id = _cache_result_output(output_file, payload, meta)
            output_file_path = None
            log_file.seek(0)
            logs = log_file.read()
    except subprocess.TimeoutExpired as error:
        if output_file_path:
            output_file_path.unlink(missing_ok=True)
        raise TimeoutError(f"查询超过 {ARGS.timeout} 秒，已终止") from error
    except Exception:
        if output_file_path:
            output_file_path.unlink(missing_ok=True)
        raise
    stderr = logs[-MAX_LOG_OUTPUT:]
    count_pattern = r"count-found:\s*(\d+)" if meta["dimension"] == "2d" else r"(?m)^\s*results:\s*(\d+)\s*\("
    count_match = re.findall(count_pattern, logs)
    statistics_timings = parse_3d_statistics_timings(stdout) if meta["dimension"] == "3d" else None
    if meta["dimension"] == "2d":
        preprocess_seconds = timing_seconds(logs, "load")
        if preprocess_seconds is None:
            # Older FARM binaries only expose per-file ``loading ... takes`` logs.
            preprocess_seconds = legacy_timing_seconds(logs, ("loading", "loaded"))
        rasterization_seconds = timing_seconds(logs, "rasterization")
    else:
        decode_seconds = timing_seconds(logs, "decode")
        prepare_seconds = timing_seconds(logs, "prepare")
        if statistics_timings:
            decode_seconds = decode_seconds if decode_seconds is not None else statistics_timings["decode"]
            prepare_seconds = prepare_seconds if prepare_seconds is not None else statistics_timings["prepare"]
        preprocess_seconds = sum_timings(decode_seconds, prepare_seconds)
        rasterization_seconds = None
    index_seconds = timing_seconds(logs, "index")
    compute_seconds = timing_seconds(logs, "compute")
    if statistics_timings:
        index_seconds = index_seconds if index_seconds is not None else statistics_timings["index"]
        compute_seconds = compute_seconds if compute_seconds is not None else statistics_timings["compute"]
    return {
        **meta,
        "ok": process.returncode == 0,
        "exit_code": process.returncode,
        "duration_seconds": elapsed,
        "total_seconds": timing_seconds(logs, "total"),
        "preprocess_seconds": preprocess_seconds,
        "index_seconds": index_seconds,
        "rasterization_seconds": rasterization_seconds,
        "compute_seconds": compute_seconds,
        "command": command,
        "command_text": shlex.join(command),
        "binary": str(binary),
        "count": int(count_match[-1]) if count_match else None,
        "output_pair_count": output_pair_count,
        "output_omitted_pairs": output_omitted_pairs,
        "result_id": result_id,
        "stdout": stdout,
        "stderr": stderr,
    }


class Handler(BaseHTTPRequestHandler):
    server_version = "TDBaseWeb/1.0"

    def log_message(self, format: str, *args: object) -> None:
        print(f"[{self.log_date_time_string()}] {format % args}")

    def send_json(self, value: dict, status: int = HTTPStatus.OK) -> None:
        data = json.dumps(value, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self) -> None:
        parsed = urlparse(self.path)
        if parsed.path == "/api/health":
            self.send_json({
                "ok": True,
                "gpu_binary": str(executable("gpu")),
                "gpu_available": executable("gpu").is_file(),
                "cpu_binary": str(executable("cpu")),
                "cpu_available": executable("cpu").is_file(),
            })
            return
        if parsed.path in {"/", "/index.html"}:
            self.serve_file(FRONTEND / "index.html", "text/html; charset=utf-8")
            return
        if parsed.path in {"/app.js", "/styles.css", "/visualization.js",
                           "/vendor/three.module.js", "/vendor/three.core.js",
                           "/vendor/OrbitControls.js"}:
            mime = "text/javascript; charset=utf-8" if parsed.path.endswith("js") else "text/css; charset=utf-8"
            self.serve_file(FRONTEND / parsed.path.lstrip("/"), mime)
            return
        self.send_json({"error": "Not found"}, HTTPStatus.NOT_FOUND)

    def serve_file(self, path: Path, mime: str) -> None:
        try:
            data = path.read_bytes()
            if path == FRONTEND / "index.html":
                data = data.decode("utf-8").replace(
                    'data-project-root=""',
                    f'data-project-root="{escape(str(ROOT), quote=True)}"',
                ).encode("utf-8")
        except FileNotFoundError:
            self.send_json({"error": "Frontend file not found"}, HTTPStatus.NOT_FOUND)
            return
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", mime)
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def do_POST(self) -> None:
        endpoint = urlparse(self.path).path
        if endpoint not in {"/api/query", "/api/visualize", "/api/dataset-stats"}:
            self.send_json({"error": "Not found"}, HTTPStatus.NOT_FOUND)
            return
        try:
            size = int(self.headers.get("Content-Length", "0"))
            if size > 100_000:
                raise ValueError("请求过大")
            payload = json.loads(self.rfile.read(size).decode("utf-8"))
            if not isinstance(payload, dict):
                raise ValueError("请求内容必须是 JSON 对象")
            if endpoint == "/api/dataset-stats":
                self.send_json({"ok": True, **DATASET_STATS.get(payload.get("dimension"), payload.get("path"))})
            elif endpoint == "/api/visualize":
                self.send_json(visualize_result(payload))
            else:
                result = run_query(payload)
                self.send_json(result, HTTPStatus.OK if result["ok"] else HTTPStatus.BAD_REQUEST)
        except (ValueError, FileNotFoundError, TimeoutError, json.JSONDecodeError) as error:
            self.send_json({"ok": False, "error": str(error)}, HTTPStatus.BAD_REQUEST)
        except struct.error as error:
            self.send_json({"ok": False, "error": f"数据文件结构无效: {error}"}, HTTPStatus.BAD_REQUEST)
        except Exception as error:  # Keep the local UI responsive on process errors.
            self.send_json({"ok": False, "error": f"服务器错误: {error}"}, HTTPStatus.INTERNAL_SERVER_ERROR)


def main() -> None:
    print(f"TDBase Web UI: http://{ARGS.host}:{ARGS.port}")
    print(f"Project root: {ROOT}")
    if BUILD_OVERRIDE:
        print(f"Build override: {BUILD_OVERRIDE}")
    if GPU_BUILD_OVERRIDE:
        print(f"GPU build override: {GPU_BUILD_OVERRIDE}")
    if CPU_BUILD_OVERRIDE:
        print(f"CPU build override: {CPU_BUILD_OVERRIDE}")
    server = ThreadingHTTPServer((ARGS.host, ARGS.port), Handler)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping TDBase Web UI")
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
