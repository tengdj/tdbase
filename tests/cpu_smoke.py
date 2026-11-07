"""Run CPU CLI smoke checks: python3 cpu_smoke.py BIN_DIR SOURCE_DIR [LOG_DIR]."""
import json
import os
from pathlib import Path
import re
import struct
import subprocess
import sys

build = Path(sys.argv[1]).resolve()
source = Path(sys.argv[2]).resolve()
work = Path(sys.argv[3]).resolve() if len(sys.argv) > 3 else build / "cpu-validation"
work.mkdir(parents=True, exist_ok=True)
env = dict(os.environ, OMP_NUM_THREADS="2")
env.pop("USE_GPU", None)
checks = []


def run(label, command, expected_count=None, expected_pair=None, require_decode=False):
    process = subprocess.run(
        list(map(str, command)), cwd=work, env=env,
        text=True, capture_output=True, timeout=90,
    )
    output = process.stdout + process.stderr
    (work / (label + ".log")).write_text(output)
    good = process.returncode == 0
    if expected_count is not None:
        counts = re.findall(r"count-found:\s*(\d+)", output)
        good = good and bool(counts) and int(counts[-1]) == expected_count
    if expected_pair is not None:
        pairs = set(re.findall(r"^(\d+)\t(\d+)$", process.stdout, re.MULTILINE))
        good = good and expected_pair in pairs
    if require_decode:
        good = good and "decode data to" in output and "distance computation" in output
    checks.append({"name": label, "passed": good, "exit_code": process.returncode})
    print(("PASS " if good else "FAIL ") + label, flush=True)
    return output


def rect(x0, y0, x1, y1):
    return f"POLYGON (({x0} {y0},{x1} {y0},{x1} {y1},{x0} {y1},{x0} {y0}))"


inputs = {
    "a": rect(0, 0, 1, 1),
    "overlap": rect(.5, .5, 1.5, 1.5),
    "inside": rect(.25, .25, .75, .75),
    "far": rect(5, 5, 6, 6),
    "near": rect(1.01, 0, 2.01, 1),
    "self": rect(0, 0, 1, 1) + "\n" + rect(1.01, 0, 2.01, 1),
}
for name, wkt in inputs.items():
    (work / (name + ".wkt")).write_text(wkt + "\n")
    run("convert_" + name, [build / "dump_polygons", work / (name + ".wkt"), work / (name + ".dt")])

for name in ("intersect_cpu", "intersection_cpu", "within_polygon_cpu"):
    run("help_" + name, [build / name, "--help"])

for bitwidth in (2, 4, 8):
    for hierarchical in (False, True):
        opts = ["-s", work / "a.dt", "-n", "2", "-v", "1", "-w", str(bitwidth)]
        if hierarchical:
            opts += ["-h"]
        suffix = f"_w{bitwidth}_h{int(hierarchical)}"
        for name, expected in (("overlap", 1), ("inside", 1), ("far", 0)):
            run("intersect_" + name + suffix,
                [build / "intersect_cpu", *opts, "-t", work / (name + ".dt")],
                expected_count=expected)
        run("area_overlap" + suffix,
            [build / "intersection_cpu", *opts, "-t", work / "overlap.dt"])
        # Coordinates are geographic; the horizontal gap is about 1.11 km.
        for threshold, expected in ((1, 0), (2, 1)):
            run(f"within_d{threshold}" + suffix,
                [build / "within_polygon_cpu", *opts, "-t", work / "near.dt", "-d", str(threshold)],
                expected_count=expected)

run("within_self",
    [build / "within_polygon_cpu", "-s", work / "self.dt", "-h", "-v", "1", "-n", "2", "-d", "2"],
    expected_count=2)
for name in ("intersect_cpu", "intersection_cpu", "within_polygon_cpu"):
    run("approx_" + name, [build / name, "-s", work / "a.dt", "-t", work / "overlap.dt", "-h", "-a", "-v", "1", "-n", "2"])

run("tdbase_dispatch", [build / "tdbase"])
tiles = [source / "data/foo_n_nv50_nu200_vs100_r30.dt", source / "data/foo_v_nv50_nu200_vs100_r30.dt"]
base = [build / "tdbase", "join", "--tile1", tiles[0], "--tile2", tiles[1],
        "--max_objects1", "1", "--max_objects2", "1", "--cn", "2", "-p"]
run("tdbase_nn", [*base, "-q", "nn", "-k", "1", "-l", "100"], expected_pair=("0", "0"))


def first_voxels(path):
    with path.open("rb") as stream:
        assert stream.read(1) == b"\0", "expected bundled compressed tile"
        payload_size = struct.unpack("<Q", stream.read(8))[0]
        stream.seek(payload_size, 1)
        count = struct.unpack("<Q", stream.read(8))[0]
        return [struct.unpack("<9f", stream.read(36)) for _ in range(count)]


# Choose a threshold between a voxel pair's box bounds so this check also
# enters the decode/packing/computation path, instead of stopping at indexing.
min_low = float("inf")
min_high = float("inf")
for a in first_voxels(tiles[0]):
    for b in first_voxels(tiles[1]):
        low = sum(max(a[i] - b[i+3], b[i] - a[i+3], 0)**2 for i in range(3))**.5
        high = sum(((a[i] + a[i+3] - b[i] - b[i+3]) / 2)**2 for i in range(3))**.5
        min_low = min(min_low, low)
        min_high = min(min_high, high)
threshold = (min_low + min_high) / 2
run("tdbase_within_progressive", [*base, "-q", "within", "-w", str(threshold),
    "-l", "20", "-l", "40", "-l", "60", "-l", "80", "-l", "100"], require_decode=True)
run("tdbase_within_lod100", [*base, "-q", "within", "-w", str(threshold), "-l", "100"], require_decode=True)

summary = {"build": str(build), "checks": checks,
           "passed": sum(c["passed"] for c in checks), "total": len(checks)}
(work / "results.json").write_text(json.dumps(summary, indent=2))
print(f"{summary['passed']}/{summary['total']} checks passed", flush=True)
sys.exit(0 if all(c["passed"] for c in checks) else 1)
