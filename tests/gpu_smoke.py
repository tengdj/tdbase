"""GPU/CPU CLI comparisons: python3 gpu_smoke.py BIN_DIR SOURCE_DIR [LOG_DIR]."""
import json
import os
from pathlib import Path
import re
import struct
import subprocess
import sys

build = Path(sys.argv[1]).resolve()
source = Path(sys.argv[2]).resolve()
work = Path(sys.argv[3]).resolve() if len(sys.argv) > 3 else build / 'gpu-validation'
work.mkdir(parents=True, exist_ok=True)
env = dict(os.environ, OMP_NUM_THREADS='2', CUDA_LAUNCH_BLOCKING='1')
checks = []


def run(label, command, expected_count=None, require_gpu=False):
    try:
        process = subprocess.run(list(map(str, command)), cwd=work, env=env,
                                 text=True, capture_output=True, timeout=120)
        output = process.stdout + process.stderr
        good = process.returncode == 0
        if expected_count is not None:
            counts = re.findall(r'count-found:\s*(\d+)', output)
            good = good and bool(counts) and int(counts[-1]) == expected_count
        if require_gpu:
            good = good and 'GPU Memory Allocation' in output
        good = good and not re.search(r'Cuda error|launch failed|illegal memory access', output)
        code = process.returncode
    except subprocess.TimeoutExpired as exc:
        output = 'TIMEOUT\n' + str(exc)
        good, code = False, None
    (work / (label + '.log')).write_text(output)
    checks.append({'name': label, 'passed': bool(good), 'exit_code': code})
    print(('PASS ' if good else 'FAIL ') + label, flush=True)
    return output


def rect(x0, y0, x1, y1):
    return f'POLYGON (({x0} {y0},{x1} {y0},{x1} {y1},{x0} {y1},{x0} {y0}))'


fixtures = {
    'a': rect(0, 0, 1, 1),
    'overlap': rect(.5, .5, 1.5, 1.5),
    'inside': rect(.25, .25, .75, .75),
    'far': rect(5, 5, 6, 6),
    'near': rect(1.01, 0, 2.01, 1),
    'self': rect(0, 0, 1, 1) + '\n' + rect(1.01, 0, 2.01, 1),
}
for name, wkt in fixtures.items():
    (work / (name + '.wkt')).write_text(wkt + '\n')
    run('convert_' + name, [build / 'dump_polygons', work / (name + '.wkt'), work / (name + '.dt')])

for name in ('intersect', 'intersection', 'within_polygon'):
    run('help_' + name, [build / name, '--help'])

for bitwidth in (2, 4, 8):
    for hierarchy in (False, True):
        opts = ['-s', work / 'a.dt', '-n', '2', '-v', '1', '-w', str(bitwidth)]
        if hierarchy:
            opts += ['-h']
        suffix = f'_w{bitwidth}_h{int(hierarchy)}'
        for name, expected in (('overlap', 1), ('inside', 1), ('far', 0)):
            run('intersect_' + name + suffix,
                [build / 'intersect', *opts, '-t', work / (name + '.dt')],
                expected_count=expected, require_gpu=name != 'far')
        # The area CLI does not expose numeric areas; this checks CUDA execution.
        run('area_overlap' + suffix,
            [build / 'intersection', *opts, '-t', work / 'overlap.dt'], require_gpu=True)
        for distance, expected in ((1, 0), (2, 1)):
            run(f'within_d{distance}' + suffix,
                [build / 'within_polygon', *opts, '-t', work / 'near.dt', '-d', str(distance)],
                expected_count=expected, require_gpu=distance == 2)

run('within_self', [build / 'within_polygon', '-s', work / 'self.dt', '-h',
                    '-v', '1', '-n', '2', '-d', '2'], expected_count=2, require_gpu=True)
for name in ('intersect', 'intersection', 'within_polygon'):
    run('approx_' + name, [build / name, '-s', work / 'a.dt', '-t', work / 'overlap.dt',
                          '-h', '-a', '-v', '1', '-n', '2'], require_gpu=True)


def first_voxels(path):
    with path.open('rb') as stream:
        assert stream.read(1) == b'\0'
        stream.seek(struct.unpack('<Q', stream.read(8))[0], 1)
        return [struct.unpack('<9f', stream.read(36))
                for _ in range(struct.unpack('<Q', stream.read(8))[0])]


tiles = [source / 'data/foo_n_nv50_nu200_vs100_r30.dt',
         source / 'data/foo_v_nv50_nu200_vs100_r30.dt']
low, high = float('inf'), float('inf')
for a in first_voxels(tiles[0]):
    for b in first_voxels(tiles[1]):
        low = min(low, sum(max(a[i] - b[i+3], b[i] - a[i+3], 0)**2 for i in range(3))**.5)
        high = min(high, sum(((a[i] + a[i+3] - b[i] - b[i+3]) / 2)**2 for i in range(3))**.5)
base = [build / 'tdbase', 'join', '--tile1', tiles[0], '--tile2', tiles[1],
        '--max_objects1', '1', '--max_objects2', '1', '--cn', '2', '-p', '-v', '2']
for label, lods in (('progressive', ['20', '40', '60', '80', '100']), ('lod100', ['100'])):
    args = ['-q', 'within', '-w', str((low + high) / 2)]
    for lod in lods:
        args += ['-l', lod]
    cpu = run('tdbase_cpu_' + label, [*base, *args])
    gpu = run('tdbase_gpu_' + label, [*base, *args, '-g'])
    pairs = lambda output: set(re.findall(r'^(\d+)\t(\d+)$', output, re.MULTILINE))
    # Allocation and the GPU computation log prove the geometry path was used.
    good = (pairs(cpu) == pairs(gpu) and 'allocated for GPU' in gpu
            and 'started to get distance' in gpu and 'GPU 0' in gpu
            and 'decode data to' in gpu)
    checks.append({'name': 'tdbase_compare_' + label, 'passed': good})
    print(('PASS ' if good else 'FAIL ') + 'tdbase_compare_' + label, flush=True)

summary = {'build': str(build), 'checks': checks,
           'passed': sum(c['passed'] for c in checks), 'total': len(checks)}
(work / 'results.json').write_text(json.dumps(summary, indent=2))
print(f"{summary['passed']}/{summary['total']} checks passed", flush=True)
sys.exit(0 if all(c['passed'] for c in checks) else 1)
