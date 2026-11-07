"""Validate the integrated executable, including running without FARM helper binaries.

python3 tests/unified_cli.py BIN_DIR SOURCE_DIR [LOG_DIR] [--gpu]
"""
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

build, source = (Path(x).resolve() for x in sys.argv[1:3])
work = Path(sys.argv[3]).resolve() if len(sys.argv) > 3 else build / 'unified-validation'
work.mkdir(parents=True, exist_ok=True)
isolated = work / 'binary-only'
isolated.mkdir(exist_ok=True)
binary = isolated / 'tdbase'
shutil.copy2(build / 'tdbase', binary)
environment = dict(os.environ, OMP_NUM_THREADS='2', CUDA_LAUNCH_BLOCKING='1', PATH='/usr/bin:/bin')
checks = []


def check(label, args, expected=0, count=None, rows=None, area=None, json_output=False):
    process = subprocess.run([str(binary), *map(str, args)], cwd=work, env=environment,
                             capture_output=True, text=True, timeout=120)
    (work / (label + '.log')).write_text(process.stdout + process.stderr)
    good = process.returncode == expected
    if count is not None:
        found = re.findall(r'count-found:\s*(\d+)', process.stderr)
        good &= bool(found) and int(found[-1]) == count
    if rows is not None or area is not None:
        values = [line.split('\t') for line in process.stdout.strip().splitlines() if line]
        if rows is not None:
            good &= {(row[0], row[1]) for row in values} == rows
        if area is not None:
            good &= len(values) == 1 and len(values[0]) == 3 and math.isclose(float(values[0][2]), area, abs_tol=1e-6)
    if json_output:
        try:
            result = json.loads(process.stdout)
            good &= result['count'] == count and bool(result['matches'])
        except (ValueError, KeyError):
            good = False
    checks.append({'name': label, 'passed': bool(good), 'exit_code': process.returncode})
    print(('PASS ' if good else 'FAIL ') + label, flush=True)
    return process


def rect(x0, y0, x1, y1):
    return f'POLYGON (({x0} {y0},{x1} {y0},{x1} {y1},{x0} {y1},{x0} {y0}))'


fixtures = {
    'a': rect(0, 0, 1, 1),
    'overlap': rect(.5, .5, 1.5, 1.5),
    'inside': rect(.25, .25, .75, .75),
    'far': rect(5, 5, 6, 6),
    'near': rect(1.01, 0, 2.01, 1),
    'self': rect(0, 0, 1, 1) + '\n' + rect(1.01, 0, 2.01, 1),
    'multi': rect(5, 5, 6, 6) + '\n' + rect(0, 0, 1, 1),
    'exponent': 'POLYGON ((+0E0 0,1E0 0,1 1E0,0 1,+0E0 0))',
}
check('help', ['--help'])
check('farm_help', ['farm', '--help'])
check('farm_alias', ['2d', '--help'])
for name, wkt in fixtures.items():
    (work / (name + '.wkt')).write_text(wkt + '\n')
    check('convert_' + name, ['farm', 'convert', work / (name + '.wkt'), work / (name + '.dt')])

base = ['-s', work / 'a.dt', '-t', work / 'overlap.dt', '-n', '2', '-v', '1', '-p']
backends = [('cpu', [])]
if '--gpu' in sys.argv:
    backends.append(('gpu', ['-g']))
for backend, device in backends:
    for hierarchy in (False, True):
        for bits in (2, 4, 8):
            suffix = f'{backend}_h{int(hierarchy)}_w{bits}'
            opts = [*base, *device, '-w', str(bits), *(['-h'] if hierarchy else [])]
            check('intersect_' + suffix, ['farm', 'intersect', *opts], count=1, rows={('0', '0')})
            check('area_' + suffix, ['farm', 'intersection', *opts], count=1, area=.25)
            opts = ['-s', work / 'a.dt', '-t', work / 'near.dt', '-n', '2', '-v', '1',
                    '-w', str(bits), '-p', *device, *(['-h'] if hierarchy else [])]
            check('within_' + suffix, ['farm', 'within', *opts, '-d', '2'], count=1, rows={('0', '0')})
    check('area_containment_' + backend, ['farm', 'intersection', '-s', work / 'a.dt', '-t', work / 'inside.dt',
                                        '-h', '-p', '-n', '2', '-v', '1', *device], count=1, area=.25)
    check('disjoint_' + backend, ['farm', 'intersect', '-s', work / 'a.dt', '-t', work / 'far.dt',
                                 '-p', '-n', '2', *device], count=0, rows=set())
    check('self_' + backend, ['farm', 'within', '-s', work / 'self.dt', '-d', '2', '-h',
                             '-p', '-n', '2', '-v', '1', *device], count=2, rows={('0', '1'), ('1', '0')})
    check('ids_' + backend, ['farm', 'intersect', '-s', work / 'multi.dt', '-t', work / 'overlap.dt',
                            '-p', '-h', '-n', '2', '-v', '1', *device], count=1, rows={('1', '0')})
    check('join_' + backend, ['join', '--dim', '2', '-q', 'intersect', '--tile1', work / 'a.dt',
                             '--tile2', work / 'overlap.dt', '-p', '-n', '2', *device], count=1, rows={('0', '0')})
    check('join_equal_' + backend, ['join', '--dim=2', '-q', 'within', '--tile1', work / 'a.dt',
                                   '--tile2', work / 'near.dt', '--within_dist', '2', '-p', '-n', '2', *device],
          count=1, rows={('0', '0')})
    for action in ('intersect', 'intersection', 'within'):
        check('approx_' + action + '_' + backend, ['farm', action, *base, '-h', '-a', *device])
    check('json_' + backend, ['farm', 'intersect', *base, '--format', 'json', *device], count=1, json_output=True)

output = work / 'areas.tsv'
check('area_file', ['farm', 'intersection', *base[:-1], '-o', output])
checks.append({'name': 'area_file_value', 'passed': math.isclose(float(output.read_text().split()[2]), .25, abs_tol=1e-6)})
check('exponent_query', ['farm', 'intersect', '-s', work / 'exponent.dt', '-t', work / 'overlap.dt', '-p', '-n', '2'],
      count=1, rows={('0', '0')})
# Cross FARM's 10 MiB load-task boundary, where workers can finish out of order.
large_count = 140000
(work / 'chunked.wkt').write_text((rect(5, 5, 6, 6) + '\n') * large_count + rect(0, 0, 1, 1) + '\n')
check('chunked_convert', ['farm', 'convert', work / 'chunked.wkt', work / 'chunked.dt'])
check('chunked_ids', ['farm', 'intersect', '-s', work / 'chunked.dt', '-t', work / 'overlap.dt', '-p', '-n', '2'],
      count=1, rows={(str(large_count), '0')})
(work / 'filter.wkt').write_text(rect(0, 0, 1, 1) + '\n' + rect(200, 0, 201, 1) + '\n')
check('filter', ['farm', 'filter', work / 'filter.wkt', work / 'filtered.dt'])
check('filter_result', ['farm', 'intersect', '-s', work / 'filtered.dt', '-t', work / 'overlap.dt', '-n', '2'], count=1)
(work / 'empty.wkt').write_text('')
check('empty_convert', ['farm', 'convert', work / 'empty.wkt', work / 'empty.dt'])
check('empty_query', ['farm', 'intersect', '-s', work / 'empty.dt', '-t', work / 'a.dt'], count=0)
(work / 'bad.wkt').write_text('POLYGON ((0 0,1 0\n')
check('bad_wkt', ['farm', 'convert', work / 'bad.wkt', work / 'bad.dt'], expected=1)
(work / 'hole.wkt').write_text('POLYGON ((0 0,2 0,2 2,0 2,0 0),(0.5 0.5,1 0.5,1 1,0.5 1,0.5 0.5))\n')
check('unsupported_inner_ring', ['farm', 'convert', work / 'hole.wkt', work / 'hole.dt'], expected=1)
check('missing_input', ['farm', 'intersect', '-s', work / 'missing.dt', '-t', work / 'a.dt'], expected=1)
check('missing_target', ['farm', 'intersect', '-s', work / 'a.dt'], expected=1)
check('wrong_dimension_data', ['farm', 'intersect', '-s', source / 'data/foo_n_nv50_nu200_vs100_r30.dt',
                              '-t', work / 'a.dt'], expected=1)
check('bad_bitwidth', ['farm', 'intersect', *base, '-w', '9'], expected=1)
check('bad_threads', ['farm', 'intersect', *base[:-3], '-n', '0'], expected=1)
check('bad_distance', ['farm', 'within', *base, '-d', '-1'], expected=1)
check('bad_query', ['join', '--dim', '2', '-q', 'nn'], expected=1)
check('bad_dimension', ['join', '--dim', '4'], expected=1)
check('missing_dimension', ['join', '--dim'], expected=1)
check('duplicate_dimension', ['join', '--dim', '2', '--dim', '3'], expected=1)
tile1 = source / 'data/foo_n_nv50_nu200_vs100_r30.dt'
tile2 = source / 'data/foo_v_nv50_nu200_vs100_r30.dt'
check('explicit_3d', ['join', '--dim', '3', '--tile1', tile1, '--tile2', tile2,
                      '--max_objects1', '1', '--max_objects2', '1', '--cn', '2', '-q', 'nn', '-k', '1', '-l', '100', '-p'],
      rows={('0', '0')})
check('unknown', ['farm', 'not_a_command'], expected=1)
check('output_same_input', ['farm', 'intersect', *base, '-o', work / './a.dt'], expected=1)
if '--gpu' not in sys.argv:
    check('gpu_in_cpu_build', ['farm', 'intersect', *base, '-g'], expected=1)

summary = {'build': str(build), 'checks': checks,
           'passed': sum(c['passed'] for c in checks), 'total': len(checks)}
(work / 'results.json').write_text(json.dumps(summary, indent=2))
print(f"{summary['passed']}/{summary['total']} checks passed", flush=True)
sys.exit(0 if all(c['passed'] for c in checks) else 1)
