#!/usr/bin/env python3
"""Rebuild the shared backends and reevaluate the n=163 Table 6 configurations."""
import argparse
from contextlib import nullcontext
import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
METRICS = ('toffoli', 'cnot', 'width', 'toffoli_depth', 'nct_depth')


def run(command, cwd=ROOT, stdout=None):
    print('Running:', ' '.join(map(str, command)), flush=True)
    with (open(stdout, 'w') if stdout else nullcontext(None)) as stream:
        subprocess.run(list(map(str, command)), cwd=cwd, stdout=stream, check=True)


def summarize(output):
    rows = {}
    point = json.loads((output/'ac_inplace_point.json').read_text())[0]
    rows['In-place+AC'] = (point['toffoli']*328, point['cnot']*328, point['qubits'],
                           point['toffoli_depth']*328, point['full_depth']*328)
    scan = list(csv.DictReader((output/'ac_balanced_scan.csv').open()))
    if {int(row['w']) for row in scan} != set(range(1, 129)):
        raise ValueError('AC Balanced scan must cover every w from 1 through 128')
    best = min(scan, key=lambda row: int(row['dw_full']))
    rows['Balanced+AC'] = tuple(int(best[k]) for k in
        ('toffoli', 'cnot', 'qubits', 'toffoli_depth', 'full_depth'))
    kara = next(csv.DictReader((output/'karatsuba_w30.csv').open()))
    rows['Balanced+Karatsuba'] = tuple(int(kara[k]) for k in
        ('Toffoli', 'CNOT', 'Qubits', 'Toffoli_Depth', 'Full_Depth'))
    text = (output/'retained_full.txt').read_text()
    def number(label):
        return int(re.search(r'^'+re.escape(label)+r'\s*\|\s*(\d+)', text, re.M)[1])
    rows['T-optimal+AC'] = tuple(number(label) for label in
        ('Toffoli Gates', 'CNOT Gates', 'Qubits', 'Toffoli Depth', 'Full Depth'))
    window = json.loads((output/'window_tdw/summary.json').read_text())['fields']['163']
    if not window['formula_minimizers_all_available']:
        raise ValueError('Window formula minimizers are incomplete')
    selected = window['selected_optimum']
    rows['Window+AC'] = tuple(selected[k] for k in
        ('toffoli', 'cnot', 'width', 'toffoli_depth', 'nct_depth_barrier'))
    for row in json.loads((output/'jsb25_shor.json').read_text()):
        name = 'JSB25 in-place' if row['mode']=='inplace' else 'JSB25 out-of-place'
        rows[name] = tuple(row[k] for k in
            ('toffoli_count', 'cnot_count', 'qubits', 'toffoli_depth', 'full_depth'))
    resources = []
    for name, actual in rows.items():
        row = {'configuration': name, 'n': 163}
        row.update(dict(zip(METRICS, actual)))
        row['dw'] = actual[2]*actual[4]
        row['tdw'] = actual[2]*actual[3]
        resources.append(row)
    fields = list(resources[0])
    with (output/'resources.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader(); writer.writerows(resources)
    summary = {
        'backend': 'nct-v2',
        'scope': 'controlled point-addition stage including the configured outer cleanup',
        'ac_balanced_scan_range': [1,128], 'ac_balanced_best_w': int(best['w']),
        'window_selected': {'s': selected['s'], 'w': selected['w']},
        'resources': resources,
        'backend_sha256': {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
            for directory in ['cpp','python'] for p in sorted((ROOT/'circuit_backend'/directory).glob('*'))
            if p.is_file()},
    }
    (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True,
                        help='a new directory for resource records, logs, and temporary builds')
    parser.add_argument('--summarize-only', action='store_true')
    args = parser.parse_args()
    out = args.output.resolve()
    if args.summarize_only:
        summarize(out); return
    if out.exists() and any(out.iterdir()):
        parser.error('use an empty output directory; existing records are never overwritten')
    out.mkdir(parents=True, exist_ok=True)
    binaries = {}
    for name, directory, binary in [
        ('balanced', 'mul_line/ac_balanced_dw_search', 'ac_balanced_estimator'),
        ('karatsuba', 'recursive_karatsuba/C++', 'karatsuba_balanced_estimator'),
        ('retained', 'mul_line/C++', 'ShorEstimator'),
    ]:
        build = out/'build'/name
        run(['cmake','-S',ROOT/directory,'-B',build,'-DCMAKE_BUILD_TYPE=Release'])
        run(['cmake','--build',build,'-j','2'])
        binaries[name] = build/binary
    run([binaries['balanced'],'scan','--n','163','--w-start','1','--w-end','128',
         '--data-root',ROOT/'mul_line/C++/data','--layer-cache',out/'ac_balanced_layers.csv',
         '--output-csv',out/'ac_balanced_scan.csv'], stdout=out/'ac_balanced_scan.log')
    run([sys.executable,'resources/point_addition_resources.py','--sizes','163',
         '--inversion','optimal_depth','--json',out/'ac_inplace_point.json'],
        ROOT/'mul_line/algorithm1_inplace_point_add')
    run([sys.executable,'shor_table_resources.py','--n','163','--mode','both',
         '--format','json','--output',out/'jsb25_shor.json'],ROOT/'jsb25')
    run([binaries['karatsuba'],'--n','163','--w','30','--output',out/'karatsuba_w30.csv'])
    run([binaries['retained'],'--n','163'], stdout=out/'retained_full.txt')
    run([sys.executable,'run.py','run','--n','163','--objective','tdw','--scope','winners',
         '--jobs','1','--memory-budget-gb','6','--verify','--output',out/'window_tdw'],
        ROOT/'windowed_qrom_toptimal',out/'window.log')
    summarize(out)


if __name__ == '__main__':
    main()
