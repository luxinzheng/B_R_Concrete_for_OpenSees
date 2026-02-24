# -*- coding: utf-8 -*-
"""Test old exe (02210959) with the NEW fixed TCL (no ModifiedNewton, NewtonLineSearch)."""
import subprocess, os, time

root = r'e:\Basic\Concrete_Model\ADINA'
exe = os.path.join(root, 'OpenSees-02210959.exe')
d = os.path.join(root, 'Multi-layer_Shell')

for fn in ['test_disp.txt', 'test_react.txt', 'ml_log.txt']:
    fp = os.path.join(d, fn)
    if os.path.isfile(fp):
        os.remove(fp)

print(f'Testing OLD exe: {os.path.basename(exe)}', flush=True)
print(f'With NEW test_br.tcl (no ModifiedNewton, NewtonLineSearch)', flush=True)
t0 = time.time()
try:
    log_fp = os.path.join(d, 'ml_log.txt')
    with open(log_fp, 'w') as logf:
        p = subprocess.run([exe, 'test_br.tcl'], cwd=d,
                           stdout=logf, stderr=subprocess.STDOUT, timeout=1800)
    elapsed = time.time() - t0
    print(f'Finished in {elapsed:.1f}s, exit={p.returncode}', flush=True)
except subprocess.TimeoutExpired:
    elapsed = time.time() - t0
    print(f'TIMEOUT after {elapsed:.1f}s', flush=True)

# Analyze
if os.path.isfile(os.path.join(d, 'ml_log.txt')):
    with open(os.path.join(d, 'ml_log.txt'), 'r', errors='replace') as f:
        lines = f.readlines()
    nan_c = sum(1 for l in lines if '-1.#IND' in l)
    warn_c = sum(1 for l in lines if 'failed to converge' in l)
    blocks = [l.rstrip() for l in lines if 'Completed block' in l]
    print(f'\nLog: {len(lines)} lines')
    print(f'NaN (-1.#IND): {nan_c}')
    print(f'Convergence warnings: {warn_c}')
    print(f'Completed blocks: {len(blocks)}')
    for b in blocks:
        print(f'  {b}')

# Check output
for fn in ['test_disp.txt', 'test_react.txt']:
    fp = os.path.join(d, fn)
    if os.path.isfile(fp) and os.path.getsize(fp) > 0:
        with open(fp, 'r') as f:
            dl = f.readlines()
        print(f'\n{fn}: {os.path.getsize(fp)} bytes, {len(dl)} lines')
        if dl:
            print(f'  First: {dl[0].rstrip()[:80]}')
            print(f'  Last:  {dl[-1].rstrip()[:80]}')
    else:
        print(f'\n{fn}: empty or not found')

# Quick reaction check
fp = os.path.join(d, 'test_react.txt')
if os.path.isfile(fp) and os.path.getsize(fp) > 0:
    print('\nReaction force magnitude check:')
    with open(fp, 'r') as f:
        for i, line in enumerate(f):
            parts = line.strip().split()
            if len(parts) >= 10:
                vals = [float(x) for x in parts[:10]]
                total = sum(vals[1:])
                if i < 5 or i % 50 == 0:
                    print(f'  row {i:3d}: shear={total/1000:.1f} kN')
