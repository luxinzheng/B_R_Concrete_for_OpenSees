# -*- coding: utf-8 -*-
"""Run Multi-layer_Shell test with wipeAnalysis-enabled TCL."""
import subprocess, os, time

root = r'e:\Basic\Concrete_Model\ADINA'
exe = os.path.join(root, 'OpenSees-02212103.exe')
d = os.path.join(root, 'Multi-layer_Shell')

for fn in ['test_disp.txt', 'test_react.txt', 'ml_log.txt']:
    fp = os.path.join(d, fn)
    if os.path.isfile(fp):
        os.remove(fp)

print(f'Running Multi-layer_Shell with {os.path.basename(exe)}...', flush=True)
t0 = time.time()
try:
    log_fp = os.path.join(d, 'ml_log.txt')
    with open(log_fp, 'w') as logf:
        p = subprocess.run([exe, 'test_br.tcl'], cwd=d,
                           stdout=logf, stderr=subprocess.STDOUT, timeout=1800)
    elapsed = time.time() - t0
    print(f'Finished in {elapsed:.1f}s, exit code={p.returncode}', flush=True)
except subprocess.TimeoutExpired:
    elapsed = time.time() - t0
    print(f'TIMEOUT after {elapsed:.1f}s', flush=True)

# Analyze log
log_fp = os.path.join(d, 'ml_log.txt')
if os.path.isfile(log_fp):
    with open(log_fp, 'r', errors='replace') as f:
        lines = f.readlines()
    print(f'\nLog: {len(lines)} lines')
    nan_c = sum(1 for l in lines if '-1.#IND' in l)
    warn_c = sum(1 for l in lines if 'failed to converge' in l)
    print(f'NaN (-1.#IND): {nan_c}')
    print(f'Convergence warnings: {warn_c}')
    for l in lines:
        lo = l.lower()
        if 'gravity' in lo or 'completed block' in lo or 'all done' in lo:
            print(f'  {l.rstrip()}')

# Output files
print('\nOutput files:', flush=True)
for fn in ['test_disp.txt', 'test_react.txt']:
    fp = os.path.join(d, fn)
    if os.path.isfile(fp):
        sz = os.path.getsize(fp)
        lc = 0
        if sz > 0:
            with open(fp, 'r') as f:
                lc = sum(1 for _ in f)
        print(f'  {fn}: {sz} bytes, {lc} lines', flush=True)
        if lc > 0:
            with open(fp, 'r') as f:
                dl = f.readlines()
            print(f'    First: {dl[0].rstrip()}')
            print(f'    Last:  {dl[-1].rstrip()}')
    else:
        print(f'  {fn}: NOT FOUND', flush=True)
