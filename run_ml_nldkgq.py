# -*- coding: utf-8 -*-
"""Run Multi-layer_Shell (NLDKGQ) with OpenSees-02210959.exe and plot results."""
import subprocess, os, time, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root = r'e:\Basic\Concrete_Model\ADINA'
exe = os.path.join(root, 'OpenSees-02210959.exe')
wdir = os.path.join(root, 'Multi-layer_Shell')
tcl = 'test_br.tcl'

# Clean old output
for f in ['test_disp.txt', 'test_react.txt']:
    p = os.path.join(wdir, f)
    if os.path.exists(p):
        os.remove(p)

# Run
logf = os.path.join(root, 'log_ml_nldkgq.txt')
print(f'Running Multi-layer_Shell (NLDKGQ) with 02210959...')
t0 = time.time()
with open(logf, 'w') as lf:
    proc = subprocess.Popen([exe, tcl], cwd=wdir, stdout=lf, stderr=subprocess.STDOUT)
    proc.wait(timeout=3600)
elapsed = time.time() - t0
print(f'Exit code: {proc.returncode}, time: {elapsed:.0f}s')

# Check output
for f in ['test_disp.txt', 'test_react.txt']:
    p = os.path.join(wdir, f)
    if os.path.exists(p):
        lc = sum(1 for _ in open(p))
        print(f'  {f}: {lc} lines')
    else:
        print(f'  {f}: NOT FOUND')

# Show last 30 lines of log
with open(logf) as lf:
    lines = lf.readlines()
print(f'\n--- Log tail ({len(lines)} total lines) ---')
for l in lines[-30:]:
    print(l.rstrip())
