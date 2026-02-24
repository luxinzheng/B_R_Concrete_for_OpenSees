# -*- coding: utf-8 -*-
"""Run Multi-layer_Shell (NLDKGQ + sp+LoadControl) with 02210959."""
import subprocess, os, time

root = r'e:\Basic\Concrete_Model\ADINA'
exe = os.path.join(root, 'OpenSees-02230818.exe')
wdir = os.path.join(root, 'Multi-layer_Shell')

for f in ['test_disp.txt', 'test_react.txt']:
    p = os.path.join(wdir, f)
    if os.path.exists(p):
        os.remove(p)

logf = os.path.join(root, 'log_ml_sp_nldkgq.txt')
print('Running Multi-layer_Shell (NLDKGQ + sp+LC)...')
t0 = time.time()
with open(logf, 'w') as lf:
    proc = subprocess.Popen([exe, 'test_br.tcl'], cwd=wdir, stdout=lf, stderr=subprocess.STDOUT)
    proc.wait(timeout=7200)
elapsed = time.time() - t0
print(f'Exit code: {proc.returncode}, time: {elapsed:.0f}s')

for f in ['test_disp.txt', 'test_react.txt']:
    p = os.path.join(wdir, f)
    if os.path.exists(p):
        lc = sum(1 for _ in open(p))
        print(f'  {f}: {lc} lines')
    else:
        print(f'  {f}: NOT FOUND')

with open(logf) as lf:
    lines = lf.readlines()
print(f'\n--- Log tail ({len(lines)} total lines) ---')
for l in lines[-20:]:
    print(l.rstrip())
