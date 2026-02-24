# -*- coding: utf-8 -*-
"""Run ML test with sp+LoadControl using 1648 exe"""
import subprocess, os, time

root = r'e:\Basic\Concrete_Model\ADINA'
ml_dir = os.path.join(root, 'Multi-layer_Shell')
exe = os.path.join(root, 'OpenSees-02221648.exe')

for f in ['test_disp.txt', 'test_react.txt']:
    p = os.path.join(ml_dir, f)
    if os.path.exists(p):
        os.remove(p)

log = os.path.join(root, 'log_1648_ML_sp.txt')
print(f'Running ML sp+LoadControl with 1648...')
t0 = time.time()
with open(log, 'w') as lf:
    proc = subprocess.Popen([exe, 'test_br.tcl'], cwd=ml_dir, stdout=lf, stderr=subprocess.STDOUT)
    try:
        proc.wait(timeout=1800)
    except subprocess.TimeoutExpired:
        proc.kill(); proc.wait()
        print('TIMEOUT')
elapsed = time.time() - t0
print(f'Exit: {proc.returncode}, Time: {elapsed:.1f}s')

for f in ['test_disp.txt', 'test_react.txt']:
    p = os.path.join(ml_dir, f)
    if os.path.exists(p) and os.path.getsize(p) > 0:
        with open(p) as fh:
            lines = fh.readlines()
        print(f'{f}: {len(lines)} lines')
        if lines:
            print(f'  first: {lines[0].strip()[:80]}')
            print(f'  last:  {lines[-1].strip()[:80]}')
    else:
        print(f'{f}: MISSING or empty')

with open(log) as lf:
    content = lf.read()
nan_ct = content.lower().count('nan')
for l in content.strip().split('\n')[-8:]:
    print(f'  {l[:120]}')
