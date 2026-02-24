# -*- coding: utf-8 -*-
"""Run sw1-1 test with OpenSees-02220808.exe (secant fix) for regression check"""
import subprocess, os, time

root = r'e:\Basic\Concrete_Model\ADINA'
sw_dir = os.path.join(root, 'sw1-1')
exe = os.path.join(root, 'OpenSees-02220808.exe')
tcl = 'br_csp3.tcl'
timeout = 1800

for f in ['1.txt', '53.txt']:
    p = os.path.join(sw_dir, f)
    if os.path.exists(p):
        os.rename(p, p + '.bak')

print(f'Running sw1-1: {exe} {tcl}')
t0 = time.time()
log_path = os.path.join(root, 'sw1_0222_log.txt')
with open(log_path, 'w') as logf:
    proc = subprocess.Popen([exe, tcl], cwd=sw_dir, stdout=logf, stderr=subprocess.STDOUT)
    try:
        proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        proc.kill(); proc.wait()
        print(f'TIMEOUT after {timeout}s')

elapsed = time.time() - t0
print(f'Exit: {proc.returncode}, Elapsed: {elapsed:.1f}s')

for f in ['1.txt', '53.txt']:
    p = os.path.join(sw_dir, f)
    if os.path.exists(p) and os.path.getsize(p) > 0:
        with open(p) as fh: lines = fh.readlines()
        print(f'{f}: {len(lines)} lines')
    else:
        print(f'{f}: MISSING or empty')

with open(log_path) as f: log = f.read()
nan_count = log.lower().count('nan')
print(f'Log: {len(log)} chars, NaN: {nan_count}')
for l in log.strip().split('\n')[-5:]:
    print(f'  {l}')
