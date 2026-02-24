# -*- coding: utf-8 -*-
"""Run ML and sw1-1 tests with OpenSees-02220901.exe"""
import subprocess, os, time

root = r'e:\Basic\Concrete_Model\ADINA'
exe = os.path.join(root, 'OpenSees-02220901.exe')

tests = [
    ('Multi-layer_Shell', 'test_br.tcl', ['test_disp.txt', 'test_react.txt']),
    ('sw1-1', 'br_csp3.tcl', ['1.txt', '53.txt']),
]

for folder, tcl, outputs in tests:
    wdir = os.path.join(root, folder)
    print(f'\n{"="*60}')
    print(f'Test: {folder} / {tcl}')
    print(f'{"="*60}')
    
    for f in outputs:
        p = os.path.join(wdir, f)
        if os.path.exists(p):
            bak = p + '.0901bak'
            if os.path.exists(bak):
                os.remove(bak)
            os.rename(p, bak)
    
    log_path = os.path.join(root, f'{folder}_0901_log.txt')
    t0 = time.time()
    with open(log_path, 'w') as logf:
        proc = subprocess.Popen([exe, tcl], cwd=wdir, stdout=logf, stderr=subprocess.STDOUT)
        try:
            proc.wait(timeout=1800)
        except subprocess.TimeoutExpired:
            proc.kill(); proc.wait()
            print(f'  TIMEOUT after 1800s')
    
    elapsed = time.time() - t0
    print(f'  Exit: {proc.returncode}, Time: {elapsed:.1f}s')
    
    for f in outputs:
        p = os.path.join(wdir, f)
        if os.path.exists(p) and os.path.getsize(p) > 0:
            with open(p) as fh:
                lines = fh.readlines()
            print(f'  {f}: {len(lines)} lines, {os.path.getsize(p)} bytes')
            if lines:
                first = lines[0].strip()[:80]
                last = lines[-1].strip()[:80]
                print(f'    first: {first}')
                print(f'    last:  {last}')
        else:
            print(f'  {f}: MISSING or empty')
    
    with open(log_path) as f:
        log = f.read()
    nan_count = log.lower().count('nan')
    print(f'  Log: {len(log)} chars, NaN: {nan_count}')
    loglines = log.strip().split('\n')
    print(f'  Last 5 lines:')
    for l in loglines[-5:]:
        print(f'    {l[:100]}')
