# -*- coding: utf-8 -*-
"""Run ML and sw1-1 tests with OpenSees-02221648.exe (STIFAC>=0.001 fix)"""
import subprocess, os, time, shutil

root = r'e:\Basic\Concrete_Model\ADINA'
exe = os.path.join(root, 'OpenSees-02221648.exe')

tests = [
    ('sw1-1',            'br_csp3.tcl',  ['53.txt', '1.txt']),
    ('Multi-layer_Shell', 'test_br.tcl', ['test_disp.txt', 'test_react.txt']),
]

for name, tcl, outputs in tests:
    wdir = os.path.join(root, name)
    print(f'\n{"="*60}')
    print(f'Testing: {name} / {tcl}')
    print(f'{"="*60}')
    
    # Clean old outputs
    for f in outputs:
        p = os.path.join(wdir, f)
        if os.path.exists(p):
            os.remove(p)
    
    logfile = os.path.join(root, f'log_1648_{name}.txt')
    t0 = time.time()
    with open(logfile, 'w') as lf:
        proc = subprocess.Popen([exe, tcl], cwd=wdir, stdout=lf, stderr=subprocess.STDOUT)
        try:
            proc.wait(timeout=1800)
        except subprocess.TimeoutExpired:
            proc.kill(); proc.wait()
            print('  TIMEOUT (30min)')
    elapsed = time.time() - t0
    print(f'  Exit: {proc.returncode}, Time: {elapsed:.1f}s')
    
    # Check outputs
    for f in outputs:
        p = os.path.join(wdir, f)
        if os.path.exists(p) and os.path.getsize(p) > 0:
            with open(p) as fh:
                lines = fh.readlines()
            print(f'  {f}: {len(lines)} lines')
            if lines:
                print(f'    first: {lines[0].strip()[:80]}')
                print(f'    last:  {lines[-1].strip()[:80]}')
        else:
            print(f'  {f}: MISSING or empty')
    
    # Check log for issues
    with open(logfile) as lf:
        log = lf.read()
    nan_ct = log.lower().count('nan')
    fail_ct = log.lower().count('failed')
    print(f'  Log: {len(log)} chars, NaN mentions: {nan_ct}, fail mentions: {fail_ct}')
    
    # Show last few lines
    loglines = log.strip().split('\n')
    for l in loglines[-5:]:
        print(f'    {l[:100]}')

print('\n\nAll tests complete.')
