# -*- coding: utf-8 -*-
"""Run sw1-1 and sw2-1 tests with the new exe."""
import subprocess, os, time

root = r'e:\Basic\Concrete_Model\ADINA'
exe = os.path.join(root, 'OpenSees-02212103.exe')

for case in ['sw1-1', 'sw2-1']:
    d = os.path.join(root, case)
    tcl = 'br_csp3.tcl'
    if case == 'sw2-1':
        tcl = 'br_csp3.tcl'
    
    for fn in ['test_disp.txt', 'test_react.txt']:
        fp = os.path.join(d, fn)
        if os.path.isfile(fp):
            os.remove(fp)
    
    print(f'\n=== {case} with {os.path.basename(exe)} ===', flush=True)
    t0 = time.time()
    log_fp = os.path.join(d, f'{case}_log.txt')
    try:
        with open(log_fp, 'w') as logf:
            p = subprocess.run([exe, tcl], cwd=d,
                               stdout=logf, stderr=subprocess.STDOUT, timeout=600)
        elapsed = time.time() - t0
        print(f'  Finished in {elapsed:.1f}s, exit={p.returncode}', flush=True)
    except subprocess.TimeoutExpired:
        elapsed = time.time() - t0
        print(f'  TIMEOUT after {elapsed:.1f}s', flush=True)
    
    # Analyze
    if os.path.isfile(log_fp):
        with open(log_fp, 'r', errors='replace') as f:
            lines = f.readlines()
        nan_c = sum(1 for l in lines if '-1.#IND' in l)
        warn_c = sum(1 for l in lines if 'failed to converge' in l)
        print(f'  Log: {len(lines)} lines, NaN={nan_c}, warnings={warn_c}')
        for l in lines:
            lo = l.lower()
            if 'gravity' in lo or 'done' in lo or 'complete' in lo:
                print(f'    {l.rstrip()}')
    
    # Output files
    for fn in ['test_disp.txt', 'test_react.txt']:
        fp = os.path.join(d, fn)
        if os.path.isfile(fp) and os.path.getsize(fp) > 0:
            with open(fp, 'r') as f:
                lc = sum(1 for _ in f)
            print(f'  {fn}: {os.path.getsize(fp)} bytes, {lc} lines')
        else:
            print(f'  {fn}: empty or not found')
