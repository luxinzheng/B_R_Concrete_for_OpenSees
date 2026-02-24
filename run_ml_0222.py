# -*- coding: utf-8 -*-
"""Run Multi-layer_Shell test with OpenSees-02220808.exe (secant fix)"""
import subprocess, os, time

root = r'e:\Basic\Concrete_Model\ADINA'
ml_dir = os.path.join(root, 'Multi-layer_Shell')
exe = os.path.join(root, 'OpenSees-02220808.exe')
tcl = 'test_br.tcl'
timeout = 1800  # 30 min

# Clean old output
for f in ['test_disp.txt', 'test_react.txt']:
    p = os.path.join(ml_dir, f)
    if os.path.exists(p):
        os.remove(p)

print(f'Running: {exe}')
print(f'TCL: {tcl}')
print(f'Dir: {ml_dir}')
print(f'Timeout: {timeout}s')

t0 = time.time()
log_path = os.path.join(root, 'ml_0222_log.txt')
with open(log_path, 'w') as logf:
    proc = subprocess.Popen(
        [exe, tcl],
        cwd=ml_dir,
        stdout=logf,
        stderr=subprocess.STDOUT
    )
    try:
        proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.wait()
        print(f'TIMEOUT after {timeout}s')

elapsed = time.time() - t0
print(f'Exit code: {proc.returncode}')
print(f'Elapsed: {elapsed:.1f}s')

# Check outputs
for f in ['test_disp.txt', 'test_react.txt']:
    p = os.path.join(ml_dir, f)
    if os.path.exists(p) and os.path.getsize(p) > 0:
        with open(p) as fh:
            lines = fh.readlines()
        print(f'{f}: {len(lines)} lines, {os.path.getsize(p)} bytes')
    else:
        print(f'{f}: MISSING or empty')

# Check log for NaN / errors
with open(log_path) as f:
    log = f.read()
nan_count = log.lower().count('nan')
print(f'Log: {len(log)} chars, NaN mentions: {nan_count}')
# Show last lines
lines = log.strip().split('\n')
print(f'Last 10 lines:')
for l in lines[-10:]:
    print(f'  {l}')
