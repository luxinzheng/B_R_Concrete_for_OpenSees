# -*- coding: utf-8 -*-
import subprocess, os, time

root = r'e:\Basic\Concrete_Model\ADINA'
ref_exe = os.path.join(root, 'OpenSees-\u53c2\u8003\u6a21\u578b.exe')
print(f'ref exe exists: {os.path.isfile(ref_exe)}', flush=True)
d = os.path.join(root, 'Multi-layer_Shell')

for fn in ['test_disp.txt', 'test_react.txt']:
    fp = os.path.join(d, fn)
    if os.path.isfile(fp):
        os.remove(fp)

t0 = time.time()
try:
    p = subprocess.run([ref_exe, 'test_br.tcl'], cwd=d,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=300)
    elapsed = time.time() - t0
    stdout = p.stdout.decode('utf-8', errors='replace')
    print(f'Done in {elapsed:.1f}s, exit={p.returncode}', flush=True)
    for line in stdout.splitlines():
        if any(k in line.lower() for k in ['gravity', 'completed', 'block', 'all done', 'fail', '-1.#ind']):
            print(f'  {line}', flush=True)
    for fn in ['test_disp.txt', 'test_react.txt']:
        fp = os.path.join(d, fn)
        sz = os.path.getsize(fp) if os.path.isfile(fp) else 0
        print(f'  {fn}: {sz} bytes', flush=True)
except subprocess.TimeoutExpired:
    print(f'TIMEOUT after {time.time()-t0:.1f}s', flush=True)
    for fn in ['test_disp.txt', 'test_react.txt']:
        fp = os.path.join(d, fn)
        sz = os.path.getsize(fp) if os.path.isfile(fp) else 0
        print(f'  {fn}: {sz} bytes', flush=True)
except Exception as e:
    print(f'ERROR: {e}')
