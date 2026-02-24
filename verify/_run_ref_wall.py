import subprocess, os, glob

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
verify = os.path.dirname(os.path.abspath(__file__))

ref_exe = None
for entry in os.scandir(root):
    if entry.is_file() and entry.name.endswith('.exe'):
        if any(ord(c) > 127 for c in entry.name) and 'OpenSees' in entry.name:
            ref_exe = entry.path
            break

if ref_exe is None:
    print("ERROR: ref exe not found")
    exit(1)

wrapper = os.path.join(verify, '_run_ref.tcl')
with open(wrapper, 'w') as f:
    f.write('set model_type "ref"\n')
    f.write(f'source {{{os.path.join(verify, "verify_wall.tcl")}}}\n')

print(f"Running: {ref_exe}")
r = subprocess.run([ref_exe, wrapper], cwd=verify, capture_output=True, text=True, timeout=600)
print(f"Exit: {r.returncode}")
for line in r.stdout.splitlines():
    if any(k in line.upper() for k in ['GRAVITY','WALL','STEP','FAIL','SKIP']):
        print(f"  {line}")
