"""
Compare reference model exe vs 1045 exe with both TCL files.
Test matrix:
  1. OpenSees-参考模型.exe + csp3-参考模型.tcl  -> reference baseline
  2. OpenSees-1045.exe     + csp3.tcl            -> our current best
  3. OpenSees-参考模型.exe + csp3.tcl            -> ref exe + our params
  4. OpenSees-1045.exe     + csp3-参考模型.tcl   -> our exe + ref params
"""
import subprocess, os, numpy as np, shutil, time

WORK_DIR = r'e:\Basic\Concrete_Model\ADINA'

tests = [
    ('OpenSees-\u53c2\u8003\u6a21\u578b.exe', 'csp3-\u53c2\u8003\u6a21\u578b.tcl',
     'ref_ref', 'Ref.exe + Ref.tcl'),
    ('OpenSees-1045.exe', 'csp3.tcl',
     '1045_ours', '1045.exe + csp3.tcl'),
    ('OpenSees-\u53c2\u8003\u6a21\u578b.exe', 'csp3.tcl',
     'ref_ours', 'Ref.exe + csp3.tcl'),
    ('OpenSees-1045.exe', 'csp3-\u53c2\u8003\u6a21\u578b.tcl',
     '1045_ref', '1045.exe + Ref.tcl'),
]

for exe_name, tcl_name, tag, label in tests:
    print(f"\n{'='*60}")
    print(f"Running: {label}")
    print(f"  exe: {exe_name}")
    print(f"  tcl: {tcl_name}")
    print(f"{'='*60}")

    exe_path = os.path.join(WORK_DIR, exe_name)
    tcl_path = os.path.join(WORK_DIR, tcl_name)

    if not os.path.exists(exe_path):
        print(f"  ERROR: exe not found!")
        continue
    if not os.path.exists(tcl_path):
        print(f"  ERROR: tcl not found!")
        continue

    # Clean output files
    for f in ['1.txt', '53.txt']:
        p = os.path.join(WORK_DIR, f)
        if os.path.exists(p):
            os.remove(p)

    t0 = time.time()
    try:
        r = subprocess.run([exe_path, tcl_path], capture_output=True, text=True,
                          timeout=900, cwd=WORK_DIR)
        elapsed = time.time() - t0
        exit_code = r.returncode
        stdout = r.stdout
        stderr = r.stderr
    except subprocess.TimeoutExpired:
        elapsed = 900
        exit_code = -999
        stdout = ""
        stderr = "TIMEOUT"

    exit_hex = f"0x{exit_code & 0xFFFFFFFF:08X}" if exit_code < 0 or exit_code > 255 else str(exit_code)
    print(f"  Exit: {exit_hex}, Time: {elapsed:.1f}s")

    # Check for gravity
    if 'gravity analyze ok' in stdout:
        print(f"  Gravity: OK")
    else:
        print(f"  Gravity: FAIL")

    # Check output files
    for fname in ['1.txt', '53.txt']:
        fpath = os.path.join(WORK_DIR, fname)
        if os.path.exists(fpath):
            nlines = sum(1 for _ in open(fpath))
            fsize = os.path.getsize(fpath)
            print(f"  {fname}: {nlines} lines, {fsize} bytes")

            # Save with tag
            dest = os.path.join(WORK_DIR, f'{tag}_{fname}')
            shutil.copy2(fpath, dest)
        else:
            print(f"  {fname}: NOT FOUND")

    # Print last few lines of stdout
    lines = stdout.strip().split('\n')
    print(f"  stdout ({len(lines)} lines):")
    for line in lines[-5:]:
        if line.strip():
            print(f"    {line.strip()}")

    # Check stderr for errors
    if stderr and stderr != "TIMEOUT":
        err_lines = stderr.strip().split('\n')
        print(f"  stderr ({len(err_lines)} lines):")
        for line in err_lines[:3]:
            if line.strip():
                print(f"    {line.strip()}")

print("\n\nDone!")
