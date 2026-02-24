"""
Test OpenSees-1045.exe (with Section 6A cracking fix) against previous versions.
Run all 4 element tests + csp3.tcl wall analysis.
"""
import subprocess, os, numpy as np, time, shutil

WORK_DIR = r'e:\Basic\Concrete_Model\ADINA'

VERSIONS = [
    'OpenSees-0929.exe',   # baseline+NaN (reached t=1.38)
    'OpenSees-0943.exe',   # previous latest
    'OpenSees-1045.exe',   # NEW: Section 6A cracking fix
]

TESTS = [
    ('test_elem_comp',  'comp',  '单轴压缩'),
    ('test_elem_tens',  'tens',  '单轴拉伸'),
    ('test_elem_shear', 'shear', '纯剪切'),
    ('test_elem_cyclic','cyclic','循环拉压'),
]

THICKNESS = 0.1  # m

def run_tcl(exe, tcl_name, prefix):
    """Run a TCL test, return (exit_code, stdout, n_points, disp_data, react_data)."""
    exe_path = os.path.join(WORK_DIR, exe)
    tcl_path = os.path.join(WORK_DIR, tcl_name + '.tcl')

    for s in ['_disp.txt', '_react.txt']:
        f = os.path.join(WORK_DIR, prefix + s)
        if os.path.exists(f):
            os.remove(f)

    t0 = time.time()
    try:
        r = subprocess.run([exe_path, tcl_path], capture_output=True, text=True,
                           timeout=120, cwd=WORK_DIR)
        elapsed = time.time() - t0
        exit_code = r.returncode
        stdout = r.stdout
    except subprocess.TimeoutExpired:
        return -999, "TIMEOUT", 0, None, None
    except Exception as e:
        return -998, str(e), 0, None, None

    disp = react = None
    df = os.path.join(WORK_DIR, prefix + '_disp.txt')
    rf = os.path.join(WORK_DIR, prefix + '_react.txt')
    try:
        if os.path.exists(df) and os.path.getsize(df) > 0:
            disp = np.loadtxt(df)
        if os.path.exists(rf) and os.path.getsize(rf) > 0:
            react = np.loadtxt(rf)
    except:
        pass

    n = len(disp) if disp is not None else 0
    return exit_code, stdout, n, disp, react

def safe_load(filepath, ncols):
    rows = []
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == ncols:
                try:
                    vals = [float(p) for p in parts]
                    rows.append(vals)
                except ValueError:
                    pass
    return np.array(rows, dtype=float) if rows else None

def stress_range(disp, react, test_prefix):
    if disp is None or react is None:
        return "NO DATA"
    n = min(len(disp), len(react))
    if test_prefix in ['comp', 'tens', 'cyclic']:
        uy = disp[:n, 2]
        Fy = react[:n, 2] + react[:n, 4]
        stress = -Fy / THICKNESS / 1e6
    else:  # shear
        ux = disp[:n, 1]
        Fx = react[:n, 1] + react[:n, 3]
        stress = -Fx / THICKNESS / 1e6
    valid = np.isfinite(stress)
    s = stress[valid]
    return f"[{s.min():.2f}, {s.max():.2f}] MPa"

# ============================================================
# Part 1: Single-element tests
# ============================================================
print("=" * 70)
print("Part 1: 单元素测试")
print("=" * 70)

for tcl_name, prefix, label in TESTS:
    print(f"\n--- {label} ({tcl_name}) ---")
    for exe in VERSIONS:
        ver = exe.replace('.exe', '')
        ec, stdout, n, disp, react = run_tcl(exe, tcl_name, prefix)
        sr = stress_range(disp, react, prefix)

        # Count fails from stdout
        n_fails = 0
        for line in stdout.split('\n'):
            if 'fails=' in line:
                import re
                m = re.search(r'fails=(\d+)', line)
                if m:
                    n_fails = int(m.group(1))

        status = f"{n}pts, {n_fails}fails" if n > 0 else "NO DATA"
        print(f"  {ver:20s}: {status:16s} stress={sr}")

# ============================================================
# Part 2: csp3.tcl wall analysis
# ============================================================
print("\n" + "=" * 70)
print("Part 2: csp3.tcl 剪力墙分析")
print("=" * 70)

for exe in VERSIONS:
    ver = exe.replace('.exe', '')
    exe_path = os.path.join(WORK_DIR, exe)
    tcl_path = os.path.join(WORK_DIR, 'csp3.tcl')

    # Clean output files
    for f in ['1.txt', '53.txt']:
        p = os.path.join(WORK_DIR, f)
        if os.path.exists(p):
            os.remove(p)

    print(f"\n  Running {ver}...", flush=True)
    t0 = time.time()
    try:
        r = subprocess.run([exe_path, tcl_path], capture_output=True, text=True,
                           timeout=600, cwd=WORK_DIR)
        elapsed = time.time() - t0
        exit_code = r.returncode
    except subprocess.TimeoutExpired:
        elapsed = 600
        exit_code = -999
        r = type('obj', (object,), {'stdout': 'TIMEOUT', 'stderr': ''})()

    # Parse key info from stdout
    lines = r.stdout.split('\n')
    last_time = "?"
    n_fails = 0
    n_skips = 0
    gravity_ok = False
    for line in lines:
        if 'gravity analyze ok' in line:
            gravity_ok = True
        if 'time=' in line or 'time =' in line:
            import re
            m = re.search(r'time[= ]+([0-9.]+)', line)
            if m:
                last_time = m.group(1)
        if 'fails=' in line:
            m = re.search(r'fails=(\d+)', line)
            if m:
                n_fails = int(m.group(1))
        if 'skips=' in line:
            m = re.search(r'skips=(\d+)', line)
            if m:
                n_skips = int(m.group(1))

    # Check output files
    f1 = os.path.join(WORK_DIR, '1.txt')
    f53 = os.path.join(WORK_DIR, '53.txt')
    n1 = 0
    n53 = 0
    if os.path.exists(f1):
        n1 = sum(1 for _ in open(f1))
    if os.path.exists(f53):
        n53 = sum(1 for _ in open(f53))

    # Save output files with version tag
    tag = ver.replace('-', '_')
    for src, name in [(f1, f'wall_{tag}_react.txt'), (f53, f'wall_{tag}_disp.txt')]:
        if os.path.exists(src) and os.path.getsize(src) > 0:
            shutil.copy2(src, os.path.join(WORK_DIR, name))

    exit_hex = f"0x{exit_code & 0xFFFFFFFF:08X}" if exit_code < 0 else str(exit_code)
    print(f"    Exit: {exit_hex}, Time: {elapsed:.1f}s")
    print(f"    Gravity: {'OK' if gravity_ok else 'FAIL'}")
    print(f"    Last time: {last_time}, fails={n_fails}, skips={n_skips}")
    print(f"    Output: 1.txt={n1} lines, 53.txt={n53} lines")

    # Print last few status lines
    status_lines = [l for l in lines if 'Step' in l or 'FAIL' in l or 'Force' in l or 'stuck' in l or 'finished' in l]
    for sl in status_lines[-5:]:
        print(f"    {sl.strip()}")

print("\nDone!")
