"""
Test OpenSees-1535.exe (tangent consistency + crack stress zeroing fix)
vs 1205 (principal strain) and 1045 (volumetric check).
"""
import subprocess, os, sys, numpy as np, time, shutil

sys.stdout.reconfigure(line_buffering=True)

WORK_DIR = r'e:\Basic\Concrete_Model\ADINA'

VERSIONS = [
    'OpenSees-1045.exe',
    'OpenSees-1205.exe',
    'OpenSees-1535.exe',
]

TESTS = [
    ('test_elem_comp',  'comp',  'Compression'),
    ('test_elem_tens',  'tens',  'Tension'),
    ('test_elem_shear', 'shear', 'Shear'),
    ('test_elem_cyclic','cyclic','Cyclic'),
]

THICKNESS = 0.1

def safe_load(filepath, ncols):
    rows = []
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= ncols:
                try:
                    vals = [float(parts[i]) for i in range(ncols)]
                    rows.append(vals)
                except ValueError:
                    pass
    return np.array(rows) if rows else None

def run_tcl(exe, tcl_name, prefix):
    exe_path = os.path.join(WORK_DIR, exe)
    tcl_path = os.path.join(WORK_DIR, tcl_name + '.tcl')
    for s in ['_disp.txt', '_react.txt']:
        f = os.path.join(WORK_DIR, prefix + s)
        if os.path.exists(f):
            os.remove(f)
    try:
        r = subprocess.run([exe_path, tcl_path], capture_output=True, text=True,
                           timeout=120, cwd=WORK_DIR)
        exit_code = r.returncode
    except subprocess.TimeoutExpired:
        return -999, 0, None, None
    except Exception as e:
        return -998, 0, None, None

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
    return exit_code, n, disp, react

def stress_range(disp, react, test_prefix):
    if disp is None or react is None:
        return "NO DATA"
    n = min(len(disp), len(react))
    if test_prefix in ['comp', 'tens', 'cyclic']:
        Fy = react[:n, 2] + react[:n, 4]
        stress = -Fy / THICKNESS / 1e6
    else:
        Fx = react[:n, 1] + react[:n, 3]
        stress = -Fx / THICKNESS / 1e6
    valid = np.isfinite(stress)
    s = stress[valid]
    if len(s) == 0:
        return "NO VALID"
    return f"[{s.min():.2f}, {s.max():.2f}] MPa"

# ============================================================
# Part 1: Single-element tests
# ============================================================
print("=" * 70)
print("Part 1: Single-element tests")
print("=" * 70)

for tcl_name, prefix, label in TESTS:
    print(f"\n--- {label} ({tcl_name}) ---")
    for exe in VERSIONS:
        ver = exe.replace('.exe', '')
        ec, n, disp, react = run_tcl(exe, tcl_name, prefix)
        sr = stress_range(disp, react, prefix)
        status = f"{n}pts" if n > 0 else "NO DATA"
        print(f"  {ver:20s}: {status:10s} stress={sr}")

# ============================================================
# Part 2: csp3.tcl wall analysis (1535 only, longer timeout)
# ============================================================
print("\n" + "=" * 70)
print("Part 2: csp3.tcl wall analysis (OpenSees-1535 only, timeout=1200s)")
print("=" * 70)

exe = 'OpenSees-1535.exe'
ver = exe.replace('.exe', '')
exe_path = os.path.join(WORK_DIR, exe)
tcl_path = os.path.join(WORK_DIR, 'csp3.tcl')

for f in ['1.txt', '53.txt']:
    p = os.path.join(WORK_DIR, f)
    if os.path.exists(p):
        os.remove(p)

print(f"\n  Running {ver}...", flush=True)
t0 = time.time()
try:
    r = subprocess.run([exe_path, tcl_path], capture_output=True, text=True,
                      timeout=1200, cwd=WORK_DIR)
    elapsed = time.time() - t0
    exit_code = r.returncode
    stdout = r.stdout
except subprocess.TimeoutExpired:
    elapsed = 1200
    exit_code = -999
    stdout = ""

exit_str = f"0x{exit_code & 0xFFFFFFFF:08X}" if exit_code < 0 or exit_code > 255 else str(exit_code)
print(f"    Exit: {exit_str}, Time: {elapsed:.1f}s")

f1 = os.path.join(WORK_DIR, '1.txt')
f53 = os.path.join(WORK_DIR, '53.txt')
n1 = n53 = 0
if os.path.exists(f1) and os.path.getsize(f1) > 0:
    n1 = sum(1 for _ in open(f1))
if os.path.exists(f53) and os.path.getsize(f53) > 0:
    n53 = sum(1 for _ in open(f53))
print(f"    Output: 1.txt={n1} lines, 53.txt={n53} lines")

# Save with version tag
tag = ver.replace('-', '_')
for src, name in [(f1, f'wall_{tag}_react.txt'), (f53, f'wall_{tag}_disp.txt')]:
    if os.path.exists(src) and os.path.getsize(src) > 0:
        shutil.copy2(src, os.path.join(WORK_DIR, name))

# Analyze results
if n53 > 0:
    d = safe_load(f53, 2)
    if d is not None and len(d) > 10:
        lat = 0
        for i in range(1, len(d)):
            if d[i, 0] < d[i-1, 0]:
                lat = i
                break
        if lat == 0: lat = 10
        disp_mm = d[lat:, 1] * 1000
        t_lat = d[lat:, 0]
        print(f"    Lateral: {len(t_lat)} pts, t=[{t_lat[0]:.3f}, {t_lat[-1]:.3f}]")
        print(f"    Disp: [{disp_mm.min():.1f}, {disp_mm.max():.1f}] mm")

        # Cycle analysis
        r_data = safe_load(f1, 6)
        if r_data is not None:
            nr = min(len(d), len(r_data))
            force = r_data[lat:nr, 1:6].sum(axis=1) / 1000
            t_f = d[lat:nr, 0]

            # Key displacement points
            targets = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
            print("\n    Key pseudo-time points:")
            for tt in targets:
                if tt <= t_f[-1]:
                    idx = np.argmin(np.abs(t_f - tt))
                    print(f"      t={tt:.1f}: d={disp_mm[idx]:+7.2f} mm, F={force[idx]:+8.1f} kN")

            # Residual force at zero disp
            z_idx = np.where(np.abs(disp_mm[:len(force)]) < 0.05)[0]
            print(f"\n    Residual force at d≈0:")
            for i in z_idx[:5]:
                print(f"      t={t_f[i]:.4f}: d={disp_mm[i]:.3f} mm, F={force[i]:.1f} kN")

            # Force stability
            peak_pos = force.max()
            peak_neg = force.min()
            print(f"\n    Force peaks: [{peak_neg:.0f}, {peak_pos:.0f}] kN")
            print(f"    |Peak| > 500kN count: {np.sum(np.abs(force) > 500)}")

# Print last stdout lines
lines = [l.strip() for l in stdout.split('\n') if l.strip()]
if lines:
    print(f"\n    Last output: {lines[-1][:120]}")

print("\nDone!")
