"""Run uniaxial comparison + wall analysis with OpenSees-2145.exe."""
import subprocess, os, sys, time, glob
import numpy as np
sys.stdout.reconfigure(line_buffering=True)

WORKDIR = r'E:\Basic\Concrete_Model\ADINA'
ADINA_EXE = os.path.join(WORKDIR, 'OpenSees-2145.exe')

# Find reference exe
REF_EXE = None
for f in os.listdir(WORKDIR):
    if f.endswith('.exe') and 'OpenSees' in f and any(ord(c) > 127 for c in f):
        REF_EXE = os.path.join(WORKDIR, f)
        break

def run(exe, tcl, label, timeout=120):
    print(f'\n  Running {label}...', flush=True)
    t0 = time.time()
    try:
        r = subprocess.run([exe, tcl], cwd=WORKDIR, capture_output=True, text=True, timeout=timeout)
        dt = time.time() - t0
        code = r.returncode
        # Print key stdout lines
        for line in (r.stdout or '').strip().split('\n'):
            if line.strip() and ('Phase' in line or 'done' in line or 'Failed' in line or 'complete' in line or 'Step' in line or 'fails' in line or 'FAILED' in line or 'gravity' in line):
                print(f'    {line.strip()}')
        print(f'  Exit: {code} ({dt:.1f}s)')
        return code
    except subprocess.TimeoutExpired:
        print(f'  TIMEOUT ({time.time()-t0:.0f}s)')
        return -999

def safe_load(filepath, ncols):
    rows = []
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == ncols:
                try: rows.append([float(x) for x in parts])
                except: pass
    return np.array(rows) if rows else np.zeros((0, ncols))

# ================================================================
# PART 1: UNIAXIAL COMPARISON
# ================================================================
print('='*60)
print('  PART 1: UNIAXIAL SINGLE-ELEMENT TEST')
print('='*60)

# Run reference
if REF_EXE:
    run(REF_EXE, os.path.join(WORKDIR, 'test_uniax_ref.tcl'), 'Reference uniaxial')
else:
    print('  Reference exe not found!')

# Run ADINA
run(ADINA_EXE, os.path.join(WORKDIR, 'test_uniax_adina.tcl'), 'ADINA-2145 uniaxial')

# Load and analyze
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

def load_uniax(prefix):
    rf = os.path.join(WORKDIR, f'uniax_react_{prefix}.txt')
    df = os.path.join(WORKDIR, f'uniax_disp_{prefix}.txt')
    if not os.path.exists(rf) or not os.path.exists(df):
        return None, None
    react = safe_load(rf, 3)
    disp = safe_load(df, 2)
    if len(react) == 0 or len(disp) == 0:
        return None, None
    n = min(len(react), len(disp))
    stress = -(react[:n, 1] + react[:n, 2]) / 0.3 / 1e6  # MPa
    strain = disp[:n, 1]
    return strain, stress

strain_ref, stress_ref = load_uniax('ref')
strain_adi, stress_adi = load_uniax('adina')

print('\n  Results:')
if strain_ref is not None:
    print(f'  Reference: {len(strain_ref)} pts, strain [{strain_ref.min():.6f},{strain_ref.max():.6f}], '
          f'stress [{stress_ref.min():.2f},{stress_ref.max():.2f}] MPa')
    # Peak
    idx = np.argmin(stress_ref)
    print(f'    Peak: {stress_ref[idx]:.2f} MPa at strain {strain_ref[idx]:.6f}')
    # Initial E (first 30 steps)
    n30 = min(30, len(strain_ref))
    if n30 > 2:
        E = np.polyfit(strain_ref[:n30], stress_ref[:n30]*1e6, 1)[0]
        print(f'    Initial E: {E/1e9:.2f} GPa')

if strain_adi is not None:
    print(f'  ADINA:     {len(strain_adi)} pts, strain [{strain_adi.min():.6f},{strain_adi.max():.6f}], '
          f'stress [{stress_adi.min():.2f},{stress_adi.max():.2f}] MPa')
    idx = np.argmin(stress_adi)
    print(f'    Peak: {stress_adi[idx]:.2f} MPa at strain {strain_adi[idx]:.6f}')
    n30 = min(30, len(strain_adi))
    if n30 > 2:
        E = np.polyfit(strain_adi[:n30], stress_adi[:n30]*1e6, 1)[0]
        print(f'    Initial E: {E/1e9:.2f} GPa')

# Saenz envelope
def saenz_curve(eps_arr, E0=21.4e9, fc=-20.7e6, epsc=-0.002, fu=-4.14e6, epsu=-0.006):
    sigma = np.zeros_like(eps_arr)
    Rp = epsu/epsc
    Es = fc/epsc
    Eu = fu/epsu
    A = (E0/Eu + (Rp-2)*Rp**2*E0/Es - (2*Rp+1)*(Rp-1)**2) / (Rp*(Rp-1)**2)
    B = 2*E0/Es - 3 - 2*A
    C = 2 - E0/Es + A
    for i, e in enumerate(eps_arr):
        if e >= 0: sigma[i] = E0*e
        else:
            x = e/epsc
            d = 1+A*x+B*x**2+C*x**3
            sigma[i] = E0*e/d if abs(d)>1e-30 else 0
    return sigma/1e6

eps_env = np.linspace(0, -0.007, 500)
sig_env = saenz_curve(eps_env)

# Plot uniaxial comparison
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

ax = axes[0]
ax.plot(eps_env*1000, sig_env, 'k--', lw=1, alpha=0.5, label='Saenz')
if strain_ref is not None:
    ax.plot(strain_ref*1000, stress_ref, 'b-', lw=1.2, label='Reference')
if strain_adi is not None:
    ax.plot(strain_adi*1000, stress_adi, 'r-', lw=1.2, label='ADINA-2145')
ax.set_xlabel('Strain (x1000)'); ax.set_ylabel('Stress (MPa)')
ax.set_title('Full Stress-Strain'); ax.legend(); ax.grid(True, alpha=0.3)

ax = axes[1]
ax.plot(eps_env*1000, sig_env, 'k--', lw=1, alpha=0.5, label='Saenz')
if strain_ref is not None:
    mask = strain_ref < 0
    ax.plot(strain_ref[mask]*1000, stress_ref[mask], 'b-', lw=1.2, label='Reference')
if strain_adi is not None:
    mask = strain_adi < 0
    ax.plot(strain_adi[mask]*1000, stress_adi[mask], 'r-', lw=1.2, label='ADINA-2145')
ax.set_xlabel('Strain (x1000)'); ax.set_ylabel('Stress (MPa)')
ax.set_xlim([-7, 0.5]); ax.set_title('Compression Detail'); ax.legend(); ax.grid(True, alpha=0.3)

ax = axes[2]
if strain_ref is not None:
    n = len(strain_ref)
    s1, s2 = min(1000, n), min(n, 1600)
    if s2 > s1: ax.plot(strain_ref[s1:s2]*1000, stress_ref[s1:s2], 'b-o', lw=1, ms=1.5, label='Reference')
if strain_adi is not None:
    n = len(strain_adi)
    s1, s2 = min(1000, n), min(n, 1600)
    if s2 > s1: ax.plot(strain_adi[s1:s2]*1000, stress_adi[s1:s2], 'r-s', lw=1, ms=1.5, label='ADINA-2022')
ax.set_xlabel('Strain (x1000)'); ax.set_ylabel('Stress (MPa)')
ax.set_title('Tension + Reload Detail'); ax.legend(); ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(WORKDIR, 'uniax_2145_comparison.png'), dpi=150)
print('\n  Saved: uniax_2145_comparison.png')

# ================================================================
# PART 2: WALL ANALYSIS
# ================================================================
print(f'\n{"="*60}')
print('  PART 2: WALL ANALYSIS (csp3.tcl)')
print('='*60)

# Clean old output
for f in ['1.txt', '53.txt']:
    p = os.path.join(WORKDIR, f)
    if os.path.exists(p): os.remove(p)

run(ADINA_EXE, os.path.join(WORKDIR, 'csp3.tcl'), 'Wall analysis (2145)', timeout=600)

# Parse wall results
react_data, disp_data = [], []
if os.path.exists(os.path.join(WORKDIR, '1.txt')):
    with open(os.path.join(WORKDIR, '1.txt')) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 6:
                try:
                    t = float(parts[0])
                    shear = sum(float(x) for x in parts[1:6]) / 1000
                    react_data.append([t, shear])
                except: pass

if os.path.exists(os.path.join(WORKDIR, '53.txt')):
    with open(os.path.join(WORKDIR, '53.txt')) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                try: disp_data.append([float(parts[0]), float(parts[1])*1000])
                except: pass

react = np.array(react_data) if react_data else np.zeros((0,2))
disp = np.array(disp_data) if disp_data else np.zeros((0,2))

print(f'\n  Displacement: {len(disp)} points')
print(f'  Reaction:     {len(react)} points')

if len(disp) > 0:
    print(f'  Disp range: [{disp[:,1].min():.2f}, {disp[:,1].max():.2f}] mm')
if len(react) > 0:
    good = np.abs(react[:,1]) < 500
    n_good = np.sum(good)
    print(f'  Good force points (|F|<500kN): {n_good} / {len(react)}')
    if n_good > 0:
        rg = react[good]
        print(f'  Force range (good): [{rg[:,1].min():.1f}, {rg[:,1].max():.1f}] kN')

# Build force-displacement
fd_pairs = []
if len(disp) > 0 and len(react) > 0:
    di, ri = 0, 0
    while di < len(disp) and ri < len(react):
        if abs(disp[di,0] - react[ri,0]) < 0.001:
            if abs(react[ri,1]) < 500:
                fd_pairs.append([disp[di,1], react[ri,1]])
            di += 1; ri += 1
        elif disp[di,0] < react[ri,0]: di += 1
        else: ri += 1

# Load reference
ref_rows = []
try:
    with open(os.path.join(WORKDIR, '参考模型.txt')) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    rd, rf = float(parts[0]), float(parts[1])
                    ref_rows.append([rd*1000, rf/1000])
                except: pass
except: pass
ref = np.array(ref_rows) if ref_rows else None

# Plot wall comparison
fig2, axes2 = plt.subplots(1, 2, figsize=(14, 6))
ax = axes2[0]
if ref is not None:
    ax.plot(ref[:,0], ref[:,1], 'b-', lw=0.8, alpha=0.5, label='Reference')
if fd_pairs:
    fd = np.array(fd_pairs)
    ax.plot(fd[:,0], fd[:,1], 'r-', lw=0.8, label='ADINA-2145')
    print(f'\n  Force-disp pairs: {len(fd)}')
    print(f'  Disp range: [{fd[:,0].min():.2f}, {fd[:,0].max():.2f}] mm')
    peaks_disp = []
    for i in range(1, len(fd)-1):
        if (fd[i,0]>fd[i-1,0] and fd[i,0]>fd[i+1,0]) or (fd[i,0]<fd[i-1,0] and fd[i,0]<fd[i+1,0]):
            peaks_disp.append(i)
    if peaks_disp:
        print('  Peaks:')
        for p in peaks_disp[:12]:
            print(f'    d={fd[p,0]:+8.2f}mm  F={fd[p,1]:+8.1f}kN')
ax.set_xlabel('Displacement (mm)'); ax.set_ylabel('Base Shear (kN)')
ax.set_title('Hysteresis: 2145 vs Reference'); ax.legend(); ax.grid(True, alpha=0.3)

# Version comparison
ax = axes2[1]
versions = {}
for ver in ['1045', '1205', '1535', '1706', '2022', '2145']:
    if ver == '2145':
        cnt = n_good if len(react) > 0 else 0
    else:
        rfile = os.path.join(WORKDIR, f'wall_OpenSees_{ver}_react.txt')
        cnt = 0
        try:
            with open(rfile) as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 6:
                        try:
                            fs = sum(float(x) for x in parts[1:6]) / 1000
                            if abs(fs) < 500: cnt += 1
                        except: pass
        except: pass
    versions[ver] = cnt

ref_cnt = len(ref) if ref is not None else 0
print(f'\n  VERSION COMPARISON:')
for v, c in versions.items():
    bar = '#' * (c // 20)
    print(f'    OpenSees-{v}: {c:5d} pts  {bar}')
print(f'    Reference:      {ref_cnt:5d} pts')

# Bar chart
labels = list(versions.keys()) + ['Ref']
counts = list(versions.values()) + [ref_cnt]
colors = ['gray','gray','gray','gray','gray','green','blue']
ax.barh(labels, counts, color=colors)
ax.set_xlabel('Good Force Points (|F|<500kN)')
ax.set_title('Version Comparison')
for i, (l, c) in enumerate(zip(labels, counts)):
    ax.text(c+5, i, str(c), va='center', fontsize=9)

plt.tight_layout()
plt.savefig(os.path.join(WORKDIR, 'wall_2145_comparison.png'), dpi=150)
print('\n  Saved: wall_2145_comparison.png')

# Save reaction data
import shutil
for fn, dst in [('1.txt', 'wall_OpenSees_2145_react.txt'), ('53.txt', 'wall_OpenSees_2145_disp.txt')]:
    src = os.path.join(WORKDIR, fn)
    if os.path.exists(src):
        shutil.copy2(src, os.path.join(WORKDIR, dst))

print('\nDone.')
