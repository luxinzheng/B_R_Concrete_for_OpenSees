"""
Run single-element uniaxial tests with both reference and ADINA models,
then generate comparison plots.
"""
import subprocess, os, sys, time, glob
import numpy as np
sys.stdout.reconfigure(line_buffering=True)

WORKDIR = r'E:\Basic\Concrete_Model\ADINA'

# Find executables by scanning directory
REF_EXE = None
ADINA_EXE = None
for f in os.listdir(WORKDIR):
    if not f.endswith('.exe') or not f.startswith('OpenSees'):
        continue
    fp = os.path.join(WORKDIR, f)
    # Reference model has non-ASCII chars in name
    if any(ord(c) > 127 for c in f):
        REF_EXE = fp
    # Find latest numbered ADINA exe
    name = f.replace('OpenSees-', '').replace('.exe', '')
    if name.isdigit():
        if ADINA_EXE is None or f > os.path.basename(ADINA_EXE):
            ADINA_EXE = fp

if REF_EXE is None:
    REF_EXE = os.path.join(WORKDIR, 'OpenSees.exe')
if ADINA_EXE is None:
    ADINA_EXE = os.path.join(WORKDIR, 'OpenSees-1706.exe')

print(f'Reference exe: {os.path.basename(REF_EXE)}')
print(f'ADINA exe:     {os.path.basename(ADINA_EXE)}')

def run_test(exe, tcl, label):
    """Run OpenSees test and return stdout."""
    print(f'\n{"="*60}')
    print(f'  Running {label}...')
    print(f'{"="*60}')
    t0 = time.time()
    try:
        r = subprocess.run([exe, tcl], cwd=WORKDIR, capture_output=True,
                          text=True, timeout=120)
        dt = time.time() - t0
        print(f'  Exit code: {r.returncode} ({dt:.1f}s)')
        if r.stdout.strip():
            for line in r.stdout.strip().split('\n'):
                print(f'  {line}')
        if r.returncode != 0 and r.stderr.strip():
            for line in r.stderr.strip().split('\n')[:5]:
                print(f'  ERR: {line}')
        return r.returncode
    except subprocess.TimeoutExpired:
        print(f'  TIMEOUT')
        return -1

def load_data(react_file, disp_file):
    """Load reaction and displacement data, compute stress-strain."""
    # Robust parsing: skip malformed lines
    def safe_load(filepath, ncols):
        rows = []
        with open(os.path.join(WORKDIR, filepath)) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) == ncols:
                    try:
                        rows.append([float(x) for x in parts])
                    except ValueError:
                        pass
        return np.array(rows) if rows else np.zeros((0, ncols))
    
    react = safe_load(react_file, 3)  # [time, R1_y, R2_y]
    disp = safe_load(disp_file, 2)    # [time, uy]
    
    # Material stress = -reaction / area  (reaction opposes internal force)
    # Area = width * thickness = 1.0 * 0.3 = 0.3 m^2
    stress_pa = -(react[:, 1] + react[:, 2]) / 0.3
    
    # Strain = uy / height = uy / 1.0
    strain_vals = disp[:, 1]
    
    # Match by time
    n = min(len(react), len(disp))
    return strain_vals[:n], stress_pa[:n] / 1e6  # strain, stress in MPa

# ============================================================
# Run tests
# ============================================================
run_test(REF_EXE, os.path.join(WORKDIR, 'test_uniax_ref.tcl'), 'Reference Model')
run_test(ADINA_EXE, os.path.join(WORKDIR, 'test_uniax_adina.tcl'), 'ADINA Model')

# ============================================================
# Load results
# ============================================================
print(f'\n{"="*60}')
print('  Loading results...')
print(f'{"="*60}')

try:
    strain_ref, stress_ref = load_data('uniax_react_ref.txt', 'uniax_disp_ref.txt')
    print(f'  Reference: {len(strain_ref)} points')
    print(f'    Strain range: [{strain_ref.min():.6f}, {strain_ref.max():.6f}]')
    print(f'    Stress range: [{stress_ref.min():.2f}, {stress_ref.max():.2f}] MPa')
    has_ref = True
except Exception as e:
    print(f'  Reference: FAILED to load ({e})')
    has_ref = False

try:
    strain_adi, stress_adi = load_data('uniax_react_adina.txt', 'uniax_disp_adina.txt')
    print(f'  ADINA: {len(strain_adi)} points')
    print(f'    Strain range: [{strain_adi.min():.6f}, {strain_adi.max():.6f}]')
    print(f'    Stress range: [{stress_adi.min():.2f}, {stress_adi.max():.2f}] MPa')
    has_adi = True
except Exception as e:
    print(f'  ADINA: FAILED to load ({e})')
    has_adi = False

# ============================================================
# Analysis
# ============================================================
print(f'\n{"="*60}')
print('  ANALYSIS')
print(f'{"="*60}')

# Saenz curve for comparison
def saenz_curve(eps, E0, fc, epsc, fu, epsu):
    """Theoretical Saenz curve."""
    eps = np.asarray(eps)
    sigma = np.zeros_like(eps)
    for i, e in enumerate(eps):
        if e >= 0:
            sigma[i] = E0 * e  # linear tension
        else:
            x = e / epsc  # normalized strain (positive for compression)
            Rp = epsu / epsc
            Es = fc / epsc
            if abs(epsu * epsc / epsc) > 1e-30:
                Eu = (fu * fc / fc) / (epsu * epsc / epsc)
            else:
                Eu = E0
            A = (E0/Eu + (Rp-2)*Rp**2*E0/Es - (2*Rp+1)*(Rp-1)**2) / (Rp*(Rp-1)**2)
            B = 2*E0/Es - 3 - 2*A
            C = 2 - E0/Es + A
            denom = 1 + A*x + B*x**2 + C*x**3
            if abs(denom) > 1e-30:
                sigma[i] = E0 * e / denom
            else:
                sigma[i] = 0
    return sigma / 1e6  # MPa

# Material parameters
E0 = 21.4e9
fc = -20.7e6
epsc = -0.002
fu = -4.14e6
epsu = -0.006

# Saenz envelope
eps_env = np.linspace(0, -0.007, 500)
sig_env = saenz_curve(eps_env, E0, fc, epsc, fu, epsu)

if has_ref:
    # Find peak compressive stress
    idx_peak = np.argmin(stress_ref)
    print(f'  Reference peak: {stress_ref[idx_peak]:.2f} MPa at strain {strain_ref[idx_peak]:.6f}')
    
    # Find initial stiffness (slope at small strain)
    mask_small = (strain_ref < -1e-6) & (strain_ref > -5e-4)
    if np.sum(mask_small) > 2:
        E_ref = np.polyfit(strain_ref[mask_small], stress_ref[mask_small]*1e6, 1)[0]
        print(f'  Reference initial E: {E_ref/1e9:.2f} GPa')

if has_adi:
    idx_peak = np.argmin(stress_adi)
    print(f'  ADINA peak:     {stress_adi[idx_peak]:.2f} MPa at strain {strain_adi[idx_peak]:.6f}')
    
    mask_small = (strain_adi < -1e-6) & (strain_adi > -5e-4)
    if np.sum(mask_small) > 2:
        E_adi = np.polyfit(strain_adi[mask_small], stress_adi[mask_small]*1e6, 1)[0]
        print(f'  ADINA initial E:     {E_adi/1e9:.2f} GPa')

# ============================================================
# Plots
# ============================================================
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# --- Plot 1: Full stress-strain comparison ---
ax = axes[0, 0]
ax.plot(eps_env*1000, sig_env, 'k--', lw=1, label='Saenz Envelope', alpha=0.5)
if has_ref:
    ax.plot(strain_ref*1000, stress_ref, 'b-', lw=1.2, label='Reference')
if has_adi:
    ax.plot(strain_adi*1000, stress_adi, 'r-', lw=1.2, label='ADINA')
ax.set_xlabel('Strain (x1000)')
ax.set_ylabel('Stress (MPa)')
ax.set_title('Uniaxial Stress-Strain: Full Protocol')
ax.legend()
ax.grid(True, alpha=0.3)

# --- Plot 2: Zoom on compression peak ---
ax = axes[0, 1]
ax.plot(eps_env*1000, sig_env, 'k--', lw=1, label='Saenz Envelope', alpha=0.5)
if has_ref:
    mask = strain_ref < 0
    ax.plot(strain_ref[mask]*1000, stress_ref[mask], 'b-', lw=1.2, label='Reference')
if has_adi:
    mask = strain_adi < 0
    ax.plot(strain_adi[mask]*1000, stress_adi[mask], 'r-', lw=1.2, label='ADINA')
ax.set_xlabel('Strain (x1000)')
ax.set_ylabel('Stress (MPa)')
ax.set_title('Compression Region (zoom)')
ax.set_xlim([-7, 0.5])
ax.legend()
ax.grid(True, alpha=0.3)

# --- Plot 3: Stress vs step (time history) ---
ax = axes[1, 0]
if has_ref:
    ax.plot(range(len(stress_ref)), stress_ref, 'b-', lw=0.8, label='Reference')
if has_adi:
    ax.plot(range(len(stress_adi)), stress_adi, 'r-', lw=0.8, label='ADINA')
# Mark phases
phase_boundaries = [0, 600, 1200, 1250, 1600]
phase_labels = ['Compress', 'Unload', 'Tension', 'Reload']
for i, (b, lbl) in enumerate(zip(phase_boundaries[:-1], phase_labels)):
    ax.axvline(b, color='gray', ls=':', alpha=0.5)
    ax.text(b+10, ax.get_ylim()[0] if i==0 else 0, lbl, fontsize=8, va='bottom')
ax.set_xlabel('Step')
ax.set_ylabel('Stress (MPa)')
ax.set_title('Stress vs Step')
ax.legend()
ax.grid(True, alpha=0.3)

# --- Plot 4: Strain vs step ---
ax = axes[1, 1]
if has_ref:
    ax.plot(range(len(strain_ref)), strain_ref*1000, 'b-', lw=0.8, label='Reference')
if has_adi:
    ax.plot(range(len(strain_adi)), strain_adi*1000, 'r-', lw=0.8, label='ADINA')
for b in phase_boundaries[:-1]:
    ax.axvline(b, color='gray', ls=':', alpha=0.5)
ax.set_xlabel('Step')
ax.set_ylabel('Strain (x1000)')
ax.set_title('Strain vs Step')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(WORKDIR, 'uniax_comparison.png'), dpi=150)
print(f'\nSaved: uniax_comparison.png')

# --- Additional plot: Unloading detail ---
fig2, axes2 = plt.subplots(1, 2, figsize=(14, 5))

# Zoom on unloading from peak
ax = axes2[0]
if has_ref:
    # Find the unloading segment (steps ~400-700: around peak to zero)
    n = len(strain_ref)
    s1, s2 = min(400, n), min(800, n)
    ax.plot(strain_ref[s1:s2]*1000, stress_ref[s1:s2], 'b-o', lw=1, ms=2, label='Reference')
if has_adi:
    n = len(strain_adi)
    s1, s2 = min(400, n), min(800, n)
    ax.plot(strain_adi[s1:s2]*1000, stress_adi[s1:s2], 'r-s', lw=1, ms=2, label='ADINA')
ax.set_xlabel('Strain (x1000)')
ax.set_ylabel('Stress (MPa)')
ax.set_title('Unloading from Post-Peak (detail)')
ax.legend()
ax.grid(True, alpha=0.3)

# Zoom on tension + reload
ax = axes2[1]
if has_ref:
    n = len(strain_ref)
    s1, s2 = min(1100, n), min(n, 1600)
    if s2 > s1:
        ax.plot(strain_ref[s1:s2]*1000, stress_ref[s1:s2], 'b-o', lw=1, ms=2, label='Reference')
if has_adi:
    n = len(strain_adi)
    s1, s2 = min(1100, n), min(n, 1600)
    if s2 > s1:
        ax.plot(strain_adi[s1:s2]*1000, stress_adi[s1:s2], 'r-s', lw=1, ms=2, label='ADINA')
ax.set_xlabel('Strain (x1000)')
ax.set_ylabel('Stress (MPa)')
ax.set_title('Tension + Reload Compression (detail)')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(WORKDIR, 'uniax_detail.png'), dpi=150)
print(f'Saved: uniax_detail.png')

print('\nDone.')
