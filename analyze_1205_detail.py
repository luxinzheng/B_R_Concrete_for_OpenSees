"""
Detailed analysis of OpenSees-1205 wall results:
1. Force spikes/jumps identification
2. Behavior at divergence point
3. Asymmetry analysis
"""
import numpy as np
import sys
sys.stdout.reconfigure(line_buffering=True)

WORK_DIR = r'e:\Basic\Concrete_Model\ADINA'

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

# Load 1205 data
d = safe_load(f'{WORK_DIR}/wall_OpenSees_1205_disp.txt', 2)
r = safe_load(f'{WORK_DIR}/wall_OpenSees_1205_react.txt', 6)
n = min(len(d), len(r))
t = d[:n, 0]
disp = d[:n, 1] * 1000
force = r[:n, 1:6].sum(axis=1) / 1000

# Skip gravity
lat = 0
for i in range(1, len(t)):
    if t[i] < t[i-1]:
        lat = i
        break
if lat == 0: lat = 10

t = t[lat:]
disp = disp[lat:]
force = force[lat:]

print(f"Total lateral pts: {len(t)}")
print(f"Time range: [{t[0]:.4f}, {t[-1]:.4f}]")
print(f"Disp range: [{disp.min():.2f}, {disp.max():.2f}] mm")
print(f"Force range: [{force.min():.1f}, {force.max():.1f}] kN")

# 1. Force jump analysis
print("\n" + "="*60)
print("1. FORCE JUMPS (|dF/dt| > 500 kN per step)")
print("="*60)
df = np.diff(force)
dt = np.diff(t)
dFdt = df / np.maximum(dt, 1e-10)
jumps = np.where(np.abs(df) > 50)[0]  # >50 kN jump in single step
print(f"Number of jumps >50kN: {len(jumps)}")
for i in jumps[:20]:
    print(f"  step {i}: t={t[i]:.4f}→{t[i+1]:.4f}, d={disp[i]:.2f}→{disp[i+1]:.2f} mm, "
          f"F={force[i]:.1f}→{force[i+1]:.1f} kN (dF={df[i]:.1f})")

# 2. Analysis near failure point
print("\n" + "="*60)
print("2. LAST 30 STEPS (near divergence)")
print("="*60)
for i in range(max(0, len(t)-30), len(t)):
    dd = disp[i] - disp[i-1] if i > 0 else 0
    dF = force[i] - force[i-1] if i > 0 else 0
    K = dF / dd if abs(dd) > 1e-6 else 0
    print(f"  t={t[i]:.5f}: d={disp[i]:+7.2f} mm, F={force[i]:+8.1f} kN, "
          f"dF={dF:+7.1f}, K_eff={K:+8.1f} kN/mm")

# 3. Asymmetry analysis
print("\n" + "="*60)
print("3. ASYMMETRY ANALYSIS")
print("="*60)
# Find envelope points at each target displacement
targets = [1, 2, 3, 4, 5, 6, 7]
for td in targets:
    pos_idx = np.where((disp > td - 0.2) & (disp < td + 0.2))[0]
    neg_idx = np.where((disp > -td - 0.2) & (disp < -td + 0.2))[0]
    if len(pos_idx) > 0 and len(neg_idx) > 0:
        f_pos = force[pos_idx]
        f_neg = force[neg_idx]
        print(f"  d=±{td}mm: F(+)={f_pos.mean():.1f} kN [{f_pos.min():.1f}, {f_pos.max():.1f}], "
              f"F(-)={f_neg.mean():.1f} kN [{f_neg.min():.1f}, {f_neg.max():.1f}]")

# 4. Cycle-by-cycle peak analysis
print("\n" + "="*60)
print("4. CYCLE PEAKS")
print("="*60)
# Cycle boundaries (approximate pseudo-time)
cycles = [(0, 1.5, '±4mm'), (1.5, 3.0, '±7mm'), (3.0, 4.5, '±11mm')]
for tc0, tc1, label in cycles:
    mask = (t >= tc0) & (t < tc1)
    if mask.sum() < 5:
        continue
    d_c = disp[mask]
    f_c = force[mask]
    i_max = np.argmax(d_c)
    i_min = np.argmin(d_c)
    print(f"  Cycle {label}: d=[{d_c.min():.1f}, {d_c.max():.1f}] mm")
    print(f"    Peak(+): d={d_c[i_max]:.1f} mm, F={f_c[i_max]:.1f} kN")
    print(f"    Peak(-): d={d_c[i_min]:.1f} mm, F={f_c[i_min]:.1f} kN")
    # Stiffness at zero-crossing
    z_cross = np.where(np.diff(np.sign(d_c)))[0]
    for zi in z_cross[:4]:
        if abs(d_c[zi+1] - d_c[zi]) > 0.01:
            K = (f_c[zi+1] - f_c[zi]) / (d_c[zi+1] - d_c[zi])
            print(f"    Zero-crossing: K={K:.1f} kN/mm at t={t[mask][zi]:.3f}")

# 5. Force at zero displacement
print("\n" + "="*60)
print("5. FORCE AT ZERO DISPLACEMENT (residual)")
print("="*60)
z_idx = np.where(np.abs(disp) < 0.05)[0]
for i in z_idx:
    print(f"  t={t[i]:.4f}: d={disp[i]:.3f} mm, F={force[i]:.1f} kN")

# 6. Reference comparison at matched displacements
print("\n" + "="*60)
print("6. REFERENCE COMPARISON")
print("="*60)
ref_rows = []
with open(f'{WORK_DIR}/参考模型.txt', 'r') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2:
            try:
                ref_rows.append([float(parts[0]), float(parts[1])])
            except ValueError:
                pass
ref = np.array(ref_rows)
ref_d = ref[10:, 0] * 1000
ref_f = ref[10:, 1] / 1000

ref_targets = [2, 3, 4, 5, -2, -3, -4, -5, -6]
for td in sorted(ref_targets):
    r_idx = np.where(np.abs(ref_d - td) < 0.3)[0]
    o_idx = np.where(np.abs(disp - td) < 0.3)[0]
    if len(r_idx) > 0 and len(o_idx) > 0:
        # Use envelope (first occurrence going outward)
        rf = ref_f[r_idx[0]]
        of = force[o_idx[0]]
        print(f"  d={td:+.0f}mm: Ref={rf:.1f} kN, 1205={of:.1f} kN, ratio={of/rf:.2f}")
