# -*- coding: utf-8 -*-
"""Plot ML DC results with b=0.01 vs reference model"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ml_dir = r'e:\Basic\Concrete_Model\ADINA\Multi-layer_Shell'

def load_data(path, expected_cols=None):
    data = []
    ncols = expected_cols
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            try:
                vals = [float(x) for x in parts]
            except ValueError:
                continue
            if ncols is None:
                ncols = len(vals)
            if len(vals) == ncols:
                data.append(vals)
    return np.array(data) if data else None

# Load new results (b=0.01, original forces, DC)
disp = load_data(f'{ml_dir}/test_disp.txt', expected_cols=2)
react = load_data(f'{ml_dir}/test_react.txt', expected_cols=10)

# Load reference
ref_disp = load_data(f'{ml_dir}/ref_disp1.txt', expected_cols=2)
ref_shear = load_data(f'{ml_dir}/ref_shearforce1.txt', expected_cols=10)

print(f'B&R b=0.01: disp={disp.shape}, react={react.shape}')
print(f'Reference:   disp={ref_disp.shape}, shear={ref_shear.shape}')

d_mm = disp[:, 1] * 1000
lam = disp[:, 0]
base_shear = np.sum(react[:, 1:], axis=1)

# Match disp-react by lambda
t_d = disp[:, 0]
t_r = react[:, 0]
md, mv, ml = [], [], []
for i in range(len(t_r)):
    idx = np.argmin(np.abs(t_d - t_r[i]))
    rel_err = np.abs(t_d[idx] - t_r[i]) / (np.abs(t_d[idx]) + 1e-10)
    if rel_err < 0.01:
        md.append(d_mm[idx])
        mv.append(base_shear[i])
        ml.append(t_r[i])
md = np.array(md)
mv = np.array(mv)
ml = np.array(ml)

# Match reference
ref_t_d = ref_disp[:, 0]
ref_t_s = ref_shear[:, 0]
ref_d_all = ref_disp[:, 1] * 1000
ref_v_all = np.sum(ref_shear[:, 1:], axis=1)
rmd, rmv = [], []
for i in range(len(ref_t_s)):
    idx = np.argmin(np.abs(ref_t_d - ref_t_s[i]))
    if np.abs(ref_t_d[idx] - ref_t_s[i]) < 1:
        rmd.append(ref_d_all[idx])
        rmv.append(ref_v_all[i])
rmd = np.array(rmd)
rmv = np.array(rmv)

print(f'B&R matched: {len(md)} pts, disp=[{md.min():.1f},{md.max():.1f}]mm, V=[{mv.min()/1000:.0f},{mv.max()/1000:.0f}]kN')
print(f'Ref matched: {len(rmd)} pts, disp=[{rmd.min():.1f},{rmd.max():.1f}]mm, V=[{rmv.min()/1000:.0f},{rmv.max()/1000:.0f}]kN')

# Filter reasonable B&R data
mask_ok = np.abs(mv) < 5e6  # < 5000 kN
n_ok = np.sum(mask_ok)
print(f'B&R reasonable (<5000kN): {n_ok}/{len(mv)} pts')

# Find contiguous OK region
first_bad = len(mv)
for i in range(len(mv)):
    if abs(mv[i]) > 5e6:
        first_bad = i
        break
print(f'B&R contiguous OK: first {first_bad} pts, disp=[{md[:first_bad].min():.1f},{md[:first_bad].max():.1f}]mm')

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# 1. Displacement history
ax = axes[0, 0]
ax.plot(range(len(d_mm)), d_mm, 'b-', lw=0.8, label=f'B&R b=0.01 ({len(d_mm)} steps)')
ax.plot(range(len(ref_d_all)), ref_d_all, 'r--', lw=0.8, label=f'Ref HEOM2D ({len(ref_d_all)} steps)')
ax.set_xlabel('Step')
ax.set_ylabel('Displacement (mm)')
ax.set_title('Top Displacement History')
ax.legend()
ax.grid(True, alpha=0.3)

# 2. Lambda history
ax = axes[0, 1]
ax.plot(range(len(lam)), lam, 'b-', lw=0.8)
ax.axhline(0, color='gray', ls=':')
ax.set_xlabel('Step')
ax.set_ylabel('Lambda')
ax.set_title('DC Load Factor History')
ax.grid(True, alpha=0.3)

# 3. Hysteresis (bounded)
ax = axes[1, 0]
if first_bad > 1:
    ax.plot(md[:first_bad], mv[:first_bad]/1000, 'b-', lw=1.2, label=f'B&R b=0.01 ({first_bad} pts)')
ax.plot(rmd, rmv/1000, 'r-', lw=0.8, alpha=0.7, label=f'Ref HEOM2D ({len(rmd)} pts)')
ax.set_xlabel('Displacement (mm)')
ax.set_ylabel('Base Shear (kN)')
ax.set_title('Hysteresis (bounded B&R region)')
ax.legend()
ax.grid(True, alpha=0.3)

# 4. Full hysteresis with reasonable filter
ax = axes[1, 1]
ok_d = md[mask_ok]
ok_v = mv[mask_ok]
if len(ok_d) > 0:
    ax.plot(ok_d, ok_v/1000, 'b.', ms=2, alpha=0.5, label=f'B&R b=0.01 OK ({len(ok_d)} pts)')
ax.plot(rmd, rmv/1000, 'r-', lw=0.8, alpha=0.7, label='Ref HEOM2D')
ax.set_xlabel('Displacement (mm)')
ax.set_ylabel('Base Shear (kN)')
ax.set_title('Hysteresis (all reasonable B&R pts)')
ax.legend()
ax.grid(True, alpha=0.3)

plt.suptitle('Multi-layer_Shell: B&R (DC, b=0.01) vs Ref HEOM2D', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(r'e:\Basic\Concrete_Model\ADINA\ml_b001_results.png', dpi=150)
print('Saved: ml_b001_results.png')
plt.close()
