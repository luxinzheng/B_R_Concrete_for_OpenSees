# -*- coding: utf-8 -*-
"""Plot ML results from DC with 1N reference force"""
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

# Load 1N DC results
disp = load_data(f'{ml_dir}/test_disp.txt', expected_cols=2)
react = load_data(f'{ml_dir}/test_react.txt', expected_cols=10)

# Load reference
ref_disp = load_data(f'{ml_dir}/ref_disp1.txt', expected_cols=2)
ref_shear = load_data(f'{ml_dir}/ref_shearforce1.txt', expected_cols=10)

print(f'DC 1N: disp={disp.shape}, react={react.shape}')
print(f'Ref:   disp={ref_disp.shape}, shear={ref_shear.shape}')

d_mm = disp[:, 1] * 1000
lam = disp[:, 0]
base_shear = np.sum(react[:, 1:], axis=1)

# Match ref
ref_t_d = ref_disp[:, 0]
ref_t_s = ref_shear[:, 0]
ref_d = ref_disp[:, 1] * 1000
ref_v_all = np.sum(ref_shear[:, 1:], axis=1)
ref_md, ref_mv = [], []
for i in range(len(ref_t_s)):
    idx = np.argmin(np.abs(ref_t_d - ref_t_s[i]))
    if np.abs(ref_t_d[idx] - ref_t_s[i]) < 1:
        ref_md.append(ref_d[idx])
        ref_mv.append(ref_v_all[i])
ref_md = np.array(ref_md)
ref_mv = np.array(ref_mv)

print(f'DC 1N disp: [{d_mm.min():.1f}, {d_mm.max():.1f}] mm')
print(f'DC 1N lambda: [{lam.min():.1f}, {lam.max():.1f}]')
print(f'DC 1N shear: [{base_shear.min()/1000:.1f}, {base_shear.max()/1000:.1f}] kN')

# Match disp and react by lambda
t_d = disp[:, 0]
t_r = react[:, 0]
md, mv = [], []
for i in range(len(t_r)):
    idx = np.argmin(np.abs(t_d - t_r[i]))
    if np.abs(t_d[idx] - t_r[i]) / (np.abs(t_d[idx]) + 1) < 0.01:
        md.append(d_mm[idx])
        mv.append(base_shear[i])
md = np.array(md)
mv = np.array(mv)
print(f'Matched: {len(md)} points')

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# 1. Displacement history
ax = axes[0, 0]
ax.plot(range(len(d_mm)), d_mm, 'b-', lw=1, label=f'B&R DC 1N ({len(d_mm)} steps)')
ax.plot(range(len(ref_d)), ref_d, 'r--', lw=1, label=f'Ref HEOM2D ({len(ref_d)} steps)')
ax.set_xlabel('Step')
ax.set_ylabel('Displacement (mm)')
ax.set_title('Top Displacement History')
ax.legend()
ax.grid(True, alpha=0.3)

# 2. Lambda history
ax = axes[0, 1]
ax.semilogy(range(len(lam)), np.abs(lam), 'b-', lw=1)
ax.set_xlabel('Step')
ax.set_ylabel('|Lambda|')
ax.set_title('DC Load Factor (Lambda)')
ax.grid(True, alpha=0.3)

# 3. Base shear (log) - reactions include lambda*1N
ax = axes[1, 0]
if len(mv) > 0:
    ax.semilogy(range(len(mv)), np.abs(mv), 'b-', lw=1, label='B&R DC 1N')
ax.semilogy(range(len(ref_mv)), np.abs(ref_mv), 'r--', lw=1, label='Ref HEOM2D')
ax.set_xlabel('Step')
ax.set_ylabel('|Base Shear| (N)')
ax.set_title('Base Shear Magnitude')
ax.legend()
ax.grid(True, alpha=0.3)

# 4. Hysteresis (if disp range is reasonable)
ax = axes[1, 1]
if len(md) > 0:
    ax.plot(md, mv/1000, 'b-', lw=1, alpha=0.5, label='B&R DC 1N')
ax.plot(ref_md, ref_mv/1000, 'r-', lw=1, label='Ref HEOM2D')
ax.set_xlabel('Displacement (mm)')
ax.set_ylabel('Base Shear (kN)')
ax.set_title('Hysteresis Comparison')
ax.legend()
ax.grid(True, alpha=0.3)

plt.suptitle('Multi-layer_Shell: DC 1N ref force (b=0.01)', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(r'e:\Basic\Concrete_Model\ADINA\ml_1N_results.png', dpi=150)
print('Saved: ml_1N_results.png')
plt.close()
