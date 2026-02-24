# -*- coding: utf-8 -*-
"""Compare Multi-layer_Shell 0901 DC results with reference model"""
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

disp_0901 = load_data(f'{ml_dir}/test_disp.txt.0901bak', expected_cols=2)
react_0901 = load_data(f'{ml_dir}/test_react.txt.0901bak', expected_cols=10)
ref_disp = load_data(f'{ml_dir}/ref_disp1.txt', expected_cols=2)
ref_shear = load_data(f'{ml_dir}/ref_shearforce1.txt', expected_cols=10)

print(f'0901 disp: {disp_0901.shape}')
print(f'0901 react: {react_0901.shape}')
print(f'ref_disp: {ref_disp.shape}')
print(f'ref_shear: {ref_shear.shape}')

# 0901: match disp and react by time (lambda factor)
t_d = disp_0901[:, 0]
t_r = react_0901[:, 0]
d_vals = disp_0901[:, 1] * 1000  # m -> mm
base_shear_0901 = np.sum(react_0901[:, 1:], axis=1)

# Build matched arrays: for each react time, find nearest disp time
matched_d = []
matched_v = []
for i in range(len(t_r)):
    idx = np.argmin(np.abs(t_d - t_r[i]))
    if np.abs(t_d[idx] - t_r[i]) < 0.01:
        matched_d.append(d_vals[idx])
        matched_v.append(base_shear_0901[i])
matched_d = np.array(matched_d)
matched_v = np.array(matched_v)
print(f'Matched 0901: {len(matched_d)} points')

# Reference: match disp and shear
ref_d = ref_disp[:, 1] * 1000  # mm
ref_base_shear = np.sum(ref_shear[:, 1:], axis=1)
ref_t_d = ref_disp[:, 0]
ref_t_s = ref_shear[:, 0]
ref_md = []
ref_mv = []
for i in range(len(ref_t_s)):
    idx = np.argmin(np.abs(ref_t_d - ref_t_s[i]))
    if np.abs(ref_t_d[idx] - ref_t_s[i]) < 0.01:
        ref_md.append(ref_d[idx])
        ref_mv.append(ref_base_shear[i])
ref_md = np.array(ref_md)
ref_mv = np.array(ref_mv)
print(f'Matched ref: {len(ref_md)} points')

# Find reasonable 0901 range
mask_ok = np.abs(matched_v) < 5e5  # < 500kN
n_ok = np.sum(mask_ok)
print(f'0901 reasonable: {n_ok}/{len(matched_v)} points')

# Stats
print(f'Ref shear range: [{ref_mv.min()/1000:.1f}, {ref_mv.max()/1000:.1f}] kN')
print(f'Ref disp range: [{ref_md.min():.1f}, {ref_md.max():.1f}] mm')
if n_ok > 0:
    ok_v = matched_v[mask_ok]
    ok_d = matched_d[mask_ok]
    print(f'0901 shear (OK): [{ok_v.min()/1000:.1f}, {ok_v.max()/1000:.1f}] kN')
    print(f'0901 disp  (OK): [{ok_d.min():.1f}, {ok_d.max():.1f}] mm')

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# 1. Disp history
ax = axes[0, 0]
ax.plot(range(len(d_vals)), d_vals, 'b-', lw=1, label=f'B&R 0901 ({len(d_vals)} steps)')
ax.plot(range(len(ref_d)), ref_d, 'r--', lw=1, label=f'Ref HEOM2D ({len(ref_d)} steps)')
ax.set_xlabel('Step')
ax.set_ylabel('Displacement (mm)')
ax.set_title('Top Displacement History')
ax.legend()
ax.grid(True, alpha=0.3)

# 2. Base shear (log)
ax = axes[0, 1]
ax.semilogy(range(len(matched_v)), np.abs(matched_v), 'b-', lw=1, label='B&R 0901')
ax.semilogy(range(len(ref_mv)), np.abs(ref_mv), 'r--', lw=1, label='Ref HEOM2D')
ax.axhline(5e5, color='gray', ls=':', label='500kN threshold')
ax.set_xlabel('Step')
ax.set_ylabel('|Base Shear| (N)')
ax.set_title('Force Magnitude')
ax.legend()
ax.grid(True, alpha=0.3)

# 3. Hysteresis (bounded)
ax = axes[1, 0]
if n_ok > 0:
    # Find contiguous OK region
    ok_idx = np.where(mask_ok)[0]
    last = ok_idx[-1] + 1
    ax.plot(matched_d[:last], matched_v[:last]/1000, 'b-', lw=1.2, label=f'B&R 0901 ({last} steps)')
ax.plot(ref_md, ref_mv/1000, 'r--', lw=1, label=f'Ref HEOM2D ({len(ref_md)} steps)')
ax.set_xlabel('Displacement (mm)')
ax.set_ylabel('Base Shear (kN)')
ax.set_title('Hysteresis (B&R bounded region)')
ax.legend()
ax.grid(True, alpha=0.3)

# 4. Hysteresis (full ref only, with annotation)
ax = axes[1, 1]
ax.plot(ref_md, ref_mv/1000, 'r-', lw=1, label='Ref HEOM2D (complete)')
if n_ok > 1:
    ax.plot(matched_d[:last], matched_v[:last]/1000, 'b-', lw=1.5, label='B&R 0901 (partial)')
ax.set_xlabel('Displacement (mm)')
ax.set_ylabel('Base Shear (kN)')
ax.set_title('Hysteresis Comparison')
ax.legend()
ax.grid(True, alpha=0.3)
ax.annotate(f'B&R: {last} of {len(matched_d)} steps physical\n'
            f'Forces diverge after step {last}',
            xy=(0.05, 0.05), xycoords='axes fraction',
            fontsize=9, color='blue',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

plt.suptitle('Multi-layer_Shell: B&R (0901 DC) vs Reference HEOM2D', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(r'e:\Basic\Concrete_Model\ADINA\compare_ml_0901.png', dpi=150)
print('Saved: compare_ml_0901.png')
plt.close()
