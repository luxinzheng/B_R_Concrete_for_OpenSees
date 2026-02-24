# -*- coding: utf-8 -*-
"""Detailed analysis of ML 0901 DC results vs reference"""
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

d_0901 = disp_0901[:, 1] * 1000
v_0901 = np.sum(react_0901[:, 1:], axis=1)
t_0901_d = disp_0901[:, 0]
t_0901_r = react_0901[:, 0]

ref_t_d = ref_disp[:, 0]
ref_t_s = ref_shear[:, 0]
ref_d_all = ref_disp[:, 1] * 1000
ref_v_all = np.sum(ref_shear[:, 1:], axis=1)

# Match ref disp and shear by time
ref_d_list, ref_v_list = [], []
for i in range(len(ref_t_s)):
    idx = np.argmin(np.abs(ref_t_d - ref_t_s[i]))
    if np.abs(ref_t_d[idx] - ref_t_s[i]) < 1.0:
        ref_d_list.append(ref_d_all[idx])
        ref_v_list.append(ref_v_all[i])
ref_d = np.array(ref_d_list)
ref_v = np.array(ref_v_list)
print(f'Ref matched: {len(ref_d)} points, disp=[{ref_d.min():.1f},{ref_d.max():.1f}]mm, V=[{ref_v.min()/1000:.0f},{ref_v.max()/1000:.0f}]kN')

# Find first contiguous block where |V| < 500 kN
threshold = 500e3
first_bad = None
for i in range(len(v_0901)):
    if abs(v_0901[i]) > threshold:
        first_bad = i
        break

if first_bad is None:
    first_bad = len(v_0901)
print(f'0901: first divergence at react step {first_bad} (lambda={t_0901_r[first_bad] if first_bad < len(t_0901_r) else "N/A"})')
print(f'0901: first {first_bad} steps have V in [{v_0901[:first_bad].min()/1000:.1f}, {v_0901[:first_bad].max()/1000:.1f}] kN')

# Match these to displacement
good_d = []
good_v = []
for i in range(first_bad):
    t = t_0901_r[i]
    idx = np.argmin(np.abs(t_0901_d - t))
    if np.abs(t_0901_d[idx] - t) < 0.1:
        good_d.append(d_0901[idx])
        good_v.append(v_0901[i])
good_d = np.array(good_d)
good_v = np.array(good_v)
print(f'Matched good points: {len(good_d)}')
print(f'  disp range: [{good_d.min():.2f}, {good_d.max():.2f}] mm')
print(f'  shear range: [{good_v.min()/1000:.1f}, {good_v.max()/1000:.1f}] kN')

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# 1. Early-phase hysteresis comparison
ax = axes[0, 0]
ax.plot(good_d, good_v/1000, 'b-o', ms=3, lw=1.2, label=f'B&R 0901 ({len(good_d)} steps)')
# Reference: show only the same displacement range
ref_mask = (ref_d >= good_d.min() - 1) & (ref_d <= good_d.max() + 1)
ax.plot(ref_d, ref_v/1000, 'r--', lw=1, alpha=0.5, label='Ref HEOM2D (full)')
ax.set_xlim(good_d.min() - 2, good_d.max() + 2)
ax.set_ylim(good_v.min()/1000 - 50, good_v.max()/1000 + 50)
ax.set_xlabel('Displacement (mm)')
ax.set_ylabel('Base Shear (kN)')
ax.set_title('Early Phase Hysteresis')
ax.legend()
ax.grid(True, alpha=0.3)

# 2. Full reference hysteresis with B&R overlay
ax = axes[0, 1]
ax.plot(ref_d, ref_v/1000, 'r-', lw=1, label='Ref HEOM2D')
ax.plot(good_d, good_v/1000, 'b-', lw=2, label='B&R 0901')
ax.set_xlabel('Displacement (mm)')
ax.set_ylabel('Base Shear (kN)')
ax.set_title('Full Hysteresis Comparison')
ax.legend()
ax.grid(True, alpha=0.3)

# 3. Lambda factor history
ax = axes[1, 0]
ax.plot(range(len(t_0901_r[:first_bad])), t_0901_r[:first_bad], 'b-o', ms=2, label='B&R 0901 lambda')
ax.set_xlabel('Step')
ax.set_ylabel('Load Factor λ')
ax.set_title('DisplacementControl Lambda (good region)')
ax.grid(True, alpha=0.3)
ax.legend()

# 4. Force magnitude for first 50 steps
ax = axes[1, 1]
n50 = min(50, first_bad, len(v_0901))
ax.bar(range(n50), np.abs(v_0901[:n50])/1000, color='blue', alpha=0.6, label='B&R 0901')
# Ref first 50
n50r = min(50, len(ref_v))
ax.bar(np.arange(n50r) + 0.3, np.abs(ref_v[:n50r])/1000, color='red', alpha=0.4, width=0.3, label='Ref HEOM2D')
ax.set_xlabel('Step')
ax.set_ylabel('|Base Shear| (kN)')
ax.set_title('First 50 Steps: Force Magnitude')
ax.legend()
ax.grid(True, alpha=0.3)

plt.suptitle('Multi-layer_Shell: B&R 0901 (DisplacementControl) Analysis', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(r'e:\Basic\Concrete_Model\ADINA\ml_detail_0901.png', dpi=150)
print('Saved: ml_detail_0901.png')
plt.close()

# Summary of model parameters affecting convergence
print('\n=== Model Parameter Analysis ===')
print('Steel02 parameters in ML:')
print('  Mat 3: fy=403.3MPa E=208.3GPa b=0.000 (ZERO hardening)')
print('  Mat 4: fy=366.7MPa E=207.1GPa b=0.000')
print('  Mat 5: fy=350.3MPa E=206.0GPa b=0.000')
print('  Mat 9: fy=366.7MPa E=207.1GPa b=0.003')
print()
print('Steel02 parameters in sw1-1 (works well):')
print('  Mat 7: fy=379MPa E=202.7GPa b=0.01')
print('  Mat 8: fy=392MPa E=200.6GPa b=0.01')
print()
print('Key differences:')
print('  ML b=0.000 → tangent=0 after yield → singular K')
print('  sw1-1 b=0.01 → tangent=E*0.01 after yield → positive K')
print()
print('Concrete STIFAC in both: 0.0001 (very low post-crack stiffness)')
print('  This makes cracked elements contribute almost nothing to K')
