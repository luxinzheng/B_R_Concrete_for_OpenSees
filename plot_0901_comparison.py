# -*- coding: utf-8 -*-
"""Detailed comparison: 0901 vs Reference model for sw1-1"""
import numpy as np
import matplotlib.pyplot as plt
import os

plt.rcParams['font.family'] = 'Microsoft YaHei'
plt.rcParams['axes.unicode_minus'] = False

root = r'e:\Basic\Concrete_Model\ADINA'

def load_data(folder, react_file, disp_file):
    rpath = os.path.join(root, folder, react_file)
    dpath = os.path.join(root, folder, disp_file)
    
    times_r, shears = [], []
    if os.path.exists(rpath):
        with open(rpath) as f:
            ncols = None
            for line in f:
                parts = line.strip().split()
                if len(parts) < 2:
                    continue
                try:
                    vals = [float(x) for x in parts]
                    if ncols is None:
                        ncols = len(vals)
                    if len(vals) == ncols:
                        times_r.append(vals[0])
                        shears.append(sum(vals[1:]))
                except:
                    pass
    
    times_d, disps = [], []
    if os.path.exists(dpath):
        with open(dpath) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        times_d.append(float(parts[0]))
                        disps.append(float(parts[1]))
                    except:
                        pass
    
    return np.array(times_r), np.array(shears), np.array(times_d), np.array(disps)

# Load 0901 data
t_r, shear_0901, t_d, disp_0901 = load_data('sw1-1', '1.txt', '53.txt')
# Load reference data
t_r_ref, shear_ref, t_d_ref, disp_ref = load_data('sw1-1', 'ref_1.txt', 'ref_53.txt')

print(f'0901: react={len(t_r)}, disp={len(t_d)}')
print(f'Ref:  react={len(t_r_ref)}, disp={len(t_d_ref)}')

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('sw1-1: OpenSees-02220901 (B&R Secant) vs Reference (HEOM2D)', fontsize=13, fontweight='bold')

# 1. Hysteresis: Disp vs Base Shear
ax = axes[0, 0]
if len(t_d) > 0 and len(t_r) > 0:
    t_common = np.sort(np.unique(np.concatenate([t_r, t_d])))
    d_i = np.interp(t_common, t_d, disp_0901)
    f_i = np.interp(t_common, t_r, shear_0901)
    ax.plot(d_i, f_i / 1000, 'b-', linewidth=0.6, alpha=0.8, label='B&R 0901')

if len(t_d_ref) > 0 and len(t_r_ref) > 0:
    t_common_ref = np.sort(np.unique(np.concatenate([t_r_ref, t_d_ref])))
    d_i_ref = np.interp(t_common_ref, t_d_ref, disp_ref)
    f_i_ref = np.interp(t_common_ref, t_r_ref, shear_ref)
    ax.plot(d_i_ref, f_i_ref / 1000, 'r-', linewidth=0.6, alpha=0.8, label='HEOM2D Ref')

ax.set_xlabel('Top Displacement (mm)')
ax.set_ylabel('Base Shear (kN)')
ax.set_title('Hysteresis Loops')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# 2. Displacement vs Time
ax = axes[0, 1]
if len(t_d) > 0:
    ax.plot(t_d, disp_0901, 'b-', linewidth=0.8, label='B&R 0901')
if len(t_d_ref) > 0:
    ax.plot(t_d_ref, disp_ref, 'r--', linewidth=0.8, label='HEOM2D Ref')
ax.set_xlabel('Time')
ax.set_ylabel('Displacement (mm)')
ax.set_title('Displacement History')
ax.legend()
ax.grid(True, alpha=0.3)

# 3. Base Shear vs Time
ax = axes[1, 0]
if len(t_r) > 0:
    ax.plot(t_r, shear_0901/1000, 'b-', linewidth=0.6, alpha=0.8, label='B&R 0901')
if len(t_r_ref) > 0:
    ax.plot(t_r_ref, shear_ref/1000, 'r-', linewidth=0.6, alpha=0.8, label='HEOM2D Ref')
ax.set_xlabel('Time')
ax.set_ylabel('Base Shear (kN)')
ax.set_title('Base Shear History')
ax.legend()
ax.grid(True, alpha=0.3)

# 4. Envelope comparison
ax = axes[1, 1]
def extract_envelope(t_d, disp, t_r, shear):
    if len(t_d) == 0 or len(t_r) == 0:
        return [], [], [], []
    t_all = np.sort(np.unique(np.concatenate([t_d, t_r])))
    d = np.interp(t_all, t_d, disp)
    f = np.interp(t_all, t_r, shear)
    
    pos_d, pos_f, neg_d, neg_f = [], [], [], []
    d_max, d_min = 0, 0
    for i in range(len(d)):
        if d[i] > d_max:
            d_max = d[i]
            pos_d.append(d[i])
            pos_f.append(f[i])
        if d[i] < d_min:
            d_min = d[i]
            neg_d.append(d[i])
            neg_f.append(f[i])
    return pos_d, pos_f, neg_d, neg_f

pd, pf, nd, nf = extract_envelope(t_d, disp_0901, t_r, shear_0901)
if pd:
    ax.plot(pd, np.array(pf)/1000, 'b-o', markersize=3, label='B&R + envelope')
if nd:
    ax.plot(nd, np.array(nf)/1000, 'b-s', markersize=3)

pd_r, pf_r, nd_r, nf_r = extract_envelope(t_d_ref, disp_ref, t_r_ref, shear_ref)
if pd_r:
    ax.plot(pd_r, np.array(pf_r)/1000, 'r-o', markersize=3, label='HEOM2D + envelope')
if nd_r:
    ax.plot(nd_r, np.array(nf_r)/1000, 'r-s', markersize=3)

ax.set_xlabel('Displacement (mm)')
ax.set_ylabel('Base Shear (kN)')
ax.set_title('Skeleton Curve')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
outpath = os.path.join(root, 'compare_sw1_0901.png')
plt.savefig(outpath, dpi=150, bbox_inches='tight')
print(f'Saved {outpath}')

# Print comparison stats
if len(t_r) > 0:
    print(f'\nB&R 0901: shear range [{shear_0901.min()/1000:.1f}, {shear_0901.max()/1000:.1f}] kN')
if len(t_r_ref) > 0:
    print(f'HEOM2D Ref: shear range [{shear_ref.min()/1000:.1f}, {shear_ref.max()/1000:.1f}] kN')
if len(t_d) > 0:
    print(f'B&R 0901: disp range [{disp_0901.min():.2f}, {disp_0901.max():.2f}] mm')
if len(t_d_ref) > 0:
    print(f'HEOM2D Ref: disp range [{disp_ref.min():.2f}, {disp_ref.max():.2f}] mm')
