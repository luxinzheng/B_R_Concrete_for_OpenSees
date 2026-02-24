# -*- coding: utf-8 -*-
"""Compare OpenSees-02220901.exe results with reference model"""
import numpy as np
import matplotlib.pyplot as plt
import os

plt.rcParams['font.family'] = 'Microsoft YaHei'
plt.rcParams['axes.unicode_minus'] = False

root = r'e:\Basic\Concrete_Model\ADINA'

def load_sw1_data(folder, react_file='1.txt', disp_file='53.txt'):
    rpath = os.path.join(root, folder, react_file)
    dpath = os.path.join(root, folder, disp_file)
    
    react_data = []
    ncols = None
    if os.path.exists(rpath):
        with open(rpath) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        vals = [float(x) for x in parts]
                        if ncols is None:
                            ncols = len(vals)
                        if len(vals) == ncols:
                            react_data.append(vals)
                    except:
                        pass
    
    disp_data = []
    if os.path.exists(dpath):
        with open(dpath) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        disp_data.append([float(parts[0]), float(parts[1])])
                    except:
                        pass
    
    return np.array(react_data) if react_data else None, np.array(disp_data) if disp_data else None

def load_ml_data(folder, react_file='test_react.txt'):
    rpath = os.path.join(root, folder, react_file)
    react_data = []
    ncols = None
    if os.path.exists(rpath):
        with open(rpath) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        vals = [float(x) for x in parts]
                        if ncols is None:
                            ncols = len(vals)
                        if len(vals) == ncols:
                            react_data.append(vals)
                    except:
                        pass
    return np.array(react_data) if react_data else None

# ===== sw1-1 data =====
sw1_react_0901, sw1_disp_0901 = load_sw1_data('sw1-1')

# Reference model data (from previous runs)
ref_react, ref_disp = None, None
ref_dir = os.path.join(root, 'sw1-1')
for suffix in ['_ref', '-ref', '_ref_backup']:
    rp = os.path.join(ref_dir, f'1{suffix}.txt')
    dp = os.path.join(ref_dir, f'53{suffix}.txt')
    if os.path.exists(rp):
        ref_react, ref_disp = load_sw1_data('sw1-1', f'1{suffix}.txt', f'53{suffix}.txt')
        break

# Also try backup files
if ref_react is None:
    for bak in ['1.txt.0901bak', '1.txt.bak']:
        rp = os.path.join(ref_dir, bak)
        dp = os.path.join(ref_dir, bak.replace('1.txt', '53.txt'))
        if os.path.exists(rp) and os.path.exists(dp):
            try:
                with open(rp) as f:
                    rdata = [[float(x) for x in l.strip().split()] for l in f if l.strip()]
                with open(dp) as f:
                    ddata = [[float(x) for x in l.strip().split()[:2]] for l in f if l.strip()]
                ref_react = np.array(rdata)
                ref_disp = np.array(ddata)
                print(f"Using backup: {bak}")
            except:
                pass
            break

# ===== ML data =====
ml_react_0901 = load_ml_data('Multi-layer_Shell')

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('OpenSees-02220901 (Secant Stress + Consistent Tangent) Results', fontsize=14, fontweight='bold')

# SW1-1 hysteresis
ax = axes[0, 0]
if sw1_disp_0901 is not None and sw1_react_0901 is not None:
    # Merge displacement and reaction by time
    t_disp = sw1_disp_0901[:, 0]
    disp = sw1_disp_0901[:, 1]
    t_react = sw1_react_0901[:, 0]
    base_shear = np.sum(sw1_react_0901[:, 1:], axis=1) if sw1_react_0901.shape[1] > 1 else sw1_react_0901[:, 1]
    
    # Interpolate to common time
    t_common = np.intersect1d(np.round(t_disp, 4), np.round(t_react, 4))
    if len(t_common) > 10:
        d_interp = np.interp(t_common, t_disp, disp)
        f_interp = np.interp(t_common, t_react, base_shear)
        ax.plot(d_interp, f_interp / 1000, 'b-', linewidth=0.8, label='B&R 0901')
    else:
        ax.plot(disp, np.zeros_like(disp), 'b.', markersize=1)
    
ax.set_xlabel('Displacement (mm)')
ax.set_ylabel('Base Shear (kN)')
ax.set_title('sw1-1: Hysteresis Loop')
ax.legend()
ax.grid(True, alpha=0.3)

# SW1-1 envelope
ax = axes[0, 1]
if sw1_disp_0901 is not None and sw1_react_0901 is not None:
    t_react = sw1_react_0901[:, 0]
    base_shear = np.sum(sw1_react_0901[:, 1:], axis=1) / 1000
    ax.plot(t_react, base_shear, 'b-', linewidth=0.8, label='B&R 0901')
ax.set_xlabel('Time')
ax.set_ylabel('Base Shear (kN)')
ax.set_title('sw1-1: Base Shear vs Time')
ax.legend()
ax.grid(True, alpha=0.3)

# ML reaction magnitude vs step
ax = axes[1, 0]
if ml_react_0901 is not None:
    load_factor = ml_react_0901[:, 0]
    react_sum = np.sum(ml_react_0901[:, 1:], axis=1)
    steps = np.arange(len(load_factor))
    
    ax.plot(steps, react_sum / 1000, 'b-', linewidth=0.8, label='Sum of reactions')
    ax.plot(steps, load_factor, 'r--', linewidth=0.8, label='Load factor (kN)')
    ax.set_xlabel('Step')
    ax.set_ylabel('Force (kN)')
    ax.set_title(f'Multi-layer_Shell: Reactions ({len(steps)} steps)')
    ax.legend()
    ax.grid(True, alpha=0.3)

# ML max reaction magnitude
ax = axes[1, 1]
if ml_react_0901 is not None:
    max_abs_react = np.max(np.abs(ml_react_0901[:, 1:]), axis=1)
    steps = np.arange(len(max_abs_react))
    ax.semilogy(steps, max_abs_react, 'b-', linewidth=0.8)
    ax.axhline(y=1e5, color='r', linestyle='--', alpha=0.5, label='100 kN threshold')
    ax.set_xlabel('Step')
    ax.set_ylabel('Max |Reaction| (N)')
    ax.set_title('Multi-layer_Shell: Max Reaction Magnitude')
    ax.legend()
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(root, 'results_0901.png'), dpi=150, bbox_inches='tight')
print('Saved results_0901.png')
plt.close()

# Summary stats
print(f'\n=== sw1-1 ===')
if sw1_react_0901 is not None:
    bs = np.sum(sw1_react_0901[:, 1:], axis=1)
    print(f'  Steps: {len(sw1_react_0901)}')
    print(f'  Base shear range: {bs.min()/1000:.1f} to {bs.max()/1000:.1f} kN')
    print(f'  Max |base shear|: {np.max(np.abs(bs))/1000:.1f} kN')

print(f'\n=== Multi-layer_Shell ===')
if ml_react_0901 is not None:
    max_abs = np.max(np.abs(ml_react_0901[:, 1:]), axis=1)
    print(f'  Steps: {len(ml_react_0901)}')
    print(f'  Load factor range: {ml_react_0901[:,0].min():.1f} to {ml_react_0901[:,0].max():.1f}')
    good = np.sum(max_abs < 1e5)
    print(f'  Steps with max |react| < 100 kN: {good}/{len(max_abs)} ({100*good/len(max_abs):.0f}%)')
    print(f'  Max |react|: {max_abs.max():.2e} N')
