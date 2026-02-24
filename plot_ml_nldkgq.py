# -*- coding: utf-8 -*-
"""Plot Multi-layer_Shell NLDKGQ hysteresis curve."""
import os, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root = r'e:\Basic\Concrete_Model\ADINA'
ml_dir = os.path.join(root, 'Multi-layer_Shell')

def load_file(path, expected_ncols=None):
    rows = []
    nc = expected_ncols
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            try:
                vals = [float(x) for x in parts]
            except ValueError:
                continue
            if nc is None:
                nc = len(vals)
            if len(vals) == nc:
                rows.append(vals)
    return np.array(rows) if rows else None

# Load NLDKGQ results
disp = load_file(os.path.join(ml_dir, 'test_disp.txt'), expected_ncols=2)
react = load_file(os.path.join(ml_dir, 'test_react.txt'), expected_ncols=10)

if disp is not None and react is not None:
    t_d, d_m = disp[:, 0], disp[:, 1]
    t_r = react[:, 0]
    base_shear = np.sum(react[:, 1:], axis=1)

    disp_mm, shear_kn = [], []
    for i in range(len(t_r)):
        idx = np.argmin(np.abs(t_d - t_r[i]))
        if np.abs(t_d[idx] - t_r[i]) < 0.5:
            disp_mm.append(d_m[idx] * 1000.0)
            shear_kn.append(base_shear[i] / 1000.0)
    disp_mm = np.array(disp_mm)
    shear_kn = np.array(shear_kn)
    print(f'NLDKGQ: {len(disp_mm)} pts, disp=[{disp_mm.min():.2f}, {disp_mm.max():.2f}]mm, V=[{shear_kn.min():.1f}, {shear_kn.max():.1f}]kN')
else:
    print('ERROR: could not load data')
    exit(1)

# Plot
fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(disp_mm, shear_kn, 'b-', lw=1.0, label='B&R (NLDKGQ)')
ax.set_xlabel('Top Displacement (mm)', fontsize=13)
ax.set_ylabel('Base Shear (kN)', fontsize=13)
ax.set_title('Multi-layer_Shell - ShellNLDKGQ Element', fontsize=15, fontweight='bold')
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
ax.axhline(y=0, color='k', lw=0.5)
ax.axvline(x=0, color='k', lw=0.5)
plt.tight_layout()
out = os.path.join(root, 'ml_nldkgq_hysteresis.png')
plt.savefig(out, dpi=150)
print(f'Saved: {out}')
plt.close()

# Also plot displacement time history
fig2, ax2 = plt.subplots(figsize=(12, 4))
ax2.plot(t_d, d_m * 1000, 'b-', lw=0.8)
ax2.set_xlabel('Time', fontsize=12)
ax2.set_ylabel('Displacement (mm)', fontsize=12)
ax2.set_title('Displacement History - NLDKGQ', fontsize=13)
ax2.grid(True, alpha=0.3)
plt.tight_layout()
out2 = os.path.join(root, 'ml_nldkgq_disp_history.png')
plt.savefig(out2, dpi=150)
print(f'Saved: {out2}')
plt.close()
