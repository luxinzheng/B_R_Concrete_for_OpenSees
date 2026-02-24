# -*- coding: utf-8 -*-
"""Plot Multi-layer_Shell: B&R vs Reference model comparison."""
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

# --- B&R model (02230818) ---
disp_br = load_file(os.path.join(ml_dir, 'test_disp.txt'), expected_ncols=2)
react_br = load_file(os.path.join(ml_dir, 'test_react.txt'), expected_ncols=10)

br_disp_mm, br_shear_kn = None, None
if disp_br is not None and react_br is not None:
    t_d, d_m = disp_br[:, 0], disp_br[:, 1]
    t_r = react_br[:, 0]
    base_shear = -np.sum(react_br[:, 1:], axis=1)

    dmm, vkn = [], []
    for i in range(len(t_r)):
        idx = np.argmin(np.abs(t_d - t_r[i]))
        if np.abs(t_d[idx] - t_r[i]) < 0.5:
            dmm.append(d_m[idx] * 1000.0)
            vkn.append(base_shear[i] / 1000.0)
    br_disp_mm = np.array(dmm)
    br_shear_kn = np.array(vkn)
    print(f'B&R:       {len(br_disp_mm)} pts, disp=[{br_disp_mm.min():.1f}, {br_disp_mm.max():.1f}] mm, '
          f'V=[{br_shear_kn.min():.1f}, {br_shear_kn.max():.1f}] kN')
else:
    print('ERROR: B&R data not found'); exit(1)

# --- Reference model ---
ref_data = load_file(os.path.join(ml_dir, 'ref_disp1.txt'), expected_ncols=2)
if ref_data is not None:
    ref_shear_kn = ref_data[:, 0]
    ref_disp_mm = ref_data[:, 1] * 1000.0
    print(f'Reference: {len(ref_disp_mm)} pts, disp=[{ref_disp_mm.min():.1f}, {ref_disp_mm.max():.1f}] mm, '
          f'V=[{ref_shear_kn.min():.1f}, {ref_shear_kn.max():.1f}] kN')
else:
    print('ERROR: Reference data not found'); exit(1)

# --- Plot comparison ---
fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(ref_disp_mm, ref_shear_kn, color='gray', lw=1.5, alpha=0.8, label='Reference')
ax.plot(br_disp_mm, br_shear_kn, 'b-', lw=1.2, label='B&R')
ax.set_xlabel('Top Displacement (mm)', fontsize=13)
ax.set_ylabel('Base Shear (kN)', fontsize=13)
ax.set_title('Multi-layer Shell - Hysteresis Comparison', fontsize=14, fontweight='bold')
ax.legend(fontsize=12, loc='lower right')
ax.grid(True, alpha=0.3)
ax.axhline(y=0, color='k', lw=0.5)
ax.axvline(x=0, color='k', lw=0.5)
plt.tight_layout()
out = os.path.join(root, 'ml_comparison.png')
plt.savefig(out, dpi=150)
print(f'Saved: {out}')
plt.close()

# --- Envelope comparison ---
def extract_envelope(disp, shear):
    """Extract peak shear at each displacement level."""
    pos_d, pos_v, neg_d, neg_v = [], [], [], []
    d_levels = np.arange(0, max(abs(disp.min()), disp.max()) + 1, 2)
    for dl in d_levels:
        mask_p = (disp >= dl - 1) & (disp <= dl + 1)
        mask_n = (disp >= -dl - 1) & (disp <= -dl + 1)
        if mask_p.any():
            pos_d.append(dl)
            pos_v.append(shear[mask_p].max())
        if mask_n.any() and dl > 0:
            neg_d.append(-dl)
            neg_v.append(shear[mask_n].min())
    return (np.array(pos_d), np.array(pos_v),
            np.array(neg_d), np.array(neg_v))

br_pd, br_pv, br_nd, br_nv = extract_envelope(br_disp_mm, br_shear_kn)
ref_pd, ref_pv, ref_nd, ref_nv = extract_envelope(ref_disp_mm, ref_shear_kn)

fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(ref_pd, ref_pv, 's-', color='gray', ms=5, lw=1.5, label='Reference (+)')
ax.plot(ref_nd, ref_nv, 's--', color='gray', ms=5, lw=1.5, label='Reference (-)')
ax.plot(br_pd, br_pv, 'o-', color='blue', ms=5, lw=1.5, label='B&R (+)')
ax.plot(br_nd, br_nv, 'o--', color='blue', ms=5, lw=1.5, label='B&R (-)')
ax.set_xlabel('Top Displacement (mm)', fontsize=13)
ax.set_ylabel('Peak Base Shear (kN)', fontsize=13)
ax.set_title('Multi-layer Shell - Envelope Comparison', fontsize=14, fontweight='bold')
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
ax.axhline(y=0, color='k', lw=0.5)
ax.axvline(x=0, color='k', lw=0.5)
plt.tight_layout()
out2 = os.path.join(root, 'ml_envelope_comparison.png')
plt.savefig(out2, dpi=150)
print(f'Saved: {out2}')
plt.close()
