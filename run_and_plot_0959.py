# -*- coding: utf-8 -*-
"""Plot sw1-1 & sw2-1: B&R (02210959) vs Reference"""
import os, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root = r'e:\Basic\Concrete_Model\ADINA'

def load_file(path, expected_ncols=None):
    """Load data skipping first 10 lines (gravity phase) and bad lines."""
    rows = []
    nc = expected_ncols
    with open(path) as f:
        all_lines = f.readlines()
    gravity_end = 10
    for line in all_lines[gravity_end:]:
        parts = line.strip().split()
        if not parts:
            continue
        try:
            vals = [float(x) for x in parts]
        except ValueError:
            continue
        if nc is None:
            nc = len(vals)
        if len(vals) != nc:
            continue
        if nc == 2 and abs(vals[1]) > 0.1:
            continue
        rows.append(vals)
    return np.array(rows) if rows else None

def build_hysteresis(disp_path, react_path, ncols_react=6):
    """Build matched (displacement_mm, base_shear_kN) arrays."""
    d = load_file(disp_path, expected_ncols=2)
    r = load_file(react_path, expected_ncols=ncols_react)
    if d is None or r is None:
        print(f'  WARN: no data from {disp_path} or {react_path}')
        return None, None

    t_d, disp_m = d[:, 0], d[:, 1]
    t_r = r[:, 0]
    base_shear_N = np.sum(r[:, 1:], axis=1)

    disp_mm_matched = []
    shear_kN_matched = []
    for i in range(len(t_r)):
        idx = np.argmin(np.abs(t_d - t_r[i]))
        if np.abs(t_d[idx] - t_r[i]) < 0.5:
            disp_mm_matched.append(disp_m[idx] * 1000.0)
            shear_kN_matched.append(base_shear_N[i] / 1000.0)

    return np.array(disp_mm_matched), np.array(shear_kN_matched)

# ── Load all data ────────────────────────────────────────────────
datasets = {}
for name, ctrl_disp, ctrl_react, ref_disp, ref_react in [
    ('sw1-1', '53.txt', '1.txt', 'ref_53.txt', 'ref_1.txt'),
    ('sw2-1', '28.txt', '1.txt', 'ref_28.txt', 'ref_1.txt'),
]:
    wdir = os.path.join(root, name)
    print(f'\n--- {name} ---')

    print(f'  Loading B&R data...')
    br_d, br_v = build_hysteresis(
        os.path.join(wdir, ctrl_disp),
        os.path.join(wdir, ctrl_react))
    if br_d is not None:
        print(f'  B&R: {len(br_d)} pts, disp=[{br_d.min():.2f}, {br_d.max():.2f}]mm, V=[{br_v.min():.1f}, {br_v.max():.1f}]kN')

    print(f'  Loading Reference data...')
    ref_d, ref_v = build_hysteresis(
        os.path.join(wdir, ref_disp),
        os.path.join(wdir, ref_react))
    if ref_d is not None:
        print(f'  Ref: {len(ref_d)} pts, disp=[{ref_d.min():.2f}, {ref_d.max():.2f}]mm, V=[{ref_v.min():.1f}, {ref_v.max():.1f}]kN')

    datasets[name] = (br_d, br_v, ref_d, ref_v)

# ── Plot ─────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for ax, name in zip(axes, ['sw1-1', 'sw2-1']):
    br_d, br_v, ref_d, ref_v = datasets[name]
    if ref_d is not None:
        ax.plot(ref_d, ref_v, 'r-', lw=1.2, alpha=0.8, label='Reference')
    if br_d is not None:
        ax.plot(br_d, br_v, 'b-', lw=1.2, alpha=0.8, label='B&R')
    ax.set_xlabel('Top Displacement (mm)', fontsize=12)
    ax.set_ylabel('Base Shear (kN)', fontsize=12)
    ax.set_title(name, fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='k', lw=0.5)
    ax.axvline(x=0, color='k', lw=0.5)

plt.tight_layout()
out_path = os.path.join(root, 'compare_sw_0959.png')
plt.savefig(out_path, dpi=150)
print(f'\nSaved: {out_path}')
plt.close()
