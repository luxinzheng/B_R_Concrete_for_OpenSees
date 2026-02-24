# -*- coding: utf-8 -*-
"""Plot comparison: OpenSees-02241821 ADINA model vs Reference for 3 test cases."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

root = r'e:\Basic\Concrete_Model\ADINA'


def load_robust(filepath, expected_cols=None):
    """Load space-delimited data, skipping lines with wrong column counts."""
    rows = []
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            try:
                vals = [float(x) for x in parts]
            except ValueError:
                continue
            if expected_cols is not None and len(vals) != expected_cols:
                continue
            if expected_cols is None:
                if not rows:
                    expected_cols = len(vals)
                elif len(vals) != expected_cols:
                    continue
            rows.append(vals)
    if rows:
        arr = np.array(rows)
        # Trim last row: OpenSees files often have truncated last lines
        if len(arr) > 1:
            arr = arr[:-1]
        return arr
    return np.empty((0, expected_cols or 1))


def skip_gravity(data):
    """Skip gravity phase: detect time going 0.1..1.0 then resetting to ~0.
    Gravity from loadConst produces: t=0.1,0.2,...,1.0,[dup 1.0],0.1,...
    Returns data from the lateral phase start."""
    if len(data) < 12:
        return data
    t = data[:, 0]
    for i in range(5, min(20, len(t))):
        if t[i] < t[i-1] - 0.05 and t[i-1] >= 0.9 and t[i] < 0.5:
            return data[i:]
    return data


# ═══════════════════════════════════════════════════
# Load & process all data
# ═══════════════════════════════════════════════════
cases = [
    {
        'name': 'SW1-1 (Cyclic Shear Wall)',
        'dir': os.path.join(root, 'sw1-1'),
        'disp': '53.txt', 'react': '1.txt',
        'ref_disp': 'ref_53.txt', 'ref_react': 'ref_1.txt',
        'dcols': 2, 'rcols': 6, 'rr': (1, 6),
    },
    {
        'name': 'SW2-1 (Cyclic Shear Wall)',
        'dir': os.path.join(root, 'sw2-1'),
        'disp': 'test_disp.txt', 'react': 'test_react.txt',
        'ref_disp': 'ref_28.txt', 'ref_react': 'ref_1.txt',
        'dcols': 2, 'rcols': 6, 'rr': (1, 6),
    },
    {
        'name': 'Multi-layer Shell (Cyclic)',
        'dir': os.path.join(root, 'Multi-layer_Shell'),
        'disp': 'test_disp.txt', 'react': 'test_react.txt',
        'ref_disp': 'ref_disp1.txt', 'ref_react': 'ref_shearforce1.txt',
        'dcols': 2, 'rcols': 10, 'rr': (1, 10),
    },
]

all_data = []
for c in cases:
    d_raw = load_robust(os.path.join(c['dir'], c['disp']), c['dcols'])
    r_raw = load_robust(os.path.join(c['dir'], c['react']), c['rcols'])
    rd_raw = load_robust(os.path.join(c['dir'], c['ref_disp']), c['dcols'])
    rr_raw = load_robust(os.path.join(c['dir'], c['ref_react']), c['rcols'])

    d = skip_gravity(d_raw)
    r = skip_gravity(r_raw)
    rd = skip_gravity(rd_raw)
    rr = skip_gravity(rr_raw)

    n = min(len(d), len(r))
    nr = min(len(rd), len(rr))

    info = {
        'name': c['name'], 'rr': c['rr'],
        'd': d[:n], 'r': r[:n],
        'rd': rd[:nr], 'rr_data': rr[:nr],
    }
    print(f"{c['name']}: ADINA={n} pts (raw d/r={len(d_raw)}/{len(r_raw)}), "
          f"Ref={nr} pts (raw d/r={len(rd_raw)}/{len(rr_raw)})")
    if n > 0:
        print(f"  ADINA disp range: [{d[:n,1].min()*1000:.2f}, {d[:n,1].max()*1000:.2f}] mm")
    if nr > 0:
        print(f"  Ref   disp range: [{rd[:nr,1].min()*1000:.2f}, {rd[:nr,1].max()*1000:.2f}] mm")
    all_data.append(info)

# ═══════════════════════════════════════════════════
# Plot
# ═══════════════════════════════════════════════════
fig, axes = plt.subplots(3, 2, figsize=(16, 15))
fig.suptitle('OpenSees-02241821  —  ADINA Concrete Model vs Reference',
             fontsize=14, fontweight='bold', y=0.995)

for row, info in enumerate(all_data):
    ax_fd = axes[row, 0]
    ax_hist = axes[row, 1]
    rng = info['rr']

    if len(info['d']) > 0:
        u = info['d'][:, 1] * 1000
        t_ad = info['d'][:, 0]
        F = np.sum(info['r'][:, rng[0]:rng[1]], axis=1) / 1000
        ax_fd.plot(u, F, 'b-', lw=0.8, label='ADINA (02241821)')
        ax_hist.plot(t_ad, u, 'b-', lw=0.8, label='ADINA (02241821)')

    if len(info['rd']) > 0:
        u_r = info['rd'][:, 1] * 1000
        t_ref = info['rd'][:, 0]
        F_r = np.sum(info['rr_data'][:, rng[0]:rng[1]], axis=1) / 1000
        ax_fd.plot(u_r, F_r, 'r--', lw=1.2, label='Reference')
        ax_hist.plot(t_ref, u_r, 'r--', lw=1.0, label='Reference')

    ax_fd.set_xlabel('Top Displacement (mm)')
    ax_fd.set_ylabel('Base Shear (kN)')
    ax_fd.set_title(f'{info["name"]} — Force vs Displacement', fontweight='bold')
    ax_fd.legend(fontsize=9, loc='best')
    ax_fd.grid(True, alpha=0.3)
    ax_fd.axhline(0, color='k', lw=0.5)
    ax_fd.axvline(0, color='k', lw=0.5)

    ax_hist.set_xlabel('Pseudo-time')
    ax_hist.set_ylabel('Top Displacement (mm)')
    ax_hist.set_title(f'{info["name"]} — Displacement History', fontweight='bold')
    ax_hist.legend(fontsize=9, loc='best')
    ax_hist.grid(True, alpha=0.3)

    n_ad = len(info['d'])
    n_ref = len(info['rd'])
    parts = [f'ADINA: {n_ad} pts']
    if n_ad > 0:
        parts[0] += f', u=[{info["d"][:,1].min()*1000:.1f}, {info["d"][:,1].max()*1000:.1f}] mm'
    parts.append(f'Ref:   {n_ref} pts')
    if n_ref > 0:
        parts[-1] += f', u=[{info["rd"][:,1].min()*1000:.1f}, {info["rd"][:,1].max()*1000:.1f}] mm'
    ax_fd.text(0.02, 0.02, '\n'.join(parts), transform=ax_fd.transAxes, fontsize=7,
               va='bottom', ha='left', bbox=dict(boxstyle='round,pad=0.3',
               facecolor='lightyellow', alpha=0.8))

plt.tight_layout(rect=[0, 0, 1, 0.97])
out_path = os.path.join(root, 'opensees_02241821_comparison.png')
plt.savefig(out_path, dpi=150, bbox_inches='tight')
print(f"\nPlot saved: {out_path}")
plt.close()
