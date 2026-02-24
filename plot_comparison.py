#!/usr/bin/env python3
"""Plot comparison of B&R vs Reference for all available test cases."""
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

ROOT = os.path.dirname(os.path.abspath(__file__))

def load_data(filepath):
    try:
        with open(filepath, 'r') as f:
            raw = f.readlines()
        if not raw:
            return None
        rows, ncols = [], None
        for line in raw:
            parts = line.split()
            if not parts:
                continue
            vals = []
            for p in parts:
                try:
                    vals.append(float(p))
                except ValueError:
                    break
            if not vals:
                continue
            if ncols is None:
                ncols = len(vals)
            if len(vals) == ncols:
                rows.append(vals)
        if not rows:
            return None
        return np.array(rows)
    except Exception:
        return None

def parse_hysteresis(disp_file, force_file):
    dd = load_data(disp_file)
    fd = load_data(force_file)
    if dd is None or fd is None:
        return None, None
    disp = dd[:, 1]
    n = min(len(disp), len(fd))
    force = np.sum(fd[:n, 1:], axis=1)
    disp = disp[:n]
    # Remove trailing crash-corrupted rows (displacement jumps > 100x median)
    if len(disp) > 10:
        med = np.median(np.abs(disp[disp != 0])) if np.any(disp != 0) else 1
        while len(disp) > 1 and abs(disp[-1]) > max(100 * med, 1.0):
            disp = disp[:-1]
            force = force[:len(disp)]
    return disp, force

CASES = [
    {
        'name': 'Multi-layer_Shell',
        'title': 'Multi-layer Shell',
        'ref_disp': os.path.join(ROOT, 'Multi-layer_Shell', 'ref_disp1.txt'),
        'ref_force': os.path.join(ROOT, 'Multi-layer_Shell', 'ref_shearforce1.txt'),
        'br_disp': os.path.join(ROOT, 'Multi-layer_Shell', 'test_disp.txt'),
        'br_force': os.path.join(ROOT, 'Multi-layer_Shell', 'test_react.txt'),
    },
    {
        'name': 'sw1-1',
        'title': 'SW1-1',
        'ref_disp': os.path.join(ROOT, 'sw1-1', 'ref_53.txt'),
        'ref_force': os.path.join(ROOT, 'sw1-1', 'ref_1.txt'),
        'br_disp': os.path.join(ROOT, 'sw1-1', 'test_disp.txt'),
        'br_force': os.path.join(ROOT, 'sw1-1', 'test_react.txt'),
    },
    {
        'name': 'sw2-1',
        'title': 'SW2-1',
        'ref_disp': os.path.join(ROOT, 'sw2-1', 'ref_28.txt'),
        'ref_force': os.path.join(ROOT, 'sw2-1', 'ref_1.txt'),
        'br_disp': os.path.join(ROOT, 'sw2-1', 'test_disp.txt'),
        'br_force': os.path.join(ROOT, 'sw2-1', 'test_react.txt'),
    },
]

results = []
for c in CASES:
    ref_d, ref_f = parse_hysteresis(c['ref_disp'], c['ref_force'])
    br_d, br_f = parse_hysteresis(c['br_disp'], c['br_force'])
    results.append({
        'name': c['name'], 'title': c['title'],
        'ref_disp': ref_d, 'ref_force': ref_f,
        'br_disp': br_d, 'br_force': br_f,
    })

# ── Summary overview (3 subplots) ──
fig, axes = plt.subplots(1, 3, figsize=(21, 6))
for idx, r in enumerate(results):
    ax = axes[idx]
    has_ref = r['ref_disp'] is not None
    has_br = r['br_disp'] is not None

    if has_ref:
        ax.plot(r['ref_disp'] * 1e3, r['ref_force'] / 1e3,
                'b-', linewidth=0.7, alpha=0.8, label='Reference')
    if has_br:
        ax.plot(r['br_disp'] * 1e3, r['br_force'] / 1e3,
                'r-', linewidth=0.7, alpha=0.8, label='B&R')

    ax.set_xlabel('Displacement (mm)')
    ax.set_ylabel('Base Shear (kN)')
    ax.set_title(r['title'])
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='k', linewidth=0.5)
    ax.axvline(x=0, color='k', linewidth=0.5)
    if has_ref or has_br:
        ax.legend(loc='best', fontsize=9)

    ref_pk = np.max(np.abs(r['ref_force'])) / 1e3 if has_ref else 0
    br_pk = np.max(np.abs(r['br_force'])) / 1e3 if has_br else 0
    ref_n = len(r['ref_disp']) if has_ref else 0
    br_n = len(r['br_disp']) if has_br else 0
    txt = f'Ref: {ref_n} steps, peak={ref_pk:.1f}kN\n'
    txt += f'B&R: {br_n} steps, peak={br_pk:.1f}kN'
    if ref_pk > 0 and br_pk > 0:
        txt += f'\nRatio: {br_pk/ref_pk:.3f}'
    ax.text(0.02, 0.02, txt, transform=ax.transAxes, fontsize=8,
            verticalalignment='bottom',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
out = os.path.join(ROOT, 'comparison_3cases.png')
fig.savefig(out, dpi=150, bbox_inches='tight')
print(f'Saved: {out}')
plt.close()

# ── Individual detailed plots ──
for r in results:
    has_ref = r['ref_disp'] is not None
    has_br = r['br_disp'] is not None
    if not has_ref and not has_br:
        continue

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    ax = axes[0]
    if has_ref:
        ax.plot(r['ref_disp'] * 1e3, r['ref_force'] / 1e3,
                'b-', linewidth=0.7, alpha=0.8, label='Reference')
    if has_br:
        ax.plot(r['br_disp'] * 1e3, r['br_force'] / 1e3,
                'r-', linewidth=0.7, alpha=0.8, label='B&R')
    ax.set_xlabel('Displacement (mm)')
    ax.set_ylabel('Base Shear (kN)')
    ax.set_title(f'{r["title"]} - Force-Displacement Hysteresis')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='k', linewidth=0.5)
    ax.axvline(x=0, color='k', linewidth=0.5)

    ax = axes[1]
    if has_ref:
        ax.plot(np.arange(len(r['ref_disp'])), r['ref_disp'] * 1e3,
                'b-', linewidth=0.7, alpha=0.8, label='Reference')
    if has_br:
        ax.plot(np.arange(len(r['br_disp'])), r['br_disp'] * 1e3,
                'r-', linewidth=0.7, alpha=0.8, label='B&R')
    ax.set_xlabel('Step')
    ax.set_ylabel('Displacement (mm)')
    ax.set_title(f'{r["title"]} - Displacement History')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    out = os.path.join(ROOT, f'comparison_{r["name"]}.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    print(f'Saved: {out}')
    plt.close()

# ── Summary table ──
print('\n' + '=' * 85)
print('SUMMARY')
print('=' * 85)
print(f'{"Case":<25} {"Ref Steps":>10} {"B&R Steps":>10} '
      f'{"Ref Peak(kN)":>14} {"B&R Peak(kN)":>14} {"Ratio":>8}')
print('-' * 85)
for r in results:
    has_ref = r['ref_disp'] is not None
    has_br = r['br_disp'] is not None
    ref_pk = np.max(np.abs(r['ref_force'])) / 1e3 if has_ref else 0
    br_pk = np.max(np.abs(r['br_force'])) / 1e3 if has_br else 0
    ref_n = len(r['ref_disp']) if has_ref else 0
    br_n = len(r['br_disp']) if has_br else 0
    ratio = br_pk / ref_pk if ref_pk > 0 and br_pk > 0 else 0
    print(f'{r["name"]:<25} {ref_n:>10d} {br_n:>10d} '
          f'{ref_pk:>14.1f} {br_pk:>14.1f} {ratio:>8.3f}')

    if has_ref and has_br:
        n = min(len(r['ref_disp']), len(r['br_disp']))
        ref_dmax = np.max(np.abs(r['ref_disp'][:n])) * 1e3
        br_dmax = np.max(np.abs(r['br_disp'][:n])) * 1e3
        overlap_rms = np.sqrt(np.mean((r['ref_force'][:n] - r['br_force'][:n])**2)) / 1e3
        print(f'  {"":25} ref_dmax={ref_dmax:.2f}mm  br_dmax={br_dmax:.2f}mm  '
              f'RMS_diff={overlap_rms:.1f}kN')

print('\nDone.')
