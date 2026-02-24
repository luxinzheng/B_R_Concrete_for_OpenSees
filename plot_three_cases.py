# -*- coding: utf-8 -*-
"""
Plot comparison of B&R vs Reference for three test cases,
using pre-existing data files.
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei']
matplotlib.rcParams['axes.unicode_minus'] = False

ROOT = os.path.dirname(os.path.abspath(__file__))


def load_data(filepath):
    if not os.path.isfile(filepath):
        print(f'  [MISS] {filepath}')
        return None
    try:
        rows = []
        ncols = None
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split()
                try:
                    vals = [float(x) for x in parts]
                except ValueError:
                    continue
                if ncols is None:
                    ncols = len(vals)
                if len(vals) == ncols:
                    rows.append(vals)
        if len(rows) < 2:
            print(f'  [EMPTY] {filepath} ({len(rows)} valid rows)')
            return None
        d = np.array(rows)
        print(f'  [OK] {filepath}  shape={d.shape}')
        return d
    except Exception as e:
        print(f'  [ERR] {filepath}: {e}')
        return None


def get_hysteresis(disp_data, force_data, disp_col, force_cols):
    if disp_data is None or force_data is None:
        return None, None

    t_d = disp_data[:, 0]
    t_f = force_data[:, 0]
    disp_raw = disp_data[:, disp_col]
    shear_raw = np.zeros(force_data.shape[0])
    for c in force_cols:
        if c < force_data.shape[1]:
            shear_raw += force_data[:, c]

    # Align by matching time steps
    i_d, i_f = 0, 0
    disp_list, shear_list = [], []
    while i_d < len(t_d) and i_f < len(t_f):
        if abs(t_d[i_d] - t_f[i_f]) < 1e-8:
            disp_list.append(disp_raw[i_d])
            shear_list.append(shear_raw[i_f])
            i_d += 1
            i_f += 1
        elif t_d[i_d] < t_f[i_f]:
            i_d += 1
        else:
            i_f += 1

    if len(disp_list) < 2:
        n = min(len(disp_raw), len(shear_raw))
        disp_list = list(disp_raw[:n])
        shear_list = list(shear_raw[:n])

    disp = np.array(disp_list)
    shear = np.array(shear_list)

    # Filter divergent points (displacement jump > 100x typical increment)
    if len(disp) > 10:
        diffs = np.abs(np.diff(disp))
        median_diff = np.median(diffs[diffs > 1e-15]) if np.any(diffs > 1e-15) else 1e-3
        threshold = max(median_diff * 100, 0.1)
        good = np.ones(len(disp), dtype=bool)
        for j in range(1, len(disp)):
            if abs(disp[j] - disp[j-1]) > threshold:
                good[j:] = False
                break
        disp = disp[good]
        shear = shear[good]

    return disp, shear


def compute_envelope(disp, force):
    if disp is None or force is None or len(disp) < 3:
        return None
    pos_d, pos_f, neg_d, neg_f = [], [], [], []
    max_d_pos, max_d_neg = 0.0, 0.0
    n = len(disp)
    for i in range(n):
        d, f = disp[i], force[i]
        if d > max_d_pos + 1e-10:
            max_d_pos = d
            pos_d.append(d)
            pos_f.append(f)
        elif d < max_d_neg - 1e-10:
            max_d_neg = d
            neg_d.append(d)
            neg_f.append(f)
    env = []
    for d, f in zip(reversed(neg_d), reversed(neg_f)):
        env.append([d, f])
    env.append([0.0, 0.0])
    for d, f in zip(pos_d, pos_f):
        env.append([d, f])
    return np.array(env) if len(env) >= 3 else None


CASES = [
    {
        'name': 'Multi-layer_Shell',
        'title': 'Multi-layer Shell Wall',
        'dir': os.path.join(ROOT, 'Multi-layer_Shell'),
        'br_disp': 'br_disp1.txt',
        'br_force': 'br_shearforce1.txt',
        'ref_disp': 'ref_disp1.txt',
        'ref_force': 'ref_shearforce1.txt',
        'disp_col': 1,
        'force_cols': list(range(1, 10)),
    },
    {
        'name': 'sw1-1',
        'title': 'Shear Wall SW1-1 (fc=20.7 MPa)',
        'dir': os.path.join(ROOT, 'sw1-1'),
        'br_disp': 'br_53.txt',
        'br_force': 'br_1.txt',
        'ref_disp': 'ref_53.txt',
        'ref_force': 'ref_1.txt',
        'disp_col': 1,
        'force_cols': list(range(1, 6)),
    },
    {
        'name': 'sw2-1',
        'title': 'Shear Wall SW2-1 (fc=30.8 MPa)',
        'dir': os.path.join(ROOT, 'sw2-1'),
        'br_disp': 'br_28.txt',
        'br_force': 'br_1.txt',
        'ref_disp': 'ref_28.txt',
        'ref_force': 'ref_1.txt',
        'disp_col': 1,
        'force_cols': list(range(1, 6)),
    },
]


def main():
    print('Loading data...\n')
    results = []
    for case in CASES:
        d = case['dir']
        print(f'=== {case["name"]} ===')
        br_d = load_data(os.path.join(d, case['br_disp']))
        br_f = load_data(os.path.join(d, case['br_force']))
        ref_d = load_data(os.path.join(d, case['ref_disp']))
        ref_f = load_data(os.path.join(d, case['ref_force']))

        br_disp, br_force = get_hysteresis(br_d, br_f, case['disp_col'], case['force_cols'])
        ref_disp, ref_force = get_hysteresis(ref_d, ref_f, case['disp_col'], case['force_cols'])
        results.append((br_disp, br_force, ref_disp, ref_force))
        print()

    # --- Individual case plots ---
    for i, case in enumerate(CASES):
        br_disp, br_force, ref_disp, ref_force = results[i]
        if br_disp is None and ref_disp is None:
            print(f'Skipping {case["name"]}: no data')
            continue

        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        fig.suptitle(case['title'], fontsize=14, fontweight='bold')

        ax = axes[0]
        if ref_disp is not None:
            ax.plot(ref_disp * 1e3, ref_force / 1e3, 'b-', lw=0.8,
                    label='Ref', alpha=0.7)
        if br_disp is not None:
            ax.plot(br_disp * 1e3, br_force / 1e3, 'r-', lw=0.8,
                    label='B&R', alpha=0.7)
        ax.set_xlabel('Displacement (mm)')
        ax.set_ylabel('Base Shear (kN)')
        ax.set_title('Hysteresis Comparison')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0, color='k', lw=0.5)
        ax.axvline(x=0, color='k', lw=0.5)

        ax = axes[1]
        if ref_disp is not None:
            env = compute_envelope(ref_disp, ref_force)
            if env is not None:
                ax.plot(env[:, 0] * 1e3, env[:, 1] / 1e3, 'b-o',
                        ms=3, lw=1.2, label='Ref')
        if br_disp is not None:
            env = compute_envelope(br_disp, br_force)
            if env is not None:
                ax.plot(env[:, 0] * 1e3, env[:, 1] / 1e3, 'r-s',
                        ms=3, lw=1.2, label='B&R')
        ax.set_xlabel('Displacement (mm)')
        ax.set_ylabel('Base Shear (kN)')
        ax.set_title('Skeleton Curve')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0, color='k', lw=0.5)
        ax.axvline(x=0, color='k', lw=0.5)

        plt.tight_layout()
        out = os.path.join(ROOT, f'compare_{case["name"]}.png')
        fig.savefig(out, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f'Saved: {out}')

    # --- Summary 3-panel plot ---
    fig, axes = plt.subplots(1, 3, figsize=(21, 6))
    fig.suptitle('B&R (02210959) vs Reference Model Comparison', fontsize=14, fontweight='bold')
    for i, case in enumerate(CASES):
        br_disp, br_force, ref_disp, ref_force = results[i]
        ax = axes[i]
        if ref_disp is not None:
            ax.plot(ref_disp * 1e3, ref_force / 1e3, 'b-', lw=0.6,
                    label='Ref', alpha=0.7)
        if br_disp is not None:
            ax.plot(br_disp * 1e3, br_force / 1e3, 'r-', lw=0.6,
                    label='B&R', alpha=0.7)
        ax.set_xlabel('Displacement (mm)')
        ax.set_ylabel('Base Shear (kN)')
        ax.set_title(case['title'])
        if br_disp is not None or ref_disp is not None:
            ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0, color='k', lw=0.5)
        ax.axvline(x=0, color='k', lw=0.5)
    plt.tight_layout()
    out = os.path.join(ROOT, 'compare_all_three.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved: {out}')

    # --- Quantitative report ---
    print('\n' + '=' * 70)
    print('QUANTITATIVE COMPARISON REPORT')
    print('=' * 70)
    for i, case in enumerate(CASES):
        br_disp, br_force, ref_disp, ref_force = results[i]
        print(f'\n--- {case["name"]} ({case["title"]}) ---')

        if br_disp is not None:
            print(f'  B&R:  {len(br_disp)} steps, '
                  f'disp=[{br_disp.min()*1e3:.3f}, {br_disp.max()*1e3:.3f}] mm, '
                  f'peak +F={br_force.max()/1e3:.2f} kN, '
                  f'peak -F={br_force.min()/1e3:.2f} kN')
        else:
            print(f'  B&R:  NO DATA (analysis failed)')

        if ref_disp is not None:
            print(f'  Ref:  {len(ref_disp)} steps, '
                  f'disp=[{ref_disp.min()*1e3:.3f}, {ref_disp.max()*1e3:.3f}] mm, '
                  f'peak +F={ref_force.max()/1e3:.2f} kN, '
                  f'peak -F={ref_force.min()/1e3:.2f} kN')
        else:
            print(f'  Ref:  NO DATA')

        if br_disp is not None and ref_disp is not None:
            r_pos = (br_force.max() - ref_force.max()) / abs(ref_force.max()) * 100
            r_neg = (abs(br_force.min()) - abs(ref_force.min())) / abs(ref_force.min()) * 100
            print(f'  Peak +F diff: {r_pos:+.1f}%')
            print(f'  Peak -F diff: {r_neg:+.1f}%')
            print(f'  B&R completion: {len(br_disp)}/{len(ref_disp)} steps '
                  f'({len(br_disp)/len(ref_disp)*100:.0f}%)')


if __name__ == '__main__':
    main()
