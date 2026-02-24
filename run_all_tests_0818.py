# -*- coding: utf-8 -*-
"""Run all test cases with OpenSees-02230818 and generate comparison plots."""
import subprocess, os, time, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root = r'e:\Basic\Concrete_Model\ADINA'
exe = os.path.join(root, 'OpenSees-02230818.exe')

def run_case(name, tcl_script, wdir, output_files, timeout=3600):
    """Run an OpenSees test case."""
    print(f'\n{"="*60}')
    print(f'  Running: {name}')
    print(f'{"="*60}')
    for f in output_files:
        p = os.path.join(wdir, f)
        if os.path.exists(p):
            os.remove(p)
    logf = os.path.join(root, f'log_{name}_0818.txt')
    t0 = time.time()
    with open(logf, 'w') as lf:
        proc = subprocess.Popen([exe, tcl_script], cwd=wdir, stdout=lf, stderr=subprocess.STDOUT)
        proc.wait(timeout=timeout)
    elapsed = time.time() - t0
    print(f'  Exit: {proc.returncode}, time: {elapsed:.0f}s')
    with open(logf) as lf:
        lines = lf.readlines()
    for l in lines[-5:]:
        print(f'  {l.rstrip()}')
    for f in output_files:
        p = os.path.join(wdir, f)
        if os.path.exists(p):
            lc = sum(1 for _ in open(p))
            print(f'  {f}: {lc} lines')
        else:
            print(f'  {f}: NOT FOUND')
    return proc.returncode

# ── Run sw1-1 ── (skip if output exists)
sw1_1txt = os.path.join(root, 'sw1-1', '53.txt')
if not os.path.exists(sw1_1txt) or os.path.getsize(sw1_1txt) == 0:
    run_case('sw1-1', 'br_csp3.tcl',
             os.path.join(root, 'sw1-1'),
             ['53.txt', '1.txt'])
else:
    print(f'\n  sw1-1: using existing results')

# ── Run sw2-1 ── (skip if output exists)
sw2_1txt = os.path.join(root, 'sw2-1', '28.txt')
if not os.path.exists(sw2_1txt) or os.path.getsize(sw2_1txt) == 0:
    run_case('sw2-1', 'br_csp3.tcl',
             os.path.join(root, 'sw2-1'),
             ['28.txt', '1.txt'])
else:
    print(f'\n  sw2-1: using existing results')

# ── Multi-layer_Shell (already run, skip if output exists) ──
ml_dir = os.path.join(root, 'Multi-layer_Shell')
ml_disp = os.path.join(ml_dir, 'test_disp.txt')
ml_react = os.path.join(ml_dir, 'test_react.txt')
if os.path.exists(ml_disp) and os.path.exists(ml_react):
    lc_d = sum(1 for _ in open(ml_disp))
    lc_r = sum(1 for _ in open(ml_react))
    if lc_d > 700:
        print(f'\n  Multi-layer_Shell: using existing results (disp={lc_d}, react={lc_r} lines)')
    else:
        run_case('Multi-layer_Shell', 'test_br.tcl', ml_dir,
                 ['test_disp.txt', 'test_react.txt'], timeout=7200)
else:
    run_case('Multi-layer_Shell', 'test_br.tcl', ml_dir,
             ['test_disp.txt', 'test_react.txt'], timeout=7200)

# ══════════════════════════════════════════════════════════════
#  Data loading utilities
# ══════════════════════════════════════════════════════════════
def load_file(path, expected_ncols=None, skip_lines=0, max_disp=None):
    rows = []
    nc = expected_ncols
    with open(path) as f:
        all_lines = f.readlines()
    for line in all_lines[skip_lines:]:
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
        if max_disp is not None and nc == 2 and abs(vals[1]) > max_disp:
            continue
        rows.append(vals)
    return np.array(rows) if rows else None

def build_hysteresis_sw(wdir, disp_file, react_file, ncols_react, skip=10):
    """Build hysteresis for sw1-1 / sw2-1 (time-based matching)."""
    d = load_file(os.path.join(wdir, disp_file), expected_ncols=2, skip_lines=skip, max_disp=0.1)
    r = load_file(os.path.join(wdir, react_file), expected_ncols=ncols_react, skip_lines=skip)
    if d is None or r is None:
        return None, None
    t_d, disp_m = d[:, 0], d[:, 1]
    t_r = r[:, 0]
    base_shear_N = np.sum(r[:, 1:], axis=1)
    dmm, vkn = [], []
    for i in range(len(t_r)):
        idx = np.argmin(np.abs(t_d - t_r[i]))
        if np.abs(t_d[idx] - t_r[i]) < 0.5:
            dmm.append(disp_m[idx] * 1000.0)
            vkn.append(base_shear_N[i] / 1000.0)
    return np.array(dmm), np.array(vkn)

def build_hysteresis_ml(wdir):
    """Build hysteresis for Multi-layer_Shell (sign-corrected)."""
    d = load_file(os.path.join(wdir, 'test_disp.txt'), expected_ncols=2)
    r = load_file(os.path.join(wdir, 'test_react.txt'), expected_ncols=10)
    if d is None or r is None:
        return None, None
    t_d, disp_m = d[:, 0], d[:, 1]
    t_r = r[:, 0]
    base_shear_N = -np.sum(r[:, 1:], axis=1)
    dmm, vkn = [], []
    for i in range(len(t_r)):
        idx = np.argmin(np.abs(t_d - t_r[i]))
        if np.abs(t_d[idx] - t_r[i]) < 0.5:
            dmm.append(disp_m[idx] * 1000.0)
            vkn.append(base_shear_N[i] / 1000.0)
    return np.array(dmm), np.array(vkn)

def load_ref_ml(wdir):
    """Load Multi-layer_Shell reference (col1=V_kN, col2=disp_m)."""
    d = load_file(os.path.join(wdir, 'ref_disp1.txt'), expected_ncols=2)
    if d is None:
        return None, None
    return d[:, 1] * 1000.0, d[:, 0]

# ══════════════════════════════════════════════════════════════
#  Load all data
# ══════════════════════════════════════════════════════════════
print('\n' + '='*60)
print('  Loading data for plots...')
print('='*60)

results = {}

# sw1-1
sw1_dir = os.path.join(root, 'sw1-1')
br_d, br_v = build_hysteresis_sw(sw1_dir, '53.txt', '1.txt', ncols_react=6)
ref_d, ref_v = build_hysteresis_sw(sw1_dir, 'ref_53.txt', 'ref_1.txt', ncols_react=6)
results['sw1-1'] = {'br': (br_d, br_v), 'ref': (ref_d, ref_v)}
if br_d is not None:
    print(f'  sw1-1 B&R: {len(br_d)} pts, disp=[{br_d.min():.1f},{br_d.max():.1f}]mm, V=[{br_v.min():.0f},{br_v.max():.0f}]kN')
if ref_d is not None:
    print(f'  sw1-1 Ref: {len(ref_d)} pts, disp=[{ref_d.min():.1f},{ref_d.max():.1f}]mm, V=[{ref_v.min():.0f},{ref_v.max():.0f}]kN')

# sw2-1
sw2_dir = os.path.join(root, 'sw2-1')
br_d, br_v = build_hysteresis_sw(sw2_dir, '28.txt', '1.txt', ncols_react=6)
ref_d, ref_v = build_hysteresis_sw(sw2_dir, 'ref_28.txt', 'ref_1.txt', ncols_react=6)
results['sw2-1'] = {'br': (br_d, br_v), 'ref': (ref_d, ref_v)}
if br_d is not None:
    print(f'  sw2-1 B&R: {len(br_d)} pts, disp=[{br_d.min():.1f},{br_d.max():.1f}]mm, V=[{br_v.min():.0f},{br_v.max():.0f}]kN')
if ref_d is not None:
    print(f'  sw2-1 Ref: {len(ref_d)} pts, disp=[{ref_d.min():.1f},{ref_d.max():.1f}]mm, V=[{ref_v.min():.0f},{ref_v.max():.0f}]kN')

# Multi-layer_Shell
br_d, br_v = build_hysteresis_ml(ml_dir)
ref_d, ref_v = load_ref_ml(ml_dir)
results['Multi-layer Shell'] = {'br': (br_d, br_v), 'ref': (ref_d, ref_v)}
if br_d is not None:
    print(f'  ML   B&R: {len(br_d)} pts, disp=[{br_d.min():.1f},{br_d.max():.1f}]mm, V=[{br_v.min():.0f},{br_v.max():.0f}]kN')
if ref_d is not None:
    print(f'  ML   Ref: {len(ref_d)} pts, disp=[{ref_d.min():.1f},{ref_d.max():.1f}]mm, V=[{ref_v.min():.0f},{ref_v.max():.0f}]kN')

# ══════════════════════════════════════════════════════════════
#  Plot: 3-panel comparison
# ══════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

for ax, name in zip(axes, ['sw1-1', 'sw2-1', 'Multi-layer Shell']):
    br_d, br_v = results[name]['br']
    ref_d, ref_v = results[name]['ref']
    if ref_d is not None:
        ax.plot(ref_d, ref_v, color='gray', lw=1.4, alpha=0.7, label='Reference')
    if br_d is not None:
        ax.plot(br_d, br_v, 'b-', lw=1.2, label='B&R')
    ax.set_xlabel('Top Displacement (mm)', fontsize=12)
    ax.set_ylabel('Base Shear (kN)', fontsize=12)
    ax.set_title(name, fontsize=14, fontweight='bold')
    ax.legend(fontsize=11, loc='lower right')
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='k', lw=0.5)
    ax.axvline(x=0, color='k', lw=0.5)

plt.suptitle('OpenSees-02230818 (B&R) vs Reference — All Test Cases',
             fontsize=15, fontweight='bold', y=1.02)
plt.tight_layout()
out = os.path.join(root, 'all_cases_comparison_0818.png')
plt.savefig(out, dpi=150, bbox_inches='tight')
print(f'\nSaved: {out}')
plt.close()

# ══════════════════════════════════════════════════════════════
#  Individual high-resolution plots
# ══════════════════════════════════════════════════════════════
for name in ['sw1-1', 'sw2-1', 'Multi-layer Shell']:
    br_d, br_v = results[name]['br']
    ref_d, ref_v = results[name]['ref']
    fig, ax = plt.subplots(figsize=(10, 7))
    if ref_d is not None:
        ax.plot(ref_d, ref_v, color='gray', lw=1.5, alpha=0.7, label='Reference')
    if br_d is not None:
        ax.plot(br_d, br_v, 'b-', lw=1.2, label='B&R')
    ax.set_xlabel('Top Displacement (mm)', fontsize=13)
    ax.set_ylabel('Base Shear (kN)', fontsize=13)
    ax.set_title(f'{name} — Hysteresis Comparison', fontsize=14, fontweight='bold')
    ax.legend(fontsize=12, loc='lower right')
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='k', lw=0.5)
    ax.axvline(x=0, color='k', lw=0.5)
    plt.tight_layout()
    safe_name = name.replace(' ', '_').replace('-', '_')
    out = os.path.join(root, f'{safe_name}_comparison_0818.png')
    plt.savefig(out, dpi=150)
    print(f'Saved: {out}')
    plt.close()

print('\n*** All tests and plots complete! ***')
