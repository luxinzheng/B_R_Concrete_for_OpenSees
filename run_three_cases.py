# -*- coding: utf-8 -*-
"""
Run three test cases with both B&R and Reference OpenSees executables,
then plot and compare results.
"""
import subprocess, os, sys, shutil, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei']
matplotlib.rcParams['axes.unicode_minus'] = False

ROOT = os.path.dirname(os.path.abspath(__file__))
BR_EXE = os.path.join(ROOT, 'OpenSees-02210959.exe')
REF_EXE = None
for f in os.listdir(ROOT):
    if f.startswith('OpenSees-') and f.endswith('.exe') and '\u53c2\u8003' in f:
        REF_EXE = os.path.join(ROOT, f)
        break
if REF_EXE is None:
    for f in os.listdir(ROOT):
        if f.startswith('OpenSees-') and f.endswith('.exe') and '\u6a21\u578b' in f:
            REF_EXE = os.path.join(ROOT, f)
            break

CASES = {
    'Multi-layer_Shell': {
        'dir': os.path.join(ROOT, 'Multi-layer_Shell'),
        'br_tcl': 'br_main.tcl',
        'ref_tcl': 'main.tcl',
        'disp_file': 'disp1.txt',
        'force_file': 'shearforce1.txt',
        'disp_col': 1,       # column index for displacement (0=time)
        'force_cols': list(range(1, 10)),  # nodes 1-9, columns 1..9
        'title': 'Multi-layer Shell Wall',
    },
    'sw1-1': {
        'dir': os.path.join(ROOT, 'sw1-1'),
        'br_tcl': 'br_csp3.tcl',
        'ref_tcl': 'csp3.tcl',
        'disp_file': '53.txt',
        'force_file': '1.txt',
        'disp_col': 1,
        'force_cols': list(range(1, 6)),  # nodes 1-5, columns 1..5
        'title': 'Shear Wall SW1-1',
    },
    'sw2-1': {
        'dir': os.path.join(ROOT, 'sw2-1'),
        'br_tcl': 'br_csp3.tcl',
        'ref_tcl': 'csp3.tcl',
        'disp_file': '28.txt',
        'force_file': '1.txt',
        'disp_col': 1,
        'force_cols': list(range(1, 6)),  # nodes 1-5, columns 1..5
        'title': 'Shear Wall SW2-1',
    },
}


def run_opensees(exe, tcl, work_dir, label=''):
    """Run an OpenSees executable with a tcl script in work_dir."""
    print(f'  Running {label}: {os.path.basename(exe)} + {tcl} ...')
    t0 = time.time()
    try:
        p = subprocess.run(
            [exe, tcl],
            cwd=work_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=600,
        )
        stdout = p.stdout.decode('utf-8', errors='replace')
        stderr = p.stderr.decode('utf-8', errors='replace')
        elapsed = time.time() - t0
        print(f'    Done in {elapsed:.1f}s, exit={p.returncode}')
        if 'FAILED' in stdout:
            for line in stdout.splitlines():
                if 'FAILED' in line or 'failed' in line.lower():
                    print(f'    WARNING: {line.strip()}')
        return p.returncode, stdout, stderr
    except subprocess.TimeoutExpired:
        print(f'    TIMEOUT after 600s!')
        return -1, '', 'timeout'


def load_data(filepath):
    """Load space-delimited data, return array or None."""
    if not os.path.isfile(filepath):
        return None
    try:
        d = np.loadtxt(filepath)
        if d.size == 0 or d.ndim < 2 or d.shape[0] < 2:
            return None
        return d
    except Exception:
        return None


def extract_hysteresis(case_cfg, prefix, work_dir):
    """Extract displacement and total base shear from output files.
    prefix is 'br_' or 'ref_' for naming saved copies.
    Returns (disp_array, force_array) or (None, None).
    """
    disp_path = os.path.join(work_dir, case_cfg['disp_file'])
    force_path = os.path.join(work_dir, case_cfg['force_file'])

    d_disp = load_data(disp_path)
    d_force = load_data(force_path)

    if d_disp is None or d_force is None:
        return None, None

    disp = d_disp[:, case_cfg['disp_col']]
    base_shear = np.zeros(d_force.shape[0])
    for c in case_cfg['force_cols']:
        if c < d_force.shape[1]:
            base_shear += d_force[:, c]

    return disp, base_shear


def save_copy(src, dst):
    """Copy src to dst if src exists."""
    if os.path.isfile(src):
        shutil.copy2(src, dst)


def run_case(name, cfg):
    """Run both B&R and Reference for a single case."""
    work_dir = cfg['dir']
    print(f'\n{"="*60}')
    print(f'Case: {name} ({cfg["title"]})')
    print(f'{"="*60}')

    br_disp, br_force = None, None
    ref_disp, ref_force = None, None

    # --- Run B&R model ---
    rc, stdout, stderr = run_opensees(BR_EXE, cfg['br_tcl'], work_dir, 'B&R')
    br_disp_raw, br_force_raw = extract_hysteresis(cfg, 'br_', work_dir)
    if br_disp_raw is not None:
        br_disp, br_force = br_disp_raw, br_force_raw
        save_copy(os.path.join(work_dir, cfg['disp_file']),
                  os.path.join(work_dir, 'br_' + cfg['disp_file']))
        save_copy(os.path.join(work_dir, cfg['force_file']),
                  os.path.join(work_dir, 'br_' + cfg['force_file']))
        print(f'    B&R: {len(br_disp)} data points')
    else:
        print(f'    B&R: NO OUTPUT DATA')

    # --- Run Reference model ---
    rc, stdout, stderr = run_opensees(REF_EXE, cfg['ref_tcl'], work_dir, 'Ref')
    ref_disp_raw, ref_force_raw = extract_hysteresis(cfg, 'ref_', work_dir)
    if ref_disp_raw is not None:
        ref_disp, ref_force = ref_disp_raw, ref_force_raw
        save_copy(os.path.join(work_dir, cfg['disp_file']),
                  os.path.join(work_dir, 'ref_' + cfg['disp_file']))
        save_copy(os.path.join(work_dir, cfg['force_file']),
                  os.path.join(work_dir, 'ref_' + cfg['force_file']))
        print(f'    Ref: {len(ref_disp)} data points')
    else:
        # Try loading pre-existing ref data
        ref_disp_path = os.path.join(work_dir, 'ref_' + cfg['disp_file'])
        ref_force_path = os.path.join(work_dir, 'ref_' + cfg['force_file'])
        d_disp = load_data(ref_disp_path)
        d_force = load_data(ref_force_path)
        if d_disp is not None and d_force is not None:
            ref_disp = d_disp[:, cfg['disp_col']]
            ref_force = np.zeros(d_force.shape[0])
            for c in cfg['force_cols']:
                if c < d_force.shape[1]:
                    ref_force += d_force[:, c]
            print(f'    Ref: {len(ref_disp)} data points (from saved ref files)')
        else:
            print(f'    Ref: NO OUTPUT DATA')

    return br_disp, br_force, ref_disp, ref_force


def analyze_and_report(name, cfg, br_disp, br_force, ref_disp, ref_force):
    """Quantitative comparison."""
    report_lines = []
    report_lines.append(f'\n--- {name} ({cfg["title"]}) ---')

    if br_disp is not None:
        report_lines.append(f'  B&R: {len(br_disp)} steps')
        report_lines.append(f'    Disp range: [{br_disp.min()*1e3:.3f}, {br_disp.max()*1e3:.3f}] mm')
        report_lines.append(f'    Peak +force: {br_force.max()/1e3:.2f} kN')
        report_lines.append(f'    Peak -force: {br_force.min()/1e3:.2f} kN')
    else:
        report_lines.append(f'  B&R: No data')

    if ref_disp is not None:
        report_lines.append(f'  Ref: {len(ref_disp)} steps')
        report_lines.append(f'    Disp range: [{ref_disp.min()*1e3:.3f}, {ref_disp.max()*1e3:.3f}] mm')
        report_lines.append(f'    Peak +force: {ref_force.max()/1e3:.2f} kN')
        report_lines.append(f'    Peak -force: {ref_force.min()/1e3:.2f} kN')
    else:
        report_lines.append(f'  Ref: No data')

    if br_disp is not None and ref_disp is not None:
        br_peak_pos = br_force.max()
        ref_peak_pos = ref_force.max()
        br_peak_neg = br_force.min()
        ref_peak_neg = ref_force.min()
        if abs(ref_peak_pos) > 1e-6:
            ratio_pos = (br_peak_pos - ref_peak_pos) / abs(ref_peak_pos) * 100
            report_lines.append(f'  Peak +force diff: {ratio_pos:+.1f}%')
        if abs(ref_peak_neg) > 1e-6:
            ratio_neg = (abs(br_peak_neg) - abs(ref_peak_neg)) / abs(ref_peak_neg) * 100
            report_lines.append(f'  Peak -force diff: {ratio_neg:+.1f}%')

    return '\n'.join(report_lines)


def plot_case(name, cfg, br_disp, br_force, ref_disp, ref_force, output_dir):
    """Plot hysteresis comparison for a case."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    fig.suptitle(cfg['title'], fontsize=14, fontweight='bold')

    # Left: overlay
    ax = axes[0]
    if ref_disp is not None:
        ax.plot(ref_disp * 1e3, ref_force / 1e3, 'b-', linewidth=0.8,
                label='Ref', alpha=0.8)
    if br_disp is not None:
        ax.plot(br_disp * 1e3, br_force / 1e3, 'r-', linewidth=0.8,
                label='B&R', alpha=0.8)
    ax.set_xlabel('Displacement (mm)')
    ax.set_ylabel('Base Shear (kN)')
    ax.set_title('Hysteresis Comparison')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='k', linewidth=0.5)
    ax.axvline(x=0, color='k', linewidth=0.5)

    # Right: envelope
    ax = axes[1]
    if ref_disp is not None:
        env_ref = compute_envelope(ref_disp, ref_force)
        if env_ref is not None:
            ax.plot(env_ref[:, 0] * 1e3, env_ref[:, 1] / 1e3, 'b-o',
                    markersize=3, linewidth=1.2, label='Ref')
    if br_disp is not None:
        env_br = compute_envelope(br_disp, br_force)
        if env_br is not None:
            ax.plot(env_br[:, 0] * 1e3, env_br[:, 1] / 1e3, 'r-s',
                    markersize=3, linewidth=1.2, label='B&R')
    ax.set_xlabel('Displacement (mm)')
    ax.set_ylabel('Base Shear (kN)')
    ax.set_title('Skeleton Curve')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='k', linewidth=0.5)
    ax.axvline(x=0, color='k', linewidth=0.5)

    plt.tight_layout()
    out_path = os.path.join(output_dir, f'compare_{name}.png')
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  Plot saved: {out_path}')
    return out_path


def compute_envelope(disp, force):
    """Compute skeleton/envelope curve from hysteresis data.
    Picks peak force at each new max/min displacement level.
    """
    if disp is None or force is None or len(disp) < 3:
        return None

    pos_peaks_d = []
    pos_peaks_f = []
    neg_peaks_d = []
    neg_peaks_f = []

    max_d_pos = 0.0
    max_d_neg = 0.0

    n = len(disp)
    i = 0
    while i < n:
        d = disp[i]
        f = force[i]
        if d > max_d_pos + 1e-10:
            max_d_pos = d
            j = i
            f_peak = f
            while j < n and disp[j] >= max_d_pos - 1e-10:
                if force[j] > f_peak:
                    f_peak = force[j]
                j += 1
            pos_peaks_d.append(max_d_pos)
            pos_peaks_f.append(f_peak)
        elif d < max_d_neg - 1e-10:
            max_d_neg = d
            j = i
            f_peak = f
            while j < n and disp[j] <= max_d_neg + 1e-10:
                if force[j] < f_peak:
                    f_peak = force[j]
                j += 1
            neg_peaks_d.append(max_d_neg)
            neg_peaks_f.append(f_peak)
        i += 1

    env = []
    for d, f in zip(reversed(neg_peaks_d), reversed(neg_peaks_f)):
        env.append([d, f])
    env.append([0.0, 0.0])
    for d, f in zip(pos_peaks_d, pos_peaks_f):
        env.append([d, f])

    if len(env) < 3:
        return None
    return np.array(env)


def main():
    print(f'B&R exe: {BR_EXE}')
    print(f'Ref exe: {REF_EXE}')
    if not os.path.isfile(BR_EXE):
        print('ERROR: B&R exe not found!')
        return
    if REF_EXE is None or not os.path.isfile(REF_EXE):
        print('ERROR: Ref exe not found!')
        return

    output_dir = ROOT
    all_reports = []

    for name, cfg in CASES.items():
        br_disp, br_force, ref_disp, ref_force = run_case(name, cfg)
        report = analyze_and_report(name, cfg, br_disp, br_force, ref_disp, ref_force)
        all_reports.append(report)
        print(report)
        plot_case(name, cfg, br_disp, br_force, ref_disp, ref_force, output_dir)

    # Summary plot: all three cases in one figure
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    fig.suptitle('Three Test Cases Comparison: B&R vs Reference', fontsize=14, fontweight='bold')

    for idx, (name, cfg) in enumerate(CASES.items()):
        ax = axes[idx]
        work_dir = cfg['dir']

        # Load saved data
        ref_disp_path = os.path.join(work_dir, 'ref_' + cfg['disp_file'])
        ref_force_path = os.path.join(work_dir, 'ref_' + cfg['force_file'])
        br_disp_path = os.path.join(work_dir, 'br_' + cfg['disp_file'])
        br_force_path = os.path.join(work_dir, 'br_' + cfg['force_file'])

        d = load_data(ref_disp_path)
        f = load_data(ref_force_path)
        if d is not None and f is not None:
            rd = d[:, cfg['disp_col']]
            rf = np.zeros(f.shape[0])
            for c in cfg['force_cols']:
                if c < f.shape[1]:
                    rf += f[:, c]
            ax.plot(rd * 1e3, rf / 1e3, 'b-', linewidth=0.6, label='Ref', alpha=0.7)

        d = load_data(br_disp_path)
        f = load_data(br_force_path)
        if d is not None and f is not None:
            bd = d[:, cfg['disp_col']]
            bf = np.zeros(f.shape[0])
            for c in cfg['force_cols']:
                if c < f.shape[1]:
                    bf += f[:, c]
            ax.plot(bd * 1e3, bf / 1e3, 'r-', linewidth=0.6, label='B&R', alpha=0.7)

        ax.set_xlabel('Displacement (mm)')
        ax.set_ylabel('Base Shear (kN)')
        ax.set_title(cfg['title'])
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0, color='k', linewidth=0.5)
        ax.axvline(x=0, color='k', linewidth=0.5)

    plt.tight_layout()
    summary_path = os.path.join(output_dir, 'compare_all_three.png')
    fig.savefig(summary_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'\nSummary plot saved: {summary_path}')

    print('\n' + '=' * 60)
    print('SUMMARY REPORT')
    print('=' * 60)
    for r in all_reports:
        print(r)


if __name__ == '__main__':
    main()
