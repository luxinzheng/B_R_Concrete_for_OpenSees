"""
Verification Suite: Run all tests with both ADINA and reference models,
then generate comparison plots.

Usage:
    python run_verification.py
    python run_verification.py --adina OpenSees-2017.exe --ref "OpenSees-ref.exe"
    python run_verification.py --tests uniaxial wall
"""
import os, sys, subprocess, argparse, glob, time, io
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR   = os.path.dirname(SCRIPT_DIR)

def find_exe(pattern):
    matches = sorted(glob.glob(os.path.join(ROOT_DIR, pattern)))
    return matches[-1] if matches else None

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--adina', default=None, help='ADINA OpenSees exe')
    p.add_argument('--ref', default=None, help='Reference OpenSees exe')
    p.add_argument('--tests', nargs='*', default=['uniaxial','wall','beam','plate'])
    return p.parse_args()

def run_opensees(exe, tcl_file, model_type, cwd, sub_test=None):
    wrapper = os.path.join(cwd, f'_run_{model_type}.tcl')
    with open(wrapper, 'w') as f:
        f.write(f'set model_type "{model_type}"\n')
        if sub_test:
            f.write(f'set sub_test "{sub_test}"\n')
        f.write(f'source {{{tcl_file}}}\n')
    cmd = [exe, wrapper]
    print(f'  Running: {os.path.basename(exe)} -> {os.path.basename(tcl_file)} [{model_type}]')
    t0 = time.time()
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, timeout=600)
    dt = time.time() - t0
    stdout = r.stdout.decode('utf-8', errors='replace') if r.stdout else ''
    stderr = r.stderr.decode('utf-8', errors='replace') if r.stderr else ''
    print(f'    exit={r.returncode}  time={dt:.1f}s')
    if r.returncode != 0:
        print(f'    STDERR: {stderr[:500]}')
    for line in stdout.splitlines():
        if any(k in line.upper() for k in ['FAIL','DONE','COMPLETE','WALL','BEAM','PLATE','GRAVITY']):
            print(f'    {line}')
    return r.returncode == 0

def load_txt(path, ncol=None):
    if not os.path.exists(path):
        print(f'  WARNING: {path} not found')
        return None
    try:
        d = np.loadtxt(path)
        if d.size == 0 or (d.ndim >= 1 and d.shape[0] == 0):
            return None
        if d.ndim == 1:
            d = d.reshape(1, -1)
        return d
    except Exception as e:
        print(f'  WARNING: {path}: {e}')
        return None

# ─── Plotting helpers ───
def setup_axes(ax, xlabel, ylabel, title):
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

def saenz_curve(E0, fc, epsc, fu, epsu, npts=200):
    eps = np.linspace(0, epsc, npts)
    R = E0 * epsc / fc
    A = (R * (epsc/epsc) - 2.0) / (1.0 - 2.0*R)  # simplified
    sig = np.zeros(npts)
    for i, e in enumerate(eps):
        x = e / epsc
        if abs(x) < 1e-12:
            sig[i] = 0.0
        else:
            sig[i] = fc * x * R / (1.0 + (R - 2.0)*x + x**2)
    return eps, sig

# ─── Plot functions ───
def plot_uniaxial(cwd):
    print('\n=== Plotting uniaxial results ===')
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    for sub, name in [(0, 'comp'), (1, 'tens'), (2, 'cyc')]:
        ax = axes[sub]
        titles = {'comp': 'Monotonic Compression', 'tens': 'Monotonic Tension', 'cyc': 'Cyclic Compression'}

        for mt, color, ls in [('adina', '#D62728', '-'), ('ref', '#1F77B4', '--')]:
            dpath = os.path.join(cwd, f'vuniax_{mt}_{name}_disp.txt')
            rpath = os.path.join(cwd, f'vuniax_{mt}_{name}_react.txt')
            dd = load_txt(dpath); rd = load_txt(rpath)
            if dd is None or rd is None:
                continue
            disp = dd[:, 1]  # node 3 uy
            force = np.sum(rd[:, 1:], axis=1)
            area = 1.0 * 0.1  # 1m width x 0.1m thick
            stress = force / area
            strain = disp  # eps = uy / L, L=1m
            label = 'B&R' if mt == 'adina' else 'Reference'
            ax.plot(strain * 1000, stress / 1e6, color=color, ls=ls, lw=1.5, label=label)

        setup_axes(ax, 'Strain (x1e-3)', 'Stress (MPa)', titles[name])
        ax.legend(fontsize=9)
        ax.axhline(0, color='k', lw=0.5)
        ax.axvline(0, color='k', lw=0.5)

    fig.tight_layout()
    out = os.path.join(cwd, 'verify_uniaxial.png')
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f'  Saved: {out}')

def plot_wall(cwd):
    print('\n=== Plotting wall results ===')
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    peak_data = {}
    for mt, color, ls in [('adina', '#D62728', '-'), ('ref', '#1F77B4', '--')]:
        mpath = os.path.join(cwd, f'vwall_{mt}_manual.txt')
        md = load_txt(mpath)
        if md is None:
            continue
        disp = md[:, 1]
        reactions = np.sum(md[:, 2:], axis=1)
        force = -reactions / 1000.0  # kN
        label = 'B&R' if mt == 'adina' else 'Reference'

        ax = axes[0]
        ax.plot(disp * 1000, force, color=color, ls=ls, lw=1.0, label=label, alpha=0.8)

        peak_data[mt] = (disp, force)

    setup_axes(axes[0], 'Displacement (mm)', 'Lateral Force (kN)', 'Shear Wall Hysteresis')
    axes[0].legend(fontsize=9)

    ax = axes[1]
    for mt, color, ls in [('adina', '#D62728', '-'), ('ref', '#1F77B4', '--')]:
        if mt not in peak_data:
            continue
        disp, force = peak_data[mt]
        pos_env_d, pos_env_f = [], []
        neg_env_d, neg_env_f = [], []
        for i in range(1, len(force)-1):
            if force[i] >= force[i-1] and force[i] >= force[i+1] and force[i] > 0:
                pos_env_d.append(disp[i]*1000); pos_env_f.append(force[i])
            if force[i] <= force[i-1] and force[i] <= force[i+1] and force[i] < 0:
                neg_env_d.append(disp[i]*1000); neg_env_f.append(force[i])
        label = 'B&R' if mt == 'adina' else 'Reference'
        if pos_env_d:
            ax.plot(pos_env_d, pos_env_f, 'o-', color=color, ms=3, lw=1.2, label=f'{label} +')
        if neg_env_d:
            ax.plot(neg_env_d, neg_env_f, 's--', color=color, ms=3, lw=1.2, label=f'{label} -')

    setup_axes(ax, 'Displacement (mm)', 'Force (kN)', 'Envelope Comparison')
    ax.legend(fontsize=8)

    fig.tight_layout()
    out = os.path.join(cwd, 'verify_wall.png')
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f'  Saved: {out}')

    for mt in peak_data:
        d, f = peak_data[mt]
        print(f'  {mt:6s}: F=[{np.min(f):.1f}, {np.max(f):.1f}] kN  d=[{np.min(d*1000):.1f}, {np.max(d*1000):.1f}] mm')

def plot_beam(cwd):
    print('\n=== Plotting beam results ===')
    fig, ax = plt.subplots(1, 1, figsize=(8, 5.5))

    for mt, color, ls in [('adina', '#D62728', '-'), ('ref', '#1F77B4', '--')]:
        mpath = os.path.join(cwd, f'vbeam_{mt}_manual.txt')
        md = load_txt(mpath)
        if md is None:
            continue
        tip_disp = md[:, 1]
        total_react = md[:, 2]
        force = -total_react / 1000.0  # kN
        label = 'B&R' if mt == 'adina' else 'Reference'
        ax.plot(tip_disp * 1000, force, color=color, ls=ls, lw=1.2, label=label)

    setup_axes(ax, 'Tip Displacement (mm)', 'Applied Force (kN)',
               'RC Cantilever Beam: Force-Displacement')
    ax.legend(fontsize=10)
    ax.axhline(0, color='k', lw=0.5)

    fig.tight_layout()
    out = os.path.join(cwd, 'verify_beam.png')
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f'  Saved: {out}')

def plot_plate(cwd):
    print('\n=== Plotting plate results ===')
    fig, ax = plt.subplots(1, 1, figsize=(8, 5.5))

    for mt, color, ls in [('adina', '#D62728', '-'), ('ref', '#1F77B4', '--')]:
        dpath = os.path.join(cwd, f'vplate_{mt}_disp.txt')
        ppath = os.path.join(cwd, f'vplate_{mt}_protocol.txt')
        dd = load_txt(dpath)
        if dd is None:
            continue
        center_w = dd[:, 1]
        pp = load_txt(ppath)
        if pp is not None:
            load_raw = pp.flatten() if pp.ndim > 1 else pp
            load_vals = load_raw[1:len(center_w)+1]  # skip initial 0, align with steps
        else:
            load_vals = dd[:, 0]
        n = min(len(center_w), len(load_vals))
        label = 'B&R' if mt == 'adina' else 'Reference'
        ax.plot(center_w[:n] * 1000, load_vals[:n] / 1000.0, color=color, ls=ls, lw=1.2, label=label)

    setup_axes(ax, 'Center Deflection (mm)', 'Applied Load (kN)', 'RC Plate: Load-Deflection')
    ax.legend(fontsize=10)
    ax.axhline(0, color='k', lw=0.5)

    fig.tight_layout()
    out = os.path.join(cwd, 'verify_plate.png')
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f'  Saved: {out}')

def print_summary(cwd):
    print('\n' + '='*60)
    print('  VERIFICATION SUMMARY')
    print('='*60)
    tests = [
        ('uniaxial', ['vuniax_{}_comp_disp.txt','vuniax_{}_tens_disp.txt','vuniax_{}_cyc_disp.txt']),
        ('wall',     ['vwall_{}_manual.txt']),
        ('beam',     ['vbeam_{}_disp.txt']),
        ('plate',    ['vplate_{}_disp.txt']),
    ]
    for tname, patterns in tests:
        for mt in ['adina', 'ref']:
            ok = all(os.path.exists(os.path.join(cwd, p.format(mt))) for p in patterns)
            status = 'OK' if ok else 'MISSING'
            nlines = 0
            for p in patterns:
                fp = os.path.join(cwd, p.format(mt))
                if os.path.exists(fp):
                    with open(fp) as f:
                        nlines += sum(1 for _ in f)
            print(f'  {tname:10s} [{mt:5s}]: {status:7s}  ({nlines} data points)')
    print('='*60)

# ─── Main ───
def main():
    args = parse_args()

    adina_exe = args.adina
    ref_exe   = args.ref

    if adina_exe is None:
        adina_exe = find_exe('OpenSees-2017.exe')
        if adina_exe is None:
            adina_exe = find_exe('OpenSees-*.exe')
    if ref_exe is None:
        try:
            candidates = []
            for entry in os.scandir(ROOT_DIR):
                if entry.is_file() and entry.name.endswith('.exe'):
                    if any(ord(c) > 127 for c in entry.name) and 'OpenSees' in entry.name:
                        candidates.append(entry.path)
            for path in candidates:
                if '参考' in path or '参考' in os.path.basename(path):
                    ref_exe = path
                    break
            if ref_exe is None and candidates:
                ref_exe = candidates[0]
        except Exception:
            pass
        if ref_exe is None:
            ref_exe = find_exe('OpenSees-ref*.exe')

    if adina_exe:
        adina_exe = os.path.abspath(adina_exe)
    if ref_exe:
        ref_exe = os.path.abspath(ref_exe)

    print('='*60)
    print('  ADINA Concrete Model Verification Suite')
    print('='*60)
    print(f'  ADINA exe: {adina_exe}')
    print(f'  Ref   exe: {ref_exe}')
    print(f'  Tests:     {args.tests}')
    print(f'  Work dir:  {SCRIPT_DIR}')
    print()

    test_map = {
        'uniaxial': [('verify_uniaxial.tcl', 'comp'),
                     ('verify_uniaxial.tcl', 'tens'),
                     ('verify_uniaxial.tcl', 'cyc')],
        'wall':     [('verify_wall.tcl', None)],
        'beam':     [('verify_beam.tcl', None)],
        'plate':    [('verify_plate.tcl', None)],
    }

    for tname in args.tests:
        entries = test_map.get(tname)
        if entries is None:
            print(f'Unknown test: {tname}')
            continue
        print(f'\n{"="*50}')
        print(f'  {tname.upper()}')
        print(f'{"="*50}')

        for tcl_file, sub_test in entries:
            tcl_path = os.path.join(SCRIPT_DIR, tcl_file)
            sub_label = f' ({sub_test})' if sub_test else ''
            print(f'\n--- {tname}{sub_label} ---')

            for mt, exe in [('adina', adina_exe), ('ref', ref_exe)]:
                if exe and os.path.isfile(exe):
                    run_opensees(exe, tcl_path, mt, SCRIPT_DIR, sub_test=sub_test)
                else:
                    print(f'  SKIP {mt} (exe not found)')

    print('\n--- Generating plots ---')
    if 'uniaxial' in args.tests: plot_uniaxial(SCRIPT_DIR)
    if 'wall' in args.tests:     plot_wall(SCRIPT_DIR)
    if 'beam' in args.tests:     plot_beam(SCRIPT_DIR)
    if 'plate' in args.tests:    plot_plate(SCRIPT_DIR)

    print_summary(SCRIPT_DIR)

if __name__ == '__main__':
    main()
