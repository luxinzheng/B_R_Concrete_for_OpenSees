"""
Run wall analysis with specified OpenSees executable and generate comparison plots.

Usage:
    python run_wall_analysis.py <exe_name>
    e.g. python run_wall_analysis.py OpenSees-2300.exe
"""
import subprocess, sys, os, shutil, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

WORKDIR = r'E:\Basic\Concrete_Model\ADINA'

def find_exe(name_hint=None):
    """Find OpenSees exe. If hint given, use it. Otherwise find latest numbered exe."""
    if name_hint and os.path.exists(os.path.join(WORKDIR, name_hint)):
        return name_hint
    exes = []
    for f in os.listdir(WORKDIR):
        if f.startswith('OpenSees-') and f.endswith('.exe'):
            try:
                num = int(f.replace('OpenSees-', '').replace('.exe', ''))
                exes.append((num, f))
            except ValueError:
                pass
    if exes:
        exes.sort(reverse=True)
        return exes[0][1]
    return None

def find_ref_exe():
    """Find reference model exe (Chinese name)."""
    for f in os.listdir(WORKDIR):
        if f.endswith('.exe') and '参考' in f:
            return f
    return None

def load_manual_wall(manual_file, force_limit=500):
    """Load from manual disp_manual.txt: time disp r1 r2 r3 r4 r5 [SKIP]"""
    if not os.path.exists(manual_file):
        return None
    fd_pairs = []
    with open(manual_file) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 7:
                continue
            if 'SKIP' in parts:
                continue  # Skip force-through data
            try:
                disp_m = float(parts[1])
                shear_n = sum(float(parts[i]) for i in range(2, 7))
                disp_mm = disp_m * 1000
                shear_kn = shear_n / 1000
                if abs(shear_kn) < force_limit:
                    fd_pairs.append([disp_mm, shear_kn])
            except (ValueError, IndexError):
                pass
    return np.array(fd_pairs) if fd_pairs else None

def load_opensees_wall(react_file, disp_file, force_limit=500):
    """Load and match reaction/displacement data."""
    react_data, disp_data = [], []
    if os.path.exists(react_file):
        with open(react_file) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 6:
                    try:
                        t = float(parts[0])
                        shear = sum(float(x) for x in parts[1:6]) / 1000
                        react_data.append([t, shear])
                    except: pass
    if os.path.exists(disp_file):
        with open(disp_file) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        disp_data.append([float(parts[0]), float(parts[1]) * 1000])
                    except: pass

    react = np.array(react_data) if react_data else np.zeros((0, 2))
    disp = np.array(disp_data) if disp_data else np.zeros((0, 2))

    fd_pairs = []
    if len(disp) > 0 and len(react) > 0:
        di, ri = 0, 0
        while di < len(disp) and ri < len(react):
            if abs(disp[di, 0] - react[ri, 0]) < 0.001:
                if abs(react[ri, 1]) < force_limit:
                    fd_pairs.append([disp[di, 1], react[ri, 1]])
                di += 1; ri += 1
            elif disp[di, 0] < react[ri, 0]:
                di += 1
            else:
                ri += 1
    return np.array(fd_pairs) if fd_pairs else None

def extract_envelope(fd):
    """Extract peak force at each displacement amplitude."""
    if fd is None or len(fd) == 0:
        return None, None
    d, f = fd[:, 0], fd[:, 1]
    pos_env, neg_env = [], []
    for i in range(1, len(d) - 1):
        if d[i] > d[i-1] and d[i] > d[i+1] and d[i] > 0.5:
            pos_env.append([d[i], f[i]])
        elif d[i] < d[i-1] and d[i] < d[i+1] and d[i] < -0.5:
            neg_env.append([d[i], f[i]])
    return (np.array(pos_env) if pos_env else None,
            np.array(neg_env) if neg_env else None)

def run_analysis(exe_name, tag):
    """Run csp3.tcl with the given executable."""
    exe_path = os.path.join(WORKDIR, exe_name)
    print(f'\n{"="*60}')
    print(f'  Running wall analysis with {exe_name}')
    print(f'{"="*60}')

    # Clean old outputs
    for f in ['1.txt', '53.txt', 'disp_manual.txt']:
        fp = os.path.join(WORKDIR, f)
        if os.path.exists(fp):
            os.remove(fp)

    t0 = time.time()
    try:
        result = subprocess.run(
            [exe_path, 'csp3.tcl'],
            cwd=WORKDIR, capture_output=True, text=True, timeout=1200)
        elapsed = time.time() - t0
        print(f'  Completed in {elapsed:.1f}s (exit code {result.returncode})')
        if result.stdout:
            for line in result.stdout.strip().split('\n')[-10:]:
                print(f'    {line}')
        if result.returncode != 0 and result.stderr:
            for line in result.stderr.strip().split('\n')[-5:]:
                print(f'    [stderr] {line}')
    except subprocess.TimeoutExpired:
        elapsed = time.time() - t0
        print(f'  TIMEOUT ({elapsed:.0f}s)')
        result = type('', (), {'returncode': -1, 'stdout': '', 'stderr': 'TIMEOUT'})()

    # Save output files with tag
    for src, suffix in [('1.txt', 'react'), ('53.txt', 'disp'), ('disp_manual.txt', 'manual')]:
        src_path = os.path.join(WORKDIR, src)
        dst_path = os.path.join(WORKDIR, f'wall_{tag}_{suffix}.txt')
        if os.path.exists(src_path):
            shutil.copy2(src_path, dst_path)
            lines = sum(1 for _ in open(dst_path))
            print(f'  {src} -> wall_{tag}_{suffix}.txt ({lines} lines)')

def main():
    # Determine exe name
    exe_hint = sys.argv[1] if len(sys.argv) > 1 else None
    exe_name = find_exe(exe_hint)
    if not exe_name:
        print('ERROR: No OpenSees executable found!')
        sys.exit(1)

    # Extract version tag
    tag = exe_name.replace('OpenSees-', '').replace('.exe', '')
    print(f'Using executable: {exe_name} (tag: {tag})')

    # Run analysis
    run_analysis(exe_name, tag)

    # Load results - try manual file first (more reliable), then standard recorders
    manual_file = os.path.join(WORKDIR, f'wall_{tag}_manual.txt')
    fd_new = load_manual_wall(manual_file)
    if fd_new is None:
        fd_new = load_opensees_wall(
            os.path.join(WORKDIR, f'wall_{tag}_react.txt'),
            os.path.join(WORKDIR, f'wall_{tag}_disp.txt'))

    # Load reference
    ref_rows = []
    try:
        with open(os.path.join(WORKDIR, '参考模型.txt')) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        ref_rows.append([float(parts[0]) * 1000, float(parts[1]) / 1000])
                    except: pass
    except: pass
    ref = np.array(ref_rows) if ref_rows else None

    # Load previous best (2145)
    fd_prev = load_opensees_wall(
        os.path.join(WORKDIR, 'wall_OpenSees_2145_react.txt'),
        os.path.join(WORKDIR, 'wall_OpenSees_2145_disp.txt'))

    print(f'\nResults:')
    print(f'  New ({tag}):    {len(fd_new) if fd_new is not None else 0} points')
    print(f'  Reference:     {len(ref) if ref is not None else 0} points')
    print(f'  Previous(2145):{len(fd_prev) if fd_prev is not None else 0} points')

    # ============================================================
    # PLOT 1: Hysteresis comparison
    # ============================================================
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    ax = axes[0]
    ax.set_title(f'Hysteresis: ADINA-{tag} vs Reference', fontsize=13, fontweight='bold')
    if ref is not None:
        ax.plot(ref[:, 0], ref[:, 1], 'b-', lw=0.7, alpha=0.6, label='Reference')
    if fd_new is not None:
        ax.plot(fd_new[:, 0], fd_new[:, 1], 'r-', lw=0.9, label=f'ADINA-{tag}')
    ax.set_xlabel('Displacement (mm)', fontsize=11)
    ax.set_ylabel('Base Shear (kN)', fontsize=11)
    ax.legend(fontsize=10, loc='upper left')
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='k', lw=0.5)
    ax.axvline(0, color='k', lw=0.5)

    ax = axes[1]
    ax.set_title(f'Detail (ADINA-{tag} range)', fontsize=13, fontweight='bold')
    if ref is not None:
        ax.plot(ref[:, 0], ref[:, 1], 'b-', lw=0.7, alpha=0.6, label='Reference')
    if fd_new is not None:
        ax.plot(fd_new[:, 0], fd_new[:, 1], 'r-', lw=0.9, label=f'ADINA-{tag}')
        xlim = max(abs(fd_new[:, 0].min()), abs(fd_new[:, 0].max())) * 1.15
        ax.set_xlim([-xlim, xlim])
    ax.set_xlabel('Displacement (mm)', fontsize=11)
    ax.set_ylabel('Base Shear (kN)', fontsize=11)
    ax.legend(fontsize=10, loc='upper left')
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='k', lw=0.5)
    ax.axvline(0, color='k', lw=0.5)

    plt.tight_layout()
    fname = os.path.join(WORKDIR, f'wall_hysteresis_{tag}.png')
    plt.savefig(fname, dpi=150)
    print(f'Saved: {fname}')

    # ============================================================
    # PLOT 2: Envelope comparison
    # ============================================================
    fig2, axes2 = plt.subplots(1, 2, figsize=(16, 7))

    ax = axes2[0]
    ax.set_title('Envelope Comparison', fontsize=13, fontweight='bold')
    for label, fd, color, marker in [
        ('Reference', ref, 'b', 'o'),
        (f'ADINA-{tag}', fd_new, 'r', 's'),
        ('ADINA-2145', fd_prev, 'gray', '^')]:
        if fd is not None:
            ax.plot(fd[:, 0], fd[:, 1], color=color, lw=0.3, alpha=0.15)
            pos, neg = extract_envelope(fd)
            if pos is not None:
                ax.plot(pos[:, 0], pos[:, 1], color=color, marker=marker, ms=5,
                        lw=1.5, label=f'{label} (+)')
            if neg is not None:
                ax.plot(neg[:, 0], neg[:, 1], color=color, marker=marker, ms=5,
                        lw=1.5, ls='--', label=f'{label} (-)')
    ax.set_xlabel('Displacement (mm)', fontsize=11)
    ax.set_ylabel('Base Shear (kN)', fontsize=11)
    ax.legend(fontsize=8, loc='best')
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='k', lw=0.5)

    ax = axes2[1]
    ax.set_title('Displacement History', fontsize=13, fontweight='bold')
    if ref is not None:
        ax.plot(np.arange(len(ref)), ref[:, 0], 'b-', lw=0.5, alpha=0.6, label='Reference')
    if fd_new is not None:
        ax.plot(np.arange(len(fd_new)), fd_new[:, 0], 'r-', lw=0.8, label=f'ADINA-{tag}')
    if fd_prev is not None:
        ax.plot(np.arange(len(fd_prev)), fd_prev[:, 0], 'gray', lw=0.5, alpha=0.5, label='ADINA-2145')
    ax.set_xlabel('Step', fontsize=11)
    ax.set_ylabel('Displacement (mm)', fontsize=11)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fname2 = os.path.join(WORKDIR, f'wall_envelope_{tag}.png')
    plt.savefig(fname2, dpi=150)
    print(f'Saved: {fname2}')

    # ============================================================
    # Statistics
    # ============================================================
    print(f'\n{"="*60}')
    print(f'  ANALYSIS SUMMARY')
    print(f'{"="*60}')
    if fd_new is not None:
        d, f = fd_new[:, 0], fd_new[:, 1]
        print(f'  ADINA-{tag}:')
        print(f'    Points:      {len(fd_new)}')
        print(f'    Disp range:  [{d.min():.2f}, {d.max():.2f}] mm')
        print(f'    Force range: [{f.min():.1f}, {f.max():.1f}] kN')
        sign_changes = sum(1 for i in range(1, len(d)) if d[i]*d[i-1] < 0)
        print(f'    Approx cycles: {sign_changes // 2}')
        pos, neg = extract_envelope(fd_new)
        if pos is not None:
            print(f'    Positive peaks:')
            for p in pos:
                print(f'      d={p[0]:+7.2f} mm  F={p[1]:+8.1f} kN')
        if neg is not None:
            print(f'    Negative peaks:')
            for p in neg:
                print(f'      d={p[0]:+7.2f} mm  F={p[1]:+8.1f} kN')
    if ref is not None:
        d, f = ref[:, 0], ref[:, 1]
        print(f'\n  Reference:')
        print(f'    Points:      {len(ref)}')
        print(f'    Disp range:  [{d.min():.2f}, {d.max():.2f}] mm')
        print(f'    Force range: [{f.min():.1f}, {f.max():.1f}] kN')

if __name__ == '__main__':
    main()
