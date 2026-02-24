"""Run wall analysis with OpenSees-1706.exe and analyze results."""
import subprocess, os, sys, time, re
import numpy as np
sys.stdout.reconfigure(line_buffering=True)

WORKDIR = r'E:\Basic\Concrete_Model\ADINA'
EXE = os.path.join(WORKDIR, 'OpenSees-1706.exe')
TCL = os.path.join(WORKDIR, 'csp3.tcl')

# Clean old output
for f in ['1.txt', '53.txt']:
    p = os.path.join(WORKDIR, f)
    if os.path.exists(p):
        os.remove(p)

# Run wall analysis (10 min timeout)
print('='*60)
print('  WALL ANALYSIS with OpenSees-1706.exe')
print('='*60)
print('  Running...', flush=True)
t0 = time.time()
try:
    result = subprocess.run([EXE, TCL], cwd=WORKDIR, capture_output=True,
                           text=True, timeout=600)
    elapsed = time.time() - t0
    code = result.returncode
    ok = (code == 0 or code == -1073740940 or code == 3221226356)
    print(f'  Exit: {"OK" if ok else "FAIL"} (code={code}, {elapsed:.1f}s)')
except subprocess.TimeoutExpired:
    elapsed = time.time() - t0
    print(f'  TIMEOUT ({elapsed:.1f}s)')

# Parse results
disp_data, react_data = [], []

if os.path.exists(os.path.join(WORKDIR, '53.txt')):
    with open(os.path.join(WORKDIR, '53.txt')) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    disp_data.append([float(parts[0]), float(parts[1])*1000])
                except: pass

if os.path.exists(os.path.join(WORKDIR, '1.txt')):
    with open(os.path.join(WORKDIR, '1.txt')) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 6:
                try:
                    t = float(parts[0])
                    shear = sum(float(x) for x in parts[1:6]) / 1000
                    react_data.append([t, shear])
                except: pass

disp = np.array(disp_data) if disp_data else np.zeros((0,2))
react = np.array(react_data) if react_data else np.zeros((0,2))

print(f'\n  Displacement: {len(disp)} points')
print(f'  Reaction:     {len(react)} points')

if len(disp) > 0:
    print(f'  Disp range: [{disp[:,1].min():.2f}, {disp[:,1].max():.2f}] mm')
    print(f'  Time range: [{disp[:,0].min():.3f}, {disp[:,0].max():.3f}]')

if len(react) > 0:
    good = np.abs(react[:,1]) < 500
    n_good = np.sum(good)
    print(f'  Good force points: {n_good} / {len(react)}')

    if n_good > 0:
        rg = react[good]
        print(f'  Force range (good): [{rg[:,1].min():.1f}, {rg[:,1].max():.1f}] kN')

    # Find divergence
    for i in range(len(react)-1):
        if abs(react[i+1,1]) > 500 and abs(react[i,1]) < 500:
            print(f'\n  Divergence at step {i} (t={react[i,0]:.5f}):')
            s = max(0, i-3)
            for j in range(s, min(i+5, len(react))):
                d_at_t = ''
                for dd in disp_data:
                    if abs(dd[0] - react[j,0]) < 0.0001:
                        d_at_t = f' d={dd[1]:.2f}mm'
                        break
                print(f'    t={react[j,0]:.5f} F={react[j,1]:+10.1f} kN{d_at_t}')
            break

# Build force-displacement
fd_pairs = []
if len(disp) > 0 and len(react) > 0:
    di, ri = 0, 0
    while di < len(disp) and ri < len(react):
        if abs(disp[di,0] - react[ri,0]) < 0.001:
            if abs(react[ri,1]) < 500:
                fd_pairs.append([disp[di,1], react[ri,1]])
            di += 1; ri += 1
        elif disp[di,0] < react[ri,0]:
            di += 1
        else:
            ri += 1

if fd_pairs:
    fd = np.array(fd_pairs)
    print(f'\n  Force-disp pairs: {len(fd)}')
    print(f'  Disp range: [{fd[:,0].min():.2f}, {fd[:,0].max():.2f}] mm')
    print(f'  Force range: [{fd[:,1].min():.1f}, {fd[:,1].max():.1f}] kN')

    # Peaks
    peaks = []
    for i in range(1, len(fd)-1):
        if (fd[i,0] > fd[i-1,0] and fd[i,0] > fd[i+1,0]) or \
           (fd[i,0] < fd[i-1,0] and fd[i,0] < fd[i+1,0]):
            peaks.append(i)
    if peaks:
        print('\n  Displacement peaks:')
        for p in peaks[:15]:
            print(f'    d={fd[p,0]:+8.2f}mm  F={fd[p,1]:+8.1f}kN')

    # Load reference
    ref_rows = []
    with open(os.path.join(WORKDIR, '参考模型.txt')) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    rd, rf = float(parts[0]), float(parts[1])
                    if abs(rd) > 1e-6 or abs(rf) > 100:
                        ref_rows.append([rd*1000, rf/1000])
                except: pass
    ref = np.array(ref_rows) if ref_rows else None

    # Plot
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    ax = axes[0]
    if ref is not None:
        ax.plot(ref[:,0], ref[:,1], 'b-', lw=0.8, label='Reference', alpha=0.5)
    ax.plot(fd[:,0], fd[:,1], 'r-', lw=0.8, label='OpenSees-1706')
    ax.set_xlabel('Displacement (mm)'); ax.set_ylabel('Base Shear (kN)')
    ax.set_title('Hysteresis: 1706 vs Reference')
    ax.legend(); ax.grid(True, alpha=0.3)

    ax = axes[1]
    rg = react[np.abs(react[:,1]) < 500]
    ax.plot(rg[:,0], rg[:,1], 'r-', lw=0.5, label='1706')
    # Load 1045 for comparison
    try:
        r1045 = []
        with open(os.path.join(WORKDIR, 'wall_OpenSees_1045_react.txt')) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 6:
                    try:
                        t = float(parts[0])
                        fs = sum(float(x) for x in parts[1:6]) / 1000
                        if abs(fs) < 500:
                            r1045.append([t, fs])
                    except: pass
        if r1045:
            r1045 = np.array(r1045)
            ax.plot(r1045[:,0], r1045[:,1], 'g-', lw=0.5, label='1045', alpha=0.5)
    except: pass
    ax.set_xlabel('Time'); ax.set_ylabel('Base Shear (kN)')
    ax.set_title('Force vs Time Comparison')
    ax.legend(); ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(WORKDIR, 'wall_1706_comparison.png'), dpi=150)
    print('\n  Saved wall_1706_comparison.png')

# Version comparison
print('\n' + '='*60)
print('  VERSION COMPARISON (good force points)')
print('='*60)
for ver in ['1045', '1205', '1535', '1640', '1706']:
    if ver == '1706':
        cnt = n_good if len(react) > 0 else 0
    else:
        rfile = os.path.join(WORKDIR, f'wall_OpenSees_{ver}_react.txt')
        cnt = 0
        try:
            with open(rfile) as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 6:
                        try:
                            fs = sum(float(x) for x in parts[1:6]) / 1000
                            if abs(fs) < 500: cnt += 1
                        except: pass
        except: pass
    print(f'  OpenSees-{ver}: {cnt:5d} points')
print(f'  Reference:      {len(ref) if ref is not None else 0:5d} points')
