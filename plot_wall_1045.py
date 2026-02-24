"""
Plot csp3.tcl wall analysis comparison: 0929 vs 0943 vs 1045
Recorder format:
  1.txt:  time Rx1 Rx2 Rx3 Rx4 Rx5  (base shear = sum of Rx)
  53.txt: time Ux53                   (top lateral displacement)
Gravity: first 10 rows (time 0.1 to 1.0), then loadConst -time 0.0
Lateral: starts from time 0.01 after gravity reset
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

WORK_DIR = r'e:\Basic\Concrete_Model\ADINA'

def safe_load(filepath, ncols):
    rows = []
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= ncols:
                try:
                    vals = [float(parts[i]) for i in range(ncols)]
                    rows.append(vals)
                except ValueError:
                    pass
    return np.array(rows) if rows else None

versions = [
    ('OpenSees_0929', 'OpenSees-0929', '#1f77b4', '--', 1.0),
    ('OpenSees_0943', 'OpenSees-0943', '#ff7f0e', '--', 1.0),
    ('OpenSees_1045', 'OpenSees-1045 (fixed)', '#2ca02c', '-', 2.0),
]

# Load displacement protocol
disp_protocol = np.loadtxt(f'{WORK_DIR}/shuju2.txt')
prot_time = np.arange(1, len(disp_protocol)+1) * 0.1
prot_disp_mm = disp_protocol * 1000

data = {}
for tag, label, color, ls, lw in versions:
    disp_file = f'{WORK_DIR}/wall_{tag}_disp.txt'
    react_file = f'{WORK_DIR}/wall_{tag}_react.txt'
    try:
        d = safe_load(disp_file, 2)
        r = safe_load(react_file, 6)
        if d is None or r is None:
            print(f"{label}: NO DATA")
            continue

        # Separate gravity and lateral phases
        # Gravity: time goes 0.1, 0.2, ..., 1.0
        # After loadConst -time 0.0, time resets and goes 0.01, 0.02, ...
        # Detect: time[i] < time[i-1] indicates the reset point
        n = min(len(d), len(r))
        t_d = d[:n, 0]
        lateral_start = 0
        for i in range(1, len(t_d)):
            if t_d[i] < t_d[i-1]:
                lateral_start = i
                break

        if lateral_start == 0:
            print(f"{label}: Could not find gravity/lateral boundary")
            lateral_start = 10  # fallback

        t_lat = d[lateral_start:n, 0]
        disp_mm = d[lateral_start:n, 1] * 1000
        base_shear_kN = r[lateral_start:n, 1:6].sum(axis=1) / 1000

        data[tag] = {
            'time': t_lat,
            'disp': disp_mm,
            'force': base_shear_kN,
            'label': label,
            'color': color,
            'ls': ls,
            'lw': lw,
            'n_gravity': lateral_start,
            'n_lateral': len(t_lat),
        }
        print(f"{label}: {len(t_lat)} lateral pts (gravity={lateral_start}), "
              f"t=[{t_lat[0]:.3f}, {t_lat[-1]:.3f}], "
              f"disp=[{disp_mm.min():.2f}, {disp_mm.max():.2f}] mm, "
              f"shear=[{base_shear_kN.min():.1f}, {base_shear_kN.max():.1f}] kN")
    except Exception as e:
        print(f"{label}: ERROR: {e}")

fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# --- Plot 1: Hysteresis (Force vs Displacement) ---
ax = axes[0, 0]
for tag in ['OpenSees_0929', 'OpenSees_0943', 'OpenSees_1045']:
    if tag in data:
        d = data[tag]
        ax.plot(d['disp'], d['force'], d['ls'], color=d['color'],
                label=d['label'], linewidth=d['lw'], alpha=0.9)
ax.set_xlabel('Top Displacement (mm)', fontsize=12)
ax.set_ylabel('Base Shear (kN)', fontsize=12)
ax.set_title('Hysteresis Curve Comparison', fontsize=14, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
ax.axhline(y=0, color='k', linewidth=0.5)
ax.axvline(x=0, color='k', linewidth=0.5)

# --- Plot 2: Target vs actual displacement ---
ax = axes[0, 1]
# Plot target displacement protocol
ax.plot(prot_time[:30], prot_disp_mm[:30], 'k-', linewidth=0.8, alpha=0.5,
        label='Target Protocol')
for tag in ['OpenSees_0929', 'OpenSees_0943', 'OpenSees_1045']:
    if tag in data:
        d = data[tag]
        ax.plot(d['time'], d['disp'], d['ls'], color=d['color'],
                label=d['label'], linewidth=d['lw'])
ax.set_xlabel('Pseudo-time', fontsize=12)
ax.set_ylabel('Top Displacement (mm)', fontsize=12)
ax.set_title('Displacement vs Target Protocol', fontsize=14, fontweight='bold')
ax.legend(fontsize=10, loc='upper left')
ax.grid(True, alpha=0.3)

# --- Plot 3: Force time history ---
ax = axes[1, 0]
for tag in ['OpenSees_0929', 'OpenSees_0943', 'OpenSees_1045']:
    if tag in data:
        d = data[tag]
        ax.plot(d['time'], d['force'], d['ls'], color=d['color'],
                label=d['label'], linewidth=d['lw'])
ax.set_xlabel('Pseudo-time', fontsize=12)
ax.set_ylabel('Base Shear (kN)', fontsize=12)
ax.set_title('Base Shear Time History', fontsize=14, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# --- Plot 4: 1045 detail hysteresis ---
ax = axes[1, 1]
if 'OpenSees_1045' in data:
    d = data['OpenSees_1045']
    ax.plot(d['disp'], d['force'], '-', color='#2ca02c', linewidth=1.5)
    ax.plot(d['disp'][0], d['force'][0], 'go', ms=8, label='Start')
    ax.plot(d['disp'][-1], d['force'][-1], 'rs', ms=8,
            label=f'End (t={d["time"][-1]:.2f})')
    ax.set_xlabel('Top Displacement (mm)', fontsize=12)
    ax.set_ylabel('Base Shear (kN)', fontsize=12)
    ax.set_title('OpenSees-1045 Hysteresis Detail', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='k', linewidth=0.5)
    ax.axvline(x=0, color='k', linewidth=0.5)

plt.suptitle('csp3.tcl Shear Wall: Section 6A Fix Comparison\n'
             '(OpenSees-1045 with volumetric strain check)',
             fontsize=16, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(f'{WORK_DIR}/wall_1045_comparison.png', dpi=150, bbox_inches='tight')
print(f"\nSaved: wall_1045_comparison.png")

# Progress summary
print("\n" + "=" * 60)
print("ANALYSIS PROGRESS COMPARISON")
print("=" * 60)
print(f"\nDisplacement protocol from shuju2.txt:")
print(f"  First cycle: 0 -> +4mm -> 0 -> -4mm -> 0")
print(f"  Second cycle: 0 -> +7mm -> 0 -> -7mm -> 0")
for tag, label, *_ in versions:
    if tag in data:
        d = data[tag]
        t_max = d['time'][-1]
        n_pts = d['n_lateral']
        d_range = f"[{d['disp'].min():.1f}, {d['disp'].max():.1f}]"
        f_range = f"[{d['force'].min():.0f}, {d['force'].max():.0f}]"
        print(f"\n  {d['label']}:")
        print(f"    Lateral pts: {n_pts}, max time: {t_max:.3f}")
        print(f"    Disp range:  {d_range} mm")
        print(f"    Force range: {f_range} kN")
        # Estimate which part of protocol reached
        if t_max < 0.4:
            stage = "First positive push (incomplete)"
        elif t_max < 0.8:
            stage = "First positive return"
        elif t_max < 1.2:
            stage = "First negative push"
        elif t_max < 1.6:
            stage = "First negative return"
        elif t_max < 2.1:
            stage = "Second positive push (0->7mm)"
        elif t_max < 3.0:
            stage = "Second positive return + negative push"
        else:
            stage = f"Advanced (t={t_max:.1f})"
        print(f"    Stage:       {stage}")
