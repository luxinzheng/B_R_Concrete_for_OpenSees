"""
Detailed wall hysteresis plot with reasonable axis limits.
Filter out force spikes from sub-stepping and focus on the quality of results.
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
    ('OpenSees_0929', 'OpenSees-0929', '#1f77b4', '--', 1.2),
    ('OpenSees_0943', 'OpenSees-0943', '#ff7f0e', '--', 1.2),
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
    d = safe_load(disp_file, 2)
    r = safe_load(react_file, 6)
    if d is None or r is None:
        continue
    n = min(len(d), len(r))
    t_d = d[:n, 0]
    lat_start = 0
    for i in range(1, len(t_d)):
        if t_d[i] < t_d[i-1]:
            lat_start = i
            break
    if lat_start == 0:
        lat_start = 10

    t_lat = d[lat_start:n, 0]
    disp_mm = d[lat_start:n, 1] * 1000
    force_kN = r[lat_start:n, 1:6].sum(axis=1) / 1000

    data[tag] = {
        'time': t_lat, 'disp': disp_mm, 'force': force_kN,
        'label': label, 'color': color, 'ls': ls, 'lw': lw,
    }

# Determine reasonable force range from 1045 data (percentile clip)
if 'OpenSees_1045' in data:
    f = data['OpenSees_1045']['force']
    f_valid = f[np.isfinite(f)]
    # Use robust percentile to find main data range
    p5, p95 = np.percentile(f_valid, [2, 98])
    f_range = max(abs(p5), abs(p95))
    f_lim = f_range * 1.5
    print(f"1045 force P2/P98: [{p5:.0f}, {p95:.0f}] kN, plot limit: +/-{f_lim:.0f} kN")
else:
    f_lim = 500

fig, axes = plt.subplots(2, 3, figsize=(20, 12))

# --- Row 1: Overview ---

# Plot 1: Displacement time history with protocol
ax = axes[0, 0]
ax.plot(prot_time[:40], prot_disp_mm[:40], 'k-', linewidth=0.8, alpha=0.4, label='Target')
for tag in ['OpenSees_0929', 'OpenSees_0943', 'OpenSees_1045']:
    if tag in data:
        d = data[tag]
        ax.plot(d['time'], d['disp'], d['ls'], color=d['color'],
                label=d['label'], linewidth=d['lw'])
ax.set_xlabel('Pseudo-time')
ax.set_ylabel('Top Displacement (mm)')
ax.set_title('Displacement Tracking', fontsize=13, fontweight='bold')
ax.legend(fontsize=9, loc='upper left')
ax.grid(True, alpha=0.3)

# Plot 2: Hysteresis all versions (clipped force)
ax = axes[0, 1]
for tag in ['OpenSees_0929', 'OpenSees_0943', 'OpenSees_1045']:
    if tag in data:
        d = data[tag]
        force_clip = np.clip(d['force'], -f_lim, f_lim)
        ax.plot(d['disp'], force_clip, d['ls'], color=d['color'],
                label=d['label'], linewidth=d['lw'], alpha=0.8)
ax.set_xlabel('Top Displacement (mm)')
ax.set_ylabel('Base Shear (kN)')
ax.set_title(f'Hysteresis (clipped +/-{f_lim:.0f} kN)', fontsize=13, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.axhline(y=0, color='k', linewidth=0.5)
ax.axvline(x=0, color='k', linewidth=0.5)

# Plot 3: Base shear time history (clipped)
ax = axes[0, 2]
for tag in ['OpenSees_0929', 'OpenSees_0943', 'OpenSees_1045']:
    if tag in data:
        d = data[tag]
        force_clip = np.clip(d['force'], -f_lim, f_lim)
        ax.plot(d['time'], force_clip, d['ls'], color=d['color'],
                label=d['label'], linewidth=d['lw'])
ax.set_xlabel('Pseudo-time')
ax.set_ylabel('Base Shear (kN)')
ax.set_title('Base Shear History (clipped)', fontsize=13, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# --- Row 2: 1045 detailed analysis ---

if 'OpenSees_1045' in data:
    d = data['OpenSees_1045']

    # Plot 4: 1045 full hysteresis
    ax = axes[1, 0]
    force_clip = np.clip(d['force'], -f_lim, f_lim)
    ax.plot(d['disp'], force_clip, '-', color='#2ca02c', linewidth=1.2)
    ax.plot(d['disp'][0], force_clip[0], 'go', ms=8, label='Start')
    ax.plot(d['disp'][-1], force_clip[-1], 'rs', ms=8,
            label=f'End (t={d["time"][-1]:.2f})')
    ax.set_xlabel('Top Displacement (mm)')
    ax.set_ylabel('Base Shear (kN)')
    ax.set_title('1045 Full Hysteresis', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='k', linewidth=0.5)
    ax.axvline(x=0, color='k', linewidth=0.5)

    # Plot 5: First cycle only (t < 1.6)
    ax = axes[1, 1]
    mask1 = d['time'] < 1.6
    ax.plot(d['disp'][mask1], force_clip[mask1], '-', color='#2ca02c', linewidth=1.5)
    ax.set_xlabel('Top Displacement (mm)')
    ax.set_ylabel('Base Shear (kN)')
    ax.set_title('1045 First Cycle (t<1.6)', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='k', linewidth=0.5)
    ax.axvline(x=0, color='k', linewidth=0.5)

    # Plot 6: Force stability analysis - sliding window stdev
    ax = axes[1, 2]
    window = 50
    force_arr = d['force']
    time_arr = d['time']
    n = len(force_arr)
    stdev = np.array([np.std(force_arr[max(0,i-window):i+1]) for i in range(n)])
    ax.semilogy(time_arr, stdev + 1, '-', color='#d62728', linewidth=1.0)
    ax.set_xlabel('Pseudo-time')
    ax.set_ylabel('Force Std (kN, sliding window=50)')
    ax.set_title('1045 Force Stability', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)

    # Find where instability starts
    threshold = np.median(stdev[stdev > 0]) * 10
    unstable = np.where(stdev > threshold)[0]
    if len(unstable) > 0:
        t_unstable = time_arr[unstable[0]]
        ax.axvline(x=t_unstable, color='r', linestyle='--', alpha=0.7,
                   label=f'Instability onset: t={t_unstable:.2f}')
        ax.legend(fontsize=9)

plt.suptitle('csp3.tcl Shear Wall: Detailed 1045 Analysis\n'
             '(Section 6A volumetric strain fix)',
             fontsize=16, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(f'{WORK_DIR}/wall_1045_detail.png', dpi=150, bbox_inches='tight')
print(f"\nSaved: wall_1045_detail.png")

# Print key force values at specific time points for 1045
if 'OpenSees_1045' in data:
    d = data['OpenSees_1045']
    print("\n--- 1045 Key Points ---")
    targets = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
    for t_target in targets:
        idx = np.argmin(np.abs(d['time'] - t_target))
        print(f"  t={d['time'][idx]:.3f}: disp={d['disp'][idx]:.2f} mm, "
              f"shear={d['force'][idx]:.1f} kN")
