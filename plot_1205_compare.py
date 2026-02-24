"""
Compare OpenSees-1205 (principal strain fix) vs 1045 vs reference model.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
sys.stdout.reconfigure(line_buffering=True)

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

def load_wall(tag):
    d = safe_load(f'{WORK_DIR}/wall_{tag}_disp.txt', 2)
    r = safe_load(f'{WORK_DIR}/wall_{tag}_react.txt', 6)
    if d is None or r is None:
        return None, None, None
    n = min(len(d), len(r))
    t = d[:n, 0]
    lat = 0
    for i in range(1, len(t)):
        if t[i] < t[i-1]:
            lat = i
            break
    if lat == 0: lat = 10
    return (d[lat:n, 0], d[lat:n, 1]*1000, r[lat:n, 1:6].sum(axis=1)/1000)

# Load reference data
ref_rows = []
with open(f'{WORK_DIR}/参考模型.txt', 'r') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2:
            try:
                ref_rows.append([float(parts[0]), float(parts[1])])
            except ValueError:
                pass
ref = np.array(ref_rows)
ref_d = ref[10:, 0] * 1000
ref_f = ref[10:, 1] / 1000

# Load wall data
t_1045, d_1045, f_1045 = load_wall('OpenSees_1045')
t_1205, d_1205, f_1205 = load_wall('OpenSees_1205')

# Load protocol
prot = np.loadtxt(f'{WORK_DIR}/shuju2.txt')
prot_t = np.arange(1, len(prot)+1) * 0.1
prot_d = prot * 1000

print(f"Reference: {len(ref_d)} pts, disp [{ref_d.min():.0f}, {ref_d.max():.0f}] mm")
if t_1045 is not None:
    print(f"1045: {len(t_1045)} pts, t=[{t_1045[-1]:.3f}], disp [{d_1045.min():.1f}, {d_1045.max():.1f}] mm")
if t_1205 is not None:
    print(f"1205: {len(t_1205)} pts, t=[{t_1205[-1]:.3f}], disp [{d_1205.min():.1f}, {d_1205.max():.1f}] mm")

fig, axes = plt.subplots(2, 3, figsize=(21, 12))

# --- Plot 1: Reference full hysteresis ---
ax = axes[0, 0]
ax.plot(ref_d, ref_f, 'b-', linewidth=1.5, label='Reference')
ax.set_xlabel('Top Displacement (mm)')
ax.set_ylabel('Base Shear (kN)')
ax.set_title('Reference Model (352 pts)', fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3)
ax.axhline(y=0, color='k', linewidth=0.5)
ax.axvline(x=0, color='k', linewidth=0.5)

# --- Plot 2: 1205 hysteresis (clipped) ---
ax = axes[0, 1]
if t_1205 is not None:
    f_clip = np.clip(f_1205, -300, 300)
    ax.plot(d_1205, f_clip, 'g-', linewidth=1.5, label='1205 (principal strain)')
    ax.set_xlabel('Top Displacement (mm)')
    ax.set_ylabel('Base Shear (kN, clipped ±300)')
    ax.set_title(f'OpenSees-1205 ({len(t_1205)} pts, t→{t_1205[-1]:.2f})', fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='k', linewidth=0.5)
    ax.axvline(x=0, color='k', linewidth=0.5)

# --- Plot 3: First cycle overlay ---
ax = axes[0, 2]
# Reference first cycle (first 15 pts)
ax.plot(ref_d[:15], ref_f[:15], 'b-o', linewidth=2, markersize=4, label='Reference')
# 1205 first cycle
if t_1205 is not None:
    mask = t_1205 < 1.6
    f_c = np.clip(f_1205[mask], -300, 300)
    ax.plot(d_1205[mask], f_c, 'g-', linewidth=1.5, label='1205')
# 1045 first cycle
if t_1045 is not None:
    mask = t_1045 < 1.6
    f_c = np.clip(f_1045[mask], -300, 300)
    ax.plot(d_1045[mask], f_c, 'r--', linewidth=1.0, alpha=0.6, label='1045')
ax.set_xlabel('Top Displacement (mm)')
ax.set_ylabel('Base Shear (kN)')
ax.set_title('First Cycle Overlay', fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3)
ax.axhline(y=0, color='k', linewidth=0.5)
ax.axvline(x=0, color='k', linewidth=0.5)

# --- Plot 4: Displacement tracking ---
ax = axes[1, 0]
ax.plot(prot_t[:50], prot_d[:50], 'k-', linewidth=0.8, alpha=0.3, label='Target')
if t_1045 is not None:
    ax.plot(t_1045, d_1045, 'r--', linewidth=1, alpha=0.6, label=f'1045 (t→{t_1045[-1]:.2f})')
if t_1205 is not None:
    ax.plot(t_1205, d_1205, 'g-', linewidth=2, label=f'1205 (t→{t_1205[-1]:.2f})')
ax.set_xlabel('Pseudo-time')
ax.set_ylabel('Top Displacement (mm)')
ax.set_title('Analysis Progress', fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_xlim(-0.1, 5.0)

# --- Plot 5: 1205 second cycle detail ---
ax = axes[1, 1]
if t_1205 is not None:
    mask2 = (t_1205 > 1.5) & (t_1205 < 3.0)
    f_c = np.clip(f_1205[mask2], -300, 300)
    ax.plot(d_1205[mask2], f_c, 'g-', linewidth=1.5, label='1205 (cycle 2)')
    # Reference second cycle (pts 15-53)
    ax.plot(ref_d[14:43], ref_f[14:43], 'b-o', linewidth=1.5, markersize=3,
            alpha=0.7, label='Reference (cycle 2)')
    ax.set_xlabel('Top Displacement (mm)')
    ax.set_ylabel('Base Shear (kN, clipped ±300)')
    ax.set_title('Second Cycle Comparison', fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='k', linewidth=0.5)
    ax.axvline(x=0, color='k', linewidth=0.5)

# --- Plot 6: Force at key points ---
ax = axes[1, 2]
if t_1205 is not None:
    # 1205 key points
    targets = [0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 3.5]
    times = []; disps = []; forces = []
    for tt in targets:
        if tt <= t_1205[-1]:
            idx = np.argmin(np.abs(t_1205 - tt))
            times.append(t_1205[idx])
            disps.append(d_1205[idx])
            forces.append(f_1205[idx])
    print("\n--- 1205 Key Points ---")
    for i in range(len(times)):
        print(f"  t={times[i]:.3f}: disp={disps[i]:.2f} mm, shear={forces[i]:.1f} kN")

    # Also show force stability
    window = 20
    stdev = np.array([np.std(f_1205[max(0,i-window):i+1])
                       for i in range(len(f_1205))])
    ax.semilogy(t_1205, stdev + 1, 'g-', linewidth=1.0, label='1205')
    if t_1045 is not None:
        stdev_1045 = np.array([np.std(f_1045[max(0,i-window):i+1])
                                for i in range(len(f_1045))])
        ax.semilogy(t_1045, stdev_1045 + 1, 'r--', linewidth=0.8, alpha=0.6, label='1045')
    ax.set_xlabel('Pseudo-time')
    ax.set_ylabel('Force Std (kN, window=20)')
    ax.set_title('Force Stability', fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)

plt.suptitle('csp3.tcl Shear Wall: 1205 (Principal Strain Fix) vs Reference',
             fontsize=16, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(f'{WORK_DIR}/wall_1205_comparison.png', dpi=150, bbox_inches='tight')
print(f"\nSaved: wall_1205_comparison.png")
