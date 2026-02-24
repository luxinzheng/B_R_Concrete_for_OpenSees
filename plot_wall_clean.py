"""
Clean hysteresis plot focusing on the first cycle quality.
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

def load_version(tag):
    d = safe_load(f'{WORK_DIR}/wall_{tag}_disp.txt', 2)
    r = safe_load(f'{WORK_DIR}/wall_{tag}_react.txt', 6)
    n = min(len(d), len(r))
    t = d[:n, 0]
    lat = 0
    for i in range(1, len(t)):
        if t[i] < t[i-1]:
            lat = i
            break
    if lat == 0: lat = 10
    return (d[lat:n, 0], d[lat:n, 1]*1000, r[lat:n, 1:6].sum(axis=1)/1000)

t_929, d_929, f_929 = load_version('OpenSees_0929')
t_943, d_943, f_943 = load_version('OpenSees_0943')
t_1045, d_1045, f_1045 = load_version('OpenSees_1045')

fig, axes = plt.subplots(1, 3, figsize=(20, 6))

# --- First cycle of 1045 (t < 1.6) ---
ax = axes[0]
mask = t_1045 < 1.6
ax.plot(d_1045[mask], f_1045[mask], '-', color='#2ca02c', linewidth=2.0,
        label='OpenSees-1045')
ax.plot(d_1045[mask][0], f_1045[mask][0], 'go', ms=10, zorder=5, label='Start')

# Mark key points
for t_target, marker, mlabel in [(0.3, 'v', 'Peak +4mm'),
                                  (1.1, '^', 'Peak -4mm')]:
    idx = np.argmin(np.abs(t_1045 - t_target))
    ax.plot(d_1045[idx], f_1045[idx], marker, color='red', ms=10, zorder=5,
            label=f'{mlabel}: {f_1045[idx]:.0f} kN')

ax.set_xlabel('Top Displacement (mm)', fontsize=13)
ax.set_ylabel('Base Shear (kN)', fontsize=13)
ax.set_title('OpenSees-1045: First Cycle Hysteresis', fontsize=14, fontweight='bold')
ax.legend(fontsize=10, loc='upper left')
ax.grid(True, alpha=0.3)
ax.axhline(y=0, color='k', linewidth=0.5)
ax.axvline(x=0, color='k', linewidth=0.5)

# --- All data of 1045 with force clipping ---
ax = axes[1]
f_clip = np.clip(f_1045, -300, 300)
ax.plot(d_1045, f_clip, '-', color='#2ca02c', linewidth=1.0, alpha=0.7)

# Highlight first cycle
mask1 = t_1045 < 1.6
ax.plot(d_1045[mask1], np.clip(f_1045[mask1], -300, 300), '-',
        color='#2ca02c', linewidth=2.0, label='First cycle')
# Second push
mask2 = t_1045 >= 1.6
ax.plot(d_1045[mask2], f_clip[mask2], '-', color='#d62728',
        linewidth=0.8, alpha=0.5, label='Second push (t>1.6)')

ax.set_xlabel('Top Displacement (mm)', fontsize=13)
ax.set_ylabel('Base Shear (kN, clipped ±300)', fontsize=13)
ax.set_title('OpenSees-1045: Full Analysis', fontsize=14, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
ax.axhline(y=0, color='k', linewidth=0.5)
ax.axvline(x=0, color='k', linewidth=0.5)

# --- Analysis progress comparison ---
ax = axes[2]
# Load protocol
prot = np.loadtxt(f'{WORK_DIR}/shuju2.txt')
prot_t = np.arange(1, len(prot)+1) * 0.1
prot_d = prot * 1000

ax.fill_between(prot_t[:30], prot_d[:30], alpha=0.1, color='gray', label='Target Protocol')
ax.plot(prot_t[:30], prot_d[:30], 'k-', linewidth=0.8, alpha=0.4)

ax.plot(t_929, d_929, '--', color='#1f77b4', linewidth=1.2, label=f'0929 (t→{t_929[-1]:.2f})')
ax.plot(t_943, d_943, '--', color='#ff7f0e', linewidth=1.2, label=f'0943 (t→{t_943[-1]:.2f})')
ax.plot(t_1045, d_1045, '-', color='#2ca02c', linewidth=2.0, label=f'1045 (t→{t_1045[-1]:.2f})')

# Add version end markers
ax.plot(t_929[-1], d_929[-1], 'x', color='#1f77b4', ms=12, mew=2)
ax.plot(t_943[-1], d_943[-1], 'x', color='#ff7f0e', ms=12, mew=2)
ax.plot(t_1045[-1], d_1045[-1], 'x', color='#2ca02c', ms=12, mew=2)

ax.set_xlabel('Pseudo-time', fontsize=13)
ax.set_ylabel('Top Displacement (mm)', fontsize=13)
ax.set_title('Analysis Progress Comparison', fontsize=14, fontweight='bold')
ax.legend(fontsize=10, loc='upper left')
ax.grid(True, alpha=0.3)
ax.set_xlim(-0.1, 4.0)

plt.suptitle('csp3.tcl Shear Wall: Section 6A Fix Results',
             fontsize=16, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(f'{WORK_DIR}/wall_1045_clean.png', dpi=150, bbox_inches='tight')
print("Saved: wall_1045_clean.png")

# Final summary
print("\n" + "=" * 60)
print("FINAL SUMMARY")
print("=" * 60)
print(f"\nOpenSees-0929: reached t={t_929[-1]:.3f} (first positive push, partial)")
print(f"OpenSees-0943: reached t={t_943[-1]:.3f} (first negative push)")
print(f"OpenSees-1045: reached t={t_1045[-1]:.3f} (second positive push, to ~5.5mm)")
print(f"\n1045 Improvement over 0943: {t_1045[-1]/t_943[-1]:.1f}x further")
print(f"1045 Improvement over 0929: {t_1045[-1]/t_929[-1]:.1f}x further")

# First cycle quality
mask = t_1045 < 1.6
f_fc = f_1045[mask]
print(f"\n1045 First Cycle (t<1.6):")
print(f"  Force range: [{f_fc.min():.1f}, {f_fc.max():.1f}] kN")
print(f"  Force StdDev: {np.std(np.diff(f_fc)):.1f} kN/step")
print(f"  Data points:  {np.sum(mask)}")

# Instability detection
large = np.abs(f_1045) > 500
idx_first_large = np.argmax(large) if np.any(large) else -1
if idx_first_large >= 0:
    print(f"\n  Force spike first appears at t={t_1045[idx_first_large]:.3f}, "
          f"disp={d_1045[idx_first_large]:.2f} mm")
