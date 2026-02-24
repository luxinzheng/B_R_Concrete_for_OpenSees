"""
Compare reference model results with OpenSees-1045 results.
Reference: 参考模型.txt (disp_m, force_N) - from csp3-参考模型.tcl with dt=0.1
Our 1045:  wall_OpenSees_1045_*.txt - from csp3.tcl with dt=0.01
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

WORK_DIR = r'e:\Basic\Concrete_Model\ADINA'

# ============================================================
# Load reference data
# ============================================================
ref_rows = []
with open(f'{WORK_DIR}/参考模型.txt', 'r') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2:
            try:
                d, force = float(parts[0]), float(parts[1])
                ref_rows.append([d, force])
            except ValueError:
                pass
ref = np.array(ref_rows)
# First 10 rows are gravity (tiny disp), skip them
ref_grav = ref[:10]
ref_lat = ref[10:]
ref_disp_mm = ref_lat[:, 0] * 1000
ref_force_kN = ref_lat[:, 1] / 1000

print(f"Reference: {len(ref_lat)} lateral points")
print(f"  Disp: [{ref_disp_mm.min():.1f}, {ref_disp_mm.max():.1f}] mm")
print(f"  Force: [{ref_force_kN.min():.1f}, {ref_force_kN.max():.1f}] kN")

# ============================================================
# Load 1045 data (from csp3.tcl run)
# ============================================================
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

d1045 = safe_load(f'{WORK_DIR}/wall_OpenSees_1045_disp.txt', 2)
r1045 = safe_load(f'{WORK_DIR}/wall_OpenSees_1045_react.txt', 6)
n = min(len(d1045), len(r1045))
t1045 = d1045[:n, 0]
# Find gravity/lateral boundary
lat_start = 0
for i in range(1, len(t1045)):
    if t1045[i] < t1045[i-1]:
        lat_start = i
        break
if lat_start == 0: lat_start = 10

our_disp_mm = d1045[lat_start:n, 1] * 1000
our_force_kN = r1045[lat_start:n, 1:6].sum(axis=1) / 1000
our_time = d1045[lat_start:n, 0]

print(f"\n1045: {len(our_disp_mm)} lateral points")
print(f"  Disp: [{our_disp_mm.min():.1f}, {our_disp_mm.max():.1f}] mm")
print(f"  Force: [{our_force_kN.min():.1f}, {our_force_kN.max():.1f}] kN")

# ============================================================
# Plot comparison
# ============================================================
fig, axes = plt.subplots(2, 3, figsize=(20, 12))

# --- Plot 1: Full reference hysteresis ---
ax = axes[0, 0]
ax.plot(ref_disp_mm, ref_force_kN, 'b-', linewidth=1.5, label='Reference Model')
ax.set_xlabel('Top Displacement (mm)')
ax.set_ylabel('Base Shear (kN)')
ax.set_title('Reference Model: Full Hysteresis', fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3)
ax.axhline(y=0, color='k', linewidth=0.5)
ax.axvline(x=0, color='k', linewidth=0.5)

# --- Plot 2: 1045 hysteresis (first cycle only, clean) ---
ax = axes[0, 1]
mask = our_time < 1.6
ax.plot(our_disp_mm[mask], our_force_kN[mask], 'g-', linewidth=1.5, label='1045 (first cycle)')
ax.set_xlabel('Top Displacement (mm)')
ax.set_ylabel('Base Shear (kN)')
ax.set_title('1045: First Cycle', fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3)
ax.axhline(y=0, color='k', linewidth=0.5)
ax.axvline(x=0, color='k', linewidth=0.5)

# --- Plot 3: Overlay first cycle ---
ax = axes[0, 2]
# Extract reference first cycle (first 15 points: 0->+4->0->-4->0)
ref_fc = ref_lat[:15]
ref_fc_disp = ref_fc[:, 0] * 1000
ref_fc_force = ref_fc[:, 1] / 1000
ax.plot(ref_fc_disp, ref_fc_force, 'b-o', linewidth=2, markersize=4, label='Reference')
ax.plot(our_disp_mm[mask], our_force_kN[mask], 'g-', linewidth=1.5, alpha=0.8, label='1045')
ax.set_xlabel('Top Displacement (mm)')
ax.set_ylabel('Base Shear (kN)')
ax.set_title('First Cycle Overlay', fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3)
ax.axhline(y=0, color='k', linewidth=0.5)
ax.axvline(x=0, color='k', linewidth=0.5)

# --- Plot 4: Force at key displacement points ---
ax = axes[1, 0]
# Extract peak forces from reference at each cycle
ref_peaks_pos = []
ref_peaks_neg = []
prev_d = 0
for i in range(len(ref_lat)):
    d = ref_lat[i, 0]
    f = ref_lat[i, 1] / 1000
    if i > 0 and ref_lat[i-1, 0] > d and d > 0:
        ref_peaks_pos.append((ref_lat[i-1, 0] * 1000, ref_lat[i-1, 1] / 1000))
    if i > 0 and ref_lat[i-1, 0] < d and d < 0:
        ref_peaks_neg.append((ref_lat[i-1, 0] * 1000, ref_lat[i-1, 1] / 1000))

# Backbone comparison
ref_disp_peaks = [p[0] for p in ref_peaks_pos]
ref_force_peaks = [p[1] for p in ref_peaks_pos]
ax.plot(ref_disp_peaks, ref_force_peaks, 'bs-', linewidth=2, label='Ref positive peaks')

ref_disp_peaks_n = [p[0] for p in ref_peaks_neg]
ref_force_peaks_n = [p[1] for p in ref_peaks_neg]
ax.plot(ref_disp_peaks_n, ref_force_peaks_n, 'bv-', linewidth=2, label='Ref negative peaks')

ax.set_xlabel('Peak Displacement (mm)')
ax.set_ylabel('Peak Force (kN)')
ax.set_title('Backbone Curve (Reference)', fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3)

# --- Plot 5: Force-displacement point comparison at integer mm ---
ax = axes[1, 1]
# Reference: find force at each integer mm displacement
ref_milestones = {}
for i in range(len(ref_lat)):
    d_mm = round(ref_lat[i, 0] * 1000, 1)
    f_kN = ref_lat[i, 1] / 1000
    key = int(d_mm)
    if key not in ref_milestones and abs(d_mm - key) < 0.01:
        ref_milestones[key] = f_kN

# 1045: find force at each integer mm displacement
our_milestones = {}
for i in range(len(our_disp_mm)):
    d_mm = round(our_disp_mm[i], 1)
    f_kN = our_force_kN[i]
    key = int(d_mm)
    if key not in our_milestones and abs(d_mm - key) < 0.5 and abs(f_kN) < 300:
        our_milestones[key] = f_kN

common = sorted(set(ref_milestones.keys()) & set(our_milestones.keys()))
if common:
    ref_f = [ref_milestones[k] for k in common]
    our_f = [our_milestones[k] for k in common]
    width = 0.35
    x = np.arange(len(common))
    ax.bar(x - width/2, ref_f, width, label='Reference', color='#1f77b4')
    ax.bar(x + width/2, our_f, width, label='1045', color='#2ca02c')
    ax.set_xticks(x)
    ax.set_xticklabels([f'{k}mm' for k in common], rotation=45)
    ax.set_ylabel('Force (kN)')
    ax.set_title('Force at Key Displacements', fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')

# --- Plot 6: Summary comparison table ---
ax = axes[1, 2]
ax.axis('off')
summary = []
summary.append(['Metric', 'Reference', '1045'])
summary.append(['Lateral pts', str(len(ref_lat)), str(len(our_disp_mm))])
summary.append(['Disp range', f'[{ref_disp_mm.min():.0f}, {ref_disp_mm.max():.0f}] mm',
                 f'[{our_disp_mm.min():.0f}, {our_disp_mm.max():.0f}] mm'])
summary.append(['Force range', f'[{ref_force_kN.min():.0f}, {ref_force_kN.max():.0f}] kN',
                 f'[{our_force_kN[np.abs(our_force_kN)<500].min():.0f}, {our_force_kN[np.abs(our_force_kN)<500].max():.0f}] kN'])

# First cycle peak comparison
summary.append(['', '', ''])
summary.append(['First Cycle', '---', '---'])
# Ref: +4mm
ref_4mm_idx = 2  # row 13 (index 2 in ref_lat)
summary.append([f'+4mm force', f'{ref_lat[2, 1]/1000:.1f} kN', ''])
# Ref: -4mm
ref_m4mm_idx = 10
summary.append([f'-4mm force', f'{ref_lat[10, 1]/1000:.1f} kN', ''])

# 1045 first peaks
mask_pos = (our_disp_mm > 3.8) & (our_disp_mm < 4.2) & (our_time < 0.4)
if np.any(mask_pos):
    our_4mm = our_force_kN[mask_pos]
    summary[6][2] = f'{our_4mm[np.argmax(np.abs(our_4mm))]:.1f} kN'

mask_neg = (our_disp_mm < -3.8) & (our_disp_mm > -4.2) & (our_time < 1.2)
if np.any(mask_neg):
    our_m4mm = our_force_kN[mask_neg]
    summary[7][2] = f'{our_m4mm[np.argmax(np.abs(our_m4mm))]:.1f} kN'

summary.append(['', '', ''])
summary.append(['dt', '0.1', '0.01'])
summary.append(['nProps', '7', '37'])
summary.append(['nStatevs', '40', '13'])

table = ax.table(cellText=summary, cellLoc='center', loc='center',
                 colWidths=[0.35, 0.35, 0.30])
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 1.5)
for (r, c), cell in table.get_celld().items():
    if r == 0:
        cell.set_facecolor('#d4e6f1')
        cell.set_text_props(weight='bold')
    elif summary[r][0] in ['First Cycle', '']:
        cell.set_facecolor('#f2f3f4')
ax.set_title('Summary Comparison', fontweight='bold', pad=20)

plt.suptitle('Reference Model vs OpenSees-1045 Comparison',
             fontsize=16, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(f'{WORK_DIR}/ref_vs_1045.png', dpi=150, bbox_inches='tight')
print(f"\nSaved: ref_vs_1045.png")

# Print force comparison at key points
print("\n" + "="*60)
print("FORCE COMPARISON AT KEY DISPLACEMENT POINTS")
print("="*60)
print(f"{'Disp(mm)':>10} {'Ref(kN)':>12} {'1045(kN)':>12} {'Ratio':>8}")
for k in sorted(common):
    r = ref_milestones[k]
    o = our_milestones[k]
    ratio = o/r if abs(r) > 0.1 else float('inf')
    print(f"{k:>10} {r:>12.1f} {o:>12.1f} {ratio:>8.2f}")
