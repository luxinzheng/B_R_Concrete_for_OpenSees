"""
Plot hysteresis curve from OpenSees csp3.tcl analysis results.
- X-axis: displacement at node 53 (top of wall)
- Y-axis: base shear (negative sum of reactions at nodes 1-5)
Handles: time reset after gravity, force-through artifacts, truncated lines
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['font.family'] = ['Microsoft YaHei', 'SimHei', 'sans-serif']

def safe_load(filepath, ncols):
    """Load text file, skipping malformed lines."""
    rows = []
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == ncols:
                try:
                    vals = [float(p) for p in parts]
                    rows.append(vals)
                except ValueError:
                    pass
    return np.array(rows, dtype=float)

# --- Load data ---
data53 = safe_load(r'e:\Basic\Concrete_Model\ADINA\53.txt', 2)
data1  = safe_load(r'e:\Basic\Concrete_Model\ADINA\1.txt', 6)

# --- Separate gravity (first 10 steps) and cyclic phases ---
n_gravity = 10
cyc_disp  = data53[n_gravity:, :]
cyc_react = data1[n_gravity:, :]
n_cyc = min(len(cyc_disp), len(cyc_react))
cyc_disp  = cyc_disp[:n_cyc]
cyc_react = cyc_react[:n_cyc]

time_cyc = cyc_disp[:, 0]
disp_mm  = cyc_disp[:, 1] * 1000.0

# Base shear = -(sum of reactions)
reactions_sum = np.sum(cyc_react[:, 1:6], axis=1)
shear_kN = -reactions_sum / 1000.0

# --- Intelligent filtering ---
# 1. Remove points with |V| > threshold (force-through artifacts)
threshold = 300.0  # kN - well above expected max (~200 kN) but below force-through spikes
good = np.abs(shear_kN) < threshold

# 2. Also remove sudden jumps in shear (gradient filter)
dV = np.abs(np.diff(shear_kN))
dV = np.append(dV, 0)
jump_threshold = 100.0  # kN per step
no_jump = dV < jump_threshold

# 3. Also check previous step (remove both sides of jump)
dV_prev = np.append(0, np.abs(np.diff(shear_kN)))
no_jump_prev = dV_prev < jump_threshold

valid = good & no_jump & no_jump_prev

time_v  = time_cyc[valid]
disp_v  = disp_mm[valid]
shear_v = shear_kN[valid]

n_removed = np.sum(~valid)
print(f"Total cyclic points: {n_cyc}")
print(f"Valid points after filtering: {len(time_v)} (removed {n_removed})")
print(f"Time range: {time_v[0]:.4f} ~ {time_v[-1]:.4f}")
print(f"Displacement range: {disp_v.min():.3f} ~ {disp_v.max():.3f} mm")
print(f"Base shear range: {shear_v.min():.1f} ~ {shear_v.max():.1f} kN")

# --- Segment the data to avoid connecting across gaps ---
# Where time difference between consecutive valid points is large, break the line
dt = np.diff(time_v)
gap_idx = np.where(dt > 0.05)[0]  # gaps larger than 0.05 time units
segments = []
start = 0
for gi in gap_idx:
    segments.append((start, gi + 1))
    start = gi + 1
segments.append((start, len(time_v)))

print(f"Number of continuous segments: {len(segments)}")

# --- Load target displacement path ---
shuju = np.loadtxt(r'e:\Basic\Concrete_Model\ADINA\shuju2.txt')
time_target = np.arange(len(shuju)) * 0.1
disp_target_mm = shuju * 1000.0
t_max_cyc = time_v[-1]

# ====== PLOT 1: Four-panel overview ======
fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# (1) Hysteresis curve
ax1 = axes[0, 0]
for i, (s, e) in enumerate(segments):
    ax1.plot(disp_v[s:e], shear_v[s:e], 'b-', linewidth=0.8, alpha=0.8)
ax1.set_xlabel('位移 (mm)', fontsize=12)
ax1.set_ylabel('底部剪力 (kN)', fontsize=12)
ax1.set_title('滞回曲线 (底部剪力 vs 顶部位移)', fontsize=13)
ax1.grid(True, alpha=0.3)
ax1.axhline(y=0, color='k', linewidth=0.5)
ax1.axvline(x=0, color='k', linewidth=0.5)

# (2) Displacement time history
ax2 = axes[0, 1]
for s, e in segments:
    ax2.plot(time_v[s:e], disp_v[s:e], 'b-', linewidth=0.8)
mask_t = time_target <= t_max_cyc + 0.5
ax2.plot(time_target[mask_t], disp_target_mm[mask_t],
         'r--', linewidth=0.8, alpha=0.6, label='目标位移')
ax2.set_xlabel('时间', fontsize=12)
ax2.set_ylabel('位移 (mm)', fontsize=12)
ax2.set_title('位移时程', fontsize=13)
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3)

# (3) Base shear time history
ax3 = axes[1, 0]
for s, e in segments:
    ax3.plot(time_v[s:e], shear_v[s:e], 'b-', linewidth=0.8)
ax3.set_xlabel('时间', fontsize=12)
ax3.set_ylabel('底部剪力 (kN)', fontsize=12)
ax3.set_title('底部剪力时程', fontsize=13)
ax3.grid(True, alpha=0.3)

# (4) Full loading path vs analysis progress
ax4 = axes[1, 1]
ax4.plot(time_target, disp_target_mm, 'r-', linewidth=0.5, alpha=0.7, label='完整加载路径')
ax4.axvline(x=t_max_cyc, color='b', linewidth=1.5, linestyle='--',
            label=f'分析终止 t={t_max_cyc:.2f}')
ax4.set_xlabel('时间', fontsize=12)
ax4.set_ylabel('位移 (mm)', fontsize=12)
ax4.set_title('完整加载路径 vs 分析进度', fontsize=13)
ax4.legend(fontsize=10)
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(r'e:\Basic\Concrete_Model\ADINA\plot_opensees_hysteresis.png', dpi=150)
plt.close()
print("\n四面板图已保存: plot_opensees_hysteresis.png")

# ====== PLOT 2: Detailed hysteresis with color-coded phases ======
fig2, ax = plt.subplots(figsize=(14, 9))

# Color phases by time
phase_defs = [
    (0.0,  0.30, '正向推至 +4mm',   'tab:blue',   1.5),
    (0.30, 0.70, '从 +4mm 卸载',    'tab:cyan',   1.2),
    (0.70, 1.10, '反向推至 -4mm',   'tab:red',    1.5),
    (1.10, 1.50, '从 -4mm 回弹',    'tab:orange', 1.2),
]

for t_start, t_end, label, color, lw in phase_defs:
    for s, e in segments:
        mask = (time_v[s:e] >= t_start) & (time_v[s:e] <= t_end)
        idx_masked = np.arange(s, e)[mask]
        if len(idx_masked) > 1:
            ax.plot(disp_v[idx_masked], shear_v[idx_masked], '-',
                    color=color, linewidth=lw, alpha=0.9)
    # Add one invisible line for legend
    ax.plot([], [], '-', color=color, linewidth=lw, label=label)

# Mark key points
idx_max_d = np.argmax(disp_v)
idx_min_d = np.argmin(disp_v)
idx_max_v = np.argmax(shear_v)
idx_min_v = np.argmin(shear_v)

ax.plot(disp_v[idx_max_d], shear_v[idx_max_d], 'ro', markersize=10, zorder=5,
        label=f'最大位移: {disp_v[idx_max_d]:.1f}mm, V={shear_v[idx_max_d]:.0f}kN')
ax.plot(disp_v[idx_min_d], shear_v[idx_min_d], 'gs', markersize=10, zorder=5,
        label=f'最小位移: {disp_v[idx_min_d]:.1f}mm, V={shear_v[idx_min_d]:.0f}kN')
ax.plot(disp_v[idx_max_v], shear_v[idx_max_v], 'm^', markersize=10, zorder=5,
        label=f'最大剪力: {shear_v[idx_max_v]:.0f}kN @ {disp_v[idx_max_v]:.1f}mm')
ax.plot(disp_v[idx_min_v], shear_v[idx_min_v], 'cv', markersize=10, zorder=5,
        label=f'最小剪力: {shear_v[idx_min_v]:.0f}kN @ {disp_v[idx_min_v]:.1f}mm')

ax.set_xlabel('顶部位移 (mm)', fontsize=14)
ax.set_ylabel('底部剪力 (kN)', fontsize=14)
ax.set_title('OpenSees csp3 - 低周往复滞回曲线', fontsize=16)
ax.grid(True, alpha=0.3)
ax.axhline(y=0, color='k', linewidth=0.5)
ax.axvline(x=0, color='k', linewidth=0.5)
ax.legend(fontsize=9, loc='lower right', ncol=2)

plt.tight_layout()
plt.savefig(r'e:\Basic\Concrete_Model\ADINA\plot_opensees_hysteresis_detail.png', dpi=150)
plt.close()
print("详细滞回图已保存: plot_opensees_hysteresis_detail.png")

# ====== PLOT 3: Backbone curve (envelope) ======
fig3, ax3 = plt.subplots(figsize=(10, 7))

# Extract backbone: for each displacement level, find max/min shear
disp_bins_pos = np.arange(0, disp_v.max() + 0.5, 0.5)
disp_bins_neg = np.arange(disp_v.min(), 0.5, 0.5)

backbone_pos_d = []
backbone_pos_v = []
backbone_neg_d = []
backbone_neg_v = []

for d in disp_bins_pos:
    mask = (disp_v >= d - 0.25) & (disp_v <= d + 0.25) & (shear_v > 0)
    if np.sum(mask) > 0:
        backbone_pos_d.append(d)
        backbone_pos_v.append(np.max(shear_v[mask]))

for d in disp_bins_neg:
    mask = (disp_v >= d - 0.25) & (disp_v <= d + 0.25) & (shear_v < 0)
    if np.sum(mask) > 0:
        backbone_neg_d.append(d)
        backbone_neg_v.append(np.min(shear_v[mask]))

# Plot hysteresis as background
for s, e in segments:
    ax3.plot(disp_v[s:e], shear_v[s:e], 'b-', linewidth=0.4, alpha=0.3)

# Plot backbone
if backbone_pos_d:
    ax3.plot(backbone_pos_d, backbone_pos_v, 'r-o', linewidth=2, markersize=5,
             label='正向骨架曲线')
if backbone_neg_d:
    ax3.plot(backbone_neg_d, backbone_neg_v, 'g-s', linewidth=2, markersize=5,
             label='负向骨架曲线')

ax3.set_xlabel('顶部位移 (mm)', fontsize=14)
ax3.set_ylabel('底部剪力 (kN)', fontsize=14)
ax3.set_title('滞回曲线与骨架曲线', fontsize=16)
ax3.grid(True, alpha=0.3)
ax3.axhline(y=0, color='k', linewidth=0.5)
ax3.axvline(x=0, color='k', linewidth=0.5)
ax3.legend(fontsize=11)

plt.tight_layout()
plt.savefig(r'e:\Basic\Concrete_Model\ADINA\plot_opensees_backbone.png', dpi=150)
plt.close()
print("骨架曲线图已保存: plot_opensees_backbone.png")

# --- Summary ---
print(f"\n{'='*50}")
print(f"分析结果摘要")
print(f"{'='*50}")
print(f"分析终止时间: t = {t_max_cyc:.4f}")
disp_at_end = np.interp(t_max_cyc, time_target, disp_target_mm)
print(f"终止时目标位移: {disp_at_end:.2f} mm")
print(f"终止时计算位移: {disp_v[-1]:.2f} mm")
print(f"最大正向位移: {disp_v.max():.2f} mm (t={time_v[idx_max_d]:.3f})")
print(f"最大负向位移: {disp_v.min():.2f} mm (t={time_v[idx_min_d]:.3f})")
print(f"最大正向剪力: {shear_v.max():.1f} kN")
print(f"最大负向剪力: {shear_v.min():.1f} kN")
print(f"\n加载历程:")
print(f"  第1半周: +4mm → 最大正向剪力 ~{shear_v[(time_v >= 0) & (time_v <= 0.3)].max():.0f} kN")
v_neg = shear_v[(time_v >= 0.7) & (time_v <= 1.1)]
if len(v_neg) > 0:
    print(f"  第2半周: -4mm → 最大负向剪力 ~{v_neg.min():.0f} kN")
print(f"  分析在从-4mm回弹到约-1.2mm时失败 (刚度矩阵奇异)")
