# -*- coding: utf-8 -*-
"""
Final comparison: OpenSees-02210959.exe (B&R ADINA, forumat.f90 restored)
vs Reference (HEOM2D 7-param) for Multi-layer_Shell
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

plt.rcParams['font.family'] = 'Microsoft YaHei'
plt.rcParams['axes.unicode_minus'] = False

root = r'e:\Basic\Concrete_Model\ADINA'
ml = os.path.join(root, 'Multi-layer_Shell')

def load_safe(path, min_cols):
    rows = []
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return None
    with open(path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= min_cols:
                try:
                    rows.append([float(x) for x in parts[:min_cols]])
                except ValueError:
                    pass
    return np.array(rows) if rows else None

ref_d = load_safe(os.path.join(ml, 'ref_disp1.txt'), 2)
ref_r = load_safe(os.path.join(ml, 'ref_shearforce1.txt'), 10)
br_d  = load_safe(os.path.join(ml, 'test_disp.txt'), 2)
br_r  = load_safe(os.path.join(ml, 'test_react.txt'), 10)

nref = min(len(ref_d), len(ref_r))
ref_disp  = ref_d[:nref, 1] * 1000.0
ref_shear = ref_r[:nref, 1:].sum(axis=1) / 1000.0

nbr = min(len(br_d), len(br_r))
br_disp  = br_d[:nbr, 1] * 1000.0
br_shear = br_r[:nbr, 1:].sum(axis=1) / 1000.0

ref_peak = max(abs(ref_shear.max()), abs(ref_shear.min()))
print(f'Ref: {nref} pts, disp [{ref_disp.min():.1f}, {ref_disp.max():.1f}] mm, '
      f'peak shear = {ref_peak:.0f} kN')
print(f'B&R: {nbr} pts, disp [{br_disp.min():.1f}, {br_disp.max():.1f}] mm')

# Classify B&R points
reasonable = np.abs(br_shear) < ref_peak * 3
n_reasonable = reasonable.sum()
idx_first_bad = np.where(~reasonable)[0][0] if (~reasonable).any() else nbr
print(f'B&R: {n_reasonable}/{nbr} pts reasonable (|shear|<{ref_peak*3:.0f} kN)')
print(f'B&R: first divergent step = {idx_first_bad}, disp={br_disp[idx_first_bad]:.1f} mm')
print(f'B&R shear at diverge: {br_shear[idx_first_bad]:.0f} kN, max: {br_shear.max():.0f}, min: {br_shear.min():.0f} kN')

fig = plt.figure(figsize=(20, 14))
fig.suptitle('Multi-layer Shell 对比: OpenSees-02210959 (B&R ADINA)  vs  参考模型 (HEOM2D)\n'
             'forumat.f90 未修改 / test_br.tcl 分析策略已修复 (24/33 blocks完成)',
             fontsize=14, fontweight='bold')

# Row 1: Hysteresis comparisons
# 1a: Reference full hysteresis
ax1 = fig.add_subplot(2, 3, 1)
ax1.plot(ref_disp, ref_shear, 'b-', lw=1.0, alpha=0.7)
ax1.set_xlabel('顶部位移 (mm)')
ax1.set_ylabel('底部剪力 (kN)')
ax1.set_title(f'参考模型 HEOM2D\n({nref} 步全部完成, ±{ref_peak:.0f} kN)')
ax1.grid(True, alpha=0.3)
ax1.axhline(0, color='k', lw=0.5); ax1.axvline(0, color='k', lw=0.5)

# 1b: B&R only reasonable points
ax2 = fig.add_subplot(2, 3, 2)
# Only plot contiguous good region
br_d_ok = br_disp[:idx_first_bad]
br_s_ok = br_shear[:idx_first_bad]
ax2.plot(br_d_ok, br_s_ok, 'r-', lw=1.0, alpha=0.7)
ax2.set_xlabel('顶部位移 (mm)')
ax2.set_ylabel('底部剪力 (kN)')
ax2.set_title(f'B&R ADINA (漂移前 {idx_first_bad} 步)\n'
              f'shear [{br_s_ok.min():.0f}, {br_s_ok.max():.0f}] kN')
ax2.grid(True, alpha=0.3)
ax2.axhline(0, color='k', lw=0.5); ax2.axvline(0, color='k', lw=0.5)

# 1c: Overlay — same Y axis forced
ax3 = fig.add_subplot(2, 3, 3)
y_max = max(ref_peak, max(abs(br_s_ok.max()), abs(br_s_ok.min()))) * 1.2
ax3.plot(ref_disp, ref_shear, 'b-', lw=1.0, alpha=0.6, label='参考 HEOM2D')
ax3.plot(br_d_ok, br_s_ok, 'r-', lw=1.0, alpha=0.7, label=f'B&R (前 {idx_first_bad} 步)')
ax3.set_ylim(-y_max, y_max)
ax3.set_xlabel('顶部位移 (mm)')
ax3.set_ylabel('底部剪力 (kN)')
ax3.set_title('两模型叠加对比')
ax3.legend(fontsize=9)
ax3.grid(True, alpha=0.3)
ax3.axhline(0, color='k', lw=0.5); ax3.axvline(0, color='k', lw=0.5)

# Row 2: Time histories
# 2a: Displacement history
ax4 = fig.add_subplot(2, 3, 4)
ax4.plot(range(nref), ref_disp, 'b-', lw=0.8, alpha=0.7, label=f'参考 ({nref} 步)')
ax4.plot(range(len(br_d)), br_d[:, 1]*1000, 'r-', lw=1.0, label=f'B&R ({len(br_d)} 步)')
ax4.axvline(idx_first_bad, color='orange', ls='--', lw=1.2, alpha=0.8,
            label=f'应力漂移起始 (step {idx_first_bad})')
ax4.set_xlabel('步号')
ax4.set_ylabel('顶部位移 (mm)')
ax4.set_title('位移加载历史')
ax4.legend(fontsize=8)
ax4.grid(True, alpha=0.3)

# 2b: Shear force history (log)
ax5 = fig.add_subplot(2, 3, 5)
ax5.semilogy(range(nref), np.abs(ref_shear)+1, 'b-', lw=0.8, alpha=0.7, label='|参考 shear| + 1')
ax5.semilogy(range(nbr), np.abs(br_shear)+1, 'r-', lw=0.8, alpha=0.7, label='|B&R shear| + 1')
ax5.axhline(ref_peak, color='b', ls=':', lw=1, alpha=0.5, label=f'参考峰值 {ref_peak:.0f} kN')
ax5.axvline(idx_first_bad, color='orange', ls='--', lw=1.2, alpha=0.8, label=f'漂移起始')
ax5.set_xlabel('步号')
ax5.set_ylabel('|底部剪力| (kN, 对数)')
ax5.set_title('剪力幅值历程 (对数坐标)')
ax5.legend(fontsize=8)
ax5.grid(True, alpha=0.3)

# 2c: B&R shear force (linear, clipped)
ax6 = fig.add_subplot(2, 3, 6)
clip_v = ref_peak * 5
ax6.plot(range(nref), ref_shear, 'b-', lw=0.8, alpha=0.7, label='参考')
ax6.plot(range(nbr), np.clip(br_shear, -clip_v, clip_v), 'r-', lw=0.8, alpha=0.7,
         label=f'B&R (clip ±{clip_v:.0f} kN)')
ax6.axvline(idx_first_bad, color='orange', ls='--', lw=1.2, alpha=0.8, label='漂移起始')
ax6.set_xlabel('步号')
ax6.set_ylabel('底部剪力 (kN)')
ax6.set_title(f'剪力历程 (线性, clip ±{clip_v:.0f} kN)')
ax6.legend(fontsize=8)
ax6.grid(True, alpha=0.3)

plt.tight_layout(rect=[0, 0, 1, 0.92])
out = os.path.join(root, 'compare_02210959_final.png')
plt.savefig(out, dpi=150)
print(f'\nSaved: {out}')

# ======== sw1-1 comparison (data may be from 02212103, behavior equivalent) ========
fig2, axes2 = plt.subplots(1, 3, figsize=(18, 5))
fig2.suptitle('sw1-1 对比: B&R ADINA vs 参考模型', fontsize=13, fontweight='bold')

sw1 = os.path.join(root, 'sw1-1')
sw1_ref_d = load_safe(os.path.join(sw1, 'ref_53.txt'), 2)
sw1_ref_r = load_safe(os.path.join(sw1, 'ref_1.txt'), 6)
sw1_br_d  = load_safe(os.path.join(sw1, '53.txt'), 2)
sw1_br_r  = load_safe(os.path.join(sw1, '1.txt'), 6)

if sw1_ref_d is not None and sw1_ref_r is not None and sw1_br_d is not None and sw1_br_r is not None:
    nref1 = min(len(sw1_ref_d), len(sw1_ref_r))
    sw1_ref_disp = sw1_ref_d[:nref1, 1]*1000
    sw1_ref_shear = sw1_ref_r[:nref1, 1:].sum(axis=1)/1000

    nbr1 = min(len(sw1_br_d), len(sw1_br_r))
    sw1_br_disp = sw1_br_d[:nbr1, 1]*1000
    sw1_br_shear = sw1_br_r[:nbr1, 1:].sum(axis=1)/1000

    ax = axes2[0]
    ax.plot(sw1_ref_disp, sw1_ref_shear, 'b-', lw=1.2, alpha=0.7, label=f'参考 ({nref1} pts)')
    ax.set_xlabel('位移 (mm)'); ax.set_ylabel('剪力 (kN)')
    ax.set_title('参考模型 滞回'); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    ax = axes2[1]
    ax.plot(sw1_br_disp, sw1_br_shear, 'r-', lw=1.0, alpha=0.7, label=f'B&R ({nbr1} pts)')
    ax.set_xlabel('位移 (mm)'); ax.set_ylabel('剪力 (kN)')
    ax.set_title('B&R ADINA 滞回'); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    ax = axes2[2]
    ax.plot(sw1_ref_disp, sw1_ref_shear, 'b-', lw=1.0, alpha=0.6, label='参考')
    ax.plot(sw1_br_disp, sw1_br_shear, 'r-', lw=1.0, alpha=0.7, label='B&R')
    ax.set_xlabel('位移 (mm)'); ax.set_ylabel('剪力 (kN)')
    ax.set_title('叠加对比'); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    print(f'\nsw1-1 Ref: {nref1} pts, shear [{sw1_ref_shear.min():.0f}, {sw1_ref_shear.max():.0f}] kN')
    print(f'sw1-1 B&R: {nbr1} pts, shear [{sw1_br_shear.min():.0f}, {sw1_br_shear.max():.0f}] kN')

plt.tight_layout()
out2 = os.path.join(root, 'compare_sw1_final.png')
plt.savefig(out2, dpi=150)
print(f'Saved: {out2}')
