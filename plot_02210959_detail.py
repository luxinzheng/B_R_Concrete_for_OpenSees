# -*- coding: utf-8 -*-
"""
Detailed comparison: OpenSees-02210959.exe (B&R ADINA, forumat.f90 unmodified)
vs Reference (HEOM2D 7-param) for Multi-layer_Shell test.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

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
ref_disp  = ref_d[:nref, 1] * 1000.0       # mm
ref_shear = ref_r[:nref, 1:].sum(axis=1) / 1000.0  # kN

nbr = min(len(br_d), len(br_r))
br_disp  = br_d[:nbr, 1] * 1000.0
br_shear = br_r[:nbr, 1:].sum(axis=1) / 1000.0

print(f'Ref: {nref} pts, disp [{ref_disp.min():.1f}, {ref_disp.max():.1f}] mm, '
      f'shear [{ref_shear.min():.0f}, {ref_shear.max():.0f}] kN')
print(f'B&R: {nbr} pts, disp [{br_disp.min():.1f}, {br_disp.max():.1f}] mm, '
      f'shear [{br_shear.min():.0f}, {br_shear.max():.0f}] kN')

# Find where B&R shear starts to go unrealistic (>200 kN or <-200 kN)
ref_shear_range = max(abs(ref_shear.max()), abs(ref_shear.min()))
threshold = ref_shear_range * 2
good_mask = np.abs(br_shear) < threshold
if good_mask.any():
    last_good = np.where(good_mask)[0][-1]
else:
    last_good = 0
print(f'Ref shear range: ±{ref_shear_range:.0f} kN')
print(f'B&R reasonable until step {last_good}/{nbr}, disp={br_disp[last_good]:.1f} mm')

fig = plt.figure(figsize=(20, 14))
fig.suptitle('Multi-layer Shell: OpenSees-02210959 (B&R ADINA)  vs  Reference (HEOM2D)\n'
             'forumat.f90 未修改 / test_br.tcl 已修复分析策略',
             fontsize=14, fontweight='bold')

# --- Row 1: Hysteresis ---
# 1a: Full range (clipped for visualization)
ax1 = fig.add_subplot(2, 3, 1)
ax1.plot(ref_disp, ref_shear, 'b-', lw=1.0, alpha=0.7, label=f'Ref HEOM2D ({nref} pts)')
clip_kn = ref_shear_range * 3
br_clip = np.clip(br_shear, -clip_kn, clip_kn)
ax1.plot(br_disp[:nbr], br_clip[:nbr], 'r-', lw=0.8, alpha=0.7, label=f'B&R (clip ±{clip_kn:.0f}kN)')
ax1.set_xlabel('Top Displacement (mm)')
ax1.set_ylabel('Base Shear (kN)')
ax1.set_title(f'Hysteresis (clipped ±{clip_kn:.0f}kN)')
ax1.legend(fontsize=8)
ax1.grid(True, alpha=0.3)
ax1.axhline(0, color='k', lw=0.5)
ax1.axvline(0, color='k', lw=0.5)

# 1b: Zoom into reasonable range only
ax2 = fig.add_subplot(2, 3, 2)
ax2.plot(ref_disp, ref_shear, 'b-', lw=1.2, alpha=0.8, label='Ref HEOM2D')
lg = last_good + 1
ax2.plot(br_disp[:lg], br_shear[:lg], 'r-', lw=1.0, alpha=0.7, label=f'B&R (first {lg} pts)')
ax2.set_xlabel('Top Displacement (mm)')
ax2.set_ylabel('Base Shear (kN)')
ax2.set_title(f'Hysteresis (B&R first {lg} reasonable pts)')
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.3)
ax2.axhline(0, color='k', lw=0.5)
ax2.axvline(0, color='k', lw=0.5)

# 1c: Just first few blocks (small displacement range, best comparison)
ax3 = fig.add_subplot(2, 3, 3)
disp_limit = 6.0  # mm
ref_mask = np.abs(ref_disp) <= disp_limit * 1.2
br_early_mask = np.abs(br_disp[:lg]) <= disp_limit * 1.2
ax3.plot(ref_disp[ref_mask], ref_shear[ref_mask], 'b-', lw=1.2, alpha=0.8, label='Ref HEOM2D')
if br_early_mask.any():
    ax3.plot(br_disp[:lg][br_early_mask], br_shear[:lg][br_early_mask], 'r-', lw=1.0, alpha=0.7, label='B&R ADINA')
ax3.set_xlabel('Top Displacement (mm)')
ax3.set_ylabel('Base Shear (kN)')
ax3.set_title(f'Hysteresis (|disp| ≤ {disp_limit} mm)')
ax3.legend(fontsize=8)
ax3.grid(True, alpha=0.3)
ax3.axhline(0, color='k', lw=0.5)
ax3.axvline(0, color='k', lw=0.5)

# --- Row 2: Time histories ---
# 2a: Displacement history
ax4 = fig.add_subplot(2, 3, 4)
ax4.plot(range(nref), ref_disp, 'b-', lw=0.8, alpha=0.7, label=f'Ref ({nref} steps)')
ax4.plot(range(len(br_d)), br_d[:, 1]*1000, 'r-', lw=1.0, label=f'B&R ({len(br_d)} steps)')
ax4.axvline(last_good, color='orange', ls='--', lw=1, alpha=0.7, label=f'B&R diverge ≈ step {last_good}')
ax4.set_xlabel('Step Number')
ax4.set_ylabel('Top Disp (mm)')
ax4.set_title('Displacement History')
ax4.legend(fontsize=8)
ax4.grid(True, alpha=0.3)

# 2b: Shear force history
ax5 = fig.add_subplot(2, 3, 5)
ax5.plot(range(nref), ref_shear, 'b-', lw=0.8, alpha=0.7, label='Ref')
br_shear_log = br_shear.copy()
ax5.plot(range(nbr), np.clip(br_shear_log, -clip_kn, clip_kn), 'r-', lw=0.8, alpha=0.7, label=f'B&R (clip ±{clip_kn:.0f}kN)')
ax5.axvline(last_good, color='orange', ls='--', lw=1, alpha=0.7, label=f'Diverge ≈ step {last_good}')
ax5.set_xlabel('Step Number')
ax5.set_ylabel('Base Shear (kN)')
ax5.set_title('Shear Force History')
ax5.legend(fontsize=8)
ax5.grid(True, alpha=0.3)

# 2c: B&R shear force in log scale (show drift magnitude)
ax6 = fig.add_subplot(2, 3, 6)
ax6.semilogy(range(nref), np.abs(ref_shear)+1, 'b-', lw=0.8, alpha=0.7, label='|Ref shear|')
ax6.semilogy(range(nbr), np.abs(br_shear)+1, 'r-', lw=0.8, alpha=0.7, label='|B&R shear|')
ax6.axhline(ref_shear_range, color='b', ls=':', lw=0.8, alpha=0.5)
ax6.axvline(last_good, color='orange', ls='--', lw=1, alpha=0.7, label=f'Diverge ≈ step {last_good}')
ax6.set_xlabel('Step Number')
ax6.set_ylabel('|Base Shear| (kN, log)')
ax6.set_title('Shear Force Magnitude (log scale)')
ax6.legend(fontsize=8)
ax6.grid(True, alpha=0.3)

plt.tight_layout(rect=[0, 0, 1, 0.93])
out = os.path.join(root, 'compare_02210959_detail.png')
plt.savefig(out, dpi=150)
print(f'\nSaved: {out}')

# Also generate sw1-1 comparison
fig2, axes2 = plt.subplots(1, 2, figsize=(14, 5))
fig2.suptitle('sw1-1: B&R ADINA vs Reference', fontsize=13, fontweight='bold')

sw1 = os.path.join(root, 'sw1-1')
sw1_ref_d = load_safe(os.path.join(sw1, 'ref_53.txt'), 2)
sw1_ref_r = load_safe(os.path.join(sw1, 'ref_1.txt'), 6)
sw1_br_d  = load_safe(os.path.join(sw1, '53.txt'), 2)
sw1_br_r  = load_safe(os.path.join(sw1, '1.txt'), 6)

ax = axes2[0]
if sw1_ref_d is not None and sw1_ref_r is not None:
    n = min(len(sw1_ref_d), len(sw1_ref_r))
    ax.plot(sw1_ref_d[:n, 1]*1000, sw1_ref_r[:n, 1:].sum(axis=1)/1000, 'b-', lw=1.2, label='Ref')
if sw1_br_d is not None and sw1_br_r is not None:
    n = min(len(sw1_br_d), len(sw1_br_r))
    sw1_shear = sw1_br_r[:n, 1:].sum(axis=1)/1000
    sw1_clip = ref_shear_range * 5 if ref_shear_range > 0 else 500
    ax.plot(sw1_br_d[:n, 1]*1000, np.clip(sw1_shear, -sw1_clip, sw1_clip), 'r-', lw=0.8, alpha=0.7, label='B&R')
ax.set_xlabel('Top Disp (mm)')
ax.set_ylabel('Base Shear (kN)')
ax.set_title('sw1-1 Hysteresis')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

ax = axes2[1]
if sw1_ref_d is not None:
    ax.plot(range(len(sw1_ref_d)), sw1_ref_d[:, 1]*1000, 'b-', lw=0.8, label=f'Ref ({len(sw1_ref_d)} pts)')
if sw1_br_d is not None:
    ax.plot(range(len(sw1_br_d)), sw1_br_d[:, 1]*1000, 'r-', lw=0.8, label=f'B&R ({len(sw1_br_d)} pts)')
ax.set_xlabel('Step')
ax.set_ylabel('Disp (mm)')
ax.set_title('sw1-1 Disp History')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

plt.tight_layout()
out2 = os.path.join(root, 'compare_sw1_02210959.png')
plt.savefig(out2, dpi=150)
print(f'Saved: {out2}')

# sw2-1 comparison
fig3, axes3 = plt.subplots(1, 2, figsize=(14, 5))
fig3.suptitle('sw2-1: B&R ADINA vs Reference', fontsize=13, fontweight='bold')

sw2 = os.path.join(root, 'sw2-1')
sw2_ref_d = load_safe(os.path.join(sw2, 'ref_28.txt'), 2)
sw2_ref_r = load_safe(os.path.join(sw2, 'ref_1.txt'), 6)
sw2_br_d  = load_safe(os.path.join(sw2, '28.txt'), 2)
sw2_br_r  = load_safe(os.path.join(sw2, '1.txt'), 6)

ax = axes3[0]
if sw2_ref_d is not None and sw2_ref_r is not None:
    n = min(len(sw2_ref_d), len(sw2_ref_r))
    ax.plot(sw2_ref_d[:n, 1]*1000, sw2_ref_r[:n, 1:].sum(axis=1)/1000, 'b-', lw=1.2, label='Ref')
if sw2_br_d is not None and sw2_br_r is not None:
    n = min(len(sw2_br_d), len(sw2_br_r))
    sw2_shear = sw2_br_r[:n, 1:].sum(axis=1)/1000
    ax.plot(sw2_br_d[:n, 1]*1000, np.clip(sw2_shear, -500, 500), 'r-', lw=0.8, alpha=0.7, label='B&R (clip ±500kN)')
ax.set_xlabel('Top Disp (mm)')
ax.set_ylabel('Base Shear (kN)')
ax.set_title('sw2-1 Hysteresis')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

ax = axes3[1]
if sw2_ref_d is not None:
    ax.plot(range(len(sw2_ref_d)), sw2_ref_d[:, 1]*1000, 'b-', lw=0.8, label=f'Ref ({len(sw2_ref_d)} pts)')
if sw2_br_d is not None:
    ax.plot(range(len(sw2_br_d)), sw2_br_d[:, 1]*1000, 'r-', lw=0.8, label=f'B&R ({len(sw2_br_d)} pts)')
ax.set_xlabel('Step')
ax.set_ylabel('Disp (mm)')
ax.set_title('sw2-1 Disp History')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

plt.tight_layout()
out3 = os.path.join(root, 'compare_sw2_02210959.png')
plt.savefig(out3, dpi=150)
print(f'Saved: {out3}')

# Print all data stats
print('\n=== Complete Data Summary ===')
for name, d, r in [
    ('ML', br_d, br_r), ('ML Ref', ref_d[:nref], ref_r[:nref]),
    ('sw1 BR', sw1_br_d, sw1_br_r), ('sw1 Ref', sw1_ref_d, sw1_ref_r),
    ('sw2 BR', sw2_br_d, sw2_br_r), ('sw2 Ref', sw2_ref_d, sw2_ref_r)]:
    if d is not None and r is not None:
        n = min(len(d), len(r))
        if r.shape[1] >= 6:
            s = r[:n, 1:].sum(axis=1)/1000
        elif r.shape[1] >= 10:
            s = r[:n, 1:].sum(axis=1)/1000
        else:
            s = r[:n, 1]/1000
        print(f'{name:10s}: {n:5d} pts, disp [{d[:n,1].min()*1000:10.1f}, {d[:n,1].max()*1000:8.1f}] mm, '
              f'shear [{s.min():12.0f}, {s.max():10.0f}] kN')
