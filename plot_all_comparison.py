# -*- coding: utf-8 -*-
"""
Comprehensive comparison: OpenSees-02210959 (B&R ADINA model) vs Reference (7-param HEOM2D)
Multi-layer_Shell + sw2-1 + verify uniaxial
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

root = r'e:\Basic\Concrete_Model\ADINA'

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

# ===================================================================
fig = plt.figure(figsize=(18, 16))
fig.suptitle('OpenSees-02210959 (B&R ADINA 37-param)  vs  Reference (HEOM2D 7-param)\n'
             'forumat.f90 unmodified; only test_br.tcl analysis strategy fixed',
             fontsize=14, fontweight='bold')

# --- 1. Multi-layer_Shell Hysteresis ---
ax1 = fig.add_subplot(3, 2, 1)

ml = os.path.join(root, 'Multi-layer_Shell')
ref_d = load_safe(os.path.join(ml, 'ref_disp1.txt'), 2)
ref_r = load_safe(os.path.join(ml, 'ref_shearforce1.txt'), 10)
br_d = load_safe(os.path.join(ml, 'test_disp.txt'), 2)
br_r = load_safe(os.path.join(ml, 'test_react.txt'), 10)

if ref_d is not None and ref_r is not None:
    n = min(len(ref_d), len(ref_r))
    ref_shear = ref_r[:n, 1:].sum(axis=1) / 1000.0
    ref_disp = ref_d[:n, 1] * 1000.0
    ax1.plot(ref_disp, ref_shear, 'b-', lw=1.0, alpha=0.7, label=f'Ref ({n} pts)')

if br_d is not None and br_r is not None:
    n = min(len(br_d), len(br_r))
    br_shear = br_r[:n, 1:].sum(axis=1) / 1000.0
    br_disp_mm = br_d[:n, 1] * 1000.0
    # Clip extreme values for visualization
    clip_val = 500  # kN
    br_shear_clip = np.clip(br_shear, -clip_val, clip_val)
    ax1.plot(br_disp_mm, br_shear_clip, 'r-', lw=1.0, alpha=0.7, label=f'B&R ({n} pts, clipped ±{clip_val}kN)')

ax1.set_xlabel('Top Displacement (mm)')
ax1.set_ylabel('Base Shear (kN)')
ax1.set_title('Multi-layer Shell: Hysteresis')
ax1.legend(fontsize=8)
ax1.grid(True, alpha=0.3)
ax1.axhline(0, color='k', lw=0.5)
ax1.axvline(0, color='k', lw=0.5)

# --- 2. Multi-layer_Shell Displacement History ---
ax2 = fig.add_subplot(3, 2, 2)
if ref_d is not None:
    ax2.plot(range(len(ref_d)), ref_d[:, 1]*1000, 'b-', lw=0.8, alpha=0.7, label=f'Ref ({len(ref_d)} steps)')
if br_d is not None:
    ax2.plot(range(len(br_d)), br_d[:, 1]*1000, 'r-', lw=1.0, label=f'B&R ({len(br_d)} steps)')
ax2.set_xlabel('Step Number')
ax2.set_ylabel('Top Disp (mm)')
ax2.set_title('Multi-layer Shell: Displacement History')
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.3)

# --- 3. sw2-1 Hysteresis ---
ax3 = fig.add_subplot(3, 2, 3)
sw2 = os.path.join(root, 'sw2-1')
sw2_ref_d = load_safe(os.path.join(sw2, 'ref_28.txt'), 2)
sw2_ref_r = load_safe(os.path.join(sw2, 'ref_1.txt'), 6)
sw2_br_d = load_safe(os.path.join(sw2, '28.txt'), 2)
sw2_br_r = load_safe(os.path.join(sw2, '1.txt'), 6)

if sw2_ref_d is not None and sw2_ref_r is not None:
    n = min(len(sw2_ref_d), len(sw2_ref_r))
    sw2_ref_shear = sw2_ref_r[:n, 1:].sum(axis=1) / 1000.0
    sw2_ref_disp = sw2_ref_d[:n, 1] * 1000.0
    ax3.plot(sw2_ref_disp, sw2_ref_shear, 'b-', lw=1.0, alpha=0.7, label=f'Ref ({n} pts)')

if sw2_br_d is not None and sw2_br_r is not None:
    n = min(len(sw2_br_d), len(sw2_br_r))
    sw2_br_shear = sw2_br_r[:n, 1:].sum(axis=1) / 1000.0
    sw2_br_disp = sw2_br_d[:n, 1] * 1000.0
    ax3.plot(sw2_br_disp, sw2_br_shear, 'r-', lw=1.0, alpha=0.7, label=f'B&R ({n} pts)')

ax3.set_xlabel('Top Displacement (mm)')
ax3.set_ylabel('Base Shear (kN)')
ax3.set_title('sw2-1: Hysteresis')
ax3.legend(fontsize=8)
ax3.grid(True, alpha=0.3)
ax3.axhline(0, color='k', lw=0.5)

# --- 4. sw2-1 Displacement History ---
ax4 = fig.add_subplot(3, 2, 4)
sw2_shuju = load_safe(os.path.join(sw2, 'shuju2.txt'), 1)
if sw2_shuju is not None:
    ax4.plot(np.arange(len(sw2_shuju))*0.1, sw2_shuju[:, 0]*1000, 'g--', lw=0.8, alpha=0.5, label='Protocol')
if sw2_ref_d is not None:
    ax4.plot(sw2_ref_d[:, 0], sw2_ref_d[:, 1]*1000, 'b-', lw=0.8, alpha=0.7, label=f'Ref ({len(sw2_ref_d)} pts)')
if sw2_br_d is not None:
    ax4.plot(sw2_br_d[:, 0], sw2_br_d[:, 1]*1000, 'r-', lw=1.0, label=f'B&R ({len(sw2_br_d)} pts)')
ax4.set_xlabel('Time')
ax4.set_ylabel('Top Disp (mm)')
ax4.set_title('sw2-1: Displacement History')
ax4.legend(fontsize=8)
ax4.grid(True, alpha=0.3)

# --- 5. Uniaxial Compression ---
ax5 = fig.add_subplot(3, 2, 5)
vdir = os.path.join(root, 'verify')
comp_ad = load_safe(os.path.join(vdir, 'vuniax_adina_comp_disp.txt'), 3)
comp_ar = load_safe(os.path.join(vdir, 'vuniax_adina_comp_react.txt'), 3)
comp_rd = load_safe(os.path.join(vdir, 'vuniax_ref_comp_disp.txt'), 3)
comp_rr = load_safe(os.path.join(vdir, 'vuniax_ref_comp_react.txt'), 3)

if comp_rd is not None and comp_rr is not None:
    ax5.plot(comp_rd[:, 1]*1000, comp_rr[:, 1]/1000, 'b-', lw=1.5, label='Ref')
if comp_ad is not None and comp_ar is not None:
    ax5.plot(comp_ad[:, 1]*1000, comp_ar[:, 1]/1000, 'r--', lw=1.5, label='B&R')
ax5.set_xlabel('Displacement (mm)')
ax5.set_ylabel('Reaction (kN)')
ax5.set_title('Uniaxial Compression')
ax5.legend(fontsize=8)
ax5.grid(True, alpha=0.3)

# --- 6. Uniaxial Cyclic ---
ax6 = fig.add_subplot(3, 2, 6)
cyc_ad = load_safe(os.path.join(vdir, 'vuniax_adina_cyc_disp.txt'), 3)
cyc_ar = load_safe(os.path.join(vdir, 'vuniax_adina_cyc_react.txt'), 3)
cyc_rd = load_safe(os.path.join(vdir, 'vuniax_ref_cyc_disp.txt'), 3)
cyc_rr = load_safe(os.path.join(vdir, 'vuniax_ref_cyc_react.txt'), 3)

if cyc_rd is not None and cyc_rr is not None:
    ax6.plot(cyc_rd[:, 1]*1000, cyc_rr[:, 1]/1000, 'b-', lw=1.0, alpha=0.7, label='Ref')
if cyc_ad is not None and cyc_ar is not None:
    ax6.plot(cyc_ad[:, 1]*1000, cyc_ar[:, 1]/1000, 'r-', lw=1.0, alpha=0.7, label='B&R')
ax6.set_xlabel('Displacement (mm)')
ax6.set_ylabel('Reaction (kN)')
ax6.set_title('Uniaxial Cyclic')
ax6.legend(fontsize=8)
ax6.grid(True, alpha=0.3)

plt.tight_layout(rect=[0, 0, 1, 0.94])
out = os.path.join(root, 'compare_all_02210959.png')
plt.savefig(out, dpi=150)
print(f'Saved: {out}')

# Print summary
print('\n=== Data Summary ===')
if br_d is not None:
    print(f'ML B&R: {len(br_d)} disp pts, disp range [{br_d[:,1].min()*1000:.1f}, {br_d[:,1].max()*1000:.1f}] mm')
if br_r is not None:
    br_s = br_r[:, 1:].sum(axis=1)/1000
    print(f'ML B&R: {len(br_r)} react pts, shear range [{br_s.min():.0f}, {br_s.max():.0f}] kN')
if ref_d is not None:
    print(f'ML Ref: {len(ref_d)} disp pts, disp range [{ref_d[:,1].min()*1000:.1f}, {ref_d[:,1].max()*1000:.1f}] mm')
if sw2_br_d is not None:
    print(f'sw2 B&R: {len(sw2_br_d)} pts, disp range [{sw2_br_d[:,1].min()*1000:.1f}, {sw2_br_d[:,1].max()*1000:.1f}] mm')
