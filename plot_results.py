"""
Plot ADINA Concrete Plane Stress Model - Uniaxial Test Results
Highlights: Saenz envelope softening + secant reload toward historical max
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.size'] = 11
import matplotlib.pyplot as plt

# ============================================================
# Read CSV data
# ============================================================
def read_csv(filename, has_step=False):
    eps_xx, sig_xx, eps_yy = [], [], []
    with open(filename, 'r') as f:
        f.readline()
        for line in f:
            line = line.strip()
            if not line:
                continue
            if has_step:
                parts = line.split(',')
                vals = parts[1].split()
            else:
                vals = line.split()
            eps_xx.append(float(vals[0]))
            sig_xx.append(float(vals[1]))
            eps_yy.append(float(vals[2]))
    return np.array(eps_xx), np.array(sig_xx), np.array(eps_yy)

e1, s1, ey1 = read_csv('test1_compression.csv')
e2, s2, ey2 = read_csv('test2_tension.csv')
e3, s3, ey3 = read_csv('test3_cyclic.csv', has_step=True)

# Saenz curve reference (for overlay)
def saenz_curve(eps, E0=30000, fc=-30, eps_c=-0.002, fu=-6, eps_u=-0.005):
    """Compute Saenz stress-strain curve"""
    RP = eps_u / eps_c
    ES = fc / eps_c
    EU = fu / eps_u
    A = (E0/EU + (RP-2)*RP**2*E0/ES - (2*RP+1)*(RP-1)**2) / (RP*(RP-1)**2)
    B = 2*E0/ES - 3 - 2*A
    C = 2 - E0/ES + A
    
    sigma = np.zeros_like(eps)
    for i, e in enumerate(eps):
        DE = e / eps_c
        if DE <= 0:
            sigma[i] = E0 * e
        else:
            denom = 1 + A*DE + B*DE**2 + C*DE**3
            sigma[i] = E0 * e / denom if abs(denom) > 1e-30 else 0
    return sigma

eps_ref = np.linspace(0, -0.006, 500)
sig_ref = saenz_curve(eps_ref)

e1_m = e1 * 1e3
e2_m = e2 * 1e3
e3_m = e3 * 1e3
eps_ref_m = eps_ref * 1e3

# ============================================================
# Figure 1: Monotonic Tests
# ============================================================
fig1, axes1 = plt.subplots(1, 2, figsize=(13, 5.5))

# --- Test 1: Monotonic Compression ---
ax = axes1[0]
ax.plot(e1_m, s1, 'b-', linewidth=2.0, label='Model response', zorder=3)
ax.plot(eps_ref_m, sig_ref, 'r--', linewidth=1.2, alpha=0.7, label='Saenz envelope', zorder=2)
ax.axhline(y=-30, color='gray', linestyle=':', linewidth=0.6, alpha=0.5)
ax.axhline(y=-6, color='gray', linestyle=':', linewidth=0.6, alpha=0.5)
ax.axvline(x=-2.0, color='gray', linestyle=':', linewidth=0.6, alpha=0.5)
ax.axvline(x=-5.0, color='gray', linestyle=':', linewidth=0.6, alpha=0.5)
ax.set_xlabel('$\\varepsilon_{xx}$ ($\\times 10^{-3}$)')
ax.set_ylabel('$\\sigma_{xx}$ (MPa)')
ax.set_title('Test 1: Uniaxial Monotonic Compression')
ax.legend(loc='lower left', fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_xlim([-5.5, 0.2])
ax.set_ylim([-35, 2])

# Annotations
ax.annotate('$f_c = -30$ MPa', xy=(-2.0, -30), xytext=(-3.8, -33),
            fontsize=8, color='red', arrowprops=dict(arrowstyle='->', color='red', lw=0.8))
ax.annotate('$f_u = -6$ MPa', xy=(-5.0, -6), xytext=(-3.5, -3),
            fontsize=8, color='red', arrowprops=dict(arrowstyle='->', color='red', lw=0.8))
ax.annotate('Softening', xy=(-3.5, -15), fontsize=9, color='blue',
            ha='center', style='italic')

# --- Test 2: Monotonic Tension ---
ax = axes1[1]
ax.plot(e2_m, s2, 'r-', linewidth=2.0, label='$\\sigma_{xx}$')
ax.axhline(y=3.0, color='gray', linestyle='--', linewidth=0.8, alpha=0.7, label='$f_t = 3$ MPa')
ax.set_xlabel('$\\varepsilon_{xx}$ ($\\times 10^{-3}$)')
ax.set_ylabel('$\\sigma_{xx}$ (MPa)')
ax.set_title('Test 2: Uniaxial Monotonic Tension')
ax.legend(loc='upper left', fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_xlim([-0.02, 1.05])
ax.set_ylim([-0.3, 4.0])

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
ax_inset = inset_axes(ax, width="45%", height="45%", loc='center right')
mask = e2_m < 0.15
ax_inset.plot(e2_m[mask], s2[mask], 'r-', linewidth=1.5)
ax_inset.axhline(y=3.0, color='gray', linestyle='--', linewidth=0.7, alpha=0.7)
ax_inset.set_xlim([-0.002, 0.14])
ax_inset.set_ylim([-0.2, 3.5])
ax_inset.set_title('Cracking detail', fontsize=8)
ax_inset.tick_params(labelsize=7)
ax_inset.grid(True, alpha=0.2)

fig1.tight_layout()
fig1.savefig('plot_monotonic.png', dpi=150, bbox_inches='tight')
print('Saved: plot_monotonic.png')

# ============================================================
# Figure 2: Cyclic Test — Stress-Strain with envelope
# ============================================================
fig2, axes2 = plt.subplots(1, 2, figsize=(14, 6))

phases = [
    (0, 100, 'Phase 1: Tension', '#1f77b4'),
    (100, 200, 'Phase 2: Unload', '#ff7f0e'),
    (200, 300, 'Phase 3: Compression', '#2ca02c'),
    (300, 400, 'Phase 4: Unload', '#d62728'),
    (400, 500, 'Phase 5: Re-tension', '#9467bd'),
    (500, 600, 'Phase 6: Re-compression', '#8c564b'),
]

# --- Left: Stress-Strain with envelope ---
ax = axes2[0]

# Saenz envelope
ax.plot(eps_ref_m, sig_ref, 'k--', linewidth=1.0, alpha=0.4, label='Saenz envelope', zorder=1)

for i0, i1, label, color in phases:
    idx = slice(i0, min(i1+1, len(e3)))
    ax.plot(e3_m[idx], s3[idx], color=color, linewidth=1.8, label=label, zorder=3)

ax.axhline(y=0, color='k', linestyle='-', linewidth=0.3, alpha=0.5)
ax.axvline(x=0, color='k', linestyle='-', linewidth=0.3, alpha=0.5)

ax.set_xlabel('$\\varepsilon_{xx}$ ($\\times 10^{-3}$)')
ax.set_ylabel('$\\sigma_{xx}$ (MPa)')
ax.set_title('Test 3: Cyclic — Stress vs Strain')
ax.legend(loc='lower left', fontsize=7.5, ncol=2)
ax.grid(True, alpha=0.3)
ax.set_xlim([-4.5, 0.8])
ax.set_ylim([-34, 4])

# Annotate secant reloading
# Phase 6 secant line from origin to (eps_max, sig_max)
eps_max_ph3 = -3.0  # x1e-3
sig_max_ph3 = -20.3
ax.plot([0, eps_max_ph3], [0, sig_max_ph3], 'k:', linewidth=0.8, alpha=0.5)
ax.annotate('Secant\nreload', xy=(-1.5, -10.2), fontsize=8, color='#8c564b',
            ha='center', style='italic')
ax.annotate('Softening\n(envelope)', xy=(-3.5, -15), fontsize=8, color='#8c564b',
            ha='center', style='italic')

# --- Right: Stress history ---
ax = axes2[1]
steps = np.arange(len(s3))

for i0, i1, label, color in phases:
    idx = slice(i0, min(i1+1, len(s3)))
    ax.plot(steps[idx], s3[idx], color=color, linewidth=1.8, label=label)

for i0, _, _, _ in phases:
    if i0 > 0:
        ax.axvline(x=i0, color='gray', linestyle=':', linewidth=0.5, alpha=0.5)

ax.axhline(y=-30, color='r', linestyle='--', linewidth=0.7, alpha=0.4, label='$f_c$')
ax.axhline(y=0, color='k', linestyle='-', linewidth=0.3, alpha=0.5)

ax.set_xlabel('Load Step')
ax.set_ylabel('$\\sigma_{xx}$ (MPa)')
ax.set_title('Test 3: Cyclic — Stress History')
ax.legend(loc='lower left', fontsize=7.5, ncol=2)
ax.grid(True, alpha=0.3)
ax.set_ylim([-34, 4])

phase_labels = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6']
phase_centers = [50, 150, 250, 350, 450, 550]
for lbl, xc in zip(phase_labels, phase_centers):
    ax.text(xc, 2, lbl, ha='center', fontsize=8, color='gray', fontweight='bold')

fig2.tight_layout()
fig2.savefig('plot_cyclic.png', dpi=150, bbox_inches='tight')
print('Saved: plot_cyclic.png')

# ============================================================
# Figure 3: Detail comparison — secant reload vs initial loading
# ============================================================
fig3, axes3 = plt.subplots(1, 2, figsize=(13, 5.5))

# --- Left: Overlay Phase 3 and Phase 6 compression paths ---
ax = axes3[0]

# Saenz envelope
ax.plot(eps_ref_m, sig_ref, 'k--', linewidth=1.0, alpha=0.4, label='Saenz envelope')

# Phase 3: initial compression
idx3 = slice(200, 301)
ax.plot(e3_m[idx3], s3[idx3], '#2ca02c', linewidth=2.0, label='Phase 3: 1st compression')

# Phase 6: secant reload + envelope continuation
idx6 = slice(500, 601)
ax.plot(e3_m[idx6], s3[idx6], '#8c564b', linewidth=2.0, label='Phase 6: Re-compression')

# Phase 4: unloading
idx4 = slice(300, 401)
ax.plot(e3_m[idx4], s3[idx4], '#d62728', linewidth=1.5, label='Phase 4: Unload', alpha=0.7)

# Secant line
ax.plot([0, -3.0], [0, -20.3], 'k:', linewidth=1.0, alpha=0.5, label='Secant line')

# E0 reference
ax.plot([0, -1.0], [0, -30.0], 'b:', linewidth=0.8, alpha=0.3, label='$E_0$ slope')

# Mark historical max point
ax.plot(-3.0, -20.3, 'ko', markersize=6, zorder=5)
ax.annotate('Historical max\n$(\\varepsilon_{max}, \\sigma_{env})$',
            xy=(-3.0, -20.3), xytext=(-1.5, -24),
            fontsize=8, arrowprops=dict(arrowstyle='->', lw=0.8))

ax.set_xlabel('$\\varepsilon_{xx}$ ($\\times 10^{-3}$)')
ax.set_ylabel('$\\sigma_{xx}$ (MPa)')
ax.set_title('Compression: Initial vs Reload Path')
ax.legend(loc='lower left', fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_xlim([-4.5, 0.5])
ax.set_ylim([-34, 2])

# --- Right: Secant modulus vs strain ---
ax = axes3[1]

# Compute secant modulus along the Saenz curve
eps_sec = np.linspace(-0.0005, -0.005, 200)
sig_sec = saenz_curve(eps_sec)
E_sec = sig_sec / eps_sec  # secant modulus

ax.plot(eps_sec*1e3, E_sec/1e3, 'b-', linewidth=1.8, label='Secant modulus $E_{sec}$')
ax.axhline(y=30, color='gray', linestyle='--', linewidth=0.8, alpha=0.5, label='$E_0 = 30$ GPa')
ax.axhline(y=15, color='r', linestyle=':', linewidth=0.8, alpha=0.5, label='$E_{sec}(\\varepsilon_c) = 15$ GPa')

# Mark the Phase 3 end point
E_sec_ph3 = -20.3 / -0.003
ax.plot(-3.0, E_sec_ph3/1e3, 'ko', markersize=6, zorder=5)
ax.annotate(f'$E_{{sec}}$ at Phase 3 end\n= {E_sec_ph3:.0f} MPa',
            xy=(-3.0, E_sec_ph3/1e3), xytext=(-1.5, 12),
            fontsize=8, arrowprops=dict(arrowstyle='->', lw=0.8))

ax.set_xlabel('$\\varepsilon_{xx}$ ($\\times 10^{-3}$)')
ax.set_ylabel('Secant Modulus (GPa)')
ax.set_title('Secant Modulus Degradation')
ax.legend(loc='upper right', fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_xlim([-5.5, 0])
ax.set_ylim([0, 35])

fig3.tight_layout()
fig3.savefig('plot_secant_detail.png', dpi=150, bbox_inches='tight')
print('Saved: plot_secant_detail.png')

# ============================================================
# Figure 4: Poisson effect
# ============================================================
fig4, axes4 = plt.subplots(1, 2, figsize=(13, 5))

ax = axes4[0]
ax.plot(e1_m, ey1*1e3, 'g-', linewidth=1.5, label='$\\varepsilon_{yy}$')
ax.plot(e1_m, -0.2*e1_m, 'k--', linewidth=0.8, alpha=0.5, label='$-\\nu \\varepsilon_{xx}$ (elastic)')
ax.set_xlabel('$\\varepsilon_{xx}$ ($\\times 10^{-3}$)')
ax.set_ylabel('$\\varepsilon_{yy}$ ($\\times 10^{-3}$)')
ax.set_title('Test 1: Poisson Effect (Compression)')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

ax = axes4[1]
for i0, i1, label, color in phases:
    idx = slice(i0, min(i1+1, len(ey3)))
    ax.plot(steps[idx], ey3[idx]*1e3, color=color, linewidth=1.5, label=label)
for i0, _, _, _ in phases:
    if i0 > 0:
        ax.axvline(x=i0, color='gray', linestyle=':', linewidth=0.5, alpha=0.5)
ax.set_xlabel('Load Step')
ax.set_ylabel('$\\varepsilon_{yy}$ ($\\times 10^{-3}$)')
ax.set_title('Test 3: Transverse Strain History')
ax.legend(loc='upper left', fontsize=8, ncol=2)
ax.grid(True, alpha=0.3)

fig4.tight_layout()
fig4.savefig('plot_poisson.png', dpi=150, bbox_inches='tight')
print('Saved: plot_poisson.png')

print('\nAll plots generated successfully.')
