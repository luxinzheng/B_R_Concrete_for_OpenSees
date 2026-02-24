"""
plot_biaxial_comparison.py
Compare biaxial test results with:
  1) Paper's theoretical prediction (Poisson coupling + envelope)
  2) Kupfer et al. (1969) experimental data
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os

rcParams['font.size'] = 10
rcParams['figure.dpi'] = 150

csv_path = os.path.join(os.path.dirname(__file__), 'out_biaxial_paths.csv')
data = np.genfromtxt(csv_path, delimiter=',', skip_header=1)

fc = -30.0
E0 = 30000.0
nu = 0.2
epsc = -0.002
ft = 2.5

path_ids = sorted(set(data[:, 0].astype(int)))
labels = {
    1: r'$-1:-1$', 2: r'$-1:-0.5$', 3: r'$-1:-0.2$', 4: r'$-1:0$',
    5: r'$-1:0.1$', 6: r'$-1:0.2$', 7: r'$-1:0.5$', 8: r'$1:1$',
}
sig_rat = {
    1: (-1,-1), 2: (-1,-0.5), 3: (-1,-0.2), 4: (-1,0),
    5: (-1,0.1), 6: (-1,0.2), 7: (-1,0.5), 8: (1,1),
}

cmap = plt.cm.tab10
colors = {pid: cmap(i / max(len(path_ids)-1, 1)) for i, pid in enumerate(path_ids)}


def find_peak(pid, data):
    """Find peak stress state for a given path."""
    mask = data[:, 0].astype(int) == pid
    s1 = data[mask, 4]
    s2 = data[mask, 5]
    e1 = data[mask, 2]
    if pid <= 7:
        idx = np.argmin(s1)
    else:
        idx = np.argmax(s1)
    return s1[idx], s2[idx], e1[idx]


def theoretical_biaxial_peak(r1, r2, fc, E0, nu, epsc):
    """
    Predict biaxial peak stress using the code's total-strain secant approach.
    For proportional strain loading with stress ratio r1:r2.
    Biaxial enhancement comes from Poisson coupling only (uniaxial Saenz curve).
    """
    if r1 < 0:
        if abs(r1) >= abs(r2):
            alpha = (r2/r1 - nu) / (1.0 - (r2/r1)*nu)
            eff = 1.0 + nu * alpha
            sig_peak = fc * eff / (1.0 - nu**2)  * (epsc) / (E0 / (1-nu**2) * epsc / fc)
            peak1 = fc / (1 - nu * alpha)
            peak2 = alpha * fc / (1 - nu * alpha) if abs(alpha) > 1e-10 else 0
        else:
            alpha = (r1/r2 - nu) / (1.0 - (r1/r2)*nu)
            peak2 = fc / (1 - nu * alpha)
            peak1 = alpha * fc / (1 - nu * alpha) if abs(alpha) > 1e-10 else 0
    else:
        peak1 = ft
        peak2 = ft
    return peak1, peak2


def saenz_stress(eps, E0, fc, epsc):
    """Uniaxial Saenz stress for a given strain (compression negative)."""
    if eps >= 0:
        return E0 * eps
    Es = fc / epsc
    eu = -0.006
    sig_u = 0.2 * fc
    Eu = sig_u / eu
    p = eu / epsc
    A = (E0/Es + (p**2 - 2*p) * E0/Es - (2*p**3 - 3*p**2 + 1))
    B = ((E0/Es - 3) - 2*A) * (-1)
    # Use simplified Saenz: σ/σc = (E0/Es)(ε/εc) / [1 + A(ε/εc) + B(ε/εc)^2 + C(ε/εc)^3]
    r = eps / epsc
    A2 = E0/Es - 2.0
    B2 = 1.0 - 2.0 * E0/Es
    C2 = E0/Es
    sig = fc * (E0/Es)*r / (1 + A2*r + B2*r**2 + C2*r**3) if abs(1 + A2*r + B2*r**2 + C2*r**3) > 1e-20 else fc
    return sig


# ---- Compute biaxial peaks from code results ----
code_peaks = {}
for pid in path_ids:
    s1p, s2p, e1p = find_peak(pid, data)
    code_peaks[pid] = (s1p/fc, s2p/fc)

# ---- Theoretical prediction: secant approach + Poisson coupling ----
# For proportional strain: ε₂ = α·ε₁
# σ₁ = E_sec/(1-ν²)(ε₁ + ν·α·ε₁) = E_sec·ε₁·(1+ν·α)/(1-ν²)
# At Saenz peak (ε₁ = εc for the dominant direction):
#   E_sec = fc/εc = Es
#   σ₁ = Es·εc·(1+ν·α)/(1-ν²) = fc·(1+ν·α)/(1-ν²)

theory_s1fc = []
theory_s2fc = []
for pid in path_ids:
    r1, r2 = sig_rat[pid]
    if r1 < 0:
        if abs(r1) >= abs(r2):
            r = r2 / r1
            alpha = (r - nu) / (1.0 - r * nu)
            s1_peak_fc = (1 + nu*alpha) / (1 - nu**2)
            s2_peak_fc = (alpha + nu) / (1 - nu**2)
        else:
            r = r1 / r2
            alpha = (r - nu) / (1.0 - r * nu)
            s2_peak_fc = (1 + nu*alpha) / (1 - nu**2)
            s1_peak_fc = (alpha + nu) / (1 - nu**2)
    else:
        s1_peak_fc = -ft/fc
        s2_peak_fc = -ft/fc
    theory_s1fc.append(s1_peak_fc)
    theory_s2fc.append(s2_peak_fc)

# ---- Kupfer et al. (1969) experimental biaxial data (normalized σ/fc) ----
# Compression-compression quadrant (approximate digitized data)
kupfer_s1fc = np.array([1.0,  1.10, 1.175, 1.25, 1.275, 1.25, 1.20, 1.16])
kupfer_s2fc = np.array([0.0,  0.20, 0.40,  0.52, 0.70,  0.90, 1.05, 1.16])
kupfer_t_s1 = np.array([-0.65, -0.55, -0.45, -0.35, -0.20, 0.0])
kupfer_t_s2 = np.array([0.0, 0.05, 0.10, 0.15, 0.25, -ft/fc])

# ======================================================================
fig, axes = plt.subplots(2, 2, figsize=(15, 12))
fig.suptitle('Biaxial Strength: Code Results vs Paper Model vs Kupfer Data',
             fontsize=14, fontweight='bold')

# ---- Panel (a): σ₁/fc vs ε₁ stress-strain curves ----
ax1 = axes[0, 0]
for pid in path_ids:
    mask = data[:, 0].astype(int) == pid
    eps1 = data[mask, 2]
    sig1_fc = data[mask, 7]
    ax1.plot(eps1 * 1000, sig1_fc, color=colors[pid], label=labels[pid], lw=1.2)
ax1.set_xlabel(r'$\varepsilon_1$ (‰)')
ax1.set_ylabel(r'$\sigma_1 / f_c$')
ax1.set_title(r'(a) $\sigma_1/f_c$ vs $\varepsilon_1$')
ax1.legend(fontsize=8, ncol=2)
ax1.grid(True, alpha=0.3)
ax1.axhline(0, color='k', lw=0.5)

# ---- Panel (b): σ₂/fc vs ε₂ stress-strain curves ----
ax2 = axes[0, 1]
for pid in path_ids:
    mask = data[:, 0].astype(int) == pid
    eps2 = data[mask, 3]
    sig2_fc = data[mask, 8]
    ax2.plot(eps2 * 1000, sig2_fc, color=colors[pid], label=labels[pid], lw=1.2)
ax2.set_xlabel(r'$\varepsilon_2$ (‰)')
ax2.set_ylabel(r'$\sigma_2 / f_c$')
ax2.set_title(r'(b) $\sigma_2/f_c$ vs $\varepsilon_2$')
ax2.legend(fontsize=8, ncol=2)
ax2.grid(True, alpha=0.3)
ax2.axhline(0, color='k', lw=0.5)

# ---- Panel (c): Biaxial Failure Envelope Comparison ----
ax3 = axes[1, 0]

# Kupfer compression-compression
ax3.plot(kupfer_s1fc, kupfer_s2fc, 'ks-', markersize=5, lw=1.5,
         label='Kupfer et al. (1969)', zorder=10)
# Mirror for symmetry
ax3.plot(kupfer_s2fc, kupfer_s1fc, 'ks-', markersize=5, lw=1.5, zorder=10)

# Code results (peak points)
code_x = [code_peaks[pid][0] for pid in path_ids]
code_y = [code_peaks[pid][1] for pid in path_ids]
ax3.plot(code_x, code_y, 'ro', markersize=8, label='Code (this test)', zorder=15)
for pid in path_ids:
    ax3.annotate(labels[pid], (code_peaks[pid][0], code_peaks[pid][1]),
                 textcoords='offset points', xytext=(6, 6), fontsize=7, color='red')

# Theoretical (Poisson-only)
ax3.plot(theory_s1fc, theory_s2fc, 'b^', markersize=7,
         label=r'Theory: $1/(1-\nu\cdot\alpha)$', zorder=12)

# Uniaxial reference lines
ax3.axhline(1.0, color='gray', ls=':', lw=0.5, alpha=0.5)
ax3.axvline(1.0, color='gray', ls=':', lw=0.5, alpha=0.5)
ax3.plot([0, 1.5], [0, 1.5], ':', color='gray', lw=0.5, alpha=0.3)

ax3.set_xlabel(r'$\sigma_1 / f_c$')
ax3.set_ylabel(r'$\sigma_2 / f_c$')
ax3.set_title('(c) Biaxial Failure Envelope Comparison')
ax3.legend(fontsize=8, loc='upper left')
ax3.grid(True, alpha=0.3)
ax3.set_xlim(-0.15, 1.5)
ax3.set_ylim(-0.15, 1.5)
ax3.set_aspect('equal')

# ---- Panel (d): Detailed comparison table as text ----
ax4 = axes[1, 1]
ax4.axis('off')

header = f'{"Path":<10} {"Code σ₁/fc":>11} {"Code σ₂/fc":>11} {"Theory σ₁/fc":>13} {"Kupfer σ₁/fc":>13}'
lines = [header, '─' * 60]

kupfer_ref = {
    1: 1.16, 2: 1.20, 3: 1.05, 4: 1.00,
    5: None, 6: None, 7: None, 8: None
}
kupfer_ref2 = {
    1: 1.16, 2: 0.60, 3: 0.21, 4: 0.00,
}

for i, pid in enumerate(path_ids):
    c1, c2 = code_peaks[pid]
    t1, t2 = theory_s1fc[i], theory_s2fc[i]
    kr = kupfer_ref.get(pid)
    kr_str = f'{kr:>13.4f}' if kr is not None else f'{"N/A":>13}'
    lines.append(f'{labels[pid]:<10} {c1:>11.4f} {c2:>11.4f} {t1:>13.4f} {kr_str}')

lines.append('')
lines.append('Key Observations:')
lines.append(f'  Equal biaxial (-1:-1): Code={code_peaks[1][0]:.3f}·fc')
lines.append(f'    Theory (Poisson only): {1/(1-nu):.3f}·fc')
lines.append(f'    Kupfer experimental:   ~1.16·fc')
lines.append(f'    Overestimation: {(1/(1-nu) - 1.16)/1.16*100:.1f}%')
lines.append('')
lines.append(f'  Poisson coupling: 1/(1-ν) = 1/{1-nu:.1f} = {1/(1-nu):.3f}')
lines.append(f'  Biaxial envelope (SP tables): NOT used for')
lines.append(f'    stress computation in Branch C (uncracked).')
lines.append(f'    Only used for failure checking.')

txt = '\n'.join(lines)
ax4.text(0.02, 0.98, txt, transform=ax4.transAxes, fontsize=8,
         verticalalignment='top', fontfamily='monospace',
         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
ax4.set_title('(d) Quantitative Comparison')

plt.tight_layout()
outfile = os.path.join(os.path.dirname(__file__), 'biaxial_comparison_report.png')
plt.savefig(outfile, dpi=200, bbox_inches='tight')
print(f'Saved: {outfile}')
plt.close()

# ======================================================================
# Print detailed analysis
# ======================================================================
print('\n' + '='*70)
print('BIAXIAL STRENGTH COMPARISON: Code vs Paper vs Experiment')
print('='*70)

print(f'\n--- Material: E0={E0} MPa, ν={nu}, fc={fc} MPa, εc={epsc} ---\n')

print('{:<12} {:>11} {:>10} {:>10} {:>12}'.format('Path', 'Code s1/fc', 'Theory', 'Kupfer', 'Code-Kupfer'))
print('-' * 58)
for i, pid in enumerate(path_ids):
    c1 = code_peaks[pid][0]
    t1 = theory_s1fc[i]
    kr = kupfer_ref.get(pid)
    if kr is not None:
        diff = '{:>+10.1f}%'.format((c1 - kr)/kr*100)
    else:
        diff = '{:>11}'.format('---')
    kr_s = '{:>10.3f}'.format(kr) if kr is not None else '{:>10}'.format('N/A')
    print('{:<12} {:>11.4f} {:>10.4f} {} {}'.format(labels[pid], c1, t1, kr_s, diff))

print('\n--- Analysis ---')
print('  Biaxial enhancement comes SOLELY from Poisson coupling:')
print('    sig1 = E_sec * eps1 * (1+nu*alpha) / (1-nu^2)')
print('')
print('  Equal biaxial (alpha=1): sig1/fc = 1/(1-nu) = {:.4f}'.format(1/(1-nu)))
print('  Kupfer experiment: sig1/fc ~ 1.16')
print('  Overestimation: {:.1f}%'.format((1/(1-nu) - 1.16)/1.16*100))
print('')
print('  Biaxial envelope (SP tables) NOT used for Saenz curve in Branch C.')
print('  For plane stress, P1=0 -> sp1=0 -> SP31I=SP32I=SP33I=1.0.')
print('  This matches paper for plane stress with default tables.')
print('')
print('  Paper uses same Poisson coupling -> same overestimation in 2D.')
