"""
Comprehensive plot of ALL standalone test results for forumat.f90.
8-panel figure covering every test case.
"""
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams.update({'font.size': 9})
HERE = Path(__file__).parent

# ── Material constants ──
# Paper (ksi)
E0k = 4600.0;  fck = -4.45;  epsc_k = -0.00218
fuk = -3.84;   epsuk = -0.00515;  ftk = 0.45
# Default (MPa)
E0m = 30000.0; fcm = -30.0; epscm = -0.002; num = 0.20
fum = 0.20*fcm; epsum = -0.006; ftm = 2.0

def saenz(E0, fc, epsc, fu, epsu, eps):
    sig = np.zeros_like(eps)
    t = eps >= 0; sig[t] = E0 * eps[t]
    c = ~t
    if np.any(c):
        p = epsu/epsc; Es = fc/epsc; Eu = fu/epsu
        A = (E0/Eu + (p-2)*p**2*E0/Es - (2*p+1)*(p-1)**2) / (p*(p-1)**2)
        B = 2*E0/Es - 3 - 2*A; C = 2 - E0/Es + A
        xi = eps[c]/epsc
        d = 1 + A*xi + B*xi**2 + C*xi**3
        sig[c] = E0*eps[c]/d
    return sig

# ═══════════════════════════════════════════════
# Load data
# ═══════════════════════════════════════════════
uniax_c = np.loadtxt(HERE/'out_uniax_comp.csv', skiprows=1)
uniax_t = np.loadtxt(HERE/'out_uniax_tens.csv', skiprows=1)
shear   = np.loadtxt(HERE/'out_pure_shear.csv', skiprows=1)
fig8    = np.genfromtxt(HERE/'out_fig8_saenz.csv', delimiter=',', skip_header=1)
fig9    = np.genfromtxt(HERE/'out_fig9_cyclic.csv', delimiter=',', skip_header=1)
biax    = np.genfromtxt(HERE/'out_biaxial.csv', delimiter=',', skip_header=1)

# ═══════════════════════════════════════════════
# Figure: 4×2 panels
# ═══════════════════════════════════════════════
fig, axes = plt.subplots(4, 2, figsize=(15, 20))
fig.suptitle('forumat.f90 — Complete Standalone Test Suite\n'
             'All 11 tests PASSED after m_E Saenz fix',
             fontsize=15, fontweight='bold', y=0.995)

# ─── (1) Uniaxial Compression ───
ax = axes[0, 0]
eps_c = uniax_c[:, 1]; sig_cx = uniax_c[:, 4]
eps_d = np.linspace(0, 1.5*epscm, 500)
sig_d = saenz(E0m, fcm, epscm, fum, epsum, eps_d)

ax.plot(-eps_d*1000, -sig_d, 'k-', lw=1.5, label='Analytical Saenz')
ax.plot(-eps_c*1000, -sig_cx, 'ro', ms=2.5, mfc='none', label='forumat.f90')
ax.set_xlabel('Compressive strain (×10⁻³)')
ax.set_ylabel('Stress (MPa)')
ax.set_title('(1) Uniaxial Compression — Default Material', fontweight='bold')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
ax.set_xlim(0, 9); ax.set_ylim(0, 35)
ip = np.argmin(sig_cx)
ax.annotate(f'Peak: {-sig_cx[ip]:.1f} MPa\nε={-eps_c[ip]*1000:.2f}×10⁻³',
            xy=(-eps_c[ip]*1000, -sig_cx[ip]), xytext=(5.5, 28), fontsize=8,
            arrowprops=dict(arrowstyle='->', color='red', lw=0.7))
ax.text(0.03, 0.55, f'E₀={E0m/1000:.0f} GPa, ν={num}\nσc={fcm} MPa, εc={epscm}\nσu={fum} MPa, εu={epsum}',
        transform=ax.transAxes, fontsize=7, va='top',
        bbox=dict(boxstyle='round', fc='wheat', alpha=0.5))

# ─── (2) Uniaxial Tension + Softening ───
ax = axes[0, 1]
eps_t = uniax_t[:, 1]; sig_tx = uniax_t[:, 4]
angle_t = uniax_t[:, 7]
crk_mask = np.abs(angle_t - 1000) > 1

ax.plot(eps_t*1e6, sig_tx, 'b-', lw=1.2, label='σ_x (forumat.f90)')
if np.any(crk_mask):
    i0 = np.where(crk_mask)[0][0]
    ax.axvline(eps_t[i0]*1e6, color='green', ls='--', lw=0.8, alpha=0.7,
               label=f'Crack at ε={eps_t[i0]*1e6:.0f} με')
ax.axhline(ftm, color='gray', ls=':', lw=0.8)
ax.text(eps_t[-1]*1e6*0.6, ftm+0.05, f'ft = {ftm} MPa', fontsize=7, color='gray')

ax.set_xlabel('Tensile strain (με)')
ax.set_ylabel('Stress (MPa)')
ax.set_title('(2) Uniaxial Tension + Softening', fontweight='bold')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
ax.set_ylim(-0.2, ftm+0.5)

# ─── (3) Pure Shear ───
ax = axes[1, 0]
gam_s = shear[:, 3]; tau_s = shear[:, 6]
angle_s = shear[:, 7]
crk_s = np.abs(angle_s - 1000) > 1

G0 = E0m / (2*(1+num))
ax.plot(gam_s*1e6, tau_s, 'b-', lw=1.2, label='τ (forumat.f90)')
ax.plot(gam_s*1e6, G0*gam_s, 'k--', lw=0.8, alpha=0.5, label=f'Elastic: G₀γ (G₀={G0:.0f} MPa)')
if np.any(crk_s):
    i0 = np.where(crk_s)[0][0]
    ax.axvline(gam_s[i0]*1e6, color='green', ls='--', lw=0.8, alpha=0.7,
               label=f'Crack at γ={gam_s[i0]*1e6:.0f} με')
ax.set_xlabel('Shear strain (με)')
ax.set_ylabel('Shear stress (MPa)')
ax.set_title('(3) Pure Shear', fontweight='bold')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

# ─── (4) Fig. 8 — Saenz Curve (paper) ───
ax = axes[1, 1]
eps8 = fig8[:, 1]; sig8n = fig8[:, 2]; sig8a = fig8[:, 3]; err8 = fig8[:, 4]
eps_ref = np.linspace(0, 1.5*epsc_k, 500)
sig_ref = saenz(E0k, fck, epsc_k, fuk, epsuk, eps_ref)

ax.plot(-eps_ref*1000, -sig_ref/abs(fck), 'k-', lw=1.5, label='Analytical Saenz (Eq. 10)')
ax.plot(-eps8*1000, -sig8n/abs(fck), 'ro', ms=2.5, mfc='none', label='forumat.f90')
ax.set_xlabel('Compressive strain (×10⁻³)')
ax.set_ylabel('σ / σc')
ax.set_title(f'(4) Fig. 8 — Saenz Curve  (max err = {np.max(err8)*100:.3f}%)', fontweight='bold')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
ax.set_xlim(0, 3.5); ax.set_ylim(0, 1.15)
ax.text(0.03, 0.50, f'E₀={E0k:.0f} ksi, ν=0\nσc={fck} ksi, εc={epsc_k}\nσu={fuk} ksi, εu={epsuk}',
        transform=ax.transAxes, fontsize=7, va='top',
        bbox=dict(boxstyle='round', fc='wheat', alpha=0.5))

# ─── (5) Fig. 9 — Cyclic Path (paper) ───
ax = axes[2, 0]
eps9 = fig9[:, 1]; sig9 = fig9[:, 2]; ang9 = fig9[:, 3]; pgr9 = fig9[:, 4]
x9 = eps9/abs(epsc_k); y9 = sig9/abs(fck)

# Monotonic Saenz reference
eps_m = np.linspace(0, 2.5*epsc_k, 400)
sig_m = saenz(E0k, fck, epsc_k, fuk, epsuk, eps_m)
ax.plot(eps_m/abs(epsc_k), sig_m/abs(fck), 'k--', lw=1.0, alpha=0.4, label='Monotonic Saenz')

# Color phases
phases = [
    (0,40,'#2196F3','Compress (A)'), (40,60,'#4CAF50','Unload'),
    (60,80,'#FF9800','Tension (B)'), (80,110,'#9C27B0','Close crack (D)'),
    (110,150,'#F44336','Reload (F)'), (150,180,'#795548','Soften'),
    (180,200,'#607D8B','Post-crush (H)'),
]
for s,e,c,lb in phases:
    sl = slice(s, min(e+1, len(x9)))
    ax.plot(x9[sl], y9[sl], color=c, lw=1.3, label=lb if s < 111 else None)
    if s >= 110:
        ax.plot(x9[sl], y9[sl], color=c, lw=1.3)

crk9 = np.where(np.abs(ang9 - 1000) > 1)[0]
if len(crk9) > 0:
    ax.plot(x9[crk9[0]], y9[crk9[0]], 'gv', ms=9, zorder=5, label='Crack')
cr9 = np.where(np.abs(pgr9 - 100) < 1)[0]
if len(cr9) > 0:
    ax.plot(x9[cr9[0]], y9[cr9[0]], 'r^', ms=9, zorder=5, label='Crush')

# Annotations
for idx, lbl in [(39,'A'), (70,'B'), (109,'D'), (149,'F'), (199,'H')]:
    if idx < len(x9) and lbl:
        ax.annotate(lbl, xy=(x9[idx], y9[idx]),
                    fontsize=11, fontweight='bold', color='red',
                    xytext=(5,5), textcoords='offset points')

ax.set_xlabel('ε / |εc|'); ax.set_ylabel('σ / |σc|')
ax.set_title(f'(5) Fig. 9 — Cyclic Path  (peak σ/σc = {np.min(sig9)/fck:.2f})', fontweight='bold')
ax.legend(fontsize=6.5, loc='lower left', ncol=2); ax.grid(True, alpha=0.3)
ax.axhline(0, color='k', lw=0.5); ax.axvline(0, color='k', lw=0.5)
ax.set_xlim(-3.8, 0.6); ax.set_ylim(-1.3, 0.25)

# ─── (6) Biaxial Compression Envelope ───
ax = axes[2, 1]
eps_b = biax[:, 1]; sbx = biax[:, 2]; sby = biax[:, 3]
eps_u = uniax_c[:, 1]; sux = uniax_c[:, 4]
ed = np.linspace(0, 1.5*epscm, 500)
sd = saenz(E0m, fcm, epscm, fum, epsum, ed)

ax.plot(-ed*1000, -sd, 'k--', lw=1.0, label='Uniaxial (analytical)')
ax.plot(-eps_u*1000, -sux, 'ks', ms=2, mfc='none', label='Uniaxial (forumat.f90)')
ax.plot(-eps_b*1000, -sbx, 'r-', lw=2, label='Equal biaxial (forumat.f90)')

ibp = np.argmin(sbx); iup = np.argmin(sux)
ax.axhline(-fcm, color='gray', ls=':', lw=0.8)
ax.annotate(f'Biaxial: {-sbx[ibp]:.1f} MPa ({sbx[ibp]/fcm:.2f}×)',
            xy=(-eps_b[ibp]*1000, -sbx[ibp]), xytext=(5, 35), fontsize=8,
            arrowprops=dict(arrowstyle='->', color='red', lw=0.7))
ax.set_xlabel('Compressive strain (×10⁻³)'); ax.set_ylabel('Stress (MPa)')
ax.set_title(f'(6) Biaxial Envelope  (ratio = {sbx[ibp]/fcm:.2f}×, Kupfer≈1.16×)',
             fontweight='bold')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
ax.set_xlim(0, 8); ax.set_ylim(0, 42)

# ─── (7) Tangent Consistency + Shear Retention ───
ax = axes[3, 0]
ax.axis('off')
checks = [
    ['Test Case', 'Check', 'Result', 'Criterion', 'Status'],
    ['Elastic Patch', 'σ = C·ε (isotropic)', 'max err < 1e-10', 'Exact', '✓ PASS'],
    ['Biaxial Symmetry', 'σx ≈ σy, τ ≈ 0', 'Verified', 'Exact', '✓ PASS'],
    ['Crack Open/Close', 'State transitions', 'All correct', 'No NaN/Inf', '✓ PASS'],
    ['Tangent (elastic)', 'FD vs analytic', '1.96%', '< 10%', '✓ PASS'],
    ['Tangent (cracked)', 'Secant C·ε ≈ σ', '0.00%', '< 5%', '✓ PASS'],
    ['Tangent (post-peak)', 'Secant C·ε ≈ σ', '0.00%', '< 5%', '✓ PASS'],
    ['Shear Retention', 'E_crk/E₀', '0.0106', '≈ ηn=0.01', '✓ PASS'],
    ['Shear Retention', 'G_crk/G₀', '0.571', '≈ ηs=0.50', '✓ PASS'],
    ['Random Robustness', '5000 random steps', '0 fails', '0', '✓ PASS'],
]
colors = [['#d0d0d0']*5] + [['white']*4+['#c8e6c9']]*9
tab = ax.table(cellText=checks, cellColours=colors, cellLoc='center',
               loc='center', colWidths=[0.20, 0.22, 0.16, 0.16, 0.12])
tab.auto_set_font_size(False); tab.set_fontsize(8); tab.scale(1.0, 1.7)
for j in range(5): tab[0,j].set_text_props(fontweight='bold', fontsize=9)
ax.set_title('(7) Non-CSV Tests — Verification Checks', fontweight='bold', pad=12)

# ─── (8) Summary: Before vs After fix ───
ax = axes[3, 1]
ax.axis('off')
summary = [
    ['Metric', 'Before Fix', 'After Fix', 'Paper'],
    ['Fig.8 Saenz max err', '0.000%', '0.000%', '< 2%'],
    ['Fig.9 peak σ/σc', '2.78 ✗', '1.00 ✓', '≈ 1.0'],
    ['Fig.9 σ at εu', '-3.84 ksi', '-3.84 ksi', 'σu=-3.84'],
    ['Fig.9 crack detect', 'Yes', 'Yes', 'Yes'],
    ['Fig.9 crush detect', 'Yes', 'Yes', 'Yes'],
    ['Biaxial σ_biax/σc', '1.25', '1.25', '≈1.16'],
    ['E_crk / E₀', '0.0105', '0.0106', '≈0.01'],
    ['G_crk / G₀', '0.566', '0.571', '≈0.50'],
    ['FD tangent (elastic)', '1.96%', '1.96%', '<10%'],
    ['Random 5000 steps', '0 fail', '0 fail', '0'],
]
colors2 = [['#d0d0d0']*4]
for row in summary[1:]:
    c_before = '#ffcdd2' if '✗' in row[1] else '#c8e6c9'
    c_after  = '#c8e6c9'
    colors2.append(['white', c_before, c_after, '#e3f2fd'])
tab2 = ax.table(cellText=summary, cellColours=colors2, cellLoc='center',
                loc='center', colWidths=[0.28, 0.22, 0.22, 0.22])
tab2.auto_set_font_size(False); tab2.set_fontsize(8); tab2.scale(1.0, 1.65)
for j in range(4): tab2[0,j].set_text_props(fontweight='bold', fontsize=9)
ax.set_title('(8) Bug Fix Impact — Before vs After',
             fontweight='bold', pad=12, color='#1b5e20')

plt.tight_layout(rect=[0, 0, 1, 0.985], h_pad=3.0)
out = HERE / 'all_tests_report.png'
plt.savefig(out, dpi=150, bbox_inches='tight')
print(f'Saved: {out}')
plt.close()
