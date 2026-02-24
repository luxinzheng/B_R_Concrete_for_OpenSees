"""
Biaxial Failure Envelope for ADINA Concrete Model

Computes the failure surface by:
  C-C region: step-by-step simulation, extract peak principal stress
  C-T region: FALSTR = ft * (1 - sigma_comp/fc) + Saenz peak
  T-T region: analytical (max principal = ft)

Plots sigma_1/|fc| vs sigma_2/|fc| (in-plane principal stresses).
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.size'] = 11
import matplotlib.pyplot as plt

# ============================================================
# Material parameters
# ============================================================
E0    = 30000.0
nu    = 0.2
fc    = -30.0
eps_c = -0.002
fu    = -6.0
eps_u = -0.005
ft    = 3.0
BETA  = 0.75
GAMA  = 1.0
STIFAC = 1e-4

SP1_T  = np.array([0.0, 0.25, 0.5, 0.75, 1.0, 1.2])
SP31_T = np.array([1.0, 1.4, 1.7, 2.2, 2.5, 2.8])
SP32_T = np.array([1.3, 1.5, 2.0, 2.3, 2.7, 3.2])
SP33_T = np.array([1.25, 1.45, 1.95, 2.25, 2.65, 3.15])

# ============================================================
# Helper functions (Saenz + biaxial envelope)
# ============================================================
def saenz_params_biaxial(SIGCP, EPSCP):
    RP = eps_u / eps_c
    ES = SIGCP / EPSCP
    EU_num = fu * SIGCP / fc
    EU_den = eps_u * EPSCP / eps_c
    EU = EU_num / EU_den if abs(EU_den) > 1e-30 else E0
    denom = RP * (RP - 1)**2
    A = (E0/EU + (RP-2)*RP**2*E0/ES - (2*RP+1)*(RP-1)**2) / denom if abs(denom) > 1e-30 else 0
    B = 2*E0/ES - 3 - 2*A
    C_ = 2 - E0/ES + A
    return A, B, C_

def saenz_stress_biaxial(e, SIGCP, EPSCP):
    A_b, B_b, C_b = saenz_params_biaxial(SIGCP, EPSCP)
    DE = e / EPSCP
    if DE <= 0:
        return E0 * e
    denom = 1 + A_b*DE + B_b*DE**2 + C_b*DE**3
    return E0 * e / denom if abs(denom) > 1e-30 else 0.0

def biaxial_envelope_py(P1, P2, P3):
    SIGCP, EPSCP, FALSTR = fc, eps_c, ft
    if P3 >= 0:
        return SIGCP, EPSCP, FALSTR
    if P1 >= 0 and P2 >= 0:
        FALSTR = ft * (1 - P3 / fc)
        return SIGCP, EPSCP, FALSTR
    if P1 < 0:
        TEMP = P1 / fc
        I = 1
        for k in range(1, 6):
            I = k + 1
            if TEMP < SP1_T[k]:
                break
        I = min(I, 5)
        J = I - 1
        DSP = SP1_T[I] - SP1_T[J]
        FRAC = (TEMP - SP1_T[J]) / DSP if abs(DSP) > 1e-15 else 0.0
        SP31I = SP31_T[J] + FRAC * (SP31_T[I] - SP31_T[J])
        SP32I = SP32_T[J] + FRAC * (SP32_T[I] - SP32_T[J])
        SP33I = SP33_T[J] + FRAC * (SP33_T[I] - SP33_T[J])
    else:
        SP31I, SP32I, SP33I = SP31_T[0], SP32_T[0], SP33_T[0]
        TEMP = 0.0
    RATIO = P2 / fc
    if RATIO > BETA * SP32I:
        denom = SP33I - BETA * SP32I
        SLOPE = (SP33I - SP32I) / denom if abs(denom) > 1e-15 else 0.0
        SIGCP = SP32I * fc + SLOPE * (P2 - BETA * SP32I * fc)
    else:
        denom = BETA * SP32I - TEMP
        SLOPE = (SP32I - SP31I) / denom if abs(denom) > 1e-15 else 0.0
        SIGCP = SP31I * fc + SLOPE * (P2 - P1)
        if P1 > 0:
            SIGCP = SP31I * fc + SLOPE * P2
    if SIGCP < fc:
        EPSCP = GAMA * eps_c * SIGCP / fc
    if P1 >= 0:
        FALSTR = ft * (1 - P2/SIGCP) * (1 - P3/SIGCP)
        if FALSTR < 0.001 * ft:
            FALSTR = ft
    else:
        FALSTR = SIGCP
    return SIGCP, EPSCP, FALSTR

# ============================================================
# Step-by-step simulation to find peak compressive stress
# ============================================================
def find_compressive_peak(ratio_x, ratio_y, nsteps=600):
    """
    Simulate biaxial loading. Returns the peak stress point
    as (sigma_1, sigma_2) in-plane principal stresses (sigma_1 >= sigma_2).
    """
    max_strain = 0.015
    deps_scale = max_strain / nsteps / max(abs(ratio_x), abs(ratio_y))
    
    old_sx, old_sy = 0.0, 0.0
    eps_max_comp = 0.0
    
    peak_sx, peak_sy = 0.0, 0.0
    peak_norm = 0.0
    
    for istep in range(1, nsteps + 1):
        ex = ratio_x * deps_scale * istep
        ey = ratio_y * deps_scale * istep
        
        if ex >= 0:
            break  # only handle compression in x
        
        # Biaxial envelope from old stress
        s_avg = (old_sx + old_sy) / 2
        s_diff = (old_sx - old_sy) / 2
        sp1 = s_avg + abs(s_diff)
        sp2 = s_avg - abs(s_diff)
        sp3 = 0.0
        P = sorted([sp1, sp2, sp3], reverse=True)
        SIGCP, EPSCP, FALSTR = biaxial_envelope_py(P[0], P[1], P[2])
        
        if ex <= eps_max_comp:
            sig_s = saenz_stress_biaxial(ex, SIGCP, EPSCP)
            E_sec = sig_s / ex if abs(ex) > 1e-15 else E0
            eps_max_comp = ex
        else:
            sig_s = saenz_stress_biaxial(eps_max_comp, SIGCP, EPSCP)
            E_sec = sig_s / eps_max_comp if abs(eps_max_comp) > 1e-15 else E0
        
        E_sec = max(E_sec, E0 * STIFAC)
        C11 = E_sec / (1 - nu**2)
        C12 = nu * C11
        sx = C11 * ex + C12 * ey
        sy = C12 * ex + C11 * ey
        
        norm = np.sqrt(sx**2 + sy**2)
        if norm > peak_norm:
            peak_norm = norm
            peak_sx, peak_sy = sx, sy
        
        old_sx, old_sy = sx, sy
    
    # Return as principal stresses (sigma_1 >= sigma_2)
    s1 = max(peak_sx, peak_sy)
    s2 = min(peak_sx, peak_sy)
    return s1, s2

# ============================================================
# Compute envelope regions
# ============================================================
print("Computing biaxial failure envelope...")

# --- C-C Region: biaxial compression ---
# Parameterize by sigma_2/sigma_1 ratio (both negative)
# Strain ratios that produce different stress ratios
# We vary the strain ratio and record the resulting stress ratio at peak
n_cc = 50
# ey/ex from 0 (eps_yy=0, produces sigma_yy>0 from Poisson) to 1 (equal biaxial)
strain_ratios = np.linspace(0, 1, n_cc + 1)
cc_s1, cc_s2 = [], []
for r in strain_ratios:
    s1, s2 = find_compressive_peak(-1.0, -r)
    cc_s1.append(s1)
    cc_s2.append(s2)

cc_s1 = np.array(cc_s1)
cc_s2 = np.array(cc_s2)

# Mirror for sigma_1 <-> sigma_2 symmetry
cc_s1_sym = np.flip(cc_s2[:-1])
cc_s2_sym = np.flip(cc_s1[:-1])
cc_s1_full = np.concatenate([cc_s1, cc_s1_sym])
cc_s2_full = np.concatenate([cc_s2, cc_s2_sym])

# Find the uniaxial point (where sigma_1 ~ 0)
idx_uni = np.argmin(np.abs(cc_s1))
print(f"  Uniaxial comp approx at ratio={strain_ratios[idx_uni]:.3f}: "
      f"sig=({cc_s1[idx_uni]:.2f}, {cc_s2[idx_uni]:.2f})")

# --- C-T Region: compression + tension ---
# The failure in C-T is when tensile principal stress reaches FALSTR.
# For plane stress: P1=sigma_t, P2=0, P3=sigma_c
# FALSTR = ft * (1 - sigma_c / fc)
# So the C-T envelope is:  sigma_t = ft * (1 - sigma_c / fc)
#   i.e., sigma_1 = ft * (1 + sigma_2 / |fc|)  (sigma_2 < 0)
ct_s2 = np.linspace(fc, 0, 100)  # from fc to 0
ct_s1 = ft * (1 - ct_s2 / fc)     # FALSTR formula

# --- T-T Region: biaxial tension ---
# For elastic material: cracking when max principal stress = ft
# With proportional loading eps_yy = r * eps_xx (both > 0):
#   sigma_1 = E/(1-nu^2) * (1+nu*r) * eps
#   sigma_2 = E/(1-nu^2) * (nu+r) * eps
# Failure at sigma_max = ft
# For r in [0, 1]:  sigma_1 >= sigma_2 if r <= 1
#   sigma_1 = ft
#   sigma_2 = ft * (nu + r) / (1 + nu*r)
tt_r = np.linspace(0, 1, 50)
tt_s1 = np.full_like(tt_r, ft)
tt_s2 = ft * (nu + tt_r) / (1 + nu * tt_r)

# By symmetry
tt_s1_sym = np.flip(tt_s2[:-1])
tt_s2_sym = np.flip(tt_s1[:-1])
tt_s1_full = np.concatenate([tt_s1, tt_s1_sym])
tt_s2_full = np.concatenate([tt_s2, tt_s2_sym])

# ============================================================
# Kupfer reference data (approximate)
# sigma_1/|fc| vs sigma_2/|fc|
# ============================================================
# C-C quadrant (Kupfer 1969)
kupfer_cc_s1fc = np.array([0.0, -0.10, -0.20, -0.35, -0.52, -0.65, -0.80, -1.0, -1.16])
kupfer_cc_s2fc = np.array([-1.0, -1.02, -1.05, -1.10, -1.15, -1.175, -1.16, -1.16, -1.16])

# C-T quadrant (Kupfer)
kupfer_ct_s1fc = np.array([0.0, 0.02, 0.05, 0.08, 0.10])
kupfer_ct_s2fc = np.array([-1.0, -0.85, -0.6, -0.35, 0.0])

# ============================================================
# Figure
# ============================================================
fig, ax = plt.subplots(figsize=(10, 10))

fca = abs(fc)

# --- C-C region ---
ax.plot(cc_s1_full/fca, cc_s2_full/fca, 'b-', linewidth=2.5,
        label='ADINA Model', zorder=3)

# --- C-T region ---
ax.plot(ct_s1/fca, ct_s2/fca, 'b-', linewidth=2.5, zorder=3)
# Symmetric C-T
ax.plot(ct_s2/fca, ct_s1/fca, 'b-', linewidth=2.5, zorder=3)

# --- T-T region ---
ax.plot(tt_s1_full/fca, tt_s2_full/fca, 'b-', linewidth=2.5, zorder=3)

# --- Kupfer reference ---
# C-C + symmetric
ax.plot(kupfer_cc_s1fc, kupfer_cc_s2fc, 'ko--', markersize=5, linewidth=1.5,
        alpha=0.6, label='Kupfer et al. (1969)', zorder=2)
ax.plot(kupfer_cc_s2fc, kupfer_cc_s1fc, 'ko--', markersize=5, linewidth=1.5,
        alpha=0.6, zorder=2)
# C-T + symmetric
ax.plot(kupfer_ct_s1fc, kupfer_ct_s2fc, 'ko--', markersize=5, linewidth=1.5,
        alpha=0.6, zorder=2)
ax.plot(kupfer_ct_s2fc, kupfer_ct_s1fc, 'ko--', markersize=5, linewidth=1.5,
        alpha=0.6, zorder=2)

# --- Key annotations ---
# Uniaxial compression
ax.plot(0, -1, 'rs', markersize=10, zorder=5)
ax.annotate('Uniaxial $f_c$\n$(0, -1)$', xy=(0, -1),
            xytext=(0.12, -0.9), fontsize=10, color='red',
            arrowprops=dict(arrowstyle='->', color='red', lw=1.2))

# Equal biaxial compression
eq_s1 = cc_s1[-1]/fca
eq_s2 = cc_s2[-1]/fca
ax.plot(eq_s1, eq_s2, 'r^', markersize=10, zorder=5)
ax.annotate(f'Equal biaxial\n$({eq_s1:.2f}, {eq_s2:.2f})$',
            xy=(eq_s1, eq_s2), xytext=(-1.0, -1.8), fontsize=10, color='red',
            arrowprops=dict(arrowstyle='->', color='red', lw=1.2))

# Uniaxial tension
ax.plot(ft/fca, 0, 'gs', markersize=10, zorder=5)
ax.annotate(f'$f_t/|f_c|={ft/fca:.2f}$', xy=(ft/fca, 0),
            xytext=(0.2, 0.05), fontsize=10, color='green',
            arrowprops=dict(arrowstyle='->', color='green', lw=1.2))

# --- Quadrant labels ---
ax.text(-0.7, -0.7, 'C-C', fontsize=16, ha='center', va='center',
        color='blue', alpha=0.25, fontweight='bold')
ax.text(0.04, -0.5, 'T-C', fontsize=14, ha='center', va='center',
        color='blue', alpha=0.25, fontweight='bold')
ax.text(-0.5, 0.04, 'C-T', fontsize=14, ha='center', va='center',
        color='blue', alpha=0.25, fontweight='bold')
ax.text(0.05, 0.05, 'T-T', fontsize=12, ha='center', va='center',
        color='green', alpha=0.35, fontweight='bold')

# --- Formatting ---
ax.axhline(y=0, color='k', linewidth=0.6)
ax.axvline(x=0, color='k', linewidth=0.6)
ax.plot([-2, 0.2], [-2, 0.2], 'k:', linewidth=0.5, alpha=0.3)  # symmetry line

ax.set_xlabel('$\\sigma_1 / |f_c|$', fontsize=14)
ax.set_ylabel('$\\sigma_2 / |f_c|$', fontsize=14)
ax.set_title('Biaxial Failure Envelope — ADINA Concrete Model', fontsize=15)
ax.legend(loc='lower left', fontsize=11, framealpha=0.9)
ax.set_aspect('equal')

ax.set_xlim([-1.75, 0.25])
ax.set_ylim([-1.75, 0.25])
ax.grid(True, alpha=0.15)

# Add tick marks at key values
ax.set_xticks([-1.5, -1.0, -0.5, 0.0])
ax.set_yticks([-1.5, -1.0, -0.5, 0.0])

fig.tight_layout()
fig.savefig('plot_envelope.png', dpi=150, bbox_inches='tight')
print('\nSaved: plot_envelope.png')

# ============================================================
# Summary
# ============================================================
print('\n' + '='*60)
print('  Biaxial Failure Envelope Summary')
print('='*60)
print(f'  fc = {fc:.1f} MPa,  ft = {ft:.1f} MPa,  nu = {nu}')
print(f'  ft/|fc| = {ft/fca:.3f}')
print()
print('  Key points on the envelope (sigma_1/|fc|, sigma_2/|fc|):')
print(f'    Uniaxial compression:    (  0.00, -1.00)')
print(f'    Equal biaxial comp:      ({eq_s1:+.3f}, {eq_s2:+.3f})')
print(f'    Biaxial comp factor:     {abs(eq_s2):.3f}')
print(f'    Uniaxial tension:        ({ft/fca:+.3f},  0.00)')
print(f'    Equal biaxial tension:   ({ft/fca:+.3f}, {ft/fca:+.3f})')
print()
print('  C-C envelope (selected points):')
for i in range(0, len(strain_ratios), 10):
    r = strain_ratios[i]
    s1, s2 = cc_s1[i]/fca, cc_s2[i]/fca
    print(f'    ey/ex={r:.2f} -> ({s1:+.3f}, {s2:+.3f})')
r = strain_ratios[-1]
s1, s2 = cc_s1[-1]/fca, cc_s2[-1]/fca
print(f'    ey/ex={r:.2f} -> ({s1:+.3f}, {s2:+.3f})')
print()
print('Done.')
