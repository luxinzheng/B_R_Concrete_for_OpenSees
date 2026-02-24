"""
Biaxial Monotonic Loading: Results vs Theoretical Comparison

Theoretical model: Isotropic secant stiffness based on uniaxial Saenz curve
  - E_sec = saenz(eps_xx) / eps_xx  for compression (eps_xx < 0)
  - sigma = C(E_sec, nu) * strain   (total-strain approach)
  - Cracking: when principal stress > ft, stress -> 0
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.size'] = 10
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

# ============================================================
# Saenz curve
# ============================================================
def saenz_params():
    RP = eps_u / eps_c
    ES = fc / eps_c
    EU = fu / eps_u
    A = (E0/EU + (RP-2)*RP**2*E0/ES - (2*RP+1)*(RP-1)**2) / (RP*(RP-1)**2)
    B = 2*E0/ES - 3 - 2*A
    C = 2 - E0/ES + A
    return A, B, C

def saenz_stress_scalar(e):
    """Uniaxial Saenz stress at strain e (< 0 for compression)."""
    A, B, C_ = saenz_params()
    DE = e / eps_c
    if DE <= 0:
        return E0 * e
    denom = 1 + A*DE + B*DE**2 + C_*DE**3
    if abs(denom) < 1e-30:
        return 0.0
    return E0 * e / denom

def saenz_E_sec(e):
    """Secant modulus from uniaxial Saenz curve at strain e."""
    if abs(e) < 1e-15:
        return E0
    return saenz_stress_scalar(e) / e

# ============================================================
# Biaxial envelope (ported from ADINA Fortran code)
# ============================================================
BETA = 0.75
GAMA = 1.0

# Spline tables (same as in Fortran)
SP1_T  = np.array([0.0, 0.25, 0.5, 0.75, 1.0, 1.2])
SP31_T = np.array([1.0, 1.4, 1.7, 2.2, 2.5, 2.8])
SP32_T = np.array([1.3, 1.5, 2.0, 2.3, 2.7, 3.2])
SP33_T = np.array([1.25, 1.45, 1.95, 2.25, 2.65, 3.15])

def biaxial_envelope_py(P1, P2, P3):
    """
    Compute biaxial-adjusted SIGCP and EPSCP.
    P1 >= P2 >= P3 (sorted principal stresses)
    Returns (SIGCP, EPSCP, FALSTR)
    """
    SIGCP = fc
    EPSCP = eps_c
    FALSTR = ft
    
    if P3 >= 0:
        return SIGCP, EPSCP, FALSTR
    
    if P1 >= 0 and P2 >= 0:
        FALSTR = ft * (1 - P3 / fc)
        return SIGCP, EPSCP, FALSTR
    
    # At least P2 and P3 compressive
    if P1 < 0:
        TEMP = P1 / fc
        # Interpolate in SP1_T
        I = 1
        for k in range(1, 6):
            I = k + 1
            if TEMP < SP1_T[k]:
                break
        I = min(I, 5)  # 0-indexed: max index 5
        J = I - 1
        DSP = SP1_T[I] - SP1_T[J]
        FRAC = (TEMP - SP1_T[J]) / DSP if abs(DSP) > 1e-15 else 0.0
        SP31I = SP31_T[J] + FRAC * (SP31_T[I] - SP31_T[J])
        SP32I = SP32_T[J] + FRAC * (SP32_T[I] - SP32_T[J])
        SP33I = SP33_T[J] + FRAC * (SP33_T[I] - SP33_T[J])
    else:
        SP31I = SP31_T[0]
        SP32I = SP32_T[0]
        SP33I = SP33_T[0]
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

def saenz_stress_biaxial(e, SIGCP, EPSCP):
    """Saenz stress with biaxial-adjusted parameters."""
    A_b, B_b, C_b = saenz_params_biaxial(SIGCP, EPSCP)
    DE = e / EPSCP
    if DE <= 0:
        return E0 * e
    denom = 1 + A_b*DE + B_b*DE**2 + C_b*DE**3
    if abs(denom) < 1e-30:
        return 0.0
    return E0 * e / denom

def saenz_params_biaxial(SIGCP, EPSCP):
    """Saenz curve params with biaxial peak."""
    RP = eps_u / eps_c
    ES = SIGCP / EPSCP
    EU_num = fu * SIGCP / fc
    EU_den = eps_u * EPSCP / eps_c
    EU = EU_num / EU_den if abs(EU_den) > 1e-30 else E0
    
    denom = RP * (RP - 1)**2
    if abs(denom) < 1e-30:
        A = 0.0
    else:
        A = (E0/EU + (RP-2)*RP**2*E0/ES - (2*RP+1)*(RP-1)**2) / denom
    B = 2*E0/ES - 3 - 2*A
    C_ = 2 - E0/ES + A
    return A, B, C_

# ============================================================
# Theoretical biaxial stress (Saenz secant + biaxial envelope)
# ============================================================
def compute_theory_biaxial(eps_xx_arr, eps_yy_arr):
    """
    Compute theoretical biaxial stresses step-by-step.
    Uses old stress state to determine biaxial envelope,
    then computes Saenz E_sec with adjusted parameters.
    """
    n = len(eps_xx_arr)
    sig_xx_th = np.zeros(n)
    sig_yy_th = np.zeros(n)
    
    cracked = False
    eps_max_comp = 0.0   # most negative eps_xx reached
    STIFAC = 1e-4
    
    # Track old stress for biaxial envelope
    old_sig_xx = 0.0
    old_sig_yy = 0.0
    
    for i in range(n):
        ex = eps_xx_arr[i]
        ey = eps_yy_arr[i]
        
        if ex < 0:
            # === Biaxial envelope from OLD stress ===
            # Principal stresses (plane stress: sigma_zz = 0)
            s_avg = (old_sig_xx + old_sig_yy) / 2
            s_diff = (old_sig_xx - old_sig_yy) / 2
            tau = 0  # no shear in our tests
            s_r = np.sqrt(s_diff**2 + tau**2)
            sp1 = s_avg + s_r
            sp2 = s_avg - s_r
            sp3 = 0.0  # sigma_zz
            # Sort: P1 >= P2 >= P3
            P = sorted([sp1, sp2, sp3], reverse=True)
            SIGCP, EPSCP, FALSTR = biaxial_envelope_py(P[0], P[1], P[2])
            
            # === Saenz E_sec with biaxial parameters ===
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
            
            sig_xx_th[i] = C11 * ex + C12 * ey
            sig_yy_th[i] = C12 * ex + C11 * ey
            
        elif ex >= 0:
            if cracked:
                sig_xx_th[i] = 0.0
                sig_yy_th[i] = 0.0
            else:
                C11 = E0 / (1 - nu**2)
                C12 = nu * C11
                sig_xx_th[i] = C11 * ex + C12 * ey
                sig_yy_th[i] = C12 * ex + C11 * ey
                
                sig_avg = (sig_xx_th[i] + sig_yy_th[i]) / 2
                sig_diff = (sig_xx_th[i] - sig_yy_th[i]) / 2
                sig_p1 = sig_avg + np.sqrt(sig_diff**2)
                if sig_p1 > ft:
                    cracked = True
        
        old_sig_xx = sig_xx_th[i]
        old_sig_yy = sig_yy_th[i]
    
    return sig_xx_th, sig_yy_th

# ============================================================
# Read CSV
# ============================================================
def read_biaxial_csv(filename):
    steps, ex, ey, sx, sy, txy = [], [], [], [], [], []
    with open(filename, 'r') as f:
        f.readline()  # header
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = [p.strip() for p in line.split(',')]
            steps.append(int(parts[0]))
            ex.append(float(parts[1]))
            ey.append(float(parts[2]))
            sx.append(float(parts[3]))
            sy.append(float(parts[4]))
            txy.append(float(parts[5]))
    return (np.array(steps), np.array(ex), np.array(ey),
            np.array(sx), np.array(sy), np.array(txy))

# ============================================================
# Case definitions
# ============================================================
cases = [
    {'file': 'test_biaxial_case1.csv', 'rx': -1.0, 'ry': -1.0,
     'name': 'Case 1: -1:-1 (biaxial compression)'},
    {'file': 'test_biaxial_case2.csv', 'rx':  1.0, 'ry':  1.0,
     'name': 'Case 2: 1:1 (biaxial tension)'},
    {'file': 'test_biaxial_case3.csv', 'rx': -1.0, 'ry': -0.5,
     'name': 'Case 3: -1:-0.5 (unequal biaxial comp)'},
    {'file': 'test_biaxial_case4.csv', 'rx': -1.0, 'ry':  0.1,
     'name': r'Case 4: -1:0.1 (comp + tension)'},
    {'file': 'test_biaxial_case5.csv', 'rx': -1.0, 'ry':  0.2,
     'name': r'Case 5: -1:0.2 ($\sigma_{yy}\approx 0$)'},
]

# ============================================================
# Saenz envelope for reference
# ============================================================
eps_env = np.linspace(0, -0.012, 500)
sig_env = np.array([saenz_stress_scalar(e) for e in eps_env])

# ============================================================
# Figure: 5 cases, each with sig_xx vs eps_xx and sig_yy vs eps_yy
# ============================================================
fig, axes = plt.subplots(5, 2, figsize=(14, 22))

for ic, case in enumerate(cases):
    steps, ex, ey, sx, sy, txy = read_biaxial_csv(case['file'])
    
    # Compute theory
    sx_th, sy_th = compute_theory_biaxial(ex, ey)
    
    # --- Left: sig_xx vs eps_xx ---
    ax = axes[ic, 0]
    ax.plot(ex*100, sx, 'b-', linewidth=1.5, label='Computed $\\sigma_{xx}$')
    ax.plot(ex*100, sx_th, 'r--', linewidth=1.2, alpha=0.7, label='Theory $\\sigma_{xx}$')
    
    # Overlay uniaxial Saenz envelope for compression cases
    if case['rx'] < 0:
        ax.plot(eps_env*100, sig_env, 'k:', linewidth=0.8, alpha=0.3,
                label='Uniaxial Saenz')
    
    ax.axhline(y=0, color='k', linewidth=0.3)
    ax.axvline(x=0, color='k', linewidth=0.3)
    ax.set_xlabel('$\\varepsilon_{xx}$ (%)')
    ax.set_ylabel('$\\sigma_{xx}$ (MPa)')
    ax.set_title(case['name'] + ' — $\\sigma_{xx}$')
    ax.legend(fontsize=8, loc='best')
    ax.grid(True, alpha=0.2)
    
    # --- Right: sig_yy vs eps_yy ---
    ax = axes[ic, 1]
    ax.plot(ey*100, sy, 'b-', linewidth=1.5, label='Computed $\\sigma_{yy}$')
    ax.plot(ey*100, sy_th, 'r--', linewidth=1.2, alpha=0.7, label='Theory $\\sigma_{yy}$')
    
    ax.axhline(y=0, color='k', linewidth=0.3)
    ax.axvline(x=0, color='k', linewidth=0.3)
    ax.set_xlabel('$\\varepsilon_{yy}$ (%)')
    ax.set_ylabel('$\\sigma_{yy}$ (MPa)')
    ax.set_title(case['name'] + ' — $\\sigma_{yy}$')
    ax.legend(fontsize=8, loc='best')
    ax.grid(True, alpha=0.2)
    
    # Print comparison
    diff_xx = np.max(np.abs(sx - sx_th))
    diff_yy = np.max(np.abs(sy - sy_th))
    print(f"{case['name']}")
    print(f"  Max |sig_xx error| = {diff_xx:.4f} MPa")
    print(f"  Max |sig_yy error| = {diff_yy:.4f} MPa")
    if diff_xx < 0.01 and diff_yy < 0.01:
        print(f"  => PERFECT MATCH")
    elif diff_xx < 1.0 and diff_yy < 1.0:
        print(f"  => GOOD (small diff)")
    else:
        print(f"  => MISMATCH (needs investigation)")

fig.tight_layout()
fig.savefig('plot_biaxial.png', dpi=150, bbox_inches='tight')
print('\nSaved: plot_biaxial.png')

# ============================================================
# Figure 2: Step-by-step comparison table (all 5 cases)
# ============================================================
fig2, axes2 = plt.subplots(5, 2, figsize=(14, 22))

for ic, case in enumerate(cases):
    steps, ex, ey, sx, sy, txy = read_biaxial_csv(case['file'])
    sx_th, sy_th = compute_theory_biaxial(ex, ey)
    
    # Left: error sig_xx
    ax = axes2[ic, 0]
    ax.plot(steps, sx - sx_th, 'g-', linewidth=1.0)
    ax.axhline(y=0, color='k', linewidth=0.3)
    ax.set_xlabel('Step')
    ax.set_ylabel('$\\sigma_{xx}^{comp} - \\sigma_{xx}^{theory}$ (MPa)')
    ax.set_title(case['name'] + ' — $\\sigma_{xx}$ error')
    ax.grid(True, alpha=0.2)
    
    # Right: error sig_yy
    ax = axes2[ic, 1]
    ax.plot(steps, sy - sy_th, 'g-', linewidth=1.0)
    ax.axhline(y=0, color='k', linewidth=0.3)
    ax.set_xlabel('Step')
    ax.set_ylabel('$\\sigma_{yy}^{comp} - \\sigma_{yy}^{theory}$ (MPa)')
    ax.set_title(case['name'] + ' — $\\sigma_{yy}$ error')
    ax.grid(True, alpha=0.2)

fig2.tight_layout()
fig2.savefig('plot_biaxial_error.png', dpi=150, bbox_inches='tight')
print('Saved: plot_biaxial_error.png')

print('\nAll biaxial plots generated.')
