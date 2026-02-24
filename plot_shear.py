"""
Pure Shear Test: Compare forumat.f90 output with theoretical ADINA model

Theory for pure shear (eps_xx=0, eps_yy=0, gamma_xy=gamma):

Phase 1 (Elastic, uncracked):
  Hypoelastic: tau_new = tau_old + G_eff * d_gamma
  G_eff = m_E / (2*(1+nu))
  m_E = weighted average of Saenz tangent for principal directions
  Cracking at sigma_1 > m_FALSTR (biaxial-adjusted tensile strength)

Phase 2 (Cracked, NUMCRK=1 at 45 deg):
  In crack coordinates: eps_n = gamma/2, eps_t = -gamma/2
  sigma_n = 0 (cracked)
  sigma_t += C_tt * d_eps_t,  C_tt = m_E_saenz / (1-nu^2)
  Global: sigma_xx = sigma_yy = sigma_t/2, tau = -sigma_t/2
  Crushing when |sigma_t| >= |fc_biaxial|

Phase 3 (Crushed): all stresses = 0
"""

import numpy as np
import matplotlib.pyplot as plt

# ========== Material Parameters ==========
E0     = 30000.0   # MPa
nu     = 0.2
ft     = 3.0       # tensile strength (positive)
fc     = -30.0     # compressive strength (negative)
eps_c  = -0.002    # peak compressive strain (negative)
fu     = -6.0      # ultimate compressive stress (negative)
eps_u  = -0.005    # ultimate compressive strain (negative)
STIFAC = 1.0e-4
SHEFAC = 0.5
RKAPA  = 0.7
ALFA   = 0.12

G0 = E0 / (2.0*(1.0+nu))  # = 12500

# ========== Saenz Curve Parameters ==========
def compute_saenz_params(E, SIGCP, EPSCP, SIGMAU, EPSU, SIGMAC, EPSC):
    """Compute Saenz curve shape parameters A, B, C"""
    RP = EPSU / EPSC
    ES = SIGCP / EPSCP
    if abs(EPSU * EPSCP / EPSC) > 1e-30:
        EU = (SIGMAU * SIGCP / SIGMAC) / (EPSU * EPSCP / EPSC)
    else:
        EU = E
    
    if abs(RP * (RP - 1.0)**2) < 1e-30:
        RAM5 = 0.0
    else:
        RAM5 = (E/EU + (RP-2)*RP**2*E/ES - (2*RP+1)*(RP-1)**2) / (RP*(RP-1)**2)
    RBM5 = 2*E/ES - 3 - 2*RAM5
    RCM5 = 2 - E/ES + RAM5
    return RAM5, RBM5, RCM5

def saenz_tangent(DE, E, RAM5, RBM5, RCM5):
    """Saenz tangent at normalized strain DE = eps/eps_c (positive for compression)"""
    if DE <= 0:
        return E
    denom = 1.0 + RAM5*DE + RBM5*DE**2 + RCM5*DE**3
    if abs(denom) < 1e-30:
        return 0.0
    return E * (1.0 - RBM5*DE**2 - 2.0*RCM5*DE**3) / denom**2

def saenz_stress(eps_val, E, EPSCP, RAM5, RBM5, RCM5):
    """Saenz stress at strain eps_val"""
    DE = eps_val / EPSCP
    if DE <= 0:
        return E * eps_val
    denom = 1.0 + RAM5*DE + RBM5*DE**2 + RCM5*DE**3
    if abs(denom) < 1e-30:
        return 0.0
    return E * eps_val / denom

# 3-point Gauss quadrature
XG3 = np.array([-0.774596669241483, 0.0, 0.774596669241483])
WG3 = np.array([0.555555555555556, 0.888888888888889, 0.555555555555556])

def gauss_averaged_tangent(eps_old, eps_new, EPSCP, E, RAM5, RBM5, RCM5):
    """3-point Gauss average of Saenz tangent over strain path [eps_old, eps_new]"""
    result = 0.0
    for L in range(3):
        E1 = XG3[L]
        DE = eps_old + (1.0 + E1) * (eps_new - eps_old) / 2.0
        DE = DE / EPSCP  # normalize
        TY = saenz_tangent(DE, E, RAM5, RBM5, RCM5)
        result += 0.5 * WG3[L] * TY
    return result


# ========== Theoretical Model ==========
def compute_EVV_full(sx, sy, txy, sz=0.0):
    """Compute Drucker-Prager equivalent stress from full 2D stress state"""
    SBAR = ((sx - sy)**2 + (sx - sz)**2 + (sy - sz)**2) / 6.0
    TMM = (sx + sy + sz) / 3.0
    return np.sqrt(SBAR + txy**2) + 3.0*ALFA*TMM

def compute_theory():
    """Step-by-step theoretical model replicating ADINA algorithm for pure shear"""
    NSTEPS = 100
    dgamma = 0.005 / NSTEPS
    
    results = []
    results.append((0, 0.0, 0.0, 0.0, 0.0, 'elastic'))
    
    # State variables
    stress_old = np.array([0.0, 0.0, 0.0, 0.0])  # sx, sy, txy, sz
    gamma_old = 0.0
    eps_old = np.array([0.0, 0.0, 0.0, 0.0])  # exx, eyy, gxy, ezz
    EVMAX = 0.0
    
    # Saenz parameters (uniaxial)
    RAM5, RBM5, RCM5 = compute_saenz_params(E0, fc, eps_c, fu, eps_u, fc, eps_c)
    
    # Phase tracking
    cracked = False
    crushed = False
    sigma_t_old_crack = 0.0  # stress in t-direction (crack coords)
    eps_t_old_crack = 0.0    # strain in t-direction (crack coords)
    eps_n_old_crack = 0.0    # strain in n-direction (crack coords)
    CRKSTR = np.zeros(3)
    m_ILFSET = 0
    
    for step in range(1, NSTEPS + 1):
        gamma_new = dgamma * step
        strain_new = np.array([0.0, 0.0, gamma_new, 0.0])  # exx, eyy, gxy, ezz
        m_ILFSET = 0  # reset each step
        
        if crushed:
            # Phase 3: fully failed (STRAIN(1) = 0 >= 0 → zero stress)
            results.append((step, gamma_new, 0.0, 0.0, 0.0, 'crushed'))
            eps_old = strain_new.copy()
            continue
        
        # Compute EVV from old stress (IKAS check)
        sx_o, sy_o, txy_o, sz_o = stress_old
        EVV = compute_EVV_full(sx_o, sy_o, txy_o, sz_o)
        IKAS = 1 if (abs(EVV - EVMAX) <= abs(EVMAX)*1e-8 or EVMAX < 1e-15) else -1
        
        if cracked:
            # Phase 2: cracked state (Section 7, NUMCRK=1 at 45 deg)
            
            # Rotate old stress to crack coords (45 deg)
            # sigma_n = (sx+sy)/2 + (sx-sy)/2*cos90 + txy*sin90 = (sx+sy)/2 + txy
            # sigma_t = (sx+sy)/2 - (sx-sy)/2*cos90 - txy*sin90 = (sx+sy)/2 - txy
            sigma_n_old = (sx_o + sy_o)/2 + txy_o
            sigma_t_old = (sx_o + sy_o)/2 - txy_o
            sigma_z_old = sz_o  # z is unrotated
            tau_nt_old = -(sx_o - sy_o)/2 * 1.0 + txy_o * 0.0  # sin90=1, cos90=0
            
            # Rotate strains to crack coords
            eps_n_old = eps_old[2] / 2.0   # gamma_old/2
            eps_t_old = -eps_old[2] / 2.0  # -gamma_old/2
            eps_n_new = gamma_new / 2.0
            eps_t_new = -gamma_new / 2.0
            
            d_eps_n = eps_n_new - eps_n_old
            d_eps_t = eps_t_new - eps_t_old
            
            # NUMCRK=1 (crack open since eps_n > 0)
            
            # Biaxial envelope from crack-coord stresses (OLD)
            sigp = np.array([sigma_n_old, sigma_t_old, sigma_z_old])
            P = sorted(sigp, reverse=True)  # P1>=P2>=P3
            
            # For this stress state (nearly uniaxial compression in t):
            SIGCP_eff = fc
            EPSCP_eff = eps_c
            
            RAM5_c, RBM5_c, RCM5_c = compute_saenz_params(E0, SIGCP_eff, EPSCP_eff,
                                                            fu, eps_u, fc, eps_c)
            
            if IKAS > 0:
                # Loading: compute Saenz tangent for each direction
                # m_YP(1) for n-direction, m_YP(2) for t-direction, m_YP(3) for z
                DENM = abs(sigma_n_old) + abs(sigma_t_old) + abs(sigma_z_old)
                
                if DENM <= 0.00001 * ft:
                    m_E_eff = E0
                else:
                    # J=1 (n-direction): sigma_n is usually ~0 or tensile → m_YP(1) = E0
                    m_YP1 = E0
                    
                    # J=2 (t-direction): compressive → compute Saenz tangent
                    if sigma_t_old >= 0.001 * SIGCP_eff:
                        m_YP2 = E0
                    else:
                        m_YP2 = gauss_averaged_tangent(eps_t_old, eps_t_new, EPSCP_eff,
                                                        E0, RAM5_c, RBM5_c, RCM5_c)
                    
                    # J=4 (z-direction): sigma_z ~0 → m_YP(3) = E0
                    m_YP3 = E0
                    
                    # Weighted average
                    m_E_eff = (abs(sigma_n_old)*m_YP1 + abs(sigma_t_old)*m_YP2 +
                               abs(sigma_z_old)*m_YP3) / DENM
            else:
                # Unloading: elastic
                m_E_eff = E0
            
            # Build constitutive for NUMCRK=1 (isotropic tangent)
            m_RK = m_E_eff / (3.0*(1.0 - 2.0*nu))
            m_G = m_E_eff / (2.0*(1.0 + nu))
            A2 = m_E_eff / (1.0 - nu**2)
            B2 = m_E_eff * nu / (1.0 - nu**2)
            
            # C(1,1) for NUMCRK=1 = STIFAC * ... ≈ 0
            # C(2,2) = A2, C(2,4) = B2, C(4,4) = A2
            # C(2,1) ≈ 0 (STIFAC scaled)
            
            # Incremental stress update in crack coords
            # sigma_t_new = sigma_t_old + A2*d_eps_t + B2*d_eps_z
            # sigma_z_new = sigma_z_old + B2*d_eps_t + A2*d_eps_z
            # d_eps_z = 0 (global eps_z = 0, z is unrotated)
            sigma_t_new = sigma_t_old + A2 * d_eps_t
            sigma_z_new = sigma_z_old + B2 * d_eps_t
            
            # Shear: sigma_nt += SHEFAC*G * d_gamma_nt
            # d_gamma_nt = 0 for pure shear at 45 deg
            tau_nt_new = tau_nt_old  # no change
            
            # Zero out cracked direction (NUMCRK=1: sigma_n = 0)
            sigma_n_new = 0.0
            
            # Apply shear retention factor (ETATAU = 1 since SHEFAC > 0.001)
            # tau_nt *= 1.0 (no change needed)
            
            # Transform to global
            # sigma_xx = sigma_n*cos²45 + sigma_t*sin²45 - 2*tau_nt*sin45*cos45
            #          = sigma_n*0.5 + sigma_t*0.5 - tau_nt
            # sigma_yy = sigma_n*0.5 + sigma_t*0.5 + tau_nt
            # tau_xy = (sigma_n - sigma_t)*sin45*cos45 + tau_nt*(cos²45 - sin²45)
            #        = (sigma_n - sigma_t)*0.5 + 0
            sx_new = sigma_n_new * 0.5 + sigma_t_new * 0.5 - tau_nt_new
            sy_new = sigma_n_new * 0.5 + sigma_t_new * 0.5 + tau_nt_new
            txy_new = (sigma_n_new - sigma_t_new) * 0.5
            sz_new = 0.0  # plane stress: sigma_z_global = 0
            
            # Check crushing: sort crack-coord principal stresses
            sigp_new = sorted([sigma_n_new, sigma_t_new, sigma_z_new])
            P3_new = sigp_new[0]  # most compressive
            if P3_new < SIGCP_eff:
                crushed = True
                results.append((step, gamma_new, txy_new, sx_new, sy_new, 'crushing'))
            else:
                results.append((step, gamma_new, txy_new, sx_new, sy_new, 'cracked'))
            
            # Update EVI and EVMAX
            EVI = compute_EVV_full(sx_new, sy_new, txy_new, sz_new)
            if EVI > EVMAX:
                EVMAX = EVI
            
            stress_old = np.array([sx_new, sy_new, txy_new, sz_new])
            eps_old = strain_new.copy()
            gamma_old = gamma_new
            continue
        
        # Phase 1: uncracked elastic with Saenz tangent
        d_gamma = gamma_new - gamma_old
        
        if IKAS > 0 and (abs(sx_o) + abs(sy_o) + abs(txy_o)) > 1e-10:
            # Loading: compute Saenz tangent (weighted average of principal directions)
            # Principal direction 2 (compressive): eps_2 = -gamma/2
            eps2_old = -gamma_old / 2.0
            eps2_new = -gamma_new / 2.0
            
            # Stress principal values
            tau_val = abs(txy_o)
            sigma1 = tau_val   # tensile
            sigma2 = -tau_val  # compressive
            
            # m_YP for tensile direction: sigma1 > 0 → m_YP1 = E0
            m_YP1 = E0
            
            # m_YP for compressive direction
            if sigma2 >= 0.001 * fc:
                m_YP2 = E0
            else:
                m_YP2 = gauss_averaged_tangent(eps2_old, eps2_new, eps_c, E0, RAM5, RBM5, RCM5)
            
            # Weighted average
            DENM = abs(sigma1) + abs(sigma2)
            if DENM > 1e-10:
                m_E = (abs(sigma1)*m_YP1 + abs(sigma2)*m_YP2) / DENM
            else:
                m_E = E0
        else:
            m_E = E0
        
        m_G = m_E / (2.0*(1.0 + nu))
        
        # Hypoelastic stress update: tau_new = tau_old + G * d_gamma
        tau_new = txy_o + m_G * d_gamma
        sx_new = 0.0
        sy_new = 0.0
        
        # Check cracking: principal stress sigma1 = tau_new
        # Biaxial-adjusted failure stress:
        m_FALSTR = ft * (1.0 - tau_new / abs(fc))
        if m_FALSTR < 0.001 * ft:
            m_FALSTR = ft
        
        if tau_new > m_FALSTR:
            # Cracking detected
            cracked = True
            m_ILFSET = 1
            results.append((step, gamma_new, tau_new, 0.0, 0.0, 'cracking'))
        else:
            results.append((step, gamma_new, tau_new, 0.0, 0.0, 'elastic'))
        
        # Update EVI and EVMAX
        EVI = compute_EVV_full(sx_new, sy_new, tau_new, 0.0)
        if IKAS == 1 and m_ILFSET == 1:
            EVMAX = EVI
        if EVI > EVMAX:
            EVMAX = EVI
        
        stress_old = np.array([sx_new, sy_new, tau_new, 0.0])
        eps_old = strain_new.copy()
        gamma_old = gamma_new
    
    return results


# ========== Read Fortran Output ==========
def read_fortran_csv(filename='test_shear.csv'):
    data = np.loadtxt(filename, delimiter=',', skiprows=1)
    return data  # columns: step, gamma, tau, sx, sy


# ========== Plotting ==========
def main():
    # Read Fortran output
    fort = read_fortran_csv()
    
    # Compute theoretical
    theory_raw = compute_theory()
    theory = np.array([(r[0], r[1], r[2], r[3], r[4]) for r in theory_raw])
    phases = [r[5] for r in theory_raw]
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Pure Shear Test: ADINA Model vs Theory', fontsize=14, fontweight='bold')
    
    # --- Plot 1: tau_xy vs gamma_xy ---
    ax1 = axes[0, 0]
    ax1.plot(fort[:, 1]*100, fort[:, 2], 'b-o', markersize=2, label='Fortran (forumat.f90)')
    ax1.plot(theory[:, 1]*100, theory[:, 2], 'r--x', markersize=3, label='Theory (ADINA algorithm)')
    # Add elastic reference
    gamma_range = np.linspace(0, 0.005, 100)
    ax1.plot(gamma_range*100, G0*gamma_range, 'g:', alpha=0.5, label=f'Elastic: G₀·γ (G₀={G0:.0f})')
    ax1.set_xlabel('γ_xy (%)')
    ax1.set_ylabel('τ_xy (MPa)')
    ax1.set_title('Shear Stress vs Shear Strain')
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)
    
    # --- Plot 2: sigma_xx vs gamma_xy ---
    ax2 = axes[0, 1]
    ax2.plot(fort[:, 1]*100, fort[:, 3], 'b-o', markersize=2, label='σ_xx (Fortran)')
    ax2.plot(theory[:, 1]*100, theory[:, 3], 'r--x', markersize=3, label='σ_xx (Theory)')
    ax2.plot(fort[:, 1]*100, fort[:, 4], 'b-.', alpha=0.5, markersize=1, label='σ_yy (Fortran)')
    ax2.set_xlabel('γ_xy (%)')
    ax2.set_ylabel('σ (MPa)')
    ax2.set_title('Normal Stresses vs Shear Strain')
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)
    
    # --- Plot 3: Error tau ---
    ax3 = axes[1, 0]
    n = min(len(fort), len(theory))
    err_tau = fort[:n, 2] - theory[:n, 2]
    err_sx = fort[:n, 3] - theory[:n, 3]
    ax3.plot(fort[:n, 1]*100, err_tau, 'b-o', markersize=2, label='Δτ_xy')
    ax3.plot(fort[:n, 1]*100, err_sx, 'r-s', markersize=2, label='Δσ_xx')
    ax3.set_xlabel('γ_xy (%)')
    ax3.set_ylabel('Error (MPa)')
    ax3.set_title('Error: Fortran - Theory')
    ax3.legend(fontsize=8)
    ax3.grid(True, alpha=0.3)
    
    # --- Plot 4: Phase diagram and key points ---
    ax4 = axes[1, 1]
    
    # Color-code by phase
    phase_colors = {'elastic': 'green', 'cracking': 'orange', 'cracked': 'blue', 
                    'crushing': 'red', 'crushed': 'gray'}
    for phase_name, color in phase_colors.items():
        mask = [p == phase_name for p in phases]
        if any(mask):
            idx = [i for i, m in enumerate(mask) if m]
            ax4.scatter(theory[idx, 1]*100, theory[idx, 2], c=color, s=15, 
                       label=f'{phase_name}', zorder=5)
    
    ax4.plot(fort[:, 1]*100, fort[:, 2], 'k-', alpha=0.3, linewidth=1)
    ax4.set_xlabel('γ_xy (%)')
    ax4.set_ylabel('τ_xy (MPa)')
    ax4.set_title('Theoretical Phases')
    ax4.legend(fontsize=8)
    ax4.grid(True, alpha=0.3)
    
    # Annotate key points
    # Find cracking step
    for i, p in enumerate(phases):
        if p == 'cracking':
            ax4.annotate(f'Crack: γ={theory[i,1]*100:.3f}%\nτ={theory[i,2]:.2f}',
                        xy=(theory[i,1]*100, theory[i,2]), fontsize=7,
                        xytext=(theory[i,1]*100+0.05, theory[i,2]+1),
                        arrowprops=dict(arrowstyle='->', color='orange'))
            break
    
    # Find crushing step
    for i, p in enumerate(phases):
        if p == 'crushing':
            ax4.annotate(f'Crush: γ={theory[i,1]*100:.3f}%\nτ={theory[i,2]:.2f}',
                        xy=(theory[i,1]*100, theory[i,2]), fontsize=7,
                        xytext=(theory[i,1]*100+0.05, theory[i,2]-2),
                        arrowprops=dict(arrowstyle='->', color='red'))
            break
    
    plt.tight_layout()
    plt.savefig('plot_shear.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: plot_shear.png")
    
    # ========== Detailed Comparison Table ==========
    print("\n" + "="*80)
    print(f"{'Step':>4}  {'gamma%':>8}  {'tau_Fort':>10}  {'tau_Theo':>10}  {'err':>10}  "
          f"{'sx_Fort':>10}  {'sx_Theo':>10}  {'err_sx':>10}  {'Phase':>10}")
    print("="*80)
    
    for i in range(min(len(fort), len(theory))):
        step = int(fort[i, 0])
        if step % 5 == 0 or step <= 10 or (step >= 65 and step <= 70):
            err_t = fort[i, 2] - theory[i, 2]
            err_s = fort[i, 3] - theory[i, 3]
            print(f"{step:4d}  {fort[i,1]*100:8.4f}  {fort[i,2]:10.4f}  {theory[i,2]:10.4f}  "
                  f"{err_t:10.4f}  {fort[i,3]:10.4f}  {theory[i,3]:10.4f}  {err_s:10.4f}  "
                  f"{phases[i]:>10}")
    
    # RMSE for non-zero steps
    mask = fort[:n, 2] != 0
    if mask.any():
        rmse_tau = np.sqrt(np.mean((fort[:n, 2][mask] - theory[:n, 2][mask])**2))
        rmse_sx = np.sqrt(np.mean((fort[:n, 3][mask] - theory[:n, 3][mask])**2))
        max_tau = np.max(np.abs(fort[:n, 2][mask] - theory[:n, 2][mask]))
        print(f"\nNon-zero region RMSE: τ_xy = {rmse_tau:.6f} MPa, σ_xx = {rmse_sx:.6f} MPa")
        print(f"Max absolute error: τ_xy = {max_tau:.6f} MPa")
    
    # Summary
    print("\n--- Key Events ---")
    for i, p in enumerate(phases):
        if p == 'cracking' and (i == 0 or phases[i-1] != 'cracking'):
            print(f"  Cracking at step {int(theory[i,0])}: γ = {theory[i,1]*100:.4f}%, "
                  f"τ = {theory[i,2]:.4f} MPa")
            print(f"    Theory: cracking when τ > ft*(1-τ/|fc|) → τ_cr = {ft*abs(fc)/(abs(fc)+ft):.4f} MPa")
        if p == 'crushing' and (i == 0 or phases[i-1] != 'crushing'):
            print(f"  Crushing at step {int(theory[i,0])}: γ = {theory[i,1]*100:.4f}%, "
                  f"τ = {theory[i,2]:.4f} MPa")
        if p == 'crushed' and (i == 0 or phases[i-1] != 'crushed'):
            print(f"  First zero-stress step: {int(theory[i,0])}")


if __name__ == '__main__':
    main()
