"""
Plot Cyclic Hysteresis Test Results vs Theoretical Curves
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.size'] = 11
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ============================================================
# Saenz curve (theoretical)
# ============================================================
def saenz_params(E0=30000, fc=-30, eps_c=-0.002, fu=-6, eps_u=-0.005):
    RP = eps_u / eps_c
    ES = fc / eps_c
    EU = fu / eps_u
    A = (E0/EU + (RP-2)*RP**2*E0/ES - (2*RP+1)*(RP-1)**2) / (RP*(RP-1)**2)
    B = 2*E0/ES - 3 - 2*A
    C = 2 - E0/ES + A
    return A, B, C

def saenz_stress(eps, E0=30000, eps_c=-0.002):
    A, B, C = saenz_params()
    sigma = np.zeros_like(eps)
    for i, e in enumerate(eps):
        DE = e / eps_c
        if DE <= 0:
            sigma[i] = E0 * e  # elastic for tension
        else:
            denom = 1 + A*DE + B*DE**2 + C*DE**3
            sigma[i] = E0 * e / denom if abs(denom) > 1e-30 else 0
    return sigma

# ============================================================
# Read CSV data
# ============================================================
def read_hysteresis_csv(filename):
    steps, eps_xx, sig_xx, eps_yy = [], [], [], []
    with open(filename, 'r') as f:
        f.readline()
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(',')
            steps.append(int(parts[0]))
            vals = parts[1].split()
            eps_xx.append(float(vals[0]))
            sig_xx.append(float(vals[1]))
            eps_yy.append(float(vals[2]))
    return np.array(steps), np.array(eps_xx), np.array(sig_xx), np.array(eps_yy)

steps, eps, sig, eps_yy = read_hysteresis_csv('test_hysteresis.csv')

# Material parameters
E0 = 30000.0
fc = -30.0
eps_c = -0.002
fu = -6.0
eps_u = -0.005
ft = 3.0
eps_t = ft / E0  # cracking strain

# ============================================================
# Compute theoretical response
# ============================================================
def compute_theoretical(eps):
    """
    Theoretical uniaxial cyclic response (matches forumat.f90):
    - ALL compression uses Saenz total-strain approach:
      * Envelope: sigma = saenz(eps) when strain exceeds historical max
      * Unload/reload: secant to historical max (E_sec = saenz(eps_max)/eps_max)
        - The secant path is a straight line from origin to max point on Saenz curve
        - Transition from inside to envelope is continuous (no discontinuity)
        - Secant modulus naturally degrades from E0 to fc/eps_c
    - Tension: elastic up to ft, then zero (cracked)
    """
    n = len(eps)
    sig_th = np.zeros(n)
    
    cracked = False
    eps_max_comp = 0.0  # most negative strain reached (always <= 0)
    
    for i in range(n):
        e = eps[i]
        
        if e >= 0:
            # Tension side
            if cracked or abs(eps_max_comp) >= abs(eps_c):
                sig_th[i] = 0.0  # no tensile capacity after cracking or post-peak
            else:
                sig_th[i] = E0 * e
                if sig_th[i] > ft:
                    sig_th[i] = ft
                    cracked = True
        else:
            # Compression side — Saenz total-strain approach
            if e <= eps_max_comp:
                # On envelope (strain exceeds historical max)
                sig_th[i] = saenz_stress(np.array([e]))[0]
                eps_max_comp = e
            else:
                # Inside envelope: secant to historical max point
                sig_at_max = saenz_stress(np.array([eps_max_comp]))[0]
                if abs(eps_max_comp) > 1e-15:
                    E_sec = sig_at_max / eps_max_comp
                else:
                    E_sec = E0
                sig_th[i] = E_sec * e
    
    return sig_th

sig_th = compute_theoretical(eps)

# Saenz envelope for overlay
eps_env = np.linspace(0, -0.012, 500)
sig_env = saenz_stress(eps_env)

# ============================================================
# Figure 1: Full hysteresis loop
# ============================================================
fig1, ax1 = plt.subplots(figsize=(10, 7))

# Saenz envelope
ax1.plot(eps_env*100, sig_env, 'r--', linewidth=1.0, alpha=0.5, label='Saenz envelope')

# Computed response
ax1.plot(eps*100, sig, 'b-', linewidth=0.8, alpha=0.8, label='Model (computed)')

# Theoretical response
ax1.plot(eps*100, sig_th, 'k--', linewidth=0.6, alpha=0.4, label='Theory (secant model)')

ax1.axhline(y=0, color='k', linewidth=0.3)
ax1.axvline(x=0, color='k', linewidth=0.3)
ax1.set_xlabel('$\\varepsilon_{xx}$ (%)')
ax1.set_ylabel('$\\sigma_{xx}$ (MPa)')
ax1.set_title('Cyclic Hysteresis: Increasing Amplitude (0.05% to 1.0%)')
ax1.legend(loc='lower left', fontsize=9)
ax1.grid(True, alpha=0.2)
ax1.set_xlim([-1.1, 1.1])
ax1.set_ylim([-35, 5])

fig1.tight_layout()
fig1.savefig('plot_hysteresis_full.png', dpi=150, bbox_inches='tight')
print('Saved: plot_hysteresis_full.png')

# ============================================================
# Figure 2: Zoomed views (pre-peak, post-peak, detail)
# ============================================================
fig2, axes2 = plt.subplots(2, 2, figsize=(14, 11))

# --- Top Left: Pre-peak cycles (before crushing) ---
ax = axes2[0, 0]
# Find the step where post-peak starts (stress first reaches fc)
idx_peak = 0
for i in range(len(sig)):
    if sig[i] < fc * 0.99:
        idx_peak = i
        break

# Show up to a bit beyond the peak
idx_show = min(idx_peak + 200, len(eps))
ax.plot(eps[:idx_show]*100, sig[:idx_show], 'b-', linewidth=1.0, label='Computed')
ax.plot(eps[:idx_show]*100, sig_th[:idx_show], 'r--', linewidth=0.8, alpha=0.6, label='Theory')
ax.plot(eps_env*100, sig_env, 'k:', linewidth=0.8, alpha=0.3, label='Saenz envelope')
ax.set_xlabel('$\\varepsilon_{xx}$ (%)')
ax.set_ylabel('$\\sigma_{xx}$ (MPa)')
ax.set_title(f'Pre-peak region (up to step {idx_show})')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.2)
ax.set_xlim([-0.35, 0.35])
ax.set_ylim([-35, 5])

# --- Top Right: Post-peak region ---
ax = axes2[0, 1]
ax.plot(eps*100, sig, 'b-', linewidth=0.8, label='Computed')
ax.plot(eps*100, sig_th, 'r--', linewidth=0.6, alpha=0.6, label='Theory')
ax.plot(eps_env*100, sig_env, 'k:', linewidth=0.8, alpha=0.3, label='Saenz envelope')
ax.set_xlabel('$\\varepsilon_{xx}$ (%)')
ax.set_ylabel('$\\sigma_{xx}$ (MPa)')
ax.set_title('Post-peak region')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.2)
ax.set_xlim([-1.1, 0.6])
ax.set_ylim([-35, 2])

# --- Bottom Left: Stress difference ---
ax = axes2[1, 0]
diff = sig - sig_th
ax.plot(steps, diff, 'g-', linewidth=0.6, alpha=0.8)
ax.axhline(y=0, color='k', linewidth=0.3)
ax.set_xlabel('Load Step')
ax.set_ylabel('$\\sigma_{comp} - \\sigma_{theory}$ (MPa)')
ax.set_title('Stress Difference (Computed - Theory)')
ax.grid(True, alpha=0.2)

# Statistics
rmse = np.sqrt(np.mean(diff**2))
max_err = np.max(np.abs(diff))
ax.text(0.95, 0.95, f'RMSE = {rmse:.3f} MPa\nMax |error| = {max_err:.3f} MPa',
        transform=ax.transAxes, ha='right', va='top', fontsize=9,
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# --- Bottom Right: Secant modulus evolution ---
ax = axes2[1, 1]

# Track secant modulus at each compressive reversal
# Find reversal points (compression peaks)
reversals_eps = []
reversals_sig = []
for i in range(1, len(eps)-1):
    if eps[i] < eps[i-1] and eps[i] < eps[i+1] and eps[i] < -0.001:
        reversals_eps.append(eps[i])
        reversals_sig.append(sig[i])

reversals_eps = np.array(reversals_eps)
reversals_sig = np.array(reversals_sig)

if len(reversals_eps) > 0:
    E_sec_comp = reversals_sig / reversals_eps  # secant modulus (positive)
    E_sec_th = saenz_stress(reversals_eps) / reversals_eps
    
    ax.plot(reversals_eps*100, E_sec_comp/1000, 'bo-', linewidth=1.2,
            markersize=4, label='Computed $E_{sec}$')
    ax.plot(reversals_eps*100, E_sec_th/1000, 'r^--', linewidth=1.0,
            markersize=4, alpha=0.7, label='Theory $E_{sec}$')
    ax.axhline(y=E0/1000, color='gray', linestyle=':', linewidth=0.8,
               alpha=0.5, label=f'$E_0 = {E0/1000:.0f}$ GPa')

ax.set_xlabel('Peak compressive strain (%)')
ax.set_ylabel('Secant Modulus (GPa)')
ax.set_title('Secant Modulus at Compression Reversals')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.2)
ax.set_ylim([0, 35])

fig2.tight_layout()
fig2.savefig('plot_hysteresis_detail.png', dpi=150, bbox_inches='tight')
print('Saved: plot_hysteresis_detail.png')

# ============================================================
# Figure 3: Stress history and strain history
# ============================================================
fig3, axes3 = plt.subplots(2, 1, figsize=(14, 8))

ax = axes3[0]
ax.plot(steps, eps*100, 'b-', linewidth=0.6)
ax.set_xlabel('Load Step')
ax.set_ylabel('$\\varepsilon_{xx}$ (%)')
ax.set_title('Strain History (Increasing Amplitude Protocol)')
ax.grid(True, alpha=0.2)
ax.axhline(y=0, color='k', linewidth=0.3)

# Mark amplitude levels
for i in range(1, 21):
    ax.axhline(y=i*0.05, color='gray', linewidth=0.3, alpha=0.3)
    ax.axhline(y=-i*0.05, color='gray', linewidth=0.3, alpha=0.3)

ax = axes3[1]
ax.plot(steps, sig, 'b-', linewidth=0.6, label='Computed')
ax.plot(steps, sig_th, 'r--', linewidth=0.5, alpha=0.5, label='Theory')
ax.axhline(y=fc, color='gray', linestyle=':', linewidth=0.8, alpha=0.5, label=f'$f_c = {fc}$ MPa')
ax.axhline(y=0, color='k', linewidth=0.3)
ax.set_xlabel('Load Step')
ax.set_ylabel('$\\sigma_{xx}$ (MPa)')
ax.set_title('Stress History')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.2)

fig3.tight_layout()
fig3.savefig('plot_hysteresis_history.png', dpi=150, bbox_inches='tight')
print('Saved: plot_hysteresis_history.png')

print('\nAll hysteresis plots generated successfully.')
