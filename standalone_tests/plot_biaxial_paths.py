"""
plot_biaxial_paths.py
Plot biaxial stress-path test results:
  Panel 1: σ₁/fc vs ε₁ for all 8 paths
  Panel 2: σ₂/fc vs ε₂ for all 8 paths
  Panel 3: Biaxial failure envelope (σ₂/fc vs σ₁/fc) with peak-stress points
  Panel 4: Stress path in σ₁-σ₂ space (full trajectories)
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os

rcParams['font.size'] = 10
rcParams['figure.dpi'] = 150

csv_path = os.path.join(os.path.dirname(__file__), 'out_biaxial_paths.csv')
data = np.genfromtxt(csv_path, delimiter=',', skip_header=1)

path_ids = sorted(set(data[:, 0].astype(int)))
labels = {
    1: r'$-1:-1$',
    2: r'$-1:-0.5$',
    3: r'$-1:-0.2$',
    4: r'$-1:0$',
    5: r'$-1:0.1$',
    6: r'$-1:0.2$',
    7: r'$-1:0.5$',
    8: r'$1:1$',
}

cmap = plt.cm.tab10
colors = {pid: cmap(i / max(len(path_ids)-1, 1)) for i, pid in enumerate(path_ids)}

fc = -30.0

fig, axes = plt.subplots(2, 2, figsize=(14, 11))
fig.suptitle('Biaxial Stress Path Test Results', fontsize=14, fontweight='bold')

# ------------------------------------------------------------------
# Panel 1: σ₁/fc  vs  ε₁
# ------------------------------------------------------------------
ax1 = axes[0, 0]
for pid in path_ids:
    mask = data[:, 0].astype(int) == pid
    eps1 = data[mask, 2]
    sig1_fc = data[mask, 7]
    ax1.plot(eps1 * 1000, sig1_fc, color=colors[pid], label=labels[pid], linewidth=1.2)
ax1.set_xlabel(r'$\varepsilon_1$ (‰)')
ax1.set_ylabel(r'$\sigma_1 / f_c$')
ax1.set_title(r'(a) $\sigma_1 / f_c$ vs $\varepsilon_1$')
ax1.legend(fontsize=8, ncol=2)
ax1.grid(True, alpha=0.3)
ax1.axhline(0, color='k', linewidth=0.5)
ax1.axvline(0, color='k', linewidth=0.5)

# ------------------------------------------------------------------
# Panel 2: σ₂/fc  vs  ε₂
# ------------------------------------------------------------------
ax2 = axes[0, 1]
for pid in path_ids:
    mask = data[:, 0].astype(int) == pid
    eps2 = data[mask, 3]
    sig2_fc = data[mask, 8]
    ax2.plot(eps2 * 1000, sig2_fc, color=colors[pid], label=labels[pid], linewidth=1.2)
ax2.set_xlabel(r'$\varepsilon_2$ (‰)')
ax2.set_ylabel(r'$\sigma_2 / f_c$')
ax2.set_title(r'(b) $\sigma_2 / f_c$ vs $\varepsilon_2$')
ax2.legend(fontsize=8, ncol=2)
ax2.grid(True, alpha=0.3)
ax2.axhline(0, color='k', linewidth=0.5)
ax2.axvline(0, color='k', linewidth=0.5)

# ------------------------------------------------------------------
# Panel 3: Biaxial Failure Envelope - peak stress points
# ------------------------------------------------------------------
ax3 = axes[1, 0]

peak_sig1_fc = []
peak_sig2_fc = []
for pid in path_ids:
    mask = data[:, 0].astype(int) == pid
    s1 = data[mask, 4]
    s2 = data[mask, 5]
    s1_fc = data[mask, 7]
    s2_fc = data[mask, 8]

    if pid <= 7:
        idx_peak = np.argmin(s1)
    else:
        idx_peak = np.argmax(s1)

    peak_sig1_fc.append(s1_fc[idx_peak])
    peak_sig2_fc.append(s2_fc[idx_peak])

    ax3.plot(s1_fc[idx_peak], s2_fc[idx_peak], 'o', color=colors[pid],
             markersize=8, zorder=5)
    ax3.annotate(labels[pid], (s1_fc[idx_peak], s2_fc[idx_peak]),
                 textcoords='offset points', xytext=(6, 6), fontsize=7,
                 color=colors[pid])

peak_sig1_fc = np.array(peak_sig1_fc)
peak_sig2_fc = np.array(peak_sig2_fc)
comp_mask = peak_sig1_fc > 0.1
if np.sum(comp_mask) > 1:
    order = np.argsort(peak_sig1_fc[comp_mask])
    ax3.plot(peak_sig1_fc[comp_mask][order], peak_sig2_fc[comp_mask][order],
             '--', color='gray', linewidth=1, alpha=0.7)

ax3.set_xlabel(r'$\sigma_1 / f_c$')
ax3.set_ylabel(r'$\sigma_2 / f_c$')
ax3.set_title('(c) Biaxial Failure Envelope (Peak Stress)')
ax3.grid(True, alpha=0.3)
ax3.axhline(0, color='k', linewidth=0.5)
ax3.axvline(0, color='k', linewidth=0.5)
ax3.set_aspect('equal')
ax3.plot([0, 1.5], [0, 1.5], ':', color='gray', linewidth=0.5, alpha=0.5)

# ------------------------------------------------------------------
# Panel 4: Full stress trajectories σ₂/fc vs σ₁/fc
# ------------------------------------------------------------------
ax4 = axes[1, 1]
for pid in path_ids:
    mask = data[:, 0].astype(int) == pid
    s1_fc = data[mask, 7]
    s2_fc = data[mask, 8]
    ax4.plot(s1_fc, s2_fc, color=colors[pid], label=labels[pid], linewidth=1.0, alpha=0.8)

ax4.set_xlabel(r'$\sigma_1 / f_c$')
ax4.set_ylabel(r'$\sigma_2 / f_c$')
ax4.set_title(r'(d) Stress Trajectories in $\sigma_1/f_c - \sigma_2/f_c$ Space')
ax4.legend(fontsize=8, ncol=2)
ax4.grid(True, alpha=0.3)
ax4.axhline(0, color='k', linewidth=0.5)
ax4.axvline(0, color='k', linewidth=0.5)

plt.tight_layout()
outfile = os.path.join(os.path.dirname(__file__), 'biaxial_paths_report.png')
plt.savefig(outfile, dpi=200, bbox_inches='tight')
print(f'Saved: {outfile}')
plt.close()

# ------------------------------------------------------------------
# Summary table
# ------------------------------------------------------------------
print('\n=== Biaxial Path Peak Stress Summary ===')
print(f'{"Path":<12} {"sig1_peak/fc":>13} {"sig2_peak/fc":>13} {"sig1_peak(MPa)":>15} {"sig2_peak(MPa)":>15}')
for i, pid in enumerate(path_ids):
    mask = data[:, 0].astype(int) == pid
    s1 = data[mask, 4]
    s2 = data[mask, 5]
    if pid <= 7:
        idx = np.argmin(s1)
    else:
        idx = np.argmax(s1)
    print(f'{labels[pid]:<12} {s1[idx]/fc:>13.4f} {s2[idx]/fc:>13.4f} '
          f'{s1[idx]:>15.4f} {s2[idx]:>15.4f}')
