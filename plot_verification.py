import matplotlib.pyplot as plt
import numpy as np
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))
plt.rcParams['font.size'] = 10

fig, axes = plt.subplots(2, 3, figsize=(16, 10))

# 1. Compression curve
d = np.genfromtxt('curve_compression.csv', delimiter=',', skip_header=1)
ax = axes[0,0]
ax.plot(d[:,0]*1e3, d[:,1]/1e6, 'b-', linewidth=1.5)
ax.set_xlabel('Strain (‰)')
ax.set_ylabel('Stress (MPa)')
ax.set_title('Uniaxial Compression')
ax.axhline(0, color='k', linewidth=0.5)
ax.grid(True, alpha=0.3)

# 2. Tension softening (detailed)
d = np.genfromtxt('curve_tension_detail.csv', delimiter=',', skip_header=1)
ax = axes[0,1]
ax.plot(d[:,0]*1e3, d[:,1]/1e6, 'r-', linewidth=1.5)
ax.set_xlabel('Strain (‰)')
ax.set_ylabel('Stress (MPa)')
ax.set_title('Tension Softening')
ax.axhline(2.07, color='gray', linestyle='--', linewidth=0.5, label='ft')
ax.grid(True, alpha=0.3)
ax.legend()

# 3. Cyclic compression
d = np.genfromtxt('curve_cyclic.csv', delimiter=',', skip_header=1)
ax = axes[0,2]
ax.plot(d[:,0]*1e3, d[:,1]/1e6, 'g-', linewidth=1.5)
ax.set_xlabel('Strain (‰)')
ax.set_ylabel('Stress (MPa)')
ax.set_title('Cyclic Compression-Tension')
ax.axhline(0, color='k', linewidth=0.5)
ax.grid(True, alpha=0.3)

# 4. Crack cyclic (crack → close → reopen)
d = np.genfromtxt('curve_crack_cyclic.csv', delimiter=',', skip_header=1)
ax = axes[1,0]
ax.plot(d[:,0]*1e3, d[:,1]/1e6, 'm-', linewidth=1.5)
ax.set_xlabel('Strain (‰)')
ax.set_ylabel('Stress (MPa)')
ax.set_title('Crack Open → Close → Reopen')
ax.axhline(2.07, color='gray', linestyle='--', linewidth=0.5, label='ft')
ax.axhline(-20.7, color='gray', linestyle=':', linewidth=0.5, label='fc')
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8)

# 5. Shear wall: tau vs gamma
d = np.genfromtxt('curve_shear_wall.csv', delimiter=',', skip_header=1)
ax = axes[1,1]
ax.plot(d[:,3]*1e3, d[:,6]/1e6, 'c-', linewidth=1.0)
ax.set_xlabel('Shear Strain γ (‰)')
ax.set_ylabel('Shear Stress τ (MPa)')
ax.set_title('Shear Wall: τ-γ Response')
ax.grid(True, alpha=0.3)

# 6. Tension softening tangent
d = np.genfromtxt('curve_tension_detail.csv', delimiter=',', skip_header=1)
ax = axes[1,2]
ax.plot(d[:,0]*1e3, d[:,3]/1e9, 'k-', linewidth=1.5)
ax.set_xlabel('Strain (‰)')
ax.set_ylabel('C11 (GPa)')
ax.set_title('Tangent C11 in Tension')
ax.axhline(0, color='r', linestyle='--', linewidth=0.5)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('verification_plots.png', dpi=150)
plt.close()
print('Saved: verification_plots.png')
