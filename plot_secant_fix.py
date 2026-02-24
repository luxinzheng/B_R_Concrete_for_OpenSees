# -*- coding: utf-8 -*-
"""Plot secant stress correction verification results"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

plt.rcParams['font.family'] = 'Microsoft YaHei'
plt.rcParams['axes.unicode_minus'] = False

root = r'e:\Basic\Concrete_Model\ADINA'

def load_csv(name):
    path = os.path.join(root, name)
    return np.genfromtxt(path, delimiter=',', skip_header=1, filling_values=0.0)

fc = -2.58e7
ft = 2.0e6

fig = plt.figure(figsize=(20, 14))
fig.suptitle('割线应力修正验证 (Secant Stress Correction)\nforumat.f90 修改后的独立材料测试',
             fontsize=14, fontweight='bold')

# 1. Uniaxial cyclic: stress-strain hysteresis
d1 = load_csv('secant_uniaxial_cyclic.csv')
ax1 = fig.add_subplot(2, 3, 1)
ax1.plot(d1[:, 1]*1000, d1[:, 2]/1e6, 'r-', lw=0.8, alpha=0.7)
ax1.axhline(fc/1e6, color='b', ls='--', lw=0.8, alpha=0.5, label=f'fc={fc/1e6:.1f} MPa')
ax1.axhline(ft/1e6, color='g', ls='--', lw=0.8, alpha=0.5, label=f'ft={ft/1e6:.1f} MPa')
ax1.set_xlabel('应变 ε_xx (×10⁻³)')
ax1.set_ylabel('应力 σ_xx (MPa)')
ax1.set_title('单轴循环 20 周期\n应力-应变滞回')
ax1.legend(fontsize=8)
ax1.grid(True, alpha=0.3)

# 2. Uniaxial cyclic: stress history
ax2 = fig.add_subplot(2, 3, 2)
ax2.plot(d1[:, 0], d1[:, 2]/1e6, 'r-', lw=0.5, alpha=0.7, label='σ_xx')
ax2.axhline(fc/1e6, color='b', ls='--', lw=0.8, alpha=0.5)
ax2.axhline(ft/1e6, color='g', ls='--', lw=0.8, alpha=0.5)
ax2.set_xlabel('步号')
ax2.set_ylabel('σ_xx (MPa)')
ax2.set_title('单轴循环: 应力历程\n20 周期无漂移')
ax2.grid(True, alpha=0.3)

# 3. Crack cycle: stress-strain
d2 = load_csv('secant_crack_cycle.csv')
ax3 = fig.add_subplot(2, 3, 3)
ax3.plot(d2[:, 1]*1000, d2[:, 2]/1e6, 'r-', lw=0.5, alpha=0.7)
ax3.axhline(fc/1e6, color='b', ls='--', lw=0.8, alpha=0.5)
ax3.axhline(ft/1e6, color='g', ls='--', lw=0.8, alpha=0.5)
ax3.set_xlabel('应变 ε_xx (×10⁻³)')
ax3.set_ylabel('σ_xx (MPa)')
ax3.set_title('裂缝开/闭/重开 50 周期\n应力有界')
ax3.grid(True, alpha=0.3)

# 4. Crack cycle: stress history
ax4 = fig.add_subplot(2, 3, 4)
ax4.plot(d2[:, 0], d2[:, 2]/1e6, 'r-', lw=0.5, alpha=0.7, label='σ_xx')
ax4.plot(d2[:, 0], d2[:, 3]/1e6, 'b-', lw=0.5, alpha=0.5, label='σ_yy')
ax4.axhline(fc/1e6, color='k', ls=':', lw=0.5, alpha=0.3)
ax4.axhline(ft/1e6, color='k', ls=':', lw=0.5, alpha=0.3)
ax4.set_xlabel('步号')
ax4.set_ylabel('应力 (MPa)')
ax4.set_title('裂缝循环: σ 历程\n50 周期无漂移')
ax4.legend(fontsize=8)
ax4.grid(True, alpha=0.3)

# 5. ML scenario: stress-strain
d3 = load_csv('secant_ml_scenario.csv')
ax5 = fig.add_subplot(2, 3, 5)
ax5.plot(d3[:, 1]*1000, d3[:, 2]/1e6, 'r-', lw=0.5, alpha=0.7)
ax5.axhline(fc/1e6, color='b', ls='--', lw=0.8, alpha=0.5, label=f'fc')
ax5.axhline(ft/1e6, color='g', ls='--', lw=0.8, alpha=0.5, label=f'ft')
ax5.set_xlabel('应变 ε_xx (×10⁻³)')
ax5.set_ylabel('σ_xx (MPa)')
ax5.set_title('ML 场景 11 块循环\n应力-应变 (有剪切)')
ax5.legend(fontsize=8)
ax5.grid(True, alpha=0.3)

# 6. ML scenario: all stress histories
ax6 = fig.add_subplot(2, 3, 6)
ax6.plot(d3[:, 0], d3[:, 2]/1e6, 'r-', lw=0.5, alpha=0.7, label='σ_xx')
ax6.plot(d3[:, 0], d3[:, 3]/1e6, 'b-', lw=0.5, alpha=0.5, label='σ_yy')
ax6.plot(d3[:, 0], d3[:, 4]/1e6, 'g-', lw=0.5, alpha=0.5, label='τ_xy')
ax6.axhline(fc/1e6, color='k', ls=':', lw=0.5, alpha=0.3)
ax6.axhline(ft/1e6, color='k', ls=':', lw=0.5, alpha=0.3)
ax6.set_xlabel('步号')
ax6.set_ylabel('应力 (MPa)')
ax6.set_title('ML 场景: 全分量历程\n11 块循环无漂移')
ax6.legend(fontsize=8)
ax6.grid(True, alpha=0.3)

plt.tight_layout(rect=[0, 0, 1, 0.93])
out = os.path.join(root, 'verify_secant_fix.png')
plt.savefig(out, dpi=150)
print(f'Saved: {out}')

# Print stats
for name, d in [('Uniaxial', d1), ('CrackCycle', d2), ('ML', d3)]:
    smax = d[:, 2].max()/1e6
    smin = d[:, 2].min()/1e6
    print(f'{name}: sig_xx range [{smin:.1f}, {smax:.1f}] MPa, {len(d)} pts')
