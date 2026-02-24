# -*- coding: utf-8 -*-
"""绘制单轴压缩应力-应变曲线"""
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'SimSun', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def main():
    base = Path(__file__).parent
    data = np.loadtxt(base / 'test_adina_uniax_compression.txt', skiprows=1)
    strain = data[:, 0]   # 应变 (无量纲)
    stress = data[:, 1]   # 应力 (Pa)

    fig, ax = plt.subplots(1, 1, figsize=(7, 5))
    ax.plot(strain * 1e3, stress / 1e6, 'r-', linewidth=2, label=r'$\sigma_{xx}$–$\varepsilon_{xx}$')
    ax.set_xlabel(r'轴向应变 $\varepsilon_{xx}$ (×10$^{-3}$)', fontsize=12)
    ax.set_ylabel(r'轴向应力 $\sigma_{xx}$ (MPa)', fontsize=12)
    ax.set_title('单轴压缩应力-应变曲线 (ADINA 平面应力混凝土)', fontsize=13)
    ax.grid(True, alpha=0.4)
    ax.axhline(0, color='k', linewidth=0.5)
    ax.axvline(0, color='k', linewidth=0.5)
    ax.legend(loc='lower left', fontsize=10)
    ax.set_xlim(left=min(strain) * 1e3 * 1.05, right=0)
    plt.tight_layout()

    out = base / 'uniax_compression_curve.png'
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print('Saved:', out)

if __name__ == '__main__':
    main()
