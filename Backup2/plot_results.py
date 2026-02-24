# -*- coding: utf-8 -*-
"""绘制 ADINA 平面应力混凝土本构测试结果"""
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# 使用支持中文的字体（Windows 常见）
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'SimSun', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def load_data(fpath):
    """读取两列数据，跳过标题行"""
    data = np.loadtxt(fpath, skiprows=1)
    return data[:, 0], data[:, 1]

def main():
    base = Path(__file__).parent
    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    axes = axes.flatten()

    # 1. 单轴拉伸
    x, y = load_data(base / 'test_adina_uniax_tension.txt')
    ax = axes[0]
    ax.plot(x * 1e3, y / 1e6, 'b-', linewidth=1.5)
    ax.set_xlabel(r'应变 $\varepsilon_{xx}$ (×10$^{-3}$)')
    ax.set_ylabel(r'应力 $\sigma_{xx}$ (MPa)')
    ax.set_title('单轴拉伸')
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='k', linewidth=0.5)
    ax.axvline(0, color='k', linewidth=0.5)

    # 2. 单轴压缩
    x, y = load_data(base / 'test_adina_uniax_compression.txt')
    ax = axes[1]
    ax.plot(x * 1e3, y / 1e6, 'r-', linewidth=1.5)
    ax.set_xlabel(r'应变 $\varepsilon_{xx}$ (×10$^{-3}$)')
    ax.set_ylabel(r'应力 $\sigma_{xx}$ (MPa)')
    ax.set_title('单轴压缩')
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='k', linewidth=0.5)
    ax.axvline(0, color='k', linewidth=0.5)

    # 3. 纯剪切
    x, y = load_data(base / 'test_adina_shear.txt')
    ax = axes[2]
    ax.plot(x * 1e3, y / 1e6, 'g-', linewidth=1.5)
    ax.set_xlabel(r'剪应变 $\gamma_{xy}$ (×10$^{-3}$)')
    ax.set_ylabel(r'剪应力 $\tau_{xy}$ (MPa)')
    ax.set_title('纯剪切')
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='k', linewidth=0.5)
    ax.axvline(0, color='k', linewidth=0.5)

    # 4. 往复单轴
    x, y = load_data(base / 'test_adina_cyclic_uniax.txt')
    ax = axes[3]
    ax.plot(x * 1e3, y / 1e6, 'b-', linewidth=1)
    ax.set_xlabel(r'应变 $\varepsilon_{xx}$ (×10$^{-3}$)')
    ax.set_ylabel(r'应力 $\sigma_{xx}$ (MPa)')
    ax.set_title('往复拉压')
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='k', linewidth=0.5)
    ax.axvline(0, color='k', linewidth=0.5)

    # 5. 往复剪切
    x, y = load_data(base / 'test_adina_cyclic_shear.txt')
    ax = axes[4]
    ax.plot(x * 1e3, y / 1e6, 'g-', linewidth=1)
    ax.set_xlabel(r'剪应变 $\gamma_{xy}$ (×10$^{-3}$)')
    ax.set_ylabel(r'剪应力 $\tau_{xy}$ (MPa)')
    ax.set_title('往复剪切')
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='k', linewidth=0.5)
    ax.axvline(0, color='k', linewidth=0.5)

    axes[5].axis('off')

    plt.tight_layout()
    out_png = base / 'adina_concrete_results.png'
    out_pdf = base / 'adina_concrete_results.pdf'
    plt.savefig(out_png, dpi=150, bbox_inches='tight')
    plt.savefig(out_pdf, bbox_inches='tight')
    print('已保存:', out_png, out_pdf)
    plt.show()

if __name__ == '__main__':
    main()
