# -*- coding: utf-8 -*-
"""
剪切与双轴：forumat 与理论对比图。
需先运行 test_shear.exe + compare_shear.py，以及 test_biax_5paths.exe + compare_biax_5paths.py。
"""
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

BASE = Path(__file__).parent
MPA = 1e-6  # Pa -> MPa


def load_shear():
    data = np.loadtxt(BASE / 'compare_shear.txt', skiprows=2)
    step = data[:, 0]
    strain_xy = data[:, 1]
    theory_xx = data[:, 2] * MPA
    theory_yy = data[:, 3] * MPA
    theory_xy = data[:, 4] * MPA
    forumat_xx = data[:, 5] * MPA
    forumat_yy = data[:, 6] * MPA
    forumat_xy = data[:, 7] * MPA
    return strain_xy, theory_xx, theory_yy, theory_xy, forumat_xx, forumat_yy, forumat_xy


def load_biax():
    data = np.loadtxt(BASE / 'compare_biax_5paths.txt', skiprows=1)
    path_id = data[:, 0].astype(int)
    strain_xx = data[:, 2]
    strain_yy = data[:, 3]
    theory_xx = data[:, 4] * MPA
    theory_yy = data[:, 5] * MPA
    forumat_xx = data[:, 6] * MPA
    forumat_yy = data[:, 7] * MPA
    return path_id, strain_xx, strain_yy, theory_xx, theory_yy, forumat_xx, forumat_yy


def plot_shear():
    strain_xy, theory_xx, theory_yy, theory_xy, forumat_xx, forumat_yy, forumat_xy = load_shear()
    # 剪应变用百分数显示
    strain_pct = strain_xy * 100

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    fig.suptitle('纯剪切 0→0.5%（100 步）：理论 vs forumat', fontsize=12)

    axes[0].plot(strain_pct, theory_xx, 'b-', lw=2, label='理论')
    axes[0].plot(strain_pct, forumat_xx, 'r--', lw=1.5, label='forumat')
    axes[0].set_xlabel(r'剪应变 $\gamma/2$ (%)')
    axes[0].set_ylabel(r'$\sigma_{xx}$ (MPa)')
    axes[0].set_title(r'$\sigma_{xx}$')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(strain_pct, theory_yy, 'b-', lw=2, label='理论')
    axes[1].plot(strain_pct, forumat_yy, 'r--', lw=1.5, label='forumat')
    axes[1].set_xlabel(r'剪应变 $\gamma/2$ (%)')
    axes[1].set_ylabel(r'$\sigma_{yy}$ (MPa)')
    axes[1].set_title(r'$\sigma_{yy}$')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(strain_pct, theory_xy, 'b-', lw=2, label='理论')
    axes[2].plot(strain_pct, forumat_xy, 'r--', lw=1.5, label='forumat')
    axes[2].set_xlabel(r'剪应变 $\gamma/2$ (%)')
    axes[2].set_ylabel(r'$\sigma_{xy}$ (MPa)')
    axes[2].set_title(r'$\sigma_{xy}$')
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    out = BASE / 'plot_shear_compare.png'
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print('Saved:', out)


def plot_biax():
    path_id, strain_xx, strain_yy, theory_xx, theory_yy, forumat_xx, forumat_yy = load_biax()
    paths = [
        (1, '-1:-1', '双轴等压'),
        (2, '1:1', '双轴等拉'),
        (3, '-1:-0.5', '双轴压 -1:-0.5'),
        (4, '-1:0.1', '压拉 -1:0.1'),
        (5, '-1:0.2', '压拉 -1:0.2'),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    axes = axes.flatten()
    fig.suptitle('双轴 5 路径（应变至 1%）：理论 vs forumat', fontsize=12)

    for idx, (pid, ratio, title) in enumerate(paths):
        mask = path_id == pid
        ex = strain_xx[mask] * 100
        ey = strain_yy[mask] * 100
        tx = theory_xx[mask]
        ty = theory_yy[mask]
        fx = forumat_xx[mask]
        fy = forumat_yy[mask]

        ax = axes[idx]
        ax.plot(ex, tx, 'b-', lw=2, label=r'理论 $\sigma_{xx}$')
        ax.plot(ex, fx, 'b--', lw=1.5, label=r'forumat $\sigma_{xx}$')
        ax.plot(ey, ty, 'g-', lw=2, label=r'理论 $\sigma_{yy}$')
        ax.plot(ey, fy, 'g--', lw=1.5, label=r'forumat $\sigma_{yy}$')
        ax.set_xlabel(r'应变 (%)')
        ax.set_ylabel(r'应力 (MPa)')
        ax.set_title(f'路径{pid} {ratio}\n{title}')
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

    # 隐藏多余子图
    axes[5].set_visible(False)
    plt.tight_layout()
    out = BASE / 'plot_biax_compare.png'
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print('Saved:', out)


def main():
    if not (BASE / 'compare_shear.txt').exists():
        print('请先运行 test_shear.exe 与 compare_shear.py 生成 compare_shear.txt')
    else:
        plot_shear()

    if not (BASE / 'compare_biax_5paths.txt').exists():
        print('请先运行 test_biax_5paths.exe 与 compare_biax_5paths.py 生成 compare_biax_5paths.txt')
    else:
        plot_biax()


if __name__ == '__main__':
    main()
