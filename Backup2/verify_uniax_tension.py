# -*- coding: utf-8 -*-
"""
单轴拉伸：对比 forumat.f90 与理论应力-应变曲线。
理论：弹性段 σ=E·ε 至 ft，软化段线性降至 0（至 10×ε_peak）。
材料参数与 test_uniax_tension_30steps.f90 一致。
"""
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'SimSun', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

E = 30.0e9   # Pa
ft = 3.0e6   # 抗拉强度
eps_peak = ft / E   # 峰值拉应变


def theory_tension_stress(eps):
    """理论单轴拉应力：弹性至 ft，随后线性软化至 0（在 10×ε_peak）"""
    eps = np.asarray(eps)
    sigma = np.zeros_like(eps)
    for i in range(len(eps)):
        e = eps[i]
        if e <= 0:
            sigma[i] = 0.0
        elif e <= eps_peak:
            sigma[i] = E * e
        elif e <= 10.0 * eps_peak:
            sigma[i] = ft * (1.0 - (e - eps_peak) / (9.0 * eps_peak))
        else:
            sigma[i] = 0.0
    return sigma


def main():
    base = Path(__file__).parent

    # 1) 加载 forumat 单轴拉伸 30 步结果
    strain_f = np.array([])
    stress_f = np.array([])
    try:
        data = np.loadtxt(base / 'uniax_tension_30steps.txt', skiprows=1)
        strain_f = data[:, 1]
        stress_f = data[:, 2]
    except Exception as e:
        print('未找到 uniax_tension_30steps.txt，请先运行: test_uniax_tension_30steps.exe')
        print(e)

    # 2) 理论曲线：光滑曲线，应变从 0 到略超过 10×ε_peak
    strain_theory = np.linspace(0, 10.5 * eps_peak, 300)
    stress_theory = theory_tension_stress(strain_theory)

    # 3) 绘图
    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    ax.plot(strain_theory * 1e6, stress_theory / 1e6, 'b-', linewidth=2, label='理论 (弹性+线性软化)')
    if len(strain_f) > 0:
        ax.plot(strain_f * 1e6, stress_f / 1e6, 'ro', markersize=6, markeredgecolor='darkred',
                markeredgewidth=0.8, label='forumat.f90 (30步)')
    ax.set_xlabel(r'轴向拉应变 $\varepsilon_{xx}$ (×10$^{-6}$)', fontsize=12)
    ax.set_ylabel(r'轴向拉应力 $\sigma_{xx}$ (MPa)', fontsize=12)
    ax.set_title('单轴拉伸：理论曲线与 forumat.f90 计算曲线对比', fontsize=13)
    ax.grid(True, alpha=0.4)
    ax.axhline(0, color='k', linewidth=0.5)
    ax.axvline(0, color='k', linewidth=0.5)
    ax.legend(loc='upper right', fontsize=10)
    ax.set_xlim(left=0, right=strain_theory[-1] * 1e6 * 1.02)
    ax.set_ylim(bottom=0, top=ft / 1e6 * 1.1)
    plt.tight_layout()

    out = base / 'verify_uniax_tension.png'
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print('已保存:', out)

    if len(strain_f) > 0:
        stress_theory_at_f = theory_tension_stress(strain_f)
        err = np.abs(stress_f - stress_theory_at_f)
        mask = np.abs(stress_theory_at_f) > 1e3
        rel = np.zeros_like(err)
        np.place(rel, mask, 100.0 * err[mask] / np.abs(stress_theory_at_f[mask]))
        print('单轴拉伸 30 步对比:')
        print('  最大绝对误差 (Pa):', np.max(err))
        print('  最大相对误差 (%):', np.max(rel) if np.any(mask) else 0)


if __name__ == '__main__':
    main()
