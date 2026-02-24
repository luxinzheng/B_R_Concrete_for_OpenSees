# -*- coding: utf-8 -*-
"""
单轴压缩：对比 forumat.f90 与 ADINA 理论 Saenz 应力-应变曲线。
ADINA 切线模量 (all_adina.for 25184-25186, forumat 318-320):
  DE = |应变|/EPSCP (受压为正)
  Et = E * (1 - RBM5*DE^2 - 2*RCM5*DE^3) / (1 + RAM5*DE + RBM5*DE^2 + RCM5*DE^3)^2
应力曲线由 Et 对应变积分得到。
材料参数需与 test_adina_concrete.f90 中单轴压缩一致。
"""
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'SimSun', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 与 test_adina_concrete.f90 完全一致
E     = 30.0e9   # Pa
nu    = 0.2
ft    = 3.0e6
fc    = 30.0e6   # 抗压强度
epsc  = 0.002    # 峰值压应变
fu    = 0.85     # 极限应力/峰值
epsu  = 0.004    # 极限压应变
# 单轴时 EPSCP=epsc, SIGCP=fc
EPSCP = epsc
SIGCP = fc
RP    = epsu / epsc
ES    = SIGCP / EPSCP
# EU = (fu*SIGCP)/(epsu*EPSCP/epsc) = fu*epsc/epsu (当 SIGCP=fc, EPSCP=epsc)
EU    = (fu * SIGCP) / (epsu * EPSCP / epsc)
RAM5  = E/EU + (RP - 2)*RP*RP*E/ES - (2*RP + 1)*(RP - 1)**2
RAM5  = RAM5 / (RP * (RP - 1)**2)
RBM5  = 2*E/ES - 3 - 2*RAM5
RCM5  = 2 - E/ES + RAM5


def saenz_tangent(DE):
    """Saenz 切线模量 Et(DE)，DE = 应变/EPSCP（受压应变取绝对值，DE>0 表示压）"""
    if DE <= 0:
        return E
    den = 1 + RAM5*DE + RBM5*DE**2 + RCM5*DE**3
    if abs(den) < 1e-30:
        return -abs(E)  # 下降段负刚度
    num = E * (1 - RBM5*DE**2 - 2*RCM5*DE**3)
    Et = num / den**2
    return max(Et, -abs(E))


def adina_uniax_curve(eps_array):
    """对应变序列积分 Et 得到 ADINA 理论单轴压应力 (eps 为负表示压)"""
    eps = np.asarray(eps_array)
    DE = -eps / EPSCP  # 压应变 -> DE > 0
    Et = np.array([saenz_tangent(d) for d in DE])
    sigma = np.zeros_like(eps)
    for i in range(1, len(eps)):
        d_eps = eps[i] - eps[i-1]
        Et_avg = 0.5 * (Et[i] + Et[i-1])
        sigma[i] = sigma[i-1] + Et_avg * d_eps
    return sigma


def main():
    base = Path(__file__).parent

    # 1) 加载 forumat 单轴压缩结果 (strain_xx, stress_xx)
    try:
        data = np.loadtxt(base / 'test_adina_uniax_compression.txt', skiprows=1)
        strain_forumat = data[:, 0]
        stress_forumat = data[:, 1]
    except Exception as e:
        print('未找到 test_adina_uniax_compression.txt，请先运行: test_adina_concrete.exe')
        print(e)
        strain_forumat = np.array([])
        stress_forumat = np.array([])

    # 1b) 加载 0→峰值 20 步离散点 (uniax_20steps.txt: step, strain_xx, stress_xx)
    strain_20 = np.array([])
    stress_20 = np.array([])
    try:
        data20 = np.loadtxt(base / 'uniax_20steps.txt', skiprows=1)
        strain_20 = data20[:, 1]
        stress_20 = data20[:, 2]
    except Exception:
        pass

    # 2) ADINA 理论曲线：从 0 开始积分，故应变序列需含 0；与 forumat 相同应变点便于逐点对比
    if len(strain_forumat) > 0:
        # 保证应变单调（加载顺序 0 -> 更负），理论积分从 0 开始
        strain_ref = np.concatenate([[0.0], np.sort(strain_forumat)[::-1]])
    else:
        strain_ref = np.linspace(0, -epsu * 1.05, 200)
    stress_ref = adina_uniax_curve(strain_ref)

    # 3) 绘图
    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    ax.plot(strain_ref * 1e3, -stress_ref / 1e6, 'b-', linewidth=2, label='ADINA 理论 (Saenz 积分)')
    if len(strain_forumat) > 0:
        ax.plot(strain_forumat * 1e3, -stress_forumat / 1e6, 'r.', markersize=3, label='forumat.f90 (PSUMAT)')
    if len(strain_20) > 0:
        ax.plot(strain_20 * 1e3, -stress_20 / 1e6, 'go', markersize=6, markeredgecolor='darkgreen', markeredgewidth=0.8, label='forumat 20步(0→峰值)')
    ax.set_xlabel(r'轴向压应变 $\varepsilon_{xx}$ (×10$^{-3}$)', fontsize=12)
    ax.set_ylabel(r'轴向压应力 $\sigma_{xx}$ (MPa)', fontsize=12)
    ax.set_title('单轴压缩：forumat.f90 与 ADINA 应力-应变曲线对比', fontsize=13)
    ax.grid(True, alpha=0.4)
    ax.axhline(0, color='k', linewidth=0.5)
    ax.axvline(0, color='k', linewidth=0.5)
    ax.legend(loc='lower left', fontsize=10)
    xmin = np.min(strain_ref)
    if len(strain_forumat) > 0:
        xmin = min(xmin, np.min(strain_forumat))
    if len(strain_20) > 0:
        xmin = min(xmin, np.min(strain_20))
    ax.set_xlim(left=xmin * 1e3 * 1.05, right=0)
    plt.tight_layout()

    out = base / 'verify_uniax_adina.png'
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print('Saved:', out)

    # 4) 若有 forumat 数据，输出最大误差 (在相同应变点插值参考曲线对比)
    if len(strain_forumat) > 0 and len(strain_ref) > 0:
        # np.interp 要求 xp 单调递增，strain_ref 为 0→负，故取反序
        stress_ref_at_forumat = np.interp(strain_forumat, strain_ref[::-1], stress_ref[::-1])
        err = np.abs(stress_forumat - stress_ref_at_forumat)
        mask = np.abs(stress_ref_at_forumat) > 1e6
        rel = np.zeros_like(err)
        np.place(rel, mask, err[mask] / np.abs(stress_ref_at_forumat[mask]))
        print('单轴压缩对比 (压应力取绝对值):')
        print('  最大绝对误差 (Pa):', np.max(err))
        print('  最大相对误差 (|σ|>1MPa):', np.max(rel) if np.any(mask) else 0)
        print('  峰值应力 ADINA 理论 (MPa):', -np.min(stress_ref) / 1e6)
        print('  峰值应力 forumat (MPa):', np.max(-stress_forumat) / 1e6)
        # 上升段 (应变 0~epsc) 相对误差
        rise = (strain_forumat >= -epsc * 1.01) & (strain_forumat <= 0)
        if np.sum(rise) > 0:
            rel_rise = rel[rise]
            rel_rise = rel_rise[np.abs(stress_ref_at_forumat[rise]) > 1e6]
            if len(rel_rise) > 0:
                print('  上升段 |ε|≤epsc 最大相对误差:', np.max(rel_rise))


if __name__ == '__main__':
    main()
