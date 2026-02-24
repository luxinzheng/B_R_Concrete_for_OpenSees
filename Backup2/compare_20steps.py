# -*- coding: utf-8 -*-
"""
单轴压缩 0→峰值 20 步：逐步对比 forumat 与 ADINA 理论应力。
- 读 uniax_20steps.txt（由 test_uniax_20steps.exe 生成）
- 对相同应变点用 Saenz 积分算理论应力
- 输出：step, strain, stress_theory, stress_forumat, err_abs, err_rel(%)，并写表到 compare_20steps.txt
"""
import numpy as np
from pathlib import Path

# 与 test_adina_concrete / test_uniax_20steps 一致
E     = 30.0e9
nu    = 0.2
ft    = 3.0e6
fc    = 30.0e6
epsc  = 0.002
fu    = 0.85
epsu  = 0.004
EPSCP = epsc
SIGCP = fc
RP    = epsu / epsc
ES    = SIGCP / EPSCP
EU    = (fu * SIGCP) / (epsu * EPSCP / epsc)
RAM5  = E/EU + (RP - 2)*RP*RP*E/ES - (2*RP + 1)*(RP - 1)**2
RAM5  = RAM5 / (RP * (RP - 1)**2)
RBM5  = 2*E/ES - 3 - 2*RAM5
RCM5  = 2 - E/ES + RAM5


def saenz_tangent(DE):
    if DE <= 0:
        return E
    den = 1 + RAM5*DE + RBM5*DE**2 + RCM5*DE**3
    if abs(den) < 1e-30:
        return -abs(E)
    num = E * (1 - RBM5*DE**2 - 2*RCM5*DE**3)
    Et = num / den**2
    return max(Et, -abs(E))


def adina_uniax_curve(eps_array):
    """对应变序列积分 Et 得到 ADINA 理论单轴压应力 (eps 为负表示压)"""
    eps = np.asarray(eps_array)
    DE = -eps / EPSCP
    Et = np.array([saenz_tangent(d) for d in DE])
    sigma = np.zeros_like(eps)
    for i in range(1, len(eps)):
        d_eps = eps[i] - eps[i-1]
        Et_avg = 0.5 * (Et[i] + Et[i-1])
        sigma[i] = sigma[i-1] + Et_avg * d_eps
    return sigma


def main():
    base = Path(__file__).parent
    nstep = 20

    # 1) 读 forumat 20 步结果
    try:
        data = np.loadtxt(base / 'uniax_20steps.txt', skiprows=1)
        steps = data[:, 0].astype(int)
        strain_f = data[:, 1]
        stress_f = data[:, 2]
    except Exception as e:
        print('请先运行: test_uniax_20steps.exe 生成 uniax_20steps.txt')
        print(e)
        return

    if len(steps) != nstep:
        print('期望 20 步，当前行数:', len(steps))
        return

    # 2) 理论应变点：0, -epsc/20, -2*epsc/20, ..., -epsc（与 Fortran 一致）
    strain_theory = np.array([0.0] + [ -epsc * i / nstep for i in range(1, nstep + 1) ])
    stress_theory = adina_uniax_curve(strain_theory)
    # 对比 forumat 的 step i 对应应变 strain_f[i-1] ≈ -i*epsc/20，对应理论应力 stress_theory[i]
    theory_at_steps = stress_theory[1 : nstep + 1]

    # 3) 误差
    err_abs = stress_f - theory_at_steps
    err_rel = np.zeros(nstep)
    for i in range(nstep):
        if abs(theory_at_steps[i]) > 1e3:
            err_rel[i] = 100.0 * err_abs[i] / abs(theory_at_steps[i])
        else:
            err_rel[i] = 0.0 if abs(err_abs[i]) < 1e-6 else 100.0

    # 4) 控制台输出
    print('单轴压缩 0→峰值 20 步对比 (应变×1000, 应力 MPa)')
    print('step  strain(1e-3)  theory(MPa)  forumat(MPa)  err_abs(Pa)   err_rel(%)')
    print('-' * 70)
    for i in range(nstep):
        print('{:4d}  {:12.6f}  {:11.4f}  {:11.4f}  {:12.2f}  {:10.2f}'.format(
            i + 1,
            strain_f[i] * 1e3,
            -theory_at_steps[i] / 1e6,
            -stress_f[i] / 1e6,
            err_abs[i],
            err_rel[i],
        ))
    print('-' * 70)
    print('最大绝对误差 (Pa):', np.max(np.abs(err_abs)))
    print('最大相对误差 (%):', np.max(np.abs(err_rel)))
    imax = np.argmax(np.abs(err_rel))
    print('最大相对误差步:', imax + 1, '  strain=', strain_f[imax]*1e3, '  rel%=', err_rel[imax])

    # 5) 写表
    out_table = base / 'compare_20steps.txt'
    with open(out_table, 'w', encoding='utf-8') as f:
        f.write('step  strain_xx       stress_theory   stress_forumat  err_abs(Pa)   err_rel_pct\n')
        for i in range(nstep):
            f.write('{:4d}  {:16.8e}  {:16.8e}  {:16.8e}  {:14.4f}  {:10.2f}\n'.format(
                i + 1, strain_f[i], theory_at_steps[i], stress_f[i], err_abs[i], err_rel[i]))
    print('已写入:', out_table)

    # 6) 简要诊断建议
    if np.argmax(np.abs(err_rel)) == 0:
        print('\n[建议] 第 1 步误差最大：检查首步主应变/主应力顺序与 ANG，单轴压轴向对应关系。')
    if np.max(np.abs(err_rel)) > 20 and np.argmax(np.abs(err_rel)) > 5:
        print('\n[建议] 中后段误差大：检查 YP/Et 是否按 Saenz、压应力限幅与开裂判断。')


if __name__ == '__main__':
    main()
