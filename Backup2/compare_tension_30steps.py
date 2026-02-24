# -*- coding: utf-8 -*-
"""
单轴拉伸 30 步：阶段1 为 0→峰值应变 10 步，阶段2 为 峰值→10×峰值 20 步。
理论：弹性段 sigma=E*eps 至 ft，软化段线性降至 0（至 10*eps_peak）。
逐步对比 forumat 与理论，输出差异表。
"""
import numpy as np
from pathlib import Path

E   = 30.0e9
ft  = 3.0e6
eps_peak = ft / E   # 峰值拉应变
n1, n2 = 10, 20


def theory_tension_stress(eps):
    """理论单轴拉应力：弹性至 ft，随后线性软化至 0（在 10*eps_peak）"""
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

    try:
        data = np.loadtxt(base / 'uniax_tension_30steps.txt', skiprows=1)
        steps = data[:, 0].astype(int)
        strain_f = data[:, 1]
        stress_f = data[:, 2]
        phase = data[:, 3].astype(int)
    except Exception as e:
        print('请先运行: test_uniax_tension_30steps.exe 生成 uniax_tension_30steps.txt')
        print(e)
        return

    nstep = n1 + n2
    if len(steps) != nstep:
        print('期望 30 步，当前行数:', len(steps))
        return

    stress_theory = theory_tension_stress(strain_f)
    err_abs = stress_f - stress_theory
    err_rel = np.zeros(nstep)
    for i in range(nstep):
        if abs(stress_theory[i]) > 1e3:
            err_rel[i] = 100.0 * err_abs[i] / abs(stress_theory[i])
        else:
            err_rel[i] = 0.0 if abs(err_abs[i]) < 1e-6 else 100.0

    print('单轴拉伸 30 步对比 (应变×1e6, 应力 MPa)')
    print('step  strain(1e-6)  theory(MPa)  forumat(MPa)  err_abs(Pa)   err_rel(%)  phase')
    print('-' * 80)
    for i in range(nstep):
        print('{:4d}  {:12.4f}  {:11.4f}  {:11.4f}  {:12.2f}  {:10.2f}  {:5d}'.format(
            steps[i],
            strain_f[i] * 1e6,
            stress_theory[i] / 1e6,
            stress_f[i] / 1e6,
            err_abs[i],
            err_rel[i],
            phase[i],
        ))
    print('-' * 80)
    print('最大绝对误差 (Pa):', np.max(np.abs(err_abs)))
    print('最大相对误差 (%):', np.max(np.abs(err_rel)))
    imax = np.argmax(np.abs(err_rel))
    print('最大相对误差步:', steps[imax], '  strain=', strain_f[imax]*1e6, '  rel%=', err_rel[imax])

    out_table = base / 'compare_tension_30steps.txt'
    with open(out_table, 'w', encoding='utf-8') as f:
        f.write('step  strain_xx       stress_theory   stress_forumat  err_abs(Pa)   err_rel_pct  phase\n')
        for i in range(nstep):
            f.write('{:4d}  {:16.8e}  {:16.8e}  {:16.8e}  {:14.4f}  {:10.2f}  {:5d}\n'.format(
                steps[i], strain_f[i], stress_theory[i], stress_f[i], err_abs[i], err_rel[i], phase[i]))
    print('已写入:', out_table)


if __name__ == '__main__':
    main()
