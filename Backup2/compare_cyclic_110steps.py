# -*- coding: utf-8 -*-
"""
单轴拉压往复 110 步：6 阶段逐步对比 forumat 与理论。
理论：压骨架 Saenz；拉骨架 弹性至 ft 后线性软化；卸载/再加载为弹性 E。
"""
import numpy as np
from pathlib import Path

# 与 test_uniax_cyclic / forumat 一致
E     = 30.0e9
fc    = 30.0e6
ft    = 3.0e6
epsc  = 0.002
eps_peak = ft / E
EPSCP = epsc
SIGCP = fc
RP    = 0.004 / epsc
ES    = SIGCP / EPSCP
EU    = (0.85 * SIGCP) / (0.004 * EPSCP / epsc)
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
    return max(num / den**2, -abs(E))


def compression_backbone(eps):
    """压骨架：eps<=0 为 Saenz 应力(负)，eps>0 为 0。"""
    if eps >= 0:
        return 0.0
    n = 200
    ep = np.linspace(0, eps, n)
    DE = -ep / EPSCP
    Et = np.array([saenz_tangent(d) for d in DE])
    s = 0.0
    for i in range(1, n):
        s = s + 0.5*(Et[i]+Et[i-1])*(ep[i]-ep[i-1])
    return s


def tension_backbone(eps):
    """拉骨架：0~eps_peak 弹性至 ft，>eps_peak 线性软化至 10*eps_peak 为 0；与 forumat 一致。"""
    if eps <= 0:
        return 0.0
    if eps <= eps_peak:
        # 弹性段：应力 = E*eps，且不超过抗拉强度 ft
        return min(E * eps, ft)
    if eps <= 10.0 * eps_peak:
        return ft * (1.0 - (eps - eps_peak) / (9.0 * eps_peak))
    return 0.0


def theory_cyclic_step_by_step(strain_arr):
    """
    往复理论：逐步计算。
    应变增加 -> 向拉；应变减少 -> 向压。
    卸载用弹性 E；触及骨架后沿骨架。
    """
    strain_arr = np.asarray(strain_arr)
    nstep = len(strain_arr)
    sigma = np.zeros(nstep)
    for i in range(nstep):
        eps = strain_arr[i]
        if i == 0:
            if eps <= 0:
                sigma[i] = compression_backbone(eps)
            else:
                sigma[i] = tension_backbone(eps)
            continue
        eps_prev = strain_arr[i-1]
        sig_prev = sigma[i-1]
        sigma_el = sig_prev + E * (eps - eps_prev)
        if eps > eps_prev:
            # 应变增加：向拉
            if sig_prev <= 0:
                # 从压卸载
                sigma[i] = sigma_el
                tb = tension_backbone(eps)
                if sigma_el > tb:
                    sigma[i] = tb
            else:
                sigma[i] = tension_backbone(eps)
        else:
            # 应变减少：向压
            if sig_prev >= 0:
                sigma[i] = sigma_el
                cb = compression_backbone(eps)
                if sigma_el < cb:
                    sigma[i] = cb
            else:
                sigma[i] = compression_backbone(eps)
    return sigma


def main():
    base = Path(__file__).parent
    try:
        data = np.loadtxt(base / 'uniax_cyclic_110steps.txt', skiprows=1)
        steps = data[:, 0].astype(int)
        strain_f = data[:, 1]
        stress_f = data[:, 2]
        phase = data[:, 3].astype(int)
    except Exception as e:
        print('请先运行: test_uniax_cyclic_110steps.exe 生成 uniax_cyclic_110steps.txt')
        print(e)
        return

    nstep = len(steps)
    stress_theory = theory_cyclic_step_by_step(strain_f)
    err_abs = stress_f - stress_theory
    err_rel = np.zeros(nstep)
    for i in range(nstep):
        if abs(stress_theory[i]) > 1e3:
            err_rel[i] = 100.0 * err_abs[i] / abs(stress_theory[i])
        else:
            err_rel[i] = 0.0 if abs(err_abs[i]) < 1.0 else 100.0

    print('单轴拉压往复 110 步对比 (应力 MPa)')
    print('step  strain(1e-6)  theory(MPa)  forumat(MPa)  err_abs(Pa)   err_rel(%)  phase')
    print('-' * 85)
    for i in range(nstep):
        print('{:4d}  {:12.4f}  {:11.4f}  {:11.4f}  {:12.2f}  {:10.2f}  {:5d}'.format(
            steps[i], strain_f[i]*1e6, stress_theory[i]/1e6, stress_f[i]/1e6,
            err_abs[i], err_rel[i], phase[i]))
    print('-' * 85)
    print('最大绝对误差 (Pa):', np.max(np.abs(err_abs)))
    print('最大相对误差 (%):', np.max(np.abs(err_rel)))
    imax = np.argmax(np.abs(err_rel))
    print('最大相对误差步:', steps[imax], '  strain=', strain_f[imax]*1e6, '  rel%=', err_rel[imax])

    out_table = base / 'compare_cyclic_110steps.txt'
    with open(out_table, 'w', encoding='utf-8') as f:
        f.write('step  strain_xx       stress_theory   stress_forumat  err_abs(Pa)   err_rel_pct  phase\n')
        for i in range(nstep):
            f.write('{:4d}  {:16.8e}  {:16.8e}  {:16.8e}  {:14.4f}  {:10.2f}  {:5d}\n'.format(
                steps[i], strain_f[i], stress_theory[i], stress_f[i], err_abs[i], err_rel[i], phase[i]))
    print('已写入:', out_table)


if __name__ == '__main__':
    main()
