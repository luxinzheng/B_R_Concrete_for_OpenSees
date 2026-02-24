# -*- coding: utf-8 -*-
"""
双轴单调 5 路径：逐步对比 forumat 与理论。
理论：主向独立单轴响应。当 strain_xy=0 时坐标轴即主向，sigma_xx = uniax(eps_xx), sigma_yy = uniax(eps_yy)。
"""
import numpy as np
from pathlib import Path

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
    if eps <= 0:
        return 0.0
    if eps <= eps_peak:
        return E * eps
    if eps <= 10.0 * eps_peak:
        return ft * (1.0 - (eps - eps_peak) / (9.0 * eps_peak))
    return 0.0


def uniaxial_sigma(eps):
    """单轴应力：压用 Saenz，拉用弹性+软化。"""
    if eps <= 0:
        return compression_backbone(eps)
    return tension_backbone(eps)


def main():
    base = Path(__file__).parent
    try:
        data = np.loadtxt(base / 'biax_5paths.txt', skiprows=1)
        path_id = data[:, 0].astype(int)
        steps = data[:, 1].astype(int)
        strain_xx = data[:, 2]
        strain_yy = data[:, 3]
        stress_xx_f = data[:, 4]
        stress_yy_f = data[:, 5]
    except Exception as e:
        print('请先运行: test_biax_5paths.exe 生成 biax_5paths.txt')
        print(e)
        return

    theory_xx = np.array([uniaxial_sigma(e) for e in strain_xx])
    theory_yy = np.array([uniaxial_sigma(e) for e in strain_yy])
    err_xx = stress_xx_f - theory_xx
    err_yy = stress_yy_f - theory_yy
    err_rel_xx = np.zeros(len(steps))
    err_rel_yy = np.zeros(len(steps))
    for i in range(len(steps)):
        if abs(theory_xx[i]) > 1e3:
            err_rel_xx[i] = 100.0 * err_xx[i] / abs(theory_xx[i])
        if abs(theory_yy[i]) > 1e3:
            err_rel_yy[i] = 100.0 * err_yy[i] / abs(theory_yy[i])

    path_names = ['-1:-1', '1:1', '-1:-0.5', '-1:0.1', '-1:0.2']
    print('双轴 5 路径逐步对比 (应力 MPa)')
    print('path  step  strain_xx   strain_yy   theory_xx  theory_yy  forumat_xx forumat_yy  err_xx(Pa)  err_yy(Pa)  rel_xx%  rel_yy%')
    print('-' * 115)
    for i in range(len(steps)):
        p = path_id[i]
        print('{:2d}  {:4d}  {:9.4e}  {:9.4e}  {:9.4f}  {:9.4f}  {:9.4f}  {:9.4f}  {:10.2f}  {:10.2f}  {:6.2f}  {:6.2f}'.format(
            p, steps[i], strain_xx[i], strain_yy[i],
            theory_xx[i]/1e6, theory_yy[i]/1e6, stress_xx_f[i]/1e6, stress_yy_f[i]/1e6,
            err_xx[i], err_yy[i], err_rel_xx[i], err_rel_yy[i]))
    print('-' * 115)
    print('最大绝对误差 stress_xx (Pa):', np.max(np.abs(err_xx)))
    print('最大绝对误差 stress_yy (Pa):', np.max(np.abs(err_yy)))
    print('最大相对误差 xx (%):', np.max(np.abs(err_rel_xx)))
    print('最大相对误差 yy (%):', np.max(np.abs(err_rel_yy)))

    out_table = base / 'compare_biax_5paths.txt'
    with open(out_table, 'w', encoding='utf-8') as f:
        f.write('path  step  strain_xx       strain_yy       theory_xx     theory_yy     forumat_xx    forumat_yy    err_xx(Pa)   err_yy(Pa)   rel_xx%  rel_yy%\n')
        for i in range(len(steps)):
            f.write('{:2d}  {:4d}  {:16.8e}  {:16.8e}  {:13.4e}  {:13.4e}  {:13.4e}  {:13.4e}  {:12.4f}  {:12.4f}  {:7.2f}  {:7.2f}\n'.format(
                path_id[i], steps[i], strain_xx[i], strain_yy[i],
                theory_xx[i], theory_yy[i], stress_xx_f[i], stress_yy_f[i],
                err_xx[i], err_yy[i], err_rel_xx[i], err_rel_yy[i]))
    print('已写入:', out_table)


if __name__ == '__main__':
    main()
