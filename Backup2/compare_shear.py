# -*- coding: utf-8 -*-
"""
纯剪切 100 步：对比 forumat 与理论。
理论：主应变 e1=strain_xy, e2=-strain_xy（45°主向），主应力 s1=uniax(e1), s2=uniax(e2)；
      sigma_xx = sigma_yy = (s1+s2)/2, sigma_xy = (s1-s2)/2.
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
    if eps <= 0:
        return compression_backbone(eps)
    return tension_backbone(eps)


def theory_shear(strain_xy):
    """纯剪切理论: e1=strain_xy, e2=-strain_xy; s1,s2; sigma_xx=sigma_yy=(s1+s2)/2, sigma_xy=(s1-s2)/2."""
    e1 = strain_xy
    e2 = -strain_xy
    s1 = uniaxial_sigma(e1)
    s2 = uniaxial_sigma(e2)
    sxx = (s1 + s2) / 2.0
    sxy = (s1 - s2) / 2.0
    return sxx, sxx, sxy


def main():
    base = Path(__file__).parent
    try:
        data = np.loadtxt(base / 'shear_100steps.txt', skiprows=1)
        steps = data[:, 0].astype(int)
        strain_xy = data[:, 1]
        stress_xx_f = data[:, 2]
        stress_yy_f = data[:, 3]
        stress_xy_f = data[:, 4]
    except Exception as e:
        print('请先运行: test_shear.exe 生成 shear_100steps.txt')
        print(e)
        return

    theory_xx = np.zeros_like(strain_xy)
    theory_yy = np.zeros_like(strain_xy)
    theory_xy = np.zeros_like(strain_xy)
    for i, exy in enumerate(strain_xy):
        theory_xx[i], theory_yy[i], theory_xy[i] = theory_shear(exy)

    err_xx = stress_xx_f - theory_xx
    err_yy = stress_yy_f - theory_yy
    err_xy = stress_xy_f - theory_xy
    rel_xx = np.where(np.abs(theory_xx) > 1e-20, 100.0 * err_xx / theory_xx, 0.0)
    rel_yy = np.where(np.abs(theory_yy) > 1e-20, 100.0 * err_yy / theory_yy, 0.0)
    rel_xy = np.where(np.abs(theory_xy) > 1e-20, 100.0 * err_xy / theory_xy, 0.0)

    out = base / 'compare_shear.txt'
    with open(out, 'w', encoding='utf-8') as f:
        f.write('step  strain_xy     theory_xx  theory_yy  theory_xy  forumat_xx forumat_yy forumat_xy  err_xx(Pa) err_yy(Pa) err_xy(Pa)  rel_xx% rel_yy% rel_xy%\n')
        f.write('-' * 120 + '\n')
        for i in range(len(steps)):
            f.write(f'{steps[i]:4d}  {strain_xy[i]:.4e}  {theory_xx[i]:10.2f} {theory_yy[i]:10.2f} {theory_xy[i]:10.2f}  '
                    f'{stress_xx_f[i]:10.2f} {stress_yy_f[i]:10.2f} {stress_xy_f[i]:10.2f}  '
                    f'{err_xx[i]:10.2f} {err_yy[i]:10.2f} {err_xy[i]:10.2f}  '
                    f'{rel_xx[i]:6.2f} {rel_yy[i]:6.2f} {rel_xy[i]:6.2f}\n')

    max_err_xx = np.max(np.abs(err_xx))
    max_err_yy = np.max(np.abs(err_yy))
    max_err_xy = np.max(np.abs(err_xy))
    max_rel_xx = np.max(np.abs(rel_xx))
    max_rel_yy = np.max(np.abs(rel_yy))
    max_rel_xy = np.max(np.abs(rel_xy))
    print('剪切 100 步对比 (应力 MPa)')
    print('最大绝对误差 (Pa): xx=', max_err_xx, ' yy=', max_err_yy, ' xy=', max_err_xy)
    print('最大相对误差 (%):  xx=', max_rel_xx, ' yy=', max_rel_yy, ' xy=', max_rel_xy)
    print('已写入:', out)
    # 打印前 15 步与后 5 步
    print('\n前 15 步:')
    for i in range(min(15, len(steps))):
        print(f'  step {steps[i]:3d}  exy={strain_xy[i]:.6f}  theory sxx={theory_xx[i]:8.2f} syy={theory_yy[i]:8.2f} sxy={theory_xy[i]:8.2f}  '
              f'forumat sxx={stress_xx_f[i]:8.2f} syy={stress_yy_f[i]:8.2f} sxy={stress_xy_f[i]:8.2f}  rel_xy%={rel_xy[i]:.2f}')
    print('...')
    print('后 5 步:')
    for i in range(max(0, len(steps)-5), len(steps)):
        print(f'  step {steps[i]:3d}  exy={strain_xy[i]:.6f}  theory sxx={theory_xx[i]:8.2f} syy={theory_yy[i]:8.2f} sxy={theory_xy[i]:8.2f}  '
              f'forumat sxx={stress_xx_f[i]:8.2f} syy={stress_yy_f[i]:8.2f} sxy={stress_xy_f[i]:8.2f}  rel_xy%={rel_xy[i]:.2f}')


if __name__ == '__main__':
    main()
