# -*- coding: utf-8 -*-
"""计算往复加载各阶段目标应变。与 verify_uniax_adina 一致的 Saenz 参数。"""
import numpy as np

E     = 30.0e9
fc    = 30.0e6
epsc  = 0.002
ft    = 3.0e6
eps_peak = ft / E
EPSCP = epsc
SIGCP = fc
RP    = 0.004 / epsc  # epsu/epsc
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
    Et = num / den**2
    return max(Et, -abs(E))

def saenz_stress_single(eps):
    """单点压应变 eps (负) 对应的 Saenz 应力 (负)。"""
    if eps >= 0:
        return 0.0
    n = 200
    ep = np.linspace(0, eps, n)
    DE = -ep / EPSCP
    Et = np.array([saenz_tangent(d) for d in DE])
    sigma = np.zeros(n)
    for i in range(1, n):
        sigma[i] = sigma[i-1] + 0.5*(Et[i]+Et[i-1])*(ep[i]-ep[i-1])
    return sigma[-1]

# 求应变使应力 = 0.5*fc (压应力 -0.5*fc)
target_sigma = -0.5 * fc
eps_lo, eps_hi = -epsc * 1.2, 0.0
for _ in range(60):
    mid = 0.5 * (eps_lo + eps_hi)
    s = saenz_stress_single(mid)
    if s > target_sigma:
        eps_hi = mid
    else:
        eps_lo = mid
eps_c50 = 0.5 * (eps_lo + eps_hi)
print("Phase1 target strain (50% fc):", eps_c50)
print("Check stress at eps_c50:", saenz_stress_single(eps_c50) / 1e6, "MPa (expect -15)")
