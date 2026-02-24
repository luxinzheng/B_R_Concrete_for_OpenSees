# -*- coding: utf-8 -*-
import numpy as np, sys

def report(fname):
    try:
        d = np.loadtxt(fname, skiprows=1)
    except Exception as e:
        print(f"  {fname}: LOAD ERROR {e}")
        return
    step = d[:,0]
    eps_x, eps_y, gamma = d[:,1], d[:,2], d[:,3]
    sig_x, sig_y, tau   = d[:,4], d[:,5], d[:,6]
    angle, pgrav         = d[:,7], d[:,8]
    print(f"  {fname}: {len(step)} steps")
    print(f"    sig_x: min={sig_x.min():.4f}  max={sig_x.max():.4f}")
    print(f"    sig_y: min={sig_y.min():.4f}  max={sig_y.max():.4f}")
    print(f"    tau:   min={tau.min():.4f}  max={tau.max():.4f}")
    print(f"    ANGLE: unique={np.unique(angle)}")
    print(f"    PGRAV: unique={np.unique(pgrav)}")
    if np.any(np.isnan(d)) or np.any(np.isinf(d)):
        print(f"    *** NaN/Inf DETECTED ***")
    if np.any(np.abs(d[:,4:7]) > 1e3):
        print(f"    *** STRESS > 1e3 MPa DETECTED ***")

for f in ['out_uniax_comp.csv', 'out_uniax_tens.csv', 'out_pure_shear.csv']:
    report(f)
