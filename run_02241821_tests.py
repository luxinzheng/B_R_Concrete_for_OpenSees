# -*- coding: utf-8 -*-
"""Run sw1-1, sw2-1, Multi-layer_Shell with OpenSees-02241821.exe and plot comparison."""
import subprocess, os, time, sys
import numpy as np

root = r'e:\Basic\Concrete_Model\ADINA'
exe = os.path.join(root, 'OpenSees-02241821.exe')

cases = [
    {'name': 'sw1-1',  'tcl': 'br_test.tcl',
     'disp': '53.txt', 'react': '1.txt',
     'ref_disp': 'ref_53.txt', 'ref_react': 'ref_1.txt',
     'disp_cols': (0,1), 'react_sum_cols': list(range(1,6)),
     'label': 'SW1-1 (Cyclic Shear Wall)'},
    {'name': 'sw2-1',  'tcl': 'test_br.tcl',
     'disp': 'test_disp.txt', 'react': 'test_react.txt',
     'ref_disp': 'ref_28.txt', 'ref_react': 'ref_1.txt',
     'disp_cols': (0,1), 'react_sum_cols': list(range(1,6)),
     'label': 'SW2-1 (Cyclic Shear Wall)'},
    {'name': 'Multi-layer_Shell', 'tcl': 'br_test.tcl',
     'disp': 'disp1.txt', 'react': 'shearforce1.txt',
     'ref_disp': 'ref_disp1.txt', 'ref_react': 'ref_shearforce1.txt',
     'disp_cols': (0,1), 'react_sum_cols': list(range(1,10)),
     'label': 'Multi-layer Shell (Cyclic)'},
]

results = {}
for c in cases:
    d = os.path.join(root, c['name'])
    for fn in [c['disp'], c['react']]:
        fp = os.path.join(d, fn)
        if os.path.isfile(fp):
            os.remove(fp)

    print(f"\n{'='*60}", flush=True)
    print(f"Running {c['name']} with {os.path.basename(exe)}...", flush=True)
    t0 = time.time()
    log_fp = os.path.join(d, f"{c['name']}_02241821_log.txt")
    try:
        with open(log_fp, 'w') as logf:
            p = subprocess.run([exe, c['tcl']], cwd=d,
                               stdout=logf, stderr=subprocess.STDOUT, timeout=600)
        elapsed = time.time() - t0
        print(f"  Finished in {elapsed:.1f}s, exit={p.returncode}", flush=True)
    except subprocess.TimeoutExpired:
        elapsed = time.time() - t0
        print(f"  TIMEOUT after {elapsed:.1f}s", flush=True)

    if os.path.isfile(log_fp):
        with open(log_fp, 'r', errors='replace') as f:
            lines = f.readlines()
        nan_c = sum(1 for l in lines if '-1.#IND' in l or 'nan' in l.lower())
        warn_c = sum(1 for l in lines if 'failed to converge' in l.lower())
        print(f"  Log: {len(lines)} lines, NaN={nan_c}, warnings={warn_c}")
        for l in lines[-5:]:
            print(f"    {l.rstrip()}")

    info = {'elapsed': elapsed, 'ok': True}
    for key in ['disp', 'react', 'ref_disp', 'ref_react']:
        fp = os.path.join(d, c[key])
        if os.path.isfile(fp) and os.path.getsize(fp) > 0:
            try:
                info[key] = np.loadtxt(fp)
                print(f"  {c[key]}: {info[key].shape}")
            except Exception as e:
                print(f"  {c[key]}: load error {e}")
                info['ok'] = False
        else:
            print(f"  {c[key]}: MISSING or empty")
            info['ok'] = False
    results[c['name']] = info

# ═══════════════════════════════════════════════════════════
# Plot comparison
# ═══════════════════════════════════════════════════════════
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10})

fig, axes = plt.subplots(len(cases), 2, figsize=(16, 5*len(cases)))
fig.suptitle(f'OpenSees-02241821.exe  —  All TCL Test Cases vs Reference',
             fontsize=14, fontweight='bold', y=0.99)

for row, c in enumerate(cases):
    info = results[c['name']]
    d = os.path.join(root, c['name'])

    # Left: Force-Displacement
    ax = axes[row, 0]
    if info['ok'] and 'disp' in info and 'react' in info:
        disp_data = info['disp']
        react_data = info['react']
        t_d = disp_data[:, 0]
        u_d = disp_data[:, 1] if disp_data.ndim > 1 else disp_data
        t_r = react_data[:, 0]
        F_r = np.sum(react_data[:, c['react_sum_cols']], axis=1)
        n = min(len(t_d), len(t_r))
        ax.plot(u_d[:n]*1000, F_r[:n]/1000, 'b-', lw=0.8, label='ADINA model (02241821)')

    if 'ref_disp' in info and 'ref_react' in info:
        rd = info['ref_disp']
        rr = info['ref_react']
        t_rd = rd[:, 0]
        u_rd = rd[:, 1] if rd.ndim > 1 else rd
        t_rr = rr[:, 0]
        F_rr = np.sum(rr[:, c['react_sum_cols']], axis=1) if rr.ndim > 1 else rr
        nr = min(len(t_rd), len(t_rr))
        ax.plot(u_rd[:nr]*1000, F_rr[:nr]/1000, 'r--', lw=1.2, label='Reference model')

    ax.set_xlabel('Displacement (mm)')
    ax.set_ylabel('Base Shear (kN)')
    ax.set_title(f'{c["label"]} — Force vs Displacement', fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='k', lw=0.5)
    ax.axvline(0, color='k', lw=0.5)

    # Right: Displacement history
    ax = axes[row, 1]
    if info['ok'] and 'disp' in info:
        disp_data = info['disp']
        ax.plot(disp_data[:, 0], disp_data[:, 1]*1000, 'b-', lw=0.8, label='ADINA (02241821)')
    if 'ref_disp' in info:
        rd = info['ref_disp']
        ax.plot(rd[:, 0], rd[:, 1]*1000, 'r--', lw=1.0, label='Reference')
    ax.set_xlabel('Pseudo-time')
    ax.set_ylabel('Displacement (mm)')
    ax.set_title(f'{c["label"]} — Displacement History', fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

plt.tight_layout(rect=[0, 0, 1, 0.97])
out_path = os.path.join(root, 'opensees_02241821_comparison.png')
plt.savefig(out_path, dpi=150, bbox_inches='tight')
print(f"\nPlot saved: {out_path}")
plt.close()
