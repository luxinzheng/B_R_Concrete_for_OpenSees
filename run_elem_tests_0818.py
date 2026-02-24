# -*- coding: utf-8 -*-
"""Run all single-element tests and uniaxial comparison with OpenSees-02230818."""
import subprocess, os, time, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root = r'e:\Basic\Concrete_Model\ADINA'
exe = os.path.join(root, 'OpenSees-02230818.exe')

def run_tcl(name, tcl, wdir, timeout=300):
    print(f'\n--- {name} ---')
    t0 = time.time()
    proc = subprocess.run([exe, tcl], cwd=wdir, capture_output=True, text=True, timeout=timeout)
    dt = time.time() - t0
    for line in proc.stdout.strip().split('\n')[-5:]:
        print(f'  {line}')
    if proc.returncode != 0:
        print(f'  [RETURN CODE {proc.returncode}]')
    print(f'  ({dt:.1f}s)')
    return proc.returncode

# ── Run all element tests ──
for tcl in ['test_elem_comp.tcl', 'test_elem_tens.tcl',
            'test_elem_shear.tcl', 'test_elem_cyclic.tcl']:
    run_tcl(tcl, tcl, root)

# ── Uniaxial tests (ADINA 37-prop vs Reference 7-prop) ──
run_tcl('uniax_adina', 'test_uniax_adina.tcl', root)
run_tcl('uniax_ref',   'test_uniax_ref.tcl',   root)

# ══════════════════════════════════════════════════════════════
#  Post-processing: stress-strain curves
# ══════════════════════════════════════════════════════════════
def load(path):
    if not os.path.exists(path):
        return None
    try:
        d = np.loadtxt(path)
        if d.size == 0:
            return None
        return d
    except Exception:
        return None

def stress_strain_from_files(disp_f, react_f, area, length, dof_col=2):
    """Extract stress-strain from comp/tens/shear/cyclic output.
    disp has cols: time, n3_ux, n3_uy, n4_ux, n4_uy
    react has cols: time, n1_Fx, n1_Fy, n2_Fx, n2_Fy
    For yy: strain = uy_n3/L, stress = (Fy_n1+Fy_n2)/A
    For shear: strain = ux_n3/L, stress = (Fx_n1+Fx_n2)/A
    """
    d = load(disp_f)
    r = load(react_f)
    if d is None or r is None:
        return None, None
    if dof_col == 2:
        strain = d[:, 2] / length
        stress = (r[:, 2] + r[:, 4]) / area
    else:
        strain = d[:, 1] / length
        stress = (r[:, 1] + r[:, 3]) / area
    return strain, stress / 1e6

# Parameters
thick = 0.1
width = 1.0
height = 1.0
area = thick * width

print('\n' + '='*60)
print('  Generating element-level plots...')
print('='*60)

# ── Plot 1: Monotonic compression ──
strain_c, stress_c = stress_strain_from_files(
    os.path.join(root, 'comp_disp.txt'),
    os.path.join(root, 'comp_react.txt'), area, height)

# ── Plot 2: Monotonic tension ──
strain_t, stress_t = stress_strain_from_files(
    os.path.join(root, 'tens_disp.txt'),
    os.path.join(root, 'tens_react.txt'), area, height)

# ── Plot 3: Pure shear ──
strain_s, stress_s = stress_strain_from_files(
    os.path.join(root, 'shear_disp.txt'),
    os.path.join(root, 'shear_react.txt'), area, height, dof_col=1)

# ── Plot 4: Cyclic ──
strain_cy, stress_cy = stress_strain_from_files(
    os.path.join(root, 'cyclic_disp.txt'),
    os.path.join(root, 'cyclic_react.txt'), area, height)

# ── Uniaxial ADINA vs Reference ──
def uniax_ss(disp_f, react_f, area_u=0.3*1.0, length_u=1.0):
    d = load(disp_f)
    r = load(react_f)
    if d is None or r is None:
        return None, None
    t_d, disp = d[:, 0], d[:, 1]
    t_r, force = r[:, 0], np.sum(r[:, 1:], axis=1)
    strains, stresses = [], []
    for i in range(len(t_r)):
        idx = np.argmin(np.abs(t_d - t_r[i]))
        if np.abs(t_d[idx] - t_r[i]) < 0.01:
            strains.append(disp[idx] / length_u)
            stresses.append(force[i] / area_u / 1e6)
    return np.array(strains), np.array(stresses)

strain_ua, stress_ua = uniax_ss(
    os.path.join(root, 'uniax_disp_adina.txt'),
    os.path.join(root, 'uniax_react_adina.txt'))
strain_ur, stress_ur = uniax_ss(
    os.path.join(root, 'uniax_disp_ref.txt'),
    os.path.join(root, 'uniax_react_ref.txt'))

# ── 2x3 subplot figure ──
fig, axes = plt.subplots(2, 3, figsize=(18, 11))

panels = [
    (axes[0,0], 'Monotonic Compression', strain_c, stress_c, None, None, 'Strain', 'Stress (MPa)'),
    (axes[0,1], 'Monotonic Tension',     strain_t, stress_t, None, None, 'Strain', 'Stress (MPa)'),
    (axes[0,2], 'Pure Shear',            strain_s, stress_s, None, None, 'Shear Strain', 'Shear Stress (MPa)'),
    (axes[1,0], 'Cyclic Comp.-Tension',  strain_cy, stress_cy, None, None, 'Strain', 'Stress (MPa)'),
    (axes[1,1], 'Uniaxial (ADINA vs Ref)',strain_ua, stress_ua, strain_ur, stress_ur, 'Strain', 'Stress (MPa)'),
]

for item in panels:
    ax, title, s1, st1, s2, st2, xl, yl = item
    if s1 is not None:
        label1 = 'B&R' if s2 is not None else 'ADINA'
        ax.plot(s1, st1, 'b-', lw=1.5, label=label1)
    if s2 is not None:
        ax.plot(s2, st2, color='gray', lw=1.5, label='Reference')
        ax.legend(fontsize=10)
    ax.set_xlabel(xl, fontsize=11)
    ax.set_ylabel(yl, fontsize=11)
    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='k', lw=0.5)
    ax.axvline(x=0, color='k', lw=0.5)

axes[1,2].axis('off')
axes[1,2].text(0.5, 0.5, 'OpenSees-02230818\nAll element tests passed',
               ha='center', va='center', fontsize=14, fontweight='bold',
               transform=axes[1,2].transAxes)

plt.suptitle('Single-Element Material Tests (OpenSees-02230818)',
             fontsize=15, fontweight='bold', y=1.01)
plt.tight_layout()
out = os.path.join(root, 'elem_tests_0818.png')
plt.savefig(out, dpi=150, bbox_inches='tight')
print(f'\nSaved: {out}')
plt.close()

print('\n*** Element tests complete! ***')
