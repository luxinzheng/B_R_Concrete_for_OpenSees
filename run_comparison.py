#!/usr/bin/env python3
"""
run_comparison.py - Compare B&R (OpenSees-02210959) vs Reference model
for three test cases: Multi-layer_Shell, sw1-1, sw2-1
"""
import os, sys, subprocess, re, shutil, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

ROOT = os.path.dirname(os.path.abspath(__file__))
BR_EXE = os.path.join(ROOT, 'OpenSees-02210959.exe')

def find_ref_exe():
    for entry in os.scandir(ROOT):
        if entry.is_file() and entry.name.endswith('.exe'):
            if '\u53c2\u8003' in entry.name:
                return entry.path
    return None

def aci_E0(fc_abs_pa):
    return 4730.0 * ((fc_abs_pa / 1e6) ** 0.5) * 1e6

EXTRA_30 = [
    1.0e-4, 0.5, 0.75, 1.0, 0.7, 0.12,
    0.0, 0.25, 0.5, 0.75, 1.0, 1.2,
    1.0, 1.4, 1.7, 2.2, 2.5, 2.8,
    1.3, 1.5, 2.0, 2.3, 2.7, 3.2,
    1.25, 1.45, 1.95, 2.25, 2.65, 3.15,
]

def br_mat_line(tag, fc_abs, ft, fcu, epsc, epscu):
    E0 = aci_E0(fc_abs)
    nu = 0.2
    params = [E0, nu, ft, -fc_abs, epsc, fcu, epscu] + EXTRA_30
    ps = '  '.join(f'{p:.6g}' for p in params)
    return f'nDMaterial PlaneStressUserMaterial    {tag}   40   37   {ps}'

# ── Robust analysis TCL snippet (same strategy as verify_wall.tcl) ──
ROBUST_LOOP_SW = r'''
# Robust lateral analysis with fallback strategies
set totalTime 50.0
set dt 0.1
set step 0; set nFail 0; set nSkip 0; set consecFail 0

while {[getTime] < $totalTime} {
    set ok [analyze 1]; incr step
    if {$ok == 0} {
        set consecFail 0
        if {[expr {$step % 50}] == 0} { record; puts "  step $step t=[getTime] f=$nFail s=$nSkip" }
        continue
    }
    algorithm ModifiedNewton; set ok [analyze 1]; algorithm KrylovNewton
    if {$ok == 0} { set consecFail 0; continue }

    set nSub 10
    integrator LoadControl [expr {$dt/$nSub}]
    test NormDispIncr 1.0e-4 500 0; algorithm Newton
    set ok 0
    for {set i 0} {$i < $nSub} {incr i} { set ok [analyze 1]; if {$ok != 0} break }
    integrator LoadControl $dt; algorithm KrylovNewton; test NormDispIncr 1.0e-5 1000 2
    if {$ok == 0} { set consecFail 0; continue }

    set nSub 100
    integrator LoadControl [expr {$dt/$nSub}]
    test NormDispIncr 1.0e-3 1000 0; algorithm Newton
    set ok 0
    for {set i 0} {$i < $nSub} {incr i} {
        set ok [analyze 1]
        if {$ok != 0} { algorithm ModifiedNewton; set ok [analyze 1]; algorithm Newton }
        if {$ok != 0} break
    }
    integrator LoadControl $dt; algorithm KrylovNewton; test NormDispIncr 1.0e-5 1000 2
    if {$ok == 0} { set consecFail 0; continue }

    set nSub 1000
    integrator LoadControl [expr {$dt/$nSub}]
    algorithm ModifiedNewton -initial; test NormDispIncr 1.0e10 1 0
    for {set i 0} {$i < $nSub} {incr i} { set ok [analyze 1]; if {$ok != 0} break }
    integrator LoadControl $dt; algorithm KrylovNewton; test NormDispIncr 1.0e-5 1000 2
    if {$ok == 0} { incr nSkip; set consecFail 0; continue }

    incr nFail; incr consecFail
    if {$consecFail >= 20} { puts "*** $consecFail consecutive fails -- stopping ***"; break }
}
record
puts "Lateral done: t=[getTime] steps=$step fails=$nFail skips=$nSkip"
'''

def robust_disp_control_snippet(ctrl_node, ctrl_dof, dincr, nsteps):
    """Generate robust displacement-controlled analysis block."""
    return f'''
# Robust displacement control: node={ctrl_node} dof={ctrl_dof} dincr={dincr} nsteps={nsteps}
set _dincr {dincr}
set _nstep {nsteps}
for {{set _s 0}} {{$_s < $_nstep}} {{incr _s}} {{
    integrator DisplacementControl {ctrl_node} {ctrl_dof} $_dincr
    set ok [analyze 1]
    if {{$ok == 0}} {{ continue }}
    algorithm ModifiedNewton
    integrator DisplacementControl {ctrl_node} {ctrl_dof} $_dincr
    set ok [analyze 1]; algorithm KrylovNewton
    if {{$ok == 0}} {{ continue }}
    set _nSub 10
    set _dsub [expr {{$_dincr / $_nSub}}]
    integrator DisplacementControl {ctrl_node} {ctrl_dof} $_dsub
    test NormDispIncr 1.0e-4 500 0; algorithm Newton
    set ok 0
    for {{set _i 0}} {{$_i < $_nSub}} {{incr _i}} {{ set ok [analyze 1]; if {{$ok != 0}} break }}
    algorithm KrylovNewton; test NormDispIncr 1.0e-4 200 2
    if {{$ok == 0}} {{ continue }}
    set _nSub 100
    set _dsub [expr {{$_dincr / $_nSub}}]
    integrator DisplacementControl {ctrl_node} {ctrl_dof} $_dsub
    test NormDispIncr 1.0e-3 1000 0; algorithm Newton
    set ok 0
    for {{set _i 0}} {{$_i < $_nSub}} {{incr _i}} {{
        set ok [analyze 1]
        if {{$ok != 0}} {{ algorithm ModifiedNewton; set ok [analyze 1]; algorithm Newton }}
        if {{$ok != 0}} break
    }}
    algorithm KrylovNewton; test NormDispIncr 1.0e-4 200 2
    if {{$ok == 0}} {{ continue }}
    puts "DC step failed: node={ctrl_node} dincr=$_dincr step=$_s"
}}
record
'''

def create_br_tcl_sw(src_tcl, dst_tcl, mat_defs):
    """Create B&R TCL for sw models with robust analysis."""
    with open(src_tcl, 'r', encoding='utf-8', errors='replace') as f:
        content = f.read()
    for tag, fc_abs, ft, fcu, epsc, epscu in mat_defs:
        pat = re.compile(
            r'nDMaterial\s+PlaneStressUserMaterial\s+' + str(tag)
            + r'\s+\d+\s+7\s+[^\n]+', re.IGNORECASE)
        content = pat.sub(br_mat_line(tag, fc_abs, ft, fcu, epsc, epscu), content)
    # Fix gravity analysis: BFGS→KrylovNewton, loosen tolerance
    content = re.sub(r'algorithm\s+BFGS\s+-count\s+\d+',
                     'algorithm KrylovNewton', content)
    content = re.sub(r'test\s+NormDispIncr\s+1\.0e-6\s+200\s*;?',
                     'test NormDispIncr 1.0e-5 500', content)
    # Replace the final analyze 500 with robust loop
    content = re.sub(r'analyze\s+500\s*$', ROBUST_LOOP_SW, content, flags=re.MULTILINE)
    with open(dst_tcl, 'w', encoding='utf-8') as f:
        f.write(content)

def create_br_tcl_multilayer(src_tcl, dst_tcl, mat_defs):
    """Create B&R TCL for Multi-layer_Shell with robust displacement control."""
    with open(src_tcl, 'r', encoding='utf-8', errors='replace') as f:
        lines = f.readlines()
    out = []
    for line in lines:
        stripped = line.strip()
        replaced = False
        for tag, fc_abs, ft, fcu, epsc, epscu in mat_defs:
            pat = (r'nDMaterial\s+PlaneStressUserMaterial\s+'
                   + str(tag) + r'\s+\d+\s+7\s+')
            if re.match(pat, stripped, re.IGNORECASE):
                out.append(br_mat_line(tag, fc_abs, ft, fcu, epsc, epscu) + '\n')
                replaced = True
                break
        if replaced:
            continue

        # Replace displacement control blocks
        m = re.match(
            r'\s*integrator\s+DisplacementControl\s+(\d+)\s+(\d+)\s+([^\s;]+)',
            stripped)
        if m:
            ctrl_node = m.group(1)
            ctrl_dof = m.group(2)
            dincr = m.group(3)
            # Look for the following analyze line
            out.append(f'# Original: {stripped}\n')
            out.append(f'set _dc_node {ctrl_node}\n')
            out.append(f'set _dc_dof {ctrl_dof}\n')
            out.append(f'set _dc_dincr {dincr}\n')
            continue

        m2 = re.match(r'\s*analyze\s+(\d+)\s*;?\s*$', stripped)
        if m2 and any('_dc_dincr' in l for l in out[-5:]):
            nsteps = m2.group(1)
            snippet = robust_disp_control_snippet(
                '$_dc_node', '$_dc_dof', '$_dc_dincr', nsteps)
            out.append(snippet + '\n')
            continue

        out.append(line)

    with open(dst_tcl, 'w', encoding='utf-8') as f:
        f.writelines(out)


def run_opensees(exe, tcl_file, cwd, timeout=1800):
    cmd = [exe, os.path.basename(tcl_file)]
    print(f'  Running: {os.path.basename(exe)} {os.path.basename(tcl_file)}')
    t0 = time.time()
    try:
        r = subprocess.run(cmd, cwd=cwd, capture_output=True, timeout=timeout)
        dt = time.time() - t0
        stdout = r.stdout.decode('utf-8', errors='replace') if r.stdout else ''
        stderr = r.stderr.decode('utf-8', errors='replace') if r.stderr else ''
        print(f'    exit={r.returncode}  time={dt:.1f}s')
        if r.returncode != 0:
            print(f'    STDERR: {stderr[:200]}')
        for line in (stdout + '\n' + stderr).splitlines():
            lup = line.upper().strip()
            if any(k in lup for k in ['FAIL', 'DONE', 'COMPLETE', 'GRAVITY',
                                       'LATERAL', 'WALL', 'STEP ', 'STOPPING']):
                if 'CTEST' not in lup and 'CONVERGENCETEST' not in lup:
                    print(f'    {line.strip()}')
        return r.returncode == 0, dt
    except subprocess.TimeoutExpired:
        print(f'    TIMEOUT after {time.time()-t0:.1f}s')
        return False, time.time() - t0
    except Exception as e:
        print(f'    ERROR: {e}')
        return False, time.time() - t0


def load_data(filepath):
    """Load whitespace-delimited data, tolerating rows with varying columns."""
    try:
        with open(filepath, 'r') as f:
            raw_lines = f.readlines()
        if not raw_lines:
            return None
        rows = []
        ncols = None
        for line in raw_lines:
            parts = line.split()
            if not parts:
                continue
            vals = []
            for p in parts:
                try:
                    vals.append(float(p))
                except ValueError:
                    break
            if not vals:
                continue
            if ncols is None:
                ncols = len(vals)
            if len(vals) == ncols:
                rows.append(vals)
        if not rows:
            return None
        return np.array(rows)
    except Exception:
        return None


def parse_hysteresis(disp_file, force_file):
    disp_data = load_data(disp_file)
    force_data = load_data(force_file)
    if disp_data is None or force_data is None:
        return None, None
    disp = disp_data[:, 1]
    n = min(len(disp), len(force_data))
    force = np.sum(force_data[:n, 1:], axis=1)
    disp = disp[:n]
    return disp, force


CASES = {
    'Multi-layer_Shell': {
        'dir': os.path.join(ROOT, 'Multi-layer_Shell'),
        'tcl': 'main.tcl',
        'mat_defs': [
            (1, 25.8e6, 2e6, -5.16e6, -0.003, -0.021),
            (11, 22.8e6, 2e6, -4.56e6, -0.002, -0.012),
        ],
        'disp_file': 'disp1.txt',
        'force_file': 'shearforce1.txt',
        'title': 'Multi-layer Shell',
        'type': 'multilayer',
    },
    'sw1-1': {
        'dir': os.path.join(ROOT, 'sw1-1'),
        'tcl': 'csp3.tcl',
        'mat_defs': [
            (1, 20.7e6, 2.07e6, -4.14e6, -0.002, -0.006),
        ],
        'disp_file': '53.txt',
        'force_file': '1.txt',
        'title': 'SW1-1',
        'type': 'sw',
    },
    'sw2-1': {
        'dir': os.path.join(ROOT, 'sw2-1'),
        'tcl': 'csp3.tcl',
        'mat_defs': [
            (1, 30.8e6, 3.08e6, -6.16e6, -0.002, -0.005),
        ],
        'disp_file': '28.txt',
        'force_file': '1.txt',
        'title': 'SW2-1',
        'type': 'sw',
    },
}


def plot_comparison(results):
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    for idx, (name, res) in enumerate(results.items()):
        ax = axes[idx]
        cfg = res['cfg']
        has_ref = res['ref_disp'] is not None
        has_br = res['br_disp'] is not None
        if has_ref:
            ax.plot(res['ref_disp'] * 1e3, res['ref_force'] / 1e3,
                    'b-', linewidth=0.7, alpha=0.8, label='Reference')
        if has_br:
            ax.plot(res['br_disp'] * 1e3, res['br_force'] / 1e3,
                    'r-', linewidth=0.7, alpha=0.8, label='B&R')
        ax.set_xlabel('Displacement (mm)')
        ax.set_ylabel('Base Shear (kN)')
        ax.set_title(cfg['title'])
        if has_ref or has_br:
            ax.legend(loc='best')
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0, color='k', linewidth=0.5)
        ax.axvline(x=0, color='k', linewidth=0.5)

        ref_pk = np.max(np.abs(res['ref_force'])) / 1e3 if has_ref else 0
        br_pk = np.max(np.abs(res['br_force'])) / 1e3 if has_br else 0
        txt = f'Ref: {res["ref_steps"]} steps, peak={ref_pk:.1f}kN\n'
        txt += f'B&R: {res["br_steps"]} steps, peak={br_pk:.1f}kN'
        if ref_pk > 0 and br_pk > 0:
            txt += f'\nRatio: {br_pk/ref_pk:.3f}'
        ax.text(0.02, 0.02, txt, transform=ax.transAxes, fontsize=8,
                verticalalignment='bottom',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    plt.tight_layout()
    out = os.path.join(ROOT, 'comparison_3cases.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    print(f'Saved: {out}')
    plt.close()


def plot_individual(name, res):
    cfg = res['cfg']
    has_ref = res['ref_disp'] is not None
    has_br = res['br_disp'] is not None
    if not has_ref and not has_br:
        return

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    ax = axes[0]
    if has_ref:
        ax.plot(res['ref_disp'] * 1e3, res['ref_force'] / 1e3,
                'b-', linewidth=0.7, alpha=0.8, label='Reference')
    if has_br:
        ax.plot(res['br_disp'] * 1e3, res['br_force'] / 1e3,
                'r-', linewidth=0.7, alpha=0.8, label='B&R')
    ax.set_xlabel('Displacement (mm)')
    ax.set_ylabel('Base Shear (kN)')
    ax.set_title(f'{cfg["title"]} - Force-Displacement Hysteresis')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='k', linewidth=0.5)
    ax.axvline(x=0, color='k', linewidth=0.5)

    ax = axes[1]
    if has_ref:
        ax.plot(np.arange(len(res['ref_disp'])), res['ref_disp'] * 1e3,
                'b-', linewidth=0.7, alpha=0.8, label='Reference')
    if has_br:
        ax.plot(np.arange(len(res['br_disp'])), res['br_disp'] * 1e3,
                'r-', linewidth=0.7, alpha=0.8, label='B&R')
    ax.set_xlabel('Step')
    ax.set_ylabel('Displacement (mm)')
    ax.set_title(f'{cfg["title"]} - Displacement History')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    out = os.path.join(ROOT, f'comparison_{name}.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    print(f'Saved: {out}')
    plt.close()


def main():
    ref_exe = find_ref_exe()
    if ref_exe is None:
        print('ERROR: Cannot find reference exe')
        sys.exit(1)
    if not os.path.isfile(BR_EXE):
        print(f'ERROR: Cannot find B&R exe: {BR_EXE}')
        sys.exit(1)

    print(f'B&R exe: {os.path.basename(BR_EXE)}')
    print(f'Ref exe: {os.path.basename(ref_exe)}')
    print('=' * 70)

    # Check for existing ref data (skip re-running if recent)
    skip_ref = '--br-only' in sys.argv

    results = {}

    for name, cfg in CASES.items():
        print(f'\n{"="*70}')
        print(f'Case: {name} ({cfg["title"]})')
        print(f'{"="*70}')

        case_dir = cfg['dir']
        tcl_file = os.path.join(case_dir, cfg['tcl'])
        br_tcl = os.path.join(case_dir, f'br_{cfg["tcl"]}')
        out_files = [cfg['disp_file'], cfg['force_file']]

        ref_files_exist = all(
            os.path.isfile(os.path.join(case_dir, f'ref_{fn}'))
            and os.path.getsize(os.path.join(case_dir, f'ref_{fn}')) > 0
            for fn in out_files)

        # --- Reference Model ---
        if skip_ref and ref_files_exist:
            print('\n[Reference Model] Using cached data')
            ref_ok, ref_time = True, 0.0
        else:
            print('\n[Reference Model]')
            ref_ok, ref_time = run_opensees(ref_exe, tcl_file, case_dir)
            for fn in out_files:
                src = os.path.join(case_dir, fn)
                dst = os.path.join(case_dir, f'ref_{fn}')
                if os.path.isfile(src):
                    shutil.copy2(src, dst)

        # --- B&R Model ---
        print('\n[B&R Model]')
        if cfg['type'] == 'multilayer':
            create_br_tcl_multilayer(tcl_file, br_tcl, cfg['mat_defs'])
        else:
            create_br_tcl_sw(tcl_file, br_tcl, cfg['mat_defs'])
        print(f'  Created: {os.path.basename(br_tcl)}')
        for tag_info in cfg['mat_defs']:
            tag, fc_abs = tag_info[0], tag_info[1]
            E0 = aci_E0(fc_abs)
            print(f'    mat {tag}: fc={fc_abs/1e6:.1f}MPa  E0={E0/1e9:.2f}GPa')

        br_ok, br_time = run_opensees(BR_EXE, br_tcl, case_dir)
        for fn in out_files:
            src = os.path.join(case_dir, fn)
            dst = os.path.join(case_dir, f'br_{fn}')
            if os.path.isfile(src):
                shutil.copy2(src, dst)

        # --- Parse ---
        ref_disp, ref_force = parse_hysteresis(
            os.path.join(case_dir, f'ref_{cfg["disp_file"]}'),
            os.path.join(case_dir, f'ref_{cfg["force_file"]}'))
        br_disp, br_force = parse_hysteresis(
            os.path.join(case_dir, f'br_{cfg["disp_file"]}'),
            os.path.join(case_dir, f'br_{cfg["force_file"]}'))

        r = {
            'ref_ok': ref_ok, 'ref_time': ref_time,
            'br_ok': br_ok, 'br_time': br_time,
            'ref_disp': ref_disp, 'ref_force': ref_force,
            'br_disp': br_disp, 'br_force': br_force,
            'ref_steps': len(ref_disp) if ref_disp is not None else 0,
            'br_steps': len(br_disp) if br_disp is not None else 0,
            'cfg': cfg,
        }
        results[name] = r

        ref_pk = np.max(np.abs(ref_force)) / 1e3 if ref_force is not None else 0
        br_pk = np.max(np.abs(br_force)) / 1e3 if br_force is not None else 0
        ref_dmax = np.max(np.abs(ref_disp)) * 1e3 if ref_disp is not None else 0
        br_dmax = np.max(np.abs(br_disp)) * 1e3 if br_disp is not None else 0
        print(f'\n  Results:')
        print(f'    Reference: steps={r["ref_steps"]}, peak={ref_pk:.1f}kN, '
              f'max_disp={ref_dmax:.2f}mm, time={ref_time:.1f}s')
        print(f'    B&R:       steps={r["br_steps"]}, peak={br_pk:.1f}kN, '
              f'max_disp={br_dmax:.2f}mm, time={br_time:.1f}s')
        if ref_pk > 0 and br_pk > 0:
            print(f'    Peak ratio (B&R/Ref): {br_pk/ref_pk:.3f}')

    # --- Plots ---
    print('\n' + '=' * 70)
    print('Generating comparison plots...')
    plot_comparison(results)
    for name, res in results.items():
        plot_individual(name, res)

    # --- Summary ---
    print('\n' + '=' * 70)
    print('SUMMARY')
    print('=' * 70)
    hdr = f'{"Case":<25} {"Ref Steps":>10} {"B&R Steps":>10} '
    hdr += f'{"Ref Peak(kN)":>14} {"B&R Peak(kN)":>14} {"Ratio":>8}'
    print(hdr)
    print('-' * 85)
    for name, res in results.items():
        ref_pk = (np.max(np.abs(res['ref_force'])) / 1e3
                  if res['ref_force'] is not None else 0)
        br_pk = (np.max(np.abs(res['br_force'])) / 1e3
                 if res['br_force'] is not None else 0)
        ratio = br_pk / ref_pk if ref_pk > 0 else 0
        print(f'{name:<25} {res["ref_steps"]:>10d} {res["br_steps"]:>10d} '
              f'{ref_pk:>14.1f} {br_pk:>14.1f} {ratio:>8.3f}')
    print('\nDone.')


if __name__ == '__main__':
    main()
