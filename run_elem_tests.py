"""
Run single-element TCL tests with all available OpenSees versions.
Compare results to identify numerical stability issues.
"""
import subprocess
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['font.family'] = ['Microsoft YaHei', 'SimHei', 'sans-serif']
import glob
import re
import time

WORK_DIR = r'e:\Basic\Concrete_Model\ADINA'
TESTS = ['test_elem_comp', 'test_elem_tens', 'test_elem_shear', 'test_elem_cyclic']
TEST_LABELS = {
    'test_elem_comp': '单轴压缩',
    'test_elem_tens': '单轴拉伸',
    'test_elem_shear': '纯剪切',
    'test_elem_cyclic': '循环拉压',
}

# Output file prefixes for each test
OUTPUT_PREFIX = {
    'test_elem_comp': 'comp',
    'test_elem_tens': 'tens',
    'test_elem_shear': 'shear',
    'test_elem_cyclic': 'cyclic',
}

# Element geometry
WIDTH = 1.0   # m
THICKNESS = 0.1  # m (2 layers x 0.05m)
AREA = WIDTH * THICKNESS

def find_opensees_versions():
    """Find all OpenSees executables."""
    exes = glob.glob(os.path.join(WORK_DIR, 'OpenSees*.exe'))
    versions = []
    for exe in sorted(exes):
        name = os.path.basename(exe)
        versions.append((name, exe))
    return versions

def run_test(exe_path, tcl_file, version_name, test_name):
    """Run a single test and return results."""
    prefix = OUTPUT_PREFIX[test_name]

    # Clean old output files
    for suffix in ['_disp.txt', '_react.txt']:
        f = os.path.join(WORK_DIR, prefix + suffix)
        if os.path.exists(f):
            os.remove(f)

    # Run OpenSees
    tcl_path = os.path.join(WORK_DIR, tcl_file + '.tcl')
    t0 = time.time()
    try:
        result = subprocess.run(
            [exe_path, tcl_path],
            capture_output=True, text=True, timeout=60,
            cwd=WORK_DIR
        )
        elapsed = time.time() - t0
        exit_code = result.returncode
        stdout = result.stdout
        stderr = result.stderr
    except subprocess.TimeoutExpired:
        elapsed = 60.0
        exit_code = -999
        stdout = "TIMEOUT"
        stderr = ""
    except Exception as e:
        elapsed = time.time() - t0
        exit_code = -998
        stdout = str(e)
        stderr = ""

    # Parse summary from stdout
    summary_line = ""
    n_fails = -1
    for line in stdout.split('\n'):
        if '_TEST:' in line:
            summary_line = line.strip()
        m = re.search(r'fails=(\d+)', line)
        if m:
            n_fails = int(m.group(1))

    # Load results
    disp_data = None
    react_data = None
    disp_file = os.path.join(WORK_DIR, prefix + '_disp.txt')
    react_file = os.path.join(WORK_DIR, prefix + '_react.txt')

    if os.path.exists(disp_file):
        try:
            disp_data = np.loadtxt(disp_file)
        except:
            pass
    if os.path.exists(react_file):
        try:
            react_data = np.loadtxt(react_file)
        except:
            pass

    return {
        'version': version_name,
        'test': test_name,
        'exit_code': exit_code,
        'elapsed': elapsed,
        'stdout': stdout,
        'stderr': stderr,
        'summary': summary_line,
        'n_fails': n_fails,
        'disp': disp_data,
        'react': react_data,
    }

def compute_stress_strain(result, test_name):
    """Compute stress and strain from displacement and reaction data."""
    disp = result['disp']
    react = result['react']
    if disp is None or react is None:
        return None, None, None, None

    time_arr = disp[:, 0]
    n = min(len(disp), len(react))
    time_arr = time_arr[:n]
    disp = disp[:n]
    react = react[:n]

    if test_name in ['test_elem_comp', 'test_elem_tens', 'test_elem_cyclic']:
        # Uniaxial y-loading
        # disp: time, ux3, uy3, ux4, uy4
        # react: time, Fx1, Fy1, Fx2, Fy2
        uy = disp[:, 2]  # uy at node 3
        strain = uy / 1.0  # eps_yy = uy / L

        # Total reaction in y at bottom (nodes 1+2)
        Fy_total = react[:, 2] + react[:, 4]  # Fy1 + Fy2
        stress = -Fy_total / AREA  # sigma_yy (negate: reaction opposes applied)

        return time_arr, strain, stress, 'eps_yy', 'sigma_yy'

    elif test_name == 'test_elem_shear':
        # Shear loading
        # disp: time, ux3, uy3, ux4, uy4
        # react: time, Fx1, Fy1, Fx2, Fy2
        ux = disp[:, 1]  # ux at node 3
        gamma = ux / 1.0  # gamma_xy = ux / L

        # Total reaction in x at bottom
        Fx_total = react[:, 1] + react[:, 3]  # Fx1 + Fx2
        tau = -Fx_total / AREA  # tau_xy

        return time_arr, gamma, tau, 'gamma_xy', 'tau_xy'

    return None, None, None, None, None

def main():
    versions = find_opensees_versions()
    print(f"Found {len(versions)} OpenSees versions:")
    for name, path in versions:
        print(f"  {name}")
    print()

    # Run all tests with all versions
    all_results = {}
    for test in TESTS:
        all_results[test] = {}
        print(f"{'='*60}")
        print(f"Test: {TEST_LABELS[test]} ({test})")
        print(f"{'='*60}")

        for ver_name, ver_path in versions:
            print(f"  Running with {ver_name}...", end='', flush=True)
            r = run_test(ver_path, test, ver_name, test)
            all_results[test][ver_name] = r
            status = "OK" if r['exit_code'] == 0 else f"EXIT={r['exit_code']}"
            n_pts = len(r['disp']) if r['disp'] is not None else 0
            fail_str = f", fails={r['n_fails']}" if r['n_fails'] >= 0 else ""
            print(f" {status}, {r['elapsed']:.1f}s, {n_pts}pts{fail_str}")
            if r['summary']:
                print(f"    {r['summary']}")
        print()

    # ========== Generate comparison plots ==========
    print("\n" + "="*60)
    print("Generating comparison plots...")
    print("="*60)

    # Use distinct colors and styles for each version
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red',
              'tab:purple', 'tab:brown', 'tab:pink']
    styles = ['-', '--', '-.', ':', '-', '--', '-.']

    # ---- Plot 1: Stress-strain for each test ----
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    ax_list = [axes[0,0], axes[0,1], axes[1,0], axes[1,1]]

    for idx, test in enumerate(TESTS):
        ax = ax_list[idx]
        has_data = False

        for vi, (ver_name, _) in enumerate(versions):
            r = all_results[test][ver_name]
            result = compute_stress_strain(r, test)
            if result[0] is None:
                continue

            t, strain, stress, xlabel, ylabel = result

            # Filter out NaN/Inf
            valid = np.isfinite(strain) & np.isfinite(stress)
            strain = strain[valid]
            stress = stress[valid]

            # Convert to engineering units
            strain_pct = strain * 100  # percent
            stress_MPa = stress / 1e6  # MPa

            label = ver_name.replace('.exe', '')
            ax.plot(strain_pct, stress_MPa,
                    color=colors[vi % len(colors)],
                    linestyle=styles[vi % len(styles)],
                    linewidth=1.2, alpha=0.8, label=label)
            has_data = True

        ax.set_xlabel(f'应变 (%)', fontsize=11)
        ax.set_ylabel(f'应力 (MPa)', fontsize=11)
        ax.set_title(TEST_LABELS[test], fontsize=13)
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0, color='k', linewidth=0.5)
        ax.axvline(x=0, color='k', linewidth=0.5)
        if has_data:
            ax.legend(fontsize=7, loc='best')

    plt.tight_layout()
    plt.savefig(os.path.join(WORK_DIR, 'elem_test_comparison.png'), dpi=150)
    plt.close()
    print("Saved: elem_test_comparison.png")

    # ---- Plot 2: Summary table ----
    fig2, ax2 = plt.subplots(figsize=(14, 6))
    ax2.axis('off')

    col_labels = ['OpenSees版本'] + [TEST_LABELS[t] for t in TESTS]
    table_data = []

    for ver_name, _ in versions:
        row = [ver_name.replace('.exe', '')]
        for test in TESTS:
            r = all_results[test][ver_name]
            n_pts = len(r['disp']) if r['disp'] is not None else 0
            has_data = n_pts > 0
            if has_data:
                cell = f"OK ({n_pts}pts)"
                if r['n_fails'] > 0:
                    cell += f"\n{r['n_fails']}fails"
                if r['exit_code'] != 0:
                    cell += f"\nexit={r['exit_code']}"
            elif r['exit_code'] == -999:
                cell = "TIMEOUT"
            else:
                cell = f"NO DATA\nexit={r['exit_code']}"
            row.append(cell)
        table_data.append(row)

    table = ax2.table(cellText=table_data, colLabels=col_labels,
                      loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 2.0)

    # Color cells
    for i, row in enumerate(table_data):
        for j in range(1, len(row)):
            cell = table[i+1, j]
            text = row[j]
            if 'OK' in text and 'fail' not in text.lower():
                cell.set_facecolor('#c8e6c9')  # green
            elif 'OK' in text:
                cell.set_facecolor('#fff9c4')  # yellow
            elif 'CRASH' in text:
                cell.set_facecolor('#ffcdd2')  # red
            elif 'TIMEOUT' in text:
                cell.set_facecolor('#e1bee7')  # purple

    ax2.set_title('各版本OpenSees单元素测试结果汇总', fontsize=14, pad=20)
    plt.tight_layout()
    plt.savefig(os.path.join(WORK_DIR, 'elem_test_summary.png'), dpi=150)
    plt.close()
    print("Saved: elem_test_summary.png")

    # ---- Plot 3: Detailed cyclic comparison ----
    fig3, axes3 = plt.subplots(1, 2, figsize=(16, 7))

    # Left: stress-strain hysteresis
    ax = axes3[0]
    for vi, (ver_name, _) in enumerate(versions):
        r = all_results['test_elem_cyclic'][ver_name]
        result = compute_stress_strain(r, 'test_elem_cyclic')
        if result[0] is None:
            continue
        t, strain, stress, _, _ = result
        valid = np.isfinite(strain) & np.isfinite(stress)
        strain_pct = strain[valid] * 100
        stress_MPa = stress[valid] / 1e6
        label = ver_name.replace('.exe', '')
        ax.plot(strain_pct, stress_MPa,
                color=colors[vi % len(colors)],
                linestyle=styles[vi % len(styles)],
                linewidth=1.0, alpha=0.8, label=label)

    ax.set_xlabel('应变 ε_yy (%)', fontsize=12)
    ax.set_ylabel('应力 σ_yy (MPa)', fontsize=12)
    ax.set_title('循环测试 - 应力应变滞回', fontsize=13)
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='k', linewidth=0.5)
    ax.axvline(x=0, color='k', linewidth=0.5)
    ax.legend(fontsize=8)

    # Right: stress time history
    ax = axes3[1]
    for vi, (ver_name, _) in enumerate(versions):
        r = all_results['test_elem_cyclic'][ver_name]
        result = compute_stress_strain(r, 'test_elem_cyclic')
        if result[0] is None:
            continue
        t, strain, stress, _, _ = result
        valid = np.isfinite(strain) & np.isfinite(stress)
        stress_MPa = stress[valid] / 1e6
        t_valid = t[valid] if len(t) == len(stress) else t[:len(stress)][valid]
        label = ver_name.replace('.exe', '')
        ax.plot(t_valid, stress_MPa,
                color=colors[vi % len(colors)],
                linestyle=styles[vi % len(styles)],
                linewidth=1.0, alpha=0.8, label=label)

    ax.set_xlabel('时间', fontsize=12)
    ax.set_ylabel('应力 σ_yy (MPa)', fontsize=12)
    ax.set_title('循环测试 - 应力时程', fontsize=13)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(os.path.join(WORK_DIR, 'elem_test_cyclic_detail.png'), dpi=150)
    plt.close()
    print("Saved: elem_test_cyclic_detail.png")

    # ---- Print detailed comparison ----
    print("\n" + "="*60)
    print("详细结果对比")
    print("="*60)

    for test in TESTS:
        print(f"\n--- {TEST_LABELS[test]} ---")
        for ver_name, _ in versions:
            r = all_results[test][ver_name]
            result = compute_stress_strain(r, test)
            ver_short = ver_name.replace('.exe', '')

            if result[0] is not None:
                t, strain, stress, _, _ = result
                valid = np.isfinite(strain) & np.isfinite(stress)
                strain_v = strain[valid]
                stress_v = stress[valid] / 1e6

                print(f"  {ver_short:20s}: exit={r['exit_code']:3d}, "
                      f"pts={len(strain_v):4d}, "
                      f"strain=[{strain_v.min():.6f}, {strain_v.max():.6f}], "
                      f"stress=[{stress_v.min():.2f}, {stress_v.max():.2f}] MPa, "
                      f"fails={r['n_fails']}")
            else:
                print(f"  {ver_short:20s}: exit={r['exit_code']:3d}, NO DATA")

    print("\nDone!")

if __name__ == '__main__':
    main()
