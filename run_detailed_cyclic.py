"""
Detailed cyclic test comparison between OpenSees versions.
Saves per-version data and creates zoomed-in comparison plots.
"""
import subprocess, os, numpy as np, matplotlib.pyplot as plt, matplotlib, shutil
matplotlib.rcParams['font.family'] = ['Microsoft YaHei', 'SimHei', 'sans-serif']

WORK_DIR = r'e:\Basic\Concrete_Model\ADINA'
VERSIONS = [
    ('OpenSees-0917.exe', '0917 (E_sec_8a fix)'),
    ('OpenSees-0929.exe', '0929 (baseline+NaN)'),
    ('OpenSees-0943.exe', '0943 (latest)'),
    ('OpenSees.exe',      'OpenSees (original)'),
]
THICKNESS = 0.1

def run_and_save(exe_name, tag):
    exe = os.path.join(WORK_DIR, exe_name)
    tcl = os.path.join(WORK_DIR, 'test_elem_cyclic.tcl')

    for f in ['cyclic_disp.txt', 'cyclic_react.txt']:
        p = os.path.join(WORK_DIR, f)
        if os.path.exists(p): os.remove(p)

    r = subprocess.run([exe, tcl], capture_output=True, text=True, timeout=60, cwd=WORK_DIR)

    disp_file = os.path.join(WORK_DIR, 'cyclic_disp.txt')
    react_file = os.path.join(WORK_DIR, 'cyclic_react.txt')

    # Save copies with version tag
    for src, suffix in [(disp_file, '_disp'), (react_file, '_react')]:
        if os.path.exists(src) and os.path.getsize(src) > 0:
            dst = os.path.join(WORK_DIR, f'cyclic_{tag}{suffix}.txt')
            shutil.copy2(src, dst)

    # Load data
    disp = np.loadtxt(disp_file) if os.path.exists(disp_file) and os.path.getsize(disp_file) > 0 else None
    react = np.loadtxt(react_file) if os.path.exists(react_file) and os.path.getsize(react_file) > 0 else None
    return disp, react, r.stdout

def compute(disp, react):
    n = min(len(disp), len(react))
    time = disp[:n, 0]
    uy = disp[:n, 2]
    strain = uy / 1.0
    Fy_total = react[:n, 2] + react[:n, 4]
    stress = -Fy_total / THICKNESS
    return time, strain, stress

# ---- Run tests ----
results = {}
for exe_name, label in VERSIONS:
    tag = exe_name.replace('.exe','').replace('-','_')
    print(f"Running {exe_name}...", end='', flush=True)
    disp, react, stdout = run_and_save(exe_name, tag)
    if disp is not None and react is not None:
        t, eps, sig = compute(disp, react)
        results[label] = (t, eps, sig)
        print(f" OK, {len(t)} pts, stress range: [{sig.min()/1e6:.2f}, {sig.max()/1e6:.2f}] MPa")
    else:
        print(" NO DATA")

# ---- Target displacement path (for reference) ----
target_vals = np.array([0, -0.001, 0, 0.00015, 0, -0.002, 0, 0.0003, 0, -0.003, 0])
target_time = np.arange(len(target_vals))

# ---- PLOT 1: Full stress-strain comparison ----
fig, axes = plt.subplots(2, 2, figsize=(18, 14))
colors = ['tab:blue', 'tab:red', 'tab:green', 'tab:purple']

# (a) Full hysteresis
ax = axes[0, 0]
for i, (label, (t, eps, sig)) in enumerate(results.items()):
    ax.plot(eps*100, sig/1e6, color=colors[i], linewidth=1.2, alpha=0.8, label=label)
ax.set_xlabel('应变 ε_yy (%)')
ax.set_ylabel('应力 σ_yy (MPa)')
ax.set_title('循环测试 - 完整应力应变曲线')
ax.grid(True, alpha=0.3); ax.legend(fontsize=9)
ax.axhline(y=0, color='k', lw=0.5); ax.axvline(x=0, color='k', lw=0.5)

# (b) Zoom into tension region
ax = axes[0, 1]
for i, (label, (t, eps, sig)) in enumerate(results.items()):
    mask = eps > -0.0002  # zoom into tension and near-zero region
    ax.plot(eps[mask]*100, sig[mask]/1e6, color=colors[i], linewidth=1.5, alpha=0.8, label=label)
ax.axhline(y=2.07, color='gray', ls='--', lw=1, label='ft = 2.07 MPa')
ax.axhline(y=0, color='k', lw=0.5)
ax.axvline(x=0, color='k', lw=0.5)
ax.set_xlabel('应变 ε_yy (%)')
ax.set_ylabel('应力 σ_yy (MPa)')
ax.set_title('循环测试 - 拉伸区域放大')
ax.grid(True, alpha=0.3); ax.legend(fontsize=9)

# (c) Stress time history
ax = axes[1, 0]
for i, (label, (t, eps, sig)) in enumerate(results.items()):
    ax.plot(t, sig/1e6, color=colors[i], linewidth=1.0, alpha=0.8, label=label)
ax.set_xlabel('时间')
ax.set_ylabel('应力 σ_yy (MPa)')
ax.set_title('循环测试 - 应力时程')
ax.grid(True, alpha=0.3); ax.legend(fontsize=9)

# Add phase labels
phase_labels = ['压-0.1%', '卸载', '拉+0.015%', '卸载', '压-0.2%',
                '卸载', '拉+0.03%', '卸载', '压-0.3%', '卸载']
for j in range(10):
    ax.axvline(x=j+1, color='gray', ls=':', lw=0.5)
    ax.text(j+0.5, ax.get_ylim()[1]*0.95, phase_labels[j], fontsize=6,
            ha='center', va='top', rotation=45)

# (d) Strain time history (sanity check)
ax = axes[1, 1]
for i, (label, (t, eps, sig)) in enumerate(results.items()):
    ax.plot(t, eps*100, color=colors[i], linewidth=1.0, alpha=0.8, label=label)
ax.plot(target_time, target_vals*100, 'k--', lw=1.5, label='目标')
ax.set_xlabel('时间')
ax.set_ylabel('应变 ε_yy (%)')
ax.set_title('循环测试 - 应变时程')
ax.grid(True, alpha=0.3); ax.legend(fontsize=9)

plt.tight_layout()
plt.savefig(os.path.join(WORK_DIR, 'cyclic_version_comparison.png'), dpi=150)
plt.close()
print("\nSaved: cyclic_version_comparison.png")

# ---- PLOT 2: Phase-by-phase stress comparison at phase boundaries ----
fig2, ax2 = plt.subplots(figsize=(14, 6))
phase_times = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
phase_names = ['comp -0.1%', 'release 1', 'tens +0.015%', 'release 2',
               'comp -0.2%', 'release 3', 'tens +0.03%', 'release 4',
               'comp -0.3%', 'release 5']

x = np.arange(len(phase_times))
width = 0.2
for i, (label, (t, eps, sig)) in enumerate(results.items()):
    stresses_at_phases = []
    for pt in phase_times:
        idx = np.argmin(np.abs(t - pt))
        stresses_at_phases.append(sig[idx] / 1e6)
    ax2.bar(x + i*width, stresses_at_phases, width, label=label, alpha=0.8)

ax2.set_xticks(x + width * (len(results)-1) / 2)
ax2.set_xticklabels(phase_names, rotation=45, ha='right', fontsize=8)
ax2.set_ylabel('应力 σ_yy (MPa)', fontsize=12)
ax2.set_title('各版本在各阶段结束时的应力对比', fontsize=14)
ax2.grid(True, alpha=0.3, axis='y')
ax2.legend(fontsize=9)
plt.tight_layout()
plt.savefig(os.path.join(WORK_DIR, 'cyclic_phase_comparison.png'), dpi=150)
plt.close()
print("Saved: cyclic_phase_comparison.png")

# ---- Print detailed per-phase comparison ----
print("\n" + "="*80)
print("各版本各阶段应力对比 (MPa)")
print("="*80)
header = f"{'阶段':<20s}"
for label in results:
    header += f"  {label:>18s}"
print(header)
print("-"*80)
for j, (pt, pn) in enumerate(zip(phase_times, phase_names)):
    row = f"{pn:<20s}"
    for label, (t, eps, sig) in results.items():
        idx = np.argmin(np.abs(t - pt))
        row += f"  {sig[idx]/1e6:>18.4f}"
    print(row)
