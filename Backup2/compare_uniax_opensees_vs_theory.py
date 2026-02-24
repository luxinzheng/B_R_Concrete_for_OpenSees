# -*- coding: utf-8 -*-
"""
单轴往复 110 步：OpenSees (TCL) 计算结果 vs 理论/Fortran 结果对比。
理论参考：uniax_cyclic_110steps.txt（Fortran forumat 单点本构输出）。
"""
import numpy as np
from pathlib import Path

def main():
    base = Path(__file__).parent
    # OpenSees 结果
    try:
        data_os = np.loadtxt(base / 'opensees_uniax_cyclic_110steps.txt', skiprows=1)
        step_os = data_os[:, 0].astype(int)
        strain_os = data_os[:, 1]
        stress_os = data_os[:, 2]
        phase_os = data_os[:, 3].astype(int)
    except Exception as e:
        print('请先运行: openSees.exe test_uniax_cyclic_opensees.tcl 生成 opensees_uniax_cyclic_110steps.txt')
        print(e)
        return

    # 理论/Fortran 结果
    try:
        data_th = np.loadtxt(base / 'uniax_cyclic_110steps.txt', skiprows=1)
        step_th = data_th[:, 0].astype(int)
        strain_th = data_th[:, 1]
        stress_th = data_th[:, 2]
        phase_th = data_th[:, 3].astype(int)
    except Exception as e:
        print('请确保存在 uniax_cyclic_110steps.txt（Fortran 单轴往复输出）')
        print(e)
        return

    n = min(len(step_os), len(step_th))
    step_os, strain_os, stress_os = step_os[:n], strain_os[:n], stress_os[:n]
    step_th, strain_th, stress_th = step_th[:n], strain_th[:n], stress_th[:n]

    err_abs = stress_os - stress_th
    err_rel = np.zeros(n)
    for i in range(n):
        if abs(stress_th[i]) > 1e3:
            err_rel[i] = 100.0 * err_abs[i] / abs(stress_th[i])
        else:
            err_rel[i] = 0.0 if abs(err_abs[i]) < 1.0 else 100.0

    print('单轴往复 110 步：OpenSees vs 理论(Fortran)')
    print('step  strain(1e-6)  theory(MPa)  OpenSees(MPa)  err_abs(Pa)   err_rel(%)  phase')
    print('-' * 90)
    for i in range(0, n, max(1, n // 25)):
        print('{:4d}  {:12.4f}  {:11.4f}  {:12.4f}  {:12.2f}  {:10.2f}  {:5d}'.format(
            step_th[i], strain_th[i]*1e6, stress_th[i]/1e6, stress_os[i]/1e6,
            err_abs[i], err_rel[i], phase_th[i]))
    print('-' * 90)
    print('最大绝对误差 (Pa):', np.max(np.abs(err_abs)))
    print('最大相对误差 (%):', np.max(np.abs(err_rel)))
    print('均方根应力误差 (Pa):', np.sqrt(np.mean(err_abs**2)))
    imax = np.argmax(np.abs(err_rel))
    print('最大相对误差步: step={}, strain={:.4e}, rel%={:.2f}'.format(
        step_th[imax], strain_th[imax], err_rel[imax]))

    # 写表
    out_table = base / 'compare_uniax_opensees_vs_theory.txt'
    with open(out_table, 'w', encoding='utf-8') as f:
        f.write('step  strain_xx       stress_theory   stress_OpenSees  err_abs(Pa)   err_rel_pct  phase\n')
        for i in range(n):
            f.write('{:4d}  {:16.8e}  {:16.8e}  {:16.8e}  {:14.4f}  {:10.2f}  {:5d}\n'.format(
                step_th[i], strain_th[i], stress_th[i], stress_os[i], err_abs[i], err_rel[i], phase_th[i]))
    print('已写入:', out_table)

    # 应力-应变对比图
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 1, figsize=(8, 5))
        ax.plot(strain_th * 1e6, stress_th / 1e6, 'b-', lw=1.5, label='理论/Fortran')
        ax.plot(strain_os * 1e6, stress_os / 1e6, 'r--', lw=1, alpha=0.8, label='OpenSees')
        ax.set_xlabel('Strain (1e-6)')
        ax.set_ylabel('Stress (MPa)')
        ax.set_title('Uniaxial cyclic: OpenSees vs Theory (Fortran)')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.axhline(0, color='k', linewidth=0.5)
        ax.axvline(0, color='k', linewidth=0.5)
        fig.tight_layout()
        fig.savefig(base / 'plot_uniax_cyclic_opensees_vs_theory.png', dpi=150)
        plt.close()
        print('已保存: plot_uniax_cyclic_opensees_vs_theory.png')
    except Exception as e:
        print('出图跳过:', e)

if __name__ == '__main__':
    main()
