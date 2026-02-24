import numpy as np

da = np.loadtxt('vwall_adina_manual.txt')
dr = np.loadtxt('vwall_ref_manual.txt')

disp_a = da[:,1]*1000; react_a = -np.sum(da[:,2:], axis=1)/1000
disp_r = dr[:,1]*1000; react_r = -np.sum(dr[:,2:], axis=1)/1000

# Key observation: at d=2mm (step 0), ratio=1.02 (nearly identical)
# But at d=-2mm (step 8), ratio=1.38 already
# And during unloading d=1mm (step 5), ratio drops to 0.35

# Check: during first unloading (from +4mm back to 0), what happens?
print("=== First unloading from +4mm to -4mm ===")
for i in range(2, 11):
    print(f"  step {i:2d}  d_BR={disp_a[i]:6.1f}  F_BR={react_a[i]:8.1f}  "
          f"d_Ref={disp_r[i]:6.1f}  F_Ref={react_r[i]:8.1f}  "
          f"dF_BR={react_a[i]-react_a[i-1]:+7.1f}  dF_Ref={react_r[i]-react_r[i-1]:+7.1f}")

print()
print("=== Unloading stiffness comparison ===")
# First unloading: step 2->6 (d=4->0mm)
dF_a = react_a[6] - react_a[2]; dd_a = disp_a[6] - disp_a[2]
dF_r = react_r[6] - react_r[2]; dd_r = disp_r[6] - disp_r[2]
print(f"  First unload (4->0mm): BR dF/dd = {dF_a/dd_a:.1f} kN/mm,  Ref dF/dd = {dF_r/dd_r:.1f} kN/mm")

# Look at single-element uniaxial data for comparison
print()
print("=== Uniaxial compression comparison ===")
dc = np.loadtxt('vuniax_adina_comp_disp.txt')
rc = np.loadtxt('vuniax_adina_comp_react.txt')
dc_r = np.loadtxt('vuniax_ref_comp_disp.txt')
rc_r = np.loadtxt('vuniax_ref_comp_react.txt')

area = 1.0 * 0.1  # width * thickness
strain_a = dc[:,1]  # node displacement = strain (L=1m)
stress_a = np.sum(rc[:,1:], axis=1) / area / 1e6
strain_r = dc_r[:,1]
stress_r = np.sum(rc_r[:,1:], axis=1) / area / 1e6

# Compare at key strain levels
for target_eps in [-0.001, -0.002, -0.003, -0.004]:
    idx_a = np.argmin(np.abs(strain_a - target_eps))
    idx_r = np.argmin(np.abs(strain_r - target_eps))
    sig_a = stress_a[idx_a]
    sig_r = stress_r[idx_r]
    ratio = sig_a/sig_r if abs(sig_r) > 0.01 else 0
    print(f"  eps={target_eps:.4f}: BR sig={sig_a:8.2f} MPa, Ref sig={sig_r:8.2f} MPa, ratio={ratio:.3f}")

# Check uniaxial tension
print()
print("=== Uniaxial tension comparison ===")
dt = np.loadtxt('vuniax_adina_tens_disp.txt')
rt = np.loadtxt('vuniax_adina_tens_react.txt')
dt_r = np.loadtxt('vuniax_ref_tens_disp.txt')
rt_r = np.loadtxt('vuniax_ref_tens_react.txt')

strain_t_a = dt[:,1]
stress_t_a = np.sum(rt[:,1:], axis=1) / area / 1e6
strain_t_r = dt_r[:,1]
stress_t_r = np.sum(rt_r[:,1:], axis=1) / area / 1e6

for target_eps in [0.0001, 0.0002, 0.0005, 0.001]:
    idx_a = np.argmin(np.abs(strain_t_a - target_eps))
    idx_r = np.argmin(np.abs(strain_t_r - target_eps))
    sig_a = stress_t_a[idx_a]
    sig_r = stress_t_r[idx_r]
    ratio = sig_a/sig_r if abs(sig_r) > 0.001 else 0
    print(f"  eps={target_eps:.5f}: BR sig={sig_a:8.4f} MPa, Ref sig={sig_r:8.4f} MPa, ratio={ratio:.3f}")
