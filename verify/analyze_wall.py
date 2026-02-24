import numpy as np

da = np.loadtxt('vwall_adina_manual.txt')
dr = np.loadtxt('vwall_ref_manual.txt')

disp_a = da[:,1]*1000; react_a = -np.sum(da[:,2:], axis=1)/1000
disp_r = dr[:,1]*1000; react_r = -np.sum(dr[:,2:], axis=1)/1000

def get_peaks(disp, force):
    pos_d, pos_f, neg_d, neg_f = [], [], [], []
    for i in range(1, len(force)-1):
        if force[i] >= force[i-1] and force[i] >= force[i+1] and force[i] > 10:
            pos_d.append(disp[i]); pos_f.append(force[i])
        if force[i] <= force[i-1] and force[i] <= force[i+1] and force[i] < -10:
            neg_d.append(disp[i]); neg_f.append(force[i])
    return pos_d, pos_f, neg_d, neg_f

ppd_a, ppf_a, npd_a, npf_a = get_peaks(disp_a, react_a)
ppd_r, ppf_r, npd_r, npf_r = get_peaks(disp_r, react_r)

print("=== Positive Peak Envelope ===")
print(f"  {'d(mm)':>8s}  {'BR(kN)':>10s}  {'Ref(kN)':>10s}  {'Ratio':>8s}")
for i in range(min(len(ppf_a), len(ppf_r))):
    ratio = ppf_a[i]/ppf_r[i] if abs(ppf_r[i])>1 else 0
    print(f"  {ppd_a[i]:8.1f}  {ppf_a[i]:10.1f}  {ppf_r[i]:10.1f}  {ratio:8.2f}")

print()
print("=== Negative Peak Envelope ===")
print(f"  {'d(mm)':>8s}  {'BR(kN)':>10s}  {'Ref(kN)':>10s}  {'Ratio':>8s}")
for i in range(min(len(npf_a), len(npf_r))):
    ratio = npf_a[i]/npf_r[i] if abs(npf_r[i])>1 else 0
    print(f"  {npd_a[i]:8.1f}  {npf_a[i]:10.1f}  {npf_r[i]:10.1f}  {ratio:8.2f}")

# Also look at initial stiffness
print()
print("=== Initial loading comparison (first 15 steps) ===")
print(f"  {'step':>4s}  {'d_BR':>8s}  {'F_BR':>10s}  {'d_Ref':>8s}  {'F_Ref':>10s}  {'F_ratio':>8s}")
for i in range(min(15, len(react_a), len(react_r))):
    ratio = react_a[i]/react_r[i] if abs(react_r[i])>1 else 0
    print(f"  {i:4d}  {disp_a[i]:8.2f}  {react_a[i]:10.2f}  {disp_r[i]:8.2f}  {react_r[i]:10.2f}  {ratio:8.2f}")

# Secant stiffness at different displacement levels
print()
print("=== Secant stiffness comparison ===")
for target_d in [2, 5, 10, 15, 20, 25]:
    idx_a = np.argmin(np.abs(disp_a - target_d))
    idx_r = np.argmin(np.abs(disp_r - target_d))
    if abs(disp_a[idx_a]) > 0.1 and abs(disp_r[idx_r]) > 0.1:
        k_a = react_a[idx_a] / disp_a[idx_a]
        k_r = react_r[idx_r] / disp_r[idx_r]
        print(f"  d~{target_d:2d}mm: BR K={k_a:7.2f} kN/mm  Ref K={k_r:7.2f} kN/mm  ratio={k_a/k_r:.2f}")
