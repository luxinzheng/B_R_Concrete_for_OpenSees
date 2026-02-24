"""Analyze the mechanism of the B&R wall capacity overshoot."""
import numpy as np

da = np.loadtxt('vwall_adina_manual.txt')
dr = np.loadtxt('vwall_ref_manual.txt')
disp_a = da[:,1]*1000; react_a = -np.sum(da[:,2:], axis=1)/1000
disp_r = dr[:,1]*1000; react_r = -np.sum(dr[:,2:], axis=1)/1000

# 1. Track peak force ratio at each cycle's peak
def get_cycle_peaks(disp, force):
    """Get the peak force at each displacement reversal."""
    peaks = []
    for i in range(1, len(force)-1):
        if (disp[i] > disp[i-1] and disp[i] > disp[i+1]) or \
           (disp[i] < disp[i-1] and disp[i] < disp[i+1]):
            if abs(force[i]) > 20:
                peaks.append((disp[i], force[i]))
    return peaks

peaks_a = get_cycle_peaks(disp_a, react_a)
peaks_r = get_cycle_peaks(disp_r, react_r)

# 2. Unloading stiffness at each reversal
print("=== Force evolution at displacement reversals ===")
print(f"{'Cycle':>5s}  {'d_BR':>7s}  {'F_BR':>8s}  {'d_Ref':>7s}  {'F_Ref':>8s}  {'Ratio':>6s}")
for i in range(min(len(peaks_a), len(peaks_r))):
    d_a, f_a = peaks_a[i]
    d_r, f_r = peaks_r[i]
    ratio = f_a/f_r if abs(f_r)>1 else 0
    print(f"  {i+1:3d}  {d_a:7.1f}  {f_a:8.1f}  {d_r:7.1f}  {f_r:8.1f}  {ratio:6.2f}")

# 3. Measure unloading stiffness at each reversal
print()
print("=== Unloading stiffness (kN/mm) at each reversal ===")
print("  (dF/dd over the first 3 steps after reversal)")
for model_name, disp_m, react_m in [("B&R", disp_a, react_a), ("Ref", disp_r, react_r)]:
    print(f"  {model_name}:")
    reversals = []
    for i in range(2, len(disp_m)-3):
        # Detect reversal: displacement direction changes
        if (disp_m[i]-disp_m[i-1]) * (disp_m[i+1]-disp_m[i]) < 0:
            dd = disp_m[i+3] - disp_m[i]
            df = react_m[i+3] - react_m[i]
            if abs(dd) > 0.1:
                k = df/dd
                reversals.append((disp_m[i], react_m[i], k))
    for j, (d, f, k) in enumerate(reversals[:10]):
        print(f"    rev {j+1:2d}: d={d:7.1f}mm F={f:8.1f}kN  K_unload={k:8.1f} kN/mm")

# 4. Check: Is the ratio growing over cycles?
print()
print("=== Ratio growth with cycle number ===")
pos_ratios = []
neg_ratios = []
for i in range(min(len(peaks_a), len(peaks_r))):
    d_a, f_a = peaks_a[i]
    d_r, f_r = peaks_r[i]
    if abs(f_r) > 10:
        ratio = f_a/f_r
        if f_a > 0:
            pos_ratios.append(ratio)
        else:
            neg_ratios.append(ratio)

if pos_ratios:
    print(f"  Positive ratios: {' -> '.join(f'{r:.2f}' for r in pos_ratios)}")
if neg_ratios:
    print(f"  Negative ratios: {' -> '.join(f'{r:.2f}' for r in neg_ratios)}")
