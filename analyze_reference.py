"""
Analyze reference model behavior characteristics for standalone test targets.
"""
import numpy as np
import sys
sys.stdout.reconfigure(line_buffering=True)

# Load reference data
ref_rows = []
with open(r'e:\Basic\Concrete_Model\ADINA\参考模型.txt', 'r') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2:
            try:
                d, force = float(parts[0]), float(parts[1])
                if abs(d) > 1e-6 or abs(force) > 100:
                    ref_rows.append([d*1000, force/1000])  # mm, kN
            except ValueError:
                pass

ref = np.array(ref_rows)
print(f"Reference data: {len(ref)} points")
print(f"Displacement range: [{ref[:,0].min():.1f}, {ref[:,0].max():.1f}] mm")
print(f"Force range: [{ref[:,1].min():.1f}, {ref[:,1].max():.1f}] kN")

# Identify cycle boundaries
print("\n" + "="*60)
print("CYCLE ANALYSIS")
print("="*60)

# Find peaks (direction changes)
peaks = [0]
for i in range(1, len(ref)-1):
    if (ref[i,0] > ref[i-1,0] and ref[i,0] > ref[i+1,0]) or \
       (ref[i,0] < ref[i-1,0] and ref[i,0] < ref[i+1,0]):
        peaks.append(i)
peaks.append(len(ref)-1)

print("\nPeak points (displacement reversals):")
for i in peaks:
    print(f"  pt {i:3d}: d={ref[i,0]:+8.1f} mm, F={ref[i,1]:+8.1f} kN")

# Key behavioral metrics
print("\n" + "="*60)
print("KEY BEHAVIORAL METRICS")
print("="*60)

# Initial stiffness
K0_pos = ref[0,1] / ref[0,0]  # kN/mm
print(f"\nInitial stiffness:")
print(f"  Push(+): d=+{ref[0,0]:.0f}mm → F={ref[0,1]:.0f} kN → K={K0_pos:.1f} kN/mm")

# Symmetry check at matched displacements
print("\nSymmetry check (|F(+d)| vs |F(-d)|):")
for d_target in [1, 2, 3, 4, 5, 6, 7]:
    pos = [abs(ref[i,1]) for i in range(len(ref)) if abs(ref[i,0] - d_target) < 0.1]
    neg = [abs(ref[i,1]) for i in range(len(ref)) if abs(ref[i,0] + d_target) < 0.1]
    if pos and neg:
        avg_pos = np.mean(pos)
        avg_neg = np.mean(neg)
        ratio = avg_neg / avg_pos if avg_pos > 0 else 0
        print(f"  d=±{d_target}mm: |F(+)|={avg_pos:.1f}, |F(-)|={avg_neg:.1f}, ratio={ratio:.3f}")

# Unloading stiffness
print("\nUnloading stiffness (from peak to zero-crossing):")
for i in range(len(ref)-1):
    d0, d1 = ref[i,0], ref[i+1,0]
    f0, f1 = ref[i,1], ref[i+1,1]
    # Check if this is an unloading from positive peak through zero
    if d0 > 3 and d1 < d0 and d1 > d0 - 1.5:
        K = (f1 - f0) / (d1 - d0)
        print(f"  d={d0:+.0f}→{d1:+.0f}mm: K={K:.1f} kN/mm (F: {f0:.0f}→{f1:.0f})")

# Residual force at zero displacement
print("\nResidual force at d≈0:")
for i in range(len(ref)):
    if abs(ref[i,0]) < 0.01:
        print(f"  pt {i}: d={ref[i,0]:.4f}mm, F={ref[i,1]:.1f} kN")

# Envelope curve
print("\nEnvelope (first occurrence at each displacement):")
seen_pos = set()
seen_neg = set()
for i in range(len(ref)):
    d = ref[i,0]
    f = ref[i,1]
    d_int = int(round(d))
    if d_int > 0 and d_int not in seen_pos:
        seen_pos.add(d_int)
        print(f"  d=+{d_int}mm: F={f:.1f} kN (first)")
    if d_int < 0 and d_int not in seen_neg:
        seen_neg.add(d_int)
        print(f"  d={d_int}mm: F={f:.1f} kN (first)")

# Stiffness degradation per cycle
print("\nStiffness degradation (secant at peak displacements):")
for i in range(len(ref)):
    d = ref[i,0]
    f = ref[i,1]
    if abs(d) > 0.5:
        K_sec = f / d
        if abs(d - 4) < 0.1 or abs(d - 7) < 0.1 or abs(d + 4) < 0.1 or abs(d + 7) < 0.1 or \
           abs(d - 11) < 0.1 or abs(d + 11) < 0.1:
            print(f"  pt {i}: d={d:+.0f}mm, F={f:+.0f} kN, K_sec={K_sec:.1f} kN/mm")
