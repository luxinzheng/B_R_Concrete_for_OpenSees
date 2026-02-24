import csv, sys
sys.stdout.reconfigure(line_buffering=True)

rows = []
with open('test_cyclic_nr.csv') as f:
    reader = csv.DictReader(f)
    for r in reader:
        rows.append(r)

print(f'Total rows: {len(rows)}')

# Find high NR iterations
print('\nSteps with NR iterations >= 10:')
for r in rows:
    it = int(r['nr_iters'])
    eps = float(r['eps_xx'])
    sig = float(r['sig_xx'])
    if it >= 10:
        print(f'  step={r["step"]:>4s} eps={eps:+.6f} sig={sig:+10.4f} iters={it}')

# Find large stress jumps
print('\nLarge stress jumps (>1 MPa):')
prev_sig = 0.0
for r in rows:
    sig = float(r['sig_xx'])
    dsig = abs(sig - prev_sig)
    if dsig > 1.0:
        eps = float(r['eps_xx'])
        print(f'  step={r["step"]:>4s} eps={eps:+.6f} sig={sig:+10.4f} dsig={dsig:.4f}')
    prev_sig = sig

# Symmetry: stress at eps = +/- values
print('\nStress at matching +/- strain:')
for target in [0.0005, 0.001, 0.0015, 0.002, 0.0025, 0.003]:
    pos_vals = [float(r['sig_xx']) for r in rows if abs(float(r['eps_xx']) - target) < 1e-6]
    neg_vals = [float(r['sig_xx']) for r in rows if abs(float(r['eps_xx']) + target) < 1e-6]
    if pos_vals and neg_vals:
        print(f'  eps=+/-{target:.4f}:')
        for i, (p, n) in enumerate(zip(pos_vals, neg_vals)):
            print(f'    cycle {i}: sig(+)={p:+.3f}  sig(-)={n:+.3f}  ratio={abs(n/p):.2f}' if abs(p) > 0.01 else f'    cycle {i}: sig(+)={p:+.6f}  sig(-)={n:+.3f}')

# Check eps_yy (Poisson strain)
print('\neps_yy at peak strains:')
for r in rows:
    eps = float(r['eps_xx'])
    if abs(abs(eps) - 0.003) < 1e-6:
        print(f'  step={r["step"]:>4s} eps_xx={eps:+.6f} sig_xx={float(r["sig_xx"]):+.3f} eps_yy={float(r["eps_yy"]):+.8f}')
