import csv, sys
sys.stdout.reconfigure(line_buffering=True)

rows = []
with open('test_transition.csv') as f:
    reader = csv.DictReader(f)
    for r in reader:
        rows.append(r)

print('Transition detail around zero crossing:')
print(f'{"step":>4s}  {"eps_xx":>12s}  {"sig_xx":>12s}  {"dsig_xx":>12s}')
for r in rows:
    sig = float(r['sig_xx'])
    dsig = float(r['dsig_xx'])
    eps = float(r['eps_xx'])
    step = int(r['step'])
    if abs(dsig) > 0.5 or (step >= 48 and step <= 65):
        print(f'{step:4d}  {eps:+12.6f}  {sig:+12.4f}  {dsig:+12.4f}')
