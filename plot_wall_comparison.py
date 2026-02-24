"""Generate wall hysteresis comparison plots from existing data files."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

WORKDIR = r'E:\Basic\Concrete_Model\ADINA'
os.chdir(WORKDIR)

def load_manual_wall(fname, force_limit=500):
    if not os.path.exists(fname):
        return None
    fd = []
    with open(fname) as f:
        for line in f:
            p = line.strip().split()
            if len(p) < 7: continue
            if 'SKIP' in p: continue
            try:
                disp_mm = float(p[1]) * 1000
                shear_kn = -sum(float(p[i]) for i in range(2, 7)) / 1000
                if abs(shear_kn) < force_limit:
                    fd.append([disp_mm, shear_kn])
            except: pass
    return np.array(fd) if fd else None

def load_ref():
    rows = []
    try:
        with open(os.path.join(WORKDIR, '参考模型.txt')) as f:
            for line in f:
                p = line.strip().split()
                if len(p) >= 2:
                    try: rows.append([float(p[0])*1000, float(p[1])/1000])
                    except: pass
    except: pass
    return np.array(rows) if rows else None

def extract_envelope(fd):
    if fd is None or len(fd) == 0: return None, None
    d, f = fd[:,0], fd[:,1]
    pos, neg = [], []
    for i in range(1, len(d)-1):
        if d[i] > d[i-1] and d[i] > d[i+1] and d[i] > 0.5:
            pos.append([d[i], f[i]])
        elif d[i] < d[i-1] and d[i] < d[i+1] and d[i] < -0.5:
            neg.append([d[i], f[i]])
    return (np.array(pos) if pos else None, np.array(neg) if neg else None)

# Load data
fd_new = load_manual_wall('wall_2025_manual.txt')
if fd_new is None:
    fd_new = load_manual_wall('disp_manual.txt')
ref = load_ref()

tag = '2025'
print(f'New ({tag}): {len(fd_new) if fd_new is not None else 0} points')
print(f'Reference:  {len(ref) if ref is not None else 0} points')

# Plot 1: Hysteresis
fig, axes = plt.subplots(1, 2, figsize=(16, 7))

ax = axes[0]
ax.set_title(f'Hysteresis: ADINA-{tag} vs Reference', fontsize=13, fontweight='bold')
if ref is not None:
    ax.plot(ref[:,0], ref[:,1], 'b-', lw=0.7, alpha=0.6, label='Reference')
if fd_new is not None:
    ax.plot(fd_new[:,0], fd_new[:,1], 'r-', lw=0.9, label=f'ADINA-{tag}')
ax.set_xlabel('Displacement (mm)', fontsize=11)
ax.set_ylabel('Base Shear (kN)', fontsize=11)
ax.legend(fontsize=10, loc='upper left')
ax.grid(True, alpha=0.3)
ax.axhline(0, color='k', lw=0.5); ax.axvline(0, color='k', lw=0.5)

ax = axes[1]
ax.set_title(f'Detail View', fontsize=13, fontweight='bold')
if ref is not None:
    ax.plot(ref[:,0], ref[:,1], 'b-', lw=0.7, alpha=0.6, label='Reference')
if fd_new is not None:
    ax.plot(fd_new[:,0], fd_new[:,1], 'r-', lw=0.9, label=f'ADINA-{tag}')
    xlim = max(abs(fd_new[:,0].min()), abs(fd_new[:,0].max())) * 1.15
    ax.set_xlim([-xlim, xlim])
ax.set_xlabel('Displacement (mm)', fontsize=11)
ax.set_ylabel('Base Shear (kN)', fontsize=11)
ax.legend(fontsize=10, loc='upper left')
ax.grid(True, alpha=0.3)
ax.axhline(0, color='k', lw=0.5); ax.axvline(0, color='k', lw=0.5)

plt.tight_layout()
plt.savefig(f'wall_hysteresis_{tag}.png', dpi=150)
print(f'Saved: wall_hysteresis_{tag}.png')

# Plot 2: Envelope
fig2, axes2 = plt.subplots(1, 2, figsize=(16, 7))

ax = axes2[0]
ax.set_title('Envelope Comparison', fontsize=13, fontweight='bold')
for label, fd, color, marker in [
    ('Reference', ref, 'b', 'o'),
    (f'ADINA-{tag}', fd_new, 'r', 's')]:
    if fd is not None:
        ax.plot(fd[:,0], fd[:,1], color=color, lw=0.3, alpha=0.15)
        pos, neg = extract_envelope(fd)
        if pos is not None:
            ax.plot(pos[:,0], pos[:,1], color=color, marker=marker, ms=5, lw=1.5, label=f'{label} (+)')
        if neg is not None:
            ax.plot(neg[:,0], neg[:,1], color=color, marker=marker, ms=5, lw=1.5, ls='--', label=f'{label} (-)')
ax.set_xlabel('Displacement (mm)', fontsize=11)
ax.set_ylabel('Base Shear (kN)', fontsize=11)
ax.legend(fontsize=8, loc='best')
ax.grid(True, alpha=0.3)
ax.axhline(0, color='k', lw=0.5)

ax = axes2[1]
ax.set_title('Displacement History', fontsize=13, fontweight='bold')
if ref is not None:
    ax.plot(np.arange(len(ref)), ref[:,0], 'b-', lw=0.5, alpha=0.6, label='Reference')
if fd_new is not None:
    ax.plot(np.arange(len(fd_new)), fd_new[:,0], 'r-', lw=0.8, label=f'ADINA-{tag}')
ax.set_xlabel('Step', fontsize=11)
ax.set_ylabel('Displacement (mm)', fontsize=11)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f'wall_envelope_{tag}.png', dpi=150)
print(f'Saved: wall_envelope_{tag}.png')

# Statistics
if fd_new is not None:
    d, f = fd_new[:,0], fd_new[:,1]
    print(f'\n{"="*60}')
    print(f'  ANALYSIS SUMMARY — ADINA-{tag}')
    print(f'{"="*60}')
    print(f'  Points:      {len(fd_new)}')
    print(f'  Disp range:  [{d.min():.2f}, {d.max():.2f}] mm')
    print(f'  Force range: [{f.min():.1f}, {f.max():.1f}] kN')
    pos, neg = extract_envelope(fd_new)
    if pos is not None:
        print(f'  Positive peaks:')
        for p in pos: print(f'    d={p[0]:+7.2f} mm  F={p[1]:+8.1f} kN')
    if neg is not None:
        print(f'  Negative peaks:')
        for p in neg: print(f'    d={p[0]:+7.2f} mm  F={p[1]:+8.1f} kN')
if ref is not None:
    d, f = ref[:,0], ref[:,1]
    print(f'\n  Reference:')
    print(f'    Points:      {len(ref)}')
    print(f'    Disp range:  [{d.min():.2f}, {d.max():.2f}] mm')
    print(f'    Force range: [{f.min():.1f}, {f.max():.1f}] kN')
