"""Plot wall hysteresis comparison: ADINA-2145 vs Reference model."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

WORKDIR = r'E:\Basic\Concrete_Model\ADINA'

def load_opensees_wall(react_file, disp_file):
    """Load OpenSees wall results (reaction + displacement)."""
    react_data, disp_data = [], []
    if os.path.exists(react_file):
        with open(react_file) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 6:
                    try:
                        t = float(parts[0])
                        shear = sum(float(x) for x in parts[1:6]) / 1000  # N -> kN
                        react_data.append([t, shear])
                    except: pass

    if os.path.exists(disp_file):
        with open(disp_file) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        disp_data.append([float(parts[0]), float(parts[1]) * 1000])  # m -> mm
                    except: pass

    react = np.array(react_data) if react_data else np.zeros((0, 2))
    disp = np.array(disp_data) if disp_data else np.zeros((0, 2))

    # Match by time
    fd_pairs = []
    if len(disp) > 0 and len(react) > 0:
        di, ri = 0, 0
        while di < len(disp) and ri < len(react):
            if abs(disp[di, 0] - react[ri, 0]) < 0.001:
                if abs(react[ri, 1]) < 500:
                    fd_pairs.append([disp[di, 1], react[ri, 1]])
                di += 1; ri += 1
            elif disp[di, 0] < react[ri, 0]:
                di += 1
            else:
                ri += 1
    return np.array(fd_pairs) if fd_pairs else None

# Load reference model
ref_rows = []
try:
    with open(os.path.join(WORKDIR, '参考模型.txt')) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    ref_rows.append([float(parts[0]) * 1000, float(parts[1]) / 1000])
                except: pass
except: pass
ref = np.array(ref_rows) if ref_rows else None

# Load ADINA-2145
fd_2145 = load_opensees_wall(
    os.path.join(WORKDIR, 'wall_OpenSees_2145_react.txt'),
    os.path.join(WORKDIR, 'wall_OpenSees_2145_disp.txt'))

# Load ADINA-1045 (best previous)
fd_1045 = load_opensees_wall(
    os.path.join(WORKDIR, 'wall_OpenSees_1045_react.txt'),
    os.path.join(WORKDIR, 'wall_OpenSees_1045_disp.txt'))

print(f'Reference: {len(ref) if ref is not None else 0} points')
print(f'ADINA-2145: {len(fd_2145) if fd_2145 is not None else 0} points')
print(f'ADINA-1045: {len(fd_1045) if fd_1045 is not None else 0} points')

# ================================================================
# FIGURE 1: Main hysteresis comparison (2145 vs Reference)
# ================================================================
fig, axes = plt.subplots(1, 2, figsize=(16, 7))

ax = axes[0]
ax.set_title('Hysteresis: ADINA-2145 vs Reference', fontsize=13, fontweight='bold')
if ref is not None:
    ax.plot(ref[:, 0], ref[:, 1], 'b-', lw=0.7, alpha=0.6, label='Reference')
if fd_2145 is not None:
    ax.plot(fd_2145[:, 0], fd_2145[:, 1], 'r-', lw=0.9, label='ADINA-2145')
ax.set_xlabel('Displacement (mm)', fontsize=11)
ax.set_ylabel('Base Shear (kN)', fontsize=11)
ax.legend(fontsize=10, loc='upper left')
ax.grid(True, alpha=0.3)
ax.axhline(0, color='k', lw=0.5)
ax.axvline(0, color='k', lw=0.5)

# Zoom to ADINA range
ax = axes[1]
ax.set_title('Hysteresis Detail (ADINA range)', fontsize=13, fontweight='bold')
if ref is not None:
    ax.plot(ref[:, 0], ref[:, 1], 'b-', lw=0.7, alpha=0.6, label='Reference')
if fd_2145 is not None:
    ax.plot(fd_2145[:, 0], fd_2145[:, 1], 'r-', lw=0.9, label='ADINA-2145')
    xlim = max(abs(fd_2145[:, 0].min()), abs(fd_2145[:, 0].max())) * 1.15
    ax.set_xlim([-xlim, xlim])
ax.set_xlabel('Displacement (mm)', fontsize=11)
ax.set_ylabel('Base Shear (kN)', fontsize=11)
ax.legend(fontsize=10, loc='upper left')
ax.grid(True, alpha=0.3)
ax.axhline(0, color='k', lw=0.5)
ax.axvline(0, color='k', lw=0.5)

plt.tight_layout()
plt.savefig(os.path.join(WORKDIR, 'wall_hysteresis_2145.png'), dpi=150)
print('Saved: wall_hysteresis_2145.png')

# ================================================================
# FIGURE 2: Envelope + version comparison
# ================================================================
fig2, axes2 = plt.subplots(1, 2, figsize=(16, 7))

# Envelope comparison
ax = axes2[0]
ax.set_title('Envelope Comparison', fontsize=13, fontweight='bold')

def extract_envelope(fd):
    """Extract peak forces at each displacement amplitude."""
    if fd is None or len(fd) == 0:
        return None, None
    pos_env, neg_env = [], []
    # Group by displacement peaks
    d = fd[:, 0]
    f = fd[:, 1]
    # Find displacement peaks
    for i in range(1, len(d) - 1):
        if d[i] > d[i-1] and d[i] > d[i+1] and d[i] > 0.5:
            pos_env.append([d[i], f[i]])
        elif d[i] < d[i-1] and d[i] < d[i+1] and d[i] < -0.5:
            neg_env.append([d[i], f[i]])
    return (np.array(pos_env) if pos_env else None,
            np.array(neg_env) if neg_env else None)

for label, fd, color, marker in [
    ('Reference', ref, 'b', 'o'),
    ('ADINA-2145', fd_2145, 'r', 's'),
    ('ADINA-1045', fd_1045, 'gray', '^')]:
    if fd is not None:
        ax.plot(fd[:, 0], fd[:, 1], color=color, lw=0.3, alpha=0.2)
        pos, neg = extract_envelope(fd)
        if pos is not None:
            ax.plot(pos[:, 0], pos[:, 1], color=color, marker=marker, ms=5,
                    lw=1.5, label=f'{label} (pos peak)')
        if neg is not None:
            ax.plot(neg[:, 0], neg[:, 1], color=color, marker=marker, ms=5,
                    lw=1.5, ls='--', label=f'{label} (neg peak)')

ax.set_xlabel('Displacement (mm)', fontsize=11)
ax.set_ylabel('Base Shear (kN)', fontsize=11)
ax.legend(fontsize=8, loc='best')
ax.grid(True, alpha=0.3)
ax.axhline(0, color='k', lw=0.5)

# Force-time history
ax = axes2[1]
ax.set_title('Displacement History', fontsize=13, fontweight='bold')
if ref is not None:
    n = len(ref)
    ax.plot(np.arange(n), ref[:, 0], 'b-', lw=0.5, alpha=0.6, label='Reference')
if fd_2145 is not None:
    n = len(fd_2145)
    ax.plot(np.arange(n), fd_2145[:, 0], 'r-', lw=0.8, label='ADINA-2145')
if fd_1045 is not None:
    n = len(fd_1045)
    ax.plot(np.arange(n), fd_1045[:, 0], 'gray', lw=0.5, alpha=0.5, label='ADINA-1045')
ax.set_xlabel('Step', fontsize=11)
ax.set_ylabel('Displacement (mm)', fontsize=11)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(WORKDIR, 'wall_envelope_2145.png'), dpi=150)
print('Saved: wall_envelope_2145.png')

# ================================================================
# Print statistics
# ================================================================
print('\n' + '='*60)
print('  ANALYSIS SUMMARY')
print('='*60)

if fd_2145 is not None:
    d = fd_2145[:, 0]
    f = fd_2145[:, 1]
    print(f'  ADINA-2145:')
    print(f'    Points:     {len(fd_2145)}')
    print(f'    Disp range: [{d.min():.2f}, {d.max():.2f}] mm')
    print(f'    Force range:[{f.min():.1f}, {f.max():.1f}] kN')
    # Count cycles
    sign_changes = 0
    for i in range(1, len(d)):
        if d[i] * d[i-1] < 0:
            sign_changes += 1
    print(f'    Approx cycles: {sign_changes // 2}')

    pos, neg = extract_envelope(fd_2145)
    if pos is not None:
        print(f'    Positive peaks:')
        for p in pos:
            print(f'      d={p[0]:+7.2f} mm  F={p[1]:+8.1f} kN')
    if neg is not None:
        print(f'    Negative peaks:')
        for p in neg:
            print(f'      d={p[0]:+7.2f} mm  F={p[1]:+8.1f} kN')

if ref is not None:
    d = ref[:, 0]
    f = ref[:, 1]
    print(f'\n  Reference:')
    print(f'    Points:     {len(ref)}')
    print(f'    Disp range: [{d.min():.2f}, {d.max():.2f}] mm')
    print(f'    Force range:[{f.min():.1f}, {f.max():.1f}] kN')

print()
