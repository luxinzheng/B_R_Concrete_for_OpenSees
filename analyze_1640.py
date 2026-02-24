"""Analyze OpenSees-1640 wall results."""
import numpy as np
import sys
sys.stdout.reconfigure(line_buffering=True)

# Read displacement data
disp_data = []
with open('53.txt', 'r') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2:
            try:
                t = float(parts[0])
                d = float(parts[1])
                disp_data.append([t, d * 1000])  # m to mm
            except:
                pass

# Read reaction data
react_data = []
with open('1.txt', 'r') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 6:
            try:
                t = float(parts[0])
                # Sum all x-reactions (nodes 1-5, dof 1)
                base_shear = sum(float(x) for x in parts[1:6])
                react_data.append([t, base_shear / 1000])  # N to kN
            except:
                pass

disp = np.array(disp_data) if disp_data else np.zeros((0,2))
react = np.array(react_data) if react_data else np.zeros((0,2))

print(f'Displacement data: {len(disp)} points')
print(f'Reaction data:     {len(react)} points')

if len(disp) > 0:
    print(f'\nDisplacement range: [{disp[:,1].min():.4f}, {disp[:,1].max():.4f}] mm')
    print(f'Time range: [{disp[:,0].min():.4f}, {disp[:,0].max():.4f}]')

if len(react) > 0:
    # Find where forces become unreasonable (> 500 kN)
    reasonable = np.abs(react[:,1]) < 500
    n_reasonable = np.sum(reasonable)
    print(f'\nReasonable force points (|F| < 500 kN): {n_reasonable} / {len(react)}')
    
    if n_reasonable > 0:
        r_good = react[reasonable]
        print(f'Force range (good): [{r_good[:,1].min():.1f}, {r_good[:,1].max():.1f}] kN')
        print(f'Time range (good):  [{r_good[:,0].min():.4f}, {r_good[:,0].max():.4f}]')
    
    # Find divergence point
    for i in range(len(react)-1):
        if abs(react[i+1,1]) > 500 and abs(react[i,1]) < 500:
            print(f'\nDivergence at step {i}:')
            start = max(0, i-5)
            for j in range(start, min(i+10, len(react))):
                print(f'  t={react[j,0]:.5f}  F={react[j,1]:+14.1f} kN')
            break

# Merge by time to get force-displacement
if len(disp) > 0 and len(react) > 0:
    # Match times
    fd_pairs = []
    di, ri = 0, 0
    while di < len(disp) and ri < len(react):
        if abs(disp[di,0] - react[ri,0]) < 0.001:
            if abs(react[ri,1]) < 500:
                fd_pairs.append([disp[di,1], react[ri,1]])
            di += 1
            ri += 1
        elif disp[di,0] < react[ri,0]:
            di += 1
        else:
            ri += 1
    
    if fd_pairs:
        fd = np.array(fd_pairs)
        print(f'\nForce-displacement pairs: {len(fd)}')
        print(f'Displacement range: [{fd[:,0].min():.2f}, {fd[:,0].max():.2f}] mm')
        print(f'Force range: [{fd[:,1].min():.1f}, {fd[:,1].max():.1f}] kN')
        
        # Compare with reference
        ref_rows = []
        with open('参考模型.txt', 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        rd, rforce = float(parts[0]), float(parts[1])
                        if abs(rd) > 1e-6 or abs(rforce) > 100:
                            ref_rows.append([rd*1000, rforce/1000])
                    except:
                        pass
        
        ref = np.array(ref_rows) if ref_rows else None
        
        # Plot
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        ax = axes[0]
        if ref is not None:
            ax.plot(ref[:,0], ref[:,1], 'b-', linewidth=0.8, label='Reference', alpha=0.5)
        ax.plot(fd[:,0], fd[:,1], 'r-', linewidth=0.8, label='OpenSees-1640')
        ax.set_xlabel('Displacement (mm)')
        ax.set_ylabel('Base Shear (kN)')
        ax.set_title('Hysteresis Comparison')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Force vs time
        ax = axes[1]
        r_good = react[np.abs(react[:,1]) < 500]
        ax.plot(r_good[:,0], r_good[:,1], 'r-', linewidth=0.5)
        ax.set_xlabel('Time')
        ax.set_ylabel('Base Shear (kN)')
        ax.set_title('Force vs Time (before divergence)')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('wall_1640_comparison.png', dpi=150)
        print('\nSaved wall_1640_comparison.png')
        
        # Key metrics
        print('\nKey displacement peaks:')
        peaks = [0]
        for i in range(1, len(fd)-1):
            if (fd[i,0] > fd[i-1,0] and fd[i,0] > fd[i+1,0]) or \
               (fd[i,0] < fd[i-1,0] and fd[i,0] < fd[i+1,0]):
                peaks.append(i)
        for p in peaks[:20]:
            print(f'  d={fd[p,0]:+8.2f}mm  F={fd[p,1]:+8.1f}kN')

# Compare with previous versions
print('\n' + '='*60)
print('VERSION COMPARISON')
print('='*60)
# Check if we can compare with OpenSees-1205
for ver in ['1045', '1205', '1535']:
    rfile = f'wall_OpenSees_{ver}_react.txt'
    try:
        with open(rfile) as f:
            lines = f.readlines()
        valid = 0
        for line in lines:
            parts = line.strip().split()
            if len(parts) >= 6:
                try:
                    fs = sum(float(x) for x in parts[1:6]) / 1000
                    if abs(fs) < 500:
                        valid += 1
                except:
                    pass
        print(f'  OpenSees-{ver}: {valid} good force points')
    except FileNotFoundError:
        print(f'  OpenSees-{ver}: no data')

print(f'  OpenSees-1640: {n_reasonable} good force points')
if ref is not None:
    print(f'  Reference:     {len(ref)} points')
