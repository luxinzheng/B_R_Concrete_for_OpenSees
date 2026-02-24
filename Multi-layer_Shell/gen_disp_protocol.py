# -*- coding: utf-8 -*-
"""Generate displacement protocol file for Multi-layer_Shell sp+LoadControl approach"""

# Original DC protocol: {dincr(mm) nsteps} in groups of 3 (push, pull-back, return)
blocks = [
    (1, 2),  (-1, 4),   (1, 2),     # ±2mm
    (1, 4),  (-1, 8),   (1, 4),     # ±4mm
    (1, 5),  (-1, 10),  (1, 5),     # ±5mm
    (1, 9),  (-1, 18),  (1, 9),     # ±9mm
    (1, 13), (-1, 25),  (1, 12),    # ±13mm
    (1, 18), (-1, 36),  (1, 18),    # ±18mm
    (1, 24), (-1, 48),  (1, 24),    # ±24mm
    (1, 30), (-1, 60),  (1, 30),    # ±30mm
    (1, 36), (-1, 72),  (1, 36),    # ±36mm
    (1, 45), (-1, 90),  (1, 45),    # ±45mm
    (1, 52), (-1, 102), (1, 50),    # ±52mm
]

disp = 0.0
disp_list = [0.0]  # start at 0

for dincr_mm, nsteps in blocks:
    for _ in range(nsteps):
        disp += dincr_mm * 1e-3  # convert mm to m... wait, units?
        disp_list.append(disp)

# Check: the original protocol uses 1e-3 as dincr for DisplacementControl.
# In the OpenSees model, the node coordinates are in meters (node.tcl uses 0.2, 0.4, etc.)
# So 1e-3 = 1mm displacement in m units
# The sp command imposes: displacement = timeSeries_value * sp_factor
# We use sp_factor = 1, so timeSeries values = actual displacement in model units (m)

# Actually, let me re-check. In sw1-1, sp 53 1 1 means impose disp at node 53, dof 1, scale=1.
# The timeSeries Path values in shuju2.txt are the actual displacement values.
# Let me check what unit system the ML model uses.

# Looking at node.tcl for ML: coordinates are like 0.0, 0.2, 0.4, ... 3.6 (meters)
# The original dincr = 1e-3 for DisplacementControl = 0.001 m = 1 mm
# So the displacement values should be in meters

# Write displacement values (in meters, since the model uses m)
with open('ml_disp_protocol.txt', 'w') as f:
    for d in disp_list:
        f.write(f'{d:.6f}\n')

print(f'Total steps: {len(disp_list) - 1}')
print(f'Disp range: [{min(disp_list)*1000:.1f}, {max(disp_list)*1000:.1f}] mm')
print(f'Final disp: {disp_list[-1]*1000:.1f} mm')
print(f'File: ml_disp_protocol.txt ({len(disp_list)} values)')

# Also show the protocol summary
print('\nProtocol summary:')
pos = 0.0
for i, (d, n) in enumerate(blocks):
    pos += d * n * 1e-3
    print(f'  Block {i+1:2d}: dincr={d:+d}mm × {n:3d} steps → position = {pos*1000:+.0f} mm')
