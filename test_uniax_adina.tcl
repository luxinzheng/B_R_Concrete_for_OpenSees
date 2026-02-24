# ==============================================================================
# Single-Element Uniaxial Test -- ADINA Model (37-prop concrete)
# 
# 1m x 1m ShellMITC4, 3 layers x 0.1m = 0.3m thick
# Stress = total_reaction_y / (1.0m * 0.3m)
# Strain = uy_top / 1.0m
# Loading: compression -> unload -> tension -> reload
# ==============================================================================
wipe
model basic -ndm 3 -ndf 6

# --- Material: ADINA concrete (37 props) ---
set fc_abs   20.7e6
set ft       2.07e6
set E0       21.4e9
set nu       0.2
set fc       [expr {-$fc_abs}]
set epsc    -0.002
set fu       [expr {-0.2*$fc_abs}]
set epsu    -0.006

nDMaterial PlaneStressUserMaterial 1 13 37 \
    $E0  $nu  $ft  $fc  $epsc  $fu  $epsu \
    1.0e-4  0.5  0.75  1.0  0.7  0.12 \
    0.0  0.25  0.5  0.75  1.0  1.2 \
    1.0  1.4   1.7  2.2   2.5  2.8 \
    1.3  1.5   2.0  2.3   2.7  3.2 \
    1.25 1.45  1.95 2.25  2.65 3.15

set G_out [expr {$E0 / (2.0*(1.0+$nu))}]
nDMaterial PlateFromPlaneStress 2 1 $G_out

# --- Section: 3 layers, total 0.3m ---
section LayeredShell 1 3 2 0.1 2 0.1 2 0.1

# --- Nodes ---
node 1  0.0  0.0  0.0
node 2  1.0  0.0  0.0
node 3  1.0  1.0  0.0
node 4  0.0  1.0  0.0

# --- Element ---
element ShellMITC4 1  1 2 3 4  1

# --- Boundary conditions ---
# True uniaxial stress: all nodes free in ux (Poisson expansion)
fix 1  1 1 1 1 1 1
fix 2  0 1 1 1 1 1
fix 3  0 0 1 1 1 1
fix 4  0 0 1 1 1 1
equalDOF 3 4 1 2

# --- Loading ---
timeSeries Linear 1
pattern Plain 1 1 {
    sp 3 2 1.0
    sp 4 2 1.0
}

# --- Recorders ---
recorder Node -file uniax_react_adina.txt -time -node 1 2 -dof 2 reaction
recorder Node -file uniax_disp_adina.txt  -time -node 3    -dof 2 disp

# --- Analysis ---
constraints Penalty 1.0e20 1.0e20
numberer RCM
system BandGeneral
test NormDispIncr 1.0e-8 200 0
algorithm Newton

# ==============================================================
# Phase 1: Compression  0 -> -0.006 (600 steps)
# ==============================================================
puts "Phase 1: Compression to -0.006..."
set dL 1.0e-5
integrator LoadControl [expr {-$dL}]
analysis Static
set ok 0
for {set i 0} {$i < 600} {incr i} {
    set ok [analyze 1]
    if {$ok != 0} {
        puts "  Failed at step $i, disp=[nodeDisp 3 2]"
        break
    }
}
puts "  Phase 1 done: disp=[nodeDisp 3 2]"

# ==============================================================
# Phase 2: Unload  -0.006 -> 0 (600 steps)
# ==============================================================
puts "Phase 2: Unload to 0..."
integrator LoadControl [expr {$dL}]
analysis Static
for {set i 0} {$i < 600} {incr i} {
    set ok [analyze 1]
    if {$ok != 0} {
        puts "  Failed at step $i, disp=[nodeDisp 3 2]"
        break
    }
}
puts "  Phase 2 done: disp=[nodeDisp 3 2]"

# ==============================================================
# Phase 3: Tension  0 -> +0.0005 (50 steps)
# ==============================================================
puts "Phase 3: Tension to +0.0005..."
integrator LoadControl [expr {$dL}]
analysis Static
for {set i 0} {$i < 50} {incr i} {
    set ok [analyze 1]
    if {$ok != 0} {
        puts "  Failed at step $i, disp=[nodeDisp 3 2]"
        break
    }
}
puts "  Phase 3 done: disp=[nodeDisp 3 2]"

# ==============================================================
# Phase 4: Reload compression  +0.0005 -> -0.003 (350 steps)
# ==============================================================
puts "Phase 4: Reload compression to -0.003..."
integrator LoadControl [expr {-$dL}]
analysis Static
for {set i 0} {$i < 350} {incr i} {
    set ok [analyze 1]
    if {$ok != 0} {
        puts "  Failed at step $i, disp=[nodeDisp 3 2]"
        break
    }
}
puts "  Phase 4 done: disp=[nodeDisp 3 2]"

puts "=== Test complete ==="
wipe
