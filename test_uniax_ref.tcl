# ==============================================================================
# Single-Element Uniaxial Test -- REFERENCE Model (built-in 7-prop concrete)
# 
# 1m x 1m ShellMITC4, 3 layers x 0.1m = 0.3m thick
# Stress = total_reaction_y / (1.0m * 0.3m)
# Strain = uy_top / 1.0m
# Loading: compression -> unload -> tension -> reload
# ==============================================================================
wipe
model basic -ndm 3 -ndf 6

# --- Material: reference concrete (7 props) ---
#                                                fc       ft      fcu    epsc0   epscu  epstu   stc
nDMaterial PlaneStressUserMaterial 1 40 7    20.7e6   2.07e6  -4.14e6  -0.002  -0.006  0.001   0.08

nDMaterial PlateFromPlaneStress 2 1 1.25e10

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
# Node 1: pin (fix ux,uy to prevent rigid body motion)
# Node 2: roller (fix uy only)
# Nodes 3,4: free ux, prescribed uy via sp
fix 1  1 1 1 1 1 1
fix 2  0 1 1 1 1 1
fix 3  0 0 1 1 1 1
fix 4  0 0 1 1 1 1
equalDOF 3 4 1 2

# --- Loading: prescribed displacement via sp + LoadControl ---
timeSeries Linear 1
pattern Plain 1 1 {
    sp 3 2 1.0
    sp 4 2 1.0
}

# --- Recorders ---
recorder Node -file uniax_react_ref.txt -time -node 1 2 -dof 2 reaction
recorder Node -file uniax_disp_ref.txt  -time -node 3    -dof 2 disp

# --- Analysis ---
constraints Penalty 1.0e20 1.0e20
numberer RCM
system BandGeneral
test NormDispIncr 1.0e-8 200 0
algorithm Newton

# ==============================================================
# Phase 1: Compression  0 -> -0.006 (600 steps, dLambda = -1e-5)
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
