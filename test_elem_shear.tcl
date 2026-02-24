# ==============================================================================
# Single-element PURE SHEAR test
# 1m x 1m ShellMITC4, 0.1m thick pure concrete
# Monotonic shear to gamma_xy = 0.005
# ==============================================================================
wipe
model basic -ndm 3 -ndf 6

# ---- Material ----
set E0    21.4e9
set nu    0.2
set ft    2.07e6
set fc    -20.7e6
set epsc  -0.002
set fu    -4.14e6
set epsu  -0.006

nDMaterial PlaneStressUserMaterial 1 13 37 \
    $E0  $nu  $ft  $fc  $epsc  $fu  $epsu \
    1.0e-4  0.5  0.75  1.0  0.7  0.12 \
    0.0  0.25  0.5  0.75  1.0  1.2 \
    1.0  1.4   1.7  2.2   2.5  2.8 \
    1.3  1.5   2.0  2.3   2.7  3.2 \
    1.25 1.45  1.95 2.25  2.65 3.15

set G_out [expr {$E0 / (2.0*(1.0+$nu))}]
nDMaterial PlateFromPlaneStress 2 1 $G_out
section LayeredShell 1 4  2 0.025  2 0.025  2 0.025  2 0.025

# ---- Nodes ----
node 1 0.0 0.0 0.0
node 2 1.0 0.0 0.0
node 3 1.0 1.0 0.0
node 4 0.0 1.0 0.0
element ShellMITC4 1  1 2 3 4  1

# ---- BCs for pure shear (eps_xx=0, eps_yy=0, gamma_xy!=0) ----
# Bottom: fix ux, uy; Top: prescribe ux, fix uy
fix 1  1 1 1 1 1 1
fix 2  1 1 1 1 1 1
fix 3  0 1 1 1 1 1
fix 4  0 1 1 1 1 1

# ---- Recorders ----
recorder Node -file shear_disp.txt  -time -node 3 4 -dof 1 2 disp
recorder Node -file shear_react.txt -time -node 1 2 -dof 1 2 reaction

# ---- Loading: shear ux = 0.005m at top (gamma = 0.005) ----
set target 0.005
set nSteps 100
set dt [expr {1.0/$nSteps}]

timeSeries Linear 1
pattern Plain 1 1 {
    sp 3 1 $target
    sp 4 1 $target
}

# ---- Analysis ----
constraints Penalty 1e20 1e20
numberer Plain
system FullGeneral
test NormDispIncr 1.0e-8 50 2
algorithm Newton
integrator LoadControl $dt
analysis Static

set nFail 0
for {set i 1} {$i <= $nSteps} {incr i} {
    set ok [analyze 1]
    if {$ok != 0} {
        algorithm ModifiedNewton
        set ok [analyze 1]
        algorithm Newton
    }
    if {$ok != 0} {
        integrator LoadControl [expr {$dt/10.0}]
        test NormDispIncr 1.0e-6 200 0
        set ok 0
        for {set j 0} {$j < 10} {incr j} {
            set ok [analyze 1]
            if {$ok != 0} break
        }
        integrator LoadControl $dt
        test NormDispIncr 1.0e-8 50 2
        algorithm Newton
    }
    if {$ok != 0} {
        incr nFail
        puts "FAIL step $i at time [getTime]"
    }
}
puts "SHEAR_TEST: time=[getTime] steps=$nSteps fails=$nFail"
remove recorders
