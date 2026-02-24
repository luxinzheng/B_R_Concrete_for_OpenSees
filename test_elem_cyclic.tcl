# ==============================================================================
# Single-element CYCLIC COMPRESSION-TENSION test
# 1m x 1m ShellMITC4, 0.1m thick pure concrete
# Cyclic path that exercises critical state transitions:
#   compress -> release -> tension (crack) -> release -> re-compress -> ...
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

# ---- BCs ----
fix 1  1 1 1 1 1 1
fix 2  0 1 1 1 1 1
fix 3  0 0 1 1 1 1
fix 4  1 0 1 1 1 1

# ---- Recorders ----
recorder Node -file cyclic_disp.txt  -time -node 3 4 -dof 1 2 disp
recorder Node -file cyclic_react.txt -time -node 1 2 -dof 1 2 reaction

# ---- Cyclic displacement path (m) ----
# Phase 1: 0 -> -0.001 (compress eps=-0.001, pre-peak)
# Phase 2: -0.001 -> 0 (release)
# Phase 3: 0 -> +0.00015 (tension, should crack at ~1e-4)
# Phase 4: +0.00015 -> 0 (release from cracked state)
# Phase 5: 0 -> -0.002 (re-compress to peak, through cracked material)
# Phase 6: -0.002 -> 0 (release)
# Phase 7: 0 -> +0.0003 (tension, deeper cracking)
# Phase 8: +0.0003 -> 0 (release)
# Phase 9: 0 -> -0.003 (compress past peak, softening)
# Phase 10: -0.003 -> 0 (release from softened state)

timeSeries Path 1 -dt 1.0 -values {
    0.0
    -0.001
    0.0
    0.00015
    0.0
    -0.002
    0.0
    0.0003
    0.0
    -0.003
    0.0
}

pattern Plain 1 1 {
    sp 3 2 1.0
    sp 4 2 1.0
}

# ---- Analysis: 50 sub-steps per phase, total 500 steps ----
set nStepsPerPhase 50
set nPhases 10
set totalSteps [expr {$nStepsPerPhase * $nPhases}]
set dt [expr {1.0/$nStepsPerPhase}]

constraints Penalty 1e20 1e20
numberer Plain
system FullGeneral
test NormDispIncr 1.0e-8 50 2
algorithm Newton
integrator LoadControl $dt
analysis Static

set nFail 0
set consecFail 0
for {set i 1} {$i <= $totalSteps} {incr i} {
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
        # 100 sub-steps with initial stiffness
        integrator LoadControl [expr {$dt/100.0}]
        algorithm ModifiedNewton -initial
        test NormDispIncr 1.0e-4 500 0
        set ok 0
        for {set j 0} {$j < 100} {incr j} {
            set ok [analyze 1]
            if {$ok != 0} break
        }
        integrator LoadControl $dt
        test NormDispIncr 1.0e-8 50 2
        algorithm Newton
    }

    if {$ok != 0} {
        incr nFail
        incr consecFail
        set phase [expr {int(($i-1)/$nStepsPerPhase) + 1}]
        puts "FAIL step $i phase=$phase at time [getTime] (consec=$consecFail)"
        if {$consecFail >= 10} {
            puts "*** 10 consecutive failures - stopping ***"
            break
        }
    } else {
        set consecFail 0
    }

    if {[expr {$i % $nStepsPerPhase}] == 0} {
        set phase [expr {$i / $nStepsPerPhase}]
        puts "  Phase $phase complete: time=[getTime], fails=$nFail"
    }
}
puts "CYCLIC_TEST: time=[getTime] steps=$i/$totalSteps fails=$nFail"
remove recorders
