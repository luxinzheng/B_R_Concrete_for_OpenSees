wipe
model basic -ndm 3 -ndf 6

# B&R material (37 params) for sw2-1: fc=30.8MPa
nDMaterial PlaneStressUserMaterial    1   40   37   2.62524e+10  0.2  3.08e+06  -3.08e+07  -0.002  -6.16e+06  -0.005  0.0001  0.5  0.75  1  0.7  0.12  0  0.25  0.5  0.75  1  1.2  1  1.4  1.7  2.2  2.5  2.8  1.3  1.5  2  2.3  2.7  3.2  1.25  1.45  1.95  2.25  2.65  3.15
nDMaterial   PlateFromPlaneStress    4        1   1.283e10

uniaxialMaterial   Steel02        7 379e6   202.7e9  0.01 18.5 0.925 0.15
uniaxialMaterial   Steel02        8 392e6   200.6e9  0.01 18.5 0.925 0.15
nDMaterial   PlateRebar          9               7     90
nDMaterial   PlateRebar         10               8     90
nDMaterial   PlateRebar         11               8     0

section   LayeredShell       1        10        4      0.0125 11 0.0003023  11  0.0004367  4 0.0246305 4 0.0246305 4 0.0246305 4 0.0246305  11  0.0004367 11 0.0003023 4 0.0125
section   LayeredShell 2 8 4 0.0125 11 0.0003023  10  0.0002356 4 0.0494621 4 0.0494621  10  0.0002356 11 0.0003023 4 0.0125

source node.tcl
source element.tcl

element truss 41 1 6 223.53e-6 7
element truss 42 6 11 223.53e-6 7
element truss 43 11 16 223.53e-6 7
element truss 44 16 21 223.53e-6 7
element truss 45 21 26 223.53e-6 7

element truss 51  2  7 223.53e-6 7
element truss 52  7 12 223.53e-6 7
element truss 53 12 17 223.53e-6 7
element truss 54 17 22 223.53e-6 7
element truss 55 22 27 223.53e-6 7

element truss 61  4  9 223.53e-6 7
element truss 62  9 14 223.53e-6 7
element truss 63 14 19 223.53e-6 7
element truss 64 19 24 223.53e-6 7
element truss 65 24 29 223.53e-6 7

element truss 71 5  10 223.53e-6 7
element truss 72 10 15 223.53e-6 7
element truss 73 15 20 223.53e-6 7
element truss 74 20 25 223.53e-6 7
element truss 75 25 30 223.53e-6 7

fixY 0.0 1 1 1 1 1 1

recorder Node -file test_disp.txt -time -node  28 -dof 1  disp
recorder Node -file test_react.txt -time -node 1 2 3 4 5 -dof 1 reaction

pattern Plain 1 Linear {
load 27 0 -493000 0 0 0 0
load 29 0 -493000 0 0 0 0
}

constraints Plain
numberer RCM
system BandGeneral
test NormDispIncr 1.0e-5 500
algorithm KrylovNewton
integrator LoadControl 0.1
analysis Static

puts "Starting gravity..."
set gok [analyze 10]
record
if {$gok != 0} { puts "Gravity FAILED" } else { puts "Gravity OK" }
loadConst -time 0.0

timeSeries Path 1 -dt 0.1 -filePath shuju2.txt
pattern Plain 2 1 { sp 28 1 1 }

constraints Penalty 1e20 1e20
numberer RCM
system BandGeneral
test NormDispIncr 1.0e-5 1000 2
algorithm KrylovNewton
integrator LoadControl 0.1
analysis Static

set step 0; set nFail 0; set consecFail 0
while {[getTime] < 50.0} {
    set ok [analyze 1]; incr step
    if {$ok == 0} {
        set consecFail 0
        if {[expr {$step % 50}] == 0} { record; puts "  step $step t=[getTime] f=$nFail" }
        continue
    }
    algorithm ModifiedNewton; set ok [analyze 1]; algorithm KrylovNewton
    if {$ok == 0} { set consecFail 0; continue }
    set nSub 10
    integrator LoadControl [expr {0.1/$nSub}]
    test NormDispIncr 1.0e-4 500 0; algorithm Newton
    set ok 0
    for {set i 0} {$i < $nSub} {incr i} { set ok [analyze 1]; if {$ok != 0} break }
    integrator LoadControl 0.1; algorithm KrylovNewton; test NormDispIncr 1.0e-5 1000 2
    if {$ok == 0} { set consecFail 0; continue }
    set nSub 100
    integrator LoadControl [expr {0.1/$nSub}]
    test NormDispIncr 1.0e-3 1000 0; algorithm Newton
    set ok 0
    for {set i 0} {$i < $nSub} {incr i} {
        set ok [analyze 1]
        if {$ok != 0} { algorithm ModifiedNewton; set ok [analyze 1]; algorithm Newton }
        if {$ok != 0} break
    }
    integrator LoadControl 0.1; algorithm KrylovNewton; test NormDispIncr 1.0e-5 1000 2
    if {$ok == 0} { set consecFail 0; continue }
    set nSub 1000
    integrator LoadControl [expr {0.1/$nSub}]
    algorithm ModifiedNewton -initial; test NormDispIncr 1.0e10 1 0
    for {set i 0} {$i < $nSub} {incr i} { set ok [analyze 1]; if {$ok != 0} break }
    integrator LoadControl 0.1; algorithm KrylovNewton; test NormDispIncr 1.0e-5 1000 2
    if {$ok == 0} { incr nFail; set consecFail 0; continue }
    incr nFail; incr consecFail
    if {$consecFail >= 20} { puts "*** $consecFail consecutive fails -- stopping ***"; break }
}
record
puts "Done: steps=$step fails=$nFail"
