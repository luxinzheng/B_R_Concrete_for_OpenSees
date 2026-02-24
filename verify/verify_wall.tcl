# ==============================================================================
# VERIFICATION 2: Shear Wall — cyclic lateral loading
#   Same geometry as csp3.tcl (1.0m wide x 2.0m tall, 0.125m thick)
#   55 nodes, 40 shell elements, truss rebar in confined regions
#   Gravity + cyclic displacement from shuju2.txt
# ==============================================================================
wipe
model basic -ndm 3 -ndf 6

if {![info exists model_type]} { set model_type "adina" }
set prefix "vwall_${model_type}"
source [file join [file dirname [info script]] verify_mat.tcl]

# Sections
section LayeredShell 1 10 \
    4 0.0125  11 0.0002403  11 0.0003676 \
    4 0.024696  4 0.024696  4 0.024696  4 0.024696 \
    11 0.0003676  11 0.0002403  4 0.0125
section LayeredShell 2 8 \
    4 0.0125  11 0.0002403  10 0.0002356 \
    4 0.0495241  4 0.0495241 \
    10 0.0002356  11 0.0002403  4 0.0125

source [file join [file dirname [info script]] .. node.tcl]
source [file join [file dirname [info script]] .. element.tcl]

# Truss rebar
foreach {eid n1 n2} {
    41 1 6  42 6 11  43 11 16  44 16 21  45 21 26
    46 26 31 47 31 36 48 36 41 49 41 46 50 46 51
    51 2 7  52 7 12  53 12 17  54 17 22  55 22 27
    56 27 32 57 32 37 58 37 42 59 42 47 60 47 52
    61 4 9  62 9 14  63 14 19  64 19 24  65 24 29
    66 29 34 67 34 39 68 39 44 69 44 49 70 49 54
    71 5 10 72 10 15 73 15 20 74 20 25 75 25 30
    76 30 35 77 35 40 78 40 45 79 45 50 80 50 55
} { element truss $eid $n1 $n2 223.53e-6 7 }

fixY 0.0 1 1 1 1 1 1

recorder Node -file ${prefix}_react.txt -time -node 1 2 3 4 5 -dof 1 reaction
recorder Node -file ${prefix}_disp.txt  -time -node 53 -dof 1 disp

# Gravity
pattern Plain 1 Linear { load 53 0 -246000 0 0 0 0 }

constraints Plain; numberer RCM; system BandGeneral
test NormDispIncr 1.0e-6 200; algorithm BFGS -count 100
integrator LoadControl 0.1; analysis Static
analyze 10
puts "gravity ok"
loadConst -time 0.0

# Cyclic loading
timeSeries Path 1 -dt 0.1 -filePath [file join [file dirname [info script]] .. shuju2.txt]
pattern Plain 2 1 { sp 53 1 1 }

constraints Penalty 1e20 1e20; numberer RCM; system BandGeneral
test NormDispIncr 1.0e-5 1000 2; algorithm KrylovNewton
integrator LoadControl 0.1; analysis Static

set fd [open "${prefix}_manual.txt" w]
proc logStep {fd} {
    set t [getTime]
    set d53 [nodeDisp 53 1]
    set r1 [nodeReaction 1 1]; set r2 [nodeReaction 2 1]
    set r3 [nodeReaction 3 1]; set r4 [nodeReaction 4 1]; set r5 [nodeReaction 5 1]
    puts $fd "$t $d53 $r1 $r2 $r3 $r4 $r5"
}

set totalTime 50.0
set dt 0.1
set step 0; set nFail 0; set nSkip 0; set consecFail 0

while {[getTime] < $totalTime} {
    set ok [analyze 1]; incr step
    if {$ok == 0} { set consecFail 0; logStep $fd
        if {[expr {$step % 50}] == 0} { flush $fd; puts "  step $step t=[getTime] f=$nFail s=$nSkip" }
        continue
    }
    algorithm ModifiedNewton; set ok [analyze 1]; algorithm KrylovNewton
    if {$ok == 0} { set consecFail 0; logStep $fd; continue }

    set nSub 10
    integrator LoadControl [expr {$dt/$nSub}]
    test NormDispIncr 1.0e-4 500 0; algorithm Newton
    set ok 0
    for {set i 0} {$i < $nSub} {incr i} { set ok [analyze 1]; if {$ok != 0} break }
    integrator LoadControl $dt; algorithm KrylovNewton; test NormDispIncr 1.0e-5 1000 2
    if {$ok == 0} { set consecFail 0; logStep $fd; continue }

    set nSub 100
    integrator LoadControl [expr {$dt/$nSub}]
    test NormDispIncr 1.0e-3 1000 0; algorithm Newton
    set ok 0
    for {set i 0} {$i < $nSub} {incr i} {
        set ok [analyze 1]
        if {$ok != 0} { algorithm ModifiedNewton; set ok [analyze 1]; algorithm Newton }
        if {$ok != 0} break
    }
    integrator LoadControl $dt; algorithm KrylovNewton; test NormDispIncr 1.0e-5 1000 2
    if {$ok == 0} { set consecFail 0; logStep $fd; continue }

    set nSub 1000
    integrator LoadControl [expr {$dt/$nSub}]
    algorithm ModifiedNewton -initial; test NormDispIncr 1.0e10 1 0
    for {set i 0} {$i < $nSub} {incr i} { set ok [analyze 1]; if {$ok != 0} break }
    integrator LoadControl $dt; algorithm KrylovNewton; test NormDispIncr 1.0e-5 1000 2
    if {$ok == 0} { incr nSkip; set consecFail 0; logStep $fd; continue }

    incr nFail; incr consecFail
    if {$consecFail >= 20} { puts "*** $consecFail consecutive fails — stopping ***"; break }
}
flush $fd; close $fd
puts "WALL ($model_type): t=[getTime] steps=$step fails=$nFail skips=$nSkip"
