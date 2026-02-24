# ==============================================================================
# VERIFICATION 3: RC Cantilever Beam — in-plane bending (XY plane)
#   Cantilever beam: L=2.0m, h=0.5m, b=0.125m (same thickness as wall)
#   8x4 mesh = 32 ShellMITC4
#   Fixed at left end (x=0), cyclic tip load at right end (x=L)
#   Uses same section/material as wall confined region
# ==============================================================================
wipe
model basic -ndm 3 -ndf 6

if {![info exists model_type]} { set model_type "adina" }
set prefix "vbeam_${model_type}"
source [file join [file dirname [info script]] verify_mat.tcl]

# Same confined section as wall
section LayeredShell 1 10 \
    4 0.0125  11 0.0002403  11 0.0003676 \
    4 0.024696  4 0.024696  4 0.024696  4 0.024696 \
    11 0.0003676  11 0.0002403  4 0.0125

set L 2.0; set H 0.5; set nx 8; set ny 4
set dx [expr {$L / $nx}]; set dy [expr {$H / $ny}]

for {set j 0} {$j <= $ny} {incr j} {
    for {set i 0} {$i <= $nx} {incr i} {
        set nid [expr {$j*($nx+1) + $i + 1}]
        node $nid [expr {$i*$dx}] [expr {$j*$dy}] 0.0
    }
}

set eid 1
for {set j 0} {$j < $ny} {incr j} {
    for {set i 0} {$i < $nx} {incr i} {
        set n1 [expr {$j*($nx+1) + $i + 1}]
        set n2 [expr {$n1 + 1}]
        set n3 [expr {$n2 + ($nx+1)}]
        set n4 [expr {$n1 + ($nx+1)}]
        element ShellMITC4 $eid $n1 $n2 $n3 $n4 1
        incr eid
    }
}

# Fix left end (x=0): all DOFs
for {set j 0} {$j <= $ny} {incr j} {
    set nid [expr {$j*($nx+1) + 1}]
    fix $nid 1 1 1 1 1 1
}

# Tip load node (top-right corner)
set tipNode [expr {$ny*($nx+1) + $nx + 1}]

recorder Node -file ${prefix}_disp.txt  -time -node $tipNode -dof 1 disp
# Reactions at left edge
set reactNodes [list]
for {set j 0} {$j <= $ny} {incr j} {
    lappend reactNodes [expr {$j*($nx+1) + 1}]
}
recorder Node -file ${prefix}_react.txt -time -node {*}$reactNodes -dof 1 reaction

# Cyclic horizontal load at tip (force-controlled)
set cyc_data [list 0.0]
foreach peak {5000 -5000 10000 -10000 20000 -20000 30000 -30000 40000 -40000 0} {
    set prev [lindex $cyc_data end]
    for {set i 1} {$i <= 10} {incr i} {
        lappend cyc_data [expr {$prev + ($peak - $prev) * $i / 10.0}]
    }
}
set nPts [llength $cyc_data]

set fd [open "${prefix}_protocol.txt" w]
foreach v $cyc_data { puts $fd $v }
close $fd

timeSeries Path 1 -dt 1.0 -filePath "${prefix}_protocol.txt"
pattern Plain 1 1 { load $tipNode 1.0 0 0 0 0 0 }

constraints Penalty 1e20 1e20; numberer RCM; system BandGeneral
test NormDispIncr 1.0e-6 200 0; algorithm KrylovNewton
integrator LoadControl 1.0; analysis Static

set nFail 0; set nSkip 0
set fd_log [open "${prefix}_manual.txt" w]
for {set i 1} {$i < $nPts} {incr i} {
    set ok [analyze 1]
    if {$ok != 0} {
        algorithm ModifiedNewton; set ok [analyze 1]; algorithm KrylovNewton
    }
    if {$ok != 0} {
        integrator LoadControl 0.1
        test NormDispIncr 1.0e-4 500 0; algorithm Newton
        set ok 0
        for {set j 0} {$j < 10} {incr j} { set ok [analyze 1]; if {$ok != 0} break }
        integrator LoadControl 1.0; test NormDispIncr 1.0e-6 200 0; algorithm KrylovNewton
    }
    if {$ok != 0} {
        integrator LoadControl 0.01
        algorithm ModifiedNewton -initial; test NormDispIncr 1.0e10 1 0
        for {set j 0} {$j < 100} {incr j} { analyze 1 }
        integrator LoadControl 1.0; algorithm KrylovNewton; test NormDispIncr 1.0e-6 200 0
        incr nSkip
    }
    set d [nodeDisp $tipNode 1]
    set rsum 0.0
    foreach rn $reactNodes { set rsum [expr {$rsum + [nodeReaction $rn 1]}] }
    puts $fd_log "[getTime] $d $rsum"
}
close $fd_log
puts "BEAM ($model_type): steps=[expr {$nPts-1}] fails=$nFail skips=$nSkip"
remove recorders
