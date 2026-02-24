# ==============================================================================
# VERIFICATION 1: Uniaxial Tests
#   sub_test = "comp" | "tens" | "cyc"   (set before sourcing)
#
# Single 1m x 1m x 0.1m ShellMITC4, pure concrete (4 equal layers)
# Uniaxial stress state: sigma_xx free, sigma_yy loaded
# ==============================================================================
wipe
model basic -ndm 3 -ndf 6

if {![info exists model_type]} { set model_type "adina" }
if {![info exists sub_test]}   { set sub_test "comp" }
set prefix "vuniax_${model_type}"
source [file join [file dirname [info script]] verify_mat.tcl]

set G_ps [expr {$E0 / (2.0*(1.0+$nu))}]
nDMaterial PlateFromPlaneStress 99 1 $G_ps
section LayeredShell 99 4  99 0.025  99 0.025  99 0.025  99 0.025

node 1 0.0 0.0 0.0; node 2 1.0 0.0 0.0
node 3 1.0 1.0 0.0; node 4 0.0 1.0 0.0
element ShellMITC4 1  1 2 3 4  99

fix 1  1 1 1 1 1 1; fix 2  0 1 1 1 1 1
fix 3  0 0 1 1 1 1; fix 4  1 0 1 1 1 1

recorder Node -file ${prefix}_${sub_test}_disp.txt  -time -node 3 4 -dof 2 disp
recorder Node -file ${prefix}_${sub_test}_react.txt -time -node 1 2 -dof 2 reaction

proc adaptive_analyze {nSteps dt} {
    set nFail 0
    for {set i 1} {$i <= $nSteps} {incr i} {
        set ok [analyze 1]
        if {$ok != 0} {
            algorithm ModifiedNewton; set ok [analyze 1]; algorithm Newton
        }
        if {$ok != 0} {
            integrator LoadControl [expr {$dt/10.0}]
            test NormDispIncr 1.0e-6 200 0
            for {set j 0} {$j < 10} {incr j} {
                set ok [analyze 1]; if {$ok != 0} break
            }
            integrator LoadControl $dt
            test NormDispIncr 1.0e-8 100 0; algorithm Newton
        }
        if {$ok != 0} { incr nFail; puts "  FAIL step $i" }
    }
    return $nFail
}

if {$sub_test eq "comp"} {
    # ── Monotonic compression to eps = -0.004 ──
    set nSteps 80; set dt [expr {1.0/$nSteps}]
    timeSeries Linear 1
    pattern Plain 1 1 { sp 3 2 -0.004; sp 4 2 -0.004 }
    constraints Penalty 1e20 1e20; numberer Plain; system FullGeneral
    test NormDispIncr 1.0e-8 100 0; algorithm Newton
    integrator LoadControl $dt; analysis Static
    set nf [adaptive_analyze $nSteps $dt]
    puts "COMP ($model_type): time=[getTime] fails=$nf"

} elseif {$sub_test eq "tens"} {
    # ── Monotonic tension to eps = +0.002 ──
    set nSteps 40; set dt [expr {1.0/$nSteps}]
    timeSeries Linear 1
    pattern Plain 1 1 { sp 3 2 0.002; sp 4 2 0.002 }
    constraints Penalty 1e20 1e20; numberer Plain; system FullGeneral
    test NormDispIncr 1.0e-8 100 0; algorithm Newton
    integrator LoadControl $dt; analysis Static
    set nf [adaptive_analyze $nSteps $dt]
    puts "TENS ($model_type): time=[getTime] fails=$nf"

} elseif {$sub_test eq "cyc"} {
    # ── Cyclic compression: 0 -> -0.002 -> 0 -> -0.003 -> 0 -> -0.004 -> 0 ──
    set cyc_data [list 0.0]
    foreach peak {-0.002 -0.003 -0.004} {
        set nHalf 30
        for {set i 1} {$i <= $nHalf} {incr i} {
            lappend cyc_data [expr {$peak * $i / double($nHalf)}]
        }
        for {set i [expr {$nHalf-1}]} {$i >= 0} {incr i -1} {
            lappend cyc_data [expr {$peak * $i / double($nHalf)}]
        }
    }
    set nPts [llength $cyc_data]
    set fd [open "${prefix}_cyc_protocol.txt" w]
    foreach v $cyc_data { puts $fd $v }
    close $fd

    timeSeries Path 1 -dt 1.0 -filePath "${prefix}_cyc_protocol.txt"
    pattern Plain 1 1 { sp 3 2 1.0; sp 4 2 1.0 }
    constraints Penalty 1e20 1e20; numberer Plain; system FullGeneral
    test NormDispIncr 1.0e-8 200 0; algorithm KrylovNewton
    integrator LoadControl 1.0; analysis Static

    set nFail 0
    for {set i 1} {$i < $nPts} {incr i} {
        set ok [analyze 1]
        if {$ok != 0} {
            algorithm ModifiedNewton; set ok [analyze 1]; algorithm KrylovNewton
        }
        if {$ok != 0} {
            integrator LoadControl 0.1
            test NormDispIncr 1.0e-6 500 0; algorithm Newton
            set ok 0
            for {set j 0} {$j < 10} {incr j} { set ok [analyze 1]; if {$ok != 0} break }
            integrator LoadControl 1.0; test NormDispIncr 1.0e-8 200 0; algorithm KrylovNewton
        }
        if {$ok != 0} { incr nFail; puts "  FAIL cyc step $i" }
    }
    puts "CYC ($model_type): time=[getTime] fails=$nFail"
}

remove recorders
puts "DONE: $sub_test $model_type"
