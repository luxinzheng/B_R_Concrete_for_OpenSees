# ==============================================================================
# VERIFICATION 4: RC Plate under cyclic out-of-plane point load
#   Square plate 1.0m x 1.0m, thickness 0.12m
#   4x4 mesh = 16 ShellMITC4 elements (25 nodes)
#   Simply supported on all 4 edges (fix uz on edges)
#   Point load at center node (force-controlled, cyclic)
#   LayeredShell: concrete + bidir rebar
# ==============================================================================
wipe
model basic -ndm 3 -ndf 6

if {![info exists model_type]} { set model_type "adina" }
set prefix "vplate_${model_type}"
source [file join [file dirname [info script]] verify_mat.tcl]

# Plate rebar (d=8mm @ 200mm spacing in both directions)
uniaxialMaterial Steel02 30 400e6 200e9 0.01 18.5 0.925 0.15
nDMaterial PlateRebar 31 30  0
nDMaterial PlateRebar 32 30 90

# Plate section: 120mm thick
set t_cover 0.020
set t_rebar [expr {50.27e-6 / 0.25}]
set t_core  [expr {(0.12 - 2.0*$t_cover - 4.0*$t_rebar) / 2.0}]
section LayeredShell 20 8 \
    4 $t_cover  31 $t_rebar  32 $t_rebar \
    4 $t_core  4 $t_core \
    32 $t_rebar  31 $t_rebar  4 $t_cover

# 5x5 node grid
set nDiv 4; set sideLen 1.0; set dx [expr {$sideLen / $nDiv}]
for {set j 0} {$j <= $nDiv} {incr j} {
    for {set i 0} {$i <= $nDiv} {incr i} {
        set nid [expr {$j*($nDiv+1) + $i + 1}]
        node $nid [expr {$i*$dx}] [expr {$j*$dx}] 0.0
    }
}

set eid 1
for {set j 0} {$j < $nDiv} {incr j} {
    for {set i 0} {$i < $nDiv} {incr i} {
        set n1 [expr {$j*($nDiv+1) + $i + 1}]
        set n2 [expr {$n1 + 1}]
        set n3 [expr {$n2 + ($nDiv+1)}]
        set n4 [expr {$n1 + ($nDiv+1)}]
        element ShellMITC4 $eid $n1 $n2 $n3 $n4 20
        incr eid
    }
}

# BCs: edges simply supported (fix uz), fix in-plane at corners
for {set j 0} {$j <= $nDiv} {incr j} {
    for {set i 0} {$i <= $nDiv} {incr i} {
        set nid [expr {$j*($nDiv+1) + $i + 1}]
        set onEdge [expr {$i==0 || $i==$nDiv || $j==0 || $j==$nDiv}]
        if {$i==0 && $j==0} {
            fix $nid 1 1 1 0 0 1
        } elseif {$i==$nDiv && $j==0} {
            fix $nid 0 1 1 0 0 1
        } elseif {$onEdge} {
            fix $nid 0 0 1 0 0 0
        }
    }
}

set centerNode [expr {($nDiv/2)*($nDiv+1) + ($nDiv/2) + 1}]

recorder Node -file ${prefix}_disp.txt  -time -node $centerNode -dof 3 disp

# Force-controlled cyclic loading
set cyc_loads [list 0.0]
foreach peak {5000 -5000 10000 -10000 15000 -15000 20000 -20000 0} {
    set prev [lindex $cyc_loads end]
    for {set i 1} {$i <= 10} {incr i} {
        lappend cyc_loads [expr {$prev + ($peak - $prev) * $i / 10.0}]
    }
}
set nPts [llength $cyc_loads]

set fd_p [open "${prefix}_protocol.txt" w]
foreach v $cyc_loads { puts $fd_p $v }
close $fd_p

timeSeries Path 1 -dt 1.0 -filePath "${prefix}_protocol.txt"
pattern Plain 1 1 { load $centerNode 0 0 1.0 0 0 0 }

constraints Penalty 1e20 1e20; numberer RCM; system BandGeneral
test NormDispIncr 1.0e-6 200 0; algorithm KrylovNewton
integrator LoadControl 1.0; analysis Static

set nFail 0
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
        incr nFail
    }
}
puts "PLATE ($model_type): steps=[expr {$nPts-1}] fails=$nFail"
remove recorders
