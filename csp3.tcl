wipe

model basic -ndm 3 -ndf 6

# ==============================================================================
# ADINA Concrete Plane-Stress Constitutive Model (forumat.f90)
#
# nDMaterial PlaneStressUserMaterial $matTag $nStatevs $nProps $Prop1...$Propn
#
# Props layout (nProps = 37):                          Units: Pa, dimensionless
#   1: E0       Young's modulus                        = 21.4e9 Pa
#   2: VNU      Poisson's ratio                        = 0.2
#   3: SIGMAT   Uniaxial tensile strength (> 0)        = 2.07e6 Pa
#   4: SIGMAC   Uniaxial compressive strength (< 0)    = -20.7e6 Pa
#   5: EPSC     Strain at peak compressive stress (< 0)= -0.002
#   6: SIGMAU   Ultimate compressive stress (<= 0)     = -4.14e6 Pa
#   7: EPSU     Ultimate compressive strain (< 0)      = -0.006
#   8: STIFAC   Post-crack stiffness factor             = 1.0e-4
#   9: SHEFAC   Shear retention factor (0~1)            = 0.5
#  10: BETA     Biaxial envelope transition param       = 0.75
#  11: GAMA     Strain ratio parameter                  = 1.0
#  12: RKAPA    Anisotropic transition ratio             = 0.7
#  13: ALFA     Drucker-Prager alpha                     = 0.12
#  14-19: SP1(6)   sigma1/fc spline table               = 0.0  0.25 0.5  0.75 1.0  1.2
#  20-25: SP31(6)  failure ratio table 1                 = 1.0  1.4  1.7  2.2  2.5  2.8
#  26-31: SP32(6)  failure ratio table 2                 = 1.3  1.5  2.0  2.3  2.7  3.2
#  32-37: SP33(6)  failure ratio table 3                 = 1.25 1.45 1.95 2.25 2.65 3.15
#
# State variables: nStatevs >= 13
#   1: EVMAX    2: ANGLE    3-5: CRKSTR    6: PGRAV
#   7: init_flag    8: eps_comp_max    9-13: spare
# ==============================================================================

# ---- Concrete material parameters (fc = 20.7 MPa) ----
set fc_abs   20.7e6;          # absolute value of compressive strength (Pa)
set ft       2.07e6;          # tensile strength (Pa)
set E0       21.4e9;          # initial Young's modulus ~ 4700*sqrt(fc_MPa) (Pa)
set nu       0.2;             # Poisson's ratio
set fc       [expr {-$fc_abs}];   # compressive strength (NEGATIVE)
set epsc    -0.002;           # peak compressive strain (NEGATIVE)
set fu       [expr {-0.2*$fc_abs}];  # ultimate stress = -4.14e6 (NEGATIVE)
set epsu    -0.006;           # ultimate compressive strain (NEGATIVE)

nDMaterial PlaneStressUserMaterial 1 13 37 \
    $E0  $nu  $ft  $fc  $epsc  $fu  $epsu \
    1.0e-4  0.5  0.75  1.0  0.7  0.12 \
    0.0  0.25  0.5  0.75  1.0  1.2 \
    1.0  1.4   1.7  2.2   2.5  2.8 \
    1.3  1.5   2.0  2.3   2.7  3.2 \
    1.25 1.45  1.95 2.25  2.65 3.15

# nDMaterial PlateFromPlaneStress $matTag $PlaneStressMatTag $OutOfPlaneShearModulus
# Match reference model: G_out = 1.25e10 Pa (same as csp3-????.tcl)
nDMaterial   PlateFromPlaneStress     4         1              1.25e10

#steel
##d=10 longitudinal reinforced steel in the confined region 
uniaxialMaterial   Steel02        7 379e6   202.7e9  0.01 18.5 0.925 0.15
##d=6  transverse reinforced steel and longitudinal reinforced steel in the middle region
uniaxialMaterial   Steel02        8 392e6   200.6e9  0.01 18.5 0.925 0.15

#angle=90 longitudinal reinforced steel
##d=10
nDMaterial   PlateRebar          9               7     90
##d=6
nDMaterial   PlateRebar         10               8     90

#angle=0 transverse reinforced steel
##d=6
nDMaterial   PlateRebar         11               8     0

                                                            
#confined region is divided into 10 layers??middle region is divided into 8 layers

#confined region
# material    absolute thickness   angle(steel)    material tag
##cover               12.5                              4
##d=6transverse     0.2403                0             11
##d=6transverse     0.3676                0             11
##core              24.696                              4
##core              24.696                              4
##core              24.696                              4
##core              24.696                              4
##d=6transverse     0.3676                0             11
##d=6transverse     0.2403                0             11
##cover               12.5                              4

# section LayeredShell $sectionTag $nLayers $matTag1 $thickness1...$matTagn $thicknessn
section   LayeredShell      1          10     4       0.0125   11 0.0002403  11  0.0003676  4 0.024696 4 0.024696 4 0.024696 4 0.024696  11  0.0003676 11 0.0002403 4 0.0125

#middle region
# material    absolute thickness   angle(steel)    material tag
##cover              12.5                               4
##d=6transverse     0.2403                0             11
##d=6longitudinal   0.2356               90             10
##core             49.5241                              4
##core             49.5241                              4
##d=6longitudinal   0.2356               90             10
##d=6transverse     0.2403                0             11
##cover              12.5                               4

section   LayeredShell 2 8 4 0.0125 11 0.0002403  10  0.0002356 4 0.0495241 4 0.0495241  10  0.0002356 11 0.0002403 4 0.0125

source node.tcl
source element.tcl


#longitudinal reinforced steel in the confined region insert in the shell elements as truss elements
element truss 41 1 6 223.53e-6 7
element truss 42 6 11 223.53e-6 7
element truss 43 11 16 223.53e-6 7
element truss 44 16 21 223.53e-6 7
element truss 45 21 26 223.53e-6 7
element truss 46 26 31 223.53e-6 7
element truss 47 31 36 223.53e-6 7
element truss 48 36 41 223.53e-6 7
element truss 49 41 46 223.53e-6 7
element truss 50 46 51 223.53e-6 7

element truss 51  2  7 223.53e-6 7
element truss 52  7 12 223.53e-6 7
element truss 53 12 17 223.53e-6 7
element truss 54 17 22 223.53e-6 7
element truss 55 22 27 223.53e-6 7
element truss 56 27 32 223.53e-6 7
element truss 57 32 37 223.53e-6 7
element truss 58 37 42 223.53e-6 7
element truss 59 42 47 223.53e-6 7
element truss 60 47 52 223.53e-6 7

element truss 61  4  9 223.53e-6 7
element truss 62  9 14 223.53e-6 7
element truss 63 14 19 223.53e-6 7
element truss 64 19 24 223.53e-6 7
element truss 65 24 29 223.53e-6 7
element truss 66 29 34 223.53e-6 7
element truss 67 34 39 223.53e-6 7
element truss 68 39 44 223.53e-6 7
element truss 69 44 49 223.53e-6 7
element truss 70 49 54 223.53e-6 7

element truss 71 5  10 223.53e-6 7
element truss 72 10 15 223.53e-6 7
element truss 73 15 20 223.53e-6 7
element truss 74 20 25 223.53e-6 7
element truss 75 25 30 223.53e-6 7
element truss 76 30 35 223.53e-6 7
element truss 77 35 40 223.53e-6 7
element truss 78 40 45 223.53e-6 7
element truss 79 45 50 223.53e-6 7
element truss 80 50 55 223.53e-6 7
fixY 0.0 1 1 1 1 1 1

recorder Node -file 1.txt -time -node 1 2 3 4 5 -dof 1 reaction 

pattern Plain 1 Linear {
load 53 0 -246000 0 0 0 0
}

recorder Node -file 53.txt -time -node  53 -dof 1  disp

constraints Plain
numberer RCM
system BandGeneral
test NormDispIncr 1.0e-6 200 ;
algorithm BFGS -count 100
integrator LoadControl 0.1;				
analysis Static				
analyze 10;					

puts "gravity analyze ok..."
loadConst -time 0.0;

timeSeries Path 1 -dt 0.1 -filePath shuju2.txt ;
pattern Plain 2 1 {
  sp 53 1 1
 }


constraints Penalty 1e20 1e20;
numberer RCM;
system BandGeneral;
test NormDispIncr 1.0e-5 1000 2;
algorithm KrylovNewton;
integrator LoadControl 0.1;
analysis Static;

# ---- Adaptive analysis with manual displacement logging ----
set totalTime 50.0
set dt       0.1
set currentTime [getTime]
set ok 0
set step 0
set nFail 0
set nSkip 0
set maxConsecFail 0
set consecFail 0

# Manual displacement logger
set fd [open "disp_manual.txt" w]

proc logStep {fd} {
    set t [getTime]
    set d53 [nodeDisp 53 1]
    set r1 [nodeReaction 1 1]; set r2 [nodeReaction 2 1]
    set r3 [nodeReaction 3 1]; set r4 [nodeReaction 4 1]; set r5 [nodeReaction 5 1]
    puts $fd "$t $d53 $r1 $r2 $r3 $r4 $r5"
}

while {$currentTime < $totalTime} {
    set ok [analyze 1]
    incr step
    
    if {$ok == 0} {
        set consecFail 0
        set currentTime [getTime]
        logStep $fd
        if {[expr {$step % 50}] == 0} {
            flush $fd
            puts "  Step $step, time = $currentTime (fails=$nFail, skips=$nSkip)"
        }
        continue
    }
    
    # --- Level 1: ModifiedNewton ---
    algorithm ModifiedNewton
    set ok [analyze 1]
    algorithm KrylovNewton
    if {$ok == 0} {
        set consecFail 0
        set currentTime [getTime]
        logStep $fd
        continue
    }
    
    # --- Level 2: 10 sub-steps ---
    set nSub 10
    integrator LoadControl [expr {$dt/$nSub}]
    test NormDispIncr 1.0e-4 500 0
    algorithm Newton
    set ok 0
    for {set i 0} {$i < $nSub} {incr i} {
        set ok [analyze 1]
        if {$ok != 0} break
    }
    integrator LoadControl $dt
    algorithm KrylovNewton
    test NormDispIncr 1.0e-5 1000 2
    if {$ok == 0} {
        set consecFail 0
        set currentTime [getTime]
        logStep $fd
        continue
    }
    
    # --- Level 3: 100 sub-steps ---
    set nSub 100
    integrator LoadControl [expr {$dt/$nSub}]
    test NormDispIncr 1.0e-3 1000 0
    algorithm Newton
    set ok 0
    for {set i 0} {$i < $nSub} {incr i} {
        set ok [analyze 1]
        if {$ok != 0} {
            algorithm ModifiedNewton
            set ok [analyze 1]
            algorithm Newton
        }
        if {$ok != 0} break
    }
    integrator LoadControl $dt
    algorithm KrylovNewton
    test NormDispIncr 1.0e-5 1000 2
    if {$ok == 0} {
        set consecFail 0
        set currentTime [getTime]
        logStep $fd
        continue
    }
    
    # --- Level 4: Force-through (1000 sub-steps) ---
    set nSub 1000
    integrator LoadControl [expr {$dt/$nSub}]
    algorithm ModifiedNewton -initial
    test NormDispIncr 1.0e10 1 0
    set ok 0
    for {set i 0} {$i < $nSub} {incr i} {
        set ok [analyze 1]
        if {$ok != 0} break
    }
    integrator LoadControl $dt
    algorithm KrylovNewton
    test NormDispIncr 1.0e-5 1000 2
    if {$ok == 0} {
        incr nSkip
        set consecFail 0
        set currentTime [getTime]
        logStep $fd
        puts "  >> Force-through at time $currentTime, step $step (skip #$nSkip)"
        continue
    }
    
    # All levels failed
    incr nFail
    incr consecFail
    if {$consecFail > $maxConsecFail} { set maxConsecFail $consecFail }
    puts "*** FAILED at time [getTime], step $step (fail #$nFail, consec=$consecFail) ***"
    
    if {$consecFail >= 20} {
        puts "*** $consecFail consecutive failures - stopping ***"
        break
    }
    set currentTime [getTime]
}

flush $fd
close $fd
puts "Analysis finished: time=$currentTime, steps=$step, fails=$nFail, skips=$nSkip, maxConsecFail=$maxConsecFail"