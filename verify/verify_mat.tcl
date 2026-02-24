# ==============================================================================
# Common material definitions for verification tests
# Set $model_type to "adina" or "ref" BEFORE sourcing this file.
# ==============================================================================

if {![info exists model_type]} { set model_type "adina" }

set E0    21.4e9
set nu    0.2
set ft    2.07e6
set fc_abs 20.7e6
set fc    [expr {-$fc_abs}]
set epsc  -0.002
set fu    [expr {-0.2*$fc_abs}]
set epsu  -0.006
set G_out [expr {$E0 / (2.0*(1.0+$nu))}]

if {$model_type eq "adina"} {
    nDMaterial PlaneStressUserMaterial 1 13 37 \
        $E0  $nu  $ft  $fc  $epsc  $fu  $epsu \
        1.0e-4  0.5  0.75  1.0  0.7  0.12 \
        0.0  0.25  0.5  0.75  1.0  1.2 \
        1.0  1.4   1.7  2.2   2.5  2.8 \
        1.3  1.5   2.0  2.3   2.7  3.2 \
        1.25 1.45  1.95 2.25  2.65 3.15
} else {
    nDMaterial PlaneStressUserMaterial 1 40 7 \
        $fc_abs  $ft  $fu  $epsc  $epsu  0.001  0.08
}

nDMaterial PlateFromPlaneStress 4 1 1.25e10

# Steel materials (shared by wall/beam tests)
uniaxialMaterial Steel02 7 379e6  202.7e9 0.01 18.5 0.925 0.15
uniaxialMaterial Steel02 8 392e6  200.6e9 0.01 18.5 0.925 0.15
nDMaterial PlateRebar  9  7 90
nDMaterial PlateRebar 10  8 90
nDMaterial PlateRebar 11  8  0
