wipe
model basic -ndm 3 -ndf 6

# B&R materials (37 params) for Multi-layer_Shell
# Material 1: fc=25.8MPa, E0=24.03GPa (ACI)
nDMaterial PlaneStressUserMaterial    1   40   37   2.40254e+10  0.2  2e+06  -2.58e+07  -0.003  -5.16e+06  -0.021  0.0001  0.5  0.75  1  0.7  0.12  0  0.25  0.5  0.75  1  1.2  1  1.4  1.7  2.2  2.5  2.8  1.3  1.5  2  2.3  2.7  3.2  1.25  1.45  1.95  2.25  2.65  3.15
nDMaterial   PlateFromPlaneStress     2        1  9.5e8

uniaxialMaterial   Steel02  3 403.3e6 2.083e11  0.000 18.5 0.925 0.15
uniaxialMaterial   Steel02  4 366.7e6 2.071e11  0.000 18.5 0.925 0.15
uniaxialMaterial   Steel02  5 350.3e6 2.060e11  0.000 18.5 0.925 0.15
uniaxialMaterial   Steel02  9 366.7e6 2.071e11  0.003 18.5 0.925 0.15

nDMaterial   PlateRebar         6     3     90
nDMaterial   PlateRebar         7     4     90
nDMaterial   PlateRebar         8     5     0
nDMaterial   PlateRebar         10     5     0

# Material 11: fc=22.8MPa, E0=22.59GPa (ACI)
nDMaterial PlaneStressUserMaterial    11   40   37   2.25854e+10  0.2  2e+06  -2.28e+07  -0.002  -4.56e+06  -0.012  0.0001  0.5  0.75  1  0.7  0.12  0  0.25  0.5  0.75  1  1.2  1  1.4  1.7  2.2  2.5  2.8  1.3  1.5  2  2.3  2.7  3.2  1.25  1.45  1.95  2.25  2.65  3.15
nDMaterial   PlateFromPlaneStress     12        11  9.5e8

section   LayeredShell 1 8 12 0.010000 10 0.000277  2 0.019317 2 0.019317 2 0.019317 2 0.019317  10 0.000277 12 0.010000
section   LayeredShell 2 8 12 0.010000  7 0.000251 12 0.019736 12 0.019736 12 0.019736 12 0.019736 7 0.000251  12 0.010000

source node.tcl
source element.tcl
element truss 145 1 10 201.061930e-6 3
element truss 146 10 19 201.061930e-6 3
element truss 147 19 28 201.061930e-6 3
element truss 148 28 37 201.061930e-6 3
element truss 149 37 46 201.061930e-6 3
element truss 150 46 55 201.061930e-6 3
element truss 151 55 64 201.061930e-6 3
element truss 152 64 73 201.061930e-6 3
element truss 153 73 82 201.061930e-6 3
element truss 154 82 91 201.061930e-6 3
element truss 155 91 100 201.061930e-6 3
element truss 156 100 109 201.061930e-6 3
element truss 157 109 118 201.061930e-6 3
element truss 158 118 127 201.061930e-6 3
element truss 159 127 136 201.061930e-6 3
element truss 160 136 145 201.061930e-6 3
element truss 161 145 154 201.061930e-6 3
element truss 162 154 163 201.061930e-6 3

element truss 163 2 11 201.061930e-6 3
element truss 164 11 20 201.061930e-6 3
element truss 165 20 29 201.061930e-6 3
element truss 166 29 38 201.061930e-6 3
element truss 167 38 47 201.061930e-6 3
element truss 168 47 56 201.061930e-6 3
element truss 169 56 65 201.061930e-6 3
element truss 170 65 74 201.061930e-6 3
element truss 171 74 83 201.061930e-6 3
element truss 172 83 92 201.061930e-6 3
element truss 173 92 101 201.061930e-6 3
element truss 174 101 110 201.061930e-6 3
element truss 175 110 119 201.061930e-6 3
element truss 176 119 128 201.061930e-6 3
element truss 177 128 137 201.061930e-6 3
element truss 178 137 146 201.061930e-6 3
element truss 179 146 155 201.061930e-6 3
element truss 180 155 164 201.061930e-6 3

element truss 181 8 17 201.061930e-6 3
element truss 182 17 26 201.061930e-6 3
element truss 183 26 35 201.061930e-6 3
element truss 184 35 44 201.061930e-6 3
element truss 185 44 53 201.061930e-6 3
element truss 186 53 62 201.061930e-6 3
element truss 187 62 71 201.061930e-6 3
element truss 188 71 80 201.061930e-6 3
element truss 189 80 89 201.061930e-6 3
element truss 190 89 98 201.061930e-6 3
element truss 191 98 107 201.061930e-6 3
element truss 192 107 116 201.061930e-6 3
element truss 193 116 125 201.061930e-6 3
element truss 194 125 134 201.061930e-6 3
element truss 195 134 143 201.061930e-6 3
element truss 196 143 152 201.061930e-6 3
element truss 197 152 161 201.061930e-6 3
element truss 198 161 170 201.061930e-6 3

element truss 199 9 18 201.061930e-6 3
element truss 200 18 27 201.061930e-6 3
element truss 201 27 36 201.061930e-6 3
element truss 202 36 45 201.061930e-6 3
element truss 203 45 54 201.061930e-6 3
element truss 204 54 63 201.061930e-6 3
element truss 205 63 72 201.061930e-6 3
element truss 206 72 81 201.061930e-6 3
element truss 207 81 90 201.061930e-6 3
element truss 208 90 99 201.061930e-6 3
element truss 209 99 108 201.061930e-6 3
element truss 210 108 117 201.061930e-6 3
element truss 211 117 126 201.061930e-6 3
element truss 212 126 135 201.061930e-6 3
element truss 213 135 144 201.061930e-6 3
element truss 214 144 153 201.061930e-6 3
element truss 215 153 162 201.061930e-6 3
element truss 216 162 171 201.061930e-6 3

fixY 0.0 1 1 1 1 1 1

# Fix drilling DOF (theta_z) for all free nodes - NLDKGQ has no drilling stiffness
for {set i 10} {$i <= 171} {incr i} {
    fix $i 0 0 0 0 0 1
}

pattern Plain 1 Linear {
load 164 0 -91428.571429 0 0 0 0
load 165 0 -91428.571429 0 0 0 0
load 166 0 -91428.571429 0 0 0 0
load 167 0 -91428.571429 0 0 0 0
load 168 0 -91428.571429 0 0 0 0
load 169 0 -91428.571429 0 0 0 0
load 170 0 -91428.571429 0 0 0 0
}

constraints Plain
numberer RCM
system BandGeneral
test NormDispIncr 1.0e-5 500
algorithm KrylovNewton
integrator LoadControl 0.1
analysis Static

puts "Gravity..."
set gok [analyze 10]
record
if {$gok != 0} { puts "Gravity FAILED" } else { puts "Gravity OK" }
loadConst -time 0.0

recorder Node -file test_disp.txt -time -node  167 -dof 1  disp
recorder Node -file test_react.txt -time -node 1 2 3 4 5 6 7 8 9 -dof 1 reaction

timeSeries Path 2 -dt 0.1 -filePath ml_disp_protocol.txt
pattern Plain 200 2 {
  sp 167 1 1
}

constraints Penalty 1e16 1e16
numberer RCM
system UmfPack
test NormDispIncr 1.0e-4 1000 2
algorithm KrylovNewton
integrator LoadControl 0.1
analysis Static

set totalTime 94.7
set dt 0.1
set step 0; set nFail 0; set nSkip 0; set consecFail 0

while {[getTime] < $totalTime} {
    set ok [analyze 1]; incr step
    if {$ok == 0} {
        set consecFail 0
        if {[expr {$step % 50}] == 0} { record; puts "  step $step t=[getTime] f=$nFail s=$nSkip" }
        continue
    }
    algorithm ModifiedNewton; set ok [analyze 1]; algorithm KrylovNewton
    if {$ok == 0} { set consecFail 0; continue }

    set nSub 10
    integrator LoadControl [expr {$dt/$nSub}]
    test NormDispIncr 1.0e-4 500 0; algorithm Newton
    set ok 0
    for {set i 0} {$i < $nSub} {incr i} { set ok [analyze 1]; if {$ok != 0} break }
    integrator LoadControl $dt; algorithm KrylovNewton; test NormDispIncr 1.0e-4 1000 2
    if {$ok == 0} { set consecFail 0; continue }

    set nSub 100
    integrator LoadControl [expr {$dt/$nSub}]
    test NormDispIncr 1.0e-3 1000 0; algorithm Newton
    set ok 0
    for {set i 0} {$i < $nSub} {incr i} {
        set ok [analyze 1]
        if {$ok != 0} { algorithm ModifiedNewton; set ok [analyze 1]; algorithm Newton }
        if {$ok != 0} break
    }
    integrator LoadControl $dt; algorithm KrylovNewton; test NormDispIncr 1.0e-4 1000 2
    if {$ok == 0} { set consecFail 0; continue }

    set nSub 1000
    integrator LoadControl [expr {$dt/$nSub}]
    algorithm ModifiedNewton -initial; test NormDispIncr 1.0e10 1 0
    for {set i 0} {$i < $nSub} {incr i} { set ok [analyze 1]; if {$ok != 0} break }
    integrator LoadControl $dt; algorithm KrylovNewton; test NormDispIncr 1.0e-4 1000 2
    if {$ok == 0} { incr nSkip; set consecFail 0; continue }

    incr nFail; incr consecFail
    if {$consecFail >= 20} { puts "*** $consecFail consecutive fails -- stopping ***"; break }
}
record
puts "Lateral done: t=[getTime] steps=$step fails=$nFail skips=$nSkip"
