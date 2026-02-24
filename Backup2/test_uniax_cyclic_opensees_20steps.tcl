# 单轴往复 20 步快速测试 (验证 forumat 嵌入)
wipe
model basic -ndm 3 -ndf 6
set E 30.0e9; set fc 30.0e6; set ft 3.0e6; set epsc 0.002
nDMaterial PlaneStressUserMaterial 1 40 13 $E 0.2 $ft $fc $epsc 0.85 0.004 0.5 2.0 0.3 0.1 1.0e-6 0.25
nDMaterial PlateFromPlaneStress 4 1 [expr $E/(2.0*1.2)]
section LayeredShell 1 3 4 0.003333 4 0.003334 4 0.003333
node 1 0 0 0; node 2 0 1 0; node 3 1 1 0; node 4 1 0 0
fix 1 1 1 1 1 1; fix 2 1 1 1 1 1; fix 3 0 1 1 1 1; fix 4 0 1 1 1 1
equalDOF 3 4 1
element ShellMITC4 1 1 2 3 4 1
# 20 步位移路径: 0 -> -0.0005 (压)
set f [open disp_20.txt w]
for {set i 1} {$i <= 20} {incr i} { puts $f [expr -0.0005*$i/20.0] }; close $f
timeSeries Path 1 -dt 1.0 -filePath disp_20.txt
pattern Plain 1 1 { sp 3 1 1.0 }
recorder Node -file disp_out_20.txt -time -node 3 -dof 1 disp
recorder Node -file react_out_20.txt -time -node 1 2 -dof 1 reaction
constraints Plain; numberer RCM; system BandGeneral
test NormDispIncr 1.0e-8 30; algorithm Newton; integrator LoadControl 1.0
analysis Static
set ok [analyze 20]
puts "analyze 20: ok = $ok"
puts "Done. Check disp_out_20.txt and react_out_20.txt"
