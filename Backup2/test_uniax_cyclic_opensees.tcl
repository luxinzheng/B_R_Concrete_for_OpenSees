# ==============================================================================
# OpenSees 单轴往复拉压测试 — 用于验证 forumat.f90 (PlaneStressUserMaterial) 嵌入
# 单层壳单元，左端固定，右端施加 x 向往复位移，得到单轴应力状态 (sigma_xx)
# 运行: openSees.exe test_uniax_cyclic_opensees.tcl
# 注意: 若进程退出时崩溃(如 -1073740940)，recorder 文件可能为空，但 "Analysis converged" 表示本构调用正常
# ==============================================================================
wipe
model basic -ndm 3 -ndf 6

# ---------- 材料 (与 forumat.f90 一致: 13 参数, nStatevs>=21) ----------
# props: E, VNU, ft, fc, epsc, fu, epsu, beta, gama, rkapa, alfa, stifac, shefac
set E       30.0e9
set fc      30.0e6
set ft      3.0e6
set epsc    0.002
nDMaterial PlaneStressUserMaterial 1 40 13 $E 0.2 $ft $fc $epsc 0.85 0.004 0.5 2.0 0.3 0.1 1.0e-6 0.25

# 平面应力 -> 壳 (出平面剪切模量)
nDMaterial PlateFromPlaneStress 4 1 [expr $E/(2.0*1.2)]

# ---------- 截面: 3 层混凝土壳 (OpenSees 2.4 LayeredShell 要求 nLayers>2), 总厚 0.01 m ----------
section LayeredShell 1 3 4 0.003333 4 0.003334 4 0.003333

# ---------- 几何: 单单元 1m x 1m, 左端固定, 右端施加 x 向位移 ----------
# 节点顺序: 1(0,0,0) 2(0,1,0) 3(1,1,0) 4(1,0,0) -> x 为加载方向, 长度 L=1
node 1 0.0 0.0 0.0
node 2 0.0 1.0 0.0
node 3 1.0 1.0 0.0
node 4 1.0 0.0 0.0

# 单轴: 左端固定, 右端仅 x 向自由 (ndf=6: ux uy uz rx ry rz)
fix 1 1 1 1 1 1
fix 2 1 1 1 1 1
fix 3 0 1 1 1 1
fix 4 0 1 1 1 1

# 右端两节点 x 向同位移 (单轴均匀应变)
equalDOF 3 4 1

# 壳单元 (4 节点, 截面 1)
element ShellMITC4 1 1 2 3 4 1

# ---------- 单轴往复: 生成位移路径文件 (与 test_uniax_cyclic 6 阶段 110 步一致) ----------
set L 1.0
set eps_peak [expr $ft/$E]
set e1 -0.00052
set e2 $eps_peak
set e3 [expr -1.5*$epsc]
set e4 [expr 2.0*$eps_peak]
set e5 [expr -2.5*$epsc]
set e6 [expr 5.0*$eps_peak]
set n1 10; set n2 20; set n3 20; set n4 20; set n5 20; set n6 20
set f [open disp_path_cyclic.txt w]
set strain_now 0.0
# Phase 1
for {set i 1} {$i <= $n1} {incr i} { set d [expr $e1*$i/$n1*$L]; puts $f $d }
set strain_now $e1
# Phase 2
for {set i 1} {$i <= $n2} {incr i} { set strain_now [expr $strain_now + ($e2-$e1)/$n2]; puts $f [expr $strain_now*$L] }
set strain_now $e2
# Phase 3
for {set i 1} {$i <= $n3} {incr i} { set strain_now [expr $strain_now + ($e3-$e2)/$n3]; puts $f [expr $strain_now*$L] }
set strain_now $e3
# Phase 4
for {set i 1} {$i <= $n4} {incr i} { set strain_now [expr $strain_now + ($e4-$e3)/$n4]; puts $f [expr $strain_now*$L] }
set strain_now $e4
# Phase 5
for {set i 1} {$i <= $n5} {incr i} { set strain_now [expr $strain_now + ($e5-$e4)/$n5]; puts $f [expr $strain_now*$L] }
set strain_now $e5
# Phase 6
for {set i 1} {$i <= $n6} {incr i} { set strain_now [expr $strain_now + ($e6-$e5)/$n6]; puts $f [expr $strain_now*$L] }
close $f
puts "Created disp_path_cyclic.txt (110 displacement values)"

# ---------- 位移由 Path 时间序列驱动: sp 参考值 1.0 * factor = 位移 ----------
timeSeries Path 1 -dt 1.0 -filePath disp_path_cyclic.txt
pattern Plain 1 1 { sp 3 1 1.0 }

# ---------- 每步写入文件并 flush，避免进程退出崩溃时 recorder 为空 ----------
set L 1.0
set area 0.01
set res [open opensees_uniax_cyclic_110steps.txt w]
puts $res "step  strain_xx       stress_xx   phase"
flush $res
# 预读位移路径作为应变来源（nodeDisp 在本约束下可能不更新）
set disp_f [open disp_path_cyclic.txt r]
set disp_list [split [read $disp_f] \n]
close $disp_f

constraints Transformation
numberer RCM
system BandGeneral
test NormDispIncr 1.0e-8 50
algorithm Newton
integrator LoadControl 1.0
analysis Static

# 阶段边界 (与 disp_path 一致): n1=10, n2=20, n3=20, n4=20, n5=20, n6=20
set phase_bound {0 10 30 50 70 90 110}
puts "Running 110 steps (6 phases)..."
set ok 0
for {set step 1} {$step <= 110} {incr step} {
    set ok [analyze 1]
    if {$ok != 0} {
        puts "WARNING: did not converge at step $step, ok = $ok"
        break
    }
    reactions
    set disp [string trim [lindex $disp_list [expr {$step - 1}]]]
    if {$disp == ""} { set disp 0.0 }
    set r1 [nodeReaction 1 1]
    set r2 [nodeReaction 2 1]
    set strain_xx [expr {$disp / $L}]
    # 反力为约束施加的力，取负得截面轴力 (压负拉正)
set stress_xx [expr {-1.0 * ($r1 + $r2) / $area}]
    set ph 1
    for {set p 1} {$p <= 6} {incr p} {
        if {$step <= [lindex $phase_bound $p]} { set ph $p; break }
    }
    puts $res [format "%4d  %16.8e  %16.8e  %4d" $step $strain_xx $stress_xx $ph]
    flush $res
}
close $res
if {$ok == 0} {
    puts "Analysis converged. Results in opensees_uniax_cyclic_110steps.txt"
} else {
    puts "Stopped at step. Partial results in opensees_uniax_cyclic_110steps.txt"
}
puts "Done. Stress_xx = (R1+R2)/0.01 Pa; strain_xx = disp/L (L=1)"
