wipe

model basic -ndm 3 -ndf 6

#===============================================================================
# ADINA平面应力混凝土本构模型定义 (修正版)
#===============================================================================

puts ""
puts "=========================================="
puts "初始化ADINA混凝土本构模型..."
puts "=========================================="

#-------------------------------------------------------------------------------
# 混凝土材料参数定义 (C30混凝土)
#-------------------------------------------------------------------------------

# 基本强度参数
set fc      30.0e6      ;# 单轴抗压强度 (Pa) = 30 MPa
set ft      2.0e6       ;# 抗拉强度 (Pa) = 2.0 MPa
set Ec      30.0e9      ;# 弹性模量 (Pa) = 30 GPa
set nu      0.2         ;# 泊松比

# 压缩软化参数
set epsc    0.002       ;# 峰值压应变
set fu      15.0e6      ;# 极限压应力 (Pa) = 15 MPa
set epsu    0.0035      ;# 极限压应变

# 本构形状参数
set beta    0.6         ;# 软化参数
set gama    3.0         ;# 形状参数
set rkapa   0.85        ;# 软化开始判据
set alfa    -0.05       ;# 多轴影响系数

# 开裂后刚度折减系数
set stifac  0.02        ;# 法向刚度折减
set shefac  0.2         ;# 剪切刚度折减

puts "混凝土参数:"
puts "  fc   = [expr $fc/1e6] MPa"
puts "  ft   = [expr $ft/1e6] MPa"
puts "  Ec   = [expr $Ec/1e9] GPa"
puts "  nu   = $nu"

#-------------------------------------------------------------------------------
# 定义平面应力混凝土材料
#-------------------------------------------------------------------------------
# 注意: OpenSees的PlaneStressUserMaterial语法
# nDMaterial PlaneStressUserMaterial $matTag $nStatev $nProps <$prop1 $prop2 ...>
#
# 重要: 对于用户材料,状态变量数nStatev在OpenSees中的处理方式:
#   - 如果材料子程序内部管理状态变量,此处应设为0
#   - 如果需要OpenSees分配空间,此处设为实际数量
#
# 根据错误信息,这里应该设置为实际需要的状态变量数

# 方案1: 让材料子程序自己管理状态变量 (推荐)
# 设置 nStatev = 0,让Fortran子程序内部管理所有状态变量
nDMaterial PlaneStressUserMaterial 1 0 13 \
    $Ec $nu $ft $fc $epsc $fu $epsu \
    $beta $gama $rkapa $alfa $stifac $shefac

puts "混凝土材料定义完成 (材料标签=1, 状态变量=内部管理, 参数=13)"

#-------------------------------------------------------------------------------
# 定义平板材料 (面外剪切)
#-------------------------------------------------------------------------------
set G_out [expr $Ec / (2.0 * (1.0 + $nu))]
puts "Out-of-plane shear modulus G = [expr $G_out/1e9] GPa"

nDMaterial PlateFromPlaneStress 4 1 $G_out

#===============================================================================
# 钢筋材料定义
#===============================================================================

puts ""
puts "定义钢筋材料..."

# 纵向受力钢筋 d=10mm (约束区)
uniaxialMaterial Steel02 7 379e6 202.7e9 0.01 18.5 0.925 0.15
puts "  材料7: d=10mm纵筋, fy=379 MPa"

# 箍筋和中部纵筋 d=6mm
uniaxialMaterial Steel02 8 392e6 200.6e9 0.01 18.5 0.925 0.15
puts "  材料8: d=6mm钢筋, fy=392 MPa"

# 钢筋层定义
nDMaterial PlateRebar 9  7  90    ;# d=10mm纵筋 (角度90)
nDMaterial PlateRebar 10 8  90    ;# d=6mm纵筋 (角度90)
nDMaterial PlateRebar 11 8  0     ;# d=6mm箍筋 (角度0)

#===============================================================================
# 分层壳截面定义
#===============================================================================

puts ""
puts "定义分层壳截面..."

# 约束区截面 (10层, 总厚125mm)
section LayeredShell 1 10 \
    4  0.0125      \
    11 0.0002403   \
    11 0.0003676   \
    4  0.024696    \
    4  0.024696    \
    4  0.024696    \
    4  0.024696    \
    11 0.0003676   \
    11 0.0002403   \
    4  0.0125
puts "  截面1: 约束区, 10层, 125mm"

# 中部区域截面 (8层, 总厚125mm)
section LayeredShell 2 8 \
    4  0.0125      \
    11 0.0002403   \
    10 0.0002356   \
    4  0.0495241   \
    4  0.0495241   \
    10 0.0002356   \
    11 0.0002403   \
    4  0.0125
puts "  截面2: 中部区, 8层, 125mm"

#===============================================================================
# 读取节点和单元定义
#===============================================================================

puts ""
puts "读取节点和单元..."

source node.tcl
puts "  节点定义完成"

source element.tcl
puts "  壳单元定义完成"

#===============================================================================
# 定义桁架单元 (约束区纵向钢筋)
#===============================================================================

puts ""
puts "定义桁架单元 (纵向钢筋)..."

# 第1列纵筋
element truss 41  1  6 223.53e-6 7
element truss 42  6 11 223.53e-6 7
element truss 43 11 16 223.53e-6 7
element truss 44 16 21 223.53e-6 7
element truss 45 21 26 223.53e-6 7
element truss 46 26 31 223.53e-6 7
element truss 47 31 36 223.53e-6 7
element truss 48 36 41 223.53e-6 7
element truss 49 41 46 223.53e-6 7
element truss 50 46 51 223.53e-6 7

# 第2列纵筋
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

# 第4列纵筋
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

# 第5列纵筋
element truss 71  5 10 223.53e-6 7
element truss 72 10 15 223.53e-6 7
element truss 73 15 20 223.53e-6 7
element truss 74 20 25 223.53e-6 7
element truss 75 25 30 223.53e-6 7
element truss 76 30 35 223.53e-6 7
element truss 77 35 40 223.53e-6 7
element truss 78 40 45 223.53e-6 7
element truss 79 45 50 223.53e-6 7
element truss 80 50 55 223.53e-6 7

puts "  桁架单元定义完成 (40根纵筋)"

#===============================================================================
# 边界条件
#===============================================================================

puts ""
puts "施加边界条件..."

fixY 0.0 1 1 1 1 1 1
puts "  底部节点全固定"

#===============================================================================
# 记录器定义
#===============================================================================

puts ""
puts "设置记录器..."

recorder Node -file 1.txt -time -node 1 2 3 4 5 -dof 1 reaction
puts "  反力记录器: 1.txt"

#===============================================================================
# 分析步骤1: 重力荷载
#===============================================================================

puts ""
puts "=========================================="
puts "阶段1: 重力荷载分析"
puts "=========================================="

pattern Plain 1 Linear {
    load 53 0 -246000 0 0 0 0
}
puts "节点53施加竖向荷载: -246 kN"

recorder Node -file 53.txt -time -node 53 -dof 1 disp
puts "位移记录器: 53.txt"

puts ""
puts "分析设置:"
constraints Plain
numberer RCM
system BandGeneral
test NormDispIncr 1.0e-6 200 0
puts "  收敛准则: NormDispIncr, tol=1e-6, maxIter=200"

algorithm BFGS -count 100
puts "  算法: BFGS"

integrator LoadControl 0.1
puts "  积分器: LoadControl, step=0.1"

analysis Static
puts "  分析类型: Static"

puts ""
puts "开始重力分析 (10步)..."

set ok [analyze 10]

if {$ok == 0} {
    puts "重力分析成功完成!"
} else {
    puts "警告: 重力分析未完全收敛,尝试调整参数..."
    
    # 尝试更小的步长
    integrator LoadControl 0.05
    set ok [analyze 10]
    
    if {$ok == 0} {
        puts "使用更小步长后分析成功!"
    } else {
        puts "错误: 重力分析失败,请检查模型!"
        exit
    }
}

loadConst -time 0.0
puts "重力荷载固定"

#===============================================================================
# 分析步骤2: 循环水平加载
#===============================================================================

puts ""
puts "=========================================="
puts "阶段2: 循环水平加载"
puts "=========================================="

timeSeries Path 1 -dt 0.1 -filePath shuju2.txt
puts "加载时程: shuju2.txt"

pattern Plain 2 1 {
    sp 53 1 1
}
puts "节点53施加位移控制"

puts ""
puts "更新分析设置:"
constraints Penalty 1e20 1e20
numberer RCM
system BandGeneral
test NormDispIncr 1.0e-5 1000 2
puts "  收敛准则: NormDispIncr, tol=1e-5, maxIter=1000"

algorithm KrylovNewton
puts "  算法: KrylovNewton"

integrator LoadControl 0.1
analysis Static

puts ""
puts "开始循环加载分析 (500步)..."

set currentStep 0
set totalSteps 500
set printInterval 50

for {set i 1} {$i <= $totalSteps} {incr i} {
    set ok [analyze 1]
    
    if {$ok != 0} {
        puts "  步骤 $i 不收敛,尝试调整..."
        
        # 尝试1: 使用NewtonLineSearch
        algorithm NewtonLineSearch
        set ok [analyze 1]
        algorithm KrylovNewton
        
        if {$ok != 0} {
            # 尝试2: 减小步长
            integrator LoadControl 0.05
            set ok [analyze 1]
            integrator LoadControl 0.1
        }
        
        if {$ok != 0} {
            # 尝试3: 使用BFGS
            algorithm BFGS
            set ok [analyze 1]
            algorithm KrylovNewton
        }
        
        if {$ok != 0} {
            puts "  错误: 步骤 $i 分析失败!"
            puts "  已完成 $i 步中的 [expr $i-1] 步"
            break
        } else {
            puts "  步骤 $i 通过调整后收敛"
        }
    }
    
    set currentStep $i
    
    # 定期输出进度
    if {[expr $i % $printInterval] == 0} {
        puts "  进度: $i / $totalSteps 步 ([expr int($i*100.0/$totalSteps)]%)"
    }
}

puts ""
if {$currentStep == $totalSteps} {
    puts "循环加载分析成功完成!"
} else {
    puts "循环加载分析部分完成: $currentStep / $totalSteps 步"
}

puts ""
puts "=========================================="
puts "分析完成"
puts "=========================================="
puts ""
puts "输出文件:"
puts "  1.txt   - 底部节点反力"
puts "  53.txt  - 顶部节点位移"
puts ""

wipe
puts "模型已清除"
puts ""
