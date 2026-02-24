!===============================================================================
!  测试 ADINA 平面应力混凝土本构 (PSUMAT)
!  工况: 单轴拉伸、单轴压缩、纯剪切、往复拉压、往复剪切
!  编译: gfortran -O2 -o test_adina_concrete.exe test_adina_concrete.f90 forumat.f90
!===============================================================================
program test_adina_concrete
    use adina_concrete_constants
    implicit none
    integer,  parameter :: nprops = 13, nstatevs = 33
    real(dp) :: props(nprops), statev(nstatevs)
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3), tangent(3,3)
    real(dp) :: E, nu, ft, fc, epsc, fu, epsu
    integer  :: i, nstep, itest
    real(dp) :: de, e0, e1, g0, g1

    ! 材料参数 (可调)
    E    = 30.0e9_dp      ! 弹性模量 Pa
    nu   = 0.2_dp         ! 泊松比
    ft   = 3.0e6_dp      ! 抗拉强度
    fc   = 30.0e6_dp     ! 抗压强度
    epsc = 0.002_dp      ! 峰值压应变
    fu   = 0.85_dp       ! 极限压应力/峰值
    epsu = 0.004_dp      ! 极限压应变
    props(1)  = E
    props(2)  = nu
    props(3)  = ft
    props(4)  = fc
    props(5)  = epsc
    props(6)  = fu
    props(7)  = epsu
    props(8)  = 0.5_dp   ! beta
    props(9)  = 2.0_dp   ! gama
    props(10) = 0.3_dp   ! rkapa
    props(11) = 0.1_dp   ! alfa
    props(12) = 1.0e-6_dp ! stifac (开裂后拉刚度)
    props(13) = 0.25_dp  ! shefac (剪切传递)

    open(10, file='test_adina_uniax_tension.txt', status='replace')
    open(11, file='test_adina_uniax_compression.txt', status='replace')
    open(12, file='test_adina_shear.txt', status='replace')
    open(13, file='test_adina_cyclic_uniax.txt', status='replace')
    open(14, file='test_adina_cyclic_shear.txt', status='replace')
    write(10,*) 'strain_xx stress_xx'
    write(11,*) 'strain_xx stress_xx'
    write(12,*) 'strain_xy stress_xy'
    write(13,*) 'strain_xx stress_xx'
    write(14,*) 'strain_xy stress_xy'

    ! ---------- 1. 单轴拉伸 (应变 xx 从 0 到 0.0005)
    statev = 0.0_dp
    stress = 0.0_dp
    strain0 = 0.0_dp
    nstep = 100
    de = 0.0005_dp / nstep
    do i = 1, nstep
        strain1(1) = strain0(1) + de
        strain1(2) = strain0(2)
        strain1(3) = strain0(3)
        dstrain(1) = de
        dstrain(2) = 0.0_dp
        dstrain(3) = 0.0_dp
        call PSUMAT(nstatevs, nprops, props, stress, strain0, strain1, dstrain, statev, tangent)
        write(10,*) strain1(1), stress(1)
        strain0 = strain1
    end do
    close(10)

    ! ---------- 2. 单轴压缩 (应变 xx 从 0 到 -0.005，步数多以便与 ADINA 理论曲线对比)
    statev = 0.0_dp
    stress = 0.0_dp
    strain0 = 0.0_dp
    nstep = 500
    de = -0.005_dp / nstep
    do i = 1, nstep
        strain1(1) = strain0(1) + de
        strain1(2) = strain0(2)
        strain1(3) = strain0(3)
        dstrain(1) = de
        dstrain(2) = 0.0_dp
        dstrain(3) = 0.0_dp
        call PSUMAT(nstatevs, nprops, props, stress, strain0, strain1, dstrain, statev, tangent)
        write(11,*) strain1(1), stress(1)
        strain0 = strain1
    end do
    close(11)

    ! ---------- 3. 纯剪切 (应变 xy 从 0 到 0.002)
    statev = 0.0_dp
    stress = 0.0_dp
    strain0 = 0.0_dp
    nstep = 100
    de = 0.002_dp / nstep
    do i = 1, nstep
        strain1(1) = strain0(1)
        strain1(2) = strain0(2)
        strain1(3) = strain0(3) + de
        dstrain(1) = 0.0_dp
        dstrain(2) = 0.0_dp
        dstrain(3) = de
        call PSUMAT(nstatevs, nprops, props, stress, strain0, strain1, dstrain, statev, tangent)
        write(12,*) strain1(3), stress(3)
        strain0 = strain1
    end do
    close(12)

    ! ---------- 4. 往复单轴 (拉-压-拉-压)
    statev = 0.0_dp
    stress = 0.0_dp
    strain0 = 0.0_dp
    nstep = 50
    de = 0.0004_dp / nstep
    do i = 1, nstep
        strain1(1) = strain0(1) + de
        strain1(2) = strain0(2)
        strain1(3) = strain0(3)
        dstrain(1) = de
        dstrain(2) = 0.0_dp
        dstrain(3) = 0.0_dp
        call PSUMAT(nstatevs, nprops, props, stress, strain0, strain1, dstrain, statev, tangent)
        write(13,*) strain1(1), stress(1)
        strain0 = strain1
    end do
    de = -0.0008_dp / nstep
    do i = 1, nstep
        strain1(1) = strain0(1) + de
        strain1(2) = strain0(2)
        strain1(3) = strain0(3)
        dstrain(1) = de
        dstrain(2) = 0.0_dp
        dstrain(3) = 0.0_dp
        call PSUMAT(nstatevs, nprops, props, stress, strain0, strain1, dstrain, statev, tangent)
        write(13,*) strain1(1), stress(1)
        strain0 = strain1
    end do
    de = 0.0008_dp / nstep
    do i = 1, nstep
        strain1(1) = strain0(1) + de
        strain1(2) = strain0(2)
        strain1(3) = strain0(3)
        dstrain(1) = de
        dstrain(2) = 0.0_dp
        dstrain(3) = 0.0_dp
        call PSUMAT(nstatevs, nprops, props, stress, strain0, strain1, dstrain, statev, tangent)
        write(13,*) strain1(1), stress(1)
        strain0 = strain1
    end do
    de = -0.0004_dp / nstep
    do i = 1, nstep
        strain1(1) = strain0(1) + de
        strain1(2) = strain0(2)
        strain1(3) = strain0(3)
        dstrain(1) = de
        dstrain(2) = 0.0_dp
        dstrain(3) = 0.0_dp
        call PSUMAT(nstatevs, nprops, props, stress, strain0, strain1, dstrain, statev, tangent)
        write(13,*) strain1(1), stress(1)
        strain0 = strain1
    end do
    close(13)

    ! ---------- 5. 往复剪切
    statev = 0.0_dp
    stress = 0.0_dp
    strain0 = 0.0_dp
    nstep = 40
    de = 0.001_dp / nstep
    do i = 1, nstep
        strain1(1) = strain0(1)
        strain1(2) = strain0(2)
        strain1(3) = strain0(3) + de
        dstrain(1) = 0.0_dp
        dstrain(2) = 0.0_dp
        dstrain(3) = de
        call PSUMAT(nstatevs, nprops, props, stress, strain0, strain1, dstrain, statev, tangent)
        write(14,*) strain1(3), stress(3)
        strain0 = strain1
    end do
    de = -0.002_dp / nstep
    do i = 1, nstep
        strain1(1) = strain0(1)
        strain1(2) = strain0(2)
        strain1(3) = strain0(3) + de
        dstrain(1) = 0.0_dp
        dstrain(2) = 0.0_dp
        dstrain(3) = de
        call PSUMAT(nstatevs, nprops, props, stress, strain0, strain1, dstrain, statev, tangent)
        write(14,*) strain1(3), stress(3)
        strain0 = strain1
    end do
    de = 0.002_dp / nstep
    do i = 1, nstep
        strain1(1) = strain0(1)
        strain1(2) = strain0(2)
        strain1(3) = strain0(3) + de
        dstrain(1) = 0.0_dp
        dstrain(2) = 0.0_dp
        dstrain(3) = de
        call PSUMAT(nstatevs, nprops, props, stress, strain0, strain1, dstrain, statev, tangent)
        write(14,*) strain1(3), stress(3)
        strain0 = strain1
    end do
    de = -0.001_dp / nstep
    do i = 1, nstep
        strain1(1) = strain0(1)
        strain1(2) = strain0(2)
        strain1(3) = strain0(3) + de
        dstrain(1) = 0.0_dp
        dstrain(2) = 0.0_dp
        dstrain(3) = de
        call PSUMAT(nstatevs, nprops, props, stress, strain0, strain1, dstrain, statev, tangent)
        write(14,*) strain1(3), stress(3)
        strain0 = strain1
    end do
    close(14)

    print *, 'Done. Output files:'
    print *, '  test_adina_uniax_tension.txt'
    print *, '  test_adina_uniax_compression.txt'
    print *, '  test_adina_shear.txt'
    print *, '  test_adina_cyclic_uniax.txt'
    print *, '  test_adina_cyclic_shear.txt'
end program test_adina_concrete
