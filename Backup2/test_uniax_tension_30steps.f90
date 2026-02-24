!===============================================================================
!  单轴拉伸：0→峰值应变 10 步，峰值→10倍峰值应变 20 步，共 30 步
!  编译: gfortran -O2 -o test_uniax_tension_30steps.exe test_uniax_tension_30steps.f90 forumat.f90
!  材料参数与 test_adina_concrete.f90 一致，峰值拉应变 = ft/E
!===============================================================================
program test_uniax_tension_30steps
    use adina_concrete_constants
    implicit none
    integer,  parameter :: nprops = 13, nstatevs = 33, n1 = 10, n2 = 20
    real(dp) :: props(nprops), statev(nstatevs)
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3), tangent(3,3)
    real(dp) :: E, nu, ft, fc, epsc, eps_peak
    integer  :: i
    real(dp) :: de

    E    = 30.0e9_dp
    nu   = 0.2_dp
    ft   = 3.0e6_dp
    fc   = 30.0e6_dp
    epsc = 0.002_dp
    props(1)  = E
    props(2)  = nu
    props(3)  = ft
    props(4)  = fc
    props(5)  = epsc
    props(6)  = 0.85_dp
    props(7)  = 0.004_dp
    props(8)  = 0.5_dp
    props(9)  = 2.0_dp
    props(10) = 0.3_dp
    props(11) = 0.1_dp
    props(12) = 1.0e-6_dp
    props(13) = 0.25_dp

    eps_peak = ft / E   ! 峰值拉应变

    statev = 0.0_dp
    stress = 0.0_dp
    strain0 = 0.0_dp

    open(30, file='uniax_tension_30steps.txt', status='replace')
    write(30, '(a)') 'step  strain_xx       stress_xx   phase'

    ! 阶段1: 0 → 峰值应变，10 步
    de = eps_peak / real(n1, dp)
    do i = 1, n1
        strain1(1) = strain0(1) + de
        strain1(2) = strain0(2)
        strain1(3) = strain0(3)
        dstrain(1) = de
        dstrain(2) = 0.0_dp
        dstrain(3) = 0.0_dp
        call PSUMAT(nstatevs, nprops, props, stress, strain0, strain1, dstrain, statev, tangent)
        write(30, '(i4,2(1x,es16.8),i4)') i, strain1(1), stress(1), 1
        strain0 = strain1
    end do

    ! 阶段2: 峰值 → 10*峰值应变，20 步
    de = (10.0_dp * eps_peak - eps_peak) / real(n2, dp)
    do i = 1, n2
        strain1(1) = strain0(1) + de
        strain1(2) = strain0(2)
        strain1(3) = strain0(3)
        dstrain(1) = de
        dstrain(2) = 0.0_dp
        dstrain(3) = 0.0_dp
        call PSUMAT(nstatevs, nprops, props, stress, strain0, strain1, dstrain, statev, tangent)
        write(30, '(i4,2(1x,es16.8),i4)') n1 + i, strain1(1), stress(1), 2
        strain0 = strain1
    end do

    close(30)
    print *, 'Done. Output: uniax_tension_30steps.txt (10 steps 0->peak + 20 steps peak->10*peak)'
end program test_uniax_tension_30steps
