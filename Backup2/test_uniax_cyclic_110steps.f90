!===============================================================================
!  单轴拉压往复：6 阶段共 110 步
!  1) 10步 压缩到峰值压应力50%； 2) 20步 拉伸到峰值拉应变；
!  3) 20步 压缩到峰值压应变150%； 4) 20步 拉伸到峰值拉应变200%；
!  5) 20步 压缩到峰值压应变250%； 6) 20步 拉伸到峰值拉应变500%
!  编译: gfortran -O2 -o test_uniax_cyclic_110steps.exe test_uniax_cyclic_110steps.f90 forumat.f90
!===============================================================================
program test_uniax_cyclic_110steps
    use adina_concrete_constants
    implicit none
    integer,  parameter :: nprops = 13, nstatevs = 33
    integer,  parameter :: n(6) = [10, 20, 20, 20, 20, 20]
    real(dp) :: props(nprops), statev(nstatevs)
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3), tangent(3,3)
    real(dp) :: E, nu, ft, fc, epsc, eps_peak
    real(dp) :: strain_end(6), de
    integer  :: i, phase, k, nstep

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

    eps_peak = ft / E
    ! 各阶段终点应变 (与 get_cyclic_targets.py 一致: 50%fc 对应约 -0.0005238)
    strain_end(1) = -0.0005238057850097551_dp   ! 50% 峰值压应力
    strain_end(2) = eps_peak                    ! 峰值拉应变
    strain_end(3) = -1.5_dp * epsc               ! 150% 峰值压应变
    strain_end(4) = 2.0_dp * eps_peak            ! 200% 峰值拉应变
    strain_end(5) = -2.5_dp * epsc               ! 250% 峰值压应变
    strain_end(6) = 5.0_dp * eps_peak            ! 500% 峰值拉应变

    statev = 0.0_dp
    stress = 0.0_dp
    strain0 = 0.0_dp
    nstep = 0

    open(30, file='uniax_cyclic_110steps.txt', status='replace')
    write(30, '(a)') 'step  strain_xx       stress_xx   phase'

    do phase = 1, 6
        de = (strain_end(phase) - strain0(1)) / real(n(phase), dp)
        do k = 1, n(phase)
            strain1(1) = strain0(1) + de
            strain1(2) = strain0(2)
            strain1(3) = strain0(3)
            dstrain(1) = de
            dstrain(2) = 0.0_dp
            dstrain(3) = 0.0_dp
            call PSUMAT(nstatevs, nprops, props, stress, strain0, strain1, dstrain, statev, tangent)
            nstep = nstep + 1
            write(30, '(i4,2(1x,es16.8),i4)') nstep, strain1(1), stress(1), phase
            strain0 = strain1
        end do
    end do

    close(30)
    print *, 'Done. Output: uniax_cyclic_110steps.txt (110 steps, 6 phases)'
end program test_uniax_cyclic_110steps
