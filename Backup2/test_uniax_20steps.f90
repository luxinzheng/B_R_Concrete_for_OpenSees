!===============================================================================
!  单轴压缩 0→峰值 恰好 20 步，输出每步应变与应力，供与理论逐点对比
!  编译: gfortran -O2 -o test_uniax_20steps.exe test_uniax_20steps.f90 forumat.f90
!  材料参数与 verify_uniax_adina.py / test_adina_concrete.f90 一致
!===============================================================================
program test_uniax_20steps
    use adina_concrete_constants
    implicit none
    integer,  parameter :: nprops = 13, nstatevs = 33, nstep = 20
    real(dp) :: props(nprops), statev(nstatevs)
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3), tangent(3,3)
    real(dp) :: E, nu, ft, fc, epsc
    integer  :: i
    real(dp) :: de, eps_end

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

    eps_end = -epsc
    de = eps_end / real(nstep, dp)

    statev = 0.0_dp
    stress = 0.0_dp
    strain0 = 0.0_dp
    open(20, file='uniax_20steps.txt', status='replace')
    write(20, '(a)') 'step  strain_xx       stress_xx'

    do i = 1, nstep
        strain1(1) = strain0(1) + de
        strain1(2) = strain0(2)
        strain1(3) = strain0(3)
        dstrain(1) = de
        dstrain(2) = 0.0_dp
        dstrain(3) = 0.0_dp
        call PSUMAT(nstatevs, nprops, props, stress, strain0, strain1, dstrain, statev, tangent)
        write(20, '(i4,2(1x,es16.8))') i, strain1(1), stress(1)
        strain0 = strain1
    end do
    close(20)
    print *, 'Done. Output: uniax_20steps.txt (20 steps 0 -> peak strain -epsc)'
end program test_uniax_20steps
