!===============================================================================
!  单轴压缩逐步调试：每步打印应变、应力及 PSUMAT 内部关键量
!  编译: gfortran -O2 -o test_debug_compression.exe test_debug_compression.f90 forumat.f90
!  运行: test_debug_compression.exe
!  输出: debug_uniax_compression.txt
!===============================================================================
program test_debug_compression
    use adina_concrete_constants
    implicit none
    integer,  parameter :: nprops = 13, nstatevs = 22
    real(dp) :: props(nprops), statev(nstatevs)
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3), tangent(3,3)
    integer  :: i, nstep
    real(dp) :: de

    props(1)  = 30.0e9_dp
    props(2)  = 0.2_dp
    props(3)  = 3.0e6_dp
    props(4)  = 30.0e6_dp
    props(5)  = 0.002_dp
    props(6)  = 0.85_dp
    props(7)  = 0.004_dp
    props(8)  = 0.5_dp
    props(9)  = 2.0_dp
    props(10) = 0.3_dp
    props(11) = 0.1_dp
    props(12) = 1.0e-6_dp
    props(13) = 0.25_dp

    open(20, file='debug_uniax_compression.txt', status='replace')
    write(20,'(a)') 'step  strain1(1)    stress(1)     statev(1:6)  statev(19:22)'
    write(20,'(a)') '--------------------------------------------------------------'

    statev = 0.0_dp
    statev(22) = 1.0_dp
    stress = 0.0_dp
    strain0 = 0.0_dp
    nstep = 60
    de = -0.003_dp / nstep

    do i = 1, nstep
        strain1(1) = strain0(1) + de
        strain1(2) = strain0(2)
        strain1(3) = strain0(3)
        dstrain(1) = de
        dstrain(2) = 0.0_dp
        dstrain(3) = 0.0_dp
        call PSUMAT(nstatevs, nprops, props, stress, strain0, strain1, dstrain, statev, tangent)
        write(20,'(i4,2(1x,es14.6),4(1x,es10.2),2(1x,es10.2),4(1x,f8.2))') &
            i, strain1(1), stress(1), statev(1), statev(2), statev(3), statev(4), statev(5), statev(6), &
            statev(19), statev(20), statev(21), statev(22)
        strain0 = strain1
    end do
    close(20)
    print *, 'Debug output: debug_uniax_compression.txt'
end program test_debug_compression
