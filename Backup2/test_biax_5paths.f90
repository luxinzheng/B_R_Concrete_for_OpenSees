!===============================================================================
!  双轴单调加载：5 条路径，每路径不少于 30 步，直至 max(|eps_xx|,|eps_yy|) 达 1%
!  路径比例 (拉为正压为负): 1)-1:-1  2)1:1  3)-1:-0.5  4)-1:0.1  5)-1:0.2
!  编译: gfortran -O2 -o test_biax_5paths.exe test_biax_5paths.f90 forumat.f90
!===============================================================================
program test_biax_5paths
    use adina_concrete_constants
    implicit none
    integer,  parameter :: nprops = 13, nstatevs = 33, nstep_per_path = 30
    integer,  parameter :: npaths = 5
    real(dp) :: props(nprops), statev(nstatevs)
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3), tangent(3,3)
    real(dp) :: E, nu, ft, fc, epsc
    real(dp) :: rxx, ryy, t_max, t, dt
    integer  :: i, path, nstep

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

    open(30, file='biax_5paths.txt', status='replace')
    write(30, '(a)') 'path  step  strain_xx       strain_yy       stress_xx       stress_yy'

    do path = 1, npaths
        statev = 0.0_dp
        stress = 0.0_dp
        strain0 = 0.0_dp

        select case (path)
        case (1)
            rxx = -1.0_dp
            ryy = -1.0_dp
        case (2)
            rxx = 1.0_dp
            ryy = 1.0_dp
        case (3)
            rxx = -1.0_dp
            ryy = -0.5_dp
        case (4)
            rxx = -1.0_dp
            ryy = 0.1_dp
        case (5)
            rxx = -1.0_dp
            ryy = 0.2_dp
        case default
            rxx = 0.0_dp
            ryy = 0.0_dp
        end select
        t_max = 0.01_dp / max(abs(rxx), abs(ryy), 1.0e-30_dp)
        dt = t_max / real(nstep_per_path, dp)

        do i = 1, nstep_per_path
            t = real(i, dp) * dt
            strain1(1) = rxx * t
            strain1(2) = ryy * t
            strain1(3) = 0.0_dp
            dstrain(1) = strain1(1) - strain0(1)
            dstrain(2) = strain1(2) - strain0(2)
            dstrain(3) = 0.0_dp
            call PSUMAT(nstatevs, nprops, props, stress, strain0, strain1, dstrain, statev, tangent)
            nstep = i
            write(30, '(i2,1x,i4,4(1x,es16.8))') path, nstep, strain1(1), strain1(2), stress(1), stress(2)
            strain0 = strain1
        end do
    end do

    close(30)
    print *, 'Done. Output: biax_5paths.txt (5 paths x 30 steps)'
end program test_biax_5paths
