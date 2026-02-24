!===============================================================================
!  纯剪切单调加载：剪应变 0 -> 0.5%（张量 strain_xy 0 -> 0.0025），100 步
!  应变: strain_xx=0, strain_yy=0, strain_xy = gamma/2, gamma=0.005
!  编译: gfortran -O2 -o test_shear.exe test_shear.f90 forumat.f90
!===============================================================================
program test_shear
    use adina_concrete_constants
    implicit none
    integer,  parameter :: nprops = 13, nstatevs = 33, nstep = 100
    real(dp) :: props(nprops), statev(nstatevs)
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3), tangent(3,3)
    real(dp) :: E, nu, ft, fc, epsc, gamma_max, dgamma, gamma
    integer  :: i

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

    gamma_max = 0.005_dp   ! 工程剪应变 0.5%
    dgamma    = gamma_max / real(nstep, dp)

    open(30, file='shear_100steps.txt', status='replace')
    write(30, '(a)') 'step  strain_xy       stress_xx       stress_yy       stress_xy'

    statev = 0.0_dp
    stress = 0.0_dp
    strain0 = 0.0_dp

    do i = 1, nstep
        gamma = real(i, dp) * dgamma
        strain1(1) = 0.0_dp
        strain1(2) = 0.0_dp
        strain1(3) = gamma / 2.0_dp   ! 张量剪应变 = 工程剪应变/2
        dstrain(1) = strain1(1) - strain0(1)
        dstrain(2) = strain1(2) - strain0(2)
        dstrain(3) = strain1(3) - strain0(3)
        call PSUMAT(nstatevs, nprops, props, stress, strain0, strain1, dstrain, statev, tangent)
        write(30, '(i4,4(1x,es16.8))') i, strain1(3), stress(1), stress(2), stress(3)
        strain0 = strain1
    end do

    close(30)
    print *, 'Done. Output: shear_100steps.txt (100 steps, shear strain to 0.5%)'
end program test_shear
