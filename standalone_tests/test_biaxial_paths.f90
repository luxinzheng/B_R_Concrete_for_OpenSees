!===============================================================================
! test_biaxial_paths.f90
!   Test ADINA concrete model under 8 proportional biaxial stress paths.
!===============================================================================
program test_biaxial_paths
    use br_test_params, only: dp, set_props_default
    implicit none

    external :: PSUMAT

    integer, parameter :: NPROPS = 37, NSTATV = 40, NSTEPS = 2000, NPATHS = 8
    real(dp) :: props(NPROPS)
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: statev(NSTATV), tangent(3,3)
    real(dp) :: E0, nu, fc, epsc

    real(dp) :: sig_rat1(NPATHS), sig_rat2(NPATHS)
    real(dp) :: alpha_s(NPATHS), eps_end(NPATHS)
    integer  :: sign_drv(NPATHS)

    real(dp) :: r, de1, e1, e2
    integer  :: ip, istep, iu
    logical  :: drv_is_1

    call set_props_default(props)
    E0   = props(1)
    nu   = props(2)
    fc   = props(4)
    epsc = props(5)

    sig_rat1 = (/ -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, &
                  -1.0_dp, -1.0_dp, -1.0_dp,  1.0_dp /)
    sig_rat2 = (/ -1.0_dp, -0.5_dp, -0.2_dp,  0.0_dp, &
                   0.1_dp,  0.2_dp,  0.5_dp,  1.0_dp /)

    do ip = 1, NPATHS
        if (abs(sig_rat1(ip)) >= abs(sig_rat2(ip))) then
            r = sig_rat2(ip) / sig_rat1(ip)
        else
            r = sig_rat1(ip) / sig_rat2(ip)
        end if
        alpha_s(ip) = (r - nu) / (1.0_dp - r * nu)

        if (sig_rat1(ip) < 0.0_dp) then
            eps_end(ip) = 3.0_dp * epsc
        else
            eps_end(ip) = 5.0_dp * props(3) / E0
        end if
    end do

    open(newunit=iu, file='out_biaxial_paths.csv', status='replace')
    write(iu, '(A)') 'path,step,eps1,eps2,sig1,sig2,sig12,sig1_fc,sig2_fc'

    do ip = 1, NPATHS
        stress  = 0.0_dp
        strain0 = 0.0_dp
        strain1 = 0.0_dp
        statev  = 0.0_dp

        de1 = eps_end(ip) / dble(NSTEPS)
        drv_is_1 = (abs(sig_rat1(ip)) >= abs(sig_rat2(ip)))

        do istep = 1, NSTEPS
            strain0 = strain1

            e1 = de1 * dble(istep)
            if (drv_is_1) then
                strain1(1) = e1
                strain1(2) = alpha_s(ip) * e1
            else
                strain1(2) = e1
                strain1(1) = alpha_s(ip) * e1
            end if
            strain1(3) = 0.0_dp

            dstrain = strain1 - strain0

            call PSUMAT(NSTATV, NPROPS, props, stress, strain0, strain1, &
                        dstrain, statev, tangent)

            write(iu, '(I2,",",I6,",",ES15.7,",",ES15.7,",",ES15.7,",",' // &
                  'ES15.7,",",ES15.7,",",ES15.7,",",ES15.7)') &
                ip, istep, strain1(1), strain1(2), &
                stress(1), stress(2), stress(3), &
                stress(1)/fc, stress(2)/fc

            if (stress(1) /= stress(1)) then
                write(*,'(A,I2,A,I6)') '  Path ', ip, ' NaN at step ', istep
                exit
            end if
        end do

        write(*,'(A,I2,A,F6.2,A,F6.2,A,F10.4,A,F10.4)') &
            '  Path ', ip, ' (', sig_rat1(ip), ':', sig_rat2(ip), &
            ')  sig_min=', minval(stress(1:2)), '  sig_max=', maxval(stress(1:2))
    end do

    close(iu)
    write(*,'(A)') 'Done. Output: out_biaxial_paths.csv'

end program test_biaxial_paths
