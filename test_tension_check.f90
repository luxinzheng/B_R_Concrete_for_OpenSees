program test_tension_check
    implicit none
    integer, parameter :: dp = selected_real_kind(15,307)

    integer, parameter :: nstatevs = 13, nprops = 37
    real(dp) :: props(nprops)
    real(dp) :: stress(3), statev(nstatevs), tangent(3,3)
    real(dp) :: strain0(3), strain1(3), dstrain(3)
    real(dp) :: E0, nu, ft, fc, epsc
    integer  :: i, nsteps
    real(dp) :: eps_yy, eps_xx, d_eps, sig_yy, sig_xx

    ! Material parameters (same as csp3.tcl)
    E0   = 21.4d9
    nu   = 0.2d0
    ft   = 2.07d6
    fc   = -20.7d6
    epsc = -0.002d0

    props = 0.0d0
    props(1)  = E0
    props(2)  = nu
    props(3)  = ft
    props(4)  = fc
    props(5)  = epsc
    props(6)  = -4.14d6
    props(7)  = -0.006d0
    props(8)  = 1.0d-4
    props(9)  = 0.5d0
    props(10) = 0.75d0
    props(11) = 1.0d0
    props(12) = 0.7d0
    props(13) = 0.12d0
    props(14:19) = (/0.0d0, 0.25d0, 0.5d0, 0.75d0, 1.0d0, 1.2d0/)
    props(20:25) = (/1.0d0, 1.4d0, 1.7d0, 2.2d0, 2.5d0, 2.8d0/)
    props(26:31) = (/1.3d0, 1.5d0, 2.0d0, 2.3d0, 2.7d0, 3.2d0/)
    props(32:37) = (/1.25d0, 1.45d0, 1.95d0, 2.25d0, 2.65d0, 3.15d0/)

    ! ========================================
    ! Test 1: Uniaxial Y-TENSION (sigma_xx = 0)
    ! eps_yy from 0 to +0.0005 in 50 steps
    ! eps_xx = -nu * eps_yy (free lateral contraction)
    ! ========================================
    nsteps = 50
    d_eps = 0.0005d0 / nsteps

    stress  = 0.0d0
    statev  = 0.0d0
    strain0 = 0.0d0
    strain1 = 0.0d0

    write(*,'(A)') '========================================'
    write(*,'(A)') 'Test: Uniaxial Y-Tension (sigma_xx=0)'
    write(*,'(A)') '  eps_yy to +0.05%, cracking at ft/E ~0.0097%'
    write(*,'(A)') '========================================'
    write(*,'(A)') '  Step   eps_yy(%)   sig_xx(MPa)   sig_yy(MPa)  ANGLE    CRKSTR1'

    do i = 1, nsteps
        eps_yy = d_eps * i
        eps_xx = -nu * eps_yy

        strain0 = strain1
        strain1(1) = eps_xx
        strain1(2) = eps_yy
        strain1(3) = 0.0d0
        dstrain = strain1 - strain0

        call PSUMAT(nstatevs, nprops, props, stress, strain0, strain1, dstrain, statev, tangent)

        sig_xx = stress(1) / 1.0d6
        sig_yy = stress(2) / 1.0d6

        if (i <= 15 .or. mod(i,10) == 0) then
            write(*,'(I6,F12.4,2F14.4,F10.1,E12.4)') &
                i, eps_yy*100, sig_xx, sig_yy, statev(2), statev(3)
        end if
    end do

    ! ========================================
    ! Test 2: Cyclic: Compress then Tension
    ! Phase A: compress eps_yy to -0.001 (50 steps)
    ! Phase B: release to 0 (50 steps)
    ! Phase C: tension eps_yy to +0.0003 (50 steps)
    ! ========================================
    write(*,'(/,A)') '========================================'
    write(*,'(A)')   'Test: Cyclic Compress-Release-Tension'
    write(*,'(A)')   '========================================'
    write(*,'(A)') '  Step  Phase   eps_yy(%)   sig_xx(MPa)   sig_yy(MPa)  ANGLE    CRKSTR1'

    stress  = 0.0d0
    statev  = 0.0d0
    strain0 = 0.0d0
    strain1 = 0.0d0

    ! Phase A: compress
    nsteps = 50
    do i = 1, nsteps
        eps_yy = -0.001d0 * dble(i) / dble(nsteps)
        eps_xx = -nu * eps_yy

        strain0 = strain1
        strain1(1) = eps_xx; strain1(2) = eps_yy; strain1(3) = 0.0d0
        dstrain = strain1 - strain0

        call PSUMAT(nstatevs, nprops, props, stress, strain0, strain1, dstrain, statev, tangent)

        if (i == 1 .or. i == nsteps) then
            write(*,'(I6,"  CompA",F12.4,2F14.4,F10.1,E12.4)') &
                i, eps_yy*100, stress(1)/1d6, stress(2)/1d6, statev(2), statev(3)
        end if
    end do

    ! Phase B: release
    do i = 1, nsteps
        eps_yy = -0.001d0 * (1.0d0 - dble(i)/dble(nsteps))
        eps_xx = -nu * eps_yy

        strain0 = strain1
        strain1(1) = eps_xx; strain1(2) = eps_yy; strain1(3) = 0.0d0
        dstrain = strain1 - strain0

        call PSUMAT(nstatevs, nprops, props, stress, strain0, strain1, dstrain, statev, tangent)

        if (i == nsteps) then
            write(*,'(I6,"  RelB ",F12.4,2F14.4,F10.1,E12.4)') &
                i, eps_yy*100, stress(1)/1d6, stress(2)/1d6, statev(2), statev(3)
        end if
    end do

    ! Phase C: tension
    do i = 1, nsteps
        eps_yy = 0.0003d0 * dble(i) / dble(nsteps)
        eps_xx = -nu * eps_yy

        strain0 = strain1
        strain1(1) = eps_xx; strain1(2) = eps_yy; strain1(3) = 0.0d0
        dstrain = strain1 - strain0

        call PSUMAT(nstatevs, nprops, props, stress, strain0, strain1, dstrain, statev, tangent)

        sig_xx = stress(1) / 1.0d6
        sig_yy = stress(2) / 1.0d6

        if (i <= 20 .or. mod(i,10) == 0) then
            write(*,'(I6,"  TensC",F12.4,2F14.4,F10.1,E12.4)') &
                i, eps_yy*100, sig_xx, sig_yy, statev(2), statev(3)
        end if
    end do

    write(*,'(/,A)') 'Done.'
end program
