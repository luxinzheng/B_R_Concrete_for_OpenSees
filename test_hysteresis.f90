!===============================================================================
!  Cyclic Hysteresis Test — Increasing Amplitude
!
!  Loading protocol:
!    0 → -0.05% → +0.05% → -0.1% → +0.1% → ... → -1.0% → +1.0%
!  Step size: 0.005% strain per step
!
!  Output: test_hysteresis.csv
!===============================================================================
program test_hysteresis
    use adina_concrete_mod, only: dp
    implicit none
    
    integer, parameter :: NPROPS   = 37
    integer, parameter :: NSTATEVS = 13
    integer, parameter :: MAX_ITER = 80
    real(dp), parameter :: TOL = 1.0d-10
    
    real(dp) :: props(NPROPS)
    real(dp) :: statev(NSTATEVS), stress(3), strain0(3)
    real(dp) :: eps_xx, eps_yy, sig_xx
    real(dp) :: target_eps, amp, deps
    real(dp), parameter :: deps_step = 5.0d-5   ! 0.005% strain per step
    integer, parameter  :: n_amp = 20            ! amplitude levels: 0.05% to 1.0%
    integer :: i, j, nsteps, step_total, iunit
    
    call setup_material(props)
    
    statev  = 0.0d0
    stress  = 0.0d0
    strain0 = 0.0d0
    eps_xx  = 0.0d0
    eps_yy  = 0.0d0
    
    iunit = 20
    open(unit=iunit, file='test_hysteresis.csv', status='replace')
    write(iunit, '(A)') 'step,eps_xx,sig_xx,eps_yy'
    
    step_total = 0
    write(iunit, '(I8,A,3ES22.12)') step_total, ',', 0.0d0, 0.0d0, 0.0d0
    
    write(*,*) '============================================'
    write(*,*) ' Cyclic Hysteresis Test (increasing amplitude)'
    write(*,*) ' Amplitude: 0.05% to 1.0% in 0.05% steps'
    write(*,*) ' Load step: 0.005% strain'
    write(*,*) '============================================'
    
    do i = 1, n_amp
        amp = 5.0d-4 * real(i, dp)   ! 0.05%, 0.1%, ..., 1.0%
        
        ! --- Compression half-cycle: current → -amp ---
        target_eps = -amp
        deps = -deps_step
        nsteps = nint(abs(target_eps - eps_xx) / deps_step)
        
        write(*,'(A,I3,A,F8.4,A,I5,A)') &
            '  Amp ', i, ': compress to ', target_eps*100.0d0, '%, ', nsteps, ' steps'
        
        do j = 1, nsteps
            step_total = step_total + 1
            eps_xx = eps_xx + deps
            call uniaxial_step(eps_xx, eps_yy, stress, strain0, statev, &
                              props, sig_xx, eps_yy)
            write(iunit, '(I8,A,3ES22.12)') step_total, ',', eps_xx, sig_xx, eps_yy
        end do
        
        ! --- Tension half-cycle: current → +amp ---
        target_eps = +amp
        deps = +deps_step
        nsteps = nint(abs(target_eps - eps_xx) / deps_step)
        
        write(*,'(A,I3,A,F8.4,A,I5,A)') &
            '  Amp ', i, ': tension  to ', target_eps*100.0d0, '%, ', nsteps, ' steps'
        
        do j = 1, nsteps
            step_total = step_total + 1
            eps_xx = eps_xx + deps
            call uniaxial_step(eps_xx, eps_yy, stress, strain0, statev, &
                              props, sig_xx, eps_yy)
            write(iunit, '(I8,A,3ES22.12)') step_total, ',', eps_xx, sig_xx, eps_yy
        end do
    end do
    
    close(iunit)
    write(*,'(A,I6,A)') '  Total steps: ', step_total, '  Output: test_hysteresis.csv'
    write(*,*) 'Done.'
    
contains

!-----------------------------------------------------------------------
subroutine setup_material(props)
    real(dp), intent(out) :: props(NPROPS)
    props = 0.0d0
    props(1)  = 30000.0d0     ! E0
    props(2)  = 0.2d0         ! nu
    props(3)  = 3.0d0         ! ft
    props(4)  = -30.0d0       ! fc
    props(5)  = -0.002d0      ! eps_c
    props(6)  = -6.0d0        ! fu
    props(7)  = -0.005d0      ! eps_u
    props(8)  = 1.0d-4        ! STIFAC
    props(9)  = 0.5d0         ! SHEFAC
    props(10) = 0.75d0; props(11) = 1.0d0; props(12) = 0.7d0; props(13) = 0.12d0
    props(14) = 0.0d0;  props(15) = 0.25d0; props(16) = 0.5d0
    props(17) = 0.75d0; props(18) = 1.0d0;  props(19) = 1.2d0
    props(20) = 1.0d0;  props(21) = 1.4d0;  props(22) = 1.7d0
    props(23) = 2.2d0;  props(24) = 2.5d0;  props(25) = 2.8d0
    props(26) = 1.3d0;  props(27) = 1.5d0;  props(28) = 2.0d0
    props(29) = 2.3d0;  props(30) = 2.7d0;  props(31) = 3.2d0
    props(32) = 1.25d0; props(33) = 1.45d0; props(34) = 1.95d0
    props(35) = 2.25d0; props(36) = 2.65d0; props(37) = 3.15d0
end subroutine

!-----------------------------------------------------------------------
subroutine uniaxial_step(eps_xx, eps_yy_init, stress, strain0, statev, &
                          props, sig_xx_out, eps_yy_out)
    real(dp), intent(in)    :: eps_xx, eps_yy_init
    real(dp), intent(inout) :: stress(3), strain0(3), statev(NSTATEVS)
    real(dp), intent(in)    :: props(NPROPS)
    real(dp), intent(out)   :: sig_xx_out, eps_yy_out
    
    real(dp) :: strain1(3), dstrain(3), tangent(3,3)
    real(dp) :: stress_trial(3), statev_trial(NSTATEVS)
    real(dp) :: eps_yy, residual, d_eps_yy, abs_res_old
    real(dp) :: max_step
    integer  :: iter
    
    eps_yy = eps_yy_init
    abs_res_old = 1.0d30
    max_step = max(abs(eps_xx) * 0.5d0, 1.0d-5)
    
    do iter = 1, MAX_ITER
        strain1(1) = eps_xx
        strain1(2) = eps_yy
        strain1(3) = 0.0d0
        dstrain = strain1 - strain0
        
        stress_trial = stress
        statev_trial = statev
        
        call PSUMAT(NSTATEVS, NPROPS, props, stress_trial, strain0, &
                    strain1, dstrain, statev_trial, tangent)
        
        residual = stress_trial(2)
        if (abs(residual) < TOL) exit
        if (abs(tangent(2,2)) < 1.0d-30) exit
        
        d_eps_yy = -residual / tangent(2,2)
        
        if (abs(d_eps_yy) > max_step) then
            d_eps_yy = sign(max_step, d_eps_yy)
        end if
        
        if (abs(residual) > abs_res_old * 1.1d0 .and. iter > 1) then
            d_eps_yy = d_eps_yy * 0.5d0
        end if
        
        abs_res_old = abs(residual)
        eps_yy = eps_yy + d_eps_yy
        if (abs(d_eps_yy) < 1.0d-15) exit
    end do
    
    strain1(1) = eps_xx
    strain1(2) = eps_yy
    strain1(3) = 0.0d0
    dstrain = strain1 - strain0
    
    call PSUMAT(NSTATEVS, NPROPS, props, stress, strain0, &
                strain1, dstrain, statev, tangent)
    
    strain0 = strain1
    sig_xx_out = stress(1)
    eps_yy_out = eps_yy
end subroutine

end program test_hysteresis
