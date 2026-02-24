!===============================================================================
!  Stand-Alone Test Program for ADINA Concrete Plane Stress Model
!  
!  Tests:
!    1. Uniaxial monotonic compression  (sigma_yy = 0, gamma_xy = 0)
!    2. Uniaxial monotonic tension      (sigma_yy = 0, gamma_xy = 0)
!    3. Uniaxial cyclic tension-compression
!
!  For uniaxial test: prescribe eps_xx, iterate to find eps_yy such that
!  sigma_yy = 0 (plane stress condition with free transverse direction).
!===============================================================================
program test_adina_concrete
    use adina_concrete_mod, only: dp
    implicit none
    
    integer, parameter :: NPROPS   = 37
    integer, parameter :: NSTATEVS = 13
    integer, parameter :: MAX_ITER = 50
    real(dp), parameter :: TOL = 1.0d-10
    
    real(dp) :: props(NPROPS)
    
    ! Set up material properties (typical concrete: C30)
    call setup_material(props)
    
    ! Test 1: Monotonic compression
    write(*,*) '============================================'
    write(*,*) ' Test 1: Uniaxial Monotonic Compression'
    write(*,*) '============================================'
    call test_monotonic_compression(props)
    
    ! Test 2: Monotonic tension
    write(*,*) ''
    write(*,*) '============================================'
    write(*,*) ' Test 2: Uniaxial Monotonic Tension'
    write(*,*) '============================================'
    call test_monotonic_tension(props)
    
    ! Test 3: Cyclic loading
    write(*,*) ''
    write(*,*) '============================================'
    write(*,*) ' Test 3: Uniaxial Cyclic Loading'
    write(*,*) '============================================'
    call test_cyclic(props)
    
    write(*,*) ''
    write(*,*) 'All tests completed. Check CSV files for results.'
    
contains

!-----------------------------------------------------------------------
! Setup default material properties for C30 concrete
!-----------------------------------------------------------------------
subroutine setup_material(props)
    real(dp), intent(out) :: props(NPROPS)
    
    props = 0.0d0
    
    ! Basic properties
    props(1)  = 30000.0d0     ! E0 = 30 GPa (in MPa)
    props(2)  = 0.2d0         ! Poisson's ratio
    props(3)  = 3.0d0         ! ft = 3 MPa (tensile strength, positive)
    props(4)  = -30.0d0       ! fc = -30 MPa (compressive strength, negative)
    props(5)  = -0.002d0      ! eps_c = -0.002 (peak compressive strain)
    props(6)  = -6.0d0        ! fu = -6 MPa (ultimate compressive stress)
    props(7)  = -0.005d0      ! eps_u = -0.005 (ultimate compressive strain)
    props(8)  = 1.0d-4        ! STIFAC (stiffness reduction factor)
    props(9)  = 0.5d0         ! SHEFAC (shear retention factor)
    
    ! Biaxial parameters
    props(10) = 0.75d0        ! BETA
    props(11) = 1.0d0         ! GAMA
    props(12) = 0.7d0         ! RKAPA
    props(13) = 0.12d0        ! ALFA
    
    ! Default biaxial failure envelope (ADINA defaults)
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
! Perform one uniaxial load step with Newton iteration for sigma_yy = 0
! Uses damped Newton with bisection fallback for robustness in softening.
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
    real(dp) :: max_step, damping
    integer  :: iter
    
    eps_yy = eps_yy_init  ! Copy to local variable (avoid aliasing)
    abs_res_old = 1.0d30
    
    ! Max step size: proportional to current strain magnitude
    max_step = max(abs(eps_xx) * 0.5d0, 1.0d-5)
    
    do iter = 1, MAX_ITER
        strain1(1) = eps_xx
        strain1(2) = eps_yy
        strain1(3) = 0.0d0
        dstrain = strain1 - strain0
        
        ! Copy state for trial (don't modify committed state)
        stress_trial = stress
        statev_trial = statev
        
        call PSUMAT(NSTATEVS, NPROPS, props, stress_trial, strain0, &
                    strain1, dstrain, statev_trial, tangent)
        
        residual = stress_trial(2)
        if (abs(residual) < TOL) exit
        if (abs(tangent(2,2)) < 1.0d-30) exit
        
        d_eps_yy = -residual / tangent(2,2)
        
        ! Step size limiting (damping)
        if (abs(d_eps_yy) > max_step) then
            d_eps_yy = sign(max_step, d_eps_yy)
        end if
        
        ! If residual increased, reduce step (line search)
        if (abs(residual) > abs_res_old * 1.1d0 .and. iter > 1) then
            d_eps_yy = d_eps_yy * 0.5d0
        end if
        
        abs_res_old = abs(residual)
        eps_yy = eps_yy + d_eps_yy
        if (abs(d_eps_yy) < 1.0d-15) exit
    end do
    
    ! Final commit call with converged eps_yy
    strain1(1) = eps_xx
    strain1(2) = eps_yy
    strain1(3) = 0.0d0
    dstrain = strain1 - strain0
    
    call PSUMAT(NSTATEVS, NPROPS, props, stress, strain0, &
                strain1, dstrain, statev, tangent)
    
    ! Update for next step
    strain0 = strain1
    sig_xx_out = stress(1)
    eps_yy_out = eps_yy
end subroutine

!-----------------------------------------------------------------------
! Test 1: Uniaxial monotonic compression
!-----------------------------------------------------------------------
subroutine test_monotonic_compression(props)
    real(dp), intent(in) :: props(NPROPS)
    
    real(dp) :: statev(NSTATEVS), stress(3), strain0(3)
    real(dp) :: eps_xx, eps_yy_guess, eps_yy, sig_xx
    real(dp) :: eps_max, deps
    integer  :: i, nsteps, iunit
    
    statev  = 0.0d0;  stress  = 0.0d0;  strain0 = 0.0d0
    
    eps_max = -0.005d0;  nsteps = 200
    deps = eps_max / real(nsteps, dp)
    
    iunit = 11
    open(unit=iunit, file='test1_compression.csv', status='replace')
    write(iunit, '(A)') 'eps_xx,sig_xx,eps_yy'
    write(iunit, '(3ES20.10)') 0.0d0, 0.0d0, 0.0d0
    
    eps_yy = 0.0d0  ! Track converged eps_yy across steps
    do i = 1, nsteps
        eps_xx = deps * real(i, dp)
        ! Use previous converged eps_yy as initial guess (much better than -nu*eps_xx)
        eps_yy_guess = eps_yy
        
        call uniaxial_step(eps_xx, eps_yy_guess, stress, strain0, statev, &
                          props, sig_xx, eps_yy)
        
        write(iunit, '(3ES20.10)') eps_xx, sig_xx, eps_yy
        if (mod(i, 50) == 0) then
            write(*,'(A,I4,A,ES12.4,A,ES12.4)') &
                '  Step ', i, ':  eps_xx = ', eps_xx, '  sig_xx = ', sig_xx
        end if
    end do
    
    close(iunit)
    write(*,*) '  Output: test1_compression.csv'
end subroutine

!-----------------------------------------------------------------------
! Test 2: Uniaxial monotonic tension
!-----------------------------------------------------------------------
subroutine test_monotonic_tension(props)
    real(dp), intent(in) :: props(NPROPS)
    
    real(dp) :: statev(NSTATEVS), stress(3), strain0(3)
    real(dp) :: eps_xx, eps_yy_guess, eps_yy, sig_xx
    real(dp) :: eps_max, deps
    integer  :: i, nsteps, iunit
    
    statev  = 0.0d0;  stress  = 0.0d0;  strain0 = 0.0d0
    
    eps_max = 0.001d0;  nsteps = 200
    deps = eps_max / real(nsteps, dp)
    
    iunit = 12
    open(unit=iunit, file='test2_tension.csv', status='replace')
    write(iunit, '(A)') 'eps_xx,sig_xx,eps_yy'
    write(iunit, '(3ES20.10)') 0.0d0, 0.0d0, 0.0d0
    
    eps_yy = 0.0d0
    do i = 1, nsteps
        eps_xx = deps * real(i, dp)
        eps_yy_guess = eps_yy
        
        call uniaxial_step(eps_xx, eps_yy_guess, stress, strain0, statev, &
                          props, sig_xx, eps_yy)
        
        write(iunit, '(3ES20.10)') eps_xx, sig_xx, eps_yy
        if (mod(i, 50) == 0) then
            write(*,'(A,I4,A,ES12.4,A,ES12.4)') &
                '  Step ', i, ':  eps_xx = ', eps_xx, '  sig_xx = ', sig_xx
        end if
    end do
    
    close(iunit)
    write(*,*) '  Output: test2_tension.csv'
end subroutine

!-----------------------------------------------------------------------
! Test 3: Uniaxial cyclic loading (tension-compression)
!-----------------------------------------------------------------------
subroutine test_cyclic(props)
    real(dp), intent(in) :: props(NPROPS)
    
    real(dp) :: statev(NSTATEVS), stress(3), strain0(3)
    real(dp) :: eps_xx, eps_yy_guess, eps_yy, sig_xx, deps
    integer  :: i, nsteps_per_phase, iunit, step_total
    
    statev  = 0.0d0;  stress  = 0.0d0;  strain0 = 0.0d0
    
    iunit = 13
    open(unit=iunit, file='test3_cyclic.csv', status='replace')
    write(iunit, '(A)') 'step,eps_xx,sig_xx,eps_yy'
    write(iunit, '(I6,A,3ES20.10)') 0, ',', 0.0d0, 0.0d0, 0.0d0
    
    nsteps_per_phase = 100
    step_total = 0
    eps_xx = 0.0d0
    eps_yy = 0.0d0  ! Track converged eps_yy across all phases
    
    ! Phase 1: Tension to eps = 0.0003 (near cracking)
    write(*,*) '  Phase 1: Tension to 0.0003'
    deps = 0.0003d0 / real(nsteps_per_phase, dp)
    do i = 1, nsteps_per_phase
        step_total = step_total + 1
        eps_xx = deps * real(i, dp)
        eps_yy_guess = eps_yy
        call uniaxial_step(eps_xx, eps_yy_guess, stress, strain0, statev, &
                          props, sig_xx, eps_yy)
        write(iunit, '(I6,A,3ES20.10)') step_total, ',', eps_xx, sig_xx, eps_yy
    end do
    write(*,'(A,ES12.4,A,ES12.4)') '    eps = ', eps_xx, '  sig = ', sig_xx
    
    ! Phase 2: Unload to eps = 0
    write(*,*) '  Phase 2: Unload to 0'
    deps = -eps_xx / real(nsteps_per_phase, dp)
    do i = 1, nsteps_per_phase
        step_total = step_total + 1
        eps_xx = eps_xx + deps
        eps_yy_guess = eps_yy
        call uniaxial_step(eps_xx, eps_yy_guess, stress, strain0, statev, &
                          props, sig_xx, eps_yy)
        write(iunit, '(I6,A,3ES20.10)') step_total, ',', eps_xx, sig_xx, eps_yy
    end do
    write(*,'(A,ES12.4,A,ES12.4)') '    eps = ', eps_xx, '  sig = ', sig_xx
    
    ! Phase 3: Compression to eps = -0.003
    write(*,*) '  Phase 3: Compression to -0.003'
    deps = (-0.003d0 - eps_xx) / real(nsteps_per_phase, dp)
    do i = 1, nsteps_per_phase
        step_total = step_total + 1
        eps_xx = eps_xx + deps
        eps_yy_guess = eps_yy
        call uniaxial_step(eps_xx, eps_yy_guess, stress, strain0, statev, &
                          props, sig_xx, eps_yy)
        write(iunit, '(I6,A,3ES20.10)') step_total, ',', eps_xx, sig_xx, eps_yy
    end do
    write(*,'(A,ES12.4,A,ES12.4)') '    eps = ', eps_xx, '  sig = ', sig_xx
    
    ! Phase 4: Unload to eps = 0
    write(*,*) '  Phase 4: Unload to 0'
    deps = -eps_xx / real(nsteps_per_phase, dp)
    do i = 1, nsteps_per_phase
        step_total = step_total + 1
        eps_xx = eps_xx + deps
        eps_yy_guess = eps_yy
        call uniaxial_step(eps_xx, eps_yy_guess, stress, strain0, statev, &
                          props, sig_xx, eps_yy)
        write(iunit, '(I6,A,3ES20.10)') step_total, ',', eps_xx, sig_xx, eps_yy
    end do
    write(*,'(A,ES12.4,A,ES12.4)') '    eps = ', eps_xx, '  sig = ', sig_xx
    
    ! Phase 5: Tension to eps = 0.0005 (post cracking)
    write(*,*) '  Phase 5: Tension to 0.0005'
    deps = (0.0005d0 - eps_xx) / real(nsteps_per_phase, dp)
    do i = 1, nsteps_per_phase
        step_total = step_total + 1
        eps_xx = eps_xx + deps
        eps_yy_guess = eps_yy
        call uniaxial_step(eps_xx, eps_yy_guess, stress, strain0, statev, &
                          props, sig_xx, eps_yy)
        write(iunit, '(I6,A,3ES20.10)') step_total, ',', eps_xx, sig_xx, eps_yy
    end do
    write(*,'(A,ES12.4,A,ES12.4)') '    eps = ', eps_xx, '  sig = ', sig_xx
    
    ! Phase 6: Compression to eps = -0.004
    write(*,*) '  Phase 6: Compression to -0.004'
    deps = (-0.004d0 - eps_xx) / real(nsteps_per_phase, dp)
    do i = 1, nsteps_per_phase
        step_total = step_total + 1
        eps_xx = eps_xx + deps
        eps_yy_guess = eps_yy
        call uniaxial_step(eps_xx, eps_yy_guess, stress, strain0, statev, &
                          props, sig_xx, eps_yy)
        write(iunit, '(I6,A,3ES20.10)') step_total, ',', eps_xx, sig_xx, eps_yy
    end do
    write(*,'(A,ES12.4,A,ES12.4)') '    eps = ', eps_xx, '  sig = ', sig_xx
    
    close(iunit)
    write(*,*) '  Output: test3_cyclic.csv'
end subroutine

end program test_adina_concrete
