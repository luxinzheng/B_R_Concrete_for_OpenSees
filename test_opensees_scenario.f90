!===============================================================================
! Standalone test replicating OpenSees shell element integration point
! behavior during cyclic pushover (csp3.tcl).
!
! Tests:
!   A. Stress continuity at eps_xx = 0 (with compression history)
!   B. Tangent consistency (numerical tangent vs. returned tangent)
!   C. NR convergence under displacement control
!   D. Multi-layer cyclic wall: tangent positive-definiteness
!===============================================================================
program test_opensees_scenario
    use adina_concrete_mod
    implicit none
    
    integer, parameter :: nstatevs = 13, nprops = 37
    real(dp) :: props(nprops)
    integer :: n_fail_total
    
    call setup_props(props)
    n_fail_total = 0
    
    write(*,'(A)') '============================================================'
    write(*,'(A)') '  OpenSees Scenario Replication Tests'
    write(*,'(A)') '============================================================'
    write(*,*)
    
    call test_A_stress_continuity(props, n_fail_total)
    call test_B_tangent_consistency(props, n_fail_total)
    call test_C_NR_convergence(props, n_fail_total)
    call test_D_cyclic_wall(props, n_fail_total)
    
    write(*,*)
    write(*,'(A)') '============================================================'
    if (n_fail_total == 0) then
        write(*,'(A)') '  ALL TESTS PASSED'
    else
        write(*,'(A,I4,A)') '  ', n_fail_total, ' FAILURES DETECTED'
    end if
    write(*,'(A)') '============================================================'

contains

subroutine setup_props(props)
    real(dp), intent(out) :: props(37)
    props = 0.0d0
    props(1)  = 21.4d9;  props(2)  = 0.2d0
    props(3)  = 2.07d6;  props(4)  = -20.7d6
    props(5)  = -0.002d0; props(6)  = -4.14d6
    props(7)  = -0.006d0; props(8)  = 1.0d-4
    props(9)  = 0.5d0;   props(10) = 0.75d0
    props(11) = 1.0d0;   props(12) = 0.7d0;  props(13) = 0.12d0
    props(14:19) = (/0.0d0, 0.25d0, 0.5d0, 0.75d0, 1.0d0, 1.2d0/)
    props(20:25) = (/1.0d0, 1.4d0, 1.7d0, 2.2d0, 2.5d0, 2.8d0/)
    props(26:31) = (/1.3d0, 1.5d0, 2.0d0, 2.3d0, 2.7d0, 3.2d0/)
    props(32:37) = (/1.25d0, 1.45d0, 1.95d0, 2.25d0, 2.65d0, 3.15d0/)
end subroutine

subroutine call_mat(props, statev, stress, s0, s1, tangent)
    real(dp), intent(in)    :: props(37), s0(3), s1(3)
    real(dp), intent(inout) :: statev(13), stress(3)
    real(dp), intent(out)   :: tangent(3,3)
    real(dp) :: ds(3)
    ds = s1 - s0
    call PSUMAT(13, 37, props, stress, s0, s1, ds, statev, tangent)
end subroutine

!-----------------------------------------------------------------------
! TEST A: Stress continuity at eps_xx = 0 transition
!-----------------------------------------------------------------------
subroutine test_A_stress_continuity(props, n_fail)
    real(dp), intent(in) :: props(37)
    integer, intent(inout) :: n_fail
    
    real(dp) :: statev(13), stress(3), s0(3), s1(3), tangent(3,3)
    real(dp) :: eps_xx, deps, eps_yy
    real(dp) :: sig_prev(3)
    real(dp) :: jump, max_jump_near_zero
    integer  :: i, nsteps
    
    write(*,'(A)') '--- Test A: Stress continuity at eps_xx = 0 ---'
    
    eps_yy = -5.0d-5
    deps = 1.0d-6
    nsteps = 600   ! from -3e-4 to +3e-4
    max_jump_near_zero = 0.0d0
    
    statev = 0.0d0; stress = 0.0d0; s0 = 0.0d0
    
    ! Apply gravity compression
    s1 = (/ -3.0d-4, eps_yy, 0.0d0 /)
    call call_mat(props, statev, stress, s0, s1, tangent)
    s0 = s1; sig_prev = stress
    
    ! Sweep eps_xx from -3e-4 to +3e-4
    do i = 1, nsteps
        eps_xx = -3.0d-4 + dble(i) * deps
        s1 = (/ eps_xx, eps_yy, 0.0d0 /)
        call call_mat(props, statev, stress, s0, s1, tangent)
        
        jump = sqrt((stress(1)-sig_prev(1))**2 + &
                    (stress(2)-sig_prev(2))**2 + &
                    (stress(3)-sig_prev(3))**2)
        
        ! Only track jumps near eps_xx = 0 (the critical transition)
        if (abs(eps_xx) < 5.0d-5) then
            if (jump > max_jump_near_zero) max_jump_near_zero = jump
        end if
        
        if (abs(eps_xx) < 2.0d-6) then
            write(*,'(A,ES12.4,A,3ES14.6,A,ES10.3)') '    eps_xx=', eps_xx, &
                '  sig=', stress(1), stress(2), stress(3), '  jump=', jump
        end if
        
        sig_prev = stress; s0 = s1
    end do
    
    ! Expected jump per step: ~E_sec * deps ≈ E0 * 1e-6 = 21400 Pa
    ! Allow 5x for Poisson coupling effects
    if (max_jump_near_zero > 5.0d0 * props(1) * deps) then
        write(*,'(A,ES12.4,A,ES12.4)') '  FAIL: max jump near zero = ', &
            max_jump_near_zero, ' > threshold ', 5.0d0 * props(1) * deps
        n_fail = n_fail + 1
    else
        write(*,'(A,ES12.4,A)') '  PASS: max jump near zero = ', &
            max_jump_near_zero, ' (smooth)'
    end if
    write(*,*)
end subroutine

!-----------------------------------------------------------------------
! TEST B: Tangent consistency (numerical vs returned tangent)
!
! CRITICAL: State must be saved BEFORE calling PSUMAT for the base case.
! Both base and perturbed calls must start from the SAME old state.
!-----------------------------------------------------------------------
subroutine test_B_tangent_consistency(props, n_fail)
    real(dp), intent(in) :: props(37)
    integer, intent(inout) :: n_fail
    
    real(dp) :: statev(13), stress(3), s0(3), s1(3), tangent(3,3)
    real(dp) :: statev_save(13), stress_save(3)  ! PRE-call state
    real(dp) :: stress_base(3)                    ! base result
    real(dp) :: statev_p(13), stress_p(3), tangent_p(3,3), s1p(3)
    real(dp) :: num_tangent(3,3), pert
    real(dp) :: eps_yy, eps_xx
    real(dp) :: err, max_err, max_err_loc
    integer  :: j, k, n_test
    real(dp) :: test_eps(10)
    
    write(*,'(A)') '--- Test B: Tangent consistency (numerical vs returned) ---'
    
    pert = 1.0d-9
    eps_yy = -5.0d-5
    max_err = 0.0d0
    max_err_loc = 0.0d0
    
    test_eps = (/ -5.0d-4, -2.0d-4, -5.0d-5, -1.0d-5, -1.0d-6, &
                   1.0d-6,  1.0d-5,  5.0d-5,  1.0d-4,  2.0d-4 /)
    
    do n_test = 1, 10
        eps_xx = test_eps(n_test)
        
        ! Initialize fresh state for each test point
        statev = 0.0d0; stress = 0.0d0; s0 = 0.0d0
        
        ! Step 1: Apply gravity (build up compression history)
        s1 = (/ -1.5d-4, eps_yy, 0.0d0 /)
        call call_mat(props, statev, stress, s0, s1, tangent)
        s0 = s1
        
        ! *** SAVE STATE BEFORE the base call ***
        statev_save = statev
        stress_save = stress
        
        ! Step 2: Base call at target strain
        s1 = (/ eps_xx, eps_yy, 0.0d0 /)
        call call_mat(props, statev, stress, s0, s1, tangent)
        stress_base = stress  ! result of base call
        
        ! Step 3: Compute numerical tangent by perturbation
        ! Each perturbed call starts from the SAME old state (statev_save, stress_save)
        num_tangent = 0.0d0
        do j = 1, 3
            s1p = s1
            s1p(j) = s1p(j) + pert
            
            statev_p = statev_save
            stress_p = stress_save
            call call_mat(props, statev_p, stress_p, s0, s1p, tangent_p)
            
            do k = 1, 3
                num_tangent(k,j) = (stress_p(k) - stress_base(k)) / pert
            end do
        end do
        
        ! Compare: max relative error
        err = 0.0d0
        do j = 1, 3
            do k = 1, 3
                if (abs(tangent(k,j)) > 1.0d3) then
                    err = max(err, abs(tangent(k,j) - num_tangent(k,j)) / abs(tangent(k,j)))
                end if
            end do
        end do
        
        if (err > max_err) then
            max_err = err
            max_err_loc = eps_xx
        end if
        
        if (err > 0.10d0) then
            write(*,'(A,ES10.3,A,ES10.3)') '    eps_xx=', eps_xx, &
                '  MISMATCH rel_err=', err
            write(*,'(A)') '    Returned tangent:'
            do k = 1, 3
                write(*,'(A,3ES14.5)') '      ', tangent(k,1), tangent(k,2), tangent(k,3)
            end do
            write(*,'(A)') '    Numerical tangent:'
            do k = 1, 3
                write(*,'(A,3ES14.5)') '      ', num_tangent(k,1), num_tangent(k,2), num_tangent(k,3)
            end do
        else
            write(*,'(A,ES10.3,A,ES10.3)') '    eps_xx=', eps_xx, &
                '  OK  rel_err=', err
        end if
    end do
    
    if (max_err > 0.50d0) then
        write(*,'(A,ES10.3,A,ES10.3)') '  FAIL: max tangent error = ', max_err, &
            ' at eps_xx=', max_err_loc
        n_fail = n_fail + 1
    else
        write(*,'(A,ES10.3)') '  PASS: max tangent error = ', max_err
    end if
    write(*,*)
end subroutine

!-----------------------------------------------------------------------
! TEST C: NR convergence (displacement-controlled uniaxial cycling)
!
! Prescribes eps_xx, iterates eps_yy to satisfy sigma_yy ≈ 0.
!-----------------------------------------------------------------------
subroutine test_C_NR_convergence(props, n_fail)
    real(dp), intent(in) :: props(37)
    integer, intent(inout) :: n_fail
    
    real(dp) :: statev(13), stress(3), s0(3), s1(3), tangent(3,3)
    real(dp) :: statev_save(13), stress_save(3)
    real(dp) :: eps_xx, eps_yy, deps_yy, residual
    integer  :: step, iter, max_iter, n_hard, max_steps
    real(dp) :: eps_path(80), dstep
    logical  :: converged
    
    write(*,'(A)') '--- Test C: NR convergence (displacement-controlled cycling) ---'
    
    dstep = 2.5d-4
    max_steps = 0
    
    ! Phase 1: 0 → +1e-3
    do step = 1, 4
        max_steps = max_steps + 1
        eps_path(max_steps) = dble(step) * dstep
    end do
    ! Phase 2: +1e-3 → -2e-3
    do step = 1, 12
        max_steps = max_steps + 1
        eps_path(max_steps) = 1.0d-3 - dble(step) * dstep
    end do
    ! Phase 3: -2e-3 → +3e-3
    do step = 1, 20
        max_steps = max_steps + 1
        eps_path(max_steps) = -2.0d-3 + dble(step) * dstep
    end do
    ! Phase 4: +3e-3 → -4e-3
    do step = 1, 28
        max_steps = max_steps + 1
        eps_path(max_steps) = 3.0d-3 - dble(step) * dstep
    end do
    
    statev = 0.0d0; stress = 0.0d0; s0 = 0.0d0
    eps_yy = 0.0d0
    max_iter = 50
    n_hard = 0
    
    open(unit=20, file='test_nr_convergence.csv', status='replace')
    write(20,'(A)') 'step,eps_xx,eps_yy,sig_xx,sig_yy,tau_xy,iterations,converged'
    
    do step = 1, max_steps
        eps_xx = eps_path(step)
        converged = .false.
        
        ! Save state before NR iteration
        statev_save = statev
        stress_save = stress
        
        do iter = 1, max_iter
            ! Reset to start-of-step state for each trial
            statev = statev_save
            stress = stress_save
            
            s1 = (/ eps_xx, eps_yy, 0.0d0 /)
            call call_mat(props, statev, stress, s0, s1, tangent)
            
            residual = stress(2)  ! sigma_yy should be 0
            
            if (abs(residual) < 1.0d0) then  ! 1 Pa tolerance
                converged = .true.
                exit
            end if
            
            if (abs(tangent(2,2)) < 1.0d0) exit  ! singular
            deps_yy = -residual / tangent(2,2)
            eps_yy = eps_yy + deps_yy
        end do
        
        if (.not. converged) then
            n_hard = n_hard + 1
            if (n_hard <= 5) then
                write(*,'(A,I4,A,ES10.3,A,I4,A,ES10.3)') '    Step ', step, &
                    ' eps_xx=', eps_xx, ' NOT converged in ', max_iter, &
                    ' iters, |res|=', abs(residual)
            end if
        else
            ! Accept the converged state
        end if
        
        write(20,'(I4,7(",",ES14.6))') step, eps_xx, eps_yy, &
            stress(1), stress(2), stress(3), dble(iter), dble(merge(1,0,converged))
        
        s0 = s1
    end do
    close(20)
    
    if (n_hard > 0) then
        write(*,'(A,I4,A,I4,A)') '  FAIL: ', n_hard, ' of ', max_steps, &
            ' steps failed to converge'
        n_fail = n_fail + 1
    else
        write(*,'(A,I4,A)') '  PASS: all ', max_steps, ' steps converged'
    end if
    write(*,'(A)') '  Output: test_nr_convergence.csv'
    write(*,*)
end subroutine

!-----------------------------------------------------------------------
! TEST D: Multi-layer cyclic wall (tangent positive-definiteness)
!-----------------------------------------------------------------------
subroutine test_D_cyclic_wall(props, n_fail)
    real(dp), intent(in) :: props(37)
    integer, intent(inout) :: n_fail
    
    integer, parameter :: nlayers = 5
    real(dp) :: statev(13,nlayers), stress(3,nlayers)
    real(dp) :: s0(3,nlayers), s1(3), tangent(3,3)
    real(dp) :: eps_xx, eps_yy_grav, curvature, y_layer
    real(dp) :: eig1, eig2, eig3, det33
    real(dp) :: min_eig
    integer  :: step, layer, n_singular, nsteps
    real(dp) :: curv_path(200), dstep
    
    write(*,'(A)') '--- Test D: Multi-layer cyclic wall simulation ---'
    
    nsteps = 0
    dstep = 2.0d-3
    
    ! Phase 1: 0 → +0.02
    do step = 1, 10
        nsteps = nsteps + 1; curv_path(nsteps) = dble(step) * dstep
    end do
    ! Phase 2: +0.02 → -0.04
    do step = 1, 30
        nsteps = nsteps + 1; curv_path(nsteps) = 0.02d0 - dble(step) * dstep
    end do
    ! Phase 3: -0.04 → +0.06
    do step = 1, 50
        nsteps = nsteps + 1; curv_path(nsteps) = -0.04d0 + dble(step) * dstep
    end do
    ! Phase 4: +0.06 → -0.06
    do step = 1, 60
        nsteps = nsteps + 1; curv_path(nsteps) = 0.06d0 - dble(step) * dstep
    end do
    
    statev = 0.0d0; stress = 0.0d0; s0 = 0.0d0
    eps_yy_grav = -2.0d-5
    n_singular = 0; min_eig = 1.0d30
    
    ! Gravity
    do layer = 1, nlayers
        s1 = (/ 0.0d0, eps_yy_grav, 0.0d0 /)
        call call_mat(props, statev(:,layer), stress(:,layer), &
                     s0(:,layer), s1, tangent)
        s0(:,layer) = s1
    end do
    
    do step = 1, nsteps
        curvature = curv_path(step)
        do layer = 1, nlayers
            y_layer = -0.075d0 + dble(layer-1) * 0.0375d0
            eps_xx = curvature * y_layer
            
            s1 = (/ eps_xx, eps_yy_grav, 0.0d0 /)
            call call_mat(props, statev(:,layer), stress(:,layer), &
                         s0(:,layer), s1, tangent)
            s0(:,layer) = s1
            
            ! Gershgorin eigenvalue bounds
            eig1 = tangent(1,1) - abs(tangent(1,2)) - abs(tangent(1,3))
            eig2 = tangent(2,2) - abs(tangent(2,1)) - abs(tangent(2,3))
            eig3 = tangent(3,3) - abs(tangent(3,1)) - abs(tangent(3,2))
            
            if (min(eig1,eig2,eig3) < min_eig) min_eig = min(eig1,eig2,eig3)
            
            if (eig1 <= 0.0d0 .or. eig2 <= 0.0d0 .or. eig3 <= 0.0d0) then
                n_singular = n_singular + 1
                if (n_singular <= 10) then
                    write(*,'(A,I4,A,I2,A,ES10.3,A,3ES10.2)') &
                        '    Step ', step, ' L', layer, &
                        ' eps=', eps_xx, ' eigs~', eig1, eig2, eig3
                    write(*,'(A,3ES13.4)') '      C:', tangent(1,:)
                    write(*,'(A,3ES13.4)') '      C:', tangent(2,:)
                    write(*,'(A,3ES13.4)') '      C:', tangent(3,:)
                end if
            end if
        end do
    end do
    
    write(*,'(A,ES12.4)') '  Min Gershgorin eigenvalue: ', min_eig
    
    if (n_singular > 0) then
        write(*,'(A,I6,A)') '  FAIL: ', n_singular, ' non-positive-definite tangents'
        n_fail = n_fail + 1
    else
        write(*,'(A)') '  PASS: all tangents are positive definite'
    end if
    write(*,*)
end subroutine

end program test_opensees_scenario
