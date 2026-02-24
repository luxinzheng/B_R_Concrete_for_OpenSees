!===============================================================================
!  Comprehensive standalone test for forumat.f90 ADINA concrete model
!===============================================================================
program test_comprehensive
    use adina_concrete_mod
    implicit none
    
    integer, parameter :: NPROPS = 37, NSTATEV = 13
    real(dp) :: props(NPROPS)
    integer :: npass, nfail
    
    npass = 0; nfail = 0
    
    call setup_props(props)
    
    write(*,'(A)') '========================================================'
    write(*,'(A)') '  ADINA Concrete Model — Comprehensive Standalone Test'
    write(*,'(A)') '========================================================'
    write(*,*)
    
    call test_elastic_response(props, npass, nfail)
    call test_uniaxial_compression(props, npass, nfail)
    call test_uniaxial_tension(props, npass, nfail)
    call test_cyclic_loading(props, npass, nfail)
    call test_compression_tension_transition(props, npass, nfail)
    call test_pure_shear(props, npass, nfail)
    call test_biaxial_compression(props, npass, nfail)
    call test_direction_independence(props, npass, nfail)
    call test_tangent_consistency(props, npass, nfail)
    call test_many_steps_stability(props, npass, nfail)
    
    write(*,*)
    write(*,'(A)') '========================================================'
    write(*,'(A,I4,A,I4,A,I4)') '  TOTAL: ', npass+nfail, ' tests, PASS=', npass, ', FAIL=', nfail
    write(*,'(A)') '========================================================'
    if (nfail > 0) then
        write(*,'(A)') '  *** SOME TESTS FAILED ***'
    else
        write(*,'(A)') '  ALL TESTS PASSED'
    end if
    
contains

subroutine setup_props(props)
    real(dp), intent(out) :: props(NPROPS)
    props = 0.0d0
    props(1)  = 21.4d9
    props(2)  = 0.2d0
    props(3)  = 2.07d6
    props(4)  = -20.7d6
    props(5)  = -0.002d0
    props(6)  = -4.14d6
    props(7)  = -0.006d0
    props(8)  = 1.0d-4
    props(9)  = 0.5d0
    props(10) = 0.75d0
    props(11) = 1.0d0
    props(12) = 0.7d0
    props(13) = 0.12d0
    props(14:19) = (/ 0.0d0, 0.25d0, 0.5d0, 0.75d0, 1.0d0, 1.2d0 /)
    props(20:25) = (/ 1.0d0, 1.1d0,  1.2d0, 1.15d0, 1.0d0, 0.9d0 /)
    props(26:31) = (/ 1.3d0, 1.5d0,  2.0d0, 2.3d0,  2.7d0, 3.2d0 /)
    props(32:37) = (/ 1.25d0, 1.45d0, 1.95d0, 2.25d0, 2.65d0, 3.15d0 /)
end subroutine

subroutine reset_state(stress, strain0, statev)
    real(dp), intent(out) :: stress(3), strain0(3), statev(NSTATEV)
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
end subroutine

subroutine do_step(props, stress, strain0, strain1, statev, tangent)
    real(dp), intent(in) :: props(NPROPS), strain1(3)
    real(dp), intent(inout) :: stress(3), strain0(3), statev(NSTATEV)
    real(dp), intent(out) :: tangent(3,3)
    real(dp) :: dstrain(3)
    dstrain = strain1 - strain0
    call PSUMAT(NSTATEV, NPROPS, props, stress, strain0, strain1, dstrain, statev, tangent)
    strain0 = strain1
end subroutine

subroutine check(test_name, cond, npass, nfail)
    character(len=*), intent(in) :: test_name
    logical, intent(in) :: cond
    integer, intent(inout) :: npass, nfail
    if (cond) then
        npass = npass + 1
    else
        nfail = nfail + 1
        write(*,'(A,A)') '  FAIL: ', test_name
    end if
end subroutine

! ===================================================================
! TEST 0: Elastic response (small strain)
! ===================================================================
subroutine test_elastic_response(props, npass, nfail)
    real(dp), intent(in) :: props(NPROPS)
    integer, intent(inout) :: npass, nfail
    real(dp) :: stress(3), strain0(3), strain1(3), statev(NSTATEV), tangent(3,3)
    real(dp) :: E0, v, C11_expected, sig_expected
    
    write(*,'(A)') 'TEST 0: Elastic response (small strain)'
    E0 = props(1); v = props(2)
    C11_expected = E0 / (1.0d0 - v*v)
    
    call reset_state(stress, strain0, statev)
    strain1 = 0.0d0
    strain1(1) = -1.0d-5
    call do_step(props, stress, strain0, strain1, statev, tangent)
    
    sig_expected = C11_expected * strain1(1)
    call check('Elastic: sig_xx ≈ C11*eps_xx', &
        abs(stress(1) - sig_expected) / abs(sig_expected) < 0.01d0, npass, nfail)
    call check('Elastic: tangent(1,1) ≈ C11', &
        abs(tangent(1,1) - C11_expected) / C11_expected < 0.01d0, npass, nfail)
    call check('Elastic: sig_yy has Poisson effect', &
        abs(stress(2)) > 0.0d0, npass, nfail)
    
    write(*,'(A,ES12.4,A,ES12.4)') '  sig_xx=', stress(1), ' expected=', sig_expected
    write(*,'(A,ES12.4,A,ES12.4)') '  C11=', tangent(1,1), ' expected=', C11_expected
    write(*,*)
end subroutine

! ===================================================================
! TEST 1: Uniaxial compression (total-strain approach)
! ===================================================================
subroutine test_uniaxial_compression(props, npass, nfail)
    real(dp), intent(in) :: props(NPROPS)
    integer, intent(inout) :: npass, nfail
    real(dp) :: stress(3), strain0(3), strain1(3), statev(NSTATEV), tangent(3,3)
    real(dp) :: sig_peak, eps_peak
    integer :: i
    
    write(*,'(A)') 'TEST 1: Uniaxial compression'
    call reset_state(stress, strain0, statev)
    
    sig_peak = 0.0d0; eps_peak = 0.0d0
    
    do i = 1, 50
        strain1 = 0.0d0
        strain1(1) = -0.0001d0 * i
        call do_step(props, stress, strain0, strain1, statev, tangent)
        
        if (stress(1) < sig_peak) then
            sig_peak = stress(1)
            eps_peak = strain1(1)
        end if
    end do
    
    call check('Comp: peak stress within 20% of fc', &
        abs(sig_peak - props(4)) / abs(props(4)) < 0.2d0, npass, nfail)
    call check('Comp: stress is negative', sig_peak < 0.0d0, npass, nfail)
    call check('Comp: no NaN in stress', .not. any(stress /= stress), npass, nfail)
    call check('Comp: no NaN in tangent', .not. any(tangent /= tangent), npass, nfail)
    call check('Comp: post-peak softening', stress(1) > sig_peak, npass, nfail)
    
    write(*,'(A,ES12.4,A,ES12.4)') '  Peak stress=', sig_peak, ' at eps=', eps_peak
    write(*,'(A,ES12.4)') '  Final stress=', stress(1)
    write(*,*)
end subroutine

! ===================================================================
! TEST 2: Uniaxial tension (cracking + softening)
! ===================================================================
subroutine test_uniaxial_tension(props, npass, nfail)
    real(dp), intent(in) :: props(NPROPS)
    integer, intent(inout) :: npass, nfail
    real(dp) :: stress(3), strain0(3), strain1(3), statev(NSTATEV), tangent(3,3)
    real(dp) :: E0, sig_max, eps_tu
    integer :: i
    logical :: found_crack, stress_decreased
    
    write(*,'(A)') 'TEST 2: Uniaxial tension (cracking + softening)'
    E0 = props(1)
    eps_tu = 10.0d0 * props(3) / E0
    call reset_state(stress, strain0, statev)
    
    sig_max = 0.0d0
    found_crack = .false.
    stress_decreased = .false.
    
    do i = 1, 200
        strain1 = 0.0d0
        strain1(1) = 1.0d-5 * i
        call do_step(props, stress, strain0, strain1, statev, tangent)
        
        if (stress(1) > sig_max .and. .not. found_crack) sig_max = stress(1)
        if (statev(2) < 361.0d0 .and. .not. found_crack) then
            found_crack = .true.
            write(*,'(A,I4,A,ES12.4)') '  Crack at step ', i, ', strain=', strain1(1)
        end if
        if (found_crack .and. stress(1) < 0.3d0 * sig_max) then
            stress_decreased = .true.
        end if
    end do
    
    call check('Tens: cracking detected', found_crack, npass, nfail)
    call check('Tens: peak stress near ft', abs(sig_max - props(3)) / props(3) < 0.3d0, npass, nfail)
    call check('Tens: stress softens after cracking', stress_decreased, npass, nfail)
    call check('Tens: final stress near 0', abs(stress(1)) < 0.5d0 * props(3), npass, nfail)
    call check('Tens: no NaN', .not. any(stress /= stress), npass, nfail)
    
    write(*,'(A,ES12.4,A,ES12.4)') '  Peak=', sig_max, ', Final=', stress(1)
    write(*,'(A,ES12.4)') '  eps_tu=', eps_tu
    write(*,*)
end subroutine

! ===================================================================
! TEST 3: Cyclic compression-tension-compression
! ===================================================================
subroutine test_cyclic_loading(props, npass, nfail)
    real(dp), intent(in) :: props(NPROPS)
    integer, intent(inout) :: npass, nfail
    real(dp) :: stress(3), strain0(3), strain1(3), statev(NSTATEV), tangent(3,3)
    real(dp) :: sig_at_peak, sig_unloaded, sig_reloaded
    integer :: i
    logical :: has_nan
    
    write(*,'(A)') 'TEST 3: Cyclic compression-unload-reload'
    call reset_state(stress, strain0, statev)
    has_nan = .false.
    
    ! Phase 1: compress to ε=-0.002 (around peak)
    do i = 1, 20
        strain1 = 0.0d0
        strain1(1) = -0.0001d0 * i
        call do_step(props, stress, strain0, strain1, statev, tangent)
        if (any(stress /= stress) .or. any(tangent /= tangent)) has_nan = .true.
    end do
    sig_at_peak = stress(1)
    write(*,'(A,ES12.4)') '  After compression to -0.002: sig=', sig_at_peak
    
    ! Phase 2: unload to ε=0
    do i = 1, 20
        strain1 = 0.0d0
        strain1(1) = -0.002d0 + 0.0001d0 * i
        call do_step(props, stress, strain0, strain1, statev, tangent)
        if (any(stress /= stress) .or. any(tangent /= tangent)) has_nan = .true.
    end do
    sig_unloaded = stress(1)
    write(*,'(A,ES12.4)') '  After unload to 0: sig=', sig_unloaded
    
    ! Phase 3: reload to ε=-0.002
    do i = 1, 20
        strain1 = 0.0d0
        strain1(1) = -0.0001d0 * i
        call do_step(props, stress, strain0, strain1, statev, tangent)
        if (any(stress /= stress) .or. any(tangent /= tangent)) has_nan = .true.
    end do
    sig_reloaded = stress(1)
    write(*,'(A,ES12.4)') '  After reload to -0.002: sig=', sig_reloaded
    
    call check('Cyclic: no NaN', .not. has_nan, npass, nfail)
    call check('Cyclic: unloaded stress near 0', abs(sig_unloaded) < 0.1d0 * abs(sig_at_peak), npass, nfail)
    call check('Cyclic: reload stress ≈ peak', &
        abs(sig_reloaded - sig_at_peak) / abs(sig_at_peak) < 0.2d0, npass, nfail)
    write(*,*)
end subroutine

! ===================================================================
! TEST 4: Compression → tension transition
! ===================================================================
subroutine test_compression_tension_transition(props, npass, nfail)
    real(dp), intent(in) :: props(NPROPS)
    integer, intent(inout) :: npass, nfail
    real(dp) :: stress(3), strain0(3), strain1(3), statev(NSTATEV), tangent(3,3)
    real(dp) :: eps_val, sig_before_crack
    integer :: i
    logical :: has_nan, cracked
    
    write(*,'(A)') 'TEST 4: Compression → tension transition'
    call reset_state(stress, strain0, statev)
    has_nan = .false.; cracked = .false.
    
    ! Compress to -0.003
    do i = 1, 30
        strain1 = 0.0d0; strain1(1) = -0.0001d0 * i
        call do_step(props, stress, strain0, strain1, statev, tangent)
    end do
    write(*,'(A,ES12.4)') '  After compression: sig=', stress(1)
    
    ! Unload through 0 to tension (+0.001)
    do i = 1, 40
        eps_val = -0.003d0 + 0.0001d0 * i
        strain1 = 0.0d0; strain1(1) = eps_val
        call do_step(props, stress, strain0, strain1, statev, tangent)
        if (any(stress /= stress) .or. any(tangent /= tangent)) has_nan = .true.
        
        if (.not. cracked .and. statev(2) < 361.0d0) then
            cracked = .true.
            write(*,'(A,ES12.4,A,ES12.4)') '  Cracked at eps=', eps_val, ', sig=', stress(1)
        end if
    end do
    
    call check('CompTens: no NaN', .not. has_nan, npass, nfail)
    call check('CompTens: cracking occurs during transition', cracked, npass, nfail)
    call check('CompTens: final stress bounded', abs(stress(1)) < 1.0d8, npass, nfail)
    
    write(*,'(A,ES12.4)') '  Final sig=', stress(1)
    write(*,*)
end subroutine

! ===================================================================
! TEST 5: Pure shear
! ===================================================================
subroutine test_pure_shear(props, npass, nfail)
    real(dp), intent(in) :: props(NPROPS)
    integer, intent(inout) :: npass, nfail
    real(dp) :: stress(3), strain0(3), strain1(3), statev(NSTATEV), tangent(3,3)
    real(dp) :: G0, tau_max, tau_elastic
    integer :: i
    logical :: has_nan
    
    write(*,'(A)') 'TEST 5: Pure shear'
    G0 = props(1) / (2.0d0 * (1.0d0 + props(2)))
    call reset_state(stress, strain0, statev)
    has_nan = .false.; tau_max = 0.0d0
    
    ! First step: check elastic shear
    strain1 = 0.0d0; strain1(3) = 1.0d-6
    call do_step(props, stress, strain0, strain1, statev, tangent)
    tau_elastic = stress(3)
    
    do i = 2, 30
        strain1 = 0.0d0; strain1(3) = 1.0d-5 * i
        call do_step(props, stress, strain0, strain1, statev, tangent)
        if (abs(stress(3)) > abs(tau_max)) tau_max = stress(3)
        if (any(stress /= stress) .or. any(tangent /= tangent)) has_nan = .true.
    end do
    
    call check('Shear: no NaN', .not. has_nan, npass, nfail)
    call check('Shear: elastic τ ≈ G*γ', abs(tau_elastic - G0*1.0d-6) / abs(G0*1.0d-6) < 0.05d0, npass, nfail)
    call check('Shear: τ_max > 0', tau_max > 0.0d0, npass, nfail)
    
    write(*,'(A,ES12.4)') '  Max shear=', tau_max
    write(*,*)
end subroutine

! ===================================================================
! TEST 6: Equi-biaxial compression
! ===================================================================
subroutine test_biaxial_compression(props, npass, nfail)
    real(dp), intent(in) :: props(NPROPS)
    integer, intent(inout) :: npass, nfail
    real(dp) :: stress(3), strain0(3), strain1(3), statev(NSTATEV), tangent(3,3)
    real(dp) :: sig_peak
    integer :: i
    logical :: has_nan
    
    write(*,'(A)') 'TEST 6: Equi-biaxial compression'
    call reset_state(stress, strain0, statev)
    has_nan = .false.; sig_peak = 0.0d0
    
    do i = 1, 30
        strain1(1) = -0.0001d0 * i
        strain1(2) = -0.0001d0 * i
        strain1(3) = 0.0d0
        call do_step(props, stress, strain0, strain1, statev, tangent)
        if (stress(1) < sig_peak) sig_peak = stress(1)
        if (any(stress /= stress) .or. any(tangent /= tangent)) has_nan = .true.
    end do
    
    call check('Biaxial: no NaN', .not. has_nan, npass, nfail)
    call check('Biaxial: peak > uniaxial fc', sig_peak < props(4), npass, nfail)
    
    write(*,'(A,ES12.4,A,ES12.4)') '  Biaxial peak=', sig_peak, ' vs fc=', props(4)
    write(*,*)
end subroutine

! ===================================================================
! TEST 7: Direction independence (compression at 0° vs 45°)
! ===================================================================
subroutine test_direction_independence(props, npass, nfail)
    real(dp), intent(in) :: props(NPROPS)
    integer, intent(inout) :: npass, nfail
    real(dp) :: stress_0(3), stress_45(3)
    real(dp) :: strain0(3), strain1(3), statev(NSTATEV), tangent(3,3)
    real(dp) :: sig_p1_0, sig_p2_0, sig_p1_45, sig_p2_45
    real(dp) :: c, s, c2, s2, cs
    integer :: i
    
    write(*,'(A)') 'TEST 7: Direction independence (0° vs 45°)'
    
    call reset_state(stress_0, strain0, statev)
    do i = 1, 15
        strain1 = 0.0d0; strain1(1) = -0.0001d0 * i
        call do_step(props, stress_0, strain0, strain1, statev, tangent)
    end do
    sig_p1_0 = 0.5d0*(stress_0(1)+stress_0(2)) + &
               sqrt((0.5d0*(stress_0(1)-stress_0(2)))**2 + stress_0(3)**2)
    sig_p2_0 = 0.5d0*(stress_0(1)+stress_0(2)) - &
               sqrt((0.5d0*(stress_0(1)-stress_0(2)))**2 + stress_0(3)**2)
    
    call reset_state(stress_45, strain0, statev)
    c = cos(45.0d0 * DEG2RAD); s = sin(45.0d0 * DEG2RAD)
    c2 = c*c; s2 = s*s; cs = c*s
    do i = 1, 15
        strain1(1) = -0.0001d0 * i * c2
        strain1(2) = -0.0001d0 * i * s2
        strain1(3) = -0.0001d0 * i * 2.0d0 * cs
        call do_step(props, stress_45, strain0, strain1, statev, tangent)
    end do
    sig_p1_45 = 0.5d0*(stress_45(1)+stress_45(2)) + &
                sqrt((0.5d0*(stress_45(1)-stress_45(2)))**2 + stress_45(3)**2)
    sig_p2_45 = 0.5d0*(stress_45(1)+stress_45(2)) - &
                sqrt((0.5d0*(stress_45(1)-stress_45(2)))**2 + stress_45(3)**2)
    
    call check('DirIndep: P1 match (<5%)', &
        abs(sig_p1_0 - sig_p1_45) / max(abs(sig_p1_0), 1.0d3) < 0.05d0, npass, nfail)
    call check('DirIndep: P2 match (<5%)', &
        abs(sig_p2_0 - sig_p2_45) / max(abs(sig_p2_0), 1.0d3) < 0.05d0, npass, nfail)
    
    write(*,'(A,2ES12.4)') '  0°:  P1,P2=', sig_p1_0, sig_p2_0
    write(*,'(A,2ES12.4)') '  45°: P1,P2=', sig_p1_45, sig_p2_45
    write(*,*)
end subroutine

! ===================================================================
! TEST 8: Newton-Raphson tangent consistency (finite difference)
! ===================================================================
subroutine test_tangent_consistency(props, npass, nfail)
    real(dp), intent(in) :: props(NPROPS)
    integer, intent(inout) :: npass, nfail
    real(dp) :: stress_base(3), stress_pert(3)
    real(dp) :: strain0(3), strain1(3), strain1_p(3)
    real(dp) :: statev_base(NSTATEV), statev_pert(NSTATEV), tangent(3,3)
    real(dp) :: dstrain(3)
    real(dp) :: h, fd_col(3), err, max_err
    integer :: j, i
    
    write(*,'(A)') 'TEST 8: Tangent consistency (finite difference)'
    
    h = 1.0d-8
    max_err = 0.0d0
    
    ! Drive to a compressive state first
    call reset_state(stress_base, strain0, statev_base)
    do i = 1, 8
        strain1 = 0.0d0; strain1(1) = -0.0001d0 * i
        call do_step(props, stress_base, strain0, strain1, statev_base, tangent)
    end do
    
    ! At this state, do one more step and check tangent
    strain0 = strain1
    strain1 = 0.0d0; strain1(1) = -0.0009d0
    dstrain = strain1 - strain0
    
    statev_pert = statev_base
    stress_base = stress_base
    call PSUMAT(NSTATEV, NPROPS, props, stress_base, strain0, strain1, dstrain, statev_base, tangent)
    
    ! FD check: perturb each strain component
    do j = 1, 3
        statev_pert = statev_base
        stress_pert = stress_base
        
        ! Reset to pre-call state
        call reset_state(stress_pert, strain1_p, statev_pert)
        do i = 1, 8
            strain1_p = 0.0d0; strain1_p(1) = -0.0001d0 * i
            call do_step(props, stress_pert, strain1_p(1:3), strain1_p, statev_pert, tangent)
        end do
        ! Apply strain1 + perturbation
        strain1_p = strain1
        strain1_p(j) = strain1_p(j) + h
        dstrain = strain1_p - strain0
        call PSUMAT(NSTATEV, NPROPS, props, stress_pert, strain0, strain1_p, dstrain, statev_pert, tangent)
        
        fd_col = (stress_pert - stress_base) / h
        do i = 1, 3
            if (abs(tangent(i,j)) > 1.0d5) then
                err = abs(fd_col(i) - tangent(i,j)) / abs(tangent(i,j))
                if (err > max_err) max_err = err
            end if
        end do
    end do
    
    call check('TanCon: max error < 20%', max_err < 0.2d0, npass, nfail)
    write(*,'(A,F8.2,A)') '  Max tangent error: ', max_err*100.0d0, '%'
    write(*,*)
end subroutine

! ===================================================================
! TEST 9: Many-step stability (cyclic strain path)
! ===================================================================
subroutine test_many_steps_stability(props, npass, nfail)
    real(dp), intent(in) :: props(NPROPS)
    integer, intent(inout) :: npass, nfail
    real(dp) :: stress(3), strain0(3), strain1(3), statev(NSTATEV), tangent(3,3)
    real(dp) :: amp, phase
    integer :: i
    logical :: has_nan, has_huge
    
    write(*,'(A)') 'TEST 9: Stability over 200 cyclic steps'
    call reset_state(stress, strain0, statev)
    has_nan = .false.; has_huge = .false.
    
    do i = 1, 200
        phase = real(i, dp) / 200.0d0 * 4.0d0 * PI_VAL
        amp = 0.003d0
        strain1(1) = amp * sin(phase)
        strain1(2) = amp * 0.3d0 * sin(phase + 1.0d0)
        strain1(3) = amp * 0.5d0 * sin(phase + 2.0d0)
        call do_step(props, stress, strain0, strain1, statev, tangent)
        
        if (any(stress /= stress) .or. any(tangent /= tangent)) has_nan = .true.
        if (any(abs(stress) > 1.0d12)) has_huge = .true.
    end do
    
    call check('Stability: no NaN in 200 steps', .not. has_nan, npass, nfail)
    call check('Stability: no huge stress', .not. has_huge, npass, nfail)
    
    write(*,'(A,3ES12.4)') '  Final stress: ', stress
    write(*,*)
end subroutine

end program test_comprehensive
