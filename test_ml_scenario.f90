!===============================================================================
!  Standalone test targeting Multi-layer_Shell failure scenarios
!  Tests: crack closing/reopening, NaN input, extreme strain, biaxial cycles
!  Uses same material parameters as Multi-layer_Shell test_br.tcl (Material 1)
!===============================================================================
program test_ml_scenario
    use adina_concrete_mod
    implicit none
    
    integer, parameter :: NP = 37, NS = 40
    real(dp) :: props(NP)
    integer :: npass, nfail
    
    npass = 0; nfail = 0
    call setup_ml_props(props)
    
    write(*,'(A)') '================================================================'
    write(*,'(A)') '  Multi-layer_Shell Scenario Tests (standalone PSUMAT calls)'
    write(*,'(A)') '================================================================'
    write(*,*)
    
    call test1_crack_close_reopen(props, npass, nfail)
    call test2_biaxial_shear_cycle(props, npass, nfail)
    call test3_extreme_strain(props, npass, nfail)
    call test4_nan_input(props, npass, nfail)
    call test5_wall_strain_path(props, npass, nfail)
    call test6_tangent_continuity(props, npass, nfail)
    call test7_rapid_reversal(props, npass, nfail)
    call test8_many_cycles(props, npass, nfail)
    
    write(*,*)
    write(*,'(A)') '================================================================'
    write(*,'(A,I4,A,I4,A,I4)') '  TOTAL: ', npass+nfail, &
        ' tests,  PASS=', npass, ',  FAIL=', nfail
    write(*,'(A)') '================================================================'
    if (nfail > 0) then
        write(*,'(A)') '  *** SOME TESTS FAILED ***'
        stop 1
    else
        write(*,'(A)') '  ALL TESTS PASSED'
    end if
    
contains

subroutine setup_ml_props(p)
    real(dp), intent(out) :: p(NP)
    ! Material 1 from test_br.tcl: fc=25.8MPa, E0=24.03GPa
    p(1)  = 2.40254d10   ! E0
    p(2)  = 0.2d0        ! nu
    p(3)  = 2.0d6        ! ft
    p(4)  = -2.58d7      ! fc
    p(5)  = -0.003d0     ! epsc
    p(6)  = -5.16d6      ! fu
    p(7)  = -0.021d0     ! epsu
    p(8)  = 1.0d-4       ! STIFAC
    p(9)  = 0.5d0        ! SHEFAC
    p(10) = 0.75d0; p(11) = 1.0d0; p(12) = 0.7d0; p(13) = 0.12d0
    p(14:19) = (/ 0.0d0, 0.25d0, 0.5d0, 0.75d0, 1.0d0, 1.2d0 /)
    p(20:25) = (/ 1.0d0, 1.4d0, 1.7d0, 2.2d0, 2.5d0, 2.8d0 /)
    p(26:31) = (/ 1.3d0, 1.5d0, 2.0d0, 2.3d0, 2.7d0, 3.2d0 /)
    p(32:37) = (/ 1.25d0, 1.45d0, 1.95d0, 2.25d0, 2.65d0, 3.15d0 /)
end subroutine

subroutine reset_state(sig, eps0, sv)
    real(dp), intent(out) :: sig(3), eps0(3), sv(NS)
    sig = 0.0d0; eps0 = 0.0d0; sv = 0.0d0
end subroutine

subroutine do_step(p, sig, eps0, eps1, sv, tg)
    real(dp), intent(in)    :: p(NP), eps1(3)
    real(dp), intent(inout) :: sig(3), eps0(3), sv(NS)
    real(dp), intent(out)   :: tg(3,3)
    real(dp) :: de(3)
    de = eps1 - eps0
    call PSUMAT(NS, NP, p, sig, eps0, eps1, de, sv, tg)
    eps0 = eps1
end subroutine

logical function has_nan3(v)
    real(dp), intent(in) :: v(3)
    has_nan3 = any(v /= v)
end function

logical function has_nan33(m)
    real(dp), intent(in) :: m(3,3)
    has_nan33 = any(m /= m)
end function

subroutine report(name, pass_flag, npass, nfail)
    character(*), intent(in) :: name
    logical, intent(in) :: pass_flag
    integer, intent(inout) :: npass, nfail
    if (pass_flag) then
        npass = npass + 1
        write(*,'(A,A,A)') '  [PASS] ', name, ''
    else
        nfail = nfail + 1
        write(*,'(A,A,A)') '  [FAIL] ', name, ''
    end if
end subroutine

!-----------------------------------------------------------------------
! Test 1: Crack opening → closing → reopening (1D tension-compression)
!-----------------------------------------------------------------------
subroutine test1_crack_close_reopen(p, npass, nfail)
    real(dp), intent(in) :: p(NP)
    integer, intent(inout) :: npass, nfail
    real(dp) :: sig(3), eps0(3), sv(NS), tg(3,3), eps1(3)
    real(dp) :: sig_prev(3), max_jump
    integer :: i
    logical :: ok_nan, ok_jump
    real(dp) :: E0, ft
    
    write(*,'(A)') ''
    write(*,'(A)') '--- Test 1: Crack close/reopen cycle ---'
    
    E0 = p(1); ft = p(3)
    call reset_state(sig, eps0, sv)
    ok_nan = .true.; ok_jump = .true.
    max_jump = 0.0d0
    
    ! Phase 1: Tension to 5*eps_cr (crack opens)
    do i = 1, 50
        eps1 = 0.0d0
        eps1(1) = dble(i) * 1.0d-5
        sig_prev = sig
        call do_step(p, sig, eps0, eps1, sv, tg)
        if (has_nan3(sig) .or. has_nan33(tg)) ok_nan = .false.
        if (i > 1) then
            max_jump = max(max_jump, abs(sig(1) - sig_prev(1)))
        end if
    end do
    write(*,'(A,ES12.4,A,ES12.4)') '    After tension: sig11=', sig(1), &
        '  eps11=', eps0(1)
    
    ! Phase 2: Unload to zero then compress (crack closes)
    do i = 49, -100, -1
        eps1 = 0.0d0
        eps1(1) = dble(i) * 1.0d-5
        sig_prev = sig
        call do_step(p, sig, eps0, eps1, sv, tg)
        if (has_nan3(sig) .or. has_nan33(tg)) ok_nan = .false.
        if (abs(sig(1) - sig_prev(1)) > max_jump * 50.0d0 .and. abs(sig(1) - sig_prev(1)) > 1.0d5) then
            ok_jump = .false.
            write(*,'(A,I4,A,ES12.4,A,ES12.4)') '    JUMP at i=',i, &
                ' dsig=', sig(1)-sig_prev(1), ' sig=', sig(1)
        end if
    end do
    write(*,'(A,ES12.4,A,ES12.4)') '    After compress: sig11=', sig(1), &
        '  eps11=', eps0(1)
    
    ! Phase 3: Reload tension (crack reopens)
    do i = -99, 50
        eps1 = 0.0d0
        eps1(1) = dble(i) * 1.0d-5
        sig_prev = sig
        call do_step(p, sig, eps0, eps1, sv, tg)
        if (has_nan3(sig) .or. has_nan33(tg)) ok_nan = .false.
    end do
    write(*,'(A,ES12.4,A,ES12.4)') '    After reopen:   sig11=', sig(1), &
        '  eps11=', eps0(1)
    
    call report('No NaN in crack cycle', ok_nan, npass, nfail)
    call report('No extreme stress jumps', ok_jump, npass, nfail)
end subroutine

!-----------------------------------------------------------------------
! Test 2: Biaxial shear cycle (eps12 cycling with constant compression)
!-----------------------------------------------------------------------
subroutine test2_biaxial_shear_cycle(p, npass, nfail)
    real(dp), intent(in) :: p(NP)
    integer, intent(inout) :: npass, nfail
    real(dp) :: sig(3), eps0(3), sv(NS), tg(3,3), eps1(3)
    integer :: i
    logical :: ok
    
    write(*,'(A)') ''
    write(*,'(A)') '--- Test 2: Biaxial shear cycle ---'
    
    call reset_state(sig, eps0, sv)
    ok = .true.
    
    ! Apply constant compression + cyclic shear
    do i = 1, 200
        eps1(1) = -5.0d-4
        eps1(2) = -2.0d-4
        if (i <= 50) then
            eps1(3) = dble(i) * 2.0d-5
        else if (i <= 150) then
            eps1(3) = (100.0d0 - dble(i)) * 2.0d-5
        else
            eps1(3) = (dble(i) - 200.0d0) * 2.0d-5
        end if
        call do_step(p, sig, eps0, eps1, sv, tg)
        if (has_nan3(sig) .or. has_nan33(tg)) then
            ok = .false.
            write(*,'(A,I4,A,3ES12.4)') '    NaN at step ', i, ' sig=', sig
        end if
    end do
    write(*,'(A,3ES12.4)') '    Final stress: ', sig
    call report('Biaxial shear cycle: no NaN', ok, npass, nfail)
end subroutine

!-----------------------------------------------------------------------
! Test 3: Extreme strain guard
!-----------------------------------------------------------------------
subroutine test3_extreme_strain(p, npass, nfail)
    real(dp), intent(in) :: p(NP)
    integer, intent(inout) :: npass, nfail
    real(dp) :: sig(3), eps0(3), sv(NS), tg(3,3), eps1(3)
    logical :: ok
    
    write(*,'(A)') ''
    write(*,'(A)') '--- Test 3: Extreme strain guard ---'
    
    call reset_state(sig, eps0, sv)
    ok = .true.
    
    ! Huge tension
    eps1 = (/ 0.1d0, 0.0d0, 0.0d0 /)
    call do_step(p, sig, eps0, eps1, sv, tg)
    if (has_nan3(sig) .or. has_nan33(tg)) ok = .false.
    write(*,'(A,3ES12.4)') '    10% tension: sig=', sig
    
    ! Huge compression
    call reset_state(sig, eps0, sv)
    eps1 = (/ -0.2d0, 0.0d0, 0.0d0 /)
    call do_step(p, sig, eps0, eps1, sv, tg)
    if (has_nan3(sig) .or. has_nan33(tg)) ok = .false.
    write(*,'(A,3ES12.4)') '    20% compress: sig=', sig
    
    ! Huge shear
    call reset_state(sig, eps0, sv)
    eps1 = (/ 0.0d0, 0.0d0, 0.5d0 /)
    call do_step(p, sig, eps0, eps1, sv, tg)
    if (has_nan3(sig) .or. has_nan33(tg)) ok = .false.
    write(*,'(A,3ES12.4)') '    50% shear: sig=', sig
    
    call report('Extreme strain: no NaN', ok, npass, nfail)
end subroutine

!-----------------------------------------------------------------------
! Test 4: NaN input handling
!-----------------------------------------------------------------------
subroutine test4_nan_input(p, npass, nfail)
    real(dp), intent(in) :: p(NP)
    integer, intent(inout) :: npass, nfail
    real(dp) :: sig(3), eps0(3), sv(NS), tg(3,3), eps1(3), de(3)
    real(dp) :: nan_val
    logical :: ok
    
    write(*,'(A)') ''
    write(*,'(A)') '--- Test 4: NaN input handling ---'
    ok = .true.
    
    ! Create NaN
    nan_val = 0.0d0
    nan_val = nan_val / nan_val
    
    ! NaN in strain1
    call reset_state(sig, eps0, sv)
    eps1 = (/ nan_val, 0.0d0, 0.0d0 /)
    de = eps1 - eps0
    call PSUMAT(NS, NP, p, sig, eps0, eps1, de, sv, tg)
    if (has_nan3(sig) .or. has_nan33(tg)) then
        ok = .false.
        write(*,'(A)') '    FAIL: NaN in strain1 → NaN output'
    else
        write(*,'(A,3ES12.4)') '    NaN strain1 → sig=', sig
    end if
    
    ! NaN in strain0 (simulating corrupted committed state)
    sig = 0.0d0; sv = 0.0d0
    eps0 = (/ nan_val, 0.0d0, 0.0d0 /)
    eps1 = (/ 1.0d-4, 0.0d0, 0.0d0 /)
    de = eps1 - eps0
    call PSUMAT(NS, NP, p, sig, eps0, eps1, de, sv, tg)
    if (has_nan3(sig) .or. has_nan33(tg)) then
        ok = .false.
        write(*,'(A)') '    FAIL: NaN in strain0 → NaN output'
    else
        write(*,'(A,3ES12.4)') '    NaN strain0 → sig=', sig
    end if
    
    ! NaN in stress (corrupted committed stress)
    eps0 = 0.0d0
    sig = (/ nan_val, 0.0d0, 0.0d0 /)
    sv = 0.0d0
    eps1 = (/ 1.0d-4, 0.0d0, 0.0d0 /)
    de = eps1 - eps0
    call PSUMAT(NS, NP, p, sig, eps0, eps1, de, sv, tg)
    if (has_nan3(sig) .or. has_nan33(tg)) then
        ok = .false.
        write(*,'(A)') '    FAIL: NaN in stress → NaN output'
    else
        write(*,'(A,3ES12.4)') '    NaN stress → sig=', sig
    end if
    
    call report('NaN input: clean output', ok, npass, nfail)
end subroutine

!-----------------------------------------------------------------------
! Test 5: Wall-like strain path (biaxial with reversal)
! Simulates a Gauss point in the Multi-layer_Shell under cyclic shear
!-----------------------------------------------------------------------
subroutine test5_wall_strain_path(p, npass, nfail)
    real(dp), intent(in) :: p(NP)
    integer, intent(inout) :: npass, nfail
    real(dp) :: sig(3), eps0(3), sv(NS), tg(3,3), eps1(3)
    integer :: i, nsteps
    logical :: ok
    real(dp) :: frac, amp
    
    write(*,'(A)') ''
    write(*,'(A)') '--- Test 5: Wall-like biaxial strain path ---'
    
    call reset_state(sig, eps0, sv)
    ok = .true.
    amp = 3.0d-3
    nsteps = 500
    
    do i = 1, nsteps
        frac = dble(i) / dble(nsteps)
        ! Cyclic: 0→+amp→-amp→0 
        if (frac < 0.25d0) then
            ! Loading right: tension left, compression right + shear
            eps1(1) = frac * 4.0d0 * amp          ! tension
            eps1(2) = -frac * 4.0d0 * amp * 0.3d0 ! compression (Poisson-like)
            eps1(3) = frac * 4.0d0 * amp * 0.5d0  ! shear
        else if (frac < 0.75d0) then
            ! Reversal: from +amp to -amp
            eps1(1) = (0.5d0 - frac) * 4.0d0 * amp
            eps1(2) = -(0.5d0 - frac) * 4.0d0 * amp * 0.3d0
            eps1(3) = (0.5d0 - frac) * 4.0d0 * amp * 0.5d0
        else
            ! Return to zero
            eps1(1) = (frac - 1.0d0) * 4.0d0 * amp
            eps1(2) = -(frac - 1.0d0) * 4.0d0 * amp * 0.3d0
            eps1(3) = (frac - 1.0d0) * 4.0d0 * amp * 0.5d0
        end if
        
        call do_step(p, sig, eps0, eps1, sv, tg)
        if (has_nan3(sig) .or. has_nan33(tg)) then
            ok = .false.
            write(*,'(A,I4,A,3ES11.3)') '    NaN at step ', i, &
                ' eps=', eps1
            exit
        end if
        ! Check tangent diagonal is positive
        if (tg(1,1) < 0.0d0 .or. tg(2,2) < 0.0d0 .or. tg(3,3) < 0.0d0) then
            write(*,'(A,I4,A,3ES11.3)') '    NEG tangent at step ', i, &
                ' diag=', tg(1,1), tg(2,2), tg(3,3)
        end if
    end do
    
    if (ok) write(*,'(A)') '    Completed 500 steps without NaN'
    call report('Wall strain path: no NaN', ok, npass, nfail)
end subroutine

!-----------------------------------------------------------------------
! Test 6: Tangent continuity during crack closing
! Fine-grained strain steps across the closing threshold
!-----------------------------------------------------------------------
subroutine test6_tangent_continuity(p, npass, nfail)
    real(dp), intent(in) :: p(NP)
    integer, intent(inout) :: npass, nfail
    real(dp) :: sig(3), eps0(3), sv(NS), tg(3,3), eps1(3)
    real(dp) :: tg_prev(3,3), max_ratio
    integer :: i
    logical :: ok_nan, ok_cont
    real(dp) :: E0, ft, eps_cr
    
    write(*,'(A)') ''
    write(*,'(A)') '--- Test 6: Tangent continuity across crack closing ---'
    
    E0 = p(1); ft = p(3)
    eps_cr = ft / E0
    
    call reset_state(sig, eps0, sv)
    ok_nan = .true.; ok_cont = .true.
    max_ratio = 1.0d0
    
    ! First: crack the material
    eps1 = (/ 3.0d0 * eps_cr, 0.0d0, 0.0d0 /)
    call do_step(p, sig, eps0, eps1, sv, tg)
    tg_prev = tg
    
    ! Now: fine steps from 3*eps_cr down to -3*eps_cr (through closing)
    do i = 1, 600
        eps1(1) = (3.0d0 - dble(i) * 0.01d0) * eps_cr
        eps1(2) = 0.0d0; eps1(3) = 0.0d0
        call do_step(p, sig, eps0, eps1, sv, tg)
        if (has_nan3(sig) .or. has_nan33(tg)) then
            ok_nan = .false.
            write(*,'(A,I4)') '    NaN at fine step ', i
        end if
        ! Check tangent(1,1) ratio change
        if (tg_prev(1,1) > 1.0d3 .and. tg(1,1) > 1.0d3) then
            if (tg(1,1) / tg_prev(1,1) > max_ratio) max_ratio = tg(1,1) / tg_prev(1,1)
            if (tg_prev(1,1) / tg(1,1) > max_ratio) max_ratio = tg_prev(1,1) / tg(1,1)
        end if
        tg_prev = tg
    end do
    
    write(*,'(A,ES12.4)') '    Max tangent ratio between steps: ', max_ratio
    if (max_ratio > 1000.0d0) then
        ok_cont = .false.
        write(*,'(A)') '    WARNING: tangent ratio > 1000x (discontinuity)'
    end if
    
    call report('Tangent: no NaN during closing', ok_nan, npass, nfail)
    call report('Tangent: continuity < 1000x', ok_cont, npass, nfail)
end subroutine

!-----------------------------------------------------------------------
! Test 7: Rapid reversal (simulates failed N-R iteration)
! After cracking, apply extreme strains in opposite direction
!-----------------------------------------------------------------------
subroutine test7_rapid_reversal(p, npass, nfail)
    real(dp), intent(in) :: p(NP)
    integer, intent(inout) :: npass, nfail
    real(dp) :: sig(3), eps0(3), sv(NS), tg(3,3), eps1(3)
    real(dp) :: sig_saved(3), eps0_saved(3), sv_saved(NS)
    integer :: i
    logical :: ok
    real(dp) :: E0, ft, eps_cr
    
    write(*,'(A)') ''
    write(*,'(A)') '--- Test 7: Rapid reversal (N-R overshoot simulation) ---'
    
    E0 = p(1); ft = p(3)
    eps_cr = ft / E0
    call reset_state(sig, eps0, sv)
    ok = .true.
    
    ! Load to cracking (5 steps)
    do i = 1, 5
        eps1 = (/ dble(i) * eps_cr, 0.0d0, 0.0d0 /)
        call do_step(p, sig, eps0, eps1, sv, tg)
    end do
    write(*,'(A,ES12.4,A,ES12.4)') '    After crack: sig11=', sig(1), &
        '  tg11=', tg(1,1)
    
    ! Save committed state
    sig_saved = sig; eps0_saved = eps0; sv_saved = sv
    
    ! Simulate N-R overshoots: trial with wildly different strains
    ! (These represent trial strains that would NOT be committed)
    do i = 1, 10
        sig = sig_saved; eps0 = eps0_saved; sv = sv_saved
        eps1(1) = -dble(i) * 0.005d0  ! increasingly compressive
        eps1(2) = dble(i) * 0.002d0
        eps1(3) = dble(i) * 0.003d0
        call do_step(p, sig, eps0, eps1, sv, tg)
        if (has_nan3(sig) .or. has_nan33(tg)) then
            ok = .false.
            write(*,'(A,I2,A,ES10.3,A,3ES11.3)') '    NaN at overshoot ', i, &
                ' eps11=', eps1(1), ' sig=', sig
        end if
    end do
    
    ! After all overshoots, do one more call from the committed state
    ! with a small reversal increment (this is what happens after wipeAnalysis)
    sig = sig_saved; eps0 = eps0_saved; sv = sv_saved
    eps1 = (/ 4.0d0*eps_cr, 0.0d0, 0.0d0 /)
    call do_step(p, sig, eps0, eps1, sv, tg)
    if (has_nan3(sig) .or. has_nan33(tg)) then
        ok = .false.
        write(*,'(A)') '    NaN in recovery step after overshoots'
    else
        write(*,'(A,ES12.4,A,ES12.4)') '    Recovery: sig11=', sig(1), &
            '  tg11=', tg(1,1)
    end if
    
    call report('Rapid reversal: no NaN', ok, npass, nfail)
end subroutine

!-----------------------------------------------------------------------
! Test 8: Many full cycles (endurance test)
!-----------------------------------------------------------------------
subroutine test8_many_cycles(p, npass, nfail)
    real(dp), intent(in) :: p(NP)
    integer, intent(inout) :: npass, nfail
    real(dp) :: sig(3), eps0(3), sv(NS), tg(3,3), eps1(3)
    integer :: cycle, i, nan_count
    real(dp) :: amp, phase
    
    write(*,'(A)') ''
    write(*,'(A)') '--- Test 8: Endurance test (20 full cycles) ---'
    
    call reset_state(sig, eps0, sv)
    nan_count = 0
    
    do cycle = 1, 20
        amp = dble(cycle) * 1.5d-4
        do i = 1, 100
            phase = dble(i) / 100.0d0 * 2.0d0 * 3.14159265d0
            eps1(1) = amp * sin(phase)
            eps1(2) = -0.3d0 * amp * sin(phase)
            eps1(3) = 0.5d0 * amp * cos(phase)
            call do_step(p, sig, eps0, eps1, sv, tg)
            if (has_nan3(sig) .or. has_nan33(tg)) nan_count = nan_count + 1
        end do
    end do
    
    write(*,'(A,I4,A)') '    Completed 2000 steps, NaN count = ', nan_count, ''
    write(*,'(A,3ES12.4)') '    Final stress: ', sig
    call report('Endurance: no NaN in 2000 steps', nan_count == 0, npass, nfail)
end subroutine

end program test_ml_scenario
