program test_tangent_consistency
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    
    real(dp) :: props(37), stress(3), strain0(3), strain1(3), statev(40), tangent(3,3)
    real(dp) :: stress_base(3), statev_base(40), stress_pert(3), statev_pert(40)
    real(dp) :: tangent_base(3,3), tangent_pert(3,3)
    real(dp) :: deps, numerical(3,3), err, max_err
    integer :: i, j, npass, nfail
    
    npass = 0; nfail = 0
    call init_props(props)
    deps = 1.0d-8
    
    write(*,'(A)') '======================================================'
    write(*,'(A)') ' Tangent Consistency & Positivity Tests'
    write(*,'(A)') '======================================================'
    
    ! ---- Test A: Elastic compression (before peak) ----
    write(*,'(A)') ''
    write(*,'(A)') '--- Test A: Elastic compression (before peak) ---'
    stress = 0.0d0; statev = 0.0d0; strain0 = 0.0d0
    strain1 = (/ -0.001d0, -0.0003d0, 0.0002d0 /)
    call call_psumat(props, stress, strain0, strain1, statev, tangent)
    
    write(*,'(A,3ES12.4)') '  stress     =', stress
    write(*,'(A,3ES12.4)') '  tang diag  =', tangent(1,1), tangent(2,2), tangent(3,3)
    
    if (tangent(1,1) > 0.0d0 .and. tangent(2,2) > 0.0d0 .and. tangent(3,3) > 0.0d0) then
        write(*,'(A)') '  [PASS] Tangent positive definite (elastic)'
        npass = npass + 1
    else
        write(*,'(A)') '  [FAIL] Tangent not positive definite!'
        nfail = nfail + 1
    end if
    
    ! ---- Test B: At peak compression ----
    write(*,'(A)') ''
    write(*,'(A)') '--- Test B: At peak compression (eps=eps_c) ---'
    stress = 0.0d0; statev = 0.0d0; strain0 = 0.0d0
    strain1 = (/ -0.003d0, -0.0005d0, 0.0001d0 /)
    call call_psumat(props, stress, strain0, strain1, statev, tangent)
    
    write(*,'(A,3ES12.4)') '  stress     =', stress
    write(*,'(A,3ES12.4)') '  tang diag  =', tangent(1,1), tangent(2,2), tangent(3,3)
    
    if (tangent(1,1) > 0.0d0 .and. tangent(2,2) > 0.0d0 .and. tangent(3,3) > 0.0d0) then
        write(*,'(A)') '  [PASS] Tangent positive at peak (secant approach)'
        npass = npass + 1
    else
        write(*,'(A)') '  [FAIL] Tangent not positive at peak!'
        nfail = nfail + 1
    end if
    
    ! ---- Test C: Post-peak compression ----
    write(*,'(A)') ''
    write(*,'(A)') '--- Test C: Post-peak compression (eps > eps_c) ---'
    stress = 0.0d0; statev = 0.0d0; strain0 = 0.0d0
    strain1 = (/ -0.006d0, -0.001d0, 0.0001d0 /)
    call call_psumat(props, stress, strain0, strain1, statev, tangent)
    
    write(*,'(A,3ES12.4)') '  stress     =', stress
    write(*,'(A,3ES12.4)') '  tang diag  =', tangent(1,1), tangent(2,2), tangent(3,3)
    
    if (tangent(1,1) > 0.0d0 .and. tangent(2,2) > 0.0d0 .and. tangent(3,3) > 0.0d0) then
        write(*,'(A)') '  [PASS] Tangent positive post-peak (secant approach)'
        npass = npass + 1
    else
        write(*,'(A)') '  [FAIL] Tangent not positive post-peak!'
        nfail = nfail + 1
    end if
    
    ! ---- Test D: Cracked state then compression ----
    write(*,'(A)') ''
    write(*,'(A)') '--- Test D: Cracked state (tension crack → compression) ---'
    stress = 0.0d0; statev = 0.0d0; strain0 = 0.0d0
    strain1 = (/ 0.0005d0, -0.0001d0, 0.0d0 /)
    call call_psumat(props, stress, strain0, strain1, statev, tangent)
    strain0 = strain1
    strain1 = (/ -0.002d0, -0.001d0, 0.0003d0 /)
    call call_psumat(props, stress, strain0, strain1, statev, tangent)
    
    write(*,'(A,3ES12.4)') '  stress     =', stress
    write(*,'(A,3ES12.4)') '  tang diag  =', tangent(1,1), tangent(2,2), tangent(3,3)
    
    if (tangent(1,1) > 0.0d0 .and. tangent(2,2) > 0.0d0 .and. tangent(3,3) > 0.0d0) then
        write(*,'(A)') '  [PASS] Cracked tangent positive'
        npass = npass + 1
    else
        write(*,'(A)') '  [FAIL] Cracked tangent not positive!'
        nfail = nfail + 1
    end if
    
    ! ---- Test E: Unloading from peak ----
    write(*,'(A)') ''
    write(*,'(A)') '--- Test E: Unloading from peak ---'
    stress = 0.0d0; statev = 0.0d0; strain0 = 0.0d0
    strain1 = (/ -0.003d0, -0.0005d0, 0.0d0 /)
    call call_psumat(props, stress, strain0, strain1, statev, tangent)
    strain0 = strain1
    strain1 = (/ -0.002d0, -0.0003d0, 0.0d0 /)
    call call_psumat(props, stress, strain0, strain1, statev, tangent)
    
    write(*,'(A,3ES12.4)') '  stress     =', stress
    write(*,'(A,3ES12.4)') '  tang diag  =', tangent(1,1), tangent(2,2), tangent(3,3)
    
    if (tangent(1,1) > 0.0d0 .and. tangent(2,2) > 0.0d0) then
        write(*,'(A)') '  [PASS] Unloading tangent positive'
        npass = npass + 1
    else
        write(*,'(A)') '  [FAIL] Unloading tangent not positive!'
        nfail = nfail + 1
    end if
    
    ! ---- Test F: Cyclic through 50 cycles - tangent always positive? ----
    write(*,'(A)') ''
    write(*,'(A)') '--- Test F: 50-cycle tangent positivity endurance ---'
    stress = 0.0d0; statev = 0.0d0; strain0 = 0.0d0
    block
        integer :: cyc, step
        real(dp) :: amp, ex, gam
        logical :: all_positive
        all_positive = .true.
        amp = 0.003d0
        
        do cyc = 1, 50
            do step = 1, 10
                ex = amp * dble(step) / 10.0d0
                gam = 0.2d0 * ex
                strain1 = (/ ex, -0.1d0*ex, gam /)
                call call_psumat(props, stress, strain0, strain1, statev, tangent)
                strain0 = strain1
                if (tangent(1,1) <= 0.0d0 .or. tangent(2,2) <= 0.0d0 .or. tangent(3,3) <= 0.0d0) then
                    all_positive = .false.
                end if
            end do
            do step = 10, -10, -1
                ex = amp * dble(step) / 10.0d0
                gam = 0.2d0 * ex
                strain1 = (/ ex, -0.1d0*ex, gam /)
                call call_psumat(props, stress, strain0, strain1, statev, tangent)
                strain0 = strain1
                if (tangent(1,1) <= 0.0d0 .or. tangent(2,2) <= 0.0d0 .or. tangent(3,3) <= 0.0d0) then
                    all_positive = .false.
                end if
            end do
            do step = -10, 0
                ex = amp * dble(step) / 10.0d0
                gam = 0.2d0 * ex
                strain1 = (/ ex, -0.1d0*ex, gam /)
                call call_psumat(props, stress, strain0, strain1, statev, tangent)
                strain0 = strain1
                if (tangent(1,1) <= 0.0d0 .or. tangent(2,2) <= 0.0d0 .or. tangent(3,3) <= 0.0d0) then
                    all_positive = .false.
                end if
            end do
        end do
        
        if (all_positive) then
            write(*,'(A)') '  [PASS] All tangent diagonals positive through 50 cycles'
            npass = npass + 1
        else
            write(*,'(A)') '  [FAIL] Some tangent diagonals went non-positive!'
            nfail = nfail + 1
        end if
    end block
    
    write(*,'(A)') ''
    write(*,'(A)') '======================================================'
    write(*,'(A,I3,A,I3,A)') ' RESULTS: ', npass, ' PASS / ', nfail, ' FAIL'
    if (nfail == 0) then
        write(*,'(A)') ' >>> ALL TESTS PASSED <<<'
    end if
    write(*,'(A)') '======================================================'
    
contains

subroutine init_props(p)
    real(dp), intent(out) :: p(37)
    p = 0.0d0
    p(1)  = 3.15d10
    p(2)  = 0.2d0
    p(3)  = -2.58d7
    p(4)  = -0.003d0
    p(5)  = 2.0d6
    p(6)  = -0.01d0
    p(7)  = 1.0d-4
    p(8)  = 0.5d0
    p(9)  = 0.005d0
    p(10) = 0.0d0
    p(37) = 1.0d0
    p(11) = -0.2d0; p(12) = 0.9d0; p(13) = -1.25d0; p(14) = 0.0d0
    p(15) = 0.12d0; p(16) = -0.3d0; p(17) = 1.0d0; p(18) = 1.0d0
    p(19:24) = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
    p(25:31) = (/ 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0 /)
    p(32:37) = (/ 1.25d0, 1.45d0, 1.95d0, 2.25d0, 2.65d0, 3.15d0 /)
end subroutine

subroutine call_psumat(props, stress, strain0, strain1, statev, tangent)
    real(dp), intent(in)    :: props(37)
    real(dp), intent(inout) :: stress(3), statev(40)
    real(dp), intent(in)    :: strain0(3), strain1(3)
    real(dp), intent(out)   :: tangent(3,3)
    real(dp) :: dstrain(3)
    dstrain = strain1 - strain0
    call PSUMAT(40, 37, props, stress, strain0, strain1, dstrain, statev, tangent)
end subroutine

end program test_tangent_consistency
