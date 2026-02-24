!===============================================================================
!  NR Transition Test: verify stress continuity at crack open/close boundary
!  Simulates what happens during NR iterations when strain crosses the
!  crack-closing threshold back and forth.
!===============================================================================
program test_nr_transition
    use adina_concrete_mod
    implicit none
    integer, parameter :: NP=37, NS=13
    real(dp) :: props(NP), stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: statev(NS), tangent(3,3), statev_save(NS), stress_save(3)
    real(dp) :: sig_open(3), sig_closed(3), jump
    integer :: i
    logical :: pass_all
    
    call setup_props(props)
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    pass_all = .true.
    
    write(*,'(A)') '================================================================'
    write(*,'(A)') '  NR TRANSITION TEST: stress continuity at crack closing'
    write(*,'(A)') '================================================================'
    
    ! Phase 1: Create a crack via tension
    do i = 1, 5
        strain1 = (/2.0d-5*i, 0.0d0, 0.0d0/)
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
    end do
    write(*,'(A,3ES12.4,A,F8.1)') '  After cracking: ', stress, '  ANGLE=', statev(2)
    
    ! Phase 2: Compress past the crack (crack closes)
    do i = 1, 10
        strain1 = (/1.0d-4 - 3.0d-5*i, 0.0d0, 0.0d0/)
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
    end do
    write(*,'(A,3ES12.4,A,F8.1)') '  After closing:  ', stress, '  ANGLE=', statev(2)
    
    ! Save this converged state (crack is closed)
    stress_save = stress; statev_save = statev
    
    ! Phase 3: Simulate NR iterations around the closing threshold
    ! Try strains that are just above and below the crack-closing threshold
    write(*,*)
    write(*,'(A)') '  NR iteration simulation: strains near crack-closing threshold'
    write(*,'(A)') '  eps_xx        sig_xx       sig_yy       tau_xy       ANGLE  jump'
    
    do i = -5, 5
        strain1(1) = statev_save(3) + 1.0d-6 * i  ! CRKSTR(1) ± small offset
        strain1(2) = 0.0d0; strain1(3) = 0.0d0
        dstrain = strain1 - strain0
        
        ! Call with the SAVED converged state (like NR would do)
        stress = stress_save; statev = statev_save
        call PSUMAT(NS,NP,props,stress,(/strain0(1),strain0(2),strain0(3)/), &
                    strain1,dstrain,statev,tangent)
        
        if (i == -1) sig_closed = stress
        if (i == 1)  sig_open = stress
        
        write(*,'(ES12.4,3ES13.4,F8.1)') strain1(1), stress(1), stress(2), stress(3), statev(2)
    end do
    
    ! Check jump at threshold
    jump = maxval(abs(sig_open - sig_closed))
    write(*,*)
    write(*,'(A,ES12.4)') '  Max stress jump at threshold: ', jump
    write(*,'(A,ES12.4)') '  ft = ', props(3)
    
    if (jump < props(3) * 0.5d0) then
        write(*,'(A)') '  PASS: Jump < 0.5*ft (NR should handle this)'
    else if (jump < props(3) * 2.0d0) then
        write(*,'(A)') '  WARN: Jump is moderate (NR might need extra iterations)'
    else
        write(*,'(A)') '  FAIL: Jump is too large (NR likely to fail)'
        pass_all = .false.
    end if
    
    ! Phase 4: Multi-cycle stress recovery test
    write(*,*)
    write(*,'(A)') '  Multi-cycle recovery (5 tension-compression cycles):'
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    
    do i = 1, 5
        ! Tension to cracking
        strain1 = (/2.0d-4, 0.0d0, 0.0d0/)
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        ! Compression
        strain1 = (/-5.0d-4, 0.0d0, 0.0d0/)
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        write(*,'(A,I2,A,ES12.4,A,F6.1)') '    Cycle ', i, ': sig_xx=', stress(1), &
            '  ANGLE=', statev(2)
    end do
    
    ! Check that stress converges toward total-strain value
    ! Virgin material at eps=-5e-4: sig = E/(1-v^2) * (-5e-4) ≈ -11.1 MPa
    write(*,'(A,ES12.4)') '  Final sig_xx: ', stress(1)
    write(*,'(A,ES12.4)') '  Virgin at eps=-5e-4: ', -5.0d-4*props(1)/(1.0d0-props(2)**2)
    
    write(*,*)
    if (pass_all) then
        write(*,'(A)') '  ALL NR TRANSITION TESTS PASSED'
    else
        write(*,'(A)') '  SOME NR TRANSITION TESTS FAILED'
    end if

contains
subroutine setup_props(props)
    real(dp), intent(out) :: props(NP)
    props = 0.0d0
    props(1)=21.4d9; props(2)=0.2d0; props(3)=2.07d6; props(4)=-20.7d6
    props(5)=-0.002d0; props(6)=-4.14d6; props(7)=-0.006d0
    props(8)=1.0d-4; props(9)=0.5d0; props(10)=0.75d0; props(11)=1.0d0
    props(12)=0.7d0; props(13)=0.12d0
    props(14:19)=(/0d0,0.25d0,0.5d0,0.75d0,1d0,1.2d0/)
    props(20:25)=(/1d0,1.1d0,1.2d0,1.15d0,1d0,0.9d0/)
    props(26:31)=(/1.3d0,1.5d0,2d0,2.3d0,2.7d0,3.2d0/)
    props(32:37)=(/1.25d0,1.45d0,1.95d0,2.25d0,2.65d0,3.15d0/)
end subroutine
end program
