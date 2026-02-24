!===============================================================================
!  Realistic wall element test: combined normal + shear strain paths
!  In the FEM, all strain components change simultaneously during cyclic push
!===============================================================================
program test_wall_realistic
    use adina_concrete_mod
    implicit none
    integer, parameter :: NP=37, NS=13
    real(dp) :: props(NP), stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: statev(NS), tangent(3,3)
    real(dp) :: peak_pos(3), peak_neg(3), t
    integer :: i, nsteps, cycle
    real(dp) :: axial_pre, disp_max, dt
    
    call setup_props(props)
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    peak_pos = 0.0d0; peak_neg = 0.0d0
    
    write(*,'(A)') '================================================================'
    write(*,'(A)') '  REALISTIC WALL ELEMENT: combined axial+shear strain path'
    write(*,'(A)') '================================================================'
    
    axial_pre = -2.0d-4
    nsteps = 20
    
    ! Phase 0: Axial precompression
    write(*,'(A)') '--- Phase 0: Axial precompression ---'
    do i = 1, 5
        strain1 = (/0.0d0, axial_pre/5.0d0*i, 0.0d0/)
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
    end do
    write(*,'(A,3ES12.4)') '  After precompression: ', stress
    
    ! Cyclic push protocol (3 half-cycles)
    disp_max = 8.0d-4
    
    do cycle = 1, 3
        write(*,'(A,I2,A)') '--- Cycle ', cycle, ' ---'
        
        if (mod(cycle,2) == 1) then
            ! Forward push: eps_xx increases (tension), eps_yy constant, gamma increases
            do i = 1, nsteps
                t = dble(i) / dble(nsteps)
                strain1(1) = 2.0d-4 * t * cycle
                strain1(2) = axial_pre - 0.5d-4 * t
                strain1(3) = disp_max * t * cycle
                dstrain = strain1 - strain0
                call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
                strain0 = strain1
                if (stress(3) > peak_pos(3)) peak_pos = stress
                if (stress(3) < peak_neg(3)) peak_neg = stress
            end do
            write(*,'(A,3ES12.4)') '  At peak:   ', stress
            write(*,'(A,F8.1,A,F5.1)') '  ANGLE=', statev(2), '  PGRAV=', statev(6)
        else
            ! Reverse push: eps_xx decreases (compression), gamma reverses
            do i = 1, nsteps
                t = dble(i) / dble(nsteps)
                strain1(1) = 2.0d-4*1 * (1.0d0-t) + (-2.0d-4) * t * (cycle-1)
                strain1(2) = axial_pre - 0.5d-4 + 0.5d-4 * t
                strain1(3) = disp_max*1 * (1.0d0-t) + (-disp_max) * t * (cycle-1)
                dstrain = strain1 - strain0
                call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
                strain0 = strain1
                if (stress(3) > peak_pos(3)) peak_pos = stress
                if (stress(3) < peak_neg(3)) peak_neg = stress
            end do
            write(*,'(A,3ES12.4)') '  At peak:   ', stress
            write(*,'(A,F8.1,A,F5.1)') '  ANGLE=', statev(2), '  PGRAV=', statev(6)
        end if
    end do
    
    write(*,*)
    write(*,'(A,ES12.4)') '  Peak +tau:  ', peak_pos(3)
    write(*,'(A,ES12.4)') '  Peak -tau:  ', peak_neg(3)
    if (abs(peak_pos(3)) > 1.0d-10) then
        write(*,'(A,F6.1,A)') '  Asymmetry |neg/pos|: ', abs(peak_neg(3)/peak_pos(3))*100.0d0, '%'
    end if
    
    ! Check tangent PD
    write(*,*)
    write(*,'(A)') '  Final tangent diagonal:'
    write(*,'(A,3ES12.4)') '    C11,C22,C33: ', tangent(1,1), tangent(2,2), tangent(3,3)
    
    if (tangent(1,1) > 0.0d0 .and. tangent(2,2) > 0.0d0 .and. tangent(3,3) > 0.0d0) then
        write(*,'(A)') '  Tangent: POSITIVE diagonal [OK]'
    else
        write(*,'(A)') '  Tangent: NEGATIVE diagonal [WARN]'
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
