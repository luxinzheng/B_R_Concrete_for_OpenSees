!===============================================================================
! Focused NR convergence debug: trace exactly what happens at step 12 of
! the cyclic loading where Test C reports NR failure.
!===============================================================================
program test_nr_debug
    use adina_concrete_mod
    implicit none
    
    integer, parameter :: nstatevs = 13, nprops = 37
    real(dp) :: props(nprops)
    real(dp) :: stress(3), s0(3), s1(3), tangent(3,3)
    real(dp) :: statev(nstatevs)
    real(dp) :: statev_save(nstatevs), stress_save(3)
    real(dp) :: eps_xx, eps_yy, deps_yy, residual
    real(dp) :: dstrain(3)
    integer  :: step, iter
    real(dp) :: dstep
    
    ! Same path as Test C
    real(dp) :: eps_path(16)
    integer  :: i
    
    call setup_props(props)
    
    dstep = 2.5d-4
    ! Steps 1-4: 0 to +1e-3
    eps_path(1) = 2.5d-4; eps_path(2) = 5d-4; eps_path(3) = 7.5d-4; eps_path(4) = 1d-3
    ! Steps 5-8: +1e-3 to 0
    eps_path(5) = 7.5d-4; eps_path(6) = 5d-4; eps_path(7) = 2.5d-4; eps_path(8) = 0d0
    ! Steps 9-12: 0 to -1e-3
    eps_path(9) = -2.5d-4; eps_path(10) = -5d-4; eps_path(11) = -7.5d-4; eps_path(12) = -1d-3
    
    statev = 0.0d0; stress = 0.0d0; s0 = 0.0d0
    eps_yy = 0.0d0
    
    write(*,'(A)') '=== NR Debug: tracing steps 1-12 ==='
    
    do step = 1, 12
        eps_xx = eps_path(step)
        
        statev_save = statev
        stress_save = stress
        
        write(*,'(A)') '----------------------------------------------------'
        write(*,'(A,I3,A,ES11.3,A,ES11.3)') 'Step ', step, &
            '  eps_xx=', eps_xx, '  initial eps_yy=', eps_yy
        write(*,'(A,3ES13.5)') '  old_stress:', stress_save(1), stress_save(2), stress_save(3)
        write(*,'(A,ES11.3,A,ES11.3,A,L2)') '  statev: eps_comp_max=', statev_save(8), &
            '  ANGLE=', statev_save(2), '  cracked=', (statev_save(2) < 361.0d0)
        
        do iter = 1, 20
            statev = statev_save
            stress = stress_save
            
            s1 = (/ eps_xx, eps_yy, 0.0d0 /)
            dstrain = s1 - s0
            call PSUMAT(nstatevs, nprops, props, stress, s0, s1, dstrain, statev, tangent)
            
            residual = stress(2)
            
            write(*,'(A,I3,A,ES11.3,A,3ES12.4,A,ES11.3,A,ES11.3)') &
                '  iter', iter, ' ey=', eps_yy, &
                ' sig=', stress(1), stress(2), stress(3), &
                ' C22=', tangent(2,2), ' ecm=', statev(8)
            
            if (abs(residual) < 1.0d0) then
                write(*,'(A,I3,A)') '  -> Converged in ', iter, ' iterations'
                exit
            end if
            
            if (abs(tangent(2,2)) < 1.0d0) then
                write(*,'(A)') '  -> SINGULAR tangent(2,2)!'
                exit
            end if
            
            deps_yy = -residual / tangent(2,2)
            eps_yy = eps_yy + deps_yy
            
            if (iter >= 10) exit
        end do
        
        s0 = s1
    end do

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

end program test_nr_debug
