!===============================================================================
!  Debug Newton iteration for uniaxial tension
!===============================================================================
program test_newton
    use adina_concrete_mod, only: dp
    implicit none
    
    integer, parameter :: NPROPS = 37, NSTATEVS = 13, MAX_ITER = 50
    real(dp), parameter :: TOL = 1.0d-10
    real(dp) :: props(NPROPS), statev(NSTATEVS)
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3), tangent(3,3)
    real(dp) :: stress_t(3), statev_t(NSTATEVS)
    real(dp) :: eps_xx, eps_yy, deps, residual, d_eps_yy
    integer  :: i, iter
    
    call setup_props(props)
    statev  = 0.0d0
    stress  = 0.0d0
    strain0 = 0.0d0
    
    deps = 5.0d-6
    
    do i = 1, 25
        eps_xx = deps * real(i, dp)
        eps_yy = -props(2) * eps_xx  ! initial guess
        
        ! Newton iteration
        do iter = 1, MAX_ITER
            strain1(1) = eps_xx;  strain1(2) = eps_yy;  strain1(3) = 0.0d0
            dstrain = strain1 - strain0
            stress_t = stress
            statev_t = statev
            call PSUMAT(NSTATEVS, NPROPS, props, stress_t, strain0, strain1, dstrain, statev_t, tangent)
            residual = stress_t(2)
            
            if (i >= 19 .and. i <= 23) then
                write(*,'(A,I3,A,I2,A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)') &
                    'Step ', i, ' iter ', iter, &
                    ': eps_yy=', eps_yy, ' sig_xx=', stress_t(1), &
                    ' sig_yy=', stress_t(2), ' T22=', tangent(2,2)
            end if
            
            if (abs(residual) < TOL) exit
            if (abs(tangent(2,2)) < 1.0d-30) then
                write(*,*) 'WARNING: tangent(2,2) near zero!'
                exit
            end if
            d_eps_yy = -residual / tangent(2,2)
            eps_yy = eps_yy + d_eps_yy
        end do
        
        ! Final commit call
        strain1(1) = eps_xx;  strain1(2) = eps_yy;  strain1(3) = 0.0d0
        dstrain = strain1 - strain0
        call PSUMAT(NSTATEVS, NPROPS, props, stress, strain0, strain1, dstrain, statev, tangent)
        strain0 = strain1
        
        write(*,'(A,I3,A,ES12.4,A,ES12.4,A,ES12.4,A,F10.2)') &
            '>>> Step ', i, ': eps_xx=', eps_xx, ' sig_xx=', stress(1), &
            ' sig_yy=', stress(2), ' ANGLE=', statev(2)
    end do

contains
subroutine setup_props(props)
    real(dp), intent(out) :: props(NPROPS)
    props = 0.0d0
    props(1)  = 30000.0d0; props(2)  = 0.2d0
    props(3)  = 3.0d0;     props(4)  = -30.0d0
    props(5)  = -0.002d0;  props(6)  = -6.0d0
    props(7)  = -0.005d0;  props(8)  = 1.0d-4
    props(9)  = 0.5d0;     props(10) = 0.75d0
    props(11) = 1.0d0;     props(12) = 0.7d0;  props(13) = 0.12d0
    props(14) = 0.0d0;  props(15) = 0.25d0; props(16) = 0.5d0
    props(17) = 0.75d0; props(18) = 1.0d0;  props(19) = 1.2d0
    props(20) = 1.0d0;  props(21) = 1.4d0;  props(22) = 1.7d0
    props(23) = 2.2d0;  props(24) = 2.5d0;  props(25) = 2.8d0
    props(26) = 1.3d0;  props(27) = 1.5d0;  props(28) = 2.0d0
    props(29) = 2.3d0;  props(30) = 2.7d0;  props(31) = 3.2d0
    props(32) = 1.25d0; props(33) = 1.45d0; props(34) = 1.95d0
    props(35) = 2.25d0; props(36) = 2.65d0; props(37) = 3.15d0
end subroutine
end program
