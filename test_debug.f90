program test_debug
    use adina_concrete_mod
    implicit none
    integer, parameter :: NP=37, NS=13
    real(dp) :: props(NP), stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: statev(NS), tangent(3,3)
    real(dp) :: eps_c, eps_r, eps_min
    integer :: i
    
    props = 0.0d0
    props(1) = 21.4d9; props(2) = 0.2d0; props(3) = 2.07d6
    props(4) = -20.7d6; props(5) = -0.002d0; props(6) = -4.14d6
    props(7) = -0.006d0; props(8) = 1.0d-4; props(9) = 0.5d0
    props(10) = 0.75d0; props(11) = 1.0d0; props(12) = 0.7d0; props(13) = 0.12d0
    props(14:19) = (/ 0.0d0, 0.25d0, 0.5d0, 0.75d0, 1.0d0, 1.2d0 /)
    props(20:25) = (/ 1.0d0, 1.1d0,  1.2d0, 1.15d0, 1.0d0, 0.9d0 /)
    props(26:31) = (/ 1.3d0, 1.5d0,  2.0d0, 2.3d0,  2.7d0, 3.2d0 /)
    props(32:37) = (/ 1.25d0, 1.45d0, 1.95d0, 2.25d0, 2.65d0, 3.15d0 /)
    
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    
    write(*,'(A)') 'Step  eps_xx        sig_xx        sig_yy        ANGLE     PGRAV  eps_comp_max  C11'
    
    do i = 1, 50
        strain1 = 0.0d0
        strain1(1) = -0.0001d0 * i
        dstrain = strain1 - strain0
        
        call PSUMAT(NS, NP, props, stress, strain0, strain1, dstrain, statev, tangent)
        strain0 = strain1
        
        write(*,'(I4,2X,ES12.4,2X,ES12.4,2X,ES12.4,2X,F8.1,2X,F5.1,2X,ES12.4,2X,ES12.4)') &
            i, strain1(1), stress(1), stress(2), statev(2), statev(6), statev(8), tangent(1,1)
    end do
end program
