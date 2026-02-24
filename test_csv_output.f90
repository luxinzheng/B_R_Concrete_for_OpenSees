!===============================================================================
!  Generate CSV data for stress-strain curve verification
!===============================================================================
program test_csv_output
    use adina_concrete_mod
    implicit none
    integer, parameter :: NP=37, NS=13
    real(dp) :: props(NP)
    
    call setup_props(props)
    call write_crack_cyclic_csv(props)
    call write_shear_wall_csv(props)
    call write_tension_softening_csv(props)
    write(*,'(A)') 'All CSV files generated.'

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

! Crack → close → reopen with full cycle
subroutine write_crack_cyclic_csv(props)
    real(dp), intent(in) :: props(NP)
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: statev(NS), tangent(3,3)
    integer :: i
    
    open(unit=30, file='curve_crack_cyclic.csv', status='replace')
    write(30,'(A)') 'eps_xx,sig_xx,sig_yy,ANGLE,NUMCRK_approx'
    
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    
    ! Tension to 5e-4
    do i = 1, 50
        strain1 = 0.0d0; strain1(1) = 1.0d-5 * i
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        write(30,'(ES15.7,4(",",ES15.7))') strain1(1), stress(1), stress(2), statev(2), statev(6)
    end do
    
    ! Compress to -3e-3
    do i = 1, 35
        strain1 = 0.0d0; strain1(1) = 5.0d-4 - 1.0d-4 * i
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        write(30,'(ES15.7,4(",",ES15.7))') strain1(1), stress(1), stress(2), statev(2), statev(6)
    end do
    
    ! Reopen to 1e-3
    do i = 1, 40
        strain1 = 0.0d0; strain1(1) = -3.0d-3 + 1.0d-4 * i
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        write(30,'(ES15.7,4(",",ES15.7))') strain1(1), stress(1), stress(2), statev(2), statev(6)
    end do
    
    close(30)
    write(*,'(A)') '  Done: curve_crack_cyclic.csv'
end subroutine

! Shear wall-like loading: axial compression + cyclic shear
subroutine write_shear_wall_csv(props)
    real(dp), intent(in) :: props(NP)
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: statev(NS), tangent(3,3)
    real(dp) :: t, shear_amp
    integer :: i, ntotal
    
    open(unit=31, file='curve_shear_wall.csv', status='replace')
    write(31,'(A)') 'step,eps_xx,eps_yy,gamma_xy,sig_xx,sig_yy,tau_xy,ANGLE,PGRAV'
    
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    ntotal = 600
    
    do i = 1, ntotal
        t = dble(i) / dble(ntotal) * 6.0d0 * 3.14159265d0
        shear_amp = 2.0d-3 * sin(t) * min(1.0d0, dble(i)/80.0d0)
        
        strain1(1) = 0.0d0
        strain1(2) = -3.0d-4
        strain1(3) = shear_amp
        
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        
        write(31,'(I6,",",ES15.7,6(",",ES15.7),",",F8.1)') &
            i, strain1(1), strain1(2), strain1(3), &
            stress(1), stress(2), stress(3), statev(2), statev(6)
    end do
    
    close(31)
    write(*,'(A)') '  Done: curve_shear_wall.csv'
end subroutine

! Full tension softening curve (fine resolution)
subroutine write_tension_softening_csv(props)
    real(dp), intent(in) :: props(NP)
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: statev(NS), tangent(3,3)
    integer :: i
    
    open(unit=32, file='curve_tension_detail.csv', status='replace')
    write(32,'(A)') 'eps_xx,sig_xx,sig_yy,C11,ANGLE'
    
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    
    do i = 1, 200
        strain1 = 0.0d0; strain1(1) = 5.0d-6 * i
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        write(32,'(ES15.7,3(",",ES15.7),",",F8.1)') &
            strain1(1), stress(1), stress(2), tangent(1,1), statev(2)
    end do
    
    close(32)
    write(*,'(A)') '  Done: curve_tension_detail.csv'
end subroutine

end program test_csv_output
