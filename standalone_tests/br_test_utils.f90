
module br_test_utils
  use br_test_params, only: dp
  implicit none
contains
  subroutine assert_close(name, a, b, rtol, atol, nfail)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: a, b, rtol, atol
    integer, intent(inout) :: nfail
    real(dp) :: err, tol
    err = abs(a-b)
    tol = max(atol, rtol*max(abs(a),abs(b)))
    if (err > tol) then
      nfail = nfail + 1
      write(*,'(A,1X,ES12.4,1X,ES12.4,1X,A,1X,ES12.4)') 'FAIL:', a, b, trim(name), err
    end if
  end subroutine

  subroutine assert_finite_vec(name, x, n, nfail)
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n)
    integer, intent(inout) :: nfail
    integer :: i
    do i=1,n
      if (.not.(x(i) == x(i)) .or. abs(x(i)) > 1.0d100) then
        nfail = nfail + 1
        write(*,'(A,1X,A,1X,I0,1X,ES12.4)') 'FAIL: non-finite', trim(name), i, x(i)
      end if
    end do
  end subroutine

  subroutine write_csv_row(unit, step, eps, sig, statev, nstatevs)
    integer, intent(in) :: unit, step, nstatevs
    real(dp), intent(in) :: eps(3), sig(3), statev(nstatevs)
    write(unit,'(I0,1x,ES16.8,1x,ES16.8,1x,ES16.8,1x,ES16.8,1x,ES16.8,1x,ES16.8,1x,ES16.8,1x,ES16.8)') &
      step, eps(1), eps(2), eps(3), sig(1), sig(2), sig(3), statev(2), statev(6)
  end subroutine
end module br_test_utils
