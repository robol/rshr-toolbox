program hr
  IMPLICIT NONE

  INTEGER :: n
  DOUBLE PRECISION, ALLOCATABLE :: A(:,:), WORK(:), tau(:)
  DOUBLE PRECISION :: x
  INTEGER :: i, j, lwork, INFO, II
  INTEGER(8) cr, t1, t2

  DO II = 1, 12
    n = 2**(II + 2)
    lwork = n

  ALLOCATE(A(N,N), WORK(LWORK), tau(n-1))

  DO i = 1, n
    DO j = 1, n
      CALL RANDOM_NUMBER(x)
      A(i,j) = x
    END DO
  END DO

  CALL system_clock (count_rate = cr)
  CALL system_clock (t1)

  CALL DGEHRD(n, 1, n, A, n, tau, work, lwork, INFO)

  CALL system_clock (t2)

  WRITE(*,*) n, (t2 - t1 * 1.d0) / cr

  DEALLOCATE(A, work, tau)
END DO

end program
