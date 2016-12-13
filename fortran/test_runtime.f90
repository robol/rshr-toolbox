program test_runtime
  IMPLICIT NONE

  DOUBLE PRECISION, ALLOCATABLE :: A(:,:), U(:,:), V(:,:), dd(:), ss(:)
  INTEGER :: n, k = 4, kk
  DOUBLE PRECISION :: x, tt(2)
  INTEGER :: i, j, lwork, INFO, II
  INTEGER(8) :: cr, t1, t2

  CALL system_clock(count_rate = cr)

  DO II = 1, 12
    n = 2**(II + 2)

    DO kk = 1, 2
      IF (kk .eq. 1) THEN
        k = 4
      ELSE
        k = 32
      END IF

    ALLOCATE(A(n, -k-1:k+1), U(n,k), V(n,k), dd(n), ss(n-1))

    DO i = 1, n
      DO j = 1, k
        CALL RANDOM_NUMBER(x)
        U(i,j) = x
        CALL RANDOM_NUMBER(x)
        V(i,j) = x
      END DO
      CALL RANDOM_NUMBER(x)
      A(i,0) = x
    END DO

    CALL system_clock(t1)
    CALL hr_impl_full_to_kh(A, n, U, V, k)
    CALL hr_impl_kh_to_h(A, n, U, V, k, dd, ss)
    CALL system_clock(t2)
 
    tt(kk) = (t2 - t1 * 1.d0) / cr   

    DEALLOCATE(A, U, V, dd, ss)
    END DO

    write(*,*) n, tt(:)

  END DO

end program
