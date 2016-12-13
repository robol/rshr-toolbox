program test_runtime
  IMPLICIT NONE

  DOUBLE PRECISION, ALLOCATABLE :: A(:,:), U(:,:), V(:,:), dd(:), ss(:)
  INTEGER :: n = 2048, k
  DOUBLE PRECISION :: x
  INTEGER :: i, j, lwork, INFO, II
  INTEGER(8) :: cr, t1, t2

  CALL system_clock(count_rate = cr)

  DO II = 1, 9
    k = 2**(II-1)
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
    write(*,*) k, (t2 - t1 * 1.d0) / cr

    DEALLOCATE(A, U, V, dd, ss)
  END DO

end program
