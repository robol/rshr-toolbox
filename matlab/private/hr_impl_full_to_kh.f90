SUBROUTINE hr_apply_rotation_to_hermitian (n, k, A, i, j, c, s)
! Apply a rotation to the Hermitian matrix A on the rows i and j. 
! Since the matrix is Hermitian, we only need to update one of them 
! and the other one only on the 2 x 2 diagonal block and then
! retrieve the result by symmetry. 
!
! PARAMETERS 
!
! n      (input) The dimension of the square matrix A
! lo     (input) An integer such that A(l,i) = A(l,j) = 0
!        for every l <= lo. This is used to make the cost of
!        a rotation O(hi - lo) < O(n). 
! hi     (input) The same of lo, but this time it guaratees that 
!        A(l,i) = A(l,j) = 0 for every l >= hi. 
! A      (input/output) The matrix A to be updated with the rotations. 
! i      (input) The first index on which the rotation acts. 
! j      (input) The second index on which the rotation acts. 
! c      (input) The cosine of the angle of the rotation. 
! s      (input) The sine of the angle of the rotation. 
!
! Author: Leonardo Robol <leonardo.robol@sns.it>

  INTEGER n, i, j
  DOUBLE PRECISION A(n,-k:k), c, s

  CALL CDS_DROT2 (A, n, k, i, c, s)
  
END SUBROUTINE hr_apply_rotation_to_hermitian

FUNCTION HR_DOT_PRODUCT (n, x, y, inc) RESULT(dp)
  DOUBLE PRECISION :: dp
  INTEGER :: n, inc, j
  DOUBLE PRECISION :: x(*), y(*)

  dp = 0.d0
  DO j = 1, n * inc, inc
     dp = dp + x(j) * y(j)
  END DO
END FUNCTION HR_DOT_PRODUCT

SUBROUTINE hr_impl_full_to_kh(A, n, U, V, k)
! Perform the reduction of the matrix A + U*V' to k-Hessenberg form
! where A is diagonal and U and V are n x k matrices. 
!
! The matrices A, U, V are overwritten with the matrices
! QAQ', QU and QV such that QAQ' has bandwidth k. 
!
! In particular, this implies that QAQ' + QUV'Q' is in
! k-Hessenberg form. 
!
! PARAMETERS 
!
! A      (input/output) The diagonal matrix A that needs to 
!        be reduced to banded form. On output, the banded matrix
!        QAQ' is stored here. 
! n      (input) The size of the square matrix A. 
! U      (input/output) The first term of the low rank correction. 
!        On output, the updated upper triangular term QU is
!        stored here. 
! V      (input/output) The second term of the low rank correction. 
!        On output, the updated term QV is stored here. 
! k      (input) The number columns of U and V. 
!
! Author: Leonardo Robol <leonardo.robol@sns.it>

  INTEGER :: n, k, i, j, s, l, jj
  DOUBLE PRECISION :: A(n,-k-1:k+1), U(n,k), V(n,k), Gc, Gs, da, db

  DO i = n, 2, -1
     DO j = 1, min(k, n-i+1)
        l = i + j - 1;

        da = U(l-1,j)
        db = U(l,j)
        CALL DROTG (da, db, Gc, Gs)

        CALL DROT (k, U(l-1,1), n, U(l,1), n, Gc, Gs)

        CALL DROT (k, V(l-1,1), n, V(l,1), n, Gc, Gs)
        
        CALL HR_APPLY_ROTATION_TO_HERMITIAN (n, k + 1, A, l-1, l, Gc, Gs)
     END DO

     ! Perform Bulge chasing if needed     
     DO s = i + k, n

        ! G = planerot (A(s-1:s,s-k-1));
        da = A(s-1,-k) ! A(s-1,s-k-1)
        db = A(s,-k-1)  ! A(s,s-k-1)
        CALL DROTG (da, db, Gc, Gs)

        CALL DROT (k, V(s-1,1), n, V(s,1), n, Gc, Gs)
        
        CALL HR_APPLY_ROTATION_TO_HERMITIAN (n, k + 1, A, s-1, s, Gc, Gs)
        
     END DO
  END DO

END SUBROUTINE HR_IMPL_FULL_TO_KH

SUBROUTINE HR_IMPL_KH_TO_H(A, n, U, V, k, dd, ss)
! Reduce a k-Hessenberg matrix represented as A + U*V'
! where A is banded with bandwidth k and U and V are n x k
! matrices with U = [ X ; 0 ], X k x k to standard Hessenberg
! form. 
!
! The low rank factors U and V are updated in place
! and the n and n-1 vectors dd and ss are used to store
! the diagonal and subdiagonal elements of the Hessenberg form. 
!
! PARAMETERS 
!
! A      (input/output) The banded matrix A. This storage is used 
!        in the computation and is meaningless at the end. 
! n      (input) The size of the square matrix A. 
! U      (input/output) The first term of the low rank correction. 
!        On output, the updated upper triangular term QU is
!        stored here. 
! V      (input/output) The second term of the low rank correction. 
!        On output, the updated term QV is stored here. 
! k      (input) The number columns of U and V. 
! dd     (output) dd The diagonal elements of the Hessenberg form.
! ss     (output) the subdiagonal elements of the Hessenberg form.    
!
! Author: Leonardo Robol <leonardo.robol@sns.it>

  INTEGER :: n, k, i, j, s, l, jj
  DOUBLE PRECISION :: A(n,-k-1:k+1), U(n,k), V(n,k), Gc, Gs
  DOUBLE PRECISION :: da, db, dd(n), ss(n-1), HR_DOT_PRODUCT

  DO j = 1, n - 2
     DO i = min(n, j + k), j + 2, -1
    
        ! G = planerot (A(i-1:i,j) + U(i-1:i,:) * V(j,:)');
        da = A(i-1,j-i+1) + HR_DOT_PRODUCT(k, U(i-1,1), V(j,1), n)   ! The entry in A was (i-1,j)
        db = A(i,j-i) + HR_DOT_PRODUCT(k, U(i,1), V(j,1), n)       ! The entry in A was (i,j)
        CALL DROTG (da, db, Gc, Gs)
        
        ! U(i-1:i,:) = G * U(i-1:i,:);
        CALL DROT (k, U(i-1,1), n, U(i,1), n, Gc, Gs)

        ! V(i-1:i,:) = G * V(i-1:i,:);
        CALL DROT (k, V(i-1,1), n, V(i,1), n, Gc, Gs)
        
        CALL HR_APPLY_ROTATION_TO_HERMITIAN (n, k + 1, A, i-1, i, Gc, Gs)

        DO s = i + k, n, k
           lo = max (1, s - k - 1)
           hi = min (n, s+k+2)
           
           da = A(s-1,-k) ! A(s-1,s-k-1)
           db = A(s,-k-1)   ! A(s,s-k-1)

           CALL DROTG (da, db, Gc, Gs)
           
           CALL DROT (k, V(s-1,1), n, V(s,1), n, Gc, Gs)

           CALL HR_APPLY_ROTATION_TO_HERMITIAN (n, k + 1, A, s-1, s, Gc, Gs)
           
        END DO

        dd(j) = A(j,0) ! A(j,j)
        ss(j) = A(j+1,-1) ! A(j+1,j)

        DO jj = 1, k
           dd(j) = dd(j) + U(j,jj) * V(j,jj)
           ss(j) = ss(j) + U(j+1,jj) * V(j,jj)
        END DO
     END DO
  END DO

  dd(n-1) = A(n-1,0) + HR_DOT_PRODUCT (k, U(n-1,1), V(n-1,1), n)
  ss(n-1) = A(n,-1) + HR_DOT_PRODUCT (k, U(n,1), V(n-1,1), n)
  dd(n)   = A(n,0) + HR_DOT_PRODUCT (k, U(n,1), V(n,1), n)

END SUBROUTINE HR_IMPL_KH_TO_H
