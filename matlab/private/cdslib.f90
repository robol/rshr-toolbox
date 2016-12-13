! CDSLIB 0.1
!
!   A Fortran library for diagonal compressed storage.
!
!   Author: Leonardo Robol <leonardo.robol@cs.kuleuven.be>
!   

subroutine cds_extract_submatrix (A, n, k, i1, i2, j1, j2, B)
  ! Implementation of a MATLAB-like slicing operation that stores
  ! the full version of A(i1:i2, j1:j2) into B.
  !
  ! INPUTS and OUTPUTS
  !
  ! A    DOUBLE PRECISION   The compressed matrix stored in CDS format. This
  !                         matrix has size n x (2*k + 1). 
  !
  ! n    INTEGER            The number of rows of A.
  !
  ! k    INTEGER            The bandwidth of the matrix A.
  !
  ! i1   INTEGER            The first row-index to extract.
  !
  ! i2   INTEGER            The last row-index to extract.
  !
  ! j1   INTEGER            The first column-index to extract.
  !
  ! j2   INTEGER            The last column-index to extract.
  !
  ! B    DOUBLE PRECISION   The matrix where the slice will be stored. It needs
  !                         to have i2 - i1 + 1 rows and j2 - j1 + 1 columns.

  implicit none

  integer :: n, k, i1, i2, j1, j2, i, j
  double precision :: A(n, -k:k), B(i1:i2, j1:j2)

  do i = i1, i2
     do j = j1, j2
        if (abs(j - i) .le. k) then
           B(i, j) = A(i, j - i)
        else
           B(i, j) = 0
        end if
     end do
  end do

end subroutine cds_extract_submatrix

subroutine cds_update_submatrix (A, n, k, i1, i2, j1, j2, B)
  ! Implementation of a MATLAB-like assignment operation that stores
  ! a sliced version of A into its compressed representation. Only the
  ! elements that would fit inside the prescribed bandwidth are read from
  ! B, all the others are ignored. 
  !
  ! INPUTS and OUTPUTS
  !
  ! A    DOUBLE PRECISION   The compressed matrix stored in CDS format. This
  !                         matrix has size n x (2*k + 1). At the end it will
  !                         be updated with the elements from B. 
  !
  ! n    INTEGER            The number of rows of A.
  !
  ! k    INTEGER            The bandwidth of the matrix A.
  !
  ! i1   INTEGER            The first row-index to update.
  !
  ! i2   INTEGER            The last row-index to update.
  !
  ! j1   INTEGER            The first column-index to update.
  !
  ! j2   INTEGER            The last column-index to update.
  !
  ! B    DOUBLE PRECISION   The matrix where the slice will be read from.
  !                         It needs to have i2 - i1 + 1 rows and
  !                         j2 - j1 + 1 columns.

  implicit none

  integer n, k, i1, i2, j1, j2, i, j
  double precision A(n, -k:k), B(i1:i2,j1:j2)

  do i = i1, i2
     do j = j1, j2
        if (abs(j - i) .le. k) then
           A(i, j - i) = B(i, j)
        end if
     end do
  end do
  
end subroutine cds_update_submatrix

subroutine cds_drot2 (A, n, k, i, c, s)
  ! Apply a rotation from both sides to the matrix A. Notice that this
  ! cannot exceed the bandwidth, so there must be enough room to accomodate
  ! the extra elements.
  !
  ! A    DOUBLE PRECISION   The compressed matrix stored in CDS format. This
  !                         matrix has size n x (2*k + 1). At the end it will
  !                         be updated using the rotation defined by c and s
  !
  ! n    INTEGER            The number of rows of A.
  !
  ! k    INTEGER            The bandwidth of the matrix A.
  !
  ! i    INTEGER            The first index where the rotation acts. 
  !
  ! c    DOUBLE PRECISION   The cosine of the rotation. 
  ! 
  ! s    DOUBLE PRECISION   The sine of the rotation.
  !
  ! work DOUBLE PRECISION   A vector of lenght 4*k

  implicit none

  integer :: n, k, i, i1, i2
  double precision :: A(n, -k:k), c, s

  i1 = max (1, i - k + 1)
  i2 = min (n, i + k)

  ! Transformation on the rows  
  call drot (i2-i1+1, A(i, i1-i), n, A(i+1, i1-i-1), n, c, s)

  ! Transformation on the columns
  call drot (i2-i1+1, A(i2,i-i2), n-1, A(i2, i+1-i2), n-1, c, s)
  
end subroutine cds_drot2

subroutine cds_drot (A, n, k, i, c, s, work)
  ! Apply a rotation from both sides to the matrix A. Notice that this
  ! cannot exceed the bandwidth, so there must be enough room to accomodate
  ! the extra elements.
  !
  ! A    DOUBLE PRECISION   The compressed matrix stored in CDS format. This
  !                         matrix has size n x (2*k + 1). At the end it will
  !                         be updated using the rotation defined by c and s
  !
  ! n    INTEGER            The number of rows of A.
  !
  ! k    INTEGER            The bandwidth of the matrix A.
  !
  ! i    INTEGER            The first index where the rotation acts. 
  !
  ! c    DOUBLE PRECISION   The cosine of the rotation. 
  ! 
  ! s    DOUBLE PRECISION   The sine of the rotation.
  !
  ! work DOUBLE PRECISION   A vector of lenght 4*k
  
  implicit none

  integer :: n, k, i, j, i1, i2, l
  double precision :: A(n, -k:k), c, s
  double precision :: work(4*k)

  ! First compute the indices that we need to extract from the
  ! banded matrix
  i1 = max (i - k + 1, 1)
  i2 = min (i + k, n)

  call cds_extract_submatrix (A, n, k, i1, i2, i, i + 1, work)
  call drot (i2 - i1 + 1, work(1), 1, work(2 + i2 - i1), 1, c, s)

  ! Apply the same rotation on the rows
  call drot (2, work(i-i1+1), i2 - i1 + 1, work(i-i1+2), i2 - i1 + 1, c, s)

  ! Store the result back
  call cds_update_submatrix (A, n, k, i1, i2, i, i + 1, work)
  call cds_update_submatrix (A, n, k, i, i, i1, i2, work(1))
  call cds_update_submatrix (A, n, k, i + 1, i + 1, i1, i2, work(2 + i2 - i1))
end subroutine cds_drot
