SUBROUTINE BLANKM (N,A,MAXA)
      
!BLANKM places zeros in all positions on the NxN matrix A

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	  IMPLICIT INTEGER (I-N)

      DIMENSION A(MAXA,MAXA)

!     A=output matrix
!     MAXA=maximum dimension of A (MUST BE SQUARE)
!     N=square dimension of matrix (i.e. 3 if A->3x3 matrix)

        DO I=1,N
            DO J=1,N

                  A(I,J)=0.0D+00

			END DO
		END DO

END SUBROUTINE BLANKM

