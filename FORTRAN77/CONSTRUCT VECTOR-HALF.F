SUBROUTINE VECH (N,A,MAXA,BB)

    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	IMPLICIT INTEGER (I-N)

    DIMENSION A(MAXA,MAXA),BB(MAXA*MAXA)

!     A=input matrix
!     AA=output vector
!     MAXA=maximum dimension of A (MUST BE SQUARE)
!     N=square dimension of matrix (i.e. 3 if A->3x3 matrix)

          II=0

            DO I=1,N
				DO J=1,I

                  II=II+1
                  BB(II)=A(I,J)

				END DO
			END DO


END SUBROUTINE VECH
