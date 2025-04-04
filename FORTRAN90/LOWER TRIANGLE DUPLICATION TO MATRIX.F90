SUBROUTINE DEVEC (N,BB,A,MAXA)
      
!DEVEC receives the vector BB and outputs the matrix A = BB

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	  IMPLICIT INTEGER (I-N)

      DIMENSION A(MAXA,MAXA),BB(MAXA*MAXA)

!     A=output matrix
!     BB=input vector
!    MAXA=maximum dimension of A (MUST BE SQUARE)
!     N=square dimension of matrix (i.e. 3 if A->3x3 matrix)

          II=0

		DO I=1,N
            DO J=1,N
        
		          II=II+1
                  A(I,J)=BB(II)
		
			END DO
		END DO

END SUBROUTINE DEVEC
