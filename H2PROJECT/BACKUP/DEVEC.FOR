C     Originallly Coded:  D. Gilmore
C           Modified By:  D. Gilmore ---> 19OCT95
C     -----------------------------------------------------
C     L. ADAMOWICZ RESEARCH / UNIVERSITY OF ARIZONA / CHEM.
C     -----------------------------------------------------
      SUBROUTINE DEVEC (N,BB,A,MAXA)
C      
C     DEVEC receives the vector BB and outputs
C         the matrix A = BB
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(MAXA,MAXA),BB(5050)
C
C     A=output matrix
C     BB=input vector
C     MAXA=maximum dimension of A (MUST BE SQUARE)
C     N=square dimension of matrix (i.e. 3 if A->3x3 matrix)
C
          II=0
              DO 1000 I=1,N
              DO 1000 J=1,N
                  II=II+1
                  A(I,J)=BB(II)
 1000         CONTINUE 
C
      RETURN
      END
