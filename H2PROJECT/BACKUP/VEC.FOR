C     Originallly Coded:  D. Gilmore
C           Modified By:  D. Gilmore ---> 19OCT95
C     -----------------------------------------------------
C     L. ADAMOWICZ RESEARCH / UNIVERSITY OF ARIZONA / CHEM.
C     -----------------------------------------------------
      SUBROUTINE VEC (N,A,MAXA,AA)
C      
C     MATRIXCOPY receives the matrix A and outputs
C         the matrix B = A
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(MAXA,MAXA),AA(5050)
C
C     A=input matrix
C     AA=output vector
C     MAXA=maximum dimension of A (MUST BE SQUARE)
C     N=square dimension of matrix (i.e. 3 if A->3x3 matrix)
C
          II=0
              DO 1000 I=1,N
              DO 1000 J=1,N
                  II=II+1
                  AA(II)=A(I,J)
 1000         CONTINUE 
C
      RETURN
      END
