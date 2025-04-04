C     Originallly Coded:  D. Gilmore
C           Modified By:  D. Gilmore ---> 19OCT95
C     -----------------------------------------------------
C     L. ADAMOWICZ RESEARCH / UNIVERSITY OF ARIZONA / CHEM.
C     -----------------------------------------------------
      SUBROUTINE TRANSPOSE (N,A,MAXA,B,MAXB)
C
C     TRANSPOSE receives the matrix A and outputs
C         the matrix B as the transpose of A (A').
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(MAXA,MAXA),B(MAXB,MAXB)
C
C     A=input matrix
C     B=output matrix = A'
C     MAXA=maximum dimension of A (MUST BE SQUARE)
C     MAXB=maximum dimension of B (MUST BE SQUARE and
C         equal to MAXA
C     N=square dimension of matrix (i.e. 3 if A->3x3 matrix)
C
              DO 1000 I=1,N
              DO 1000 J=1,N
                  B(I,J)=A(J,I)
 1000         CONTINUE 
C
      RETURN
      END
