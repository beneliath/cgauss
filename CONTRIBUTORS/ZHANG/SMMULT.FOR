C     Originallly Coded:  D. Gilmore
C           Modified By:  D. Gilmore ---> 19OCT95
C     -----------------------------------------------------
C     L. ADAMOWICZ RESEARCH / UNIVERSITY OF ARIZONA / CHEM.
C     -----------------------------------------------------
      SUBROUTINE SMMULT (N,A,B,MAXB,C,MAXC)
C
C     SMMULT multiplies scalar A into the matrix B and
C         passes the result back as matrix C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION B(MAXB,MAXB),C(MAXC,MAXC)
C
C     A=input scalar
C     B=input matrix
C     MAXB=maximum dimension of B (MUST BE SQUARE)
C     MAXC=maximum dimension of C (MUST BE SQUARE and
C         equal to MAXB
C     N=square dimension of matrix (i.e. 3 if A->3x3 matrix)
C
              DO 1000 I=1,N
              DO 1000 J=1,N
                  C(I,J)=A*B(I,J)
 1000         CONTINUE 
C
      RETURN
      END
