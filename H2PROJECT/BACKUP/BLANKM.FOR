C     Originallly Coded:  D. Gilmore
C           Modified By:  D. Gilmore ---> 19OCT95
C     -----------------------------------------------------
C     L. ADAMOWICZ RESEARCH / UNIVERSITY OF ARIZONA / CHEM.
C     -----------------------------------------------------
      SUBROUTINE BLANKM (N,A,MAXA)
C      
C     BLANKM places zeros in all positions on the NxN matrix A
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(MAXA,MAXA)
C
C     A=output matrix
C     MAXA=maximum dimension of A (MUST BE SQUARE)
C     N=square dimension of matrix (i.e. 3 if A->3x3 matrix)
C
              DO 1000 I=1,N
              DO 1000 J=1,N
                  A(I,J)=0.0D+00
 1000         CONTINUE 
C
      RETURN
      END
