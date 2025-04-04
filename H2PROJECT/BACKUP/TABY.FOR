C-----------------------------------
      SUBROUTINE TABY(AX,N,M,NNN,MMM)
C-----------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)   
      DIMENSION AX(NNN,MMM)
C
C     --- PRINT A MATRIX STORED IN RECTANGULAR FORM
C
      MM=M/10
      IF(MM.EQ.0) GO TO 6
      DO II=1,MM
          JP=(II-1)*10+1
          JK=II*10
          WRITE(7,11)
          WRITE(6,12)(I,I=JP,JK)
C             ...This comment was left by Zhang
          DO I=1,N
              WRITE(7,1)I,(AX(I,J),J=JP,JK)
              WRITE(6,1)I,(AX(I,J),J=JP,JK)
          END DO
      END DO
6     CONTINUE 
      MA=MM*10+1
      IF(MA.GT.M) RETURN
      WRITE(7,12)(I,I=MA,M) 
      WRITE(6,12)(I,I=MA,M)
      DO I=1,N
          WRITE(7,1) I,(AX(I,J),J=MA,M)
          WRITE(6,1) I,(AX(I,J),J=MA,M)
      END DO
C
C
C     -----------------------------------------------------
C     --- OutPut FORMAT(s) ---
C     -----------------------------------------------------
    1 FORMAT(1X,I3,3X,*)
   11 FORMAT(/)
   12 FORMAT(6X,10(10X,I3,11X))
C     -----------------------------------------------------
C     -----------------------------------------------------
C
C
      RETURN
      END
