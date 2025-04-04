SUBROUTINE TAB(AX,N,M,NNN,MMM)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	  IMPLICIT INTEGER (I-N)

      DIMENSION AX(NNN,MMM)

!     --- PRINT A MATRIX STORED IN RECTANGULAR FORM

      MM=M/10
      IF(MM.EQ.0) GO TO 6
      DO II=1,MM
          JP=(II-1)*10+1
          JK=II*10
          WRITE(7,11)
          WRITE(6,12)(I,I=JP,JK)
!             ...This comment was left by Zhang

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


!--- OutPut FORMAT(s) ---
!-----------------------------------------------------

    1 FORMAT(1X,I3,3X,20F23.15)
   11 FORMAT(/)
   12 FORMAT(6X,10(10X,I3,11X))


END SUBROUTINE TAB
