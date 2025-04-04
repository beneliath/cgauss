      SUBROUTINE F(IN)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /T1/ ND,NDMAX,FFF,RR(180)
      COMMON /VALUE/ M,EMIN
      X=RR(1)
      Y=RR(2)
      FFF=(X-1174D+101)**2 + (Y+2325)**2
      WRITE(6,*)X,Y,FFF
C      F=FFF
      RETURN
      END
