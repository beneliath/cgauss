C*********************
C    PROGRAM MAIN
C*********************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C--MMAX - MAX. LENGTH OF EXPANSION--
C
      PARAMATER (MMAX=50)
      DIMENSION
     1   AA1(MMAX),RX1(MMAX),RY1(MMAX),RZ1(MMAX),
     2   AA2(MMAX),RX2(MMAX),RY2(MMAX),RZ2(MMAX),
     3   B12(MMAX)
      DIMENSION S(MMAX,MMAX),H(MMAX,MMAX)CO(MMAX),G(5*MMAX)
    
        COMMON /VALUE/ EMIN
       
C--DATA--
C
      CALL GUESS(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12)
    
    
C
      CALL HAMILTON(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,S,H,CO)
    
C
      CALL GRAD(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,SH,CO,G)
    
    
    
      STOP
      END
