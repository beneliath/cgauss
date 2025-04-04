C************************************************************
      SUBROUTINE OVERLP(AX,AY,AZ,A1,BX,BY,BZ,A2,                      
     1                  CX,CY,CZ,A3,DX,DY,DZ,A4,
     2                                          A5,A6,OV,KLM) 
C************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C ***                                      ***
C ***     FOUR CENTER OVERLAP INTEGRAL     ***
C ***                                      *** 
      PARAMETER(ZERO=0.0D+00) 
      PI=3.1415926535898D+00
      ALK1=A1+A2                                                       
      ALK2=A3+A4                                                       
      BLK12=A5+A6  
      DETAB=ALK1*ALK2+(ALK1+ALK2)*BLK12
      IF(DETAB.LE.ZERO) GOTO 300	
C
C --- ONE CENTER OVERLAP INTEGRAL ---
C
      OV=PI**3/(DETAB*DSQRT(DETAB)) 
C
C --- CENTERS ---
C
      CALL KFACT(AX,AY,AZ,A1,BX,BY,BZ,A2,                      
     1           CX,CY,CZ,A3,DX,DY,DZ,A4,
     2                                  A5,A6,XK) 
C
C --- FOUR-CENTER OVERLAP INTEGRAL ---
C
      OV=OV*XK
C
      GOTO 1500
300   KLM=5
1500  RETURN                                                            
      END      
