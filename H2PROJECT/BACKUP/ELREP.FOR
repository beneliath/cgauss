


C**********************************************
      SUBROUTINE ELREP(AX,AY,AZ,A1,BX,BY,BZ,A2,                      
     1                 CX,CY,CZ,A3,DX,DY,DZ,A4,
     2                 A5,A6,ER,klm)                       
C**********************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C ***                                                 ***
C ***     FOUR CENTER ELECTRON REPULSION INTEGRAL     ***
C ***                                                 ***
      DATA ONE,TWO,FOUR,SIX,EIGHT,HALF/1.,2.,4.,6.,8.,.5/               
      PI=3.1415926535898D+00
      PI52=PI**2.5D+00
      ALK1=A1+A2                                                       
      ALK2=A3+A4                                                       
      BLK12=A5+A6  
      DETAB=ALK1*ALK2+(ALK1+ALK2)*BLK12
C
C --- ONE CENTERELECTRON REPULSION INTEGRAL ---
C
      ER=TWO*PI52/(DETAB*DSQRT(ALK1+ALK2))
C
C --- CENTERS ---
C
      CALL KFACT(AX,AY,AZ,A1,BX,BY,BZ,A2,                      
     1           CX,CY,CZ,A3,DX,DY,DZ,A4,
     2                               A5,A6,XK) 
C
      RX1=(A1*AX+A2*BX)/ALK1
      RX2=(A3*CX+A4*DX)/ALK2
      RY1=(A1*AY+A2*BY)/ALK1
      RY2=(A3*CY+A4*DY)/ALK2
      RZ1=(A1*AZ+A2*BZ)/ALK1
      RZ2=(A3*CZ+A4*DZ)/ALK2
      DD=(RX1-RX2)**2+(RY1-RY2)**2+(RZ1-RZ2)**2
      ZZ=EXP(-(ALK1*ALK2*BLK12)*DD/DETAB)
C
C --- FOUR-CENTER ELECTRON REPULSION INTEGRAL ---
C
      WW=ALK1*ALK1*ALK2*ALK2*DD/((ALK1+ALK2)*DETAB) 
      ER=ER*XK*F0(WW) 
C
      GOTO 1500
300   KLM=5      
1500  RETURN                                                            
      END 
