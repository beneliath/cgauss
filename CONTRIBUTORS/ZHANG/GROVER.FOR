C**************************************************
      SUBROUTINE grover(L,K,R1LZ,AL1,R1KZ,Ak1,                      
     1              R2LZ,AL2,R2KZ,AK2,BL,BK,
     2              DSA1,DSA2,DSB,DSR1,DSR2) 
C**************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DELTA(100,100)
      DATA ONE,TWO,FOUR,SIX,EIGHT,HALF/1.,2.,4.,6.,8.,.5/               
      PI=3.1415926535898D+00

C ***                                                              ***
C ***    SUBROUTINE  CALCULATES FIRST DERIVATIVES OF overlap       ***
C ***                   
     
      DO 10 iL=1,100
      DO 10 iK=1,100
      DELTA(IL,IK)=0.0
10    CONTINUE
      DO 11 II=1,100
      DELTA(II,II)=1.0 
11    CONTINUE
      ALK1=AL1+Ak1                                                       
      ALK2=AL2+Ak2                                                       
      BLK12=BL+BK  
      DETAB=ALK1*ALK2+(ALK1+ALK2)*BLK12
	DETA1=(1.0+DELTA(L,K))*(ALK2 +BLK12)
	DETA2=(1.0+DELTA(L,K))*(ALK1 +BLK12)
	DETB=(1.0+DELTA(L,K))*(ALK1 +ALK2)
C
      OVTEM=PI**3/(DETAB*DSQRT(DETAB)) 
C     --- CENTERS ---
C
      XK1=EXP(-((AL1*AK1)/ALK1)*((R1LZ-R1KZ)**2))
      XK2=EXP(-((AL2*AK2)/ALK2)*((R2LZ-R2KZ)**2))
C
      RZ1=(AL1*R1LZ+Ak1*R1KZ)/ALK1
      RZ2=(AL2*R2LZ+Ak2*R2KZ)/ALK2
      DD=(RZ1-RZ2)**2
      XK3=EXP(-(ALK1*ALK2*BLK12)*DD/DETAB)
      XK=XK1*XK2*XK3
C
      OV=OVTEM*XK

      CALL DK(L,K,R1LZ,AL1,R1KZ,Ak1,                      
     1              R2LZ,AL2,R2KZ,AK2,BL,BK,
     2              DKA1,DKA2,DKB,DKR1,DKR2) 
C
C
      DSA1=(-1.5)*(OVTEM/DETAB)*DETA1*XK + OVTEM*DKA1
      DSA2=(-1.5)*(OVTEM/DETAB)*DETA2*XK + OVTEM*DKA2
      DSB =(-1.5)*(OVTEM/DETAB)*DETB *XK + OVTEM*DKB
      DSR1=OVTEM*DKR1
      DSR2=OVTEM*DKR2


      RETURN
      END

