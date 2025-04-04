C**************************************************
      SUBROUTINE grele(L,K,R1LZ,AL1,R1KZ,Ak1,                      
     1              R2LZ,AL2,R2KZ,AK2,BL,BK,
     2              DERA1,DERA2,DERB,DERR1,DERR2) 
C**************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DELTA(100,100)
      DATA ONE,TWO,FOUR,SIX,EIGHT,HALF/1.,2.,4.,6.,8.,.5/               
      PI=3.1415926535898D+00

C ***                                                              ***
C ***    SUBROUTINE  CALCULATES FIRST DERIVATIVES OF V_ee          ***
C ***                   
     
      DO 10 iL=1,100
      DO 10 iK=1,100
      DELTA(iL,iK)=0.0
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
c
      G1=2.0*(PI**2.5)/(DETAB*DSQRT(ALK1+ALK2))
      DG1A1=-2.0*(PI**2.5)*( DETA1*DSQRT(ALK1+ALK2) 
     1      +DETAB*0.5*(1.0+DELTA(L,K))/DSQRT(ALK1+ALK2) )
     2      /(DETAB*DSQRT(ALK1+ALK2))**2.0
      DG1A2=-2.0*(PI**2.5)*( DETA2*DSQRT(ALK1+ALK2) 
     1      +DETAB*0.5*(1.0+DELTA(L,K))/DSQRT(ALK1+ALK2) )
     2      /(DETAB*DSQRT(ALK1+ALK2))**2.0
      DG1B=-2.0*(PI**2.5)/DSQRT(ALK1+ALK2)*( DETB/DETAB**2.0 )

C     --- CENTERS ---
C
      XK1=EXP(-((AL1*AK1)/ALK1)*((R1LZ-R1KZ)**2))
      XK2=EXP(-((AL2*AK2)/ALK2)*((R2LZ-R2KZ)**2))
C
      RZ1=(AL1*R1LZ+Ak1*R1KZ)/ALK1
      RZ1A1= ( DELTA(L,K)*R1LZ+R1KZ )/ALK1
     1     -( AL1*R1LZ+Ak1*R1KZ )*(1.0+DELTA(L,K))/ALK1**2.0
      RZ1R1=(AL1*DELTA(L,K) + AK1 )/ALK1	
C
      RZ2=(AL2*R2LZ+Ak2*R2KZ)/ALK2
      RZ2A2=( DELTA(L,K)*R2LZ+R2KZ )/ALK2
     1     -( AL2*R2LZ+Ak2*R2KZ )*(1.0+DELTA(L,K))/ALK2**2.0
      RZ2R2=(AL2*DELTA(L,K) + AK2 )/ALK2	
C
      DD=(RZ1-RZ2)**2
      DDA1=  2*(RZ1-RZ2)*RZ1A1
      DDA2=- 2*(RZ1-RZ2)*RZ2A2
      DDR1= 2*(RZ1-RZ2)*RZ1R1	
      DDR2=- 2*(RZ1-RZ2)*RZ2R2
C
      XK3=EXP(-(ALK1*ALK2*BLK12)*DD/DETAB)
      XK=XK1*XK2*XK3
      CALL DK(L,K,R1LZ,AL1,R1KZ,Ak1,                      
     1              R2LZ,AL2,R2KZ,AK2,BL,BK,
     2              DKA1,DKA2,DKB,DKR1,DKR2) 
C
      WWN=ALK1*ALK1*ALK2*ALK2*DD
      WWNA1=2.0*ALK1*(1.0+DELTA(L,K))*ALK2*ALK2*DD
     1     +ALK1*ALK1*ALK2*ALK2*DDA1
      WWNA2=2.0*ALK2*(1.0+DELTA(L,K))*ALK1*ALK1*DD
     1     +ALK1*ALK1*ALK2*ALK2*DDA2
C
      WWD=(ALK1+ALK2)*DETAB
      WWDA1=(1.0+DELTA(L,K))*DETAB +(ALK1+ALK2)*DETA1
      WWDA2=(1.0+DELTA(L,K))*DETAB +(ALK1+ALK2)*DETA2
      WW=WWN/WWD
C     THE DERIVATIVES OF WW
      WWA1=(WWNA1*WWD-WWDA1*WWN)/WWD**2.0
      WWA2=(WWNA2*WWD-WWDA2*WWN)/WWD**2.0
      WWB =ALK1*ALK1*ALK2*ALK2*DD/(ALK1+ALK2) *(-DETB/DETAB**2.0 ) 
      WWR1=ALK1*ALK1*ALK2*ALK2*DDR1/((ALK1+ALK2)*DETAB) 
      WWR2=ALK1*ALK1*ALK2*ALK2*DDR2/((ALK1+ALK2)*DETAB) 

C      THE DERIVATIVES OF F0
c      F0A1=-F1(WW) *WWA1
c      F0A2=-F1(WW) *WWA2
c      F0B =-F1(WW) *WWB
c      F0R1=-F1(WW) *WWR1
c      F0R2=-F1(WW) *WWR2


C -------- THIS IS A TEST-------------------------
      ER=G1*XK*WW 
      DERA1=DG1A1*XK*WW + G1*DKA1*WW +G1*XK*WWA1
      DERA2=DG1A2*XK*WW + G1*DKA2*WW +G1*XK*WWA2
      DERB =DG1B *XK*WW + G1*DKB *WW +G1*XK*WWB
      DERR1=G1*DKR1*WW +G1*XK*WWR1
      DERR2=G1*DKR2*WW +G1*XK*WWR2
C ------- END OF TEST-------------- 


c      ER=G1*XK*F0(WW) 
c      DERA1=DG1A1*XK*F0(WW) + G1*DKA1*F0(WW) +G1*XK*F0A1
c      DERA2=DG1A2*XK*F0(WW) + G1*DKA2*F0(WW) +G1*XK*F0A2
c      DERB =DG1B *XK*F0(WW) + G1*DKB *F0(WW) +G1*XK*F0B
c      DERR1=G1*DKR1*F0(WW) +G1*XK*F0R1
c      DERR2=G1*DKR2*F0(WW) +G1*XK*F0R2
	RETURN
        END
