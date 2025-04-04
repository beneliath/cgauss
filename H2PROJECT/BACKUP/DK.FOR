C**************************************************
      SUBROUTINE DK(L,K,R1LZ,AL1,R1KZ,Ak1,                      
     1              R2LZ,AL2,R2KZ,AK2,BL,BK,
     2              DKA1,DKA2,DKB,DKR1,DKR2) 
C**************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DELTA(100,100)
      DATA ONE,TWO,FOUR,SIX,EIGHT,HALF/1.,2.,4.,6.,8.,.5/               
C ***                                                              ***
C ***    SUBROUTINE DK CALCULATES FIRST DERIVATIVES OF K-FACTOR    ***
C ***                   
     
      DO 10 IL=1,100
      DO 10 IK=1,100
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
      XK3=EXP(-(ALK1*ALK2*BLK12)*DD/DETAB)
      XK=XK1*XK2*XK3
C
C     --- FIRST DERIVATIVE WITH RESPECT TO AK1
C
	XK1A1=XK1*(-(AL1+DELTA(L,K)*AK1)/ALK1 
     1       + AL1*AK1*(1.0+DELTA(L,K))/ALK1**2.0 )*((R1LZ-R1KZ)**2)
	XK3A1=XK3*(-1)
     1  * ( ((1.0+DELTA(L,K))*ALK2*BLK12*DD+ALK1*ALK2*BLK12*DDA1)/DETAB 
     2     -ALK1*ALK2*BLK12*DD*DETA1/DETAB**2.0 )
	DKA1=XK1A1*XK2*XK3 +XK3A1*XK2*XK1
C
C     --- FIRST DERIVATIVE WITH RESPECT TO Alpha_2
C
	XK2A2=XK2*(-( AL2+DELTA(L,K)*AK2)/ALK2 
     1       + AL2*AK2*(1.0+DELTA(L,K))/ALK2**2.0 )*((R2LZ-R2KZ)**2)
	XK3A2=XK3*(-1)
     1  * ( ((1.0+DELTA(L,K))*ALK1*BLK12*DD+ALK1*ALK2*BLK12*DDA2)/DETAB
     2     -ALK1*ALK2*BLK12*DD*DETA2/DETAB**2.0 )
	DKA2=XK2A2*XK1*XK3 +XK3A2*XK1*XK2      
C
C     --- FIRST DERIVATIVE WITH RESPECT TO Beta_12
C
	XK3B=XK3*(-1)
     1  * ( (1.0+DELTA(L,K))*ALK1*ALK2*DD/DETAB
     2     -ALK1*ALK2*BLK12*DD*DETB/DETAB**2.0 )
        DKB=XK3B*XK1*XK2
C
C     --- FIRST DERIVATIVE WITH RESPECT TO R1
C
      XK1R1=XK1*(-AL1*AK1/ALK1)*(R1LZ-R1KZ)*2*( DELTA(L,K)-1.0 )
	XK3R1=XK3*(-(ALK1*ALK2*BLK12)*DDR1/DETAB)
      DKR1=XK1R1*XK2*XK3 + XK3R1*XK1*XK2

C
C     --- FIRST DERIVATIVE WITH RESPECT TO R2
C
      XK2R2=XK2*(-AL2*AK2/ALK2)*(R2LZ-R2KZ)*2*( DELTA(L,K)-1.0 )
	XK3R2=XK3*(-(ALK1*ALK2*BLK12)*DDR2/DETAB)
      DKR2=XK2R2*XK1*XK3 + XK3R2*XK1*XK2


      RETURN                                                            
      END     
