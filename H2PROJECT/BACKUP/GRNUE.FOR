C**************************************************
      SUBROUTINE grNUe(L,K,R1LZ,AL1,R1KZ,Ak1,                      
     1              R2LZ,AL2,R2KZ,AK2,BL,BK,EZ,
     2              DENA1,DENA2,DENB,DENR1,DENR2) 
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
c
      G1=2.0*PI**2.5*(
     1 (1.0/(DETAB*DSQRT(ALK2+BLK12))) +
     2 (1.0/(DETAB*DSQRT(ALK1+BLK12))) )
      G1A1=2.0*PI**2.5*(
     1 -( DETA1*DSQRT(ALK1+BLK12) 
     2   +DETAB*0.5*(1.0+DELTA(L,K))/DSQRT(ALK1+BLK12) )
     3 /(DETAB**2.0*(ALK1+BLK12) ) 
     4 - DETA1/(DETAB**2.0*DSQRT(ALK2+BLK12))   )
C
      G1A2=2.0*PI**2.5*(
     1 -( DETA2*DSQRT(ALK2+BLK12) 
     2   +DETAB*0.5*(1.0+DELTA(L,K))/DSQRT(ALK2+BLK12) )
     3 /(DETAB**2.0*(ALK2+BLK12) ) 
     4 - DETA2/(DSQRT(ALK1+BLK12)*DETAB**2.0)   )
C
       G1B=2.0*PI**2.5*(
     1 -( DETB*DSQRT(ALK1+BLK12) 
     2   +DETAB*0.5*(1.0+DELTA(L,K))/DSQRT(ALK1+BLK12) )
     3 /(DETAB**2.0*(ALK1+BLK12) ) 
     4 -( DETB*DSQRT(ALK2+BLK12) 
     5   +DETAB*0.5*(1.0+DELTA(L,K))/DSQRT(ALK2+BLK12) )
     6 /(DETAB**2.0*(ALK2+BLK12) ) )
    
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
C      RZ1=(AL1*R1LZ+Ak1*R1KZ)/ALK1
C      RZ2=(AL2*R2LZ+Ak2*R2KZ)/ALK2
C
      XK3=EXP(-(ALK1*ALK2*BLK12)*DD/DETAB)
      XK=XK1*XK2*XK3
      CALL DK(L,K,R1LZ,AL1,R1KZ,Ak1,                      
     1              R2LZ,AL2,R2KZ,AK2,BL,BK,
     2              DKA1,DKA2,DKB,DKR1,DKR2) 
C
C     THIS IS FOR f0 FUNCTION
C
C ***  FOR THE FIRST ELECTRON
      XX=ALK1*(ALK2+BLK12)*RZ1+ALK2*BLK12*RZ2
      XXA1=(1.0+DELTA(L,K))*(ALK2+BLK12)*RZ1 +ALK1*(ALK2+BLK12)*RZ1A1
      XXA2=ALK1*(1.0+DELTA(L,K))*RZ1+(1.0+DELTA(L,K))*BLK12*RZ2
     1    +ALK2*BLK12*RZ2A2
      XXB=ALK1*(1.0+DELTA(L,K))*RZ1+ALK2*(1.0+DELTA(L,K))*RZ2
      XXR1=ALK1*(ALK2+BLK12)*RZ1R1
      XXR2=ALK2*BLK12*RZ2R2
C
      RRZ=XX/DETAB 
      RRZA1=(XXA1*DETAB-XX*DETA1)/DETAB**2.0
      RRZA2=(XXA2*DETAB-XX*DETA2)/DETAB**2.0     
      RRZB=(XXB*DETAB-XX*DETB)/DETAB**2.0
      RRZR1=XXR1/DETAB 
      RRZR2=XXR2/DETAB 
C
      RENO=(RRZ-EZ)**2
      RENOA1=2.0*(RRZ-EZ)*RRZA1
      RENOA2=2.0*(RRZ-EZ)*RRZA2
      RENOB =2.0*(RRZ-EZ)*RRZB
      RENOR1 =2.0*(RRZ-EZ)*RRZR1
      RENOR2 =2.0*(RRZ-EZ)*RRZR2
C
      DEL=(ALK2+BLK12)/DETAB
      DELA1=-(ALK2+BLK12)*DETA1/DETAB**2.0
      DELA2=( (1.0+DELTA(L,K))*DETAB-(ALK2+BLK12)*DETA2 )/DETAB**2.0
      DELB= ( (1.0+DELTA(L,K))*DETAB-(ALK2+BLK12)*DETB )/DETAB**2.0
C
      RE1=RENO/DEL
      RE1A1=( RENOA1*DEL-RENO*DELA1 )/DEL**2.0
      RE1A2=( RENOA2*DEL-RENO*DELA2 )/DEL**2.0
      RE1B =( RENOB*DEL-RENO*DELB )/DEL**2.0
      RE1R1=RENOR1/DEL
      RE1R2=RENOR2/DEL
C ***  FOR THE SECOND ELECTRON
      XX=ALK2*(ALK1+BLK12)*RZ2+ALK1*BLK12*RZ1
      XXA2=(1.0+DELTA(L,K))*(ALK1+BLK12)*RZ2 +ALK2*(ALK1+BLK12)*RZ2A2
      XXA1=ALK2*(1.0+DELTA(L,K))*RZ2+(1.0+DELTA(L,K))*BLK12*RZ1
     1    +ALK1*BLK12*RZ1A1
      XXB=ALK2*(1.0+DELTA(L,K))*RZ2+ALK1*(1.0+DELTA(L,K))*RZ1
      XXR2=ALK2*(ALK1+BLK12)*RZ2R2
      XXR1=ALK1*BLK12*RZ1R1
C
      RRZ=XX/DETAB 
      RRZA1=(XXA1*DETAB-XX*DETA1)/DETAB**2.0
      RRZA2=(XXA2*DETAB-XX*DETA2)/DETAB**2.0     
      RRZB=(XXB*DETAB-XX*DETB)/DETAB**2.0
      RRZR1=XXR1/DETAB 
      RRZR2=XXR2/DETAB 
C
      RENO=(RRZ-EZ)**2
      RENOA1=2.0*(RRZ-EZ)*RRZA1
      RENOA2=2.0*(RRZ-EZ)*RRZA2
      RENOB =2.0*(RRZ-EZ)*RRZB
      RENOR1 =2.0*(RRZ-EZ)*RRZR1
      RENOR2 =2.0*(RRZ-EZ)*RRZR2
C
      DEL=(ALK1+BLK12)/DETAB
      DELA2=-(ALK1+BLK12)*DETA2/DETAB**2.0
      DELA1=( (1.0+DELTA(L,K))*DETAB-(ALK1+BLK12)*DETA1 )/DETAB**2.0
      DELB= ( (1.0+DELTA(L,K))*DETAB-(ALK1+BLK12)*DETB )/DETAB**2.0
C
      RE2=RENO/DEL
      RE2A1=( RENOA1*DEL-RENO*DELA1 )/DEL**2.0
      RE2A2=( RENOA2*DEL-RENO*DELA2 )/DEL**2.0
      RE2B =( RENOB*DEL-RENO*DELB )/DEL**2.0
      RE2R1=RENOR1/DEL
      RE2R2=RENOR2/DEL
C
      RE=RE1+RE2
      REA1=RE1A1+RE2A1
      REA2=RE1A2+RE2A2
      REB=RE1B+RE2B
      RER1=RE1R1+RE2R1
      RER2=RE1R2+RE2R2
C
c      F0A1=-F1(RE) *REA1
c      F0A2=-F1(RE) *REA2
c      F0B =-F1(RE) *REB
c      F0R1=-F1(RE) *RER1
c      F0R2=-F1(RE) *RER2
C
 
C------THIS IS A TEST--------------------------------
      EN=G1*XK*RE
      DENA1=G1A1*XK*RE +G1*DKA1*RE +G1*XK*REA1
      DENA2=G1A2*XK*RE +G1*DKA2*RE +G1*XK*REA2
      DENB=G1B*XK*RE +G1*DKB*RE +G1*XK*REB
      DENR1=G1*DKR1*RE +G1*XK*RER1
      DENR2=G1*DKR2*RE +G1*XK*RER2
C--------------END OF TEST------------------------------


c      EN=G1*XK*F0(RE)
c      DENA1=G1A1*XK*F0(RE) +G1*DKA1*F0(RE) +G1*XK*F0A1
c      DENA2=G1A2*XK*F0(RE) +G1*DKA2*F0(RE) +G1*XK*F0A2
c      DENB=G1B*XK*F0(RE) +G1*DKB*F0(RE) +G1*XK*F0B
c      DENR1=G1*DKR1*F0(RE) +G1*XK*F0R1
c      DENR2=G1*DKR2*F0(RE) +G1*XK*F0R2
      

	RETURN
        END
