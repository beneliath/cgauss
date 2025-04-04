C****************************************************
      SUBROUTINE ENUCL(AX,AY,AZ,A1,BX,BY,BZ,A2,                      
     1                 CX,CY,CZ,A3,DX,DY,DZ,A4,A5,A6,
     2                 EX,EY,EZ,EN,klm)                       
C****************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C ***                                                 ***
C ***     FOUR CENTER NUCLEAR ATTRACTION INTEGRAL     ***
C ***                                                 ***
      DATA ONE,TWO,FOUR,SIX,EIGHT,HALF/1.,2.,4.,6.,8.,.5/               
      PI=3.1415926535898D+00
      PI52=PI**2.5D+00
C *** --------------------------------------------------   
C *** P.M.Kozlowski & L. Adamowicz J.C.P 95, 6681 (1991)    
C *** --------------------------------------------------   
      ALK1=A1+A2                                                       
      ALK2=A3+A4                                                       
      BLK12=A5+A6  
      ALKXX=ALK2+BLK12	
      DETAB=ALK1*ALK2+(ALK1+ALK2)*BLK12
      IF(ALKXX.LE.0.0D+00) GOTO 300
      IF(DETAB.LE.0.0D+00) GOTO 300
C
C --- ONE CENTER NUCLEAR INTEGRAL ---
C
      EN=TWO*PI52/(DETAB*DSQRT(ALK2+BLK12))
C
C --- CENTERS ---
C
      XX1=(AX-BX)**2+(AY-BY)**2+(AZ-BZ)**2
      XK1=EXP(-((A1*A2)/ALK1)*XX1) 
      XX2=(CX-DX)**2+(CY-DY)**2+(CZ-DZ)**2
      XK2=EXP(-((A3*A4)/ALK2)*XX2) 
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
C     WRITE(6,*) ' ONE CENTER NUC. ATTRACTION INTEGRAL = ',EN
      EN=EN*XK1*XK2*ZZ
      RRX=(ALK1*(ALK2+BLK12)*RX1+ALK2*BLK12*RX2)/DETAB 
      RRY=(ALK1*(ALK2+BLK12)*RY1+ALK2*BLK12*RY2)/DETAB 
      RRZ=(ALK1*(ALK2+BLK12)*RZ1+ALK2*BLK12*RZ2)/DETAB 
      RE2=(RRX-EX)**2+(RRY-EY)**2+(RRZ-EZ)**2
      DELTA2=(ALK2+BLK12)/(FOUR*DETAB)
      RE2=RE2/(FOUR*DELTA2)
      EN=EN*F0(RE2) 
C
	goto 1500
300   KLM=5

1500   RETURN                                                            
      END 
