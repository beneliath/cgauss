C     L. ADAMOWICZ RESEARCH / UNIVERSITY OF ARIZONA / CHEM.
C     -----------------------------------------------------
      SUBROUTINE OVERLP(M,DFLAG,AX,AY,AZ,A1,BX,BY,BZ,A2,                      
     1                  CX,CY,CZ,A3,DX,DY,DZ,A4,
     2                  A5,A6,OV,KLM) 
C     ----------------------------------------------------- 
      IMPLICIT REAL*8 (A-H,O-Z)
	REAL*8 DETSL,DETSK,CCL,CCK
	INTEGER M,DFLAG
	COMMON /DOV/ DOVDA1,DOVDA2,DOVDA5,DOVDAZ,DOVDBZ

C     -----------------------------------------------------
      DATA ONE,TWO,THREE,FOUR,FIVE,SIX,EIGHT/1.,2.,3.,4.,5.,6.,8./
	DATA HALF/0.5/
C     -----------------------------------------------------
		ALK1=2.0D+00*A1
		ALK2=2.0D+00*A3                                                       
		BLK12=2.0D+00*A5  
		DETSL=(ALK1*ALK2)+((ALK1+ALK2)*BLK12)
	    CCL=8*A1*A3*A5*((AZ-CZ)**2)/DETSL
C
		ALK1=2.0D+00*A2
		ALK2=2.0D+00*A4                                                       
		BLK12=2.0D+00*A6 
		DETSK=(ALK1*ALK2)+((ALK1+ALK2)*BLK12)
	    CCK=8*A2*A4*A6*((BZ-DZ)**2)/DETSK
C     -----------------------------------------------------
      ALK1=A1+A2
      ALK2=A3+A4                                                       
      BLK12=A5+A6  
      DETAB=ALK1*ALK2+(ALK1+ALK2)*BLK12
C
      IF(DETAB.LE.ZERO) GOTO 300
C     -----------------------------------------------------
	OV=(DABS(DSQRT(DABS(DSQRT(DETSL**3.0D+00)))) *
	-	DABS(DSQRT(DABS(DSQRT(DETSK**3.0D+00))))) /
	-	(DETAB*DSQRT(DETAB))
C     -----------------------------------------------------
	AB2=(AX-BX)**2 + (AY-BY)**2 + (AZ-BZ)**2
	CD2=(CX-DX)**2 + (CY-DY)**2 + (CZ-DZ)**2
      RX1=(A1*AX+A2*BX)/ALK1
      RX2=(A3*CX+A4*DX)/ALK2
      RY1=(A1*AY+A2*BY)/ALK1
      RY2=(A3*CY+A4*DY)/ALK2
      RZ1=(A1*AZ+A2*BZ)/ALK1
      RZ2=(A3*CZ+A4*DZ)/ALK2
C
      DD=((RX1-RX2)**2)+((RY1-RY2)**2)+((RZ1-RZ2)**2)
C
	CKL=(((A1*A2)/ALK1)*((AX-BX)**2+(AY-BY)**2+(AZ-BZ)**2)+((A3*A4
	-)/ALK2)*((CX-DX)**2+(CY-DY)**2+(CZ-DZ)**2)+((ALK1*ALK2*BLK12)*DD/D
     -ETAB))
C
	XK=EXP(-CKL+((CCK+CCL)/2.0D+00))
C     -----------------------------------------------------
      OV=OV*XK
C     -----------------------------------------------------
C	 DERIVATIVES OF THE OVERLAP:
C     -----------------------------------------------------
	IF (DFLAG .EQ. 1) THEN
C
      DOVDA1=OV*(-(AB2*A2/ALK1)+AB2*A1*A2/ALK1**2.0D+00+ALK1*ALK2*BLK12*
	-(ALK2+BLK12)*DD/DETAB**2.0D+00-ALK2*BLK12*DD/DETAB+(-8.0D+00*A1*A3
     -*A5*(4.0D+00*A3+4.0D+00*A5)*(AZ-CZ)**2.0D+00/DETSL**2.0D+00+8.0D
     -+00*A3*A5*(AZ-CZ)**2.0D+00/DETSL)/2.0D+00-ALK1*ALK2*BLK12*(2.0D+
     -00*(AX/ALK1-(AX*A1+A2*BX)/ALK1**2.0D+00)*(RX1-RX2)+2.0D+00*(AY/ALK
     -1-(AY*A1+A2*BY)/ALK1**2.0D+00)*(RY1-RY2)+2.0D+00*(AZ/ALK1-(AZ*A1+A
     -2*BZ)/ALK1**2.0D+00)*(RZ1-RZ2))/DETAB)+3.0D+00*(4.0D+00*A3+4.0D+00
     -*A5)*DETSK**(3.0D+00/4.0D+00)*XK/(4.0D+00*DETAB**(3.0D+00/2.0D+0
     -0)*DETSL**(1.0D+00/4.0D+00))-3.0D+00*(ALK2+BLK12)*DETSK**(3.0D
     -+00/4.0D+00)*DETSL**(3.0D+00/4.0D+00)*XK/(2.0D+00*DETAB**(5.0D+0
     -0/2.0D+00))
C
	DOVDA2=OV*(-(A4*CD2/ALK2)+A3*A4*CD2/ALK2**2.0D+00+ALK1*ALK2*BLK12*
	-(ALK1+BLK12)*DD/DETAB**2.0D+00-ALK1*BLK12*DD/DETAB+(-8.0D+00*A1*A3
     -*A5*(4.0D+00*A1+4.0D+00*A5)*(AZ-CZ)**2.0D+00/DETSL**2.0D+00+8.0D
     -+00*A1*A5*(AZ-CZ)**2.0D+00/DETSL)/2.0D+00-ALK1*ALK2*BLK12*(2.0D+
     -00*(-(CX/ALK2)+(A3*CX+A4*DX)/ALK2**2.0D+00)*(RX1-RX2)+2.0D+00*(-(C
     -Y/ALK2)+(A3*CY+A4*DY)/ALK2**2.0D+00)*(RY1-RY2)+2.0D+00*(-(CZ/ALK2)
     -+(A3*CZ+A4*DZ)/ALK2**2.0D+00)*(RZ1-RZ2))/DETAB)+3.0D+00*(4.0D+00*A
     -1+4.0D+00*A5)*DETSK**(3.0D+00/4.0D+00)*XK/(4.0D+00*DETAB**(3.0D+
     -00/2.0D+00)*DETSL**(1.0D+00/4.0D+00))-3.0D+00*(ALK1+BLK12)*DETSK
     -**(3.0D+00/4.0D+00)*DETSL**(3.0D+00/4.0D+00)*XK/(2.0D+00*DETAB
     -**(5.0D+00/2.0D+00))
C
	DOVDA5=(ALK1*ALK2*(ALK1+ALK2)*BLK12*DD/DETAB**2.0D+00-ALK1*ALK2*DD
	-/DETAB+(-16.0D+00*A1*A3*(2.0D+00*A1+2.0D+00*A3)*A5*(AZ-CZ)**2.0D+0
     -0/DETSL**2.0D+00+8.0D+00*A1*A3*(AZ-CZ)**2.0D+00/DETSL)/2.0D+00
     -)*OV+3.0D+00*(2.0D+00*A1+2.0D+00*A3)*DETSK**(3.0D+00/4.0D+00)*XK
     -/(2.0D+00*DETAB**(3.0D+00/2.0D+00)*DETSL**(1.0D+00/4.0D+00))-3.0
     -D+00*(ALK1+ALK2)*DETSK**(3.0D+00/4.0D+00)*DETSL**(3.0D+00/4.0D
     -+00)*XK/(2.0D+00*DETAB**(5.0D+00/2.0D+00))
C
	DOVDAZ=OV*(-2.0D+00*A1*A2*(AZ-BZ)/ALK1+8.0D+00*A1*A3*A5*(AZ-CZ)/DE
	-TSL-2.0D+00*ALK2*A1*BLK12*(RZ1-RZ2)/DETAB)
C
	DOVDBZ=OV*(-8.0D+00*A1*A3*A5*(AZ-CZ)/DETSL-2.0D+00*A3*A4*(CZ-DZ)
	-/ALK2+2.0D+00*ALK1*A3*BLK12*(RZ1-RZ2)/DETAB)
C     -----------------------------------------------------
      GOTO 1500
C     -----------------------------------------------------
300   KLM=5
C     -----------------------------------------------------
1500  END IF


c	WRITE(6,*) OV

	RETURN                                                            
      END