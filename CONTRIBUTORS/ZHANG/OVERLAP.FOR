      SUBROUTINE OVERLAP(AX,AY,AZ,A1,BX,BY,BZ,A2,                      
     -                  CX,CY,CZ,A3,DX,DY,DZ,A4,
     -                  A5,A6,OV) 
C     ----------------------------------------------------- 
      IMPLICIT REAL*8 (A-H,O-Z)
	REAL*8 DETSL,DETSK,CCL,CCK
C     -----------------------------------------------------
	PARAMETER(ONE=1.00D+00)
	PARAMETER(TWO=2.00D+00)
	PARAMETER(THREE=3.00D+00)
	PARAMETER(FOUR=4.00D+00)
	PARAMETER(FIVE=5.00D+00)
	PARAMETER(SIX=6.00D+00)
	PARAMETER(EIGHT=8.00D+00)
	PARAMETER(HALF=0.50D+00)
C     -----------------------------------------------------
		ALK1=TWO*A1
		ALK2=TWO*A3                                                       
		BLK12=TWO*A5  
		DETSL=(ALK1*ALK2)+((ALK1+ALK2)*BLK12)
	    CCL=EIGHT*A1*A3*A5*((AZ-CZ)**TWO)/DETSL
C
		ALK1=TWO*A2
		ALK2=TWO*A4                                                       
		BLK12=TWO*A6 
		DETSK=(ALK1*ALK2)+((ALK1+ALK2)*BLK12)
	    CCK=EIGHT*A2*A4*A6*((BZ-DZ)**TWO)/DETSK
C     -----------------------------------------------------
      ALK1=A1+A2
      ALK2=A3+A4                                                       
      BLK12=A5+A6  
      DETAB=ALK1*ALK2+(ALK1+ALK2)*BLK12
C     -----------------------------------------------------
	OV=(DABS(DSQRT(DABS(DSQRT(DETSL**THREE)))) *
	-	DABS(DSQRT(DABS(DSQRT(DETSK**THREE))))) /
	-	(DETAB*DSQRT(DETAB))
C     -----------------------------------------------------
	AB2=(AX-BX)**TWO + (AY-BY)**TWO + (AZ-BZ)**TWO
	CD2=(CX-DX)**TWO + (CY-DY)**TWO + (CZ-DZ)**TWO
      RX1=(A1*AX+A2*BX)/ALK1
      RX2=(A3*CX+A4*DX)/ALK2
      RY1=(A1*AY+A2*BY)/ALK1
      RY2=(A3*CY+A4*DY)/ALK2
      RZ1=(A1*AZ+A2*BZ)/ALK1
      RZ2=(A3*CZ+A4*DZ)/ALK2
C
      DD=((RX1-RX2)**TWO)+((RY1-RY2)**TWO)+((RZ1-RZ2)**TWO)
C
	CKL=(((A1*A2)/ALK1)*((AX-BX)**TWO+(AY-BY)**TWO+(AZ-BZ)**TWO)+((A3*
	-A4)/ALK2)*((CX-DX)**TWO+(CY-DY)**TWO+(CZ-DZ)**TWO)+((ALK1*ALK2*
     -BLK12)*DD/DETAB))
C
	XK=EXP(-CKL+((CCK+CCL)/TWO))
C     -----------------------------------------------------
      OV=OV*XK
C     -----------------------------------------------------
	RETURN                                                            
      END