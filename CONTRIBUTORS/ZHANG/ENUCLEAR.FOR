      SUBROUTINE ENUCLEAR(OVERLAP,AX,AY,AZ,A1,BX,BY,BZ,A2,                      
     -                 CX,CY,CZ,A3,DX,DY,DZ,A4,A5,A6,
     -                 EX,EY,EZ,EN)                       
C     -----------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C     -----------------------------------------------------
	PARAMETER(C1=0.564691197D-04)
	PARAMETER(C2=0.758433197D-03)
	PARAMETER(C3=0.769838037D-02)
	PARAMETER(C4=0.629344460D-01)
	PARAMETER(C5=0.213271302D+00)
	PARAMETER(C6=0.720266520D-04)
	PARAMETER(C7=0.955528842D-03)
	PARAMETER(C8=0.101431553D-01)
	PARAMETER(C9=0.738522953D-01)
	PARAMETER(C10=0.338450368D+00)
	PARAMETER(C11=0.879937801D+00)
	PARAMETER(ONE=1.00D+00)
	PARAMETER(TWO=2.00D+00)
	PARAMETER(FOUR=4.00D+00)
	PARAMETER(SIX=6.00D+00)
	PARAMETER(EIGHT=8.00D+00)
	PARAMETER(HALF=0.50D+00)
	PARAMETER(PI=3.1415926535898D+00)
C
C	...THE FOLLOWING PARAMETERS ARE ONLY FOR THE CONDITION OF
C		CYLINDRICAL SYMMETRY ABOUT THE Z-AXIS...	
C
	PARAMETER(RX1=0.)
	PARAMETER(RX2=0.)
	PARAMETER(RY1=0.)
	PARAMETER(RY2=0.)
	PARAMETER(RRX=0.)
	PARAMETER(RRY=0.)
C     -----------------------------------------------------
C     *** FOUR CENTER NUCLEAR ATTRACTION INTEGRAL ***
C     -----------------------------------------------------
      ALK1=A1+A2                                                       
      ALK2=A3+A4                                                       
      BLK12=A5+A6  
      ALKXX=ALK2+BLK12
      DETAB=ALK1*ALK2+(ALK1+ALK2)*BLK12
C     -----------------------------------------------------
C     --- ONE CENTER NUCLEAR INTEGRAL ---
C     -----------------------------------------------------
C      RX1=(A1*AX+A2*BX)/ALK1
C      RX2=(A3*CX+A4*DX)/ALK2
C      RY1=(A1*AY+A2*BY)/ALK1
C      RY2=(A3*CY+A4*DY)/ALK2
      RZ1=(A1*AZ+A2*BZ)/ALK1
      RZ2=(A3*CZ+A4*DZ)/ALK2
C     -----------------------------------------------------
C     --- FOUR-CENTER ELECTRON REPULSION INTEGRAL ---
C     -----------------------------------------------------
	EN=TWO*OVERLAP*DSQRT(DETAB/(PI*(ALK2+BLK12)))
C
C      RRX=(ALK1*(ALK2+BLK12)*RX1+ALK2*BLK12*RX2)/DETAB 
C      RRY=(ALK1*(ALK2+BLK12)*RY1+ALK2*BLK12*RY2)/DETAB 
      RRZ=(ALK1*(ALK2+BLK12)*RZ1+ALK2*BLK12*RZ2)/DETAB 
C
C      RE2=(RRX-EX)**TWO+(RRY-EY)**TWO+(RRZ-EZ)**TWO
      RE2=(RRZ-EZ)**TWO
C
      DELTA2=(ALK2+BLK12)/(FOUR*DETAB)
      RE2=RE2/(FOUR*DELTA2)
C
      EN=EN*F0(RE2) 
C
      RETURN                                                            
      END