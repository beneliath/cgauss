SUBROUTINE dOVERLAP_dA1(AX,AY,AZ,A1,BX,BY,BZ,A2,	&
					CX,CY,CZ,A3,DX,DY,DZ,A4,		&
                    A5,A6,DOV1,DOV2,DOV3,DOV4) 


	USE scalars_to_be_shared			  


	DOUBLE PRECISION :: AX,AY,AZ,A1,BX,BY,BZ,A2,	&
						CX,CY,CZ,A3,DX,DY,DZ,A4,	&
						A5,A6,DOV1,DOV2,DOV3,DOV4,	&
						DETSL,DETSK,CCL,CCK,CKL,	&
						ALK1,ALK2,BLK12,DETAB,		&
						AB2,CD2,DD,XK,				&
						RRX1,RRX2,RRY1,				&
						RRY2,RRZ1,RRZ2


		ALK1=TWO*A1
		ALK2=TWO*A3                                                       
		BLK12=TWO*A5  
		DETSL=(ALK1*ALK2)+((ALK1+ALK2)*BLK12)
	    CCL=EIGHT*A1*A3*A5*((AZ-CZ)**TWO)/DETSL

		ALK1=TWO*A2
		ALK2=TWO*A4                                                       
		BLK12=TWO*A6 
		DETSK=(ALK1*ALK2)+((ALK1+ALK2)*BLK12)
	    CCK=EIGHT*A2*A4*A6*((BZ-DZ)**TWO)/DETSK
!-----------------------------------------------------

      ALK1=A1+A2
      ALK2=A3+A4                                                       
      BLK12=A5+A6  
      DETAB=ALK1*ALK2+(ALK1+ALK2)*BLK12
!-----------------------------------------------------

	OV=(DABS(DSQRT(DABS(DSQRT(DETSL**THREE)))) *		&
		DABS(DSQRT(DABS(DSQRT(DETSK**THREE))))) /		&
		(DETAB*DSQRT(DETAB))
!-----------------------------------------------------

!	AB2=(AX-BX)**TWO + (AY-BY)**TWO + (AZ-BZ)**TWO
!	CD2=(CX-DX)**TWO + (CY-DY)**TWO + (CZ-DZ)**TWO

	AB2=(AZ-BZ)**TWO
	CD2=(CZ-DZ)**TWO

!      RRX1=(A1*AX+A2*BX)/ALK1
!      RRX2=(A3*CX+A4*DX)/ALK2
!      RRY1=(A1*AY+A2*BY)/ALK1
!      RRY2=(A3*CY+A4*DY)/ALK2

!!!      RRX1=0.00D+00
!!!      RRX2=0.00D+00
!!!      RRY1=0.00D+00
!!!      RRY2=0.00D+00

      RRZ1=(A1*AZ+A2*BZ)/ALK1
      RRZ2=(A3*CZ+A4*DZ)/ALK2

!      DD=((RRX1-RRX2)**TWO)+((RRY1-RRY2)**TWO)+((RRZ1-RRZ2)**TWO)

      DD=(RRZ1-RRZ2)**TWO

!	CKL=(((A1*A2)/ALK1)*((AX-BX)**TWO+(AY-BY)**TWO+(AZ-BZ)**TWO)+((A3*	&
!		A4)/ALK2)*((CX-DX)**TWO+(CY-DY)**TWO+(CZ-DZ)**TWO)+((ALK1*ALK2*	&
!		BLK12)*DD/DETAB))

	CKL=(((A1*A2)/ALK1)*((AZ-BZ)**TWO)+((A3*	&
		A4)/ALK2)*((CZ-DZ)**TWO)+((ALK1*ALK2*	&
		BLK12)*DD/DETAB))

	XK=EXP(-CKL+((CCK+CCL)/TWO))
!-----------------------------------------------------

      OV=OV*XK

!      DOVDA1=OV*(-(AB2*A2/ALK1)+AB2*A1*A2/ALK1**2.0D+00+ALK1*ALK2*BLK12*
!	-(ALK2+BLK12)*DD/DETAB**2.0D+00-ALK2*BLK12*DD/DETAB+(-8.0D+00*A1*A3
!     -*A5*(4.0D+00*A3+4.0D+00*A5)*(AZ-CZ)**2.0D+00/DETSL**2.0D+00+8.0D
!     -+00*A3*A5*(AZ-CZ)**2.0D+00/DETSL)/2.0D+00-ALK1*ALK2*BLK12*(2.0D+
!     -00*(AX/ALK1-(AX*A1+A2*BX)/ALK1**2.0D+00)*(RX1-RX2)+2.0D+00*(AY/ALK
!     -1-(AY*A1+A2*BY)/ALK1**2.0D+00)*(RY1-RY2)+2.0D+00*(AZ/ALK1-(AZ*A1+A
!     -2*BZ)/ALK1**2.0D+00)*(RZ1-RZ2))/DETAB)+3.0D+00*(4.0D+00*A3+4.0D+00
!     -*A5)*DETSK**(3.0D+00/4.0D+00)*XK/(4.0D+00*DETAB**(3.0D+00/2.0D+0
!     -0)*DETSL**(1.0D+00/4.0D+00))-3.0D+00*(ALK2+BLK12)*DETSK**(3.0D
!     -+00/4.0D+00)*DETSL**(3.0D+00/4.0D+00)*XK/(2.0D+00*DETAB**(5.0D+0
!     -0/2.0D+00))


END SUBROUTINE dOVERLAP_dA1