SUBROUTINE OVERLAP(	AZ,A1,BZ,A2,CZ,A3,DZ,A4,A5,A6,OV) 


	USE scalars_to_be_shared			  


	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	IMPLICIT          INTEGER     (I-N)

		ACZ=AZ-CZ
		BDZ=BZ-DZ

		ABZ=AZ-BZ
		CDZ=CZ-DZ

		ALK1=TWO*A1
		ALK2=TWO*A3                                                       
		BLK12=TWO*A5  
		DETSL=(ALK1*ALK2)+((ALK1+ALK2)*BLK12)
	    CCL=EIGHT*A1*A3*A5*(ACZ*ACZ)/DETSL

		ALK1=TWO*A2
		ALK2=TWO*A4                                                       
		BLK12=TWO*A6 
		DETSK=(ALK1*ALK2)+((ALK1+ALK2)*BLK12)
	    CCK=EIGHT*A2*A4*A6*(BDZ*BDZ)/DETSK
!-----------------------------------------------------

      ALK1=A1+A2
      ALK2=A3+A4                                                       
      BLK12=A5+A6  
      DETAB=ALK1*ALK2+(ALK1+ALK2)*BLK12
!-----------------------------------------------------

	OV=(DSQRT(DSQRT(DABS(DETSL*DETSL*DETSL))) *		&
		DSQRT(DSQRT(DABS(DETSK*DETSK*DETSK)))) /		&
		(DETAB*DSQRT(DETAB))
!-----------------------------------------------------

!	AB2=(AX-BX)**TWO + (AY-BY)**TWO + (AZ-BZ)**TWO
!	CD2=(CX-DX)**TWO + (CY-DY)**TWO + (CZ-DZ)**TWO

	AB2=ABZ*ABZ
	CD2=CDZ*CDZ

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

	RRZ12=RRZ1-RRZ2

!      DD=((RRX1-RRX2)**TWO)+((RRY1-RRY2)**TWO)+((RRZ1-RRZ2)**TWO)

      DD=RRZ12*RRZ12

!	CKL=(((A1*A2)/ALK1)*((AX-BX)**TWO+(AY-BY)**TWO+(AZ-BZ)**TWO)+((A3*	&
!		A4)/ALK2)*((CX-DX)**TWO+(CY-DY)**TWO+(CZ-DZ)**TWO)+((ALK1*ALK2*	&
!		BLK12)*DD/DETAB))

	CKL=(((A1*A2)/ALK1)*(ABZ*ABZ)+((A3*		&
		A4)/ALK2)*(CDZ*CDZ)+((ALK1*ALK2*	&
		BLK12)*DD/DETAB))

	XK=EXP(-CKL+(HALF*(CCK+CCL)))
!-----------------------------------------------------

      OV=OV*XK


END SUBROUTINE OVERLAP