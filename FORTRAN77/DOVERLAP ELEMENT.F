!      PROJECT TITLE:  correlated gaussian geminals for molecular hydrogen
!MAIN EXECUTION FILE:  cgauss.F90
!
!   Program FILENAME:  cgauss.exe
!        Lead Author:  D. Gilmore
!    Other Author(s):  L. Adamowicz, D. Kinghorn
!
!        LAST EDITED:  22may96:wed/20:05
!

!NOTES FOR FUTURE EDITS:
!===================================================================================================
!	1.  unnecessary sharing of modules should be eliminated
!
!	2.	this subroutine --PASSES-- ALL 5 finite difference tests!
!===================================================================================================
!
!

!-------------------------------------------------------&   !=======================================
SUBROUTINE dOVERLAP(AZ,A1,BZ,A2,CZ,A3,DZ,A4,		&
                    A5,A6,							&
					DOV1,DOV2,DOV3,DOV4,DOV5,OV) 


	USE scalars_to_be_shared


	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	IMPLICIT INTEGER (I-N)


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

		ALK1=A1+A2
		ALK2=A3+A4                                                       
		BLK12=A5+A6  
		DETAB=ALK1*ALK2+(ALK1+ALK2)*BLK12

!-----------------------------------------------------

	OV=(DSQRT(DSQRT(DABS(DETSL*DETSL*DETSL))) *		&
		DSQRT(DSQRT(DABS(DETSK*DETSK*DETSK)))) /		&
		(DETAB*DSQRT(DETAB))

!-----------------------------------------------------

	AB2=ABZ*ABZ
	CD2=CDZ*CDZ

      RZ1=(A1*AZ+A2*BZ)/ALK1
      RZ2=(A3*CZ+A4*DZ)/ALK2

	RZ12=RZ1-RZ2

      DD=RZ12*RZ12

	CKL=(((A1*A2)/ALK1)*ABZ*ABZ+((A3*A4)/ALK2)*CDZ*CDZ	&
	+((ALK1*ALK2*BLK12)*DD/DETAB))

	XK=EXP(-CKL+(HALF*(CCK+CCL)))

!-----------------------------------------------------

      OV=OV*XK


      t2 = ONE/ALK1
      t6 = T2/ALK1
      t10 = ALK1*ALK2*BLK12
      t11 = ALK2+BLK12
      t13 = DETAB*DETAB
      t18 = ONE/DETAB
      t24 = A3+A5
      t27 = ACZ*ACZ
      t28 = DETSL*DETSL
      t60 = DSQRT(DABS(DETAB*DETAB*DETAB))
      t57 = DSQRT(DSQRT(DABS(DETSK*DETSK*DETSK)))
      t64 = DSQRT(DSQRT(DABS(DETSL)))
      t70 = DSQRT(DSQRT(DABS(DETSL*DETSL*DETSL)))
      t73 = DSQRT(DABS(DETAB*DETAB*DETAB*DETAB*DETAB))

      DOV1=OV*(-AB2*A2*t2+AB2*A1*A2*t6+t10*t11*DD/t13-ALK2*BLK12*DD*t18+		&
	  (-EIGHT*A1*A3*A5*FOUR*t24*t27/t28+EIGHT*A3*A5*t27/DETSL)*HALF-t10*TWO*	&
	  (AZ*t2-(AZ*A1+A2*BZ)*t6)*RZ12*t18)+THREE*t24*t57*XK/t60/t64-				&
	  THREE*t11*t57*t70*XK*HALF/t73

      tt2 = ONE/ALK2
      tt6 = tt2/ALK2
      tt11 = ALK1+BLK12
      t21 = EIGHT*A1
      t25 = FOUR*(A1+A5)

      DOV2 = OV*(-A4*CD2*tt2+A3*A4*CD2*tt6+t10*tt11*DD/t13-ALK1*BLK12*DD*t18+	&
	  (-t21*A3*A5*t25*t27/t28+t21*A5*t27/DETSL)*HALF-t10*TWO*(-CZ*tt2+(A3*CZ+	&
	  A4*DZ)*tt6)*RZ12*t18)+THREE*t25*t57*XK*QUARTER/t60/t64-THREE*tt11*		&
	  t57*t70*XK*HALF/t73


      t1 = ALK1*ALK2
      ttt2 = ALK1+ALK2
      t16 = TWO*(A1+A3)

      DOV3 = (t1*ttt2*BLK12*DD/t13-t1*DD/DETAB+(-SIXTEEN*A1*A3*t16*A5*t27/t28+		&
	  EIGHT*A1*A3*t27/DETSL)*HALF)*OV+THREE*t16*t57*XK*HALF/t60/t64-THREE*ttt2*t57*	&
	  t70*XK*HALF/t73


	  DOV4=OV*(-t2*TWO*A1*A2*ABZ+EIGHT*A1*A3*A5*ACZ/DETSL-t18*TWO*ALK2*A1*BLK12*RZ12)


	  DOV5=OV*(-EIGHT*A1*A3*A5*ACZ/DETSL-tt2*TWO*A3*A4*CDZ+t18*TWO*ALK1*A3*BLK12*RZ12)


END SUBROUTINE dOVERLAP
