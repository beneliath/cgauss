!      PROJECT TITLE:  correlated gaussian geminals for molecular hydrogen
!MAIN EXECUTION FILE:  cgauss.F90
!
!   Program FILENAME:  cgauss.exe
!        Lead Author:  D. Gilmore
!    Other Author(s):  L. Adamowicz, D. Kinghorn
!
!        LAST EDITED:  23may96:wed/08:07
!

!NOTES FOR FUTURE EDITS:
!===================================================================================================
!	1.  unnecessary sharing of modules should be eliminated
!===================================================================================================
!
!

!-------------------------------------------------------&   !=======================================
SUBROUTINE Calculate_Rayleigh_Quotient
	USE scalars_to_be_shared
	USE vectors_to_be_shared
	USE matrices_to_be_shared

	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	IMPLICIT          INTEGER    (I-N)
!-----------------------------------------------------------
	HE = ZERO
	SE = ZERO

	DO I=1,M
		DO J=1,M
			HE = HE + (C(I)*H(I,J)*C(J))
			SE = SE + (C(I)*S(I,J)*C(J))
		END DO
	END DO

	R_ENER=HE/SE

END SUBROUTINE Calculate_Rayleigh_Quotient
