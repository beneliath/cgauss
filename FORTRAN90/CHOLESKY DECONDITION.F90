!      PROJECT TITLE:  correlated gaussian geminals for molecular hydrogen
!MAIN EXECUTION FILE:  cgauss.F90
!
!   Program FILENAME:  cgauss.exe
!        Lead Author:  D. Gilmore
!    Other Author(s):  L. Adamowicz, D. Kinghorn
!
!        LAST EDITED:  22may96:wed/18:21
!

!NOTES FOR FUTURE EDITS:
!===================================================================================================
!	1.  unnecessary sharing of modules should be eliminated
!===================================================================================================
!
!

!-------------------------------------------------------&   !=======================================
SUBROUTINE Cholesky_Deconditioning							!alphas & betas -> cholesky factors
!-----------------------------------------------------------
	USE scalars_to_be_shared
	USE vectors_to_be_shared
!-----------------------------------------------------------

	DO I=1,M
			T1(I) = DSQRT(AA1(I)+B12(I))	
			T2(I) = DSQRT(AA2(I)+B12(I)- ((B12(I)*B12(I))/(AA1(I)+B12(I))))
			T3(I) =-B12(I)/DSQRT(AA1(I)+B12(I))
	END DO

END SUBROUTINE Cholesky_Deconditioning