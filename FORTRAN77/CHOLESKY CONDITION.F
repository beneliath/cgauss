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
!===================================================================================================
!
!

!-------------------------------------------------------&   !=======================================
SUBROUTINE Cholesky_Condition_Parameters					!cholesky factors -> alphas & betas
!-----------------------------------------------------------
   USE scalars_to_be_shared
   USE vectors_to_be_shared
!-----------------------------------------------------------
	DO I=1,M
		AA1(I)= T1(I)*T1(I)               + T1(I)*T3(I)
        AA2(I)= T2(I)*T2(I) + T3(I)*T3(I) + T1(I)*T3(I)
        B12(I)=                           - T1(I)*T3(I)
	END DO

END SUBROUTINE Cholesky_Condition_Parameters
