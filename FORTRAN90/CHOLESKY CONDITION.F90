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
	OPEN (5,FILE="POINT1.DAT")
	OPEN (8,FILE="POINT2.DAT")
	DO I=1,M
		AA1(I)= T1(I)*T1(I)               + T1(I)*T3(I)
        AA2(I)= T2(I)*T2(I) + T3(I)*T3(I) + T1(I)*T3(I)
        B12(I)=                           - T1(I)*T3(I)

		WRITE(5,*)AA1(I),RZ1(I)
		WRITE(8,*)AA2(I),RZ2(I)

	END DO
	CLOSE(5)
	CLOSE(8)

END SUBROUTINE Cholesky_Condition_Parameters