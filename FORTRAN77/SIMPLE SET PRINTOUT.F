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
SUBROUTINE Simple_Wave_Function_Printout					!print wave fcn in simple format
!-----------------------------------------------------------
	USE scalars_to_be_shared
	USE vectors_to_be_shared
!-----------------------------------------------------------! ---:END OF FORMATTED COMMENTS:---

	CALL BIOUT('  AA1         AA2         B12         C           RZ1         RZ2')
    CALL BIOUT('  ---------------------------------------------------------------------')

	DO I=1,M
				WRITE(6,709)AA1(I),AA2(I),B12(I),C(I),RZ1(I),RZ2(I)
				WRITE(7,709)AA1(I),AA2(I),B12(I),C(I),RZ1(I),RZ2(I)
	END DO

 	CALL BIOUT('  =====================================================================')
	CALL BIOUT(' ')


!FORMAT(s):
!-----------------------------------------------------------
	709  FORMAT(6E12.4)

END SUBROUTINE Simple_Wave_Function_Printout
