!      PROJECT TITLE:  correlated gaussian geminals for molecular hydrogen
!MAIN EXECUTION FILE:  cgauss.F90
!
!   Program FILENAME:  cgauss.exe
!        Lead Author:  D. Gilmore
!    Other Author(s):  L. Adamowicz, D. Kinghorn
!
!        LAST EDITED:  22may96:wed/17:29
!

!NOTES FOR FUTURE EDITS:
!===================================================================================================
!	1.  unnecessary sharing of modules should be eliminated
!
!	2.	eliminate unnecessary vars in this module!
!===================================================================================================
!
!

!-------------------------------------------------------&   !=======================================
MODULE	program_settings_to_be_shared

	INTEGER ::  IMODE,METHODOPT,METHODGEN,RESTART,		&
				STGROW,PREOPT,ISTEND,STGROWPER,TEST

END MODULE	program_settings_to_be_shared
