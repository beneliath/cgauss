!      PROJECT TITLE:  correlated gaussian geminals for molecular hydrogen
!MAIN EXECUTION FILE:  cgauss.F90
!
!   Program FILENAME:  cgauss.exe
!        Lead Author:  D. Gilmore
!    Other Author(s):  L. Adamowicz, D. Kinghorn
!
!        LAST EDITED:  22may96:wed/17:28
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
MODULE	strings_to_be_shared

	CHARACTER(LEN=66) ::	FILENAME,RESTARTNAME,		&
							RESTARTOUTPUT,OPTIMIZEDWAVE,&
							METHODNAME,GENNAME

END MODULE strings_to_be_shared