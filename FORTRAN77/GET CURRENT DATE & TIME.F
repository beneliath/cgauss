!      PROJECT TITLE:  correlated gaussian geminals for molecular hydrogen
!MAIN EXECUTION FILE:  cgauss.F90
!
!   Program FILENAME:  cgauss.exe
!        Lead Author:  D. Gilmore
!    Other Author(s):  L. Adamowicz, D. Kinghorn
!
!        LAST EDITED:  22may96:wed/18:00
!

!NOTES FOR FUTURE EDITS:
!===================================================================================================
!	1.  unnecessary sharing of modules should be eliminated
!===================================================================================================
!
!

!-------------------------------------------------------&   !=======================================
SUBROUTINE Get_Current_Time(TD,TT,TZ)

	INTRINSIC DATE_AND_TIME

	CHARACTER(LEN=10) :: TT
	CHARACTER(LEN=5)  :: TZ
	CHARACTER(LEN=8)  :: TD

	CALL DATE_AND_TIME (DATE=TD,TIME=TT,ZONE=TZ )

END SUBROUTINE Get_Current_Time
