!      PROJECT TITLE:  correlated gaussian geminals for molecular hydrogen
!MAIN EXECUTION FILE:  cgauss.F90
!
!   Program FILENAME:  cgauss.exe
!        Lead Author:  D. Gilmore
!    Other Author(s):  L. Adamowicz, D. Kinghorn
!
!        LAST EDITED:  22may96:wed/17:59
!

!NOTES FOR FUTURE EDITS:
!===================================================================================================
!
!===================================================================================================
!
!

!-------------------------------------------------------&   !=======================================
SUBROUTINE BIOUT(STRING)

      CHARACTER*(*) STRING
      WRITE(6,*) STRING
      WRITE(7,*) STRING

END SUBROUTINE BIOUT