!      PROJECT TITLE:  correlated gaussian geminals for molecular hydrogen
!MAIN EXECUTION FILE:  cgauss.F90
!
!   Program FILENAME:  cgauss.exe
!        Lead Author:  D. Gilmore
!    Other Author(s):  L. Adamowicz, D. Kinghorn
!
!        LAST EDITED:  22may96:wed/17:58
!

!NOTES FOR FUTURE EDITS:
!===================================================================================================
!	1.  unnecessary sharing of modules should be eliminated
!===================================================================================================
!
!

!-------------------------------------------------------&   !=======================================
SUBROUTINE Print_Preliminary_Matter							!PRINTS HEADERs TO SCREEN & file.out...
!-----------------------------------------------------------
	USE          scalars_to_be_shared						!shared modules...
	USE	         vectors_to_be_shared			  
	USE	        matrices_to_be_shared			  
	USE	         strings_to_be_shared			  
	USE	program_settings_to_be_shared

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)						!implicit real typing
	IMPLICIT          INTEGER     (I-N)						!implicit int  typing
!-----------------------------------------------------------
    CALL BIOUT('CGAUSS:  H-H')
    CALL BIOUT('===================================='	&
				//'===============')
    CALL BIOUT('ADAMOWICZ RESEARCH                 '	&
				//'[UofA chemistry]')
	CALL BIOUT('Programming:[D. Gilmore, L. Adamowicz,'	&
				//' D. Kinghorn]') 
    CALL BIOUT('===================================='	&
				//'===============')

	METHODNAME='IMSL:DUMING'

	IF(RESTART.EQ.1) THEN
		GENNAME='RESTART'
		M=0
		ELSE
			GENNAME='RANDOM'
	END IF

	IF(TEST .EQ. 1) GENNAME='INTERNAL TEST PARAMETERS'

	DO I=6,7
		WRITE(I,*)'NUMBER OF FUNCTIONS =',M
		WRITE(I,*)'FUNCTION GENERATION =',GENNAME
		WRITE(I,*)'OPTIMIZATION METHOD =',METHODNAME
		WRITE(I,*)'OUTPUT FILENAME     =',FILENAME
	    WRITE(I,*)'NEW RESTART OUTPUT  =',RESTARTOUTPUT
	
		IF(RESTART .EQ. 1)								&
			WRITE(I,*)'       RESTART FROM=',RESTARTNAME
	END DO

END SUBROUTINE Print_Preliminary_Matter