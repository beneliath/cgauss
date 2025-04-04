!      PROJECT TITLE:  correlated gaussian geminals for molecular hydrogen
!MAIN EXECUTION FILE:  cgauss.F90
!
!   Program FILENAME:  cgauss.exe
!        Lead Author:  D. Gilmore
!    Other Author(s):  L. Adamowicz, D. Kinghorn
!
!        LAST EDITED:  22may96:wed/17:48
!

!NOTES FOR FUTURE EDITS:
!===================================================================================================
!	1.  unnecessary sharing of modules should be eliminated
!===================================================================================================
!
!

!-------------------------------------------------------&   !=======================================
SUBROUTINE ProcessPipedInput								!PIPED INPUT -> PROGRAM SETTINGS

	USE program_settings_to_be_shared
	USE strings_to_be_shared

	CHARACTER(LEN=66) :: rSTR

	READ(*,*) rSTR

		DO WHILE (rSTR .NE. ":BEGIN:")						!wait for :BEGIN: tag...
			READ(*,*) rSTR
		END DO

		DO WHILE (rSTR .NE. ':END:')						!loop until end of section (:END:)...
			READ(*,*) rSTR

			IF (rSTR .EQ. 'NUMBER_OF_FUNCTIONS=')		&	!read number of basis fcns...
				READ(*,*) M

			IF (rSTR .EQ.								&	!name of file for opt wave fcn...
				'OPTIMIZED_WAVEFUNCTION_OUTPUT_FILE=') 	&
				READ(*,*) OPTIMIZEDWAVE

			IF (rSTR .EQ. 'OUTPUT_FILE=')				&	!general file output (file.out)...
				READ(*,*) FILENAME

			IF (rSTR .EQ. 'RESTART_FILE=')				&	!name for restart output...
				READ(*,*) RESTARTOUTPUT

			IF (rSTR .EQ. 'START_FROM_RESTART=') THEN		!start from restart?...
				READ(*,*) rSTR
				IF (rSTR .EQ.  'true') THEN
					RESTART=1
					ELSE
						RESTART=0
				END IF
			END IF

			IF (rSTR .EQ. 'RESTART_FROM=')				&	!name of file to restart from...
				READ(*,*) RESTARTNAME

			IF (rSTR .EQ.								&	!optimize before stoch growth?...
			'PRE-OPTIMIZE_BEFORE_STOCHASTIC_GROWTH=')	&
			THEN
				READ(*,*) rSTR
				IF (rSTR .EQ.  'true') THEN
					PREOPT=1
					ELSE
						PREOPT=0
				END IF
			END IF

			IF (rSTR .EQ. 'STOCHASTIC_END=')			&	!stochaistically grow to...
				READ(*,*) ISTEND

		END DO


END SUBROUTINE ProcessPipedInput
