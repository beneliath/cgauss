!      PROJECT TITLE:  correlated gaussian geminals for molecular hydrogen
!MAIN EXECUTION FILE:  cgauss.F90
!
!   Program FILENAME:  cgauss.exe
!        Lead Author:  D. Gilmore
!    Other Author(s):  L. Adamowicz, D. Kinghorn
!
!        LAST EDITED:  22may96:wed/17:40
!

!NOTES FOR FUTURE EDITS:
!===================================================================================================
!	1.  unnecessary sharing of modules should be eliminated
!===================================================================================================
!
!

!-------------------------------------------------------&   !=======================================
SUBROUTINE Default_Program_Settings
!-----------------------------------------------------------
	USE scalars_to_be_shared
	USE	strings_to_be_shared			  
	USE	program_settings_to_be_shared
!-----------------------------------------------------------
		TEST   = 0											!DEFAULT PROGRAM SETTINGS...
		IMODE  = 0									
	    M      = 1											!number of basis functions to use
		FILENAME      = 'H2_PROJECT.OUT'				
		RESTARTNAME   = 'restart.dat'					
		RESTARTOUTPUT = 'RESTART_OUT.DAT'			
		OPTIMIZEDWAVE = 'OPTIMAL_WAVEFUNCTION.DAT'
		RESTART= 1								
		PREOPT = 1								
		ISTEND = 40								
		BD     = 1.4011D+00									!experimental equillibrium bond distance
															!of H2

END	SUBROUTINE Default_Program_Settings