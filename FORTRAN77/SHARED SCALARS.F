!      PROJECT TITLE:  correlated gaussian geminals for molecular hydrogen
!MAIN EXECUTION FILE:  cgauss.F90
!
!   Program FILENAME:  cgauss.exe
!        Lead Author:  D. Gilmore
!    Other Author(s):  L. Adamowicz, D. Kinghorn
!
!        LAST EDITED:  22may96:wed/17:04
!

!NOTES FOR FUTURE EDITS:
!===================================================================================================
!	1.  unnecessary sharing of modules should be eliminated
!
!	2.	nuclear geometry variables(EX1-EZ1,EX2-EZ2) can be reduced under cylindrical symmetry
!===================================================================================================
!
!

!-------------------------------------------------------&   !=======================================
MODULE	scalars_to_be_shared

	DOUBLE PRECISION, PARAMETER :: ZERO    = 0.00D+00		!REAL NUMERICAL CONSTANTS...
	DOUBLE PRECISION, PARAMETER :: HALF    = 0.50D+00
	DOUBLE PRECISION, PARAMETER :: QUARTER = 0.25D+00
	DOUBLE PRECISION, PARAMETER :: EIGHTH  = QUARTER*HALF
	DOUBLE PRECISION, PARAMETER :: ONE     = 1.00D+00
	DOUBLE PRECISION, PARAMETER :: TWO     = 2.00D+00
	DOUBLE PRECISION, PARAMETER :: THREE   = 3.00D+00
	DOUBLE PRECISION, PARAMETER :: FOUR    = 4.00D+00
	DOUBLE PRECISION, PARAMETER :: FIVE    = 5.00D+00
	DOUBLE PRECISION, PARAMETER :: SIX     = 6.00D+00
	DOUBLE PRECISION, PARAMETER :: SEVEN   = 7.00D+00
	DOUBLE PRECISION, PARAMETER :: EIGHT   = 8.00D+00
	DOUBLE PRECISION, PARAMETER :: SIXTEEN = 16.00D+00

	DOUBLE PRECISION, PARAMETER :: C1  = 0.564691197D-04	!REAL PADE' APPROXIMATION CONSTANTS...
	DOUBLE PRECISION, PARAMETER :: C2  = 0.758433197D-03	!	...for the error fnc approximation...
	DOUBLE PRECISION, PARAMETER :: C3  = 0.769838037D-02
	DOUBLE PRECISION, PARAMETER :: C4  = 0.629344460D-01
	DOUBLE PRECISION, PARAMETER :: C5  = 0.213271302D+00
	DOUBLE PRECISION, PARAMETER :: C6  = 0.720266520D-04
	DOUBLE PRECISION, PARAMETER :: C7  = 0.955528842D-03
	DOUBLE PRECISION, PARAMETER :: C8  = 0.101431553D-01
	DOUBLE PRECISION, PARAMETER :: C9  = 0.738522953D-01
	DOUBLE PRECISION, PARAMETER :: C10 = 0.338450368D+00
	DOUBLE PRECISION, PARAMETER :: C11 = 0.879937801D+00


	DOUBLE PRECISION, PARAMETER :: PI=3.1415926535898D+00	!f(circle) -> (circumference/diameter)

	INTEGER, PARAMETER :: NMAX = 5							!NUMBER OF NONLINEAR PARAMETERS

	DOUBLE PRECISION :: BD,R_ENER,EMIN,SQNORMGRAD			!BOND DISTANCE/ENERGIES/|GRADIENT|**2
	DOUBLE PRECISION :: EZ1,EZ2								!NUCLEAR GEOMETRY VARIABLES


!MICROSOFT PROPRIETARY OBJECT!
!-----------------------------------------------------------
!...this block must be removed for generalized porting to
!	other platforms...
!-----------------------------------------------------------
	DOUBLE PRECISION :: ELAPSED_TIME
!-----------------------------------------------------------
!                              MICROSOFT PROPRIETARY OBJECT!

	INTEGER ::	M											!NUMBER OF BASIS FCNS

END MODULE scalars_to_be_shared
