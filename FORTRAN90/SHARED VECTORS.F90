!      PROJECT TITLE:  correlated gaussian geminals for molecular hydrogen
!MAIN EXECUTION FILE:  cgauss.F90
!
!   Program FILENAME:  cgauss.exe
!        Lead Author:  D. Gilmore
!    Other Author(s):  L. Adamowicz, D. Kinghorn
!
!        LAST EDITED:  22may96:wed/17:10
!

!NOTES FOR FUTURE EDITS:
!===================================================================================================
!	1.  unnecessary sharing of modules should be eliminated
!
!	2.	under cylindrical symmetry, vectors RX1,RY1,RX2,RY2 can be eliminated
!
!	3.	elimination of redundant linear coefficient vector CO
!===================================================================================================
!
!

!-------------------------------------------------------&   !=======================================
MODULE	vectors_to_be_shared

	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::		&
							AA1,RX1,RY1,RZ1,			&	!alpha1,...
							AA2,RX2,RY2,RZ2,B12,C,		&	!alpha2,...,correlation,linear coeff
							GRAD,CO,					&	!gradient vector,another lin coeff
							T1,T2,T3,GRADA,GRADC,		&	!storage for cholesky params & grad vec
							TSTORE1T,TSTORE2T,TSTORE3T,	&	!temporary storage for cholesky params...
							RZ1STORET,RZ2STORET,CSTORET		!coords & lin coeff.

END MODULE vectors_to_be_shared