!      PROJECT TITLE:  correlated gaussian geminals for molecular hydrogen
!MAIN EXECUTION FILE:  cgauss.F90
!
!   Program FILENAME:  cgauss.exe
!        Lead Author:  D. Gilmore
!    Other Author(s):  L. Adamowicz, D. Kinghorn
!
!        LAST EDITED:  22may96:wed/17:26
!

!NOTES FOR FUTURE EDITS:
!===================================================================================================
!	1.  unnecessary sharing of modules should be eliminated
!
!	2.	elimination of unnecessary matrix storage space (or reworking)...
!		Kinetic,Potential,OVV1-OVV4
!
!	3.	find out what TSTORE?,?STORE,DSS,Q,Y does...
!===================================================================================================
!
!

!-------------------------------------------------------&   !=======================================
MODULE	matrices_to_be_shared

	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::	&
							S,H,CSTORE,TSTORE1,TSTORE2,	&
							TSTORE3,RZ1STORE,RZ2STORE,	&
							Kinetic,Potential,			&
							OVV1,OVV2,OVV3,OVV4
														

	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::	&
							DS,DH,DSS,GRADM
	
	DOUBLE PRECISION, DIMENSION(2,2) :: Q,Y

END MODULE matrices_to_be_shared
