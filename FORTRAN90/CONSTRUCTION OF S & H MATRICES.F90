!      PROJECT TITLE:  correlated gaussian geminals for molecular hydrogen
!MAIN EXECUTION FILE:  cgauss.F90
!
!   Program FILENAME:  cgauss.exe
!        Lead Author:  D. Gilmore
!    Other Author(s):  L. Adamowicz, D. Kinghorn
!
!        LAST EDITED:  23may96:wed/07:31
!

!NOTES FOR FUTURE EDITS:
!===================================================================================================
!	1.  unnecessary sharing of modules should be eliminated
!===================================================================================================
!
!

!-------------------------------------------------------&   !=======================================
SUBROUTINE Construct_S_and_H_Matrices
!-----------------------------------------------------------
	USE scalars_to_be_shared
	USE	matrices_to_be_shared

	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	IMPLICIT          INTEGER    (I-N)
!-----------------------------------------------------------
		CALL Calculate_Hamiltonian							!only the lower triangles are calcd
															!because S & H are hermitian...
		DO L=1,M
			DO K=1,L										!...so load upper triangle from calcd
					H(K,L)        = H(L,K)					!lower triangle
					S(K,L)        = S(L,K)
			END DO
		END DO
!Do i=1,M
!do j=1,M
!		write(6,*) i,S(i,i)
!End do
!end do

!stop
!	OPEN (5,FILE='S.dat')

!Do i=1,M
!		write(5,*) (S(i,j),j=1,M)
!write(5,*)''
!		write(6,*) (S(i,j),j=1,M)
!write(6,*)'line break'
!End do
!Close(5)

!Stop
!-----------------------------------------------------------!TESTING AREA   :BEGIN:
!		CALL BIOUT('S-MATRIX:')								
!			CALL TAB(S,M,M,M,M)								!printout overlap matrix
!		CALL BIOUT(' ')
!			                  
!		CALL BIOUT('H-MATRIX:')
!			CALL TAB(H,M,M,M,M)                  			!printout hamiltonian matrix
!		CALL BIOUT(' ')
!STOP
!-----------------------------------------------------------!TESTING AREA     :END:

END SUBROUTINE Construct_S_and_H_Matrices