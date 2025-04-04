!      PROJECT TITLE:  correlated gaussian geminals for molecular hydrogen
!MAIN EXECUTION FILE:  cgauss.F90
!
!   Program FILENAME:  cgauss.exe
!        Lead Author:  D. Gilmore
!    Other Author(s):  L. Adamowicz, D. Kinghorn
!
!        LAST EDITED:  23may96:wed/07:58
!

!NOTES FOR FUTURE EDITS:
!===================================================================================================
!	1.  unnecessary sharing of modules should be eliminated
!===================================================================================================
!
!

!-------------------------------------------------------&   !=======================================
SUBROUTINE Solve_Eigen_Value_Problem
!-----------------------------------------------------------
	USE scalars_to_be_shared			  
	USE	vectors_to_be_shared			  
	USE	matrices_to_be_shared			  

    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	IMPLICIT INTEGER (I-N)

	DOUBLE PRECISION, DIMENSION  (:), ALLOCATABLE ::		&
		AA,BIG,JB,TP

	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::	&
		UT,U,TEM,SS,SHALF,SNHALF,HH,XY,XT,W,CA
!-----------------------------------------------------------
	ALLOCATE( AA(M*M),BIG(M),JB(M),TP(M*M),UT(M,M),		&
		U(M,M),TEM(M,M),SS(M,M),SHALF(M,M),SNHALF(M,M),	&
		HH(M,M),XY(M,M),XT(M,M),W(M,M),CA(M,M)	)

	MMAX=M
	MNMAX=M*NMAX
	NPMAX=NMAX*(NMAX-1)/2

    NNP1D2=M*(M+1)/2

		CALL VECH (M,S,MMAX,AA)								!DIAGNOLIZE THE S MATRIX
		CALL YACOBI(AA,TP,M,NNP1D2,BIG,JB)

		CALL DEVEC (M,TP,UT,MMAX)							!transpose matrix UT -> the vector TP

		U=TRANSPOSE(UT)										!rotation matrix -> transpose (UT)


		TEM=MATMUL(UT,S)									!diagonalization of S...
		SS =MATMUL(TEM,U)
															!DIAGONALIZE THE H MATRIX:  PART I
		CALL BLANKM(M,SHALF,MMAX)							!Create S**(1/2) & S**(-1/2) matrices...

			DO II=1,M
				SHALF(II,II)=DSQRT(SS(II,II))
			END DO

		TEM=MATMUL(U,SHALF)
		SHALF=MATMUL(TEM,UT)

		CALL BLANKM(M,SNHALF,MMAX)							!Blank values in matrix SNHALF -> ZERO

		DO II=1,M											!Construct the S**(-1/2) matrix...
			SNHALF(II,II)=DSQRT(SS(II,II))
			SNHALF(II,II)=1.0D+00/SNHALF(II,II)
		END DO

		TEM=MATMUL(U,SNHALF)
		SNHALF=MATMUL(TEM,UT)
															!CREATE MATRIX HH=S^(-1/2)*H*S^(-1/2)...
		TEM=MATMUL(SNHALF,H)
		HH=MATMUL(TEM,SNHALF)
															!DIAGNOLIZE HH: eigenval W & eigenvec C

		CALL VECH (M,HH,MMAX,AA)							!1-D array AA(II) = HH matrix for YABOCI

		CALL YACOBI(AA,TP,M,NNP1D2,BIG,JB)

		CALL DEVEC (M,TP,XT,MMAX)							!transpose XT <- vector TP

		XY=TRANSPOSE(XT)									!rotation matrix from its transpose (XT)

															!Perform operations to diagonalize HH
		TEM=MATMUL(XT,HH)
		W=MATMUL(TEM,XY)

		EMIN=W(1,1)											!Ground State Energy <- Eigenvalues...
		IX=1

		DO I=1,M
			CC=ZERO
			CC=W(I,I)
			IF(CC .LT. EMIN) THEN
				EMIN=CC
				IX=I
			ENDIF
		END DO
															!***   FIND EIGENVECTOR C   ***
		CA=MATMUL(SNHALF,XY)

		DO I=1,M
			DO J=1,M
				CA(I,J)=ZERO
					DO K=1,M
						CA(I,J)=CA(I,J)+SNHALF(I,K)*	&
								XY(K,J)
					END DO
			END DO
		END DO

		DO I=1,M
			CO(I) = CA(I,IX)
		END DO

	DEALLOCATE( AA,BIG,JB,TP,UT,U,TEM,SS,SHALF,SNHALF,	&
		HH,XY,XT,W,CA )

END SUBROUTINE Solve_Eigen_Value_Problem
