SUBROUTINE Solve_2_by_2(MMAX,H,S,E)

	USE scalars_to_be_shared
!	USE vectors_to_be_shared			  

    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	IMPLICIT INTEGER (I-N)

	DOUBLE PRECISION, DIMENSION(2,2) :: H,S
	DOUBLE PRECISION, DIMENSION(4) :: AA,TP

	DOUBLE PRECISION, DIMENSION(2) ::	BIG,JB

	DOUBLE PRECISION, DIMENSION(2,2):: UT,U,TEM,SS,TEM1,	&
									TEM2,TEM3,SHALF,SNHALF,		&
									CHECKS,HH,XY,XT,W,CT,		&
									CHECKH,CA


!	MMAX=2
	MNMAX=MMAX*NMAX
	NPMAX=NMAX*(NMAX-1)/2


      NNP1D2=MMAX*(MMAX+1)/2


!DIAGNOLIZE THE S MATRIX
!---------------------------------------------------------------
      CALL VECH (MMAX,S,MMAX,AA)

      CALL YACOBI(AA,TP,2,NNP1D2,BIG,JB)

!Construct the transpose matrix UT FROM the vector TP...
      CALL DEVEC (2,TP,UT,MMAX)

!PRINT-OUT THE TRANSPOSE ROTATION MATRIX:
!===============================================================
!      CALL BIOUT(' ')
!      CALL BIOUT(' TRANSPOSE matrix ')
!      CALL BIOUT(' ---------------------- ')
!      CALL TAB(UT,M,M,100,100)
!===============================================================
!      STOP


!Construct the regular rotation matrix from its transpose (UT)...
!      CALL TRANSPOSE(M,UT,256,U,256)
		U=TRANSPOSE(UT)


!SHOW THAT THE ROTATION IS A UNITARY TRANSFORMATION:
!===============================================================
!      CALL BIOUT(' ')
!      CALL BIOUT(' Unitary Transformation?: ')
!      CALL BIOUT(' ------------------------ ' )
!      CALL MATRIXM (M,U,100,UT,100,TEM3,100)
!      CALL TAB(TEM3,M,M,100,100)
!      CALL BIOUT(' ')
!===============================================================
!      STOP


!Perform the operations to diagonalize S
!      CALL MATRIXM(M,UT,256,S,MMAX,TEM,256)
!      CALL MATRIXM(M,TEM,256,U,256,SS,256)

	TEM=MATMUL(UT,S)
	SS =MATMUL(TEM,U)


!PRINT-OUT THE DIAGONALIZED S MATRIX:
!===============================================================
!      CALL BIOUT(' ')
!      CALL BIOUT(' Diagonalized (S) matrix: ')
!      CALL BIOUT(' ----------------------- ' )
!      CALL TAB(SS,M,M,100,100)
!      CALL BIOUT(' ')
!===============================================================
!      STOP


!DIAGONALIZE THE H MATRIX:  PART I
!---------------------------------------------------------------
!Create S**(1/2) & S**(-1/2) matrices...

      CALL BLANKM(2,SHALF,MMAX)

!      DO I=1,M
!          DO J=1,M
!              SHALF(I,J)=0.0D+00
!		  END DO
!	  END DO

			DO II=1,2

              SHALF(II,II)=DSQRT(SS(II,II))

			END DO
!      CALL MATRIXM(M,U,256,SHALF,256,TEM,256)

		TEM=MATMUL(U,SHALF)

!      CALL MATRIXM(M,TEM,256,UT,256,SHALF,256)

		SHALF=MATMUL(TEM,UT)

!PRINT-OUT THE S**(1/2) MATRIX:
!===============================================================
!      CALL BIOUT(' ')
!      CALL BIOUT(' S**(1/2) matrix: ')
!      CALL BIOUT(' ---------------- ')
!      CALL TAB(SHALF,M,M,100,100)
!      CALL BIOUT(' ')
!===============================================================
!      STOP


!Blank the values in matrix SNHALF -->0.0D+00
      CALL BLANKM(2,SNHALF,MMAX)


!Construct the S**(-1/2) matrix...
      DO II=1,2
          SNHALF(II,II)=DSQRT(SS(II,II))
          SNHALF(II,II)=1.0D+00/SNHALF(II,II)
	  END DO

!      CALL MATRIXM(M,U,256,SNHALF,256,TEM,256)

		TEM=MATMUL(U,SNHALF)

!      CALL MATRIXM(M,TEM,256,UT,256,SNHALF,256)

		SNHALF=MATMUL(TEM,UT)

!PRINT-OUT THE S**(-1/2) MATRIX:
!===============================================================
!      CALL BIOUT(' ')
!      CALL BIOUT(' S**(-1/2) matrix: ')
!      CALL BIOUT(' ---------------- ')
!      CALL TAB(SNHALF,M,M,100,100)
!      CALL BIOUT(' ')
!===============================================================
!      STOP


!CHECK FOR IDENTITY MATRIX:
!===============================================================
!      CALL MATRIXM(M,SNHALF,100,SHALF,100,CHECKS,100)
!      CALL BIOUT(' ')
!      CALL BIOUT(' S**(-1/2) * S**(1/2) = I CHECK:')
!      CALL BIOUT(' -------------------------------')
!      CALL TAB(CHECKS,M,M,100,100)
!===============================================================
!      STOP 


!CREATE A NEW MATRIX HH=S**(-1/2) * H * S**(-1/2)
      
!	  CALL MATRIXM(M,SNHALF,256,H,MMAX,TEM,256)

		TEM=MATMUL(SNHALF,H)

!      CALL MATRIXM(M,TEM,256,SNHALF,256,HH,256)

		HH=MATMUL(TEM,SNHALF)

!PRINT-OUT THE HH MATRIX:
!===============================================================
!      CALL BIOUT(' ')
!      CALL BIOUT(' HH matrix: ')
!      CALL BIOUT(' ---------- ')
!      CALL TAB(HH,M,M,100,100)
!      CALL BIOUT(' ')
!===============================================================
!      STOP


!DIAGNOLIZE HH MATRIX TO GET EIGENVALUE W AND EIGENVECTOR C-----


!Define the one dimensional array AA(II) as the HH matrix for YABOCI...
      CALL VECH (2,HH,MMAX,AA)


      CALL YACOBI(AA,TP,2,NNP1D2,BIG,JB)


!Construct the transpose matrix XT FROM the vector TP...
      CALL DEVEC (2,TP,XT,MMAX)


!PRINT-OUT THE TRANSPOSE ROTATION MATRIX:
!===============================================================
!      CALL BIOUT(' ')
!      CALL BIOUT(' TRANSPOSE matrix ')
!      CALL BIOUT(' ---------------------- ')
!         CALL TAB(XT,M,M,100,100)
!===============================================================
!      STOP


!Construct the regular rotation matrix from its transpose (XT)...

!      CALL TRANSPOSE(M,XT,256,XY,256)

	XY=TRANSPOSE(XT)

!SHOW THAT THE ROTATION IS A UNITARY TRANSFORMATION:
!===============================================================
!      CALL BIOUT(' ')
!      CALL BIOUT(' Unitary Transformation?: ')
!      CALL BIOUT(' ------------------------ ' )
!          CALL MATRIXM (M,XY,100,XT,100,TEM3,100)
!          CALL TAB(TEM3,M,M,100,100)
!      CALL BIOUT(' ')
!===============================================================
!      STOP


!CHECK FOR IDENTITY MATRIX:
!===============================================================
!      CALL MATRIXM(M,SNHALF,100,XY,100,TEM,100)
!      CALL TRANSPOSE(M,TEM,100,TEM1,100)
!      CALL MATRIXM(M,TEM1,100,S,MMAX,TEM2,100)
!      CALL MATRIXM(M,TEM2,100,TEM,100,TEM3,100)
!      CALL BIOUT(' ')
!      CALL BIOUT(' (UT * ST**(-1/2)) * S * (S**(-1/2) * U) = I CHECK:')
!      CALL BIOUT(' -------------------------------------------------')
!      CALL TAB(TEM3,M,M,100,100)
!===============================================================
!      STOP 


!Perform the operations to diagonalize HH

!      CALL MATRIXM(M,XT,256,HH,256,TEM,256)

			TEM=MATMUL(XT,HH)

!      CALL MATRIXM(M,TEM,256,XY,256,W,256)

			W=MATMUL(TEM,XY)

!PRINT-OUT THE DIAGONALIZED W MATRIX:
!===============================================================
!      CALL BIOUT(' ')
!      CALL BIOUT(' Diagonalized (W) matrix: ')
!      CALL BIOUT(' ----------------------- ' )
!      CALL TAB(W,M,M,100,100)
!      CALL BIOUT(' ')
!===============================================================
!      STOP


!Pull out the Ground State Energy from the Eigenvalues...
      EMIN=W(1,1)
      IX=1

    DO I=1,2

          CC=ZERO

          CC=W(I,I)
!		  write(6,*)cc

          IF(CC .LT. EMIN) THEN
              EMIN=CC
              IX=I
          ENDIF

	END DO

!WRITE(6,*)'E=',E



!***   FIND EIGENVECTOR C

!      CALL MATRIXM(M,SNHALF,256,XY,256,CA,256)

			CA=MATMUL(SNHALF,XY)

	DO I=1,2
		DO J=1,2

      CA(I,J)=ZERO

		DO K=1,2
			CA(I,J)=CA(I,J)+SNHALF(I,K)*XY(K,J)
		END DO

		END DO
	END DO


!		DO I=1,2
!			CO(I) = CA(I,IX)
!		END DO


E=EMIN

END SUBROUTINE Solve_2_by_2
