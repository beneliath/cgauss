C     Originallly Coded:  Z. Zhang (...poorly)
C           Modified By:  D. Gilmore ---> 19OCT95
C     -----------------------------------------------------
C     L. ADAMOWICZ RESEARCH / UNIVERSITY OF ARIZONA / CHEM.
C     -----------------------------------------------------
      SUBROUTINE ENERGY(M,H,S,EMIN,CO,KLM)
      IMPLICIT REAL*8 (A-H,O-Z)
	INTEGER M
C 
C     N=number of particles
C     M= number of functions used in the expansion
C     IPER=number of permutations
C     NP = N*(N-1)/2 = number of interacting pairs
C     MN = M*N = number of orbital exponents
C     MNP = M*N*(N-1)/2 = number of correlation exponents
C
      PARAMETER ( NMAX = 5 )
C      PARAMETER ( MMAX = 3 )
C      PARAMETER ( NPMAX = NMAX*(NMAX-1)/2 )
C      PARAMETER ( MNMAX = MMAX*NMAX )
C 
      REAL*8, DIMENSION(M) ::	CO
	REAL*8, DIMENSION(M,M) :: S,H
	REAL*8 SNORM[ALLOCATABLE](:,:)

C      DIMENSION  SNORM(MMAX,MMAX)
      DIMENSION AA(66000),BIG(256),JB(256),TEMP(256),TP(66000)
      DIMENSION UT(256,256),U(256,256),TEM(256,256),SS(256,256)
      DIMENSION TEM1(256,256),TEM2(256,256),TEM3(256,256)
      DIMENSION SHALF(256,256),SNHALF(256,256),CHECKS(256,256)
      DIMENSION HH(256,256),XY(256,256),XT(256,256),W(256,256)
      DIMENSION CT(256,256),CHECKH(256,256),C(256,256),AB(256)
C
C      COMMON /MOLPARAM/ N, M, NP 
C      COMMON /FFF/ S, H
c      COMMON /NORM/ SNORM
	ALLOCATE(SNORM(M,M))
	MMAX=M
	MNMAX=M*NMAX
	NPMAX=NMAX*(NMAX-1)/2
C
      IF(KLM.EQ.5) GOTO 1414
C
C     PRINT-OUT THE S & H MATRICES (UNNORMALIZED):
C     Uncomment all code between the (===) lines to execute...
C     ===============================================================
C      CALL BIOUT(' ')
C      CALL BIOUT(' Overlap (S) matrix: ')
C      CALL BIOUT(' ------------------- ')
C          CALL TAB(S,M,M,MMAX,MMAX)
C      CALL BIOUT(' ')
C      CALL BIOUT(' Hamiltonian (H) matrix: ')
C      CALL BIOUT(' ----------------------- ')
C          CALL TAB(H,M,M,MMAX,MMAX)
C      CALL BIOUT(' ')
C     ===============================================================
C      STOP
C
C     ---------------------------------------------------------------
C     NORMALIZATION OF S & H MATRICES:
C     ---------------------------------------------------------------
C     Copy the contents of matrix S into matrix SNORM...
c      CALL MATRIXCOPY(M,S,MMAX,SNORM,MMAX)
C
C     Normalization of the S & H matrix elements...
c      DO 59 I=1,M
c      DO 59 J=1,M
c      S(I,J) = S(I,J)/( DSQRT(SNORM(I,I))*DSQRT(SNORM(J,J)) )
c      H(I,J) = H(I,J)/( DSQRT(SNORM(I,I))*DSQRT(SNORM(J,J)) )
C      S(I,J) = S(I,J)/( (SNORM(I,I)**(1/2))*(SNORM(J,J)**(1/2)) )
C      H(I,J) = H(I,J)/( (SNORM(I,I)**(1/2))*(SNORM(J,J)**(1/2)) )
c 59   CONTINUE
C
C
C     PRINT-OUTS OF NORMALIZED S & H MATRICES:
C     Uncomment all code between the (===) lines to execute...
C     ===============================================================
C      CALL BIOUT(' ')
C      CALL BIOUT(' Normalized Overlap (S) Matrix: ')
C      CALL BIOUT(' ------------------------------ ')
C         CALL TAB(S,M,M,MMAX,MMAX)
C      CALL BIOUT(' ')
C      CALL BIOUT(' Normalized Hamiltonian (H) Matrix: ')
C      CALL BIOUT(' ---------------------------------- ')
C         CALL TAB(H,M,M,MMAX,MMAX)
C      CALL BIOUT(' ')
C     ===============================================================
C      STOP
C
C
C     NNP1D2=number of lower triangular elements for YACOBI
      NNP1D2=M*(M+1)/2
C
C     ---------------------------------------------------------------
C     DIAGNOLIZE THE S MATRIX
C     ---------------------------------------------------------------
C     Define the one dimensional array AA(II) as the S matrix for YACOBI...
      CALL VECH (M,S,MMAX,AA)
C
C
C     Call YACOBI to get the TRANSPOSE of the ROTATION MATRIX...
      CALL YACOBI(AA,TP,M,NNP1D2,BIG,JB)
C
C     Construct the transpose matrix UT FROM the vector TP...
      CALL DEVEC (M,TP,UT,256)
C
C
C     PRINT-OUT THE TRANSPOSE ROTATION MATRIX:
C     Uncomment all code between the (===) lines to execute...
C     ===============================================================
C      CALL BIOUT(' ')
C      CALL BIOUT(' TRANSPOSE matrix ')
C      CALL BIOUT(' ---------------------- ')
C         CALL TAB(UT,M,M,100,100)
C     ===============================================================
C      STOP
C
C
C     Construct the regular rotation matrix from its transpose (UT)...
      CALL TRANSPOSE(M,UT,256,U,256)
C
C
C     SHOW THAT THE ROTATION IS A UNITARY TRANSFORMATION:
C     Uncomment all code between the (===) lines to execute...
C     ===============================================================
C      CALL BIOUT(' ')
C      CALL BIOUT(' Unitary Transformation?: ')
C      CALL BIOUT(' ------------------------ ' )
C          CALL MATRIXM (M,U,100,UT,100,TEM3,100)
C          CALL TAB(TEM3,M,M,100,100)
C      CALL BIOUT(' ')
C     ===============================================================
C      STOP
C
C
C     Perform the operations to diagonalize S
      CALL MATRIXM(M,UT,256,S,MMAX,TEM,256)
      CALL MATRIXM(M,TEM,256,U,256,SS,256)

C
C
C     PRINT-OUT THE DIAGONALIZED S MATRIX:
C     Uncomment all code between the (===) lines to execute...
C     ===============================================================
C      CALL BIOUT(' ')
C      CALL BIOUT(' Diagonalized (S) matrix: ')
C      CALL BIOUT(' ----------------------- ' )
C          CALL TAB(SS,M,M,100,100)
C      CALL BIOUT(' ')
C     ===============================================================
C      STOP
C
C
C     ---------------------------------------------------------------
C     DIAGONALIZE THE H MATRIX:  PART I
C     ---------------------------------------------------------------
C     Create S**(1/2) & S**(-1/2) matrices...
C
C
      DO 3130 I=1,M
          DO 3130 J=1,M
              SHALF(I,J)=0.0D+00
3130          CONTINUE
          DO 3140  II=1,M
C      WRITE(6,*) 'SS (',II,',',II,') = ',SS(II,II)
C      WRITE(7,*) 'SS (',II,',',II,') = ',SS(II,II)
              SHALF(II,II)=DSQRT(SS(II,II))
C              SHALF(II,II)=(SS(II,II)**(1/2))
3140  CONTINUE
C
      CALL MATRIXM(M,U,256,SHALF,256,TEM,256)
C
      CALL MATRIXM(M,TEM,256,UT,256,SHALF,256)
C
C
C     PRINT-OUT THE S**(1/2) MATRIX:
C     Uncomment all code between the (===) lines to execute...
C     ===============================================================
C      CALL BIOUT(' ')
C      CALL BIOUT(' S**(1/2) matrix: ')
C      CALL BIOUT(' ---------------- ')
C          CALL TAB(SHALF,M,M,100,100)
C      CALL BIOUT(' ')
C     ===============================================================
C      STOP
C
C
C     Blank the values in matrix SNHALF -->0.0D+00
      CALL BLANKM(M,SNHALF,256)
C
C
C     Construct the S**(-1/2) matrix...
      DO 3160 II=1,M
          SNHALF(II,II)=DSQRT(SS(II,II))
C          SNHALF(II,II)=(SS(II,II)**(1/2))
          SNHALF(II,II)=1.0D+00/SNHALF(II,II)
3160  CONTINUE
C
      CALL MATRIXM(M,U,256,SNHALF,256,TEM,256)
C
C
      CALL MATRIXM(M,TEM,256,UT,256,SNHALF,256)
C
C
C     PRINT-OUT THE S**(-1/2) MATRIX:
C     Uncomment all code between the (===) lines to execute...
C     ===============================================================
C      CALL BIOUT(' ')
C      CALL BIOUT(' S**(-1/2) matrix: ')
C      CALL BIOUT(' ---------------- ')
C          CALL TAB(SNHALF,M,M,100,100)
C      CALL BIOUT(' ')
C     ===============================================================
C      STOP
C
C
C     CHECK FOR IDENTITY MATRIX:
C     Uncomment all code between the (===) lines to execute...
C     ===============================================================
C      CALL MATRIXM(M,SNHALF,100,SHALF,100,CHECKS,100)
C      CALL BIOUT(' ')
C      CALL BIOUT(' S**(-1/2) * S**(1/2) = I CHECK:')
C      CALL BIOUT(' -------------------------------')
C      CALL TAB(CHECKS,M,M,100,100)
C     ===============================================================
C      STOP 
C
C
C     CREATE A NEW MATRIX HH=S**(-1/2) * H * S**(-1/2)
C
C
      CALL MATRIXM(M,SNHALF,256,H,MMAX,TEM,256)
C
C
      CALL MATRIXM(M,TEM,256,SNHALF,256,HH,256)
C
C
C     PRINT-OUT THE HH MATRIX:
C     Uncomment all code between the (===) lines to execute...
C     ===============================================================
C      CALL BIOUT(' ')
C      CALL BIOUT(' HH matrix: ')
C      CALL BIOUT(' ---------- ')
C          CALL TAB(HH,M,M,100,100)
C      CALL BIOUT(' ')
C     ===============================================================
C      STOP
C
C
C     DIAGNOLIZE HH MATRIX TO GET EIGENVALUE W AND EIGENVECTOR C-----
C
C
C     Define the one dimensional array AA(II) as the HH matrix for YABOCI...
      CALL VECH (M,HH,256,AA)
C
C
C     Call YACOBI to get the TRANSPOSE of the ROTATION MATRIX...
      CALL YACOBI(AA,TP,M,NNP1D2,BIG,JB)
C
C
C     Construct the transpose matrix XT FROM the vector TP...
      CALL DEVEC (M,TP,XT,256)
C
C
C     PRINT-OUT THE TRANSPOSE ROTATION MATRIX:
C     Uncomment all code between the (===) lines to execute...
C     ===============================================================
C      CALL BIOUT(' ')
C      CALL BIOUT(' TRANSPOSE matrix ')
C      CALL BIOUT(' ---------------------- ')
C         CALL TAB(XT,M,M,100,100)
C     ===============================================================
C      STOP
C
C
C     Construct the regular rotation matrix from its transpose (XT)...
      CALL TRANSPOSE(M,XT,256,XY,256)
C
C
C     SHOW THAT THE ROTATION IS A UNITARY TRANSFORMATION:
C     Uncomment all code between the (===) lines to execute...
C     ===============================================================
C      CALL BIOUT(' ')
C      CALL BIOUT(' Unitary Transformation?: ')
C      CALL BIOUT(' ------------------------ ' )
C          CALL MATRIXM (M,XY,100,XT,100,TEM3,100)
C          CALL TAB(TEM3,M,M,100,100)
C      CALL BIOUT(' ')
C     ===============================================================
C      STOP
C
C
C     CHECK FOR IDENTITY MATRIX:
C     Uncomment all code between the (===) lines to execute...
C     ===============================================================
C      CALL MATRIXM(M,SNHALF,100,XY,100,TEM,100)
C      CALL TRANSPOSE(M,TEM,100,TEM1,100)
C      CALL MATRIXM(M,TEM1,100,S,MMAX,TEM2,100)
C      CALL MATRIXM(M,TEM2,100,TEM,100,TEM3,100)
C      CALL BIOUT(' ')
C      CALL BIOUT(' (UT * ST**(-1/2)) * S * (S**(-1/2) * U) = I CHECK:')
C      CALL BIOUT(' -------------------------------------------------')
C      CALL TAB(TEM3,M,M,100,100)
C     ===============================================================
C      STOP 
C
C
C     Perform the operations to diagonalize HH
      CALL MATRIXM(M,XT,256,HH,256,TEM,256)
      CALL MATRIXM(M,TEM,256,XY,256,W,256)
C
C
C     PRINT-OUT THE DIAGONALIZED W MATRIX:
C     Uncomment all code between the (===) lines to execute...
C     ===============================================================
C      CALL BIOUT(' ')
C      CALL BIOUT(' Diagonalized (W) matrix: ')
C      CALL BIOUT(' ----------------------- ' )
C          CALL TAB(W,M,M,100,100)
C      CALL BIOUT(' ')
C     ===============================================================
C      STOP
C
C
C     Pull out the Ground State Energy from the Eigenvalues...
      EMIN=W(1,1)
      IX=1
C
      DO 3211 I=1,M
          CC=0.0D+00
          CC=W(I,I)
          IF(CC .LT. EMIN) THEN
              EMIN=CC
              IX=I
          ENDIF
3211  CONTINUE
C
C
C     Print out the lowest eigenvalue = Ground State Energy
C      CALL BIOUT(' ')
C      WRITE(7,1856) EMIN
C      WRITE(6,1856) EMIN
C      RJUICE=EMIN
C
C
C      STOP 
C















C ***   FIND EIGENVECTOR C
      CALL MATRIXM(M,SNHALF,256,XY,256,C,256)
      DO 600  I=1,M
      DO 600  J=1,M
      C(I,J)=0.0D+00
      DO 650  K=1,M
      C(I,J)=C(I,J)+SNHALF(I,K)*XY(K,J)
 650   CONTINUE
 600   CONTINUE
C      open(unit=6,access='append',form='formatted')
C      CALL BIOUT('CHECK THE EIGENVECTOR C')
C      DO 220  I=1,M
C      DO 220  J=1,M
C      WRITE(6,*) I,J,C(I,J)
C      WRITE(7,*) I,J,C(I,J)
C 220   CONTINUE
C     CLOSE(6)
      DO 226  I=1,M
      CO(I) = C(I,IX)
 226   CONTINUE
C      DO 230 I=1,M
C      DO 230 J=1,M
C      CT(I,J)=C(J,I)
C 230   CONTINUE
C      CALL MATRIXM(M,CT,100,S,100,TEM,100)
C      DO 700  I=1,M
C      DO 700  J=1,M
C      TEM(I,J)=0
C      DO 750  K=1,M
C      TEM(I,J)=TEM(I,J)+CT(I,K)*S(K,J)
C 750   CONTINUE
C 700   CONTINUE
C      CALL MATRIXM(M,TEM,100,C,100,CHECKS,100)
C      DO 100  I=1,M
C      DO 100  J=1,M
C      CHECKS(I,J)=0
C      DO 150  K=1,M
C      CHECKS(I,J)=CHECKS(I,J)+TEM(I,K)*C(K,J)
C 150   CONTINUE
C 100   CONTINUE
C      CALL TAB(CHECKS,M,M,100,100)
C     STOP
C       open(unit=6,access='append',form='formatted')
C     WRITE(6,*) 'CHECK S'
C     DO 240 I=1,M
C     DO 240 J=1,M
C     WRITE(6,*) CHECKS(I,J)
C240  CONTINUE
C     close(6)
C      DO 90 L=1,M
C      DO 90 K=1,M
C      SSS = SSS + CO(L)*CO(K)*S(L,K)
C90    CONTINUE 
C *** 
C      open (unit=6,access='append',form='formatted')
C      WRITE(6,*)' NORMALIZATION = ', SSS
C      WRITE(7,*)' NORMALIZATION = ', SSS
C      CALL TAB (C,M,M,100,100)
C      CALL BIOUT('******************************************')
C      CALL BIOUT('          EXPANSION COEFFICIENTS')
C      CALL BIOUT('==========================================')
C      DO 33 J = 1, M
	
C      WRITE(6,698) J, CO(J)
C      WRITE(7,698) J, CO(J)
C698   FORMAT(10X, 'CO(', I2,') =' D17.10)
C33    CONTINUE 
C      GO TO 1515
 1414  CONTINUE
C      EMIN = 10.0D+00
C      open (unit=6,access='append',form='formatted')
C      write (7,*) ' emin = 10.00 '
C      close (6)
1515  DEALLOCATE(SNORM)
	RETURN
C
C     FORMATS:
C     --------------------------------------------
 11   FORMAT(/,/)
1855  FORMAT (/)
1856  FORMAT (' GroundState Energy (CALCULATED) = ',F18.12,' au')
C     --------------------------------------------
      END
