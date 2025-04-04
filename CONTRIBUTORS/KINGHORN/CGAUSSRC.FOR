                 PROGRAM CGAUSSRC
c
c --- Donald B. Kinghorn
c     Dept of Chemistry
c     Washington State University
c     April 24 1995
c
c =====================================================================
c Purpose
c =======
c
c Test and debug program for the ENGRAD subroutine
c The energy and gradient are computed after the hamiltonian
c overlap and gradient component matrices are built.
c
c =====================================================================
	IMPLICIT NONE
c     ..
c     .. Parameters ..
*********************************************************************
*********************************************************************
c     These are the parameters that limit the size of problem
c     that can be solved with this code.
c     CHANGE THESE TO INCREASE PROBLEM SIZE LIMITS
c     #############LOOK AT THE COMMON BLOCKS TOO!###############
c     ******THESE PARAMETERS NEED TO BE SET IN ENGRAD.FOR ALSO*******
*********************************************************************
c     NL is the maximum number of pseudo particles i.e. the 
c     maximum size of the exponent matrices.
      INTEGER             NL
      PARAMETER          (NL = 2)
c     NNMX is NL*(NL+1)/2      
      INTEGER             NNMX
      PARAMETER          (NNMX = 3)
c     NBMX is the maximum number of basis functions      
      INTEGER             NBMX
      PARAMETER          (NBMX = 256)
c     NXMX is the maximum length of the exponent vector NXMX=(NNMX*NBMX)      
      INTEGER             NXMX
      PARAMETER          (NXMX = 768)
c     NXACMX is the maximum length of the exponent + evec vector NXACMX=(NNMX*NBMX+NBMX)
      INTEGER             NXACMX
      PARAMETER          (NXACMX = 1024)
c     NMXSYMM is maximum number of terms in symmetry projector      
      INTEGER             NMXSYMM
      PARAMETER          (NMXSYMM = 2)
*********************************************************************
*********************************************************************    
c     ..
c     .. Parameter..
c     LW is length of workspace for TN  NLW = 14*NXACMX = 35840
      INTEGER             NLW
      PARAMETER          (NLW = 14336) 
	
      DOUBLE PRECISION   ZERO
      PARAMETER          (ZERO = 0.0D0)
	DOUBLE PRECISION   ONE
      PARAMETER          (ONE = 1.0D0)

	           
c     ..
c     .. Local Scalars ..
      INTEGER             NX, NXAC, I, J, K, IERROR, LW, CNT,
     &                    STOPNB, STEPNB, NEWNX, MAXIT, MAXFUN, MSGLVL 
      DOUBLE PRECISION    ENG, OLDENG, ALPHA, OLDX(NXMX)
      DOUBLE PRECISION    ETA, ACCRCY, XTOL, STEPMX, DSQRT,
     &                    VIRIAL, SUMT, SUMV  
      DOUBLE PRECISION	RN
c     ..
c     .. Local Arrays .. 
      DOUBLE PRECISION    X(NXMX), XAC(NXACMX), GRAD(NXACMX), W(NLW)
     
c     ..
c     .. External Functions ..     

c     ..
c     .. External Subroutines ..
      EXTERNAL            ENGRADRC, DRNNOA
c     ..
c     .. Common blocks .. 
*********************************************************************
**********************COMMON BLOCKS**********************************
*********************************************************************
c N       number of particles
c NB      number of basis functions
c NN      N*(N+1)/2
c NSYMM   number of terms in symmetry projection
c WEV     "which eigen value"
c SYMMC   symmetry coefs
c SYMM    symmetry projection matrices
c MASS    mass matrix for kinetic energy
c CHRG    charge vector for potential energy
c ABSTOL  convergence tolarence for the eigenvalue
      INTEGER             N, NB, NN, NSYMM
      DOUBLE PRECISION    SYMMC(NMXSYMM), SYMM(NL,NL,NMXSYMM),
     &                    MASS(NL,NL), CHRG(NNMX),
     &                    EVEC(NBMX),
     &                    T(NBMX,NBMX), V(NBMX,NBMX), SCALE
      COMMON /DATA1/ N, NB, NN, NSYMM
	COMMON /DATA2/ MASS, CHRG, SYMMC,SYMM
      COMMON /EIG/ SCALE, EVEC, T, V
*********************************************************************
*********************************************************************
********************************************************************* 
c     ..
c     .. Executable Statements .. 

      LW = NLW
c                           
c INPUT:
c
c     2 blank reads for header and N label
      READ(*,*)
      READ(*,*)
c     read N
      READ(*,*) N
c     Print NB
      WRITE(*,*) 'N IS:'
      WRITE(*,*) N
                  
c     balnk read for NB label
      READ(*,*)
c     read NB
      READ(*,*) NB
c     Print NB
      WRITE(*,*) 'NB IS:'
      WRITE(*,*) NB
            
c     blank read for NN label
      READ(*,*)
c     read NN
      READ(*,*) NN
c     Print NN
      WRITE(*,*) 'NN IS:'
      WRITE(*,*) NN
            
c     blank read for WEV label ***NOT USED***
      READ(*,*)
c     read WEV
      READ(*,*) 
c     Print WEV
C      WRITE(*,*) 'WEV IS:'
C      WRITE(*,*) WEV
            
c     blank read for ABSTOL label	***NOT USED***
      READ(*,*)
c     read ABSTOL
      READ(*,*) 
c     Print ABSTOL
C      WRITE(*,*) 'ABSTOL IS:'
C      WRITE(*,*) ABSTOL
      
c     blank read for STEPMX label
      READ(*,*)
c     read STEPMX
      READ(*,*) STEPMX
c     Print STEPMX
      WRITE(*,*) 'STEPMX IS:'
      WRITE(*,*) STEPMX         

c     blank read for ENG label
      READ(*,*)
c     read ENG
      READ(*,*) ENG
c     Print ENG
C      WRITE(*,*) 'INITIAL GUESS FOR ENERGY IS:'
C      WRITE(*,*) ENG                                
c
c     compute NX
c
      NX = NB * NN
c     Print NX
      WRITE(*,*) 'NX IS:'
      WRITE(*,*) NX        
     
c     blank read for X label
      READ(*,*)
c     read X
      READ(*,*) ( X(I), I=1,NX )
c     Print X
C      WRITE(*,*) 'INITIAL GUESS FOR EXPONENTS, X IS:'
C      WRITE(*,*) ( X(I), I=1,NX )
C      WRITE(*,*)
C      WRITE(*,*)      
c
c     Read in the charge product vector CHRG
c
c     blank read for CHRG label
      READ(*,*)
c     read CHRG
      READ(*,*) ( CHRG(I), I=1,NN )
c     Print CHRG
C      WRITE(*,*) 'CHRG IS:'
C      WRITE(*,9901) ( CHRG(I), I=1,NN )
C      WRITE(*,*)
C      WRITE(*,*)       
c 
c     Read in the mass matrix MASS
c
c     blank read for MASS label
      READ(*,*)
      READ(*,*) ( ( MASS(I,J), J=1,N), I=1,N)
c     Print MASS
C      WRITE(*,*) 'MASS IS:'
C      DO I=1,N
C          DO J=1,N
C              WRITE(*,9901)  MASS(I,J)
C          END DO
C          WRITE(*,*)
C      END DO
C      WRITE(*,*)      
c
c     Read in symmetry information
c
c     balnk read for NSYMM label
      READ(*,*)
c     read NSYMM
      READ(*,*) NSYMM
c     Print NSYMM
      WRITE(*,*) 'NSYMM IS:'
      WRITE(*,*) NSYMM
                  
c     blank read for SYMMC label
      READ(*,*) 
c     Read in SYMMC
      READ(*,*) ( SYMMC(I), I=1,NSYMM )
      DO I=1,NSYMM
          SYMMC(I) = SYMMC(I)/NSYMM
      END DO
c     Print SYMMC
c      WRITE(*,*) 'SYMMC IS:'
c      WRITE(*,9901) ( SYMMC(I), I=1,NSYMM )
c      WRITE(*,*)
c      WRITE(*,*)             
c
c     blank read for SYMM label
      READ(*,*)      
c     Read each permutation matrix seperated by a blank read
      DO K=1,NSYMM
          READ(*,*)
          READ(*,*) ( ( SYMM(I,J,K), J=1,N), I=1,N)
      END DO         
c     Print SYMM
c      WRITE(*,*) 'SYMM IS:'
c      DO K=1,NSYMM
c          DO I=1,N
c              DO J=1,N
c                  WRITE(*,9901)  SYMM(I,J,K)
c              END DO
c              WRITE(*,*)
c          END DO
c          WRITE(*,*)
c      END DO
c      WRITE(*,*)  

c     blank read for EVEC label
      READ(*,*)
c     read EVEC
      READ(*,*) ( EVEC(I), I=1,NB )
C	DO I=1,NB
C		EVEC(I) = ONE/DBLE(NB)
C	END DO

c     Make XAC
      NXAC = NX + NB
      DO I=1,NX
		XAC(I) = X(I)
	END DO
	
      DO I=1,NB
		XAC(NX+I) = EVEC(I)
	END DO	 

C     CALL LMQN FOR THE MINIMIZATION
C
C ETA    - SEVERITY OF THE LINESEARCH
C MAXFUN - MAXIMUM ALLOWABLE NUMBER OF FUNCTION EVALUATIONS
C XTOL   - DESIRED ACCURACY FOR THE SOLUTION X*
C STEPMX - MAXIMUM ALLOWABLE STEP IN THE LINESEARCH
C ACCRCY - ACCURACY OF COMPUTED FUNCTION VALUES
C MSGLVL - DETERMINES QUANTITY OF PRINTED OUTPUT
C          0 = NONE, 1 = ONE LINE PER MAJOR ITERATION.
C MAXIT  - MAXIMUM NUMBER OF INNER ITERATIONS PER STEP
C
      STOPNB = 0
      STEPNB = 0

1113  CONTINUE
      
      MAXIT  = 25
      MAXFUN = 5000
      ETA    = .25D0
c      STEPMX = 1.D+1 MOVED THIS TO INPUT FILE
      ACCRCY = 2.220446049250313D-016
      XTOL   = DSQRT(ACCRCY)
      MSGLVL = 1      
            
      CALL LMQN (IERROR, NXAC, XAC, ENG, GRAD, W, LW, ENGRADRC, 
     &     MSGLVL, MAXIT, MAXFUN, ETA, STEPMX, ACCRCY, XTOL)

   
      
C      DO I=1,NB
C          IF(ABS(EVEC(I,1)) .GT. 1.0D0) THEN
C              WRITE(*,*) 'RESTART', I
C              DO J=1,NN                  
C                   CALL DRNNOA (1, RN)
C                  X((I-1)*NN+J) = X((I-1)*NN+J)*(RN+.75D0)
C                   X((I-1)*NN+J) = RN
C              END DO
C              GOTO 1113
C          END IF
C      END DO

C
C PRINT THE RESULTS
C
      IF (IERROR .NE. 0) WRITE(*,800) IERROR
C      IF (MSGLVL .GE. 1) WRITE(*,810)
C      IF (MSGLVL .GE. 1) WRITE(*,820) (I,X(I),I=1,NX)
c      
c     Print the energy, X,and eigenvector     
c

c     call engrad one more time to make sure eigenvector is correct
C      CALL ENGRAD(NX, X, ENG, GRAD)                      
c     Print the energy
     
      WRITE(*,*) 'THE ENERGY IS:'
      WRITE(*,*) ENG
c     Print GRAD
      WRITE(*,*) 'GRAD IS:'
      WRITE(*,9901) ( GRAD(I), I=1,NXAC )
      WRITE(*,*)
C      WRITE(*,*)
      
c     Print X
      WRITE(*,*) 'X IS:'
      WRITE(*,9901) ( XAC(I), I=1,NX )
      WRITE(*,*)
      WRITE(*,*)
      
c     Print EVEC 
      WRITE(*,*) 'EVEC IS:'
      WRITE(*,9901) ( XAC(I+NX), I=1,NB )
      WRITE(*,*)
      WRITE(*,*)

C*************************************************
C     COMPUTE THE SCALED ENERGY AND VIRIAL COEF
C
C     COMPLETE UPPER TRIANGLE OF T,V
      DO I=1,NB-1
          DO J=I+1,NB
              T(I,J) = T(J,I)
              V(I,J) = V(J,I)
          END DO
      END DO
C     COMPUTE C'TC, C'VC 
        
      DO I=1,NB-1
          DO J=I+1,NB
              T(I,J) = T(J,I)
              V(I,J) = V(J,I)
          END DO
      END DO
      
      SUMT = ZERO
      SUMV = ZERO
	DO I=1,NB
		DO J=1,NB
			SUMT = SUMT + EVEC(I)*T(I,J)*EVEC(J)
			SUMV = SUMV + EVEC(I)*V(I,J)*EVEC(J)
		END DO
	END DO
	SUMT = SUMT * SCALE
	SUMV = SUMV * SCALE

C     VIRIAL COEF
      VIRIAL = -SUMV/(2.0D0*SUMT)
C     SCALED ENERGY
      ENG = VIRIAL*(SUMV/2.0D0)
C      ENG = SCALE*ENG      
      
      WRITE(*,*) 'VIRIAL COEF IS: ', VIRIAL
      WRITE(*,*) 'SCALE cSc IS:', SCALE
      WRITE(*,*) 'SCALED ENERGY IS: ', ENG
      WRITE(*,*) 'KINETIC ENERGY IS: ', SUMT 
      WRITE(*,*) 'POTENTIAL ENERGY IS: ', SUMV 
C ***********************************************

C     ADD STEPNB NEW BASIS FCNS TO LAST RESULT AND RESTART

c     Print EVEC 
C      WRITE(*,*) 'EVEC IS:'
C      WRITE(*,9901) ( EVEC(I,1), I=1,NB )
C      WRITE(*,*)
C      WRITE(*,*)
            
      NB = NB + STEPNB
      WRITE(*,*) 'NB IS: ', NB
      IF ( NB .LE. STOPNB ) THEN         
          NEWNX = NX + NN*STEPNB
          DO I=(NX+1),NEWNX
			CALL DRNNOA (1, RN)
              X(I) =   RN
          END DO
          NX = NEWNX
          
          GOTO 1113
      END IF

C     DO A LITTLE LEVELING ON THE EXPONENTS      
C      WRITE(*,*)
C	WRITE(*,*) 'Opt on Rayleigh Quotient'
C	    MAXIT  = 10
C          MAXFUN = 25
C          ETA    = .25D0
C          ACCRCY = 1.D-16
C          XTOL   = DSQRT(ACCRCY)
C          MSGLVL = 1
C          CALL LMQN (IERROR, NX, X, ENG, GRAD, W, LW, EGFIXEV, 
C     &    MSGLVL, MAXIT, MAXFUN, ETA, STEPMX, ACCRCY, XTOL)
C	GOTO 1113
              
      
c     Printing Formats
 800  FORMAT(//,' ERROR CODE =', I3,/)
 810  FORMAT(10X, 'CURRENT SOLUTION IS ',/14X, 'I', 11X, 'X(I)')
 820  FORMAT(10X, I5, 2X, 1PD22.15)
 9901 FORMAT( F22.15)         
c
c     END OF CGAUSSRC
c
      END                       
