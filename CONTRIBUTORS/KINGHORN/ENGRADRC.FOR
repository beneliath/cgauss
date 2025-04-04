      SUBROUTINE ENGRADRC( NXAC, XAC, ENG, GRAD)
c
c --- Donald B. Kinghorn
c     Dept of Chemistry
c     University of Arizona
c     Feb 24 1996
c
c     Last modified Feb 24 1996 DBK
      IMPLICIT REAL*8 (A-H,O-Z)

c
	COMMON /GRADIENT2/ EMIN,GRAND
      COMMON /PASSIT/ M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,EVEC,S,H,CO
	PARAMETER(MMAX=50)
C

      DIMENSION 
     1         AA1(MMAX),RX1(MMAX),RY1(MMAX),RZ1(MMAX),
     2         AA2(MMAX),RX2(MMAX),RY2(MMAX),RZ2(MMAX),
     3         B12(MMAX),EVEC(MMAX)
      DIMENSION S(MMAX,MMAX),H(MMAX,MMAX),CO(MMAX)

c     ..Scalar Arguments..    
      INTEGER NXAC
      REAL*8 ENG
c     ..Array Arguments..
      REAL*8 XAC(*),GRAD(*)
	REAL*8 GRAND(MMAX*5+MMAX)
	EXTERNAL HAMILTON
	
c =====================================================================
c Purpose
c =======
c
c Correlated Gaussian energy and gradient calculation using the Rayliegh Quotient
c This subroutine is set up to be called from the optimization 
c subroutine TN
c 
c This subroutine computes the energy and gradient at a given point XAC  
c 
c Arguments
c =========
c
c NXAC      (input) INTEGER
c         Size of vectors XAC and GRAD  (NXAC = #of basis functions,(NB) times
c                                             length of vech[LK], (NN) + NB)     
c
c XAC       (input) DOUBLE PRECISION vector, dimension (NXAC)
c         XAC = {vech(LK): K=1..NB : EVEC}
c         Independent variables for energy functional
c
c ENG     (output) DOUBLE PRECISION energy at point X
c         ENG is the Rayleigh qoutient
c
c GRAD    (output) DOUBLE PRECISION gradient at the point XAC 
c              
c
c =====================================================================
c 
c     ..
c     .. Parameters ..
*********************************************************************
*********************************************************************
c     These are the parameters that limit the size of problem
c     that can be solved with this code.
c     CHANGE THESE TO INCREASE PROBLEM SIZE LIMITS
c     ********LOOK IN THE CALLING ROUTIEN FOR THESE PARAMETERS TOO***
c     #############LOOK AT THE COMMON BLOCKS TOO!###############
*********************************************************************
	K=1
	DO I=1,M
		AA1(I)=XAC(K)
		AA2(I)=XAC(K+1)
		B12(I)=XAC(K+2)
		RZ1(I)=XAC(K+3)
		RZ2(I)=XAC(K+4)
	K=K+5
	END DO
	DO I=1,M
		EVEC(I)=XAC(M*5+I)
	END DO
      CALL HAMILTON(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,EVEC,S,H,CO)
	ENG=EMIN
	DO I=1,M*5+M
		GRAD(I)=GRAND(I)
	END DO
c
c     END OF ENGRADRC
c
      END           
      
