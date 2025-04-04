	SUBROUTINE TESTDOVER(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,
     1                    S,H,CO)
C     -----------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /VALUE/ EMIN
	COMMON /OVER/ DETS(50),L,K,CC(50)
	COMMON /LC/RC(50)
      PARAMETER(MMAX=50)
      DIMENSION 
     1         AA1(MMAX),RX1(MMAX),RY1(MMAX),RZ1(MMAX),
     2         AA2(MMAX),RX2(MMAX),RY2(MMAX),RZ2(MMAX),
     3         B12(MMAX)
      DIMENSION S(MMAX,MMAX),H(MMAX,MMAX),CO(MMAX) 
      DATA ONE,TWO,FOUR,SIX,EIGHT,HALF/1.,2.,4.,6.,8.,.5/
C     -----------------------------------------------------	  
      CALL GEOM(EX1,EY1,EZ1,EX2,EY2,EZ2)
C     -----------------------------------------------------	  
C	DO L=1,M
C		ALK1=AA1(L)+AA1(L)
C		ALK2=AA2(L)+AA2(L)                                                       
C		BLK12=B12(L)+B12(L)  
C		DETS(L)=ALK1*ALK2+((ALK1+ALK2)*BLK12)
C	    CC(L)=8*AA1(L)*AA2(L)*B12(L)*((RZ1(L)-RZ2(L))**2)/DETS(L)
C	END DO
C
C      DO L=1,M
C          DO K=1,M
C     -----------------------------------------------------
C     --- Calculate OVERLAP INTEGRAL ---
C     -----------------------------------------------------

	
              CALL OVERLP(RAx,RAy,RAz,AA1,
     1            RCx,RCy,RCz,AA2,                      
     2            RBx,RBy,RBz,AA3,
     3            RDx,RDy,RDz,AA4,
     4            AA5,AA6,OV1,KLM)

	AA1=AA1+1.0D-09

              CALL OVERLP(RAx,RAy,RAz,AA1,
     1            RCx,RCy,RCz,AA2,                      
     2            RBx,RBy,RBz,AA3,
     3            RDx,RDy,RDz,AA4,
     4            AA5,AA6,OV2,KLM)

	OV3=(OV2-OV1)/1.0D-09
	WRITE(6,*)OV3
	STOP
	RETURN
	END
