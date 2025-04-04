      SUBROUTINE GRADEN(N,X,G)
C
      IMPLICIT REAL*8(A-H,O-Z)
	INTEGER DFLAG
C
      PARAMETER (MMAX=500)
C
      DIMENSION 
     1         AA1(MMAX),RX1(MMAX),RY1(MMAX),RZ1(MMAX),
     2         AA2(MMAX),RX2(MMAX),RY2(MMAX),RZ2(MMAX),
     3         B12(MMAX),EVEC(MMAX),S(MMAX,MMAX),H(MMAX,MMAX),CO(MMAX)  
	DIMENSION X(MMAX*5+MMAX),G(MMAX*5+MMAX)
	REAL*8 GRAD(MMAX*5+MMAX),T1(MMAX),T2(MMAX),T3(MMAX)
	CHARACTER RESTARTNAME*66
C
	COMMON /FINALGRAD/ GG
C
C
C	REWIND (8)
	NB=(N/6)
C
C
C
	K=1
	DO I=1,NB
			T1(I)=X(K)	
			T2(I)=X(K+1)	
			T3(I)=X(K+2)	
			RZ1(I)=X(K+3)	
			RZ2(I)=X(K+4)	
		K=K+5
	END DO
C
	DO I=1,NB
			EVEC(I)=X(NB*5+I)
	END DO
C
C
	M=NB
      DO I=1,M
	    AA1(I)=T1(I)**2 + T1(I)*T3(I)
          AA2(I)=T2(I)**2 + T3(I)**2 + T1(I)*T3(I)
          B12(I)=-T1(I)*T3(I)
	END DO
C
      DFLAG=1
	CALL HAMILTON(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,EVEC,S,H,
	-				CO,EMIN,GRAD,T1,T2,T3,DFLAG)
C
      DO I=1,M
	    AA1(I)=T1(I)**2 + T1(I)*T3(I)
          AA2(I)=T2(I)**2 + T3(I)**2 + T1(I)*T3(I)
          B12(I)=-T1(I)*T3(I)
	END DO
C
	DO I=1,NB*5+NB
		G(I)=GRAD(I)
	END DO
	GG=0.00D+00
	DO I=1,NB*5+NB
		GG=GG+GRAD(I)*GRAD(I)
	END DO
	WRITE(6,*)'ENERGY=',EMIN,'   --->   ||g||^2 = ',GG
	WRITE(7,*)'ENERGY=',EMIN,'   --->   ||g||^2 = ',GG
C
C
      J=1
	OPEN (8, FILE = 'Restart_Output.DAT', ACCESS = 'SEQUENTIAL',
     -		STATUS = 'UNKNOWN')
      WRITE(8,*) '     LinCoeff:.   AA1:        AA2:        B12:        
	-RZ1:       RZ2:'
      WRITE(8,*) 'FNC: =================================================
     -====================='
 706  IF (J .LE. M) THEN
          WRITE(8,1111) J,EVEC(J),AA1(J),AA2(J),B12(J),RZ1(J),RZ2(J)
          J=J+1

      GOTO 706
      END IF
      WRITE(8,*) '------------------------------------------------------
     1---------------------'
	WRITE(8,*)'ENERGY=',EMIN,'   --->   ||g||^2 = ',GG
	WRITE(8,*)
C
      J=1
      WRITE(8,*) ':BEGIN:'
	WRITE(8,*) M
 606  IF (J .LE. M) THEN
          WRITE(8,*) J,EVEC(J),T1(J),T2(J),T3(J),RZ1(J),RZ2(J)
          J=J+1

      GOTO 606
      END IF
      WRITE(8,*) ':END:'
	WRITE(8,*)
	CLOSE(8)
C
1111  FORMAT(I3,6D12.4)
C
C
      RETURN
      END
