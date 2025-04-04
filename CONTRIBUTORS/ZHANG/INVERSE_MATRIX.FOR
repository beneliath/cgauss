      SUBROUTINE FCN(N,X,F)
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
	DIMENSION X(MMAX*5+MMAX),T1(MMAX),T2(MMAX),T3(MMAX)
C
C
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
C
      DO I=1,M
	    AA1(I)=T1(I)**2 + T1(I)*T3(I)
          AA2(I)=T2(I)**2 + T3(I)**2 + T1(I)*T3(I)
          B12(I)=-T1(I)*T3(I)
	END DO

C
	DFLAG=0
      CALL HAMILTON(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,EVEC,S,H,
	-				CO,EMIN,GRAD,T1,T2,T3,DFLAG)
C
C	J=1
C      CALL BIOUT('     LinCoeff:.   AA1:        AA2:        B12:        
C	1RZ1:       RZ2:')
C      CALL BIOUT('FNC: =================================================
C     1=====================')
C 606  IF (J .LE. M) THEN
C 1111 FORMAT(I3,6D12.4)
C          WRITE(6,1111) J,EVEC(J),AA1(J),AA2(J),B12(J),RZ1(J),RZ2(J)
C          WRITE(7,1111) J,EVEC(J),AA1(J),AA2(J),B12(J),RZ1(J),RZ2(J)
C          J=J+1
C      GOTO 606
C      END IF
C      CALL BIOUT('------------------------------------------------------
C     1---------------------')
C
C	DO I=1,N
C	 WRITE(6,*)INFO(I)
C	END DO

      F=EMIN
C
C
C	J=1
C      CALL BIOUT('     LinCoeff:.   AA1:        AA2:        B12:        
C	-RZ1:       RZ2:')
C      CALL BIOUT('FNC: =================================================
C     -=====================')
C 606  IF (J .LE. M) THEN
C          WRITE(6,1111) J,EVEC(J),AA1(J),AA2(J),B12(J),RZ1(J),RZ2(J)
C          WRITE(7,1111) J,EVEC(J),AA1(J),AA2(J),B12(J),RZ1(J),RZ2(J)
C          J=J+1
C
C      GOTO 606
C      END IF
C      CALL BIOUT('------------------------------------------------------
C     1---------------------')
CC
C1111  FORMAT(I3,6D12.4)
C
C
      RETURN
      END
