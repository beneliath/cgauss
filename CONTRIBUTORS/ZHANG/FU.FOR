      SUBROUTINE FU(I)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (MMAX=50)
      COMMON /PASSIT/ M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,EVEC,S,H,CO
	COMMON/T1/ N,NMAX,F,R(500),DR(500),CURVE,ALPHA,NOISE,EXACT,
     1           ALPH(500),RNOT(500), V(500,500), CURV(500),
     2		   RESTART,JSTART, FNOT, TIMLIM,TOLF,TOLR       
      COMMON /VALUE/ EMIN
	COMMON /LC/RC(50)
      INTEGER I,J,M
      REAL*8 
     1         AA1(MMAX),RX1(MMAX),RY1(MMAX),RZ1(MMAX),
     2         AA2(MMAX),RX2(MMAX),RY2(MMAX),RZ2(MMAX),
     3         B12(MMAX),EVEC(MMAX),S(MMAX,MMAX),H(MMAX,MMAX),CO(MMAX)  
      I=0
C	 REWIND 7
C      J=1
C      DO I=1,M
C		TEMP1=(r(J)**2)+(r(J)*r(J+2))
C		TEMP2=(r(J+1)**2)+(r(J+2)**2) + (r(J)*r(J+2))
C		TEMP3=-r(J)*r(J+2)
C	AA1(I)=TEMP1
C	AA2(I)=TEMP2
C	B12(I)=TEMP3
C	    RZ1(I)=r(J+3)
C          RZ2(I)=r(J+4)
C          J=J+5
CC		 WRITE(7,747)AA1(I),AA2(I),B12(I),RZ1(I),RZ2(I),I
C      END DO
 747	FORMAT ('A1=',F18.12, '/ A2=',F18.12, '/ B12=',F18.12,'/ Z1=',
     *F18.12, '/ Z2=',F18.12,'   function='I3)
C
      CALL HAMILTON(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,EVEC,S,H,CO)
C
C      CALL BIOUT (' ')          
C      I=1
C      J=1
C      CALL BIOUT('        AA1               AA2               B12')
C      CALL BIOUT('        ----------------------------------------------
C     1--')
C 606  IF (J .LE. M) THEN
C          WRITE(6,1111) I,P(I),P(I+1),P(I+2),P(I+3),P(I+4)
C          WRITE(7,1111) I,AA1(I),AA2(I),B12(I)
C          I=I+5
C          J=J+1
C      GOTO 606
C      END IF
C 1111 FORMAT(I3,5D10.4)
C
      F=EMIN
C
      RETURN
      END
