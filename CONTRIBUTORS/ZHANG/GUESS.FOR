C     L. ADAMOWICZ RESEARCH / UNIVERSITY OF ARIZONA / CHEM.
C     -----------------------------------------------------
      SUBROUTINE GUESS(METHOD,M,A1,RX1,RY1,RZ1,A2,RX2,RY2,RZ2,B12,EVEC)
C     -----------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
	INTEGER METHOD,M,ERROR,I,MMAX,MF
      REAL*8, DIMENSION(M) ::	A1,RX1,RY1,RZ1,
	-						A2,RX2,RY2,RZ2,
     -						B12,EVEC
      REAL*4 TEMP
      REAL*8 TEMP1,TEMP2,BD
C     -----------------------------------------------------
      PARAMETER (BD=1.4D+00)
C     -----------------------------------------------------
C     --- Data Follows ---
C     -----------------------------------------------------
      INCLUDE 'FLIB.FD'
C     -----------------------------------------------------
C     General Form of WaveFunction:
C
C     psi=Exp(A1*|rA1|**2)Exp(A2*|rB2|**2)Exp(B12*|r12|**2)
C     -----------------------------------------------------
C     =============================================
C     READ FROM 'h2.in':
C     =============================================
C      CALL BIOUT('READING (h2.in) Gaussian Geminal Data...')
C      CALL BIOUT(' ')
C      CALL BIOUT('     AA1               AA2               B12'
C     1)
C      CALL BIOUT('     -------------------------------------------------
C     1-----')
C      CALL BIOUT(' ')
C      OPEN(5, FILE='h2.in')
C      I=1
C 5207  IF (I .LE. M) THEN
C            READ(5,*) TEMP1
C            RX1(I)=TEMP1
C            READ(5,*) TEMP1
C            RY1(I)=TEMP1
C            READ(5,*) TEMP1
C            RZ1(I)=TEMP1
CC
C            READ(5,*) TEMP1
C            RX2(I)=TEMP1
C            READ(5,*) TEMP1
C            RY2(I)=TEMP1
C            READ(5,*) TEMP1
C            RZ2(I)=TEMP1
CC
C            READ(5,*) TEMP1
C            READ(5,*) TEMP2
C            A1(I)=TEMP1
C            A2(I)=TEMP1
C            B12(I)=TEMP2
C          WRITE(7,709) A1(I),A2(I),B12(I)
C          WRITE(6,709) A1(I),A2(I),B12(I)
C          I=I+1
C      GOTO 5207
C      END IF
C      CALL BIOUT(' ')
C
C
C
C     =============================================
C     RANDOMLY GENERATED PARAMETERS:
C     =============================================
	IF (METHOD.EQ.1) THEN
     		CALL BIOUT('GENERATING RANDOM GAUSSIAN GEMINAL DATA...')
		CALL BIOUT(' ')
		I=1
		CALL SEED(RND$TIMESEED)
	CALL BIOUT('      T1                T2               T12
	-    LinCoeff.')
      CALL BIOUT('     -------------------------------------------------
     --------------------------')
C
	POINT=1.00D-12
      DO I=1,M
	    CALL RANDOM(TEMP)
			TEMP1=DBLE(TEMP)*1.0D-01
          CALL RANDOM(TEMP)
			TEMP2=DBLE(TEMP)*1.0D-01
          CALL RANDOM(TEMP)
			TEMP3=DBLE(TEMP)*1.0D-12
c	A1(I)=TEMP1*10**((8.089714D+00*EXP(-POINT**2 /2/(2.3/3)**2)-
c	-3.089714)/3)
c	A2(I)=TEMP2*10**((8.089714D+00*EXP(-POINT**2 /2/(2.3/3)**2)-
c	-3.089714)/3)
	A1=TEMP1
	A2=TEMP2
C				B12(I)=1.00D-8
				B12(I)=TEMP3
C          CALL RANDOM(TEMP)
C			TEMP1=DBLE(TEMP)*1.0D-02
C				A2(I)=I+(TEMP1/I)
		RZ1(I)=0.700D+00
		RZ2(I)=-0.700D+00
		POINT=POINT+1.00D+00/M
				RX1(I)=0.0D+00
				RX2(I)=0.0D+00
				RY1(I)=0.0D+00
				RY2(I)=0.0D+00
		CALL RANDOM(TEMP)
			EVEC(I)=DBLE(TEMP)*1.0D-02
				WRITE(6,709)A1(I),A2(I),B12(I),EVEC(I),RZ1(I),RZ2(I)
				WRITE(7,709)A1(I),A2(I),B12(I),EVEC(I),RZ1(I),RZ2(I)
      END DO
	CALL BIOUT(' ')
	END IF
C     =============================================
C     HARD-CODED PARAMETERS:
C     =============================================
C      A1(1)=.160235037618
C      RX1(1)=0.0D+00
C      RY1(1)=0.0D+00
C      RZ1(1)=1.011036314027
C      A2(1)=.350880170955
C      RX2(1)=0.0D+00
C      RY2(1)=0.0D+00
C      RZ2(1)=-.063892060873
C      B12(1)=-.018048129548
C
CC
C      A1(2)=.296378877724
C      RX1(2)=0.0D+00 
C      RY1(2)=0.0D+00
C      RZ1(2)=.051541492960
C      A2(2)=1.151304709732
C      RX2(2)=0.0D+00 
C      RY2(2)=0.0D+00
C      RZ2(2)=-.598756422434
C      B12(2)= -.056849064138
C
CC
C      A1(3)=.929110466028
C      RX1(3)=0.0D+00 
C      RY1(3)=0.0D+00
C      RZ1(3)=.661673616984
C      A2(3)=3.457894696749
C      RX2(3)=0.0D+00
C      RY2(3)=0.0D+00
C      RZ2(3)=-.698358442171
C      B12(3)=-.032480513174
C
C	EVEC(1)=0.1234D+00
C	EVEC(2)=0.2345D+00
C	EVEC(3)=0.3456D+00
CC
C      A1(4)=.089838780145
C      RX1(4)=0.0D+00
C      RY1(4)=0.0D+00
C      RZ1(4)=.904358319677
C      A2(4)=8.592120547377
C      RX2(4)=0.0D+00 
C      RY2(4)=0.0D+00
C      RZ2(4)=-.687430054083
C      B12(4)=.088393827404
CC
C      A1(5)=-.160198996335
C      RX1(5)=0.0D+00
C      RY1(5)=0.0D+00
C      RZ1(5)=-.105970239485
C      A2(5)=131.528440491213
C      RX2(5)=0.0D+00
C      RY2(5)=0.0D+00
C      RZ2(5)=-.684271209383
C      B12(5)=.611556470151
CC
C      A1(6)=.21187121D+00
C      RX1(6)=0.0D+00
C      RY1(6)=0.0D+00
C      RZ1(6)=0.7D+00
C      A2(6)=.21187121D+00
C      RX2(6)=0.0D+00
C      RY2(6)=0.0D+00
C      RZ2(6)=-0.7D+00
C      B12(6)=.14772309D+00
CC
C      A1(7)=0.01D+00
C      RX1(7)=0.0 D+00
C      RY1(7)=0.0D+00
C      RZ1(7)=0.7D+00
C      A2(7)=0.01D+00
C      RX2(7)=0.0D+00
C      RY2(7)=0.0D+00
C      RZ2(7)=-0.7D+00
C      B12(7)=0.0008D+00
CC
C      A1(8)=27.0D+00
C      RX1(8)=0.0 D+00
C      RY1(8)=0.0D+00
C      RZ1(8)=0.7D+00
C      A2(8)=27.0D+00
C      RX2(8)=0.0D+00
C      RY2(8)=0.0D+00
C      RZ2(8)=-0.7D+00
C      B12(8)=0.2777D+00
CC
C      A1(9)=3.0D+00
C      RX1(9)=0.0 D+00
C      RY1(9)=0.0D+00
C      RZ1(9)=0.7D+00
C      A2(9)=3.0D+00
C      RX2(9)=0.0D+00
C      RY2(9)=0.0D+00
C      RZ2(9)=-0.7D+00
C      B12(9)=0.0003D+00
CC
C      A1(10)=4000.0D+00
C      RX1(10)=0.0 D+00
C      RY1(10)=0.0D+00
C      RZ1(10)=0.7D+00
C      A2(10)=4000.0D+00
C      RX2(10)=0.0D+00
C      RY2(10)=0.0D+00
C      RZ2(10)=-0.7D+00
C      B12(10)=12.47D+00
C     =============================================
C
C     -----------------------------------------------------
C     FORMAT(s):
C     -----------------------------------------------------
 709  FORMAT(6E18.8)
C
      RETURN
      END