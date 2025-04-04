C     Originallly Coded:  Z. Zhang
C           Modified By:  D. Gilmore ---> (24 May 95)
C     -----------------------------------------------------
C     L. ADAMOWICZ RESEARCH / UNIVERSITY OF ARIZONA / CHEM.
C     -----------------------------------------------------
      SUBROUTINE GUESS(M,A1,RX1,RY1,RZ1,A2,RX2,RY2,RZ2,B12)
C     -----------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
	PARAMETER(MMAX=50)
      PARAMETER(MF=5)
      REAL*8 A1(MMAX),RX1(MMAX),RY1(MMAX),RZ1(MMAX),
     1          A2(MMAX),RX2(MMAX),RY2(MMAX),RZ2(MMAX),
     2          B12(MMAX)
C     -----------------------------------------------------
C     --- Data Follows ---
C     -----------------------------------------------------
C
C     --- BD = Bond Distance for H2 in Atomic Units,
C             (Bohr Radii) ---
C
C      INCLUDE 'FLIB.FI'
      INCLUDE 'FLIB.FD'

      REAL*4 TEMP
      REAL*8 TEMP1,TEMP2
      INTEGER I,M

      PARAMETER(BD=1.4D+00)
      M=MF
C
C     -----------------------------------------------------
C     M = Number of Gaussian Geminals
C     A1 = Exponent Alpha _sub(1)
C     RX1 = x-Component of Position for Gaussian of A1
C     RY1 = y-Component of Position for Gaussian of A1
C     RZ1 = z-Component of Position for Gaussian of A1
C     A2 = Exponent Alpha _sub(2)
C     RX2 = x-Component of Position for Gaussian of A2
C     RY2 = y-Component of Position for Gaussian of A2
C     RZ2 = z-Component of Position for Gaussian of A2
C     B12 = Correlation Exponent Beta _sub(1,2)
C     -----------------------------------------------------
C     General Form of WaveFunction:
C
C     psi=Exp(A1*|rA1|**2)Exp(A2*|rB2|**2)Exp(B12*|r12|**2)
C
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
C     =============================================
C     RANDOMLY GENERATED PARAMETERS:
C     =============================================
C      CALL BIOUT('GENERATING RANDOM Gaussian Geminal Data...')
C      CALL BIOUT(' ')
C      I=1
C      CALL SEED(RND$TIMESEED)
C      CALL BIOUT('     AA1               AA2               B12'
C     1)
C      CALL BIOUT('     -------------------------------------------------
C     1-----')
C 207  IF (I .LE. M) THEN
C          CALL RANDOM(TEMP)
C           TEMP1=DBLE(TEMP)*1.0D+00
C          CALL RANDOM(TEMP)
C           TEMP2=DBLE(TEMP)*1.0D+00
C          A1(I)=TEMP1
CC          CALL RANDOM(TEMP)
CC           TEMP1=DBLE(TEMP)*1.0D-02
CC          A2(I)=I+(TEMP1/I)
C          A2(I)=TEMP2
C            RZ1(I)=0.7D+00
C            RZ2(I)=-0.7D+00
C             RX1(I)=0.0D+00
C             RX2(I)=0.0D+00
C             RY1(I)=0.0D+00
C             RY2(I)=0.0D+00
C          CALL RANDOM(TEMP)
C          TEMP1=DBLE(TEMP)*1.0D-3
C          B12(I)=TEMP1
C          WRITE(6,709)A1(I),A2(I),B12(I)
C          WRITE(7,709)A1(I),A2(I),B12(I)
C          I=I+1
C      GOTO 207
C      END IF
C 709  FORMAT(E18.8,E18.8,E18.8)
C      CALL BIOUT(' ')
C     =============================================
C     HARD-CODED PARAMETERS:
C     =============================================
      A1(1)=1.488413547506
      RX1(1)=0.0D+00
      RY1(1)=0.0D+00
      RZ1(1)=0.7D+00
      A2(1)=1.098578067090
      RX2(1)=0.0D+00
      RY2(1)=0.0D+00
      RZ2(1)=-0.7D+00
      B12(1)=0.0D+00

C
      A1(2)=2.019446053312
      RX1(2)=0.0D+00 
      RY1(2)=0.0D+00
      RZ1(2)=0.7D+00
      A2(2)=2.287955039062
      RX2(2)=0.0D+00 
      RY2(2)=0.0D+00
      RZ2(2)=-0.7D+00
      B12(2)= 0.0D+00

C
      A1(3)=3.477626118318
      RX1(3)=0.0D+00 
      RY1(3)=0.0D+00
      RZ1(3)=0.7D+00
      A2(3)=3.261059340400
      RX2(3)=0.0D+00
      RY2(3)=0.0D+00
      RZ2(3)=-0.7D+00
      B12(3)=0.0D+00

C
      A1(4)=4.406191682323
      RX1(4)=0.0D+00
      RY1(4)=0.0D+00
      RZ1(4)=0.7D+00
      A2(4)=4.153396206550
      RX2(4)=0.0D+00 
      RY2(4)=0.0D+00
      RZ2(4)=-0.7D+00
      B12(4)=0.0D+00

C
      A1(5)=5.765668067546
      RX1(5)=0.0D+00
      RY1(5)=0.0D+00
      RZ1(5)=0.7D+00
      A2(5)=5.394665357549
      RX2(5)=0.0D+00
      RY2(5)=0.0D+00
      RZ2(5)=-0.7D+00
      B12(5)=0.0D+00

C
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
      RETURN
      END
