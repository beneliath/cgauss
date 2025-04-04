      SUBROUTINE YACOBI (A,U,N,NNP1D2,BIG,JB)
C     L(I)=(I*(I-1))/2
CSW      1
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NNP1D2),U(N*N),BIG(N),JB(N)
      DATA ZERO/0.0D+00/
      DATA EPS/1.0D-14/
      DATA ONE/1.0D+00/
CSW      2
C      ABS(X)=DABS(X)
C      SQRT(X)=DSQRT(X)
      L(I)=(I*(I-1))/2.0D+00
C ***
C      WRITE(6,*)'YACOBI: N,NNP1D2', N,NNP1D2
      U(1)=ONE
      IF (N.EQ.1) GO TO 999
      J=0
      K=0
      N2=N*N
      DO 10 I=1,N2
   10 U(I)=ZERO
      DO 30 I=1,N
      J=J+I
      JJ=K+I
      K=K+N
      U(JJ)=ONE
C     SECTION TO LOCATE BIGGEST OFF-DIAGONAL ELEMENT IN ROW I
   15 IM1=I-1
      IF (IM1.EQ.0) GO TO 30
   16 BIG(I)=ZERO
      DO 20 JJ=1,IM1
      IJ=L(I)+JJ
      B=ABS(A(IJ))
      IF (B.LE.BIG(I)) GO TO 20
   18 BIG(I)=B
      JB(I)=JJ
   20 CONTINUE
   30 CONTINUE
C     LOCATE ROTATION PIVOT
   31 CONTINUE
      IBIG=2
      DO 35 I=2,N
      IF (BIG(I).GT.BIG(IBIG)) IBIG=I
   35 CONTINUE
      BIGGST=BIG(IBIG)
      JBIG=JB(IBIG)
      IF (BIGGST.LE.EPS) GO TO 999
C     BEGIN ROTATION
   40 II=L(IBIG)+IBIG
      IJ=L(IBIG)+JBIG
      JJ=L(JBIG)+JBIG
      AB=A(II)-A(JJ)
      AA=A(IJ)+A(IJ)
      D=SQRT(AB*AB+AA*AA)
      IF (AB.LT.ZERO) D=-D
      T=AA/(AB+D)
      TSQ=T*T
      CSQ=ONE/(ONE+TSQ)
      C=SQRT(CSQ)
      S=C*T
      AB=AA*T
      D=(A(II)+AB+TSQ*A(JJ))*CSQ
      A(JJ)=(A(JJ)-AB+TSQ*A(II))*CSQ
      A(II)=D
      A(IJ)=ZERO
      III=4
      KIX=(IBIG-1)*N
      KJX=(JBIG-1)*N
      DO 60 I=1,N
      KI=KIX+I
      KJ=KJX+I
      D=C*U(KI)+S*U(KJ)
      U(KJ)=C*U(KJ)-S*U(KI)
      U(KI)=D
      IF (I-JBIG) 45,54,46
   45 KJ=L(JBIG)+I
      KI=L(IBIG)+I
      GO TO 50
   46 KJ=L(I)+JBIG
      IF (I-IBIG)47,55,48
   47 KI=L(IBIG)+I
      III=2
      IF (JB(I).EQ.JBIG) III=3
      GO TO 50
   48 KI=L(I)+IBIG
      III=1
      IF (JB(I).EQ.IBIG.OR.JB(I).EQ.JBIG) III=3
   50 D=C*A(KI)+S*A(KJ)
      A(KJ)=C*A(KJ)-S*A(KI)
      A(KI)=D
      GO TO (51,52,55,60),III
   51 B=ABS(A(KI))
      IF (B.LE.BIG(I)) GO TO 52
      BIG(I)=B
      JB(I)=IBIG
   52 B=ABS(A(KJ))
      IF (B.LE.BIG(I)) GO TO 60
   53 BIG(I)=B
      JB(I)=JBIG
      GO TO 60
   54 IF (I.EQ.1) GO TO 60
   55 IM1=I-1
      BIG(I)=ZERO
      DO 58 J=1,IM1
      IJ=L(I)+J
      B=ABS(A(IJ))
      IF (B.LE.BIG(I)) GO TO 58
      BIG(I)=B
      JB(I)=J
   58 CONTINUE
   60 CONTINUE
C     END OF ROTATION
      GO TO 31
  999 CONTINUE
C LOCK EIGENVECTORS SO THAT FIRST COEFFICIENT IS ALWAYS
C POSITIVE.
      IJ=1
      DO 501 I=1,N
      IF (U(IJ).GE.ZERO) GO TO 503
      DO 502 J=1,N
      U(IJ+J-1)=-U(IJ+J-1)
  502 CONTINUE
  503 IJ=IJ+N
  501 CONTINUE
      RETURN
      END
