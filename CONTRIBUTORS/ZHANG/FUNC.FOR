      REAL*8 FUNCTION FUNC(p)
      PARAMETER (NMAX=50)
      PARAMETER (MMAX=50)
      COMMON /PASSIT/ M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,S,H,CO
      COMMON /SMALLONE/ RJUICE
      REAL*8 RJUICE
      INTEGER M,N,I,J
      REAL*8 p(NMAX),F
      REAL*8 
     1         AA1(MMAX),RX1(MMAX),RY1(MMAX),RZ1(MMAX),
     2         AA2(MMAX),RX2(MMAX),RY2(MMAX),RZ2(MMAX),
     3         B12(MMAX)
      REAL*8 S(MMAX,MMAX),H(MMAX,MMAX),CO(MMAX)
C          F=(AA(1)-2)**2 + (AA(2)-7)**2
C      AA(1)=ABS(AA(1))
C      AA(2)=ABS(AA(2))
C      AA(3)=ABS(AA(3))
C      AA(4)=ABS(AA(4))
C      AA(5)=ABS(AA(5))
C      AA(6)=ABS(AA(6))
      N=M*5
      J=1
 1070 IF (J .LE. N) THEN
c          p(J)=DABS(p(J))
c          p(J+1)=DABS(p(J+1))
          p(J+2)=DABS(p(J+2))
          J=J+5
          GOTO 1070
      ENDIF

      J=1
      DO I=1,M

          AA1(I)=(p(J)**2) + (p(J)*p(J+2))
          AA2(I)=(p(J+1)**2) + (p(J+2)**2) + (p(J)*p(J+2))
          B12(I)=-p(J)*p(J+2)
          RZ1(I)=p(J+3)
          RZ2(I)=p(J+4)
          J=J+5
      END DO

      CALL HAMILTON(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,EVEC,S,H,CO)
c      CALL BIOUT (' ')          
c      I=1
c      J=1
c      CALL BIOUT('        AA1               AA2               B12')
c      CALL BIOUT('        ----------------------------------------------
c     1--')
c 606  IF (J .LE. M) THEN
c 1111 FORMAT(I3,5D10.4)
c          WRITE(6,1111) I,P(I),P(I+1),P(I+2),P(I+3),P(I+4)
C          WRITE(7,1111) I,AA1(I),AA2(I),B12(I)
c          I=I+5
c          J=J+1
c      GOTO 606
c      END IF
      FUNC=RJUICE
c      WRITE(6,*) 'RJUICE = ',RJUICE
c      pause 'energy above...'
      RETURN
      END
