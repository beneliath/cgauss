      SUBROUTINE VECTORMIN
      PARAMETER (NMAX=50)
      PARAMETER (MMAX=50)
      REAL*8 
     1         AA1(MMAX),RX1(MMAX),RY1(MMAX),RZ1(MMAX),
     2         AA2(MMAX),RX2(MMAX),RY2(MMAX),RZ2(MMAX),
     3         B12(MMAX)
      REAL*8 S(MMAX,MMAX),H(MMAX,MMAX),CO(MMAX),TEMP1,TEMP2
      INTEGER N,M,I,J
      REAL*8 POINT(NMAX),DIRECTION(NMAX),COMPARE,TOLERANCE,FRET
      COMMON /PASSIT/ M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,S,H,CO
      COMMON /VALUE/ EMIN
      N=M*2
      TOLERANCE=1.0D-06
      WRITE(6,777) TOLERANCE
      WRITE(7,777) TOLERANCE
      WRITE(6,*) '    Number of Geminals = ',M
      WRITE(7,*) '    Number of Geminals = ',M
      CALL BIOUT(' ')
 777  FORMAT('     Convergence Target = ',E18.4,/)
C
C
      I=1
      J=1
 707  IF (I .LE. M) THEN
          POINT(J)=AA1(I)
          POINT(J+1)=B12(I)
          I=I+1
          J=J+2
          GOTO 707
      ENDIF
C      POINT(1)=AA1(1)
C      POINT(2)=B12(1)
C      POINT(3)=AA1(2)
C      POINT(4)=B12(2)
C      POINT(5)=AA1(3)
C      POINT(6)=B12(3)
C
C
 1001 I=1
c      J=1
 1002 CALL BLANKDIRECTION (N,DIRECTION,NMAX)
C      IF (J .LE. N) THEN
C          DIRECTION(J)=0.0D+00
C          J=J+1
C          GOTO 1002
C      ENDIF
C      WRITE(6,*) I,J
      IF (I .LE. N) THEN
              DIRECTION(I)=1.0D+00
              I=I+1
              CALL LINMIN(POINT,DIRECTION,N,FRET)
          GOTO 1002
      ENDIF
C      
      COMPARE=FRET
C
      I=1
C      J=1
 1003 CALL BLANKDIRECTION (N,DIRECTION,NMAX)
C      IF (J .LE. N) THEN
C          DIRECTION(J)=0.0D+00
C          J=J+1
C          GOTO 1003
C      ENDIF
C      WRITE(6,*) I,J
          IF (I .LE. N) THEN
              DIRECTION(I)=1.0D+00
              I=I+1
              CALL LINMIN(POINT,DIRECTION,N,FRET)
          GOTO 1003
      ENDIF
C
C
      TEMP1=FRET
      TEMP2=COMPARE
      WRITE(6,1856) FRET,ABS(TEMP1-TEMP2)
      WRITE(7,1856) FRET,ABS(TEMP1-TEMP2)
1856  FORMAT (' E_o = ',F18.12,'  /  Iteration Difference = ',E18.4)

      IF (ABS(TEMP1-TEMP2) .GT. TOLERANCE) THEN
          GOTO 1001
      END IF
      CALL BIOUT (' ')
           CALL BIOUT ('   OPTIMUM GAUSSIAN EXPONENTS:')
      CALL BIOUT (' ')          
      I=1
      J=1
      CALL BIOUT('        AA1               AA2               B12')
      CALL BIOUT('        ----------------------------------------------
     1--')
 606  IF (I .LE. M) THEN
 1111 FORMAT(I3,D18.8,D18.8,D18.8)
          WRITE(6,1111) I,AA1(I),AA2(I),B12(I)
          WRITE(7,1111) I,AA1(I),AA2(I),B12(I)
C          WRITE(7,1111) I,POINT(J),POINT(J),POINT(J+1)
          I=I+1
          J=J+2
      GOTO 606
      END IF
      CALL BIOUT (' ')
      WRITE(6,*) '   FUNCTION VALUE AT MINIMUM= ',FRET
      WRITE(7,*) '   FUNCTION VALUE AT MINIMUM= ',FRET
      CALL BIOUT (' ')
      WRITE(6,*) '   Convergence Index = ',ABS(TEMP1-TEMP2)
      WRITE(7,*) '   Convergence Index = ',ABS(TEMP1-TEMP2)
      RETURN
      END
