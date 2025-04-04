      SUBROUTINE FU(I)
      implicit real*8(a-h,o-z)
      PARAMETER (MMAX=50)
      COMMON /PASSIT/ M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,S,H,CO
C	 COMMON /PASSIT/ M
C      COMMON /SMALLONE/ RJUICE
	COMMON/T1/ N,NMAX,F,R(500),DR(500),CURVE,ALPHA,NOISE,EXACT,
     *           ALPH(500)
     1,RNOT(500), V(500,500), CURV(500),RESTART, 
     1 JSTART, FNOT, TIMLIM,      
     2 TOLF,TOLR       
      COMMON /VALUE/ EMIN
C      REAL*8 RJUICE
      INTEGER I,J,M
C      REAL*8 p(NMAX),F
      REAL*8 
     1         AA1(MMAX),RX1(MMAX),RY1(MMAX),RZ1(MMAX),
     2         AA2(MMAX),RX2(MMAX),RY2(MMAX),RZ2(MMAX),
     3         B12(MMAX)
      REAL*8 S(MMAX,MMAX),H(MMAX,MMAX),CO(MMAX)  
      I=0
C          F=(AA(1)-2)**2 + (AA(2)-7)**2
C      AA(1)=ABS(AA(1))
C      AA(2)=ABS(AA(2))
C      AA(3)=ABS(AA(3))
C      AA(4)=ABS(AA(4))
C      AA(5)=ABS(AA(5))
C      AA(6)=ABS(AA(6))
C      N=M*5
C      J=1
C 1070 IF (J .LE. N) THEN
C          r(J)=DABS(r(J))
C          r(J+1)=DABS(r(J+1))
C          r(J+2)=DABS(r(J+2))
C          J=J+5
C          GOTO 1070
C      ENDIF
C
c	 REWIND 7

      J=1
      DO I=1,M

C      write(6,*)I,J,M
         AA1(I)=(r(J)**2) + (r(J)*r(J+2))
         AA2(I)=(r(J+1)**2) + (r(J+2)**2) + (r(J)*r(J+2))
         B12(I)=-r(J)*r(J+2)
C
c      IF (AA1(I) .le. B12(I)) then
c          AA1(I)=AA1(I)+B12(I)
c      END IF
c      IF (AA2(I) .le. B12(I)) then
c          AA2(I)=AA2(I)+B12(I)
c      END IF



C      IF (r(J) .le. r(J+2)) then
C          r(J)=r(J)+r(J+2)
C      END IF
C      IF (r(J+1) .le. r(J+2)) then
C          r(J+1)=r(J+1)+r(J+2)
C      END IF
c          AA1(I)=DABS(r(J))
c          AA2(I)=DABS(r(J+1))
c          B12(I)=DABS(r(J+2))
          RZ1(I)=r(J+3)
          RZ2(I)=r(J+4)


C          WRITE(6,*)I,AA1(I),AA2(I),B12(I),RZ1(I),RZ2(I)
      
C          WRITE(7,*)I,AA1(I),AA2(I),B12(I),RZ1(I),RZ2(I)
          J=J+5
c      WRITE(7,747)AA1(I),AA2(I),B12(I),RZ1(I),RZ2(I),I
	
     
     
     
                     
      END DO


 747	FORMAT ('A1=',F18.12, '/ A2=',F18.12, '/ B12=',F18.12,'/ Z1=',
     *F18.12, '/ Z2=',F18.12,'   function='I3)



c     CALL BIOUT(' ')

      CALL HAMILTON(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,S,H,CO)
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
      F=EMIN
c      WRITE(6,*) 'RJUICE = ',RJUICE
c      pause 'energy above...'
      RETURN
      END
