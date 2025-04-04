C     Originallly Coded:  Z. Zhang
C           Modified By:  D. Gilmore ---> (24 May 95)
C     -----------------------------------------------------
C     L. ADAMOWICZ RESEARCH / UNIVERSITY OF ARIZONA / CHEM.
C     -----------------------------------------------------
C     --- MAIM PROGRAM CORE:  H2.F ---
C     -----------------------------------------------------   
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(MMAX=50)
C     MMAX = Maximum Number of Gaussians that may be used
C
      DIMENSION 
     1         AA1(MMAX),RX1(MMAX),RY1(MMAX),RZ1(MMAX),
     2         AA2(MMAX),RX2(MMAX),RY2(MMAX),RZ2(MMAX),
     3         B12(MMAX)
C      REAL*8 XI(MMAX,MMAX),P(MMAX),FRET
C     -----------------------------------------------------
C     AA1 = Exponent Alpha _sub(1)
C     RX1 = x-Component of Position for Gaussian of A1
C     RY1 = y-Component of Position for Gaussian of A1  
C     RZ1 = z-Component of Position for Gaussian of A1
C     AA2 = Exponent Alpha _sub(2)
C     RX2 = x-Component of Position for Gaussian of A2
C     RY2 = y-Component of Position for Gaussian of A2  
C     RZ2 = z-Component of Position for Gaussian of A2
C     B12 = Correlation Exponent Beta _sub(1,2)
C     -----------------------------------------------------
      DIMENSION S(MMAX,MMAX),H(MMAX,MMAX),CO(MMAX),G(5*MMAX)
C         ...Creates the S-Matrix and H-Matrix arrays
C     -----------------------------------------------------
C     S = M by M, S-Matrix / "Overlap Matrix"
C     H = M by M, H-Matrix / "Hamiltonian Matrix"
C     CO = Expansion Coefficients
C     G = ?
C     -----------------------------------------------------
      REAL*8 NOISE
      LOGICAL EXACT,RESTAR
      COMMON /T1/ ND,NMAX,FFF,RR(500),DR(500),CURVE,ALPHA,NOISE,EXACT,
     *           ALPH(500),RNOT(500),V(500,500),CURV(500),RESTAR,JSTART,
     *           FNOT,TIMLIM,TOLF,TOLR
      COMMON /VALUE/ EMIN
      COMMON /PASSIT/ M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,S,H,CO
      EXTERNAL FU
C
C
C     ---Open Output File: "OUTPUT.TXT"---
C     -----------------------------------------------------
      INCLUDE   'FLIB.FD'
      INTEGER*2 tmpday, tmpmonth, tmpyear
      INTEGER*2 tmphour, tmpminute, tmpsecond, tmphund
      OPEN (7, FILE = 'OUTPUT.TXT', ACCESS = 'SEQUENTIAL',
     1      STATUS = 'UNKNOWN')
C     -----------------------------------------------------
C     OUTPUT PROGRAM TITLE/CREDITS/& DATE STAMP:
C     -----------------------------------------------------
      WRITE(6,*)'Correlated Gaussian Geminal Output:  H_2'
      WRITE(6,*)'========================================'
      WRITE(6,*)'Adamowicz Research / [Univ. of AZ, CHEM]'
      WRITE(6,*)'Programming:      [Z. Zhang, D. Gilmore]' 
      
      WRITE (6,*) ' Current date and time...'
      CALL GETDAT(tmpyear, tmpmonth, tmpday)
      CALL GETTIM(tmphour, tmpminute, tmpsecond, tmphund)  
      WRITE (6, 9001) tmpmonth, tmpday, tmpyear
      WRITE (6, 9002) tmphour, tmpminute, tmpsecond, tmphund
      WRITE(6,*)'Output stream begins...'
      WRITE(7,*)'Correlated Gaussian Geminal Output:  H_2'
      WRITE(7,*)'========================================'
      WRITE(7,*)'Adamowicz Research / [Univ. of AZ, CHEM]'
      WRITE(7,*)'Programming:      [Z. Zhang, D. Gilmore]' 
      
      WRITE (7,*) ' PROGRAM BEGINS:  TIME STAMP...'
      CALL GETDAT(tmpyear, tmpmonth, tmpday)
      CALL GETTIM(tmphour, tmpminute, tmpsecond, tmphund)  
      WRITE (7, 9003) tmpmonth, tmpday, tmpyear
      WRITE (7, 9004) tmphour, tmpminute, tmpsecond, tmphund
      WRITE(7,*)'Output stream begins...'
C     -----------------------------------------------------
C
C     LOAD GEMINAL DATA:
C     -----------------------------------------------------
      CALL GUESS(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12)
C     -----------------------------------------------------
C

C      DO I=1,M
C      write(6,*)I,J,M
C          AA1(I)=(AA1(I)**2) + (AA1(I)*B12(I))
C           AA2(I)=(AA2(I)**2) + (B12(I)**2) + (AA1(I)*B12(I))
C         B12(I)=-AA1(I)*B12(I)
C
C      END DO




C
C     HAMILTONIAN SUBROUTINE:  Calcs. overlap & energy...
C     -----------------------------------------------------
      CALL HAMILTON(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,S,H,CO)
C     -----------------------------------------------------
      WRITE(6,*) '           UNOPTIMIZED GROUND STATE ENERGY = ',EMIN
      WRITE(7,*) '           UNOPTIMIZED GROUND STATE ENERGY = ',EMIN
      CALL BIOUT(' ')
C
c      STOP
C
      WRITE (7,*) ' UNOPTIMIZED ENERGY ROUTINE FINISHES:  TIME STAMP...'
      WRITE (6,*) ' UNOPTIMIZED ENERGY ROUTINE FINISHES:  TIME STAMP...'
C
      CALL GETDAT(tmpyear, tmpmonth, tmpday)
      CALL GETTIM(tmphour, tmpminute, tmpsecond, tmphund)  
C
      WRITE (6, 9003) tmpmonth, tmpday, tmpyear
      WRITE (6, 9004) tmphour, tmpminute, tmpsecond, tmphund
C
      WRITE (7, 9003) tmpmonth, tmpday, tmpyear
      WRITE (7, 9004) tmphour, tmpminute, tmpsecond, tmphund
C
C
C      CALL BIOUT (' ')
C      WRITE(6,*) '1-D VECTOR OPTIMIZATION BEGINS:'
C      WRITE(7,*) '1-D VECTOR OPTIMIZATION BEGINS:'
C      CALL BIOUT (' ')
C
C
C     NON-GRADIENT MINIMIZATION ROUTINE:
C     -----------------------------------------------------
C      CALL VECTORMIN
C      CALL BIOUT(' ')
C      WRITE (7,*) ' 1-D VECTOR MINIMIZATION FINISHES:  TIME STAMP...'
C      WRITE (6,*) ' 1-D VECTOR MINIMIZATION FINISHES:  TIME STAMP...'
C
C      CALL GETDAT(tmpyear, tmpmonth, tmpday)
C      CALL GETTIM(tmphour, tmpminute, tmpsecond, tmphund)  
C
C      WRITE (6, 9003) tmpmonth, tmpday, tmpyear
C      WRITE (6, 9004) tmphour, tmpminute, tmpsecond, tmphund
C
C      WRITE (7, 9003) tmpmonth, tmpday, tmpyear
C      WRITE (7, 9004) tmphour, tmpminute, tmpsecond, tmphund
C      CALL PLODDING(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,S,H)
C      STOP
C     -----------------------------------------------------
C
C      CALL TAB(xi,M,M,MMAX,MMAX)
C      CALL BIOUT(' ')
C      
C      FRETTEMP=0.0d+00
      CALL BIOUT (' ')          
C      j=1
C      DO i=1,M
C           p(j)=AA1(i)
C           p(j+1)=AA2(i)
C           p(j+2)=B12(i)
C           p(j+3)=RZ1(i)
C           p(j+4)=RZ2(i)
C           j=j+5
C      END DO

      j=1
      DO i=1,M
           rr(j)=AA1(i)
           rr(j+1)=AA2(i)
           rr(j+2)=B12(i)
           rr(j+3)=RZ1(i)
           rr(j+4)=RZ2(i)
           j=j+5
      END DO


      I=1
      J=1
      CALL BIOUT('     AA1         AA2         B12         RZ1         R
     1Z2')
      CALL BIOUT('     -------------------------------------------------
     1---------')

 606  IF (J .LE. M) THEN
 1111 FORMAT(I3,5D12.4)
          WRITE(6,1111) J,rr(I),rr(I+1),rr(I+2),rr(I+3),rr(I+4)
          I=I+5
          J=J+1
      GOTO 606
      END IF
      CALL BIOUT(' ')
C 1
C      DO i=1,M*5
C          DO j=1,M*5
C              xi(i,j)=0.0D+00
C              IF (i .eq. j) THEN
C                  xi(i,j)=1.00D+00
C              ENDIF
C          END DO
C      END DO   
C      DO i=1,M*5
C      WRITE(6,*) p(i)
C      END DO
C      CALL powell(p,xi,M*5,M*5,1.0D-12,100,fret)
c      CALL BIOUT(' ')
C
C 515     WRITE(6,*) '    POWELL OPTIMIZED GROUND STATE ENERGY = ',EMIN
C      WRITE(7,*) '    POWELL OPTIMIZED GROUND STATE ENERGY = ',EMIN
C

      TOLF = 1.0D-06
      TOLR = 1.0D-06
      NOISE = 1.0D-08
      ND = M*5
      NMAX = 180
      KKKK=-1
      FFF = EMIN

      CALL BIOUT(' ')
      WRITE (7,*) ' TRUDGE Optimization Begins:...'
      WRITE (6,*) ' TRUDGE Optimization Begins:...'
      CALL BIOUT(' ')

      CALL TRUDGE(KKKK,FU)

       CALL BIOUT (' ')          
      I=1
      J=1
      CALL BIOUT('     AA1         AA2         B12         RZ1         R
     1Z2')
      CALL BIOUT('     -------------------------------------------------
     1---------')
 607  IF (J .LE. M) THEN

          WRITE(6,1111) J,rr(I),rr(I+1),rr(I+2),rr(I+3),rr(I+4)
          I=I+5
          J=J+1
      GOTO 607
      END IF


C
C     GRADIENT MINIMIZATION ROUTINE:  (vestigial from Zhang)
C     -----------------------------------------------------
C      CALL GRAD(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,S,H,CO,G)
c      STOP
C     -----------------------------------------------------
C
C
C     FORMAT(s):
C     -----------------------------------------------------
1007  FORMAT(/) 
9001  FORMAT ('     ',I2, '/', I2.2, '/', I4.4)  
9002  FORMAT ('     ',I2, ':', I2.2, ':', I2.2, ':', I2.2, ' ', A, 'm')
9003  FORMAT ('     ',I2, '/', I2.2, '/', I4.4)  
9004  FORMAT ('     ',I2, ':', I2.2, ':', I2.2, ':', I2.2, ' ', A, 'm')
C     -----------------------------------------------------
C
C     TERMINATE PROGRAM EXECUTION!
      CALL BIOUT(' ')
      WRITE (7,*) ' PROGRAM FINISHES:  TIME STAMP...'
      WRITE (6,*) ' PROGRAM FINISHES:  TIME STAMP...'

      CALL GETDAT(tmpyear, tmpmonth, tmpday)
      CALL GETTIM(tmphour, tmpminute, tmpsecond, tmphund)  

      WRITE (6, 9003) tmpmonth, tmpday, tmpyear
      WRITE (6, 9004) tmphour, tmpminute, tmpsecond, tmphund

      WRITE (7, 9003) tmpmonth, tmpday, tmpyear
      WRITE (7, 9004) tmphour, tmpminute, tmpsecond, tmphund
C      WRITE (6,*) EMIN
      END
