C     Originallly Coded:  Z. Zhang
C           Modified By:  D. Gilmore ---> (24 May 95)
C     -----------------------------------------------------
C     L. ADAMOWICZ RESEARCH / UNIVERSITY OF ARIZONA / CHEM.
C     -----------------------------------------------------
      SUBROUTINE OVERLP(AX,AY,AZ,A1,BX,BY,BZ,A2,                      
     1                  CX,CY,CZ,A3,DX,DY,DZ,A4,
     2                  A5,A6,OV,KLM) 
C     ----------------------------------------------------- 
C         ...OVERLAP INTEGRAL = < psi(L) | psi {k) >
C     -----------------------------------------------------
C     Reference For This Calculation Algorythm:
C         "Gaussian Correlation Functions:  Two-Electron
C          Systems", Lester & Krauss, Jorn. Chem. Phys.,
C          Vol.41,No.5,pp.1407-1413,(sp. p.1411)
C     -----------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C     -----------------------------------------------------
C     A1 = Exponent Alpha _sub1(L)
C     AX = x-Component of Position for Gaussian of A1
C     AY = y-Component of Position for Gaussian of A1  
C     AZ = z-Component of Position for Gaussian of A1
C     A2 = Exponent Alpha _sub1(K)
C     BX = x-Component of Position for Gaussian of A2
C     BY = y-Component of Position for Gaussian of A2  
C     BZ = z-Component of Position for Gaussian of A2
C     A3 = Exponent Alpha _sub2(L)
C     CX = x-Component of Position for Gaussian of A3
C     CY = y-Component of Position for Gaussian of A3  
C     CZ = z-Component of Position for Gaussian of A3
C     A4 = Exponent Alpha _sub2(K)
C     DX = x-Component of Position for Gaussian of A4
C     DY = y-Component of Position for Gaussian of A4  
C     DZ = z-Component of Position for Gaussian of A4
C     A5 = Correlation Exponent Beta _sub1,2(L)
C     A6 = Correlation Exponent Beta _sub1,2(K)
C     OV = ? 
C     KLM = ?
C     -----------------------------------------------------
C     --- FOUR CENTER OVERLAP INTEGRAL ---
C     -----------------------------------------------------
      PARAMETER(ZERO=0.0D+00)
C         ..."ZERO" is a fixed parameter of Double Precision
C             and has value of 0.0E00 
C
      PI=3.1415926535898D+00
C         ..."PI" is a Double Precision variable of
C             appropriate value.
C
C     From HAMILTON.F:   DO L=1,M
C                             Do K=1,M ...
C


C      REWIND 7
c      WRITE(7,747)AZ,A1,BZ,A2,A5
c      WRITE(7,757)CZ,A3,DZ,A4,A6
c      WRITE(7,*)
c 747	FORMAT ('1,1 Z=',F18.12, '/ 1,1 A=',F18.12, '/ 1,2 Z=',F18.12,'/ 1
c     *,2 A=',F18.12, '/ 1B12=',F18.12)
c 757	FORMAT ('2,1 Z=',F18.12, '/ 2,1 A=',F18.12, '/ 2,2 Z=',F18.12,'/ 2
c     *,2 A=',F18.12, '/ 2B12=',F18.12)


      ALK1=A1+A2
C         ...ALK1 = Alpha _sub1(L) + Alpha _sub1(K)
C                                                       
      ALK2=A3+A4                                                       
C         ...ALK2 = Alpha _sub2(L) + Alpha _sub2(K)
C
      BLK12=A5+A6  
C         ...BLK12 = BETA12(L) + BETA12(K)
C
      DETAB=ALK1*ALK2+(ALK1+ALK2)*BLK12
C     Modification (dgilmore) follows: 
C      DETAB=(ALK1*ALK2)+((ALK1+ALK2)*(2.0D0*BLK12))
C         ...DETAB = Determinant of AB
C
      IF(DETAB.LE.ZERO) GOTO 300
C         ...If DETAB is LESS THAN OR EQUAL TO ZERO,
C            then GOTO 300.
C
C     -----------------------------------------------------
C     --- ONE CENTER OVERLAP INTEGRAL ---
C     -----------------------------------------------------
C      OV=(PI*PI*PI)/(DETAB**(1.5D0))
CCC      OV=PI**3/(DETAB*DSQRT(DETAB))
	OV=1.0D+00/(DETAB*DSQRT(DETAB))
C         ...OV = Magnitude of Overlap (???)
C
C     -----------------------------------------------------
C     --- CENTERS ---
C     -----------------------------------------------------
      CALL KFACT(AX,AY,AZ,A1,BX,BY,BZ,A2,                      
     1           CX,CY,CZ,A3,DX,DY,DZ,A4,
     2           A5,A6,XK) 
C         ...KFACT Obtains the K-Factor.  The Product of 
C             any two gaussians is proportional to a third
C             gaussian centered between them.  K-Factor is
C             a multiplicative factor on this product
C             gaussian.  (ref. Szabo & Ostlund, Appendix A)
C                         "Modern Quantum Chemistry", 1982.
C
C     -----------------------------------------------------
C     --- FOUR-CENTER OVERLAP INTEGRAL ---
C     -----------------------------------------------------
      OV=OV*XK
C         ...OV = Magnitude of Overlap      
C            Apparently the total overlap is simply the K-
C            Factor (XK) multiplied into the simple core
C            overlap expression above.   
C
C     -----------------------------------------------------
      GOTO 1500
C     -----------------------------------------------------
300   KLM=5
C         ...KLM = ?
C            If the Determinant of AB is less than or equal
C            to zero, then set KLM equal to 5.
C            Possibly some kind of flag ???
C
C     -----------------------------------------------------
1500  RETURN                                                            
      END      
