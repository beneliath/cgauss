C     Originallly Coded:  Z. Zhang
C           Modified By:  D. Gilmore ---> (24 May 95)
C     -----------------------------------------------------
C     L. ADAMOWICZ RESEARCH / UNIVERSITY OF ARIZONA / CHEM.
C     -----------------------------------------------------
      SUBROUTINE GEOM(X1,Y1,Z1,X2,Y2,Z2)
C     -----------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C     -----------------------------------------------------
C     --- GEOMETRY (INPUT DATA) --- 
C     -----------------------------------------------------
      PARAMETER(BD=1.4011)
C         ...BD = H2 Bond Distance in Atomic Units (Bohrs)
C
      X1= 0.00D+00
      Y1= 0.00D+00
      Z1= BD/2.0D+00
      X2= 0.00D+00
      Y2= 0.00D+00
      Z2= -BD/2.0D+00
C      rab= 1.4D+00
C     -----------------------------------------------------
C     INTERNUCLEAR DISTANCE:  BOND = R = 
C                             Sqr(x**2 + y**2 + z**2)
C     -----------------------------------------------------
C      BOND=((X1-X2)**2 +(Y1-Y2)**2 +(Z1-Z2)**2)**(0.5)
C      WRITE(6,707) BOND
C      WRITE(6,709)
C      WRITE(7,707) BOND
C      WRITE(7,709)
C     -----------------------------------------------------
C         ...DEVICE = 6 ---> CRT Screen
C            DEVICE = 7 ---> OUTPUT Text File
C
C
C     -----------------------------------------------------
C     --- OutPut FORMAT(s) ---
C     -----------------------------------------------------
707   FORMAT(/,
     1    ' Absolute Bond Distance (PARAMATER) =',F13.9, 'au')
709   FORMAT(/)   
C     -----------------------------------------------------
C     -----------------------------------------------------
C
C
      RETURN 
      END
