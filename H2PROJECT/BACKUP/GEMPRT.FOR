C     Originallly Coded:  Z. Zhang
C           Modified By:  D. Gilmore ---> (24 May 95)
C     -----------------------------------------------------
C     L. ADAMOWICZ RESEARCH / UNIVERSITY OF ARIZONA / CHEM.
C     -----------------------------------------------------
      SUBROUTINE GEMPRT(NN,
     1                  AA,X1,X2,X3,
     2                  BB,Y1,Y2,Y3,
     3                  CC)
C     -----------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL PRINT /.TRUE./
C         ...The Variable "PRINT" is of type LOGICAL and is
C            set equal to "TRUE".
C
      WRITE(6,999) NN,AA,X1,X2,X3,BB,Y1,Y2,Y3,CC
      WRITE(7,999) NN,AA,X1,X2,X3,BB,Y1,Y2,Y3,CC
C         ...DEVICE = 6 ---> CRT Screen
C            DEVICE = 7 ---> OUTPUT Text File
C
C     NN = Gaussian Geminal Number
C     AA = Exponent Alpha _sub(1)
C     X1, X2, X3 = Cart. Coordinates for Gaussian of AA
C     BB = Exponent Alpha _sub(2)
C     Y1, Y2, Y3 = Cart. Coordinates for Gaussian of BB
C     CC = Correleation Exponent
C
C
C     -----------------------------------------------------
C     --- OutPut FORMAT(s) ---
C     -----------------------------------------------------
999   FORMAT(' ----------------------------------------------',
     1/,' GEMINAL NO.=',I3,/,
C         ...I3 = Three Digit Integer (here using NN)
C
     2' ----------------------------------------------',/,
     3' Orbital exponent ALPHA1 =',F13.9,/,
C         ...F13.9 = Thirteen Digit Floating Point =
C                     '####.#########'
C 
     4' Rx1 = ',F13.9,'au',/,' Ry1 = ',F13.9,'au',/,' Rz1 = ',
     5 F13.9,'au',/,/, 
     6' Orbital exponent ALPHA2 =',F13.9,/, 
     7' Rx2 = ',F13.9,'au',/,' Ry2 = ',F13.9,'au',/,' Rz2 = ',
     8 F13.9,'au',/,/,' Correl. exponent BETA12 =',F13.9, 
     9 /,' ==============================================',/)
C
C     ...SEE OUTPUT Text File to see how this prints out...
C     -----------------------------------------------------
C     -----------------------------------------------------
C
C
      RETURN
      END