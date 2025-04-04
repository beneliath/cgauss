C     Originallly Coded:  Z. Zhang
C           Modified By:  D. Gilmore ---> (24 May 95)
C     -----------------------------------------------------
C     L. ADAMOWICZ RESEARCH / UNIVERSITY OF ARIZONA / CHEM.
C     -----------------------------------------------------
      SUBROUTINE KFACT(AX,AY,AZ,A1,BX,BY,BZ,A2,                      
     1                 CX,CY,CZ,A3,DX,DY,DZ,A4,
     2                 A5,A6,XK) 
C     -----------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
	COMMON /OVER/ DETS(50),L,K,CC(50)

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
C     XK = ?
C     -----------------------------------------------------
C     --- SUBROUTINE KFACT CALCULATES K-FACTOR ---
C     -----------------------------------------------------
      ALK1=A1+A2                                                       
      ALK2=A3+A4                                                       
      BLK12=A5+A6  
      DETAB=ALK1*ALK2+((ALK1+ALK2)*BLK12)
C      DETAB=(ALK1*ALK2)+((ALK1+ALK2)*(2.0D0*BLK12))
C
C     -----------------------------------------------------
C     --- CENTERS ---
C      -----------------------------------------------------
C
c      XK1=EXP(-((A1*A2)/ALK1)*((AX-BX)**2+(AY-BY)**2+(AZ-BZ)**2))
c      XK2=EXP(-((A3*A4)/ALK2)*((CX-DX)**2+(CY-DY)**2+(CZ-DZ)**2))
C
      RX1=(A1*AX+A2*BX)/ALK1
      RX2=(A3*CX+A4*DX)/ALK2
      RY1=(A1*AY+A2*BY)/ALK1
      RY2=(A3*CY+A4*DY)/ALK2
      RZ1=(A1*AZ+A2*BZ)/ALK1
      RZ2=(A3*CZ+A4*DZ)/ALK2
      DD=((RX1-RX2)**2)+((RY1-RY2)**2)+((RZ1-RZ2)**2)
C      XK=XK1*XK2*EXP(-((ALK1*ALK2*BLK12)/DETAB)*DD)
C	WRITE(6,*)DETAB,DD,ALK1,ALK2,BLK12,XK1,XK2
c       XK=XK1*XK2*EXP(-((ALK1*ALK2*BLK12)*DD/DETAB)+((CC(K)+CC
c	1(L))/2.0D+00))
	CKL=(((A1*A2)/ALK1)*((AX-BX)**2+(AY-BY)**2+(AZ-BZ)**2)+((A3*A4
	1)/ALK2)*((CX-DX)**2+(CY-DY)**2+(CZ-DZ)**2)+((ALK1*ALK2*BLK12)*DD/D
     2ETAB))

	XK=EXP(-CKL+((CC(K)+CC(L))/2.0D+00))

C         ...Need to Verify the EXP (...) Part Here.
C
C     -----------------------------------------------------
C
      RETURN                                                            
      END       
