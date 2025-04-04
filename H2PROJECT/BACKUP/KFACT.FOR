C*******************************************************
      SUBROUTINE KFACT(AX,AY,AZ,A1,BX,BY,BZ,A2,                      
     1                 CX,CY,CZ,A3,DX,DY,DZ,A4,
     2                                         A5,A6,XK) 
C*******************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C ***                                        ***
C ***    SUBROUTINE K CALCULATES K-FACTOR    ***
C ***                                        *** 
      ALK1=A1+A2                                                       
      ALK2=A3+A4                                                       
      BLK12=A5+A6  
      DETAB=ALK1*ALK2+(ALK1+ALK2)*BLK12
C
C     --- CENTERS ---
C
      XK1=EXP(-((A1*A2)/ALK1)*((AX-BX)**2+(AY-BY)**2+(AZ-BZ)**2))
      XK2=EXP(-((A3*A4)/ALK2)*((CX-DX)**2+(CY-DY)**2+(CZ-DZ)**2))
C
      RX1=(A1*AX+A2*BX)/ALK1
      RX2=(A3*CX+A4*DX)/ALK2
      RY1=(A1*AY+A2*BY)/ALK1
      RY2=(A3*CY+A4*DY)/ALK2
      RZ1=(A1*AZ+A2*BZ)/ALK1
      RZ2=(A3*CZ+A4*DZ)/ALK2
      DD=(RX1-RX2)**2+(RY1-RY2)**2+(RZ1-RZ2)**2
      XK=XK1*XK2*EXP(-(ALK1*ALK2*BLK12)*DD/DETAB)
      RETURN                                                            
      END       
