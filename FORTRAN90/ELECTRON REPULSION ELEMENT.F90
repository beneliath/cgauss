SUBROUTINE ELECREP(OV,AZ,A1,BZ,A2,CZ,A3,DZ,A4,A5,A6,ER)                       

	USE scalars_to_be_shared

    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	IMPLICIT INTEGER (I-N)


!*** FOUR CENTER ELECTRON REPULSION INTEGRAL ***
!-----------------------------------------------------

      ALK1=A1+A2                                                       
      ALK2=A3+A4                                                       
      BLK12=A5+A6
      ALKXX=ALK2+BLK12
      DETAB=ALK1*ALK2+(ALK1+ALK2)*BLK12


!--- ONE CENTER ELECTRON REPULSION INTEGRAL ---
!-----------------------------------------------------

	ER=TWO*OV*DSQRT(DETAB/(PI*(ALK1+ALK2)))


!--- CENTERS ---
!-----------------------------------------------------

!     RX1=(A1*AX+A2*BX)/ALK1
!     RX2=(A3*CX+A4*DX)/ALK2
!     RY1=(A1*AY+A2*BY)/ALK1
!     RY2=(A3*CY+A4*DY)/ALK2

!      RX1=0.00D+00
!      RX2=0.00D+00
!      RY1=0.00D+00
!      RY2=0.00D+00

      RZ1=(A1*AZ+A2*BZ)/ALK1
      RZ2=(A3*CZ+A4*DZ)/ALK2

!      DD=(RX1-RX2)**2+(RY1-RY2)**2+(RZ1-RZ2)**2

      D0=(RZ1-RZ2)
	  DD=D0*D0


!--- FOUR-CENTER ELECTRON REPULSION INTEGRAL ---
!-----------------------------------------------------

      WW=ALK1*ALK1*ALK2*ALK2*DD/((ALK1+ALK2)*DETAB) 


	ER=ER*F0(WW)
	  

END SUBROUTINE ELECREP
