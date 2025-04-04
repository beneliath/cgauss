SUBROUTINE ENUCLEAR(OV,AZ,A1,BZ,A2,CZ,A3,DZ,A4,A5,A6,	&
                      EZ,EN)                       

	USE scalars_to_be_shared

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	IMPLICIT INTEGER (I-N)


!*** FOUR CENTER NUCLEAR ATTRACTION INTEGRAL ***
!-----------------------------------------------------

      ALK1=A1+A2                                                       
      ALK2=A3+A4                                                       
      BLK12=A5+A6  
      ALKXX=ALK2+BLK12
      DETAB=ALK1*ALK2+(ALK1+ALK2)*BLK12


!--- ONE CENTER NUCLEAR INTEGRAL ---
!-----------------------------------------------------

!      RX1=(A1*AX+A2*BX)/ALK1
!      RX2=(A3*CX+A4*DX)/ALK2
!      RY1=(A1*AY+A2*BY)/ALK1
!      RY2=(A3*CY+A4*DY)/ALK2

      RZ1=(A1*AZ+A2*BZ)/ALK1
      RZ2=(A3*CZ+A4*DZ)/ALK2


!--- FOUR-CENTER ELECTRON REPULSION INTEGRAL ---
!-----------------------------------------------------

	EN=TWO*OV*DSQRT(DETAB/(PI*(ALK2+BLK12)))


!      RRX=(ALK1*(ALK2+BLK12)*RX1+ALK2*BLK12*RX2)/DETAB 
!      RRY=(ALK1*(ALK2+BLK12)*RY1+ALK2*BLK12*RY2)/DETAB 

      RRZ=(ALK1*(ALK2+BLK12)*RZ1+ALK2*BLK12*RZ2)/DETAB
	   

!      RE2=(RRX-EX)**TWO +(RRY-EY)**TWO +(RRZ-EZ)**TWO

      RE2=(RRZ-EZ)*(RRZ-EZ)

      DELTA2=(ALK2+BLK12)/(FOUR*DETAB)

      RE2=RE2/(FOUR*DELTA2)

      EN=EN*F0(RE2) 


END SUBROUTINE ENUCLEAR