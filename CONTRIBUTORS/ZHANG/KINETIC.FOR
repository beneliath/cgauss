      SUBROUTINE KINETIC(ISYMM,M,OVV1,OVV2,OVV3,OVV4,
	-					AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,
	-					KINET)
C     -----------------------------------------------------
      IMPLICIT NONE
	INTEGER L,K,M,ISYMM
	REAL*8 OV1,OV2,OV3,OV4
	REAL*8 T1,T2,T3,T4,T5,T6,T7,T8
      REAL*8, DIMENSION(M) ::	AA1,RX1,RY1,RZ1,
	-						AA2,RX2,RY2,RZ2,B12
	REAL*8, DIMENSION(M,M) :: KINET
	REAL*8, DIMENSION(M,M) :: OVV1,OVV2,OVV3,OVV4
C
	PARAMETER(TWO=2.00D+00)
C
	DO L=1,M
		DO K=1,L
			OV1=OVV1(L,K)
			OV2=OVV2(L,K)
			OV3=OVV3(L,K)
			OV4=OVV4(L,K)
C
C     -----------------------------------------------------
C     --- Calculate KINETIC ENERGY OF ELECTRON 1 ---
C     -----------------------------------------------------
              CALL TKINET(OV1,
     -            RX1(L),RY1(L),RZ1(L),AA1(L),
     -            RX1(K),RY1(K),RZ1(K),AA1(K),                      
     -            RX2(L),RY2(L),RZ2(L),AA2(L),
     -            RX2(K),RY2(K),RZ2(K),AA2(K),
     -            B12(L),B12(K),T1)  
C
	IF (ISYMM) THEN
C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C
              CALL TKINET(OV2,
     -            RX1(L),RY1(L),RZ1(L),AA1(L),
     -            RX2(K),RY2(K),RZ2(K),AA2(K),                      
     -            RX2(L),RY2(L),RZ2(L),AA2(L),
     -            RX1(K),RY1(K),RZ1(K),AA1(K),
     -            B12(L),B12(K),T2)  
C
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
C
              CALL TKINET(OV3,
     -            RX1(L),RY1(L),RZ1(L),AA1(L),
     -            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
     -            RX2(L),RY2(L),RZ2(L),AA2(L),
     -            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
     -            B12(L),B12(K),T3)
C
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL TKINET(OV4,
     -            RX1(L),RY1(L),RZ1(L),AA1(L),
     -            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
     -            RX2(L),RY2(L),RZ2(L),AA2(L),
     -            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
     -            B12(L),B12(K),T4)
	END IF
C     -----------------------------------------------------
C     --- Calculate KINETIC ENERGY OF ELECTRON 2 ---
C     -----------------------------------------------------
             CALL TKINET(OV1,
     -            RX2(L),RY2(L),RZ2(L),AA2(L),
     -            RX2(K),RY2(K),RZ2(K),AA2(K),                      
     -            RX1(L),RY1(L),RZ1(L),AA1(L),
     -            RX1(K),RY1(K),RZ1(K),AA1(K),
     -            B12(L),B12(K),T5)  
C
	IF (ISYMM) THEN
C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C
              CALL TKINET(OV2,
     -            RX2(L),RY2(L),RZ2(L),AA2(L),
     -            RX1(K),RY1(K),RZ1(K),AA1(K),                      
     -            RX1(L),RY1(L),RZ1(L),AA1(L),
     -            RX2(K),RY2(K),RZ2(K),AA2(K),
     -            B12(L),B12(K),T6) 
C
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
C
              CALL TKINET(OV3,
     -            RX2(L),RY2(L),RZ2(L),AA2(L),
     -            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
     -            RX1(L),RY1(L),RZ1(L),AA1(L),
     -            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
     -            B12(L),B12(K),T7)
C
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL TKINET(OV4,
     -            RX2(L),RY2(L),RZ2(L),AA2(L),
     -            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
     -            RX1(L),RY1(L),RZ1(L),AA1(L),
     -            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
     -            B12(L),B12(K),T8)
	END IF
C
	IF (ISYMM) THEN
            KINET(L,K)=(T1+T2+T3+T4+T5+T6+T7+T8)/TWO
		ELSE
		  KINET(L,K)=(T1+T5)/TWO
	END IF
C				
	IF (L .NE. K) KINET(K,L)=KINET(L,K)
C
		END DO
	END DO
C
	RETURN
	END