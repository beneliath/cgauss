      SUBROUTINE POTENTIAL(ISYMM,M,OVV1,OVV2,OVV3,OVV4,S
	-					EX1,EY1,EZ1,EX2,EY2,EZ2,
	-					AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,
	-					POTEN)
C     -----------------------------------------------------
      IMPLICIT NONE
	INTEGER L,K,M,ISYMM
	REAL*8 OV1,OV2,OV3,OV4
	REAL*8 EN1,EN2,EN3,EN4,EN5,EN6,EN7,EN8,EN9,EN10,EN11,EN12,EN13,
	-		EN14,EN15,EN16
      REAL*8 ENN1
      REAL*8, DIMENSION(M) ::	AA1,RX1,RY1,RZ1,
	-						AA2,RX2,RY2,RZ2,B12
	REAL*8, DIMENSION(M,M) :: POTEN,S
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
C     --- Calculate Nuc./Elec. Attractn. FOR ELECTRON 1 ---
C     -----------------------------------------------------
              CALL ENUCLEAR(OV1,RX1(L),RY1(L),RZ1(L),AA1(L),
     -           RX1(K),RY1(K),RZ1(K),AA1(K),                      
     -           RX2(L),RY2(L),RZ2(L),AA2(L),
     -           RX2(K),RY2(K),RZ2(K),AA2(K),
     -           B12(L),B12(K),
     -           EX1,EY1,EZ1,EN1)
C
	IF (ISYMM) THEN
C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C
              CALL ENUCLEAR(OV2,RX1(L),RY1(L),RZ1(L),AA1(L),
     -           RX2(K),RY2(K),RZ2(K),AA2(K),                      
     -           RX2(L),RY2(L),RZ2(L),AA2(L),
     -           RX1(K),RY1(K),RZ1(K),AA1(K),
     -           B12(L),B12(K),
     -           EX1,EY1,EZ1,EN2)                       
C
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
C
              CALL ENUCLEAR(OV3,RX1(L),RY1(L),RZ1(L),AA1(L),
     -            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
     -            RX2(L),RY2(L),RZ2(L),AA2(L),
     -            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
     -            B12(L),B12(K),
     -            EX1,EY1,EZ1,EN3)
C
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL ENUCLEAR(OV4,RX1(L),RY1(L),RZ1(L),AA1(L),
     -            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
     -            RX2(L),RY2(L),RZ2(L),AA2(L),
     -            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
     -            B12(L),B12(K),
     -            EX1,EY1,EZ1,EN4)
	END IF
C  
              CALL ENUCLEAR(OV1,RX1(L),RY1(L),RZ1(L),AA1(L),
     -           RX1(K),RY1(K),RZ1(K),AA1(K),                      
     -           RX2(L),RY2(L),RZ2(L),AA2(L),
     -           RX2(K),RY2(K),RZ2(K),AA2(K),
     -           B12(L),B12(K),
     -           EX2,EY2,EZ2,EN5)
C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C
	IF (ISYMM) THEN
              CALL ENUCLEAR(OV2,RX1(L),RY1(L),RZ1(L),AA1(L),
     -           RX2(K),RY2(K),RZ2(K),AA2(K),                      
     -           RX2(L),RY2(L),RZ2(L),AA2(L),
     -           RX1(K),RY1(K),RZ1(K),AA1(K),
     -           B12(L),B12(K),
     -           EX2,EY2,EZ2,EN6)                       
C
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
C
              CALL ENUCLEAR(OV3,RX1(L),RY1(L),RZ1(L),AA1(L),
     -            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
     -            RX2(L),RY2(L),RZ2(L),AA2(L),
     -            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
     -            B12(L),B12(K),
     -            EX2,EY2,EZ2,EN7)
C
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL ENUCLEAR(OV4,RX1(L),RY1(L),RZ1(L),AA1(L),
     -            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
     -            RX2(L),RY2(L),RZ2(L),AA2(L),
     -            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
     -            B12(L),B12(K),
     -            EX2,EY2,EZ2,EN8)
	END IF
C 
C     -----------------------------------------------------
C     --- Calculate Nuc./Elec. Attractn. FOR ELECTRON 2 ---
C     -----------------------------------------------------
              CALL ENUCLEAR(OV1,RX2(L),RY2(L),RZ2(L),AA2(L),
     -           RX2(K),RY2(K),RZ2(K),AA2(K),                      
     -           RX1(L),RY1(L),RZ1(L),AA1(L),
     -           RX1(K),RY1(K),RZ1(K),AA1(K),
     -           B12(L),B12(K),
     -           EX1,EY1,EZ1,EN9)                       
C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C
	IF (ISYMM) THEN
              CALL ENUCLEAR(OV2,RX2(L),RY2(L),RZ2(L),AA2(L),
     -           RX1(K),RY1(K),RZ1(K),AA1(K),                      
     -           RX1(L),RY1(L),RZ1(L),AA1(L),
     -           RX2(K),RY2(K),RZ2(K),AA2(K),
     -           B12(L),B12(K),
     -           EX1,EY1,EZ1,EN10)                       
C
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
C
              CALL ENUCLEAR(OV3,RX2(L),RY2(L),RZ2(L),AA2(L),
     -            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
     -            RX1(L),RY1(L),RZ1(L),AA1(L),
     -            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
     -            B12(L),B12(K),
     -            EX1,EY1,EZ1,EN11)
C
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL ENUCLEAR(OV4,RX2(L),RY2(L),RZ2(L),AA2(L),
     -            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
     -            RX1(L),RY1(L),RZ1(L),AA1(L),
     -            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
     -            B12(L),B12(K),
     -            EX1,EY1,EZ1,EN12)
	END IF
C 
              CALL ENUCLEAR(OV1,RX2(L),RY2(L),RZ2(L),AA2(L),
     -           RX2(K),RY2(K),RZ2(K),AA2(K),                      
     -           RX1(L),RY1(L),RZ1(L),AA1(L),
     -           RX1(K),RY1(K),RZ1(K),AA1(K),
     -           B12(L),B12(K),
     -           EX2,EY2,EZ2,EN13)                       
C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C
	IF (ISYMM) THEN
              CALL ENUCLEAR(OV2,RX2(L),RY2(L),RZ2(L),AA2(L),
     -           RX1(K),RY1(K),RZ1(K),AA1(K),                      
     -           RX1(L),RY1(L),RZ1(L),AA1(L),
     -           RX2(K),RY2(K),RZ2(K),AA2(K),
     -           B12(L),B12(K),
     -           EX2,EY2,EZ2,EN14)  
C
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
C
              CALL ENUCLEAR(OV3,RX2(L),RY2(L),RZ2(L),AA2(L),
     -            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
     -            RX1(L),RY1(L),RZ1(L),AA1(L),
     -            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
     -            B12(L),B12(K),
     -            EX2,EY2,EZ2,EN15)
C
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL ENUCLEAR(OV4,RX2(L),RY2(L),RZ2(L),AA2(L),
     -            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
     -            RX1(L),RY1(L),RZ1(L),AA1(L),
     -            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
     -            B12(L),B12(K),
     -            EX2,EY2,EZ2,EN16)
	END IF
C 
C
	IF (ISYMM) THEN
              ENN1=EN1+EN2+EN3+EN4+EN5+EN6+EN7+EN8
     -            +EN9+EN10+EN11+EN12+EN13+EN14+EN15+EN16
	ELSE
              ENN1=EN1+EN5+EN9+EN13
	END IF
C     -----------------------------------------------------
C     --- Calculate the Elec./Elec. Repulsion ---
C     -----------------------------------------------------
              CALL ELREP(DFLAG,OV1,RX1(L),RY1(L),RZ1(L),AA1(L),
     1           RX1(K),RY1(K),RZ1(K),AA1(K),                      
     2           RX2(L),RY2(L),RZ2(L),AA2(L),
     3           RX2(K),RY2(K),RZ2(K),AA2(K),
     4           B12(L),B12(K),ER1,KLM)  
C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C
	IF (ISYMM) THEN
              CALL ELREP(DFLAG,OV2,RX1(L),RY1(L),RZ1(L),AA1(L),
     1           RX2(K),RY2(K),RZ2(K),AA2(K),                      
     2           RX2(L),RY2(L),RZ2(L),AA2(L),
     3           RX1(K),RY1(K),RZ1(K),AA1(K),
     4           B12(L),B12(K),ER2,KLM)  
C
C     ...Nuclear Symmetrization Overlap Integral: P_[H1,H2]
C
              CALL ELREP(DFLAG,OV3,RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
     4            B12(L),B12(K),ER3,KLM)
C
C     ...Electro/Nuclear Symmetrization Overlap Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL ELREP(DFLAG,OV4,RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
     4            B12(L),B12(K),ER4,KLM)
	END IF

C              
	IF (ISYMM) THEN
		EER1=ER1+ER2+ER3+ER4
	ELSE
		EER1=ER1
	END IF
C
C
	IF (ISYMM) THEN
              EN=EN1+EN2+EN3+EN4+EN5+EN6+EN7+EN8
     1            +EN9+EN10+EN11+EN12+EN13+EN14+EN15+EN16
	ELSE
              EN=EN1+EN5+EN9+EN13
	END IF
C
	IF(IPERM.EQ.1) THEN
              EN=EN1+EN5+EN9+EN13
	END IF
	IF(IPERM.EQ.2) THEN
              EN=EN2+EN6+EN10+EN14
	END IF
	IF(IPERM.EQ.3) THEN
              EN=EN3+EN7+EN11+EN15
	END IF
	IF(IPERM.EQ.4) THEN
              EN=EN4+EN8+EN12+EN16
	END IF
C
	IF (ISYMM) THEN
              ER=ER1+ER2+ER3+ER4
	ELSE
			ER=ER1
	END IF
C
	IF(IPERM.EQ.1) THEN
			ER=ER1
	END IF
	IF(IPERM.EQ.2) THEN
			ER=ER2
	END IF
	IF(IPERM.EQ.3) THEN
			ER=ER3
	END IF
	IF(IPERM.EQ.4) THEN
			ER=ER4
	END IF
C
		END DO
	END DO
C
	RETURN
	END