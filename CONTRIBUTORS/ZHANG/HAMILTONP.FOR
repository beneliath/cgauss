C     L. ADAMOWICZ RESEARCH / UNIVERSITY OF ARIZONA / CHEM.
C     -----------------------------------------------------
      SUBROUTINE HAMILTONp(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,
     1                    EVEC,S,H,EMIN,XAA1,XRX1,XRY1,XRZ1,
     2					XAA2,XRX2,XRY2,XRZ2,XB12)
C     -----------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER M,ERROR,DFLAG
	REAL*8 Snp1,Hnp1
      REAL*8, DIMENSION(M) ::	AA1,RX1,RY1,RZ1,
	-						AA2,RX2,RY2,RZ2,
     -						B12,EVEC
	REAL*8, DIMENSION(M,M) :: S,H
	REAL*8, DIMENSION(M+1) :: CO
	REAL*8 HH(2,2),SS(2,2),U(2,2),BIG(2),JB(2),VH(3)
C     -----------------------------------------------------
	COMMON /OVER/ L,K
	COMMON /TNOVER/ NXAC,XAC
C     -----------------------------------------------------
	PARAMETER(ISYMM=1)
	PARAMETER(IPERM=0)
C     -----------------------------------------------------
      MMAX=M
	DFLAG=0
C     -----------------------------------------------------
      DATA ONE,TWO,FOUR,SIX,EIGHT,HALF/1.,2.,4.,6.,8.,.5/
C     -----------------------------------------------------
      CALL GEOM(EX1,EY1,EZ1,EX2,EY2,EZ2)
	L=1
	K=1
C     -----------------------------------------------------	  
C     -----------------------------------------------------
C     --- Calculate OVERLAP INTEGRAL ---
C     -----------------------------------------------------
             CALL OVERLP(M,DFLAG,XRX1,XRY1,XRZ1,XAA1,
     1            XRX1,XRY1,XRZ1,XAA1,                      
     2            XRX2,XRY2,XRZ2,XAA2,
     3            XRX2,XRY2,XRZ2,XAA2,
     4            XB12,XB12,OV1,KLM)
C
C
	IF(ISYMM) THEN
		CALL OVERLP(M,DFLAG,XRX1,XRY1,XRZ1,XAA1,
     1            XRX2,XRY2,XRZ2,XAA2,                      
     2            XRX2,XRY2,XRZ2,XAA2,
     3            XRX1,XRY1,XRZ1,XAA1,
     4            XB12,XB12,OV2,KLM)
C
C
C     ...Nuclear Symmetrization Overlap Integral: P_[H1,H2]
C
              CALL OVERLP(M,DFLAG,XRX1,XRY1,XRZ1,XAA1,
     1            -XRX1,-XRY1,-XRZ1,XAA1,                      
     2            XRX2,XRY2,XRZ2,XAA2,
     3            -XRX2,-XRY2,-XRZ2,XAA2,
     4            XB12,XB12,OV3,KLM)
C
C
C     ...Electro/Nuclear Symmetrization Overlap Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL OVERLP(M,DFLAG,XRX1,XRY1,XRZ1,XAA1,
     1            -XRX2,-XRY2,-XRZ2,XAA2,                      
     2            XRX2,XRY2,XRZ2,XAA2,
     3            -XRX1,-XRY1,-XRZ1,XAA1,
     4            XB12,XB12,OV4,KLM)
C              
	END IF
C     -----------------------------------------------------
	IF(ISYMM) THEN
              Snp1=OV1+OV2+OV3+OV4
		ELSE
              Snp1=OV1
	END IF
	IF(IPERM.EQ.1) THEN
		Snp1=OV1
	END IF
	IF(IPERM.EQ.2) THEN
		Snp1=OV2
	END IF
	IF(IPERM.EQ.3) THEN
		Snp1=OV3
	END IF
	IF(IPERM.EQ.4) THEN
		Snp1=OV4
	END IF
C
C     -----------------------------------------------------
C     --- Calculate KINETIC ENERGY OF ELECTRON 1 ---
C     -----------------------------------------------------
              CALL TK(DFLAG,OV1,
     1            XRX1,XRY1,XRZ1,XAA1,
     2            XRX1,XRY1,XRZ1,XAA1,                      
     3            XRX2,XRY2,XRZ2,XAA2,
     4            XRX2,XRY2,XRZ2,XAA2,
     5            XB12,XB12,T1)  
C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
	IF (ISYMM) THEN
              CALL TK(DFLAG,OV2,
     1            XRX1,XRY1,XRZ1,XAA1,
     2            XRX2,XRY2,XRZ2,XAA2,                      
     3            XRX2,XRY2,XRZ2,XAA2,
     4            XRX1,XRY1,XRZ1,XAA1,
     5            XB12,XB12,T2)  
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
              CALL TK(DFLAG,OV3,
     1            XRX1,XRY1,XRZ1,XAA1,
     2            -XRX1,-XRY1,-XRZ1,XAA1,                      
     3            XRX2,XRY2,XRZ2,XAA2,
     4            -XRX2,-XRY2,-XRZ2,XAA2,
     5            XB12,XB12,T3)
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
              CALL TK(DFLAG,OV4,
     1            XRX1,XRY1,XRZ1,XAA1,
     2            -XRX2,-XRY2,-XRZ2,XAA2,                      
     3            XRX2,XRY2,XRZ2,XAA2,
     4            -XRX1,-XRY1,-XRZ1,XAA1,
     5            XB12,XB12,T4)
	END IF
C     -----------------------------------------------------
C     --- Calculate KINETIC ENERGY OF ELECTRON 2 ---
C     -----------------------------------------------------
             CALL TK(DFLAG,OV1,
     1            XRX2,XRY2,XRZ2,XAA2,
     2            XRX2,XRY2,XRZ2,XAA2,                      
     3            XRX1,XRY1,XRZ1,XAA1,
     4            XRX1,XRY1,XRZ1,XAA1,
     5            XB12,XB12,T5)  
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
	IF (ISYMM) THEN
              CALL TK(DFLAG,OV2,
     1            XRX2,XRY2,XRZ2,XAA2,
     2            XRX1,XRY1,XRZ1,XAA1,                      
     3            XRX1,XRY1,XRZ1,XAA1,
     4            XRX2,XRY2,XRZ2,XAA2,
     5            XB12,XB12,T6) 
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
              CALL TK(DFLAG,OV3,
     1            XRX2,XRY2,XRZ2,XAA2,
     2            -XRX2,-XRY2,-XRZ2,XAA2,                      
     3            XRX1,XRY1,XRZ1,XAA1,
     4            -XRX1,-XRY1,-XRZ1,XAA1,
     5            XB12,XB12,T7)
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
              CALL TK(DFLAG,OV4,
     1            XRX2,XRY2,XRZ2,XAA2,
     2            -XRX1,-XRY1,-XRZ1,XAA1,                      
     3            XRX1,XRY1,XRZ1,XAA1,
     4            -XRX2,-XRY2,-XRZ2,XAA2,
     5            XB12,XB12,T8)
	END IF
C
	IF (ISYMM) THEN
            T=(T1+T2+T3+T4+T5+T6+T7+T8)/2.0D+00
	ELSE
		  T=(T1+T5)/2.0D+00
	END IF
	IF(IPERM.EQ.1) THEN
		T=(T1+T5)/2.0D+00
	END IF
	IF(IPERM.EQ.2) THEN
		T=(T2+T6)/2.0D+00
	END IF
	IF(IPERM.EQ.3) THEN
		T=(T3+T7)/2.0D+00
	END IF
	IF(IPERM.EQ.4) THEN
		T=(T4+T8)/2.0D+00
	END IF
C     -----------------------------------------------------
C     --- Calculate Nuc./Elec. Attractn. FOR ELECTRON 1 ---
C     -----------------------------------------------------
              CALL ENUCL(DFLAG,OV1,XRX1,XRY1,XRZ1,XAA1,
     1           XRX1,XRY1,XRZ1,XAA1,                      
     2           XRX2,XRY2,XRZ2,XAA2,
     3           XRX2,XRY2,XRZ2,XAA2,
     4           XB12,XB12,
     5           EX1,EY1,EZ1,EN1,KLM)
C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C
	IF (ISYMM) THEN
              CALL ENUCL(DFLAG,OV2,XRX1,XRY1,XRZ1,XAA1,
     1           XRX2,XRY2,XRZ2,XAA2,                      
     2           XRX2,XRY2,XRZ2,XAA2,
     3           XRX1,XRY1,XRZ1,XAA1,
     4           XB12,XB12,
     5           EX1,EY1,EZ1,EN2,KLM)                       
C
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
C
              CALL ENUCL(DFLAG,OV3,XRX1,XRY1,XRZ1,XAA1,
     1            -XRX1,-XRY1,-XRZ1,XAA1,                      
     2            XRX2,XRY2,XRZ2,XAA2,
     3            -XRX2,-XRY2,-XRZ2,XAA2,
     4            XB12,XB12,
     5            EX1,EY1,EZ1,EN3,KLM)
C
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL ENUCL(DFLAG,OV4,XRX1,XRY1,XRZ1,XAA1,
     1            -XRX2,-XRY2,-XRZ2,XAA2,                      
     2            XRX2,XRY2,XRZ2,XAA2,
     3            -XRX1,-XRY1,-XRZ1,XAA1,
     4            XB12,XB12,
     5            EX1,EY1,EZ1,EN4,KLM)
	END IF
C  
              CALL ENUCL(DFLAG,OV1,XRX1,XRY1,XRZ1,XAA1,
     1           XRX1,XRY1,XRZ1,XAA1,                      
     2           XRX2,XRY2,XRZ2,XAA2,
     3           XRX2,XRY2,XRZ2,XAA2,
     4           XB12,XB12,
     5           EX2,EY2,EZ2,EN5,KLM)
C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C
	IF (ISYMM) THEN
              CALL ENUCL(DFLAG,OV2,XRX1,XRY1,XRZ1,XAA1,
     1           XRX2,XRY2,XRZ2,XAA2,                      
     2           XRX2,XRY2,XRZ2,XAA2,
     3           XRX1,XRY1,XRZ1,XAA1,
     4           XB12,XB12,
     5           EX2,EY2,EZ2,EN6,KLM)                       
C
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
C
              CALL ENUCL(DFLAG,OV3,XRX1,XRY1,XRZ1,XAA1,
     1            -XRX1,-XRY1,-XRZ1,XAA1,                      
     2            XRX2,XRY2,XRZ2,XAA2,
     3            -XRX2,-XRY2,-XRZ2,XAA2,
     4            XB12,XB12,
     5            EX2,EY2,EZ2,EN7,KLM)
C
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL ENUCL(DFLAG,OV4,XRX1,XRY1,XRZ1,XAA1,
     1            -XRX2,-XRY2,-XRZ2,XAA2,                      
     2            XRX2,XRY2,XRZ2,XAA2,
     3            -XRX1,-XRY1,-XRZ1,XAA1,
     4            XB12,XB12,
     5            EX2,EY2,EZ2,EN8,KLM)
	END IF
C 
C     -----------------------------------------------------
C     --- Calculate Nuc./Elec. Attractn. FOR ELECTRON 2 ---
C     -----------------------------------------------------
              CALL ENUCL(DFLAG,OV1,XRX2,XRY2,XRZ2,XAA2,
     1           XRX2,XRY2,XRZ2,XAA2,                      
     2           XRX1,XRY1,XRZ1,XAA1,
     3           XRX1,XRY1,XRZ1,XAA1,
     4           XB12,XB12,
     5           EX1,EY1,EZ1,EN9,KLM)                       
C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C
	IF (ISYMM) THEN
              CALL ENUCL(DFLAG,OV2,XRX2,XRY2,XRZ2,XAA2,
     1           XRX1,XRY1,XRZ1,XAA1,                      
     2           XRX1,XRY1,XRZ1,XAA1,
     3           XRX2,XRY2,XRZ2,XAA2,
     4           XB12,XB12,
     5           EX1,EY1,EZ1,EN10,KLM)                       
C
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
C
              CALL ENUCL(DFLAG,OV3,XRX2,XRY2,XRZ2,XAA2,
     1            -XRX2,-XRY2,-XRZ2,XAA2,                      
     2            XRX1,XRY1,XRZ1,XAA1,
     3            -XRX1,-XRY1,-XRZ1,XAA1,
     4            XB12,XB12,
     5            EX1,EY1,EZ1,EN11,KLM)
C
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL ENUCL(DFLAG,OV4,XRX2,XRY2,XRZ2,XAA2,
     1            -XRX1,-XRY1,-XRZ1,XAA1,                      
     2            XRX1,XRY1,XRZ1,XAA1,
     3            -XRX2,-XRY2,-XRZ2,XAA2,
     4            XB12,XB12,
     5            EX1,EY1,EZ1,EN12,KLM)
	END IF
C 
              CALL ENUCL(DFLAG,OV1,XRX2,XRY2,XRZ2,XAA2,
     1           XRX2,XRY2,XRZ2,XAA2,                      
     2           XRX1,XRY1,XRZ1,XAA1,
     3           XRX1,XRY1,XRZ1,XAA1,
     4           XB12,XB12,
     5           EX2,EY2,EZ2,EN13,KLM)                       
C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C
	IF (ISYMM) THEN
              CALL ENUCL(DFLAG,OV2,XRX2,XRY2,XRZ2,XAA2,
     1           XRX1,XRY1,XRZ1,XAA1,                      
     2           XRX1,XRY1,XRZ1,XAA1,
     3           XRX2,XRY2,XRZ2,XAA2,
     4           XB12,XB12,
     5           EX2,EY2,EZ2,EN14,KLM)  
C
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
C
              CALL ENUCL(DFLAG,OV3,XRX2,XRY2,XRZ2,XAA2,
     1            -XRX2,-XRY2,-XRZ2,XAA2,                      
     2            XRX1,XRY1,XRZ1,XAA1,
     3            -XRX1,-XRY1,-XRZ1,XAA1,
     4            XB12,XB12,
     5            EX2,EY2,EZ2,EN15,KLM)
C
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL ENUCL(DFLAG,OV4,XRX2,XRY2,XRZ2,XAA2,
     1            -XRX1,-XRY1,-XRZ1,XAA1,                      
     2            XRX1,XRY1,XRZ1,XAA1,
     3            -XRX2,-XRY2,-XRZ2,XAA2,
     4            XB12,XB12,
     5            EX2,EY2,EZ2,EN16,KLM)
	END IF
C 
C
	IF (ISYMM) THEN
              ENN1=EN1+EN2+EN3+EN4+EN5+EN6+EN7+EN8
     1            +EN9+EN10+EN11+EN12+EN13+EN14+EN15+EN16
	ELSE
              ENN1=EN1+EN5+EN9+EN13
	END IF
C     -----------------------------------------------------
C     --- Calculate the Elec./Elec. Repulsion ---
C     -----------------------------------------------------
              CALL ELREP(DFLAG,OV1,XRX1,XRY1,XRZ1,XAA1,
     1           XRX1,XRY1,XRZ1,XAA1,                      
     2           XRX2,XRY2,XRZ2,XAA2,
     3           XRX2,XRY2,XRZ2,XAA2,
     4           XB12,XB12,ER1,KLM)  
C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C
	IF (ISYMM) THEN
              CALL ELREP(DFLAG,OV2,XRX1,XRY1,XRZ1,XAA1,
     1           XRX2,XRY2,XRZ2,XAA2,                      
     2           XRX2,XRY2,XRZ2,XAA2,
     3           XRX1,XRY1,XRZ1,XAA1,
     4           XB12,XB12,ER2,KLM)  
C
C     ...Nuclear Symmetrization Overlap Integral: P_[H1,H2]
C
              CALL ELREP(DFLAG,OV3,XRX1,XRY1,XRZ1,XAA1,
     1            -XRX1,-XRY1,-XRZ1,XAA1,                      
     2            XRX2,XRY2,XRZ2,XAA2,
     3            -XRX2,-XRY2,-XRZ2,XAA2,
     4            XB12,XB12,ER3,KLM)
C
C     ...Electro/Nuclear Symmetrization Overlap Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL ELREP(DFLAG,OV4,XRX1,XRY1,XRZ1,XAA1,
     1            -XRX2,-XRY2,-XRZ2,XAA2,                      
     2            XRX2,XRY2,XRZ2,XAA2,
     3            -XRX1,-XRY1,-XRZ1,XAA1,
     4            XB12,XB12,ER4,KLM)
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
              Hnp1=ER+T-EN+(1.0D+00/1.4011D+00)*Snp1
C     ----------------------------------------------------- 
C      CALL BIOUT(' ')
C	CALL ENERGY(M,H,S,EMIN,CO,KIM)
C	DO J=1,M
C		EVEC(J)=CO(J)
C	END DO
C	WRITE(6,*)CO(1),CO(2)
C      J=1
C	CALL BIOUT(' ')
C      CALL BIOUT('     LinCoeff:.   AA1:        AA2:        B12:        
C	-RZ1:       RZ2:')
C      CALL BIOUT('FNC: =================================================
C     -=====================')
C	DO J=1,M
C          WRITE(6,1111) J,CO(J),AA1(J),AA2(J),B12(J),RZ1(J),RZ2(J)
C          WRITE(7,1111) J,CO(J),AA1(J),AA2(J),B12(J),RZ1(J),RZ2(J)
C	END DO

	SS(1,1)=1.00D+00
	SS(2,2)=Snp1/Snp1
C	HH(1,1)=EMIN
      HH(2,2)=Hnp1/Snp1
	HH(1,2)=0.00D+00
	SS(1,2)=0.00D+00
C     ----------------------------------------------------- 
      IF (M.GE.2) K=2
	DO L=1,M
C     -----------------------------------------------------
C     --- Calculate OVERLAP INTEGRAL ---
C     -----------------------------------------------------
             CALL OVERLP(M,DFLAG,RX1(L),RY1(L),RZ1(L),AA1(L),
     1            XRX1,XRY1,XRZ1,XAA1,                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            XRX2,XRY2,XRZ2,XAA2,
     4            B12(L),XB12,OV1,KLM)
C
C
	IF(ISYMM) THEN
		CALL OVERLP(M,DFLAG,RX1(L),RY1(L),RZ1(L),AA1(L),
     1            XRX2,XRY2,XRZ2,XAA2,                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            XRX1,XRY1,XRZ1,XAA1,
     4            B12(L),XB12,OV2,KLM)
C
C
C     ...Nuclear Symmetrization Overlap Integral: P_[H1,H2]
C
              CALL OVERLP(M,DFLAG,RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -XRX1,-XRY1,-XRZ1,XAA1,                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -XRX2,-XRY2,-XRZ2,XAA2,
     4            B12(L),XB12,OV3,KLM)
C
C
C     ...Electro/Nuclear Symmetrization Overlap Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL OVERLP(M,DFLAG,RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -XRX2,-XRY2,-XRZ2,XAA2,                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -XRX1,-XRY1,-XRZ1,XAA1,
     4            B12(L),XB12,OV4,KLM)
C              
	END IF
C     -----------------------------------------------------
              Stemp=OV1+OV2+OV3+OV4
C
C	WRITE(6,*)'Stemp(',L,')=',Stemp
C     -----------------------------------------------------
C     --- Calculate KINETIC ENERGY OF ELECTRON 1 ---
C     -----------------------------------------------------
              CALL TK(DFLAG,OV1,
     1            RX1(L),RY1(L),RZ1(L),AA1(L),
     2            XRX1,XRY1,XRZ1,XAA1,                      
     3            RX2(L),RY2(L),RZ2(L),AA2(L),
     4            XRX2,XRY2,XRZ2,XAA2,
     5            B12(L),XB12,T1)  
C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
	IF (ISYMM) THEN
              CALL TK(DFLAG,OV2,
     1            RX1(L),RY1(L),RZ1(L),AA1(L),
     2            XRX2,XRY2,XRZ2,XAA2,                      
     3            RX2(L),RY2(L),RZ2(L),AA2(L),
     4            XRX1,XRY1,XRZ1,XAA1,
     5            B12(L),XB12,T2)  
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
              CALL TK(DFLAG,OV3,
     1            RX1(L),RY1(L),RZ1(L),AA1(L),
     2            -XRX1,-XRY1,-XRZ1,XAA1,                      
     3            RX2(L),RY2(L),RZ2(L),AA2(L),
     4            -XRX2,-XRY2,-XRZ2,XAA2,
     5            B12(L),XB12,T3)
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
              CALL TK(DFLAG,OV4,
     1            RX1(L),RY1(L),RZ1(L),AA1(L),
     2            -XRX2,-XRY2,-XRZ2,XAA2,                      
     3            RX2(L),RY2(L),RZ2(L),AA2(L),
     4            -XRX1,-XRY1,-XRZ1,XAA1,
     5            B12(L),XB12,T4)
	END IF
C     -----------------------------------------------------
C     --- Calculate KINETIC ENERGY OF ELECTRON 2 ---
C     -----------------------------------------------------
             CALL TK(DFLAG,OV1,
     1            RX2(L),RY2(L),RZ2(L),AA2(L),
     2            XRX2,XRY2,XRZ2,XAA2,                      
     3            RX1(L),RY1(L),RZ1(L),AA1(L),
     4            XRX1,XRY1,XRZ1,XAA1,
     5            B12(L),XB12,T5)  
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
	IF (ISYMM) THEN
              CALL TK(DFLAG,OV2,
     1            RX2(L),RY2(L),RZ2(L),AA2(L),
     2            XRX1,XRY1,XRZ1,XAA1,                      
     3            RX1(L),RY1(L),RZ1(L),AA1(L),
     4            XRX2,XRY2,XRZ2,XAA2,
     5            B12(L),XB12,T6) 
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
              CALL TK(DFLAG,OV3,
     1            RX2(L),RY2(L),RZ2(L),AA2(L),
     2            -XRX2,-XRY2,-XRZ2,XAA2,                      
     3            RX1(L),RY1(L),RZ1(L),AA1(L),
     4            -XRX1,-XRY1,-XRZ1,XAA1,
     5            B12(L),XB12,T7)
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
              CALL TK(DFLAG,OV4,
     1            RX2(L),RY2(L),RZ2(L),AA2(L),
     2            -XRX1,-XRY1,-XRZ1,XAA1,                      
     3            RX1(L),RY1(L),RZ1(L),AA1(L),
     4            -XRX2,-XRY2,-XRZ2,XAA2,
     5            B12(L),XB12,T8)
	END IF
C
	IF (ISYMM) THEN
            T=(T1+T2+T3+T4+T5+T6+T7+T8)/2.0D+00
	ELSE
		  T=(T1+T5)/2.0D+00
	END IF
C     -----------------------------------------------------
C     --- Calculate Nuc./Elec. Attractn. FOR ELECTRON 1 ---
C     -----------------------------------------------------
              CALL ENUCL(DFLAG,OV1,RX1(L),RY1(L),RZ1(L),AA1(L),
     1           XRX1,XRY1,XRZ1,XAA1,                      
     2           RX2(L),RY2(L),RZ2(L),AA2(L),
     3           XRX2,XRY2,XRZ2,XAA2,
     4           B12(L),XB12,
     5           EX1,EY1,EZ1,EN1,KLM)
C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C
	IF (ISYMM) THEN
              CALL ENUCL(DFLAG,OV2,RX1(L),RY1(L),RZ1(L),AA1(L),
     1           XRX2,XRY2,XRZ2,XAA2,                      
     2           RX2(L),RY2(L),RZ2(L),AA2(L),
     3           XRX1,XRY1,XRZ1,XAA1,
     4           B12(L),XB12,
     5           EX1,EY1,EZ1,EN2,KLM)                       
C
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
C
              CALL ENUCL(DFLAG,OV3,RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -XRX1,-XRY1,-XRZ1,XAA1,                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -XRX2,-XRY2,-XRZ2,XAA2,
     4            B12(L),XB12,
     5            EX1,EY1,EZ1,EN3,KLM)
C
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL ENUCL(DFLAG,OV4,RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -XRX2,-XRY2,-XRZ2,XAA2,                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -XRX1,-XRY1,-XRZ1,XAA1,
     4            B12(L),XB12,
     5            EX1,EY1,EZ1,EN4,KLM)
	END IF
C  
              CALL ENUCL(DFLAG,OV1,RX1(L),RY1(L),RZ1(L),AA1(L),
     1           XRX1,XRY1,XRZ1,XAA1,                      
     2           RX2(L),RY2(L),RZ2(L),AA2(L),
     3           XRX2,XRY2,XRZ2,XAA2,
     4           B12(L),XB12,
     5           EX2,EY2,EZ2,EN5,KLM)
C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C
	IF (ISYMM) THEN
              CALL ENUCL(DFLAG,OV2,RX1(L),RY1(L),RZ1(L),AA1(L),
     1           XRX2,XRY2,XRZ2,XAA2,                      
     2           RX2(L),RY2(L),RZ2(L),AA2(L),
     3           XRX1,XRY1,XRZ1,XAA1,
     4           B12(L),XB12,
     5           EX2,EY2,EZ2,EN6,KLM)                       
C
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
C
              CALL ENUCL(DFLAG,OV3,RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -XRX1,-XRY1,-XRZ1,XAA1,                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -XRX2,-XRY2,-XRZ2,XAA2,
     4            B12(L),XB12,
     5            EX2,EY2,EZ2,EN7,KLM)
C
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL ENUCL(DFLAG,OV4,RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -XRX2,-XRY2,-XRZ2,XAA2,                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -XRX1,-XRY1,-XRZ1,XAA1,
     4            B12(L),XB12,
     5            EX2,EY2,EZ2,EN8,KLM)
	END IF
C 
C     -----------------------------------------------------
C     --- Calculate Nuc./Elec. Attractn. FOR ELECTRON 2 ---
C     -----------------------------------------------------
              CALL ENUCL(DFLAG,OV1,RX2(L),RY2(L),RZ2(L),AA2(L),
     1           XRX2,XRY2,XRZ2,XAA2,                      
     2           RX1(L),RY1(L),RZ1(L),AA1(L),
     3           XRX1,XRY1,XRZ1,XAA1,
     4           B12(L),XB12,
     5           EX1,EY1,EZ1,EN9,KLM)                       
C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C
	IF (ISYMM) THEN
              CALL ENUCL(DFLAG,OV2,RX2(L),RY2(L),RZ2(L),AA2(L),
     1           XRX1,XRY1,XRZ1,XAA1,                      
     2           RX1(L),RY1(L),RZ1(L),AA1(L),
     3           XRX2,XRY2,XRZ2,XAA2,
     4           B12(L),XB12,
     5           EX1,EY1,EZ1,EN10,KLM)                       
C
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
C
              CALL ENUCL(DFLAG,OV3,RX2(L),RY2(L),RZ2(L),AA2(L),
     1            -XRX2,-XRY2,-XRZ2,XAA2,                      
     2            RX1(L),RY1(L),RZ1(L),AA1(L),
     3            -XRX1,-XRY1,-XRZ1,XAA1,
     4            B12(L),XB12,
     5            EX1,EY1,EZ1,EN11,KLM)
C
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL ENUCL(DFLAG,OV4,RX2(L),RY2(L),RZ2(L),AA2(L),
     1            -XRX1,-XRY1,-XRZ1,XAA1,                      
     2            RX1(L),RY1(L),RZ1(L),AA1(L),
     3            -XRX2,-XRY2,-XRZ2,XAA2,
     4            B12(L),XB12,
     5            EX1,EY1,EZ1,EN12,KLM)
	END IF
C 
              CALL ENUCL(DFLAG,OV1,RX2(L),RY2(L),RZ2(L),AA2(L),
     1           XRX2,XRY2,XRZ2,XAA2,                      
     2           RX1(L),RY1(L),RZ1(L),AA1(L),
     3           XRX1,XRY1,XRZ1,XAA1,
     4           B12(L),XB12,
     5           EX2,EY2,EZ2,EN13,KLM)                       
C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C
	IF (ISYMM) THEN
              CALL ENUCL(DFLAG,OV2,RX2(L),RY2(L),RZ2(L),AA2(L),
     1           XRX1,XRY1,XRZ1,XAA1,                      
     2           RX1(L),RY1(L),RZ1(L),AA1(L),
     3           XRX2,XRY2,XRZ2,XAA2,
     4           B12(L),XB12,
     5           EX2,EY2,EZ2,EN14,KLM)  
C
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
C
              CALL ENUCL(DFLAG,OV3,RX2(L),RY2(L),RZ2(L),AA2(L),
     1            -XRX2,-XRY2,-XRZ2,XAA2,                      
     2            RX1(L),RY1(L),RZ1(L),AA1(L),
     3            -XRX1,-XRY1,-XRZ1,XAA1,
     4            B12(L),XB12,
     5            EX2,EY2,EZ2,EN15,KLM)
C
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL ENUCL(DFLAG,OV4,RX2(L),RY2(L),RZ2(L),AA2(L),
     1            -XRX1,-XRY1,-XRZ1,XAA1,                      
     2            RX1(L),RY1(L),RZ1(L),AA1(L),
     3            -XRX2,-XRY2,-XRZ2,XAA2,
     4            B12(L),XB12,
     5            EX2,EY2,EZ2,EN16,KLM)
	END IF
C 
C
	IF (ISYMM) THEN
              ENN1=EN1+EN2+EN3+EN4+EN5+EN6+EN7+EN8
     1            +EN9+EN10+EN11+EN12+EN13+EN14+EN15+EN16
	ELSE
              ENN1=EN1+EN5+EN9+EN13
	END IF
C     -----------------------------------------------------
C     --- Calculate the Elec./Elec. Repulsion ---
C     -----------------------------------------------------
              CALL ELREP(DFLAG,OV1,RX1(L),RY1(L),RZ1(L),AA1(L),
     1           XRX1,XRY1,XRZ1,XAA1,                      
     2           RX2(L),RY2(L),RZ2(L),AA2(L),
     3           XRX2,XRY2,XRZ2,XAA2,
     4           B12(L),XB12,ER1,KLM)  
C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C
	IF (ISYMM) THEN
              CALL ELREP(DFLAG,OV2,RX1(L),RY1(L),RZ1(L),AA1(L),
     1           XRX2,XRY2,XRZ2,XAA2,                      
     2           RX2(L),RY2(L),RZ2(L),AA2(L),
     3           XRX1,XRY1,XRZ1,XAA1,
     4           B12(L),XB12,ER2,KLM)  
C
C     ...Nuclear Symmetrization Overlap Integral: P_[H1,H2]
C
              CALL ELREP(DFLAG,OV3,RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -XRX1,-XRY1,-XRZ1,XAA1,                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -XRX2,-XRY2,-XRZ2,XAA2,
     4            B12(L),XB12,ER3,KLM)
C
C     ...Electro/Nuclear Symmetrization Overlap Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL ELREP(DFLAG,OV4,RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -XRX2,-XRY2,-XRZ2,XAA2,                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -XRX1,-XRY1,-XRZ1,XAA1,
     4            B12(L),XB12,ER4,KLM)
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
C
      HH(1,2)=HH(1,2)+ EVEC(L)*(ER+T-EN+(1.0D+00/1.4011D+00)*Stemp)
	SS(1,2)=SS(1,2)+ (EVEC(L)*Stemp)
c	WRITE(6,*)'c',L,'*H(',L,',3)=',EVEC(L)*(ER+T-EN+(1.0D+00/
c	-1.4011D+00)*Stemp)
C     ----------------------------------------------------- 
      END DO
C	WRITE(6,*)'HH(1,2)=',HH(1,2)
C     -----------------------------------------------------
	Stemp=0.00D+00
	HH(1,1)=0.00D+00
      DO L=1,M
		DO K=1,M
			HH(1,1)=HH(1,1)+(EVEC(L)*EVEC(K)*H(L,K))
			Stemp=Stemp + (EVEC(L)*EVEC(K)*S(L,K))
		END DO
	END DO
	HH(1,1)=HH(1,1)/Stemp
C	WRITE(6,*)'Stemp=',Stemp
C	WRITE(6,*)'Snp1=',Snp1
C
	SS(2,2)=Snp1/Snp1
C	STOP
	HH(1,2)=HH(1,2)/DSQRT(Stemp*Snp1)
C	WRITE(6,*)'HH(1,2)=',HH(1,2)
	HH(2,1)=HH(1,2)
	SS(1,2)=SS(1,2)/DSQRT(Stemp*Snp1)
	SS(2,1)=SS(1,2)
	VH(1)=HH(1,1)
	VH(2)=HH(2,1)
	VH(3)=HH(2,2)				 
C     ----------------------------------------------------- 
C	CALL BIOUT(' ')
C      CALL BIOUT('BEFORE ROTATION...')
C	CALL BIOUT('HH MATRIX')
C	CALL TAB(HH,2,2,2,2)
C	CALL BIOUT(' ')
C	CALL BIOUT('SS MATRIX')
C	CALL TAB(SS,2,2,2,2)
C	CALL BIOUT(' ')
C	CALL YACOBI(VH,U,2,3,BIG,JB)
C	HH(1,1)=VH(1)
C	HH(1,2)=VH(2)
C	HH(2,1)=VH(2)
C	HH(2,2)=VH(3)
C	CALL BIOUT('AFTER ROTATION...')
C	CALL BIOUT('HH MATRIX')
C	CALL TAB(HH,2,2,2,2)
C	CALL BIOUT(' ')
C	CALL BIOUT('SS MATRIX')
C	CALL TAB(SS,2,2,2,2)
C	CALL BIOUT(' ')
C     ----------------------------------------------------- 
C      CALL TAB(S,M,M,MMAX,MMAX)                  
C      CALL TAB(H,M,M,MMAX,MMAX)
      CALL ENERGY(2,HH,SS,EMIN,CO,KIM)
C      J=1
C	CALL BIOUT(' ')
C      CALL BIOUT('     LinCoeff:.   AA1:        AA2:        B12:        
C	-RZ1:       RZ2:')
C      CALL BIOUT('FNC: =================================================
C     -=====================')
C	DO J=1,M
C          WRITE(6,1111) J,EVEC(J),AA1(J),AA2(J),B12(J),RZ1(J),RZ2(J)
C          WRITE(7,1111) J,EVEC(J),AA1(J),AA2(J),B12(J),RZ1(J),RZ2(J)
C	END DO
C          WRITE(6,1111) M+1,CO(2),XAA1,XAA2,XB12,XRZ1,XRZ2
C          WRITE(7,1111) M+1,CO(2),XAA1,XAA2,XB12,XRZ1,XRZ2
C		STOP
C     -----------------------------------------------------      
C	WRITE(6,*)'ENERGY DIAGNOSTIC INDEX FOR ADDITION OF THE FUNCTION=',
C	-EMIN-HH(1,1)
c      SUMH=0.00D+00
c	SUMS=0.00D+00
c	DO J=1,M
c		DO I=1,M
c			SUMH=SUMH+EVEC(I)*H(I,J)*EVEC(J)
c			SUMS=SUMS+EVEC(I)*S(I,J)*EVEC(J)
c		END DO
c	END DO
cC
c	ENG=SUMH/SUMS
c	EMIN=ENG
cC
C     -----------------------------------------------------      
1111  FORMAT(I3,6D12.4)
      RETURN
      END