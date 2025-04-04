C     L. ADAMOWICZ RESEARCH / UNIVERSITY OF ARIZONA / CHEM.
C     -----------------------------------------------------
      SUBROUTINE HAMILTON(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,
     -                    EVEC,S,H,CO,EMIN,GRAD,TTT1,TTT2,TTT3,DFLAG)
C     -----------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER M,ERROR,DFLAG
      REAL*8, DIMENSION(M) ::	AA1,RX1,RY1,RZ1,
	-						AA2,RX2,RY2,RZ2,
     -						B12,EVEC,CO,TTT1,TTT2,TTT3
	REAL*8, DIMENSION(M*6) :: GRAD
	REAL*8, DIMENSION(M,M) :: S,H
C
	REAL*8  TEMP[ALLOCATABLE](:)

	REAL*8	DHDA1[ALLOCATABLE](:,:),DHDA2[ALLOCATABLE](:,:),
	-		DHDA5[ALLOCATABLE](:,:),DHDAZ[ALLOCATABLE](:,:),
     -		DHDBZ[ALLOCATABLE](:,:)
	REAL*8	DSDA1[ALLOCATABLE](:,:),DSDA2[ALLOCATABLE](:,:),
	-		DSDA5[ALLOCATABLE](:,:),DSDAZ[ALLOCATABLE](:,:),
     -		DSDBZ[ALLOCATABLE](:,:)
	REAL*8	GRADM[ALLOCATABLE](:,:),DS[ALLOCATABLE](:,:),
	-		DH[ALLOCATABLE](:,:),GRADA[ALLOCATABLE](:),
     -		GRADC[ALLOCATABLE](:)
C     -----------------------------------------------------
	COMMON /OVER/ L,K
	COMMON /TNOVER/ NXAC,XAC
	COMMON /DOV/ DOVDA1,DOVDA2,DOVDA5,DOVDAZ,DOVDBZ
	COMMON /DOVn/ DOVnDA1,DOVnDA2,DOVnDA5,DOVnDAZ,DOVnDBZ
	COMMON /DTK/ DTKDA1,DTKDA2,DTKDA5,DTKDAZ,DTKDBZ
	COMMON /DEN/ DENDA1,DENDA2,DENDA5,DENDAZ,DENDBZ,
	1			 DENDA1P,DENDA2P,DENDA5P,DENDAZP,DENDBZP
	COMMON /DER/ DERDA1,DERDA2,DERDA5,DERDAZ,DERDBZ
C     -----------------------------------------------------
	PARAMETER(ISYMM=1)
	PARAMETER(IPERM=0)
C     -----------------------------------------------------
      MMAX=M
      ALLOCATE(	TEMP(MMAX),DHDA1(MMAX,MMAX),DHDA2(MMAX,MMAX),
	-			DHDA5(MMAX,MMAX),DHDAZ(MMAX,MMAX),DHDBZ(MMAX,MMAX),
     -			DSDA1(MMAX,MMAX),DSDA2(MMAX,MMAX),DSDA5(MMAX,MMAX),
     -			DSDAZ(MMAX,MMAX),DSDBZ(MMAX,MMAX),GRADM(MMAX,MMAX*5),
     -			DS(MMAX*5,MMAX),DH(MMAX*5,MMAX),GRADA(MMAX*5),
     -			GRADC(MMAX),STAT=ERROR)
C
	IF (ERROR .NE. 0) STOP '!NOT ENOUGH MEMORY SPACE AVAILABLE ON THIS
	- BOX! ...ABORTING!'
C     -----------------------------------------------------
      DATA ONE,TWO,FOUR,SIX,EIGHT,HALF/1.,2.,4.,6.,8.,.5/
C     -----------------------------------------------------
      CALL GEOM(EX1,EY1,EZ1,EX2,EY2,EZ2)
C     -----------------------------------------------------	  
	IF (DFLAG .EQ. 1) THEN
      DO L=1,M
          DO K=1,M
C     -----------------------------------------------------
C     --- Calculate OVERLAP INTEGRAL ---
C     -----------------------------------------------------
             CALL OVERLP(M,DFLAG,RX1(L),RY1(L),RZ1(L),AA1(L),
     1            RX1(K),RY1(K),RZ1(K),AA1(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            RX2(K),RY2(K),RZ2(K),AA2(K),
     4            B12(L),B12(K),OV1,KLM)
C
	DOV1DA1=DOVDA1
C	WRITE(6,*)'DSDA1(1)=',DOV1DA1
C	STOP
	DOV1DA2=DOVDA2
	DOV1DA5=DOVDA5
	DOV1DAZ=DOVDAZ
	DOV1DBZ=DOVDBZ
C
	IF(ISYMM) THEN
		CALL OVERLP(M,DFLAG,RX1(L),RY1(L),RZ1(L),AA1(L),
     1            RX2(K),RY2(K),RZ2(K),AA2(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            RX1(K),RY1(K),RZ1(K),AA1(K),
     4            B12(L),B12(K),OV2,KLM)
C
	DOV2DA1=DOVDA1
C	WRITE(6,*)'DSDA1(2)=',DOV2DA1
C	STOP
	DOV2DA2=DOVDA2
	DOV2DA5=DOVDA5
	DOV2DAZ=DOVDAZ
	DOV2DBZ=DOVDBZ
C
C     ...Nuclear Symmetrization Overlap Integral: P_[H1,H2]
C
              CALL OVERLP(M,DFLAG,RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
     4            B12(L),B12(K),OV3,KLM)
C
	DOV3DA1=DOVDA1
C	WRITE(6,*)'DSDA1(3)=',DOV3DA1
C	STOP
	DOV3DA2=DOVDA2
	DOV3DA5=DOVDA5
	DOV3DAZ=DOVDAZ
	DOV3DBZ=DOVDBZ
C
C     ...Electro/Nuclear Symmetrization Overlap Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL OVERLP(M,DFLAG,RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
     4            B12(L),B12(K),OV4,KLM)
C              
	DOV4DA1=DOVDA1
C	WRITE(6,*)'DSDA1(4)=',DOV4DA1
C	STOP
	DOV4DA2=DOVDA2
	DOV4DA5=DOVDA5
	DOV4DAZ=DOVDAZ
	DOV4DBZ=DOVDBZ
	END IF
C     -----------------------------------------------------
	IF(ISYMM) THEN
              S(L,K)=OV1+OV2+OV3+OV4
		ELSE
              S(L,K)=OV1
	END IF
	IF(IPERM.EQ.1) THEN
		S(L,K)=OV1
	END IF
	IF(IPERM.EQ.2) THEN
		S(L,K)=OV2
	END IF
	IF(IPERM.EQ.3) THEN
		S(L,K)=OV3
	END IF
	IF(IPERM.EQ.4) THEN
		S(L,K)=OV4
	END IF
C
	IF(ISYMM) THEN
	DSDA1(L,K)=DOV1DA1+DOV2DA1+DOV3DA1+DOV4DA1
	DSDA2(L,K)=DOV1DA2+DOV2DA2+DOV3DA2+DOV4DA2
	DSDA5(L,K)=DOV1DA5+DOV2DA5+DOV3DA5+DOV4DA5
	DSDAZ(L,K)=DOV1DAZ+DOV2DAZ+DOV3DAZ+DOV4DAZ
	DSDBZ(L,K)=DOV1DBZ+DOV2DBZ+DOV3DBZ+DOV4DBZ
	ELSE
	DSDA1(L,K)=DOV1DA1
	DSDA2(L,K)=DOV1DA2
	DSDA5(L,K)=DOV1DA5
	DSDAZ(L,K)=DOV1DAZ
	DSDBZ(L,K)=DOV1DBZ
	END IF		
	IF(IPERM.EQ.1) THEN
	DSDA1(L,K)=DOV1DA1
	DSDA2(L,K)=DOV1DA2
	DSDA5(L,K)=DOV1DA5
	DSDAZ(L,K)=DOV1DAZ
	DSDBZ(L,K)=DOV1DBZ
	END IF
	IF(IPERM.EQ.2) THEN
	DSDA1(L,K)=DOV2DA1
	DSDA2(L,K)=DOV2DA2
	DSDA5(L,K)=DOV2DA5
	DSDAZ(L,K)=DOV2DAZ
	DSDBZ(L,K)=DOV2DBZ
	END IF
	IF(IPERM.EQ.3) THEN
	DSDA1(L,K)=DOV3DA1
	DSDA2(L,K)=DOV3DA2
	DSDA5(L,K)=DOV3DA5
	DSDAZ(L,K)=DOV3DAZ
	DSDBZ(L,K)=DOV3DBZ
	END IF
	IF(IPERM.EQ.4) THEN
	DSDA1(L,K)=DOV4DA1
	DSDA2(L,K)=DOV4DA2
	DSDA5(L,K)=DOV4DA5
	DSDAZ(L,K)=DOV4DAZ
	DSDBZ(L,K)=DOV4DBZ
	END IF
C     -----------------------------------------------------
C     --- Calculate KINETIC ENERGY OF ELECTRON 1 ---
C     -----------------------------------------------------
	DOVnDA1=DOV1DA1
	DOVnDA2=DOV1DA2
	DOVnDA5=DOV1DA5
	DOVnDAZ=DOV1DAZ
	DOVnDBZ=DOV1DBZ
              CALL TK(DFLAG,OV1,
     1            RX1(L),RY1(L),RZ1(L),AA1(L),
     2            RX1(K),RY1(K),RZ1(K),AA1(K),                      
     3            RX2(L),RY2(L),RZ2(L),AA2(L),
     4            RX2(K),RY2(K),RZ2(K),AA2(K),
     5            B12(L),B12(K),T1)  
	D1T1DA1=DTKDA1
	D1T1DA2=DTKDA2
	D1T1DA5=DTKDA5
	D1T1DAZ=DTKDAZ
	D1T1DBZ=DTKDBZ
C





C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
	IF (ISYMM) THEN
	DOVnDA1=DOV2DA1
	DOVnDA2=DOV2DA2
	DOVnDA5=DOV2DA5
	DOVnDAZ=DOV2DAZ
	DOVnDBZ=DOV2DBZ
              CALL TK(DFLAG,OV2,
     1            RX1(L),RY1(L),RZ1(L),AA1(L),
     2            RX2(K),RY2(K),RZ2(K),AA2(K),                      
     3            RX2(L),RY2(L),RZ2(L),AA2(L),
     4            RX1(K),RY1(K),RZ1(K),AA1(K),
     5            B12(L),B12(K),T2)  
	D2T1DA1=DTKDA1
	D2T1DA2=DTKDA2
	D2T1DA5=DTKDA5
	D2T1DAZ=DTKDAZ
	D2T1DBZ=DTKDBZ



C     ...Nuclear Symmetrization Integral: P_[H1,H2]
	DOVnDA1=DOV3DA1
	DOVnDA2=DOV3DA2
	DOVnDA5=DOV3DA5
	DOVnDAZ=DOV3DAZ
	DOVnDBZ=DOV3DBZ
              CALL TK(DFLAG,OV3,
     1            RX1(L),RY1(L),RZ1(L),AA1(L),
     2            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
     3            RX2(L),RY2(L),RZ2(L),AA2(L),
     4            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
     5            B12(L),B12(K),T3)
	D3T1DA1=DTKDA1
	D3T1DA2=DTKDA2
	D3T1DA5=DTKDA5
	D3T1DAZ=DTKDAZ
	D3T1DBZ=DTKDBZ



C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
	DOVnDA1=DOV4DA1
	DOVnDA2=DOV4DA2
	DOVnDA5=DOV4DA5
	DOVnDAZ=DOV4DAZ
	DOVnDBZ=DOV4DBZ
              CALL TK(DFLAG,OV4,
     1            RX1(L),RY1(L),RZ1(L),AA1(L),
     2            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
     3            RX2(L),RY2(L),RZ2(L),AA2(L),
     4            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
     5            B12(L),B12(K),T4)
	D4T1DA1=DTKDA1
	D4T1DA2=DTKDA2
	D4T1DA5=DTKDA5
	D4T1DAZ=DTKDAZ
	D4T1DBZ=DTKDBZ
	END IF



C     -----------------------------------------------------
C     --- Calculate KINETIC ENERGY OF ELECTRON 2 ---
C     -----------------------------------------------------
	DOVnDA1=DOV1DA2
	DOVnDA2=DOV1DA1
	DOVnDA5=DOV1DA5
	DOVnDAZ=DOV1DBZ
	DOVnDBZ=DOV1DAZ
             CALL TK(DFLAG,OV1,
     1            RX2(L),RY2(L),RZ2(L),AA2(L),
     2            RX2(K),RY2(K),RZ2(K),AA2(K),                      
     3            RX1(L),RY1(L),RZ1(L),AA1(L),
     4            RX1(K),RY1(K),RZ1(K),AA1(K),
     5            B12(L),B12(K),T5)  
	D1T2DA1=DTKDA2
	D1T2DA2=DTKDA1
	D1T2DA5=DTKDA5
	D1T2DAZ=DTKDBZ
	D1T2DBZ=DTKDAZ




C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
	IF (ISYMM) THEN
	DOVnDA1=DOV2DA2
	DOVnDA2=DOV2DA1
	DOVnDA5=DOV2DA5
	DOVnDAZ=DOV2DBZ
	DOVnDBZ=DOV2DAZ
              CALL TK(DFLAG,OV2,
     1            RX2(L),RY2(L),RZ2(L),AA2(L),
     2            RX1(K),RY1(K),RZ1(K),AA1(K),                      
     3            RX1(L),RY1(L),RZ1(L),AA1(L),
     4            RX2(K),RY2(K),RZ2(K),AA2(K),
     5            B12(L),B12(K),T6) 
	D2T2DA1=DTKDA2
	D2T2DA2=DTKDA1
	D2T2DA5=DTKDA5
	D2T2DAZ=DTKDBZ
	D2T2DBZ=DTKDAZ



C     ...Nuclear Symmetrization Integral: P_[H1,H2]
	DOVnDA1=DOV3DA2
	DOVnDA2=DOV3DA1
	DOVnDA5=DOV3DA5
	DOVnDAZ=DOV3DBZ
	DOVnDBZ=DOV3DAZ
              CALL TK(DFLAG,OV3,
     1            RX2(L),RY2(L),RZ2(L),AA2(L),
     2            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
     3            RX1(L),RY1(L),RZ1(L),AA1(L),
     4            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
     5            B12(L),B12(K),T7)
	D3T2DA1=DTKDA2
	D3T2DA2=DTKDA1
	D3T2DA5=DTKDA5
	D3T2DAZ=DTKDBZ
	D3T2DBZ=DTKDAZ



C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
	DOVnDA1=DOV4DA2
	DOVnDA2=DOV4DA1
	DOVnDA5=DOV4DA5
	DOVnDAZ=DOV4DBZ
	DOVnDBZ=DOV4DAZ
              CALL TK(DFLAG,OV4,
     1            RX2(L),RY2(L),RZ2(L),AA2(L),
     2            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
     3            RX1(L),RY1(L),RZ1(L),AA1(L),
     4            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
     5            B12(L),B12(K),T8)
	D4T2DA1=DTKDA2
	D4T2DA2=DTKDA1
	D4T2DA5=DTKDA5
	D4T2DAZ=DTKDBZ
	D4T2DBZ=DTKDAZ
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
C
	IF (ISYMM) THEN
	DTDA1=(D1T1DA1+D2T1DA1+D3T1DA1+D4T1DA1+D1T2DA1+D2T2DA1+D3T2DA1+D4T
	-2DA1)/2.00D+00
	DTDA2=(D1T1DA2+D2T1DA2+D3T1DA2+D4T1DA2+D1T2DA2+D2T2DA2+D3T2DA2+D4T
	-2DA2)/2.00D+00
	DTDA5=(D1T1DA5+D2T1DA5+D3T1DA5+D4T1DA5+D1T2DA5+D2T2DA5+D3T2DA5+D4T
	-2DA5)/2.00D+00
	DTDAZ=(D1T1DAZ+D2T1DAZ+D3T1DAZ+D4T1DAZ+D1T2DAZ+D2T2DAZ+D3T2DAZ+D4T
	-2DAZ)/2.00D+00
	DTDBZ=(D1T1DBZ+D2T1DBZ+D3T1DBZ+D4T1DBZ+D1T2DBZ+D2T2DBZ+D3T2DBZ+D4T
	-2DBZ)/2.00D+00
	ELSE
	DTDA1=(D1T1DA1+D1T2DA1)/2.00D+00
	DTDA2=(D1T1DA2+D1T2DA2)/2.00D+00
	DTDA5=(D1T1DA5+D1T2DA5)/2.00D+00
	DTDAZ=(D1T1DAZ+D1T2DAZ)/2.00D+00
	DTDBZ=(D1T1DBZ+D1T2DBZ)/2.00D+00
	END IF
	IF(IPERM.EQ.1) THEN
	DTDA1=(D1T1DA1+D1T2DA1)/2.00D+00
	DTDA2=(D1T1DA2+D1T2DA2)/2.00D+00
	DTDA5=(D1T1DA5+D1T2DA5)/2.00D+00
	DTDAZ=(D1T1DAZ+D1T2DAZ)/2.00D+00
	DTDBZ=(D1T1DBZ+D1T2DBZ)/2.00D+00
	END IF
	IF(IPERM.EQ.2) THEN
	DTDA1=(D2T1DA1+D2T2DA1)/2.00D+00
	DTDA2=(D2T1DA2+D2T2DA2)/2.00D+00
	DTDA5=(D2T1DA5+D2T2DA5)/2.00D+00
	DTDAZ=(D2T1DAZ+D2T2DAZ)/2.00D+00
	DTDBZ=(D2T1DBZ+D2T2DBZ)/2.00D+00
	END IF
	IF(IPERM.EQ.3) THEN
	DTDA1=(D3T1DA1+D3T2DA1)/2.00D+00
	DTDA2=(D3T1DA2+D3T2DA2)/2.00D+00
	DTDA5=(D1T1DA5+D3T2DA5)/2.00D+00
	DTDAZ=(D3T1DAZ+D3T2DAZ)/2.00D+00
	DTDBZ=(D3T1DBZ+D3T2DBZ)/2.00D+00
	END IF
	IF(IPERM.EQ.4) THEN
	DTDA1=(D4T1DA1+D4T2DA1)/2.00D+00
	DTDA2=(D4T1DA2+D4T2DA2)/2.00D+00
	DTDA5=(D4T1DA5+D4T2DA5)/2.00D+00
	DTDAZ=(D4T1DAZ+D4T2DAZ)/2.00D+00
	DTDBZ=(D4T1DBZ+D4T2DBZ)/2.00D+00
	END IF
C     -----------------------------------------------------
C     --- Calculate Nuc./Elec. Attractn. FOR ELECTRON 1 ---
C     -----------------------------------------------------
	DOVnDA1=DOV1DA1
	DOVnDA2=DOV1DA2
	DOVnDA5=DOV1DA5
	DOVnDAZ=DOV1DAZ
	DOVnDBZ=DOV1DBZ
              CALL ENUCL(DFLAG,OV1,RX1(L),RY1(L),RZ1(L),AA1(L),
     1           RX1(K),RY1(K),RZ1(K),AA1(K),                      
     2           RX2(L),RY2(L),RZ2(L),AA2(L),
     3           RX2(K),RY2(K),RZ2(K),AA2(K),
     4           B12(L),B12(K),
     5           EX1,EY1,EZ1,EN1,KLM)
	D1EN1DA1=DENDA1
	D1EN1DA2=DENDA2
	D1EN1DA5=DENDA5
	D1EN1DAZ=DENDAZ
	D1EN1DBZ=DENDBZ


C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C
	IF (ISYMM) THEN
	DOVnDA1=DOV2DA1
	DOVnDA2=DOV2DA2
	DOVnDA5=DOV2DA5
	DOVnDAZ=DOV2DAZ
	DOVnDBZ=DOV2DBZ
              CALL ENUCL(DFLAG,OV2,RX1(L),RY1(L),RZ1(L),AA1(L),
     1           RX2(K),RY2(K),RZ2(K),AA2(K),                      
     2           RX2(L),RY2(L),RZ2(L),AA2(L),
     3           RX1(K),RY1(K),RZ1(K),AA1(K),
     4           B12(L),B12(K),
     5           EX1,EY1,EZ1,EN2,KLM)                       
	D1EN2DA1=DENDA1
	D1EN2DA2=DENDA2
	D1EN2DA5=DENDA5
	D1EN2DAZ=DENDAZ
	D1EN2DBZ=DENDBZ


C
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
C
	DOVnDA1=DOV3DA1
	DOVnDA2=DOV3DA2
	DOVnDA5=DOV3DA5
	DOVnDAZ=DOV3DAZ
	DOVnDBZ=DOV3DBZ
              CALL ENUCL(DFLAG,OV3,RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
     4            B12(L),B12(K),
     5            EX1,EY1,EZ1,EN3,KLM)
	D1EN3DA1=DENDA1
	D1EN3DA2=DENDA2
	D1EN3DA5=DENDA5
	D1EN3DAZ=DENDAZ
	D1EN3DBZ=DENDBZ



C
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
	DOVnDA1=DOV4DA1
	DOVnDA2=DOV4DA2
	DOVnDA5=DOV4DA5
	DOVnDAZ=DOV4DAZ
	DOVnDBZ=DOV4DBZ
              CALL ENUCL(DFLAG,OV4,RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
     4            B12(L),B12(K),
     5            EX1,EY1,EZ1,EN4,KLM)
	D1EN4DA1=DENDA1
	D1EN4DA2=DENDA2
	D1EN4DA5=DENDA5
	D1EN4DAZ=DENDAZ
	D1EN4DBZ=DENDBZ



	END IF
C  

	DOVnDA1=DOV1DA1
	DOVnDA2=DOV1DA2
	DOVnDA5=DOV1DA5
	DOVnDAZ=DOV1DAZ
	DOVnDBZ=DOV1DBZ
              CALL ENUCL(DFLAG,OV1,RX1(L),RY1(L),RZ1(L),AA1(L),
     1           RX1(K),RY1(K),RZ1(K),AA1(K),                      
     2           RX2(L),RY2(L),RZ2(L),AA2(L),
     3           RX2(K),RY2(K),RZ2(K),AA2(K),
     4           B12(L),B12(K),
     5           EX2,EY2,EZ2,EN5,KLM)
	D1EN5DA1=DENDA1P
	D1EN5DA2=DENDA2P
	D1EN5DA5=DENDA5P
	D1EN5DAZ=DENDAZP
	D1EN5DBZ=DENDBZP

C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C
	IF (ISYMM) THEN
	DOVnDA1=DOV2DA1
	DOVnDA2=DOV2DA2
	DOVnDA5=DOV2DA5
	DOVnDAZ=DOV2DAZ
	DOVnDBZ=DOV2DBZ
              CALL ENUCL(DFLAG,OV2,RX1(L),RY1(L),RZ1(L),AA1(L),
     1           RX2(K),RY2(K),RZ2(K),AA2(K),                      
     2           RX2(L),RY2(L),RZ2(L),AA2(L),
     3           RX1(K),RY1(K),RZ1(K),AA1(K),
     4           B12(L),B12(K),
     5           EX2,EY2,EZ2,EN6,KLM)                       
	D1EN6DA1=DENDA1P
	D1EN6DA2=DENDA2P
	D1EN6DA5=DENDA5P
	D1EN6DAZ=DENDAZP
	D1EN6DBZ=DENDBZP


C
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
C
	DOVnDA1=DOV3DA1
	DOVnDA2=DOV3DA2
	DOVnDA5=DOV3DA5
	DOVnDAZ=DOV3DAZ
	DOVnDBZ=DOV3DBZ
              CALL ENUCL(DFLAG,OV3,RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
     4            B12(L),B12(K),
     5            EX2,EY2,EZ2,EN7,KLM)
	D1EN7DA1=DENDA1P
	D1EN7DA2=DENDA2P
	D1EN7DA5=DENDA5P
	D1EN7DAZ=DENDAZP
	D1EN7DBZ=DENDBZP


C
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
	DOVnDA1=DOV4DA1
	DOVnDA2=DOV4DA2
	DOVnDA5=DOV4DA5
	DOVnDAZ=DOV4DAZ
	DOVnDBZ=DOV4DBZ
              CALL ENUCL(DFLAG,OV4,RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
     4            B12(L),B12(K),
     5            EX2,EY2,EZ2,EN8,KLM)
	D1EN8DA1=DENDA1P
	D1EN8DA2=DENDA2P
	D1EN8DA5=DENDA5P
	D1EN8DAZ=DENDAZP
	D1EN8DBZ=DENDBZP



	END IF
C 
C     -----------------------------------------------------
C     --- Calculate Nuc./Elec. Attractn. FOR ELECTRON 2 ---
C     -----------------------------------------------------
	DOVnDA1=DOV1DA2
	DOVnDA2=DOV1DA1
	DOVnDA5=DOV1DA5
	DOVnDAZ=DOV1DBZ
	DOVnDBZ=DOV1DAZ
              CALL ENUCL(DFLAG,OV1,RX2(L),RY2(L),RZ2(L),AA2(L),
     1           RX2(K),RY2(K),RZ2(K),AA2(K),                      
     2           RX1(L),RY1(L),RZ1(L),AA1(L),
     3           RX1(K),RY1(K),RZ1(K),AA1(K),
     4           B12(L),B12(K),
     5           EX1,EY1,EZ1,EN9,KLM)                       
	D2EN1DA1=DENDA2
	D2EN1DA2=DENDA1
	D2EN1DA5=DENDA5
	D2EN1DAZ=DENDBZ
	D2EN1DBZ=DENDAZ



C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C
	IF (ISYMM) THEN
	DOVnDA1=DOV2DA2
	DOVnDA2=DOV2DA1
	DOVnDA5=DOV2DA5
	DOVnDAZ=DOV2DBZ
	DOVnDBZ=DOV2DAZ
              CALL ENUCL(DFLAG,OV2,RX2(L),RY2(L),RZ2(L),AA2(L),
     1           RX1(K),RY1(K),RZ1(K),AA1(K),                      
     2           RX1(L),RY1(L),RZ1(L),AA1(L),
     3           RX2(K),RY2(K),RZ2(K),AA2(K),
     4           B12(L),B12(K),
     5           EX1,EY1,EZ1,EN10,KLM)                       
	D2EN2DA1=DENDA2
	D2EN2DA2=DENDA1
	D2EN2DA5=DENDA5
	D2EN2DAZ=DENDBZ
	D2EN2DBZ=DENDAZ

C
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
C
	DOVnDA1=DOV3DA2
	DOVnDA2=DOV3DA1
	DOVnDA5=DOV3DA5
	DOVnDAZ=DOV3DBZ
	DOVnDBZ=DOV3DAZ
              CALL ENUCL(DFLAG,OV3,RX2(L),RY2(L),RZ2(L),AA2(L),
     1            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
     2            RX1(L),RY1(L),RZ1(L),AA1(L),
     3            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
     4            B12(L),B12(K),
     5            EX1,EY1,EZ1,EN11,KLM)
	D2EN3DA1=DENDA2
	D2EN3DA2=DENDA1
	D2EN3DA5=DENDA5
	D2EN3DAZ=DENDBZ
	D2EN3DBZ=DENDAZ


C
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
	DOVnDA1=DOV4DA2
	DOVnDA2=DOV4DA1
	DOVnDA5=DOV4DA5
	DOVnDAZ=DOV4DBZ
	DOVnDBZ=DOV4DAZ
              CALL ENUCL(DFLAG,OV4,RX2(L),RY2(L),RZ2(L),AA2(L),
     1            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
     2            RX1(L),RY1(L),RZ1(L),AA1(L),
     3            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
     4            B12(L),B12(K),
     5            EX1,EY1,EZ1,EN12,KLM)
	D2EN4DA1=DENDA2
	D2EN4DA2=DENDA1
	D2EN4DA5=DENDA5
	D2EN4DAZ=DENDBZ
	D2EN4DBZ=DENDAZ


	END IF
C 

 
	DOVnDA1=DOV1DA2
	DOVnDA2=DOV1DA1
	DOVnDA5=DOV1DA5
	DOVnDAZ=DOV1DBZ
	DOVnDBZ=DOV1DAZ
              CALL ENUCL(DFLAG,OV1,RX2(L),RY2(L),RZ2(L),AA2(L),
     1           RX2(K),RY2(K),RZ2(K),AA2(K),                      
     2           RX1(L),RY1(L),RZ1(L),AA1(L),
     3           RX1(K),RY1(K),RZ1(K),AA1(K),
     4           B12(L),B12(K),
     5           EX2,EY2,EZ2,EN13,KLM)                       
	D2EN5DA1=DENDA2P
	D2EN5DA2=DENDA1P
	D2EN5DA5=DENDA5P
	D2EN5DAZ=DENDBZP
	D2EN5DBZ=DENDAZP


C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C
	IF (ISYMM) THEN
	DOVnDA1=DOV2DA2
	DOVnDA2=DOV2DA1
	DOVnDA5=DOV2DA5
	DOVnDAZ=DOV2DBZ
	DOVnDBZ=DOV2DAZ
              CALL ENUCL(DFLAG,OV2,RX2(L),RY2(L),RZ2(L),AA2(L),
     1           RX1(K),RY1(K),RZ1(K),AA1(K),                      
     2           RX1(L),RY1(L),RZ1(L),AA1(L),
     3           RX2(K),RY2(K),RZ2(K),AA2(K),
     4           B12(L),B12(K),
     5           EX2,EY2,EZ2,EN14,KLM)  
	D2EN6DA1=DENDA2P
	D2EN6DA2=DENDA1P
	D2EN6DA5=DENDA5P
	D2EN6DAZ=DENDBZP
	D2EN6DBZ=DENDAZP


C
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
C
	DOVnDA1=DOV3DA2
	DOVnDA2=DOV3DA1
	DOVnDA5=DOV3DA5
	DOVnDAZ=DOV3DBZ
	DOVnDBZ=DOV3DAZ
              CALL ENUCL(DFLAG,OV3,RX2(L),RY2(L),RZ2(L),AA2(L),
     1            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
     2            RX1(L),RY1(L),RZ1(L),AA1(L),
     3            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
     4            B12(L),B12(K),
     5            EX2,EY2,EZ2,EN15,KLM)
	D2EN7DA1=DENDA2P
	D2EN7DA2=DENDA1P
	D2EN7DA5=DENDA5P
	D2EN7DAZ=DENDBZP
	D2EN7DBZ=DENDAZP


C
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
	DOVnDA1=DOV4DA2
	DOVnDA2=DOV4DA1
	DOVnDA5=DOV4DA5
	DOVnDAZ=DOV4DBZ
	DOVnDBZ=DOV4DAZ
              CALL ENUCL(DFLAG,OV4,RX2(L),RY2(L),RZ2(L),AA2(L),
     1            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
     2            RX1(L),RY1(L),RZ1(L),AA1(L),
     3            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
     4            B12(L),B12(K),
     5            EX2,EY2,EZ2,EN16,KLM)
	D2EN8DA1=DENDA2P
	D2EN8DA2=DENDA1P
	D2EN8DA5=DENDA5P
	D2EN8DAZ=DENDBZP
	D2EN8DBZ=DENDAZP



	END IF
C 
C
	IF (ISYMM) THEN
              ENN1=EN1+EN2+EN3+EN4+EN5+EN6+EN7+EN8
     1            +EN9+EN10+EN11+EN12+EN13+EN14+EN15+EN16
	ELSE
              ENN1=EN1+EN5+EN9+EN13
	END IF
C
	IF (ISYMM) THEN
	DENADA1=D1EN1DA1+D1EN2DA1+D1EN3DA1+D1EN4DA1+D1EN5DA1+D1EN6DA1+D1EN
	-7DA1+D1EN8DA1+D2EN1DA1+D2EN2DA1+D2EN3DA1+D2EN4DA1+D2EN5DA1+D2EN6DA
     -1+D2EN7DA1+D2EN8DA1
	DENADA2=D1EN1DA2+D1EN2DA2+D1EN3DA2+D1EN4DA2+D1EN5DA2+D1EN6DA2+D1EN
	-7DA2+D1EN8DA2+D2EN1DA2+D2EN2DA2+D2EN3DA2+D2EN4DA2+D2EN5DA2+D2EN6DA
     -2+D2EN7DA2+D2EN8DA2
	DENADA5=D1EN1DA5+D1EN2DA5+D1EN3DA5+D1EN4DA5+D1EN5DA5+D1EN6DA5+D1EN
	-7DA5+D1EN8DA5+D2EN1DA5+D2EN2DA5+D2EN3DA5+D2EN4DA5+D2EN5DA5+D2EN6DA
     -5+D2EN7DA5+D2EN8DA5
	DENADAZ=D1EN1DAZ+D1EN2DAZ+D1EN3DAZ+D1EN4DAZ+D1EN5DAZ+D1EN6DAZ+D1EN
	-7DAZ+D1EN8DAZ+D2EN1DAZ+D2EN2DAZ+D2EN3DAZ+D2EN4DAZ+D2EN5DAZ+D2EN6DA
     -Z+D2EN7DAZ+D2EN8DAZ
	DENADBZ=D1EN1DBZ+D1EN2DBZ+D1EN3DBZ+D1EN4DBZ+D1EN5DBZ+D1EN6DBZ+D1EN
	-7DBZ+D1EN8DBZ+D2EN1DBZ+D2EN2DBZ+D2EN3DBZ+D2EN4DBZ+D2EN5DBZ+D2EN6DB
     -Z+D2EN7DBZ+D2EN8DBZ
	ELSE
	DENADA1=D1EN1DA1+D1EN5DA1+D2EN1DA1+D2EN5DA1
	DENADA2=D1EN1DA2+D1EN5DA2+D2EN1DA2+D2EN5DA2
	DENADA5=D1EN1DA5+D1EN5DA5+D2EN1DA5+D2EN5DA5
	DENADAZ=D1EN1DAZ+D1EN5DAZ+D2EN1DAZ+D2EN5DAZ
	DENADBZ=D1EN1DBZ+D1EN5DBZ+D2EN1DBZ+D2EN5DBZ
	END IF
	IF(IPERM.EQ.1) THEN
	DENADA1=D1EN1DA1+D1EN5DA1+D2EN1DA1+D2EN5DA1
	DENADA2=D1EN1DA2+D1EN5DA2+D2EN1DA2+D2EN5DA2
	DENADA5=D1EN1DA5+D1EN5DA5+D2EN1DA5+D2EN5DA5
	DENADAZ=D1EN1DAZ+D1EN5DAZ+D2EN1DAZ+D2EN5DAZ
	DENADBZ=D1EN1DBZ+D1EN5DBZ+D2EN1DBZ+D2EN5DBZ
	END IF
	IF(IPERM.EQ.2) THEN
	DENADA1=D1EN2DA1+D1EN6DA1+D2EN2DA1+D2EN6DA1
	DENADA2=D1EN2DA2+D1EN6DA2+D2EN2DA2+D2EN6DA2
	DENADA5=D1EN2DA5+D1EN6DA5+D2EN2DA5+D2EN6DA5
	DENADAZ=D1EN2DAZ+D1EN6DAZ+D2EN2DAZ+D2EN6DAZ
	DENADBZ=D1EN2DBZ+D1EN6DBZ+D2EN2DBZ+D2EN6DBZ
	END IF
	IF(IPERM.EQ.3) THEN
	DENADA1=D1EN3DA1+D1EN7DA1+D2EN3DA1+D2EN7DA1
	DENADA2=D1EN3DA2+D1EN7DA2+D2EN3DA2+D2EN7DA2
	DENADA5=D1EN3DA5+D1EN7DA5+D2EN3DA5+D2EN7DA5
	DENADAZ=D1EN3DAZ+D1EN7DAZ+D2EN3DAZ+D2EN7DAZ
	DENADBZ=D1EN3DBZ+D1EN7DBZ+D2EN3DBZ+D2EN7DBZ
	END IF
	IF(IPERM.EQ.4) THEN
	DENADA1=D1EN4DA1+D1EN8DA1+D2EN4DA1+D2EN8DA1
	DENADA2=D1EN4DA2+D1EN8DA2+D2EN4DA2+D2EN8DA2
	DENADA5=D1EN4DA5+D1EN8DA5+D2EN4DA5+D2EN8DA5
	DENADAZ=D1EN4DAZ+D1EN8DAZ+D2EN4DAZ+D2EN8DAZ
	DENADBZ=D1EN4DBZ+D1EN8DBZ+D2EN4DBZ+D2EN8DBZ
	END IF
C     -----------------------------------------------------
C     --- Calculate the Elec./Elec. Repulsion ---
C     -----------------------------------------------------
	DOVnDA1=DOV1DA1
	DOVnDA2=DOV1DA2
	DOVnDA5=DOV1DA5
	DOVnDAZ=DOV1DAZ
	DOVnDBZ=DOV1DBZ
              CALL ELREP(DFLAG,OV1,RX1(L),RY1(L),RZ1(L),AA1(L),
     1           RX1(K),RY1(K),RZ1(K),AA1(K),                      
     2           RX2(L),RY2(L),RZ2(L),AA2(L),
     3           RX2(K),RY2(K),RZ2(K),AA2(K),
     4           B12(L),B12(K),ER1,KLM)  
	DER1DA1=DERDA1
	DER1DA2=DERDA2
	DER1DA5=DERDA5
	DER1DAZ=DERDAZ
	DER1DBZ=DERDBZ



C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C
	IF (ISYMM) THEN
	DOVnDA1=DOV2DA1
	DOVnDA2=DOV2DA2
	DOVnDA5=DOV2DA5
	DOVnDAZ=DOV2DAZ
	DOVnDBZ=DOV2DBZ
              CALL ELREP(DFLAG,OV2,RX1(L),RY1(L),RZ1(L),AA1(L),
     1           RX2(K),RY2(K),RZ2(K),AA2(K),                      
     2           RX2(L),RY2(L),RZ2(L),AA2(L),
     3           RX1(K),RY1(K),RZ1(K),AA1(K),
     4           B12(L),B12(K),ER2,KLM)  
	DER2DA1=DERDA1
	DER2DA2=DERDA2
	DER2DA5=DERDA5
	DER2DAZ=DERDAZ
	DER2DBZ=DERDBZ

C
C     ...Nuclear Symmetrization Overlap Integral: P_[H1,H2]
C
	DOVnDA1=DOV3DA1
	DOVnDA2=DOV3DA2
	DOVnDA5=DOV3DA5
	DOVnDAZ=DOV3DAZ
	DOVnDBZ=DOV3DBZ
              CALL ELREP(DFLAG,OV3,RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
     4            B12(L),B12(K),ER3,KLM)
	DER3DA1=DERDA1
	DER3DA2=DERDA2
	DER3DA5=DERDA5
	DER3DAZ=DERDAZ
	DER3DBZ=DERDBZ


C
C     ...Electro/Nuclear Symmetrization Overlap Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
	DOVnDA1=DOV4DA1
	DOVnDA2=DOV4DA2
	DOVnDA5=DOV4DA5
	DOVnDAZ=DOV4DAZ
	DOVnDBZ=DOV4DBZ
              CALL ELREP(DFLAG,OV4,RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
     4            B12(L),B12(K),ER4,KLM)
	DER4DA1=DERDA1
	DER4DA2=DERDA2
	DER4DA5=DERDA5
	DER4DAZ=DERDAZ
	DER4DBZ=DERDBZ



	END IF

C              
	IF (ISYMM) THEN
		EER1=ER1+ER2+ER3+ER4
	ELSE
		EER1=ER1
	END IF
C
	IF (ISYMM) THEN
	DEERDA1=DER1DA1+DER2DA1+DER3DA1+DER4DA1
	DEERDA2=DER1DA2+DER2DA2+DER3DA2+DER4DA2
	DEERDA5=DER1DA5+DER2DA5+DER3DA5+DER4DA5
	DEERDAZ=DER1DAZ+DER2DAZ+DER3DAZ+DER4DAZ
	DEERDBZ=DER1DBZ+DER2DBZ+DER3DBZ+DER4DBZ
	ELSE
	DEERDA1=DER1DA1
	DEERDA2=DER1DA2
	DEERDA5=DER1DA5
	DEERDAZ=DER1DAZ
	DEERDBZ=DER1DBZ
	END IF
	IF(IPERM.EQ.1) THEN
	DEERDA1=DER1DA1
	DEERDA2=DER1DA2
	DEERDA5=DER1DA5
	DEERDAZ=DER1DAZ
	DEERDBZ=DER1DBZ
	END IF
	IF(IPERM.EQ.2) THEN
	DEERDA1=DER2DA1
	DEERDA2=DER2DA2
	DEERDA5=DER2DA5
	DEERDAZ=DER2DAZ
	DEERDBZ=DER2DBZ
	END IF
	IF(IPERM.EQ.3) THEN
	DEERDA1=DER3DA1
	DEERDA2=DER3DA2
	DEERDA5=DER3DA5
	DEERDAZ=DER3DAZ
	DEERDBZ=DER3DBZ
	END IF
	IF(IPERM.EQ.4) THEN
	DEERDA1=DER4DA1
	DEERDA2=DER4DA2
	DEERDA5=DER4DA5
	DEERDAZ=DER4DAZ
	DEERDBZ=DER4DBZ
	END IF
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
              H(L,K)=ER+T-EN+(1.0D+00/1.4011D+00)*S(L,K)
C
	DHDA1(L,K)=DEERDA1+DTDA1-DENADA1+(1.0D+00/1.4011D+00)*DSDA1(L,K)
	DHDA2(L,K)=DEERDA2+DTDA2-DENADA2+(1.0D+00/1.4011D+00)*DSDA2(L,K)
	DHDA5(L,K)=DEERDA5+DTDA5-DENADA5+(1.0D+00/1.4011D+00)*DSDA5(L,K)
	DHDAZ(L,K)=DEERDAZ+DTDAZ-DENADAZ+(1.0D+00/1.4011D+00)*DSDAZ(L,K)
	DHDBZ(L,K)=DEERDBZ+DTDBZ-DENADBZ+(1.0D+00/1.4011D+00)*DSDBZ(L,K)

C
	IF (L.EQ.K) THEN
		DHDA1(L,K)=2.00D+00*DHDA1(L,K)
		DHDA2(L,K)=2.00D+00*DHDA2(L,K)
		DHDA5(L,K)=2.00D+00*DHDA5(L,K)
		DHDAZ(L,K)=2.00D+00*DHDAZ(L,K)
		DHDBZ(L,K)=2.00D+00*DHDBZ(L,K)
C
		DSDA1(L,K)=2.00D+00*DSDA1(L,K)
		DSDA2(L,K)=2.00D+00*DSDA2(L,K)
		DSDA5(L,K)=2.00D+00*DSDA5(L,K)
		DSDAZ(L,K)=2.00D+00*DSDAZ(L,K)
		DSDBZ(L,K)=2.00D+00*DSDBZ(L,K)
	END IF
C     ----------------------------------------------------- 
          END DO 
      END DO

C	STOP 
	
C     ----------------------------------------------------- 
C      CALL TAB(S,M,M,MMAX,MMAX)                  
C      CALL TAB(H,M,M,MMAX,MMAX)
C      CALL ENERGY(M,H,S,EMIN,CO,KIM)
C		STOP
C     -----------------------------------------------------      
      SUMH=0.00D+00
	SUMS=0.00D+00
	DO J=1,M
		DO I=1,M
			SUMH=SUMH+EVEC(I)*H(I,J)*EVEC(J)
			SUMS=SUMS+EVEC(I)*S(I,J)*EVEC(J)
		END DO
	END DO
C
	ENG=SUMH/SUMS
	EMIN=ENG
C
	DO I=1,M
		K=1
		DO J=1,M
			DH(K,I)=DHDA1(J,I)
			DH(K+1,I)=DHDA2(J,I)
			DH(K+2,I)=DHDA5(J,I)
			DH(K+3,I)=DHDAZ(J,I)
			DH(K+4,I)=DHDBZ(J,I)
C
			DS(K,I)=DSDA1(J,I)
			DS(K+1,I)=DSDA2(J,I)
			DS(K+2,I)=DSDA5(J,I)
			DS(K+3,I)=DSDAZ(J,I)
			DS(K+4,I)=DSDBZ(J,I)
C
			K=K+5
		END DO
	END DO
C
C      CALL TAB(DH,M*5,M,MMAX*5,MMAX)                  
C	CALL BIOUT(' ')
C     CALL TAB(DS,M*5,M,MMAX*5,MMAX)                  

C	STOP
C
	SCALE=1.00D+00/SUMS
C
	DO I=1,M
		DO J=1,M*5
			GRADM(I,J)=DH(J,I)-ENG*DS(J,I)
			GRADM(I,J)=SCALE*GRADM(I,J)
		END DO
		END DO
C
c	CALL BIOUT('GRADM-MATRIX:')
c		CALL TAB(GRADM,M,M*5,MMAX,MMAX*5)
c	CALL BIOUT(' ')

c	STOP


	NN=5
	NB=M
	DO J=1,NB
		DO K=1,NN
			IC=(J-1)*NN+K
			GRADA(IC)=0.00D+00
			DO I=1,NB
				IF(I.NE.J) THEN
					GRADA(IC)=GRADA(IC)+2.00D+00*EVEC(I)*EVEC(J)*
	-						GRADM(I,IC)
				ELSE
					GRADA(IC)=GRADA(IC)+EVEC(I)*EVEC(J)*GRADM(I,IC)
				END IF
			END DO
		END DO
	END DO
C
      DO I=1,NB
		GRADC(I) = 0.00D+00
	END DO
	
      DO J=1,NB
		DO I=1,NB
		  GRADC(I)=GRADC(I)+2.00D+00*SCALE*( H(I,J)-ENG*S(I,J) )*
	-				EVEC(J)
		END DO
	END DO

	K=1
	DO I=1,NB*5,5
		DA1=GRADA(I)
		DA2=GRADA(I+1)
		DB=GRADA(I+2)
	GRADA(I)=DA1*(2.00D+00*TTT1(K)+TTT3(K))+DA2*TTT3(K)-DB*TTT3(K)
	GRADA(I+1)=2.00D+00*TTT2(K)*DA2
	GRADA(I+2)=DA1*TTT1(K)+(2.00D+00*TTT3(K)+TTT1(K))*DA2-DB*TTT1(K)
	K=K+1
      END DO
C
      NX=M*5
	DO I=1,NX
		GRAD(I) = GRADA(I)
	END DO
	
      DO I=1,NB
		GRAD(I+NX) = GRADC(I)
	END DO
C     -----------------------------------------------------      
      DEALLOCATE(	TEMP,DHDA1,DHDA2,DHDA5,DHDAZ,DHDBZ,DSDA1,DSDA2,
	-			DSDA5,DSDAZ,DSDBZ,GRADM,DS,DH,GRADA,GRADC)
C
C
C
	ELSE
C
C
C
      DO L=1,M
          DO K=1,L
C     -----------------------------------------------------
C     --- Calculate OVERLAP INTEGRAL ---
C     -----------------------------------------------------
             CALL OVERLP(M,DFLAG,RX1(L),RY1(L),RZ1(L),AA1(L),
     1            RX1(K),RY1(K),RZ1(K),AA1(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            RX2(K),RY2(K),RZ2(K),AA2(K),
     4            B12(L),B12(K),OV1,KLM)
C
C
	IF(ISYMM) THEN
		CALL OVERLP(M,DFLAG,RX1(L),RY1(L),RZ1(L),AA1(L),
     1            RX2(K),RY2(K),RZ2(K),AA2(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            RX1(K),RY1(K),RZ1(K),AA1(K),
     4            B12(L),B12(K),OV2,KLM)
C
C
C     ...Nuclear Symmetrization Overlap Integral: P_[H1,H2]
C
              CALL OVERLP(M,DFLAG,RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
     4            B12(L),B12(K),OV3,KLM)
C
C
C     ...Electro/Nuclear Symmetrization Overlap Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL OVERLP(M,DFLAG,RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
     4            B12(L),B12(K),OV4,KLM)
C              
	END IF
C     -----------------------------------------------------
	IF(ISYMM) THEN
              S(L,K)=OV1+OV2+OV3+OV4
		ELSE
              S(L,K)=OV1
	END IF
	IF(IPERM.EQ.1) THEN
		S(L,K)=OV1
	END IF
	IF(IPERM.EQ.2) THEN
		S(L,K)=OV2
	END IF
	IF(IPERM.EQ.3) THEN
		S(L,K)=OV3
	END IF
	IF(IPERM.EQ.4) THEN
		S(L,K)=OV4
	END IF
C
C	IF (K.EQ.1) THEN
C		WRITE(6,*)'S(',L,',1)=',S(L,K)
C	END IF
C     -----------------------------------------------------
C     --- Calculate KINETIC ENERGY OF ELECTRON 1 ---
C     -----------------------------------------------------
              CALL TK(DFLAG,OV1,
     1            RX1(L),RY1(L),RZ1(L),AA1(L),
     2            RX1(K),RY1(K),RZ1(K),AA1(K),                      
     3            RX2(L),RY2(L),RZ2(L),AA2(L),
     4            RX2(K),RY2(K),RZ2(K),AA2(K),
     5            B12(L),B12(K),T1)  
C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
	IF (ISYMM) THEN
              CALL TK(DFLAG,OV2,
     1            RX1(L),RY1(L),RZ1(L),AA1(L),
     2            RX2(K),RY2(K),RZ2(K),AA2(K),                      
     3            RX2(L),RY2(L),RZ2(L),AA2(L),
     4            RX1(K),RY1(K),RZ1(K),AA1(K),
     5            B12(L),B12(K),T2)  
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
              CALL TK(DFLAG,OV3,
     1            RX1(L),RY1(L),RZ1(L),AA1(L),
     2            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
     3            RX2(L),RY2(L),RZ2(L),AA2(L),
     4            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
     5            B12(L),B12(K),T3)
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
              CALL TK(DFLAG,OV4,
     1            RX1(L),RY1(L),RZ1(L),AA1(L),
     2            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
     3            RX2(L),RY2(L),RZ2(L),AA2(L),
     4            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
     5            B12(L),B12(K),T4)
	END IF
C     -----------------------------------------------------
C     --- Calculate KINETIC ENERGY OF ELECTRON 2 ---
C     -----------------------------------------------------
             CALL TK(DFLAG,OV1,
     1            RX2(L),RY2(L),RZ2(L),AA2(L),
     2            RX2(K),RY2(K),RZ2(K),AA2(K),                      
     3            RX1(L),RY1(L),RZ1(L),AA1(L),
     4            RX1(K),RY1(K),RZ1(K),AA1(K),
     5            B12(L),B12(K),T5)  
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
	IF (ISYMM) THEN
              CALL TK(DFLAG,OV2,
     1            RX2(L),RY2(L),RZ2(L),AA2(L),
     2            RX1(K),RY1(K),RZ1(K),AA1(K),                      
     3            RX1(L),RY1(L),RZ1(L),AA1(L),
     4            RX2(K),RY2(K),RZ2(K),AA2(K),
     5            B12(L),B12(K),T6) 
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
              CALL TK(DFLAG,OV3,
     1            RX2(L),RY2(L),RZ2(L),AA2(L),
     2            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
     3            RX1(L),RY1(L),RZ1(L),AA1(L),
     4            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
     5            B12(L),B12(K),T7)
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
              CALL TK(DFLAG,OV4,
     1            RX2(L),RY2(L),RZ2(L),AA2(L),
     2            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
     3            RX1(L),RY1(L),RZ1(L),AA1(L),
     4            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
     5            B12(L),B12(K),T8)
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

C	write(6,*)L,K,T
C     -----------------------------------------------------
C     --- Calculate Nuc./Elec. Attractn. FOR ELECTRON 1 ---
C     -----------------------------------------------------
              CALL ENUCL(DFLAG,OV1,RX1(L),RY1(L),RZ1(L),AA1(L),
     1           RX1(K),RY1(K),RZ1(K),AA1(K),                      
     2           RX2(L),RY2(L),RZ2(L),AA2(L),
     3           RX2(K),RY2(K),RZ2(K),AA2(K),
     4           B12(L),B12(K),
     5           EX1,EY1,EZ1,EN1,KLM)
C

c	WRITE(6,*)'EN1=',EN1

C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C
	IF (ISYMM) THEN
              CALL ENUCL(DFLAG,OV2,RX1(L),RY1(L),RZ1(L),AA1(L),
     1           RX2(K),RY2(K),RZ2(K),AA2(K),                      
     2           RX2(L),RY2(L),RZ2(L),AA2(L),
     3           RX1(K),RY1(K),RZ1(K),AA1(K),
     4           B12(L),B12(K),
     5           EX1,EY1,EZ1,EN2,KLM)                       
C
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
C
              CALL ENUCL(DFLAG,OV3,RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
     4            B12(L),B12(K),
     5            EX1,EY1,EZ1,EN3,KLM)
C
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL ENUCL(DFLAG,OV4,RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
     4            B12(L),B12(K),
     5            EX1,EY1,EZ1,EN4,KLM)
	END IF
C  
              CALL ENUCL(DFLAG,OV1,RX1(L),RY1(L),RZ1(L),AA1(L),
     1           RX1(K),RY1(K),RZ1(K),AA1(K),                      
     2           RX2(L),RY2(L),RZ2(L),AA2(L),
     3           RX2(K),RY2(K),RZ2(K),AA2(K),
     4           B12(L),B12(K),
     5           EX2,EY2,EZ2,EN5,KLM)
C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C
	IF (ISYMM) THEN
              CALL ENUCL(DFLAG,OV2,RX1(L),RY1(L),RZ1(L),AA1(L),
     1           RX2(K),RY2(K),RZ2(K),AA2(K),                      
     2           RX2(L),RY2(L),RZ2(L),AA2(L),
     3           RX1(K),RY1(K),RZ1(K),AA1(K),
     4           B12(L),B12(K),
     5           EX2,EY2,EZ2,EN6,KLM)                       
C
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
C
              CALL ENUCL(DFLAG,OV3,RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
     4            B12(L),B12(K),
     5            EX2,EY2,EZ2,EN7,KLM)
C
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL ENUCL(DFLAG,OV4,RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
     4            B12(L),B12(K),
     5            EX2,EY2,EZ2,EN8,KLM)
	END IF
C 
C     -----------------------------------------------------
C     --- Calculate Nuc./Elec. Attractn. FOR ELECTRON 2 ---
C     -----------------------------------------------------
              CALL ENUCL(DFLAG,OV1,RX2(L),RY2(L),RZ2(L),AA2(L),
     1           RX2(K),RY2(K),RZ2(K),AA2(K),                      
     2           RX1(L),RY1(L),RZ1(L),AA1(L),
     3           RX1(K),RY1(K),RZ1(K),AA1(K),
     4           B12(L),B12(K),
     5           EX1,EY1,EZ1,EN9,KLM)                       
C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C
	IF (ISYMM) THEN
              CALL ENUCL(DFLAG,OV2,RX2(L),RY2(L),RZ2(L),AA2(L),
     1           RX1(K),RY1(K),RZ1(K),AA1(K),                      
     2           RX1(L),RY1(L),RZ1(L),AA1(L),
     3           RX2(K),RY2(K),RZ2(K),AA2(K),
     4           B12(L),B12(K),
     5           EX1,EY1,EZ1,EN10,KLM)                       
C
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
C
              CALL ENUCL(DFLAG,OV3,RX2(L),RY2(L),RZ2(L),AA2(L),
     1            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
     2            RX1(L),RY1(L),RZ1(L),AA1(L),
     3            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
     4            B12(L),B12(K),
     5            EX1,EY1,EZ1,EN11,KLM)
C
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL ENUCL(DFLAG,OV4,RX2(L),RY2(L),RZ2(L),AA2(L),
     1            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
     2            RX1(L),RY1(L),RZ1(L),AA1(L),
     3            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
     4            B12(L),B12(K),
     5            EX1,EY1,EZ1,EN12,KLM)
	END IF
C 
              CALL ENUCL(DFLAG,OV1,RX2(L),RY2(L),RZ2(L),AA2(L),
     1           RX2(K),RY2(K),RZ2(K),AA2(K),                      
     2           RX1(L),RY1(L),RZ1(L),AA1(L),
     3           RX1(K),RY1(K),RZ1(K),AA1(K),
     4           B12(L),B12(K),
     5           EX2,EY2,EZ2,EN13,KLM)                       
C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C
	IF (ISYMM) THEN
              CALL ENUCL(DFLAG,OV2,RX2(L),RY2(L),RZ2(L),AA2(L),
     1           RX1(K),RY1(K),RZ1(K),AA1(K),                      
     2           RX1(L),RY1(L),RZ1(L),AA1(L),
     3           RX2(K),RY2(K),RZ2(K),AA2(K),
     4           B12(L),B12(K),
     5           EX2,EY2,EZ2,EN14,KLM)  
C
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
C
              CALL ENUCL(DFLAG,OV3,RX2(L),RY2(L),RZ2(L),AA2(L),
     1            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
     2            RX1(L),RY1(L),RZ1(L),AA1(L),
     3            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
     4            B12(L),B12(K),
     5            EX2,EY2,EZ2,EN15,KLM)
C
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL ENUCL(DFLAG,OV4,RX2(L),RY2(L),RZ2(L),AA2(L),
     1            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
     2            RX1(L),RY1(L),RZ1(L),AA1(L),
     3            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
     4            B12(L),B12(K),
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
              H(L,K)=ER+T-EN+(1.0D+00/1.4011D+00)*S(L,K)
C
c	WRITE(6,*)'ER=',ER,L,K
c	WRITE(6,*)'EN=',EN,L,K
C     ----------------------------------------------------- 
          END DO 
      END DO
	DO L=1,M
		DO K=1,L
			IF (K .NE. L) THEN
				S(K,L)=S(L,K)
				H(K,L)=H(L,K)
			END IF
		END DO
	END DO
C     ----------------------------------------------------- 
C      CALL TAB(S,M,M,MMAX,MMAX)                  
C      CALL TAB(H,M,M,MMAX,MMAX)
	IF (DFLAG .NE. 3) CALL ENERGY(M,H,S,EMIN,CO,KIM)
C		STOP
C     -----------------------------------------------------      
	IF (DFLAG .EQ. 3) THEN
      SUMH=0.00D+00
	SUMS=0.00D+00
	DO J=1,M
		DO I=1,M
			SUMH=SUMH+EVEC(I)*H(I,J)*EVEC(J)
			SUMS=SUMS+EVEC(I)*S(I,J)*EVEC(J)
		END DO
	END DO
C
	ENG=SUMH/SUMS
	EMIN=ENG
	END IF
C
C     -----------------------------------------------------      
      DEALLOCATE(	TEMP,DHDA1,DHDA2,DHDA5,DHDAZ,DHDBZ,DSDA1,DSDA2,
	-			DSDA5,DSDAZ,DSDBZ,GRADM,DS,DH,GRADA,GRADC)
	END IF

      RETURN
      END