      SUBROUTINE SOVERLAP(ISYMM,M,
	-					AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,
	-					OV1,OV2,OV3,OV4,S)
C     -----------------------------------------------------
      IMPLICIT NONE
	INTEGER L,K,M,ISYMM
	REAL*8 OV
      REAL*8, DIMENSION(M) ::	AA1,RX1,RY1,RZ1,
	-						AA2,RX2,RY2,RZ2,B12
	REAL*8, DIMENSION(M,M) :: S
	REAL*8, DIMENSION(M,M) :: OV1,OV2,OV3,OV4
C
	DO L=1,M
		DO K=1,L
             CALL OVERLAP(RX1(L),RY1(L),RZ1(L),AA1(L),
     -					RX1(K),RY1(K),RZ1(K),AA1(K),                      
     -					RX2(L),RY2(L),RZ2(L),AA2(L),
     -					RX2(K),RY2(K),RZ2(K),AA2(K),
     -					B12(L),B12(K),OV)
		OV1(L,K)=OV
C
	IF(ISYMM) THEN
C
C     ...Electron Symmetrization Overlap Integral: P_[e1,e2]
C
		CALL OVERLAP(RX1(L),RY1(L),RZ1(L),AA1(L),
     -            RX2(K),RY2(K),RZ2(K),AA2(K),                      
     -            RX2(L),RY2(L),RZ2(L),AA2(L),
     -            RX1(K),RY1(K),RZ1(K),AA1(K),
     -            B12(L),B12(K),OV)
		OV2(L,K)=OV
C
C     ...Nuclear Symmetrization Overlap Integral: P_[H1,H2]
C
              CALL OVERLAP(RX1(L),RY1(L),RZ1(L),AA1(L),
     -            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
     -            RX2(L),RY2(L),RZ2(L),AA2(L),
     -            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
     -            B12(L),B12(K),OV)
		OV3(L,K)=OV
C
C     ...Electro/Nuclear Symmetrization Overlap Integral:
C         P_[e1,e2]*P_[H1,H2]
C
              CALL OVERLAP(RX1(L),RY1(L),RZ1(L),AA1(L),
     -            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
     -            RX2(L),RY2(L),RZ2(L),AA2(L),
     -            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
     -            B12(L),B12(K),OV)
		OV4(L,K)=OV
C
	END IF
C
	IF(ISYMM) THEN
              S(L,K)=OV1(L,K)+OV2(L,K)+OV3(L,K)+OV4(L,K)
		ELSE
              S(L,K)=OV1(L,K)
	END IF
C
	IF (L .NE. K) THEN
		S(K,L)=S(L,K)
		OV1(K,L)=OV1(L,K)
		OV2(K,L)=OV2(L,K)
		OV3(K,L)=OV3(L,K)
		OV4(K,L)=OV4(L,K)
	END IF
C				
		END DO
	END DO
C
	RETURN
	END