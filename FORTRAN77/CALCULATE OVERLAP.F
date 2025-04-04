SUBROUTINE Calculate_Overlap

	USE scalars_to_be_shared			  
	USE	vectors_to_be_shared			  
	USE	matrices_to_be_shared			  


	DOUBLE PRECISION :: OV

	DO L=1,M
		DO K=1,L
             CALL OVERLAP(RX1(L),RY1(L),RZ1(L),AA1(L),	&
						RX1(K),RY1(K),RZ1(K),AA1(K),    &                  
     					RX2(L),RY2(L),RZ2(L),AA2(L),	&
     					RX2(K),RY2(K),RZ2(K),AA2(K),	&
    					B12(L),B12(K),OV)
		OVV1(L,K)=OV

!     ...Electron Symmetrization Overlap Integral: P_[e1,e2]

		CALL OVERLAP(RX1(L),RY1(L),RZ1(L),AA1(L),		&
                RX2(K),RY2(K),RZ2(K),AA2(K),            &           
				RX2(L),RY2(L),RZ2(L),AA2(L),			&
                RX1(K),RY1(K),RZ1(K),AA1(K),			&
                B12(L),B12(K),OV)

		OVV2(L,K)=OV

!     ...Nuclear Symmetrization Overlap Integral: P_[H1,H2]


		OVV3(L,K)=OV

!     ...Electro/Nuclear Symmetrization Overlap Integral:
!         P_[e1,e2]*P_[H1,H2]


		OVV4(L,K)=OV



              S(L,K)=OVV1(L,K)+OVV2(L,K)+OVV3(L,K)+OVV4(L,K)


	IF (L .NE. K) THEN
		S(K,L)=S(L,K)
		OVV1(K,L)=OVV1(L,K)
		OVV2(K,L)=OVV2(L,K)
		OVV3(K,L)=OVV3(L,K)
		OVV4(K,L)=OVV4(L,K)
	END IF

				
		END DO
	END DO



END SUBROUTINE Calculate_Overlap

