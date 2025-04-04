SUBROUTINE Calculate_Overlap_Gradient

	USE scalars_to_be_shared			  
	USE	vectors_to_be_shared			  
	USE	matrices_to_be_shared
	
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	IMPLICIT INTEGER (I-N)

	DO K=1,M

		J=1

		DO L=1,M

IF (L .NE. K) THEN
             CALL dOVERLAP(RX1(L),RY1(L),RZ1(L),AA1(L),			&
						RX1(K),RY1(K),RZ1(K),AA1(K),			&                  
     					RX2(L),RY2(L),RZ2(L),AA2(L),			&
     					RX2(K),RY2(K),RZ2(K),AA2(K),			&
    					B12(L),B12(K),							&
						D1,D2,D3,D4,D5)
ELSE
	D1=ZERO
	D2=ZERO
	D3=ZERO
	D4=ZERO
	D5=ZERO
END IF
			DSDA1=D1
			DSDA2=D2
			DSDA5=D3
			DSDAZ=D4
			DSDBZ=D5


!     ...Electron Symmetrization Overlap Integral: P_[e1,e2]

		CALL dOVERLAP(RX1(L),RY1(L),RZ1(L),AA1(L),				&
						RX2(K),RY2(K),RZ2(K),AA2(K),            &           
						RX2(L),RY2(L),RZ2(L),AA2(L),			&
						RX1(K),RY1(K),RZ1(K),AA1(K),			&
						B12(L),B12(K),							&
						D1,D2,D3,D4,D5)

			DSDA1=DSDA1+D1
			DSDA2=DSDA2+D2
			DSDA5=DSDA5+D3
			DSDAZ=DSDAZ+D4
			DSDBZ=DSDBZ+D5

!     ...Nuclear Symmetrization Overlap Integral: P_[H1,H2]

		CALL dOVERLAP(RX1(L),RY1(L),RZ1(L),AA1(L),				&
						-RX1(K),-RY1(K),-RZ1(K),AA1(K),         &             
						RX2(L),RY2(L),RZ2(L),AA2(L),			&
						-RX2(K),-RY2(K),-RZ2(K),AA2(K),			&
						B12(L),B12(K),							&
						D1,D2,D3,D4,D5)

			DSDA1=DSDA1+D1
			DSDA2=DSDA2+D2
			DSDA5=DSDA5+D3
			DSDAZ=DSDAZ+D4
			DSDBZ=DSDBZ+D5


!     ...Electro/Nuclear Symmetrization Overlap Integral:
!         P_[e1,e2]*P_[H1,H2]

		CALL dOVERLAP(RX1(L),RY1(L),RZ1(L),AA1(L),				&
						-RX2(K),-RY2(K),-RZ2(K),AA2(K),         &             
						RX2(L),RY2(L),RZ2(L),AA2(L),			&
						-RX1(K),-RY1(K),-RZ1(K),AA1(K),			&
						B12(L),B12(K),							&
						D1,D2,D3,D4,D5)

			DSDA1=DSDA1+D1
			DSDA2=DSDA2+D2
			DSDA5=DSDA5+D3
			DSDAZ=DSDAZ+D4
			DSDBZ=DSDBZ+D5


		IF (L.EQ.K) THEN
			DSDA1=TWO*DSDA1
			DSDA2=TWO*DSDA2
			DSDA5=TWO*DSDA5
			DSDAZ=TWO*DSDAZ
			DSDBZ=TWO*DSDBZ
		END IF

			DS(J,K)  =DSDA1
			DS(J+1,K)=DSDA2
			DS(J+2,K)=DSDA5
			DS(J+3,K)=DSDAZ
			DS(J+4,K)=DSDBZ

			J=J+5


		END DO
	END DO



END SUBROUTINE Calculate_Overlap_Gradient