SUBROUTINE Gradient_Projection_1

	USE scalars_to_be_shared			  
	USE	vectors_to_be_shared			  
	USE	matrices_to_be_shared
	USE projection_1_module
	
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	IMPLICIT INTEGER (I-N)

	DO L=1,M

		J=1

		DO K=1,M


		CALL dOVERLAP(RX1(L),RY1(L),RZ1(L),AA1(L),				&
						RX2(K),RY2(K),RZ2(K),AA2(K),            &           
						RX2(L),RY2(L),RZ2(L),AA2(L),			&
						RX1(K),RY1(K),RZ1(K),AA1(K),			&
						B12(L),B12(K),							&
						D1,D2,D3,D4,D5)

			DSDA1=D1
			DSDA2=D2
			DSDA5=D3
			DSDAZ=D4
			DSDBZ=D5

              CALL DTKINET(OVV2(L,K),					&
				D1,D2,D3,D4,D5,							&
				RX1(L),RY1(L),RZ1(L),AA1(L),			&
                RX2(K),RY2(K),RZ2(K),AA2(K),			&                   
                RX2(L),RY2(L),RZ2(L),AA2(L),			&
                RX1(K),RY1(K),RZ1(K),AA1(K),			&
                B12(L),B12(K),							&
				DT1,DT2,DT3,DT4,DT5)  

		DTDA1=DT1
		DTDA2=DT2
		DTDA5=DT3
		DTDAZ=DT4
		DTDBZ=DT5

              CALL DTKINET(OVV2(L,K),					&
				D2,D1,D3,D5,D4,							&
                RX2(L),RY2(L),RZ2(L),AA2(L),			&
                RX1(K),RY1(K),RZ1(K),AA1(K),			&                   
                RX1(L),RY1(L),RZ1(L),AA1(L),			&
                RX2(K),RY2(K),RZ2(K),AA2(K),			&
                B12(L),B12(K),							&
				DT1,DT2,DT3,DT4,DT5) 

		DTDA1=DTDA1+DT2
		DTDA2=DTDA2+DT1
		DTDA5=DTDA5+DT3
		DTDAZ=DTDAZ+DT5
		DTDBZ=DTDBZ+DT4

              CALL dELECREP(OVV2(L,K),									&
				D1,D2,D3,D4,D5,											&
				RX1(L),RY1(L),RZ1(L),AA1(L),							&
                RX2(K),RY2(K),RZ2(K),AA2(K),                      		&
                RX2(L),RY2(L),RZ2(L),AA2(L),							&
                RX1(K),RY1(K),RZ1(K),AA1(K),							&
                B12(L),B12(K),											&
				DER1,DER2,DER3,DER4,DER5)  

		DERDA1=DER1
		DERDA2=DER2
		DERDA5=DER3
		DERDAZ=DER4
		DERDBZ=DER5



              CALL dENUCLEAR(OVV2(L,K),				&
				D1,D2,D3,D4,D5,						&
				RX1(L),RY1(L),RZ1(L),AA1(L),		&
                RX2(K),RY2(K),RZ2(K),AA2(K),        &
                RX2(L),RY2(L),RZ2(L),AA2(L),		&
                RX1(K),RY1(K),RZ1(K),AA1(K),		&
                B12(L),B12(K),						&
                EX1,EY1,EZ1,-ONE,					&
				DEN1,DEN2,DEN3,DEN4,DEN5)                       

		DENDA1=DEN1
		DENDA2=DEN2
		DENDA5=DEN3
		DENDAZ=DEN4
		DENDBZ=DEN5

              CALL dENUCLEAR(OVV2(L,K),			&
				D1,D2,D3,D4,D5,					&
				RX1(L),RY1(L),RZ1(L),AA1(L),	&
                RX2(K),RY2(K),RZ2(K),AA2(K),    &
                RX2(L),RY2(L),RZ2(L),AA2(L),	&
                RX1(K),RY1(K),RZ1(K),AA1(K),	&
                B12(L),B12(K),					&
                EX2,EY2,EZ2,ONE,				&
				DEN1,DEN2,DEN3,DEN4,DEN5)                       

		DENDA1=DENDA1+DEN1
		DENDA2=DENDA2+DEN2
		DENDA5=DENDA5+DEN3
		DENDAZ=DENDAZ+DEN4
		DENDBZ=DENDBZ+DEN5

              CALL dENUCLEAR(OVV2(L,K),				&
				D2,D1,D3,D5,D4,						&
				RX2(L),RY2(L),RZ2(L),AA2(L),		&
                RX1(K),RY1(K),RZ1(K),AA1(K),        &
                RX1(L),RY1(L),RZ1(L),AA1(L),		&
                RX2(K),RY2(K),RZ2(K),AA2(K),		&
                B12(L),B12(K),						&
                EX1,EY1,EZ1,-ONE,					&
				DEN1,DEN2,DEN3,DEN4,DEN5)                       

		DENDA1=DENDA1+DEN2
		DENDA2=DENDA2+DEN1
		DENDA5=DENDA5+DEN3
		DENDAZ=DENDAZ+DEN5
		DENDBZ=DENDBZ+DEN4

              CALL dENUCLEAR(OVV2(L,K),				&
				D2,D1,D3,D5,D4,						&
				RX2(L),RY2(L),RZ2(L),AA2(L),		&
                RX1(K),RY1(K),RZ1(K),AA1(K),        &
                RX1(L),RY1(L),RZ1(L),AA1(L),		&
                RX2(K),RY2(K),RZ2(K),AA2(K),		&
                B12(L),B12(K),						&
                EX2,EY2,EZ2,ONE,					&
				DEN1,DEN2,DEN3,DEN4,DEN5)  

		DENDA1=DENDA1+DEN2
		DENDA2=DENDA2+DEN1
		DENDA5=DENDA5+DEN3
		DENDAZ=DENDAZ+DEN5
		DENDBZ=DENDBZ+DEN4


	DHDA1=HALF*DTDA1+DERDA1+(ONE/BD)*DSDA1-DENDA1
	DHDA2=HALF*DTDA2+DERDA2+(ONE/BD)*DSDA2-DENDA2
	DHDA5=HALF*DTDA5+DERDA5+(ONE/BD)*DSDA5-DENDA5
	DHDAZ=HALF*DTDAZ+DERDAZ+(ONE/BD)*DSDAZ-DENDAZ
	DHDBZ=HALF*DTDBZ+DERDBZ+(ONE/BD)*DSDBZ-DENDBZ

	IF (L.EQ.K) THEN
		DHDA1=TWO*DHDA1
		DHDA2=TWO*DHDA2
		DHDA5=TWO*DHDA5
		DHDAZ=TWO*DHDAZ
		DHDBZ=TWO*DHDBZ
		DSDA1=TWO*DSDA1
		DSDA2=TWO*DSDA2
		DSDA5=TWO*DSDA5
		DSDAZ=TWO*DSDAZ
		DSDBZ=TWO*DSDBZ

	END IF


			DH1(J,L)  =DHDA1
			DH1(J+1,L)=DHDA2
			DH1(J+2,L)=DHDA5
			DH1(J+3,L)=DHDAZ
			DH1(J+4,L)=DHDBZ

			DS1(J,L)  =DSDA1
			DS1(J+1,L)=DSDA2
			DS1(J+2,L)=DSDA5
			DS1(J+3,L)=DSDAZ
			DS1(J+4,L)=DSDBZ

			J=J+5


		END DO

	END DO


END SUBROUTINE Gradient_Projection_1