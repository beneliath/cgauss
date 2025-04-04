!      PROJECT TITLE:  correlated gaussian geminals for molecular hydrogen
!MAIN EXECUTION FILE:  cgauss.F90
!
!   Program FILENAME:  cgauss.exe
!        Lead Author:  D. Gilmore
!    Other Author(s):  L. Adamowicz, D. Kinghorn
!
!        LAST EDITED:  22may96:wed/18:21
!

!NOTES FOR FUTURE EDITS:
!===================================================================================================
!	1.  unnecessary sharing of modules should be eliminated
!===================================================================================================
!
!

!-------------------------------------------------------&   !=======================================
SUBROUTINE Calculate_Total_Gradient
!-----------------------------------------------------------
	USE  scalars_to_be_shared			  
	USE	 vectors_to_be_shared			  
	USE	matrices_to_be_shared

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	IMPLICIT          INTEGER     (I-N)
!-----------------------------------------------------------
!DO L=1,M
!DO K=1,M
!RMODIFY=1.00D-8
!			CALL OVERLAP(RZ1(L),AA1(L),RZ1(K),AA1(K),	&                  
!  					RZ2(L),AA2(L),RZ2(K),AA2(K),B12(L),	&
!					B12(K),OV1)
!AA1(L)=AA1(L)+RMODIFY
!			CALL OVERLAP(RZ1(L),AA1(L),RZ1(K),AA1(K),	&                  
! 					RZ2(L),AA2(L),RZ2(K),AA2(K),B12(L),	&
!					B12(K),OV2)
!FINITE=(OV2-OV1)/RMODIFY
!WRITE(6,*)L,K,'FINITE    dA2=',FINITE
!AA1(L)=AA1(L)-RMODIFY
!END DO
!END DO
!CALL BIOUT(' ')

	DO L=1,M
		J=1
		DO K=1,M

             CALL dOVERLAP(RZ1(L),AA1(L),RZ1(K),AA1(K),	&
					RZ2(L),AA2(L),RZ2(K),AA2(K),B12(L),	&
					B12(K),D1,D2,D3,D4,D5,OV1)
		DSDA1=D1
		DSDA2=D2
		DSDA5=D3
		DSDAZ=D4
		DSDBZ=D5

!IF (L .EQ. K) D1=D1*TWO
!WRITE(6,*)L,K,'ANALYTIC (R1)=',D1

              CALL DTKINET(OV1,							&
					D1,D2,D3,D4,D5,						&
					RZ1(L),AA1(L),RZ1(K),AA1(K),		&                   
					RZ2(L),AA2(L),RZ2(K),AA2(K),		&
					B12(L),B12(K),						&
					DT1,DT2,DT3,DT4,DT5) 
		DTDA1=DT1
		DTDA2=DT2
		DTDA5=DT3
		DTDAZ=DT4
		DTDBZ=DT5

             CALL DTKINET(OV1,							&
					D2,D1,D3,D5,D4,						&
					RZ2(L),AA2(L),RZ2(K),AA2(K),		&                   
					RZ1(L),AA1(L),RZ1(K),AA1(K),		&
					B12(L),B12(K),						&
					DT1,DT2,DT3,DT4,DT5)  
		DTDA1=DTDA1+DT2
		DTDA2=DTDA2+DT1
		DTDA5=DTDA5+DT3
		DTDAZ=DTDAZ+DT5
		DTDBZ=DTDBZ+DT4

              CALL dELECREP(OV1,						&
				D1,D2,D3,D4,D5,							&
				RZ1(L),AA1(L),RZ1(K),AA1(K),            &
                RZ2(L),AA2(L),RZ2(K),AA2(K),			&
                B12(L),B12(K),							&
				DER1,DER2,DER3,DER4,DER5)  
		DERDA1=DER1
		DERDA2=DER2
		DERDA5=DER3
		DERDAZ=DER4
		DERDBZ=DER5

              CALL dENUCLEAR(OV1,						&
				D1,D2,D3,D4,D5,							&
				RZ1(L),AA1(L),RZ1(K),AA1(K),			&
                RZ2(L),AA2(L),RZ2(K),AA2(K),			&
                B12(L),B12(K),							&
				EZ1,-ONE,								&
				DEN1,DEN2,DEN3,DEN4,DEN5)
		DENDA1=DEN1
		DENDA2=DEN2
		DENDA5=DEN3
		DENDAZ=DEN4
		DENDBZ=DEN5

              CALL dENUCLEAR(OV1,						&
				D1,D2,D3,D4,D5,							&
				RZ1(L),AA1(L),RZ1(K),AA1(K),			&
                RZ2(L),AA2(L),RZ2(K),AA2(K),			&
                B12(L),B12(K),							&
				EZ2,ONE,								&
				DEN1,DEN2,DEN3,DEN4,DEN5)
		DENDA1=DENDA1+DEN1
		DENDA2=DENDA2+DEN2
		DENDA5=DENDA5+DEN3
		DENDAZ=DENDAZ+DEN4
		DENDBZ=DENDBZ+DEN5

              CALL dENUCLEAR(OV1,						&
				D2,D1,D3,D5,D4,							&
				RZ2(L),AA2(L),RZ2(K),AA2(K),			&
                RZ1(L),AA1(L),RZ1(K),AA1(K),			&
                B12(L),B12(K),							&
				EZ1,-ONE,								&
				DEN1,DEN2,DEN3,DEN4,DEN5)                       
		DENDA1=DENDA1+DEN2
		DENDA2=DENDA2+DEN1
		DENDA5=DENDA5+DEN3
		DENDAZ=DENDAZ+DEN5
		DENDBZ=DENDBZ+DEN4

              CALL dENUCLEAR(OV1,						&
				D2,D1,D3,D5,D4,							&
				RZ2(L),AA2(L),RZ2(K),AA2(K),            &
                RZ1(L),AA1(L),RZ1(K),AA1(K),			&
                B12(L),B12(K),							&
                EZ2,ONE,								&
				DEN1,DEN2,DEN3,DEN4,DEN5)                       
		DENDA1=DENDA1+DEN2
		DENDA2=DENDA2+DEN1
		DENDA5=DENDA5+DEN3
		DENDAZ=DENDAZ+DEN5
		DENDBZ=DENDBZ+DEN4

!-----------------------------------------------------
		CALL dOVERLAP(RZ1(L),AA1(L),RZ2(K),AA2(K),      &           
				RZ2(L),AA2(L),RZ1(K),AA1(K),B12(L),		&
				B12(K),D1,D2,D3,D4,D5,OV2)
			DSDA1=DSDA1+D1
			DSDA2=DSDA2+D2
			DSDA5=DSDA5+D3
			DSDAZ=DSDAZ+D4
			DSDBZ=DSDBZ+D5


              CALL DTKINET(OV2,							&
				D1,D2,D3,D4,D5,							&
				RZ1(L),AA1(L),RZ2(K),AA2(K),			&                   
                RZ2(L),AA2(L),RZ1(K),AA1(K),			&
                B12(L),B12(K),DT1,DT2,DT3,DT4,DT5)  
		DTDA1=DTDA1+DT1
		DTDA2=DTDA2+DT2
		DTDA5=DTDA5+DT3
		DTDAZ=DTDAZ+DT4
		DTDBZ=DTDBZ+DT5

              CALL DTKINET(OV2,							&
				D2,D1,D3,D5,D4,							&
                RZ2(L),AA2(L),RZ1(K),AA1(K),			&                   
                RZ1(L),AA1(L),RZ2(K),AA2(K),			&
                B12(L),B12(K),DT1,DT2,DT3,DT4,DT5) 
		DTDA1=DTDA1+DT2
		DTDA2=DTDA2+DT1
		DTDA5=DTDA5+DT3
		DTDAZ=DTDAZ+DT5
		DTDBZ=DTDBZ+DT4

              CALL dELECREP(OV2,						&
				D1,D2,D3,D4,D5,							&
				RZ1(L),AA1(L),RZ2(K),AA2(K),       		&
                RZ2(L),AA2(L),RZ1(K),AA1(K),			&
                B12(L),B12(K),DER1,DER2,DER3,DER4,DER5)  
		DERDA1=DERDA1+DER1
		DERDA2=DERDA2+DER2
		DERDA5=DERDA5+DER3
		DERDAZ=DERDAZ+DER4
		DERDBZ=DERDBZ+DER5

              CALL dENUCLEAR(OV2,						&
				D1,D2,D3,D4,D5,							&
				RZ1(L),AA1(L),RZ2(K),AA2(K),			&
                RZ2(L),AA2(L),RZ1(K),AA1(K),			&
                B12(L),B12(K),EZ1,-ONE,					&
				DEN1,DEN2,DEN3,DEN4,DEN5)                       
		DENDA1=DENDA1+DEN1
		DENDA2=DENDA2+DEN2
		DENDA5=DENDA5+DEN3
		DENDAZ=DENDAZ+DEN4
		DENDBZ=DENDBZ+DEN5

              CALL dENUCLEAR(OV2,						&
				D1,D2,D3,D4,D5,							&
				RZ1(L),AA1(L),RZ2(K),AA2(K),		    &
                RZ2(L),AA2(L),RZ1(K),AA1(K),			&
                B12(L),B12(K),EZ2,ONE,					&
				DEN1,DEN2,DEN3,DEN4,DEN5)                       
		DENDA1=DENDA1+DEN1
		DENDA2=DENDA2+DEN2
		DENDA5=DENDA5+DEN3
		DENDAZ=DENDAZ+DEN4
		DENDBZ=DENDBZ+DEN5

              CALL dENUCLEAR(OV2,						&
				D2,D1,D3,D5,D4,							&
				RZ2(L),AA2(L),RZ1(K),AA1(K),			&
                RZ1(L),AA1(L),RZ2(K),AA2(K),			&	
                B12(L),B12(K),EZ1,-ONE,					&
				DEN1,DEN2,DEN3,DEN4,DEN5)                       
		DENDA1=DENDA1+DEN2
		DENDA2=DENDA2+DEN1
		DENDA5=DENDA5+DEN3
		DENDAZ=DENDAZ+DEN5
		DENDBZ=DENDBZ+DEN4

              CALL dENUCLEAR(OV2,						&
				D2,D1,D3,D5,D4,							&
				RZ2(L),AA2(L),RZ1(K),AA1(K),			&
                RZ1(L),AA1(L),RZ2(K),AA2(K),			&
                B12(L),B12(K),EZ2,ONE,					&
				DEN1,DEN2,DEN3,DEN4,DEN5)  
		DENDA1=DENDA1+DEN2
		DENDA2=DENDA2+DEN1
		DENDA5=DENDA5+DEN3
		DENDAZ=DENDAZ+DEN5
		DENDBZ=DENDBZ+DEN4

!-----------------------------------------------------
		CALL dOVERLAP(RZ1(L),AA1(L),-RZ1(K),AA1(K),     &             
				RZ2(L),AA2(L),-RZ2(K),AA2(K),B12(L),	&
				B12(K),D1,D2,D3,D4,D5,OV3)
			DSDA1=DSDA1+D1
			DSDA2=DSDA2+D2
			DSDA5=DSDA5+D3
			DSDAZ=DSDAZ+D4
			DSDBZ=DSDBZ+D5

              CALL DTKINET(OV3,							&
				D1,D2,D3,D4,D5,							&
				RZ1(L),AA1(L),-RZ1(K),AA1(K),			&                      
		        RZ2(L),AA2(L),-RZ2(K),AA2(K),			&
                B12(L),B12(K),DT1,DT2,DT3,DT4,DT5)
		DTDA1=DTDA1+DT1
		DTDA2=DTDA2+DT2
		DTDA5=DTDA5+DT3
		DTDAZ=DTDAZ+DT4
		DTDBZ=DTDBZ+DT5

              CALL DTKINET(OV3,							&
				D2,D1,D3,D5,D4,							&
				RZ2(L),AA2(L),-RZ2(K),AA2(K),			&
				RZ1(L),AA1(L),-RZ1(K),AA1(K),			&
				B12(L),B12(K),DT1,DT2,DT3,DT4,DT5)
		DTDA1=DTDA1+DT2
		DTDA2=DTDA2+DT1
		DTDA5=DTDA5+DT3
		DTDAZ=DTDAZ+DT5
		DTDBZ=DTDBZ+DT4

              CALL dELECREP(OV3,						&
				D1,D2,D3,D4,D5,							&
				RZ1(L),AA1(L),-RZ1(K),AA1(K),		    &
                RZ2(L),AA2(L),-RZ2(K),AA2(K),			&
                B12(L),B12(K),DER1,DER2,DER3,DER4,DER5)
		DERDA1=DERDA1+DER1
		DERDA2=DERDA2+DER2
		DERDA5=DERDA5+DER3
		DERDAZ=DERDAZ+DER4
		DERDBZ=DERDBZ+DER5

              CALL dENUCLEAR(OV3,						&
				D1,D2,D3,D4,D5,							&
				RZ1(L),AA1(L),-RZ1(K),AA1(K),			&
                RZ2(L),AA2(L),-RZ2(K),AA2(K),			&
                B12(L),B12(K),EZ1,-ONE,					&
				DEN1,DEN2,DEN3,DEN4,DEN5)
		DENDA1=DENDA1+DEN1
		DENDA2=DENDA2+DEN2
		DENDA5=DENDA5+DEN3
		DENDAZ=DENDAZ+DEN4
		DENDBZ=DENDBZ+DEN5

             CALL dENUCLEAR(OV3,						&
				D1,D2,D3,D4,D5,							&
				RZ1(L),AA1(L),-RZ1(K),AA1(K),			&
                RZ2(L),AA2(L),-RZ2(K),AA2(K),			&
                B12(L),B12(K),EZ2,ONE,					&
				DEN1,DEN2,DEN3,DEN4,DEN5)
		DENDA1=DENDA1+DEN1
		DENDA2=DENDA2+DEN2
		DENDA5=DENDA5+DEN3
		DENDAZ=DENDAZ+DEN4
		DENDBZ=DENDBZ+DEN5

              CALL dENUCLEAR(OV3,						&
				D2,D1,D3,D5,D4,							&
				RZ2(L),AA2(L),-RZ2(K),AA2(K),			&
                RZ1(L),AA1(L),-RZ1(K),AA1(K),			&
                B12(L),B12(K),EZ1,-ONE,					&
				DEN1,DEN2,DEN3,DEN4,DEN5)
		DENDA1=DENDA1+DEN2
		DENDA2=DENDA2+DEN1
		DENDA5=DENDA5+DEN3
		DENDAZ=DENDAZ+DEN5
		DENDBZ=DENDBZ+DEN4

              CALL dENUCLEAR(OV3,						&
				D2,D1,D3,D5,D4,							&
				RZ2(L),AA2(L),-RZ2(K),AA2(K),			&
                RZ1(L),AA1(L),-RZ1(K),AA1(K),			&
                B12(L),B12(K),EZ2,ONE,					&
				DEN1,DEN2,DEN3,DEN4,DEN5)
		DENDA1=DENDA1+DEN2
		DENDA2=DENDA2+DEN1
		DENDA5=DENDA5+DEN3
		DENDAZ=DENDAZ+DEN5
		DENDBZ=DENDBZ+DEN4

!-----------------------------------------------------
		CALL dOVERLAP(RZ1(L),AA1(L),-RZ2(K),AA2(K),		&             
				RZ2(L),AA2(L),-RZ1(K),AA1(K),			&
				B12(L),B12(K),D1,D2,D3,D4,D5,OV4)
			DSDA1=DSDA1+D1
			DSDA2=DSDA2+D2
			DSDA5=DSDA5+D3
			DSDAZ=DSDAZ+D4
			DSDBZ=DSDBZ+D5

              CALL DTKINET(OV4,							&
				D1,D2,D3,D4,D5,							&
                RZ1(L),AA1(L),-RZ2(K),AA2(K),			&                      
                RZ2(L),AA2(L),-RZ1(K),AA1(K),			&
                B12(L),B12(K),DT1,DT2,DT3,DT4,DT5)
		DTDA1=DTDA1+DT1
		DTDA2=DTDA2+DT2
		DTDA5=DTDA5+DT3
		DTDAZ=DTDAZ+DT4
		DTDBZ=DTDBZ+DT5

              CALL DTKINET(OV4,							&
				D2,D1,D3,D5,D4,							&
                RZ2(L),AA2(L),-RZ1(K),AA1(K),			&
                RZ1(L),AA1(L),-RZ2(K),AA2(K),			&
                B12(L),B12(K),DT1,DT2,DT3,DT4,DT5)
		DTDA1=DTDA1+DT2
		DTDA2=DTDA2+DT1
		DTDA5=DTDA5+DT3
		DTDAZ=DTDAZ+DT5
		DTDBZ=DTDBZ+DT4

              CALL dELECREP(OV4,						&
				D1,D2,D3,D4,D5,							&
				RZ1(L),AA1(L),-RZ2(K),AA2(K),			&
                RZ2(L),AA2(L),-RZ1(K),AA1(K),			&
                B12(L),B12(K),DER1,DER2,DER3,DER4,DER5)
		DERDA1=DERDA1+DER1
		DERDA2=DERDA2+DER2
		DERDA5=DERDA5+DER3
		DERDAZ=DERDAZ+DER4
		DERDBZ=DERDBZ+DER5

              CALL dENUCLEAR(OV4,						&
				D1,D2,D3,D4,D5,							&
				RZ1(L),AA1(L),-RZ2(K),AA2(K),			&
                RZ2(L),AA2(L),-RZ1(K),AA1(K),			&
                B12(L),B12(K),EZ1,-ONE,					&
				DEN1,DEN2,DEN3,DEN4,DEN5)
		DENDA1=DENDA1+DEN1
		DENDA2=DENDA2+DEN2
		DENDA5=DENDA5+DEN3
		DENDAZ=DENDAZ+DEN4
		DENDBZ=DENDBZ+DEN5

              CALL dENUCLEAR(OV4,						&
				D1,D2,D3,D4,D5,							&
				RZ1(L),AA1(L),-RZ2(K),AA2(K),			&
                RZ2(L),AA2(L),-RZ1(K),AA1(K),			&
                B12(L),B12(K),EZ2,ONE,					&
				DEN1,DEN2,DEN3,DEN4,DEN5)
		DENDA1=DENDA1+DEN1
		DENDA2=DENDA2+DEN2
		DENDA5=DENDA5+DEN3
		DENDAZ=DENDAZ+DEN4
		DENDBZ=DENDBZ+DEN5

              CALL dENUCLEAR(OV4,						&
				D2,D1,D3,D5,D4,							&
				RZ2(L),AA2(L),-RZ1(K),AA1(K),			&
                RZ1(L),AA1(L),-RZ2(K),AA2(K),			&
                B12(L),B12(K),EZ1,-ONE,					&
				DEN1,DEN2,DEN3,DEN4,DEN5)
		DENDA1=DENDA1+DEN2
		DENDA2=DENDA2+DEN1
		DENDA5=DENDA5+DEN3
		DENDAZ=DENDAZ+DEN5
		DENDBZ=DENDBZ+DEN4

              CALL dENUCLEAR(OV4,						&
				D2,D1,D3,D5,D4,							&
				RZ2(L),AA2(L),-RZ1(K),AA1(K),			&
                RZ1(L),AA1(L),-RZ2(K),AA2(K),			&
                B12(L),B12(K),EZ2,ONE,					&
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

			DH((L-K)*5+J,K)  =DHDA1
			DH((L-K)*5+J+1,K)=DHDA2
			DH((L-K)*5+J+2,K)=DHDA5
			DH((L-K)*5+J+3,K)=DHDAZ
			DH((L-K)*5+J+4,K)=DHDBZ

			DS((L-K)*5+J,K)  =DSDA1
			DS((L-K)*5+J+1,K)=DSDA2
			DS((L-K)*5+J+2,K)=DSDA5
			DS((L-K)*5+J+3,K)=DSDAZ
			DS((L-K)*5+J+4,K)=DSDBZ

			J=J+5


		END DO
	END DO


END SUBROUTINE Calculate_Total_Gradient