!      PROJECT TITLE:  correlated gaussian geminals for molecular hydrogen
!MAIN EXECUTION FILE:  cgauss.F90
!
!   Program FILENAME:  cgauss.exe
!        Lead Author:  D. Gilmore
!    Other Author(s):  L. Adamowicz, D. Kinghorn
!
!        LAST EDITED:  23may96:wed/08:05
!

!NOTES FOR FUTURE EDITS:
!===================================================================================================
!	1.  unnecessary sharing of modules should be eliminated
!===================================================================================================
!
!

!-------------------------------------------------------&   !=======================================
SUBROUTINE Energy_Diagnosis(XAA1,XRZ1,XAA2,XRZ2,XB12,E)

	USE scalars_to_be_shared
	USE vectors_to_be_shared
	USE matrices_to_be_shared

    IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	IMPLICIT INTEGER(I-N)

	DOUBLE PRECISION ::  Snp1,Hnp1

	DOUBLE PRECISION HH(2,2),SS(2,2),U(2,2),BIG(2),JB(2),VH(3)

    MMAX=M
	L=1
	K=1

        CALL OVERLAP(XRZ1,XAA1,	&
				XRZ1,XAA1,		&
                XRZ2,XAA2,		&
                XRZ2,XAA2,		&
                XB12,XB12,OV1)

		CALL OVERLAP(XRZ1,XAA1,	&
                XRZ2,XAA2,        &
                XRZ2,XAA2,		&
                XRZ1,XAA1,		&
                XB12,XB12,OV2)

        CALL OVERLAP(XRZ1,XAA1,	&
                -XRZ1,XAA1,     &                 
                XRZ2,XAA2,		&
                -XRZ2,XAA2,		&
                XB12,XB12,OV3)

        CALL OVERLAP(XRZ1,XAA1,	&
                -XRZ2,XAA2,     &                 
                XRZ2,XAA2,		&
                -XRZ1,XAA1,		&
                XB12,XB12,OV4)

              Snp1=OV1+OV2+OV3+OV4

              CALL TKINET(OV1,				&
                 XRZ1,XAA1,		&
                 XRZ1,XAA1,      &                
                 XRZ2,XAA2,		&
                 XRZ2,XAA2,		&
                 XB12,XB12,TK1)  

              CALL TKINET(OV2,				&
                 XRZ1,XAA1,		&
                 XRZ2,XAA2,      &                
                 XRZ2,XAA2,		&
                 XRZ1,XAA1,		&
                 XB12,XB12,TK2)  

              CALL TKINET(OV3,				&
                 XRZ1,XAA1,		&
                 -XRZ1,XAA1,   &                   
                 XRZ2,XAA2,		&
                 -XRZ2,XAA2,	&
                 XB12,XB12,TK3)

              CALL TKINET(OV4,				&
                 XRZ1,XAA1,		&
                 -XRZ2,XAA2,   &                   
                 XRZ2,XAA2,		&
                 -XRZ1,XAA1,	&
                 XB12,XB12,TK4)

             CALL TKINET(OV1,				&
                 XRZ2,XAA2,		&
                 XRZ2,XAA2,      &                
                 XRZ1,XAA1,		&
                 XRZ1,XAA1,		&
                 XB12,XB12,TK5)  

              CALL TKINET(OV2,				&
                 XRZ2,XAA2,		&
                 XRZ1,XAA1,      &                
                 XRZ1,XAA1,		&
                 XRZ2,XAA2,		&
                 XB12,XB12,TK6) 

              CALL TKINET(OV3,				&
                 XRZ2,XAA2,		&
                 -XRZ2,XAA2,   &                   
                 XRZ1,XAA1,		&
                 -XRZ1,XAA1,	&
                 XB12,XB12,TK7)

              CALL TKINET(OV4,				&
                 XRZ2,XAA2,		&
                 -XRZ1,XAA1,   &                   
                 XRZ1,XAA1,		&
                 -XRZ2,XAA2,	&
                 XB12,XB12,TK8)				

            T=(TK1+TK2+TK3+TK4+TK5+TK6+TK7+TK8)/TWO

              CALL ENUCLEAR(OV1,XRZ1,XAA1,    &
                XRZ1,XAA1,                       &
                XRZ2,XAA2,					   &
                XRZ2,XAA2,					   &
                XB12,XB12,								   &
                EZ1,EN1)

              CALL ENUCLEAR(OV2,XRZ1,XAA1,	   &
                XRZ2,XAA2,                       &
                XRZ2,XAA2,					   &
                XRZ1,XAA1,					   &
                XB12,XB12,								   &
                EZ1,EN2)                       

              CALL ENUCLEAR(OV3,XRZ1,XAA1,	   &
                 -XRZ1,XAA1,                   &   
                 XRZ2,XAA2,					   &
                 -XRZ2,XAA2,				   &
                 XB12,XB12,								   &
                 EZ1,EN3)

              CALL ENUCLEAR(OV4,XRZ1,XAA1,	   &
                 -XRZ2,XAA2,                   &   
                 XRZ2,XAA2,					   &
                 -XRZ1,XAA1,				   &
                 XB12,XB12,								   &
                 EZ1,EN4)

              CALL ENUCLEAR(OV1,XRZ1,XAA1,	   &
                XRZ1,XAA1,                       &
                XRZ2,XAA2,					   &
                XRZ2,XAA2,					   &
                XB12,XB12,								   &
                EZ2,EN5)

              CALL ENUCLEAR(OV2,XRZ1,XAA1,	   &
                XRZ2,XAA2,                       &
                XRZ2,XAA2,					   &
                XRZ1,XAA1,					   &
                XB12,XB12,								   &
                EZ2,EN6)                       

              CALL ENUCLEAR(OV3,XRZ1,XAA1,	   &
                 -XRZ1,XAA1,                   &   
                 XRZ2,XAA2,					   &
                 -XRZ2,XAA2,				   &
                 XB12,XB12,								   &
                 EZ2,EN7)

              CALL ENUCLEAR(OV4,XRZ1,XAA1,	   &
                 -XRZ2,XAA2,                   &   
                 XRZ2,XAA2,					   &
                 -XRZ1,XAA1,				   &
                 XB12,XB12,								   &
                 EZ2,EN8)

              CALL ENUCLEAR(OV1,XRZ2,XAA2,	   &
                XRZ2,XAA2,                       &
                XRZ1,XAA1,					   &
                XRZ1,XAA1,					   &
                XB12,XB12,								   &
                EZ1,EN9)                       

              CALL ENUCLEAR(OV2,XRZ2,XAA2,	   &
                XRZ1,XAA1,                       &
                XRZ1,XAA1,					   &
                XRZ2,XAA2,					   &
                XB12,XB12,								   &
                EZ1,EN10)                       

              CALL ENUCLEAR(OV3,XRZ2,XAA2,	   &
                 -XRZ2,XAA2,                   &   
                 XRZ1,XAA1,					   &
                 -XRZ1,XAA1,				   &
                 XB12,XB12,								   &
                 EZ1,EN11)

              CALL ENUCLEAR(OV4,XRZ2,XAA2,	   &
                 -XRZ1,XAA1,                   &   
                 XRZ1,XAA1,					   &
                 -XRZ2,XAA2,				   &
                 XB12,XB12,								   &
                 EZ1,EN12)					   

              CALL ENUCLEAR(OV1,XRZ2,XAA2,	   &
                XRZ2,XAA2,                       &
                XRZ1,XAA1,					   &
                XRZ1,XAA1,					   &
                XB12,XB12,								   &
                EZ2,EN13)                       


              CALL ENUCLEAR(OV2,XRZ2,XAA2,	   &
                XRZ1,XAA1,                       &
                XRZ1,XAA1,					   &
                XRZ2,XAA2,					   &
                XB12,XB12,								   &
                EZ2,EN14)  

              CALL ENUCLEAR(OV3,XRZ2,XAA2,	   &
                 -XRZ2,XAA2,                   &   
                 XRZ1,XAA1,					   &
                 -XRZ1,XAA1,				   &
                 XB12,XB12,								   &
                 EZ2,EN15)

              CALL ENUCLEAR(OV4,XRZ2,XAA2,	   &
                 -XRZ1,XAA1,                   &   
                 XRZ1,XAA1,					   &
                 -XRZ2,XAA2,				   &
                 XB12,XB12,								   &
                 EZ2,EN16)

              CALL ELECREP(OV1,XRZ1,XAA1,		&
                XRZ1,XAA1,                    &  
                XRZ2,XAA2,					&
                XRZ2,XAA2,					&
                XB12,XB12,ER1)  

              CALL ELECREP(OV2,XRZ1,XAA1,		&
                XRZ2,XAA2,                    &  
                XRZ2,XAA2,					&
                XRZ1,XAA1,					&
                XB12,XB12,ER2)  

              CALL ELECREP(OV3,XRZ1,XAA1,		&
                 -XRZ1,XAA1,                &      
                 XRZ2,XAA2,					&
                 -XRZ2,XAA2,				&
                 XB12,XB12,ER3)

              CALL ELECREP(OV4,XRZ1,XAA1,		&
                 -XRZ2,XAA2,                &      
                 XRZ2,XAA2,					&
                 -XRZ1,XAA1,				&
                 XB12,XB12,ER4)

              EN=EN1+EN2+EN3+EN4+EN5+EN6+EN7+EN8+EN9+EN10+EN11+EN12+EN13+EN14+EN15+EN16
              ER=ER1+ER2+ER3+ER4

              Hnp1=ER+T-EN+(ONE/BD)*Snp1


	SS(1,1)=ONE
	SS(2,2)=Snp1/Snp1
    HH(2,2)=Hnp1/Snp1
	HH(1,2)=ZERO
	SS(1,2)=ZERO

!    IF (M.GE.2) K=2
	DO L=1,M
		CALL OVERLAP(RZ1(L),AA1(L),	&
				XRZ1,XAA1,                &      
                RZ2(L),AA2(L),		&
                XRZ2,XAA2,				&
                B12(L),XB12,OV1)

		CALL OVERLAP(RZ1(L),AA1(L),	&
                XRZ2,XAA2,                &      
                RZ2(L),AA2(L),		&
                XRZ1,XAA1,				&
                B12(L),XB12,OV2)

		CALL OVERLAP(RZ1(L),AA1(L),	&
                -XRZ1,XAA1,             &         
                RZ2(L),AA2(L),		&
                -XRZ2,XAA2,				&
                B12(L),XB12,OV3)

		CALL OVERLAP(RZ1(L),AA1(L),	&
                -XRZ2,XAA2,             &         
                RZ2(L),AA2(L),		&
                -XRZ1,XAA1,				&
                B12(L),XB12,OV4)					

		Stemp=OV1+OV2+OV3+OV4

              CALL TKINET(OV1,					  &
                RZ1(L),AA1(L),	  &
                XRZ1,XAA1,            &          
                RZ2(L),AA2(L),	  &
                XRZ2,XAA2,			  &
                 B12(L),XB12,TK1)  

              CALL TKINET(OV2,					  &
                RZ1(L),AA1(L),	  &
                XRZ2,XAA2,            &          
                RZ2(L),AA2(L),	  &
                XRZ1,XAA1,			  &
                 B12(L),XB12,TK2)  

              CALL TKINET(OV3,					  &
                RZ1(L),AA1(L),	  &
                -XRZ1,XAA1,         &             
                RZ2(L),AA2(L),	  &
                -XRZ2,XAA2,		  &
                 B12(L),XB12,TK3)

              CALL TKINET(OV4,					  &
                RZ1(L),AA1(L),	  &
                -XRZ2,XAA2,         &             
                RZ2(L),AA2(L),	  &
                -XRZ1,XAA1,		  &
                 B12(L),XB12,TK4)

             CALL TKINET(OV1,					  &
                RZ2(L),AA2(L),	  &
                XRZ2,XAA2,            &          
                RZ1(L),AA1(L),	  &
                XRZ1,XAA1,			  &
                 B12(L),XB12,TK5)  

              CALL TKINET(OV2,					  &
                RZ2(L),AA2(L),	  &
                XRZ1,XAA1,            &          
                RZ1(L),AA1(L),	  &
                XRZ2,XAA2,			  &
                 B12(L),XB12,TK6) 

              CALL TKINET(OV3,					  &
                RZ2(L),AA2(L),	  &
                -XRZ2,XAA2,         &             
                RZ1(L),AA1(L),	  &
                -XRZ1,XAA1,		  &
                 B12(L),XB12,TK7)

              CALL TKINET(OV4,					  &
                RZ2(L),AA2(L),	  &
                -XRZ1,XAA1,         &             
                RZ1(L),AA1(L),	  &
                -XRZ2,XAA2,		  &
                 B12(L),XB12,TK8)

            T=(TK1+TK2+TK3+TK4+TK5+TK6+TK7+TK8)/TWO

              CALL ENUCLEAR(OV1,RZ1(L),AA1(L),	&
                XRZ1,XAA1,                      		&
                RZ2(L),AA2(L),					&
                XRZ2,XAA2,							&
                B12(L),XB12,									&
                EZ1,EN1)

              CALL ENUCLEAR(OV2,RZ1(L),AA1(L),	&
                XRZ2,XAA2,                      		&
                RZ2(L),AA2(L),					&
                XRZ1,XAA1,							&
                B12(L),XB12,									&
                EZ1,EN2)                       

              CALL ENUCLEAR(OV3,RZ1(L),AA1(L),	&
                 -XRZ1,XAA1,                      	&
                 RZ2(L),AA2(L),					&
                 -XRZ2,XAA2,						&
                 B12(L),XB12,									&
                 EZ1,EN3)

              CALL ENUCLEAR(OV4,RZ1(L),AA1(L),	&
                 -XRZ2,XAA2,                      	&
                 RZ2(L),AA2(L),					&
                 -XRZ1,XAA1,						&
                 B12(L),XB12,									&
                 EZ1,EN4)

              CALL ENUCLEAR(OV1,RZ1(L),AA1(L),	&
                XRZ1,XAA1,                      		&
                RZ2(L),AA2(L),					&
                XRZ2,XAA2,							&
                B12(L),XB12,									&
                EZ2,EN5)

              CALL ENUCLEAR(OV2,RZ1(L),AA1(L),	&
                XRZ2,XAA2,                      		&
                RZ2(L),AA2(L),					&
                XRZ1,XAA1,							&
                B12(L),XB12,									&
                EZ2,EN6)                       

              CALL ENUCLEAR(OV3,RZ1(L),AA1(L),	&
                 -XRZ1,XAA1,                      	&
                 RZ2(L),AA2(L),					&
                 -XRZ2,XAA2,						&
                 B12(L),XB12,									&
                 EZ2,EN7)

              CALL ENUCLEAR(OV4,RZ1(L),AA1(L),	&
                 -XRZ2,XAA2,                      	&
                 RZ2(L),AA2(L),					&
                 -XRZ1,XAA1,						&
                 B12(L),XB12,									&
                 EZ2,EN8)

              CALL ENUCLEAR(OV1,RZ2(L),AA2(L),	&
                XRZ2,XAA2,                      		&
                RZ1(L),AA1(L),					&
                XRZ1,XAA1,							&
                B12(L),XB12,									&
                EZ1,EN9)                       

              CALL ENUCLEAR(OV2,RZ2(L),AA2(L),	&
                XRZ1,XAA1,                      		&
                RZ1(L),AA1(L),					&
                XRZ2,XAA2,							&
                B12(L),XB12,									&
                EZ1,EN10)                       

              CALL ENUCLEAR(OV3,RZ2(L),AA2(L),	&
                 -XRZ2,XAA2,                      	&
                 RZ1(L),AA1(L),					&
                 -XRZ1,XAA1,						&
                 B12(L),XB12,									&
                 EZ1,EN11)

              CALL ENUCLEAR(OV4,RZ2(L),AA2(L),	&
                 -XRZ1,XAA1,                      	&
                 RZ1(L),AA1(L),					&
                 -XRZ2,XAA2,						&
                 B12(L),XB12,									&
                 EZ1,EN12)

              CALL ENUCLEAR(OV1,RZ2(L),AA2(L),	&
                XRZ2,XAA2,                      		&
                RZ1(L),AA1(L),					&
                XRZ1,XAA1,							&
                B12(L),XB12,									&
                EZ2,EN13)                       

              CALL ENUCLEAR(OV2,RZ2(L),AA2(L),	&
                XRZ1,XAA1,                      		&
                RZ1(L),AA1(L),					&
                XRZ2,XAA2,							&
                B12(L),XB12,									&
                EZ2,EN14)  

              CALL ENUCLEAR(OV3,RZ2(L),AA2(L),	&
                 -XRZ2,XAA2,                      	&
                 RZ1(L),AA1(L),					&
                 -XRZ1,XAA1,						&
                 B12(L),XB12,									&
                 EZ2,EN15)
	 
	               CALL ENUCLEAR(OV4,RZ2(L),AA2(L), &
                 -XRZ1,XAA1,                      	  &
                 RZ1(L),AA1(L),					  &
                 -XRZ2,XAA2,						  &
                 B12(L),XB12,									  &
                 EZ2,EN16)

              CALL ELECREP(OV1,RZ1(L),AA1(L),   &
                XRZ1,XAA1,                      	  &
                RZ2(L),AA2(L),				  &
                XRZ2,XAA2,						  &
                B12(L),XB12,ER1)  

              CALL ELECREP(OV2,RZ1(L),AA1(L),	  &
                XRZ2,XAA2,                      	  &
                RZ2(L),AA2(L),				  &
                XRZ1,XAA1,						  &
                B12(L),XB12,ER2)  

              CALL ELECREP(OV3,RZ1(L),AA1(L),	  &
                 -XRZ1,XAA1,                      &
                 RZ2(L),AA2(L),				  &
                 -XRZ2,XAA2,					  &
                 B12(L),XB12,ER3)

              CALL ELECREP(OV4,RZ1(L),AA1(L),	  &
                 -XRZ2,XAA2,                      &
                 RZ2(L),AA2(L),				  &
                 -XRZ1,XAA1,					  &
                 B12(L),XB12,ER4)

              EN=EN1+EN2+EN3+EN4+EN5+EN6+EN7+EN8+EN9+EN10+EN11+EN12+EN13+EN14+EN15+EN16

              ER=ER1+ER2+ER3+ER4

    HH(1,2)=HH(1,2)+ C(L)*(ER+T-EN+(ONE/BD)*Stemp)
	SS(1,2)=SS(1,2)+ (C(L)*Stemp)

      END DO

	Stemp=ZERO
	HH(1,1)=ZERO
      DO L=1,M
		DO K=1,M
			HH(1,1)=HH(1,1)+(C(L)*C(K)*H(L,K))
			Stemp=Stemp + (C(L)*C(K)*S(L,K))
		END DO
	END DO
	HH(1,1)=HH(1,1)/Stemp
	SS(2,2)=Snp1/Snp1

	HH(1,2)=HH(1,2)/DSQRT(Stemp*Snp1)
	HH(2,1)=HH(1,2)
	SS(1,2)=SS(1,2)/DSQRT(Stemp*Snp1)
	SS(2,1)=SS(1,2)
	VH(1)=HH(1,1)
	VH(2)=HH(2,1)
	VH(3)=HH(2,2)				 

!      CALL ENERGY(2,HH,SS,EMIN,CO,KIM)
	CALL Solve_2_by_2(2,HH,SS,EMIN)
E=EMIN

END SUBROUTINE Energy_Diagnosis