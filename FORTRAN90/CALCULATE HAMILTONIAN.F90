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
SUBROUTINE Calculate_Hamiltonian

	USE scalars_to_be_shared			  
	USE	vectors_to_be_shared			  
	USE	matrices_to_be_shared			  

	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	IMPLICIT          INTEGER    (I-N)

	DO L=1,M
		DO K=1,L

			CALL OVERLAP(RZ1(L),AA1(L),RZ1(K),AA1(K),	&                  
   					RZ2(L),AA2(L),RZ2(K),AA2(K),B12(L),	&
					B12(K),OV1)

            CALL ENUCLEAR(OV1,RZ1(L),AA1(L),RZ1(K),		&
					AA1(K),RZ2(L),AA2(L),RZ2(K),AA2(K),	&
					B12(L),B12(K),EZ1,EN1)
            CALL ENUCLEAR(OV1,RZ1(L),AA1(L),RZ1(K),		&
					AA1(K),RZ2(L),AA2(L),RZ2(K),AA2(K),	&
					B12(L),B12(K),EZ2,EN5)
            CALL ENUCLEAR(OV1,RZ2(L),AA2(L),RZ2(K),		&
					AA2(K),RZ1(L),AA1(L),RZ1(K),AA1(K),	&
					B12(L),B12(K),EZ1,EN9)                       
			CALL ENUCLEAR(OV1,RZ2(L),AA2(L),RZ2(K),		&
					AA2(K),RZ1(L),AA1(L),RZ1(K),AA1(K),	&
	                B12(L),B12(K),EZ2,EN13)                       

            CALL ELECREP(OV1,RZ1(L),AA1(L),RZ1(K),		&
					AA1(K),RZ2(L),AA2(L),RZ2(K),AA2(K),	&
	                B12(L),B12(K),ER1)  

            CALL TKINET(OV1,RZ1(L),AA1(L),RZ1(K),		&
					AA1(K),RZ2(L),AA2(L),RZ2(K),AA2(K),	&
					B12(L),B12(K),TK1)  
            CALL TKINET(OV1,RZ2(L),AA2(L),RZ2(K),		&
					AA2(K),RZ1(L),AA1(L),RZ1(K),AA1(K),	&
	                B12(L),B12(K),TK5)  





			CALL OVERLAP(RZ1(L),AA1(L),RZ2(K),AA2(K),	&
					RZ2(L),AA2(L),RZ1(K),AA1(K),B12(L),	&
					B12(K),OV2)

            CALL ENUCLEAR(OV2,RZ1(L),AA1(L),RZ2(K),		&
					AA2(K),RZ2(L),AA2(L),RZ1(K),AA1(K),	&
					B12(L),B12(K),EZ1,EN2)                       
            CALL ENUCLEAR(OV2,RZ1(L),AA1(L),RZ2(K),		&
					AA2(K),RZ2(L),AA2(L),RZ1(K),AA1(K),	&
	                B12(L),B12(K),EZ2,EN6)                       
            CALL ENUCLEAR(OV2,RZ2(L),AA2(L),RZ1(K),		&
					AA1(K),RZ1(L),AA1(L),RZ2(K),AA2(K),	&
	                B12(L),B12(K),EZ1,EN10)                       
            CALL ENUCLEAR(OV2,RZ2(L),AA2(L),RZ1(K),		&
					AA1(K),RZ1(L),AA1(L),RZ2(K),AA2(K),	&
	                B12(L),B12(K),EZ2,EN14)  

            CALL ELECREP(OV2,RZ1(L),AA1(L),RZ2(K),		&
					AA2(K),RZ2(L),AA2(L),RZ1(K),AA1(K),	&
	                B12(L),B12(K),ER2)  

            CALL TKINET(OV2,RZ1(L),AA1(L),RZ2(K),AA2(K),&                   
					RZ2(L),AA2(L),RZ1(K),AA1(K),B12(L),	&
					B12(K),TK2)  
            CALL TKINET(OV2,RZ2(L),AA2(L),RZ1(K),AA1(K),&                   
					RZ1(L),AA1(L),RZ2(K),AA2(K),B12(L),	&
					B12(K),TK6) 




			CALL OVERLAP(RZ1(L),AA1(L),-RZ1(K),AA1(K),  &             
					RZ2(L),AA2(L),-RZ2(K),AA2(K),B12(L),&
					B12(K),OV3)

            CALL ENUCLEAR(OV3,RZ1(L),AA1(L),-RZ1(K),	&
					AA1(K),RZ2(L),AA2(L),-RZ2(K),AA2(K),&
					B12(L),B12(K),EZ1,EN3)
            CALL ENUCLEAR(OV3,RZ1(L),AA1(L),-RZ1(K),	&
					AA1(K),RZ2(L),AA2(L),-RZ2(K),AA2(K),&
					B12(L),B12(K),EZ2,EN7)
            CALL ENUCLEAR(OV3,RZ2(L),AA2(L),-RZ2(K),	&
					AA2(K),RZ1(L),AA1(L),-RZ1(K),AA1(K),&
					B12(L),B12(K),EZ1,EN11)
            CALL ENUCLEAR(OV3,RZ2(L),AA2(L),-RZ2(K),	&
					AA2(K),RZ1(L),AA1(L),-RZ1(K),AA1(K),&
					B12(L),B12(K),EZ2,EN15)

            CALL ELECREP(OV3,RZ1(L),AA1(L),-RZ1(K),		&
					AA1(K),RZ2(L),AA2(L),-RZ2(K),AA2(K),&
					B12(L),B12(K),ER3)

            CALL TKINET(OV3,RZ1(L),AA1(L),-RZ1(K),		&
					AA1(K),RZ2(L),AA2(L),-RZ2(K),AA2(K),&
					B12(L),B12(K),TK3)
            CALL TKINET(OV3,RZ2(L),AA2(L),-RZ2(K),		&
					AA2(K),RZ1(L),AA1(L),-RZ1(K),AA1(K),&
					B12(L),B12(K),TK7)



			CALL OVERLAP(RZ1(L),AA1(L),-RZ2(K),AA2(K),  &             
					RZ2(L),AA2(L),-RZ1(K),AA1(K),B12(L),&
					B12(K),OV4)

            CALL ENUCLEAR(OV4,RZ1(L),AA1(L),-RZ2(K),	&
					AA2(K),RZ2(L),AA2(L),-RZ1(K),AA1(K),&
					B12(L),B12(K),EZ1,EN4)
            CALL ENUCLEAR(OV4,RZ1(L),AA1(L),-RZ2(K),	&
					AA2(K),RZ2(L),AA2(L),-RZ1(K),AA1(K),&
					B12(L),B12(K),EZ2,EN8)
            CALL ENUCLEAR(OV4,RZ2(L),AA2(L),-RZ1(K),	&
					AA1(K),RZ1(L),AA1(L),-RZ2(K),AA2(K),&
					B12(L),B12(K),EZ1,EN12)
            CALL ENUCLEAR(OV4,RZ2(L),AA2(L),-RZ1(K),	&
					AA1(K),RZ1(L),AA1(L),-RZ2(K),AA2(K),&
					B12(L),B12(K),EZ2,EN16)

            CALL ELECREP(OV4,RZ1(L),AA1(L),-RZ2(K),		&
					AA2(K),RZ2(L),AA2(L),-RZ1(K),AA1(K),&
					B12(L),B12(K),ER4)

            CALL TKINET(OV4,RZ1(L),AA1(L),-RZ2(K),		&
					AA2(K),RZ2(L),AA2(L),-RZ1(K),AA1(K),&
					B12(L),B12(K),TK4)
            CALL TKINET(OV4,RZ2(L),AA2(L),-RZ1(K),		&
					AA1(K),RZ1(L),AA1(L),-RZ2(K),AA2(K),&
					B12(L),B12(K),TK8)



            OV=OV1+OV2+OV3+OV4

			EN=EN1+EN2+EN3+EN4+EN5+EN6+EN7+EN8+EN9+EN10+EN11+EN12+EN13+EN14+EN15+EN16

			ER=ER1+ER2+ER3+ER4

			TK=HALF*(TK1+TK2+TK3+TK4+TK5+TK6+TK7+TK8)


		H(L,K)=TK+ER-EN+(ONE/BD)*OV
		S(L,K)=OV

		END DO
	END DO


END SUBROUTINE Calculate_Hamiltonian