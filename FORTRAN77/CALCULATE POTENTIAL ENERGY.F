SUBROUTINE Calculate_Potential_Energy

	USE scalars_to_be_shared			  
	USE	vectors_to_be_shared			  
	USE	matrices_to_be_shared			  

	INTEGER :: L,K

	DOUBLE PRECISION EN1,EN2,EN3,EN4,EN5,EN6,EN7,EN8,			&
					 EN9,EN10,EN11,EN12,EN13,EN14,EN15,EN16

	DOUBLE PRECISION ER1,ER2,ER3,ER4
	
	DOUBLE PRECISION EN,ER

	DO L=1,M
		DO K=1,L


!--- Calculate Nuc./Elec. Attractn. FOR ELECTRON 1 ---
!-----------------------------------------------------

              CALL ENUCLEAR(OVV1(L,K),RX1(L),RY1(L),RZ1(L),AA1(L),		&
               RX1(K),RY1(K),RZ1(K),AA1(K),                      		&
                RX2(L),RY2(L),RZ2(L),AA2(L),							&
                RX2(K),RY2(K),RZ2(K),AA2(K),							&
                B12(L),B12(K),											&
                EX1,EY1,EZ1,EN1)

!...Electron Symmetrization Integral: P_[elec.1,elec.2]

              CALL ENUCLEAR(OVV2(L,K),RX1(L),RY1(L),RZ1(L),AA1(L),		&
                RX2(K),RY2(K),RZ2(K),AA2(K),                      		&
                RX2(L),RY2(L),RZ2(L),AA2(L),							&
                RX1(K),RY1(K),RZ1(K),AA1(K),							&
                B12(L),B12(K),											&
                EX1,EY1,EZ1,EN2)                       

!...Nuclear Symmetrization Integral: P_[H1,H2]

              CALL ENUCLEAR(OVV3(L,K),RX1(L),RY1(L),RZ1(L),AA1(L),		&
                 -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      	&
                 RX2(L),RY2(L),RZ2(L),AA2(L),							&
                 -RX2(K),-RY2(K),-RZ2(K),AA2(K),						&
                 B12(L),B12(K),											&
                 EX1,EY1,EZ1,EN3)

!...Electro/Nuclear Symmetrization Integral:
!         P_[elec.1,elec.2]*P_[H1,H2]

              CALL ENUCLEAR(OVV4(L,K),RX1(L),RY1(L),RZ1(L),AA1(L),		&
                 -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      	&
                 RX2(L),RY2(L),RZ2(L),AA2(L),							&
                 -RX1(K),-RY1(K),-RZ1(K),AA1(K),						&
                 B12(L),B12(K),											&
                 EX1,EY1,EZ1,EN4)

  
              CALL ENUCLEAR(OVV1(L,K),RX1(L),RY1(L),RZ1(L),AA1(L),		&
                RX1(K),RY1(K),RZ1(K),AA1(K),                      		&
                RX2(L),RY2(L),RZ2(L),AA2(L),							&
                RX2(K),RY2(K),RZ2(K),AA2(K),							&
                B12(L),B12(K),											&
                EX2,EY2,EZ2,EN5)

!...Electron Symmetrization Integral: P_[elec.1,elec.2]

              CALL ENUCLEAR(OVV2(L,K),RX1(L),RY1(L),RZ1(L),AA1(L),		&
                RX2(K),RY2(K),RZ2(K),AA2(K),                      		&
                RX2(L),RY2(L),RZ2(L),AA2(L),							&
                RX1(K),RY1(K),RZ1(K),AA1(K),							&
                B12(L),B12(K),											&
                EX2,EY2,EZ2,EN6)                       

!...Nuclear Symmetrization Integral: P_[H1,H2]

             CALL ENUCLEAR(OVV3(L,K),RX1(L),RY1(L),RZ1(L),AA1(L),		&
                -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      	&
                RX2(L),RY2(L),RZ2(L),AA2(L),							&
                -RX2(K),-RY2(K),-RZ2(K),AA2(K),							&
                B12(L),B12(K),											&
                EX2,EY2,EZ2,EN7)

!...Electro/Nuclear Symmetrization Integral:
!         P_[elec.1,elec.2]*P_[H1,H2]

              CALL ENUCLEAR(OVV4(L,K),RX1(L),RY1(L),RZ1(L),AA1(L),		&
                 -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      	&
                 RX2(L),RY2(L),RZ2(L),AA2(L),							&
                 -RX1(K),-RY1(K),-RZ1(K),AA1(K),						&
                 B12(L),B12(K),											&
                 EX2,EY2,EZ2,EN8)

!--- Calculate Nuc./Elec. Attractn. FOR ELECTRON 2 ---
!-----------------------------------------------------

              CALL ENUCLEAR(OVV1(L,K),RX2(L),RY2(L),RZ2(L),AA2(L),		&
                RX2(K),RY2(K),RZ2(K),AA2(K),                      		&
                RX1(L),RY1(L),RZ1(L),AA1(L),							&
                RX1(K),RY1(K),RZ1(K),AA1(K),							&
                B12(L),B12(K),											&
                EX1,EY1,EZ1,EN9)                       

!...Electron Symmetrization Integral: P_[elec.1,elec.2]


              CALL ENUCLEAR(OVV2(L,K),RX2(L),RY2(L),RZ2(L),AA2(L),		&
                RX1(K),RY1(K),RZ1(K),AA1(K),                      		&
                RX1(L),RY1(L),RZ1(L),AA1(L),							&
                RX2(K),RY2(K),RZ2(K),AA2(K),							&
                B12(L),B12(K),											&
                EX1,EY1,EZ1,EN10)                       

!...Nuclear Symmetrization Integral: P_[H1,H2]

              CALL ENUCLEAR(OVV3(L,K),RX2(L),RY2(L),RZ2(L),AA2(L),		&
                 -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      	&
                 RX1(L),RY1(L),RZ1(L),AA1(L),							&
                 -RX1(K),-RY1(K),-RZ1(K),AA1(K),						&
                 B12(L),B12(K),											&
                 EX1,EY1,EZ1,EN11)

!...Electro/Nuclear Symmetrization Integral:
!         P_[elec.1,elec.2]*P_[H1,H2]

              CALL ENUCLEAR(OVV4(L,K),RX2(L),RY2(L),RZ2(L),AA2(L),		&
                 -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      	&
                 RX1(L),RY1(L),RZ1(L),AA1(L),							&
                 -RX2(K),-RY2(K),-RZ2(K),AA2(K),						&
                 B12(L),B12(K),											&
                 EX1,EY1,EZ1,EN12)

 
              CALL ENUCLEAR(OVV1(L,K),RX2(L),RY2(L),RZ2(L),AA2(L),		&
                RX2(K),RY2(K),RZ2(K),AA2(K),                      		&
                RX1(L),RY1(L),RZ1(L),AA1(L),							&
                RX1(K),RY1(K),RZ1(K),AA1(K),							&
                B12(L),B12(K),											&
                EX2,EY2,EZ2,EN13)                       

!...Electron Symmetrization Integral: P_[elec.1,elec.2]


              CALL ENUCLEAR(OVV2(L,K),RX2(L),RY2(L),RZ2(L),AA2(L),		&
                RX1(K),RY1(K),RZ1(K),AA1(K),                      		&
                RX1(L),RY1(L),RZ1(L),AA1(L),							&
                RX2(K),RY2(K),RZ2(K),AA2(K),							&
                B12(L),B12(K),											&
                EX2,EY2,EZ2,EN14)  

!...Nuclear Symmetrization Integral: P_[H1,H2]

              CALL ENUCLEAR(OVV3(L,K),RX2(L),RY2(L),RZ2(L),AA2(L),		&
                 -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      	&
                 RX1(L),RY1(L),RZ1(L),AA1(L),							&
                 -RX1(K),-RY1(K),-RZ1(K),AA1(K),						&
                 B12(L),B12(K),											&
                 EX2,EY2,EZ2,EN15)

!...Electro/Nuclear Symmetrization Integral:
!         P_[elec.1,elec.2]*P_[H1,H2]

              CALL ENUCLEAR(OVV4(L,K),RX2(L),RY2(L),RZ2(L),AA2(L),		&
                 -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      	&
                 RX1(L),RY1(L),RZ1(L),AA1(L),							&
                 -RX2(K),-RY2(K),-RZ2(K),AA2(K),						&
                 B12(L),B12(K),											&
                 EX2,EY2,EZ2,EN16)



!--- Calculate the Elec./Elec. Repulsion ---
!-----------------------------------------------------

              CALL ELECREP(OVV1(L,K),RX1(L),RY1(L),RZ1(L),AA1(L),		&
                RX1(K),RY1(K),RZ1(K),AA1(K),                      		&
                RX2(L),RY2(L),RZ2(L),AA2(L),							&
                RX2(K),RY2(K),RZ2(K),AA2(K),							&
                B12(L),B12(K),ER1)  

!...Electron Symmetrization Integral: P_[elec.1,elec.2]


              CALL ELECREP(OVV2(L,K),RX1(L),RY1(L),RZ1(L),AA1(L),		&
                RX2(K),RY2(K),RZ2(K),AA2(K),                      		&
                RX2(L),RY2(L),RZ2(L),AA2(L),							&
                RX1(K),RY1(K),RZ1(K),AA1(K),							&
                B12(L),B12(K),ER2)  

!...Nuclear Symmetrization Overlap Integral: P_[H1,H2]

              CALL ELECREP(OVV3(L,K),RX1(L),RY1(L),RZ1(L),AA1(L),		&
                 -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      	&
                 RX2(L),RY2(L),RZ2(L),AA2(L),							&
                -RX2(K),-RY2(K),-RZ2(K),AA2(K),							&
                 B12(L),B12(K),ER3)

!...Electro/Nuclear Symmetrization Overlap Integral:
!         P_[elec.1,elec.2]*P_[H1,H2]

              CALL ELECREP(OVV4(L,K),RX1(L),RY1(L),RZ1(L),AA1(L),		&
                 -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      	&
                 RX2(L),RY2(L),RZ2(L),AA2(L),							&
                 -RX1(K),-RY1(K),-RZ1(K),AA1(K),						&
                 B12(L),B12(K),ER4)





				EN=EN1+EN2+EN3+EN4+EN5+EN6+EN7+EN8		&
					+EN9+EN10+EN11+EN12+EN13+EN14+EN15+EN16


				ER=ER1+ER2+ER3+ER4

		Potential(L,K)=ER-EN+(1.00D+00/BD)*S(L,K)

		IF (L .NE. K) Potential(K,L)=Potential(L,K)


		END DO
	END DO


END SUBROUTINE Calculate_Potential_Energy
