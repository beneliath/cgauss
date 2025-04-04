SUBROUTINE Calculate_Kinetic_Energy

	USE scalars_to_be_shared			  
	USE	vectors_to_be_shared			  
	USE	matrices_to_be_shared			  

	DOUBLE PRECISION TK1,TK2,TK3,TK4,TK5,TK6,TK7,TK8

	INTEGER :: L,K


	DO L=1,M
		DO K=1,L


!--- Calculate KINETIC ENERGY OF ELECTRON 1 ---
!-----------------------------------------------------


!...Electron Symmetrization Integral: P_[elec.1,elec.2]


!...Nuclear Symmetrization Integral: P_[H1,H2]


!...Electro/Nuclear Symmetrization Integral:
!         P_[elec.1,elec.2]*P_[H1,H2]


!--- Calculate KINETIC ENERGY OF ELECTRON 2 ---
!-----------------------------------------------------



!...Electron Symmetrization Integral: P_[elec.1,elec.2]


!...Nuclear Symmetrization Integral: P_[H1,H2]


!...Electro/Nuclear Symmetrization Integral:
!         P_[elec.1,elec.2]*P_[H1,H2]


	

	Kinetic(L,K)=(TK1+TK2+TK3+TK4+TK5+TK6+TK7+TK8)/TWO
					
	IF (L .NE. K) Kinetic(K,L)=Kinetic(L,K)

		END DO
	END DO


END SUBROUTINE Calculate_Kinetic_Energy
