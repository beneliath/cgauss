SUBROUTINE GRADEN(N,X,G)

	USE scalars_to_be_shared
	USE vectors_to_be_shared
	USE matrices_to_be_shared
	USE strings_to_be_shared


    IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	IMPLICIT INTEGER(I-N)

	DIMENSION X(N),G(N)


	K=1
	DO I=1,M
			T1(I)=X(K)	
			T2(I)=X(K+1)	
			T3(I)=X(K+2)	
			RZ1(I)=X(K+3)	
			RZ2(I)=X(K+4)	
		K=K+5
	END DO

	DO I=1,M
			C(I)=X(M*5+I)
	END DO


	CALL Cholesky_Condition_Parameters
!	CALL Construct_S_and_H_Matrices
	CALL Construct_Gradient

!	CALL Calculate_Rayleigh_Quotient



	DO I=1,N
		G(I)=GRAD(I)
	END DO

!	SQNORMGRAD=ZERO
!	DO I=1,N
!		SQNORMGRAD=SQNORMGRAD+GRAD(I)*GRAD(I)
!	END DO

!	WRITE(6,*)'ENERGY=',R_ENER
!	WRITE(7,*)'ENERGY=',R_ENER,'   --->   ||g||^2 = ',SQNORMGRAD

 


END SUBROUTINE GRADEN
