SUBROUTINE FCN(NN,X,F)


	USE scalars_to_be_shared
	USE vectors_to_be_shared
	USE matrices_to_be_shared


    IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	IMPLICIT INTEGER (I-N)
	
	DIMENSION X(M*6)

	NB=(NN/6)

!	WRITE(6,*)'HIT FUNCTION'

	K=1
	DO I=1,NB
			T1(I) =X(K)	
			T2(I) =X(K+1)	
			T3(I) =X(K+2)	
			RZ1(I)=X(K+3)	
			RZ2(I)=X(K+4)	
		K=K+5
	END DO


	DO I=1,NB
			C(I)=X(NB*5+I)
	END DO


	CALL Cholesky_Condition_Parameters
	CALL Construct_S_and_H_Matrices

!	CALL Solve_Eigen_Value_Problem
	CALL Calculate_Rayleigh_Quotient

!	F=EMIN
	F=R_ENER

!WRITE(6,*)'FINISHED FUNCTION'


END SUBROUTINE FCN

