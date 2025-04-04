SUBROUTINE GRADEN(N,X,G)

	USE scalars_to_be_shared
	USE vectors_to_be_shared
	USE matrices_to_be_shared
	USE strings_to_be_shared


    IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	IMPLICIT INTEGER(I-N)

	DIMENSION X(M*6),G(M*6)

	NB=(N/6)

	K=1
	DO I=1,NB
			T1(I)=X(K)	
			T2(I)=X(K+1)	
			T3(I)=X(K+2)	
			RZ1(I)=X(K+3)	
			RZ2(I)=X(K+4)	
		K=K+5
	END DO

	DO I=1,NB
			C(I)=X(NB*5+I)
	END DO


	CALL Cholesky_Condition_Parameters

	CALL Construct_Gradient


	DO I=1,NB*5+NB
		G(I)=GRAD(I)
	END DO

	SQNORMGRAD=ZERO
	DO I=1,NB*5+NB
		SQNORMGRAD=SQNORMGRAD+GRAD(I)*GRAD(I)
	END DO

!	WRITE(6,*)'ENERGY=',R_ENER
!	WRITE(6,*)'ENERGY=',R_ENER,'   --->   ||g||^2 = ',SQNORMGRAD
!	WRITE(7,*)'ENERGY=',R_ENER,'   --->   ||g||^2 = ',SQNORMGRAD

   OPEN (8, FILE = 'Restart_Output.DAT', ACCESS = 'SEQUENTIAL',STATUS = 'UNKNOWN')
   WRITE(8,*) '     LinCoeff:.   AA1:        AA2:        B12:        RZ1:       RZ2:'
   WRITE(8,*) 'FNC: ======================================================================'

	DO J=1,M
    WRITE(8,1111) J,C(J),AA1(J),AA2(J),B12(J),RZ1(J),RZ2(J)
 	END DO

   WRITE(8,*) '---------------------------------------------------------------------------'
   WRITE(8,*)'ENERGY=',R_ENER,'   --->   ||g||^2 = ',SQNORMGRAD
   WRITE(8,*)

   WRITE(8,*) ':BEGIN:'
   WRITE(8,*) M

	DO J=1,M
         WRITE(8,*) J,C(J),T1(J),T2(J),T3(J),RZ1(J),RZ2(J)
	END DO

   WRITE(8,*) ':END:'
	WRITE(8,*)
	CLOSE(8)

1111  FORMAT(I3,6D12.4)


END SUBROUTINE GRADEN

