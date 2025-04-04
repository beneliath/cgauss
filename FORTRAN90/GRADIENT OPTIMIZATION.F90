SUBROUTINE Gradient_Optimization

	USE scalars_to_be_shared
	USE vectors_to_be_shared
	USE matrices_to_be_shared

    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	IMPLICIT INTEGER (I-N)


	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::	X,XSCALE,XGUESS,G

	DIMENSION IPARAM(7),RPARAM(7)

	EXTERNAL FCN
	EXTERNAL GRADEN
	EXTERNAL DU4INF

	ALLOCATE( G(M*6),X(M*6),XSCALE(M*6),XGUESS(M*6))

	NB=M
	N=NB*5+NB

	K=1
	DO I=1,NB
			X(K)  = T1(I)	
			X(K+1)= T2(I)	
			X(K+2)= T3(I)	
			X(K+3)= RZ1(I)	
			X(K+4)= RZ2(I)	 

			XGUESS(K)  = T1(I)	
			XGUESS(K+1)= T2(I)	
			XGUESS(K+2)= T3(I)	
			XGUESS(K+3)= RZ1(I)	
			XGUESS(K+4)= RZ2(I)	
		K=K+5
	END DO

	DO I=1,NB
			X(NB*5+I)=C(I)

			XGUESS(NB*5+I)=C(I)
	END DO


	DO I=1,N
		XSCALE(I)=ONE
	END DO

	FSCALE=ONE

	CALL DU4INF(IPARAM,RPARAM)
	IPARAM(3)=10000
	IPARAM(4)=40000
	IPARAM(5)=40000
!	IPARAM(6)=ZERO
!	RPARAM(1)=1.00D-6
!	RPARAM(2)=1.00D-10

!WRITE(6,*)'RPARAM2=',RPARAM(2)
!STOP

!	CALL BIOUT('IMSL:DUMING Optimization Begins...')
!	CALL BIOUT('----------------------------------')


	CALL DUMING(FCN,GRADEN,N,XGUESS,XSCALE,FSCALE,IPARAM, &
			RPARAM,X,FVALUE)

!	CALL BIOUT(' ')     
!     WRITE (6,99999)(IPARAM(L),L=3,5)
!      WRITE (7,99999)(IPARAM(L),L=3,5)
!	CALL BIOUT(' ')

99999 FORMAT ('  The number of iterations is ',						&
            10X, I4, /, '  The number of function evaluations is ',	&
            I4, /, '  The number of gradient evaluations is ', I4)

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



	DEALLOCATE (G,X,XSCALE,XGUESS)


END SUBROUTINE Gradient_Optimization