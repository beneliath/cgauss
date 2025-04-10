      SUBROUTINE GRADOPT(M,T1,RX1,RY1,RZ1,T2,RX2,RY2,RZ2,T3,EVEC,GRAD)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INTEGER M,ERROR
C
      REAL*8, DIMENSION(M) ::	RX1,RY1,RZ1,RX2,RY2,RZ2,
     -						B12,EVEC,T1,T2,T3
C
	REAL*8, DIMENSION(M*6) :: GRAD
C
	REAL*8	X[ALLOCATABLE](:),XSCALE[ALLOCATABLE](:),
	-		XGUESS[ALLOCATABLE](:)
C
	DIMENSION IPARAM(7),RPARAM(7)
C
	EXTERNAL FCN
	EXTERNAL GRADEN
	EXTERNAL DU4INF
C
	ALLOCATE( X(M*6),XSCALE(M*6),XGUESS(M*6),STAT=ERROR)
C
	IF (ERROR .NE. 0) STOP '!NOT ENOUGH MEMORY SPACE AVAILABLE ON THIS
	- BOX! ...ABORTING!'
C
	NB=M
	N=NB*5+NB
C
	K=1
	DO I=1,NB
			X(K)  =T1(I)	
			X(K+1)=T2(I)	
			X(K+2)=T3(I)	
			X(K+3)=RZ1(I)	
			X(K+4)=RZ2(I)	
C
			XGUESS(K)  =T1(I)	
			XGUESS(K+1)=T2(I)	
			XGUESS(K+2)=T3(I)	
			XGUESS(K+3)=RZ1(I)	
			XGUESS(K+4)=RZ2(I)	
		K=K+5
	END DO
C
	DO I=1,NB
			X(NB*5+I)=EVEC(I)
C
			XGUESS(NB*5+I)=EVEC(I)
	END DO
C
C
	DO I=1,N
		XSCALE(I)=1.0D+00
	END DO
C
	FSCALE=1.0D+00
	CALL DU4INF(IPARAM,RPARAM)
	IPARAM(3)=5000
	IPARAM(4)=5000
	IPARAM(5)=5000
C	RPARAM(1)=1.00D-07
C
C	DO I=1,N
C		INFO(I)=1
C	END DO
C	
C	CALL DCHGRD(FCN,GRADEN,N,X,INFO)
C	STOP
	CALL BIOUT('IMSL:DUMING Optimization Begins...')
	CALL BIOUT('----------------------------------')
	CALL DUMING(FCN,GRADEN,N,XGUESS,XSCALE,FSCALE,IPARAM,
	-		RPARAM,X,FVALUE)
C
	CALL BIOUT(' ')     
      WRITE (6,99999)(IPARAM(L),L=3,5)
      WRITE (7,99999)(IPARAM(L),L=3,5)
C
99999 FORMAT ('  The number of iterations is ',
     &       10X, I4, /, '  The number of function evaluations is ',
     &       I4, /, '  The number of gradient evaluations is ', I4)
C
	K=1
	DO I=1,NB
			T1(I)=X(K)
			T2(I)=X(K+1)
			T3(I)=X(K+2)	
			RZ1(I)=X(K+3)
			RZ2(I)=X(K+4)
C
		K=K+5
	END DO
	DO I=1,NB
			EVEC(I)=X(NB*5+I)
	END DO
C
C
C
C
	DEALLOCATE (X,XSCALE,XGUESS)
      RETURN
      END

      SUBROUTINE NGRADOPT(M,T1,RX1,RY1,RZ1,T2,RX2,RY2,RZ2,T3,EVEC,GRAD)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INTEGER M,ERROR
C
      REAL*8, DIMENSION(M) ::	RX1,RY1,RZ1,RX2,RY2,RZ2,
     -						B12,EVEC,T1,T2,T3
C
	REAL*8, DIMENSION(M*6) :: GRAD
C
	REAL*8	X[ALLOCATABLE](:),XSCALE[ALLOCATABLE](:),
	-		XGUESS[ALLOCATABLE](:)
C
	DIMENSION IPARAM(7),RPARAM(7)
C
	EXTERNAL FCN
	EXTERNAL GRADEN
	EXTERNAL DU4INF
C
	ALLOCATE( X(M*6),XSCALE(M*6),XGUESS(M*6),STAT=ERROR)
C
	IF (ERROR .NE. 0) STOP '!NOT ENOUGH MEMORY SPACE AVAILABLE ON THIS
	- BOX! ...ABORTING!'
C
	NB=M
	N=NB*5+NB
C
	K=1
	DO I=1,NB
			X(K)  =T1(I)	
			X(K+1)=T2(I)	
			X(K+2)=T3(I)	
			X(K+3)=RZ1(I)	
			X(K+4)=RZ2(I)	
C
			XGUESS(K)  =T1(I)	
			XGUESS(K+1)=T2(I)	
			XGUESS(K+2)=T3(I)	
			XGUESS(K+3)=RZ1(I)	
			XGUESS(K+4)=RZ2(I)	
		K=K+5
	END DO
C
	DO I=1,NB
			X(NB*5+I)=EVEC(I)
C
			XGUESS(NB*5+I)=EVEC(I)
	END DO
C
C
	DO I=1,N
		XSCALE(I)=1.0D+00
	END DO
C
	FSCALE=1.0D+00
	CALL DU4INF(IPARAM,RPARAM)
	IPARAM(1)=1
	IPARAM(3)=5000
	IPARAM(4)=5000
	IPARAM(5)=5000
	CALL DU4INF(IPARAM,RPARAM)
C
C	DO I=1,N
C		INFO(I)=1
C	END DO
C	
C	CALL DCHGRD(FCN,GRADEN,N,X,INFO)
C	STOP
	CALL BIOUT('IMSL:DUMINF Optimization Begins...')
	CALL BIOUT('----------------------------------')
	CALL DUMINF(FCN,N,XGUESS,XSCALE,FSCALE,IPARAM,RPARAM,X,FVALUE)
C     
	K=1
	DO I=1,NB
			T1(I)=X(K)
			T2(I)=X(K+1)
			T3(I)=X(K+2)	
			RZ1(I)=X(K+3)
			RZ2(I)=X(K+4)
C
		K=K+5
	END DO
	DO I=1,NB
			EVEC(I)=X(NB*5+I)
	END DO
C
C
C
	DEALLOCATE( X,XSCALE,XGUESS)
C
      RETURN
      END
