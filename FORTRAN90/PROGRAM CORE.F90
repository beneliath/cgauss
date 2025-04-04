!PROJECT TITLE:  correlated gaussian geminals for molecular hydrogen
!
!PROJECT SUMMARY:
!	...incorporates 4-fold symmetry projection and optimization of wavefunction based on an
!	imsl gradient (TN) algorythm utilizing the derivatives of the molecular integrals.  Basis set
!	is stochastically increased followed by gradient optimization (iteratively).
!
!  Program FILENAME:  cgauss.exe
!	    Lead Author:  D. Gilmore
!   Other Author(s):  L. Adamowicz, D. Kinghorn
!
!LANGUAGE SPECIFICATION:  Fortran90
!  DEVELOPMENT PLATFORM:  Microsoft Fortran Powerstation (professional) w/ IMSL
!
!   CURRENT EXECUTABLE SIZE:  ?MB (22may96)
!CURRENT MEMORY REQUIREMENT:  ?MB (22may96)
!
!LAST EDITED:  22may96:wed/15:30
!

!NOTES FOR FUTURE EDITS:
!===================================================================================================
!	1.	suspect that there could be a problem in the gradient...
!		entire gradient should be finite-difference-tested for many points
!
!	2.	code needs to be further optimized for lower memory usage:  TOO MANY TEMP ARRAYS/MATRICES!
!		(i.e.  OVV1-OVV4,Kinetic,Potential,...maybe more)
!
!	3.	haven't checked the modified derivative of the error function yet:  F0(z) -> dF0(z)
!
!	4.  find out the memory requirement for the IMSL duming: routine
!
!	5.  unnecessary sharing of modules should be eliminated
!
!	6.	nuclear geometry variables(EX1-EZ1,EX2-EZ2) can be reduced under cylindrical symmetry
!
!	7.	under cylindrical symmetry, vectors RX1,RY1,RX2,RY2 can be eliminated
!
!	8.	change ALL IMPLICIT typing to explicit using IMPLICIT NONE
!===================================================================================================
!
!

!-------------------------------------------------------&   !=======================================
PROGRAM	cgauss
!-----------------------------------------------------------
	USE          scalars_to_be_shared						!SHARED MODULES...
	USE	         vectors_to_be_shared			  
	USE	        matrices_to_be_shared			  
	USE	         strings_to_be_shared			  
	USE	program_settings_to_be_shared

		IMPLICIT DOUBLE PRECISION(A-H,O-Z)					!implicit real typing
		IMPLICIT          INTEGER    (I-N)					!implicit int  typing

			CHARACTER (LEN=10) :: TT						!timestamp var character typing...
			CHARACTER (LEN=5)  :: TZ
			CHARACTER (LEN=8)  :: TD

			CHARACTER (LEN=80) :: CHARSTRING				!general character storage for BIOUT

				INTRINSIC RANDOM_NUMBER						!definition of intrinsic functions...
				INTRINSIC RANDOM_SEED
!-----------------------------------------------------------
	CALL Default_Program_Settings							!INITIALIZE PROGRAM SETTINGS

	IF (IMODE.EQ. 1) CALL ProcessPipedInput					!PROCESS PIPED INPUT

	OPEN (7, FILE = FILENAME, ACCESS = 'SEQUENTIAL')		!OPEN FILE:  file.out

	CALL Print_Preliminary_Matter							!PRINT PROGRAM HEADER DATA...
															!---------------------------------------
		CALL BIOUT(':PROGRAM BEGINS:')						!	BIOUT -> screen & file.out
		CALL Get_Current_Time(TD,TT,TZ)
		CHARSTRING='Date: '//TD//' | Time: '//TT//' / '	&
										    //TZ//' gmt'
		CALL BIOUT(CHARSTRING)
															!---------------------------------------
	IF (TEST  .EQ. 1) M=2									!TEST CASE IS A SIMPLE 2x2
	
		ALLOCATE(AA1(M),RZ1(M),AA2(M),RZ2(M),			&	!ALLOCATE NECESSARY MEMORY...
				B12(M),C(M),S(M,M),H(M,M),				&
				GRAD(M*6),T1(M),T2(M),T3(M),CO(M),		&
				DS(M*5,M),DH(M*5,M) )
															
	IF (RESTART .EQ. 0 .AND. TEST .NE. 1)				&	!NORESTART: GUESS AT INITIAL POINT
		CALL Guess_Random_Starting_Data						!

	IF (RESTART .EQ. 1 .AND. TEST .NE. 1)				&	!  RESTART: LOAD RESTART POINT
		CALL Load_Restart_Point								!

	IF (TEST    .EQ. 1                  )				&	!     TEST: LOAD TESTING PARAMETERS
		CALL Load_Default_Testing_Parameters				!


!MICROSOFT PROPRIETARY OBJECT!
!-----------------------------------------------------------
!...this block must be removed for generalized porting to
!	other platforms...
!-----------------------------------------------------------
	ELAPSED_TIME = TIMEF()
!-----------------------------------------------------------
!                              MICROSOFT PROPRIETARY OBJECT!





	CALL BIOUT('Defining nuclear geometry...')
	CALL Define_Geometry									!DEFINE NUCLEAR GEOMETRY

	CALL Cholesky_Condition_Parameters						!CHOLESKY PRECONDITIONING
	CALL BIOUT('Cholesky preconditioning complete...')
!	stop

	CALL BIOUT('constructing S & H matrices...')
	CALL Construct_S_and_H_Matrices							
	CALL BIOUT('S & H matrix constructions complete...')
	
	CALL BIOUT('solving eigen-problem -> energy...')
	CALL Solve_Eigen_Value_Problem							!CALCULATE PRE-OPTIMIZED ENERGY

	CALL BIOUT('pre-optimization stats:')
	CALL Print_Wave_Function_Stats


!MICROSOFT PROPRIETARY OBJECT!
!-----------------------------------------------------------
!...this block must be removed for generalized porting to
!	other platforms...
!-----------------------------------------------------------
	ELAPSED_TIME = TIMEF()

	DO I=6,7
	WRITE(I,900) ELAPSED_TIME,         ' sec =',			&
				 ELAPSED_TIME/60,	   ' min =',			&
				 ELAPSED_TIME/60/60,   ' hrs =',			&
				 ELAPSED_TIME/60/60/24,' dys'
	END DO
!-----------------------------------------------------------
!                              MICROSOFT PROPRIETARY OBJECT!

!stop

	IF (PREOPT .EQ. 1) THEN										!PRE-OPT BEFORE STOCHASTIC GROWTH?

		DO I=1,M												!load eigenprob lin coeff solution
			C(I)=CO(I)											!into standard lin coeff vector.
		END DO

!DELTA=1.00D-09

!	CALL Construct_S_and_H_Matrices							
!	CALL Calculate_Rayleigh_Quotient
!	R1=R_ENER

	CALL Construct_Gradient

!	C(1)=C(1)+DELTA
!	CALL Cholesky_Condition_Parameters						!CHOLESKY PRECONDITIONING
!	CALL Construct_S_and_H_Matrices							
!	CALL Calculate_Rayleigh_Quotient
!	R2=R_ENER

!WRITE(6,*)'FINITE DIFFERENCE=',(R2-R1)/DELTA
!WRITE(6,*)'ANALYTIC GRADIENT=',GRAD(11)





!stop
		CALL Gradient_Optimization

		CALL Cholesky_Condition_Parameters

		CALL Construct_S_and_H_Matrices
		CALL Solve_Eigen_Value_Problem

	END IF

	ENER=EMIN

!OPEN (8, FILE='convergence.dat', ACCESS = 'SEQUENTIAL',STATUS = 'UNKNOWN')
!write(8,*) M,EMIN
!CLOSE(8)



	CALL Print_Wave_Function_Stats
	CALL Simple_Wave_Function_Printout							!printout wavefunction
	


	RUP=20.
	RDOWN=.1

	RSTEP=(RUP-RDOWN)/100.
	IRIGHT=INT((RUP-BD)/RSTEP)
	ILEFT=INT((BD-RDOWN)/RSTEP)

	OPEN(5,FILE="RIGHT.DAT")
	write(5,*)'   '
	CLOSE(5)

	DO I=1,IRIGHT

	CALL Define_Geometry									!DEFINE NUCLEAR GEOMETRY

	CALL Construct_S_and_H_Matrices							

	CALL Calculate_Rayleigh_Quotient

	CALL Construct_Gradient
		CALL Gradient_Optimization
		CALL Cholesky_Condition_Parameters
		CALL Construct_S_and_H_Matrices
		CALL Calculate_Rayleigh_Quotient

	OPEN(5,FILE="RIGHT.DAT",ACCESS='APPEND')
	WRITE(6,*)R_ENER, BD
	WRITE(5,*)R_ENER, BD
	CLOSE(5)

	BD=BD+RSTEP

	END DO


stop
!	CALL BIOUT('Iterative Stochaistic Optimization :BEGINS:')
	CALL RANDOM_SEED()

    FMAX = 10.00D+00
	FMIN = -7.00D+00
	ZMAX =  1.70D+00
	BL0  =  0.70D+00
	BL1  =  1.30D+00
	BLURBEGIN = FIVE
	BLUREND   =-SEVEN
	BETAMAX   = FIVE
	BETAMIN   =-FIVE

	IVIEW = 4000
	ZNUC  = BD/TWO

	A=FMAX-FMIN
	B=FMIN

	CALL RANDOM_NUMBER(TEMP1)

	K =-DLOG(0.0001)/(ZMAX-ZNUC)
	KA=-DLOG(0.0001)/(ZNUC)

	DO JJJ=M,ISTEND

	DO I=1,IVIEW
	TEMP1=.5
	IF (TEMP1.GE..5) THEN

			CALL RANDOM_NUMBER(TEMP1)
 4141		CALL RANDOM_NUMBER(TEMP2)
	IF (TEMP2.LT..5) GOTO 4141
	TOF1=(FMAX-FMIN)*TEMP1+FMIN
	XAA1=TEMP2*10**TOF1

			CALL RANDOM_NUMBER(TEMP1)
 5151		CALL RANDOM_NUMBER(TEMP2)
	IF (TEMP2.LT..5) GOTO 5151
	TOF2=(FMAX-FMIN)*TEMP1+FMIN
	XAA2=TEMP2*10**TOF2
	XRZ1=-DLOG((TOF1+DABS(FMIN))/(FMAX-FMIN))/K + ZNUC
	XRZ2=-DLOG((TOF2+DABS(FMIN))/(FMAX-FMIN))/K + ZNUC
	CALL RANDOM_NUMBER(TEMP1)
	IF (TEMP1.GE..5) XRZ1=-1.*XRZ1
	CALL RANDOM_NUMBER(TEMP1)
	IF (TEMP1.GE..5) XRZ2=-1*XRZ2
	IF (BL0.GT.0.) THEN
		CALL RANDOM_NUMBER(TEMP1)
		XC=1
		IF (TEMP1.GE..5) XC=-1.
		CALL RANDOM_NUMBER(TEMP1)
		IF (TOF1 .LE. BLURBEGIN) THEN
		XRZ1=XRZ1+XC*DSQRT(-2.*BL0**2 * DLOG(TEMP1))* (BLURBEGIN-TOF1)/(BLURBEGIN-BLUREND)
		END IF	

		CALL RANDOM_NUMBER(TEMP1)
		XC=1
		IF (TEMP1.GE..5) XC=-1.
		CALL RANDOM_NUMBER(TEMP1)
		IF (TOF2 .LE. BLURBEGIN) THEN
		XRZ2=XRZ2+XC*DSQRT(-2.*BL0**2 * DLOG(TEMP1))* (BLURBEGIN-TOF2)/(BLURBEGIN-BLUREND)
		END IF

	END IF


	ELSE
			CALL RANDOM_NUMBER(TEMP1)
 2121		CALL RANDOM_NUMBER(TEMP2)
	IF (TEMP2.LT..5) GOTO 2121
	TOF1=(FMAX-FMIN)*TEMP1+FMIN
	XAA1=TEMP2*10**TOF1

			CALL RANDOM_NUMBER(TEMP1)
 3131		CALL RANDOM_NUMBER(TEMP2)
	IF (TEMP2.LT..5) GOTO 3131
	TOF2=(FMAX-FMIN)*TEMP1+FMIN
	XAA2=TEMP2*10**TOF2

	XRZ1=-DLOG((TOF1+DABS(FMIN))/(FMAX-FMIN))/KA - ZNUC
	XRZ2=-DLOG((TOF2+DABS(FMIN))/(FMAX-FMIN))/KA - ZNUC
	CALL RANDOM_NUMBER(TEMP1)
	IF (TEMP1.GE..5) XRZ1=-1.*XRZ1
	CALL RANDOM_NUMBER(TEMP1)
	IF (TEMP1.GE..5) XRZ2=-1*XRZ2

	IF (BL1.GT.0.) THEN
		CALL RANDOM_NUMBER(TEMP1)
		XC=1
		IF (TEMP1.GE..5) XC=-1.
		CALL RANDOM_NUMBER(TEMP1)
		IF (TOF1 .LE. BLURBEGIN) THEN
		XRZ1=XRZ1+XC*DSQRT(-2.*BL1**2 * DLOG(TEMP1))* (BLURBEGIN-TOF1)/(BLURBEGIN-BLUREND)
		END IF
	

		CALL RANDOM_NUMBER(TEMP1)
		XC=1
		IF (TEMP1.GE..5) XC=-1.
		CALL RANDOM_NUMBER(TEMP1)
		IF (TOF2 .LE. BLURBEGIN) THEN
		XRZ2=XRZ2+XC*DSQRT(-2.*BL1**2 * DLOG(TEMP1))* (BLURBEGIN-TOF2)/(BLURBEGIN-BLUREND)
		END IF

	END IF


	END IF
	CALL RANDOM_NUMBER(TEMP1)
	CALL RANDOM_NUMBER(TEMP2)
	XB12=DABS(TEMP1*10.**(((BETAMAX+1.)-(BETAMIN+1))*TEMP2+(BETAMIN+1)))

!	CALL HAMILTONp(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,
!     1                    EVEC,S,H,EMIN,XAA1,XRX1,XRY1,XRZ1,
!     2					XAA2,XRX2,XRY2,XRZ2,XB12)


CALL Energy_Diagnosis(XAA1,XRZ1,XAA2,XRZ2,XB12,EMIN)


	IF (I .EQ. 1) ENER=EMIN

	IF (EMIN .LT. ENER) THEN 
		ENER=EMIN
		TA1=XAA1
		TA2=XAA2
		TB12=XB12
		TRZ1=XRZ1
		TRZ2=XRZ2
	END IF
!		TA1= XAA1*XAA1             + XAA1*XB12
!        TA2= XAA2*XAA2 + XB12*XB12 + XAA1*XB12


      END DO


! 	CALL BIOUT(' ')

!	WRITE(6,*)'BEST DIAGNOSTIC VALUE IN ',IVIEW,' GUESSES=',ENER
!pause
!	WRITE(6,*)'***************************'
!	WRITE(6,*)'ADDING THIS FUNCTION TO Psi'
!	WRITE(6,*)'***************************'
	ALLOCATE(CSTORET(M),TSTORE1T(M),TSTORE2T(M),TSTORE3T(M),RZ1STORET(M),RZ2STORET(M))
	DO III=1,M
		CSTORET(III)=C(III)
		TSTORE1T(III)=AA1(III)
		TSTORE2T(III)=AA2(III)
		TSTORE3T(III)=B12(III)
		RZ1STORET(III)=RZ1(III)
		RZ2STORET(III)=RZ2(III)
	END DO
	M=M+1
		DEALLOCATE(AA1,RZ1,AA2,RZ2,B12,C,S,H,GRAD,T1,	&
					T2,T3,CO,DS,DH )

		ALLOCATE(AA1(M),RZ1(M),AA2(M),RZ2(M),			&	!ALLOCATE NECESSARY MEMORY...
				B12(M),C(M),S(M,M),H(M,M),				&
				GRAD(M*6),T1(M),T2(M),T3(M),CO(M),		&
				DS(M*5,M),DH(M*5,M) )

	DO III=1,M-1
	C(III)=CSTORET(III)
	AA1(III)=TSTORE1T(III)
	AA2(III)=TSTORE2T(III)
	B12(III)=TSTORE3T(III)
	RZ1(III)=RZ1STORET(III)
	RZ2(III)=RZ2STORET(III)
	END DO
	AA1(M)=TA1
	AA2(M)=TA2
	B12(M)=TB12
	RZ1(M)=TRZ1
	RZ2(M)=TRZ2

	DEALLOCATE(CSTORET,TSTORE1T,TSTORE2T,TSTORE3T,RZ1STORET,RZ2STORET)


	CALL Construct_S_and_H_Matrices

	CALL Solve_Eigen_Value_Problem

	DO NNN=1,M
		C(NNN)=CO(NNN)
	END DO

!	WRITE(6,*)'NOW WE HAVE',M,' FUNCTIONS...'
!	WRITE(6,*)'EXACT ENERGY=',EMIN

!	WRITE(6,*)'********************************'
!	WRITE(6,*)'PERFORMING GRADIENT OPTIMIZATION'
!	WRITE(6,*)'********************************'

	DO NNN=1,M
		T1(NNN)=DSQRT(AA1(NNN)+B12(NNN))
		T2(NNN)=DSQRT(AA2(NNN)+B12(NNN)- ((B12(NNN)**TWO)/(AA1(NNN)+B12(NNN))))
		T3(NNN)=-B12(NNN)/DSQRT(AA1(NNN)+B12(NNN))
	END DO

	CALL Construct_Gradient
!	WRITE(6,*)'PRE-OPTIMIZATION ENERGY --->',EMIN

 	CALL Gradient_Optimization

	CALL Cholesky_Condition_Parameters

	CALL Construct_S_and_H_Matrices
	CALL Solve_Eigen_Value_Problem

	CALL Print_Wave_Function_Stats
!MICROSOFT PROPRIETARY OBJECT!
!-----------------------------------------------------------
!...this block must be removed for generalized porting to
!	other platforms...
!-----------------------------------------------------------
	ELAPSED_TIME = TIMEF()

	DO I=6,7
	WRITE(I,900) ELAPSED_TIME,         ' sec =',			&
				 ELAPSED_TIME/60,	   ' min =',			&
				 ELAPSED_TIME/60/60,   ' hrs =',			&
				 ELAPSED_TIME/60/60/24,' dys'
	END DO
!-----------------------------------------------------------
!                              MICROSOFT PROPRIETARY OBJECT!


	OPEN (8, FILE = 'convergence.dat', ACCESS = 'APPEND',STATUS = 'UNKNOWN')
	write(8,*) M,EMIN
	CLOSE(8)

  OPEN (8, FILE = 'Restart_Output.DAT', ACCESS = 'SEQUENTIAL',STATUS = 'UNKNOWN')
   WRITE(8,*) '     LinCoeff:.   AA1:        AA2:        B12:        RZ1:       RZ2:'
   WRITE(8,*) 'FNC: ======================================================================'

	DO J=1,M
    WRITE(8,1111) J,C(J),AA1(J),AA2(J),B12(J),RZ1(J),RZ2(J)
 	END DO

SQNORMGRAD=ZERO
DO I=1,M*6
	SQNORMGRAD=SQNORMGRAD+(GRAD(I)*GRAD(I))
END DO

   WRITE(8,*) '---------------------------------------------------------------------------'
   WRITE(8,*)'ENERGY=',EMIN,'   --->   ||g||^2 = ',SQNORMGRAD
   WRITE(8,*)

   WRITE(8,*) ':BEGIN:'
   WRITE(8,*) M

	DO J=1,M
         WRITE(8,*) J,C(J),T1(J),T2(J),T3(J),RZ1(J),RZ2(J)
	END DO

   WRITE(8,*) ':END:'
	WRITE(8,*)
	CLOSE(8)

	
	DO J=1,M
		C(J)=CO(J)
	END DO

	END DO



!CLEAN-UP...
!=====================================================

	CALL BIOUT('PROGRAM TERMINATES!')
	CALL Get_Current_Time(TD,TT,TZ)
	WRITE(6,*)'TIME STAMP==> Date: ',TD,' | Time: ',TT,' / ',TZ,'GMT'		
	WRITE(7,*)'TIME STAMP==> Date: ',TD,' | Time: ',TT,' / ',TZ,'GMT'		

	CALL BIOUT(' ')

!MICROSOFT PROPRIETARY STUFF!
	ELAPSED_TIME = TIMEF()

	WRITE(6,9900)ELAPSED_TIME,' sec =',ELAPSED_TIME/60,' min =',		&
		ELAPSED_TIME/60/60,' hrs =',ELAPSED_TIME/60/60/24,' dys'
	WRITE(7,9900)ELAPSED_TIME,' sec =',ELAPSED_TIME/60,' min =',		&
		ELAPSED_TIME/60/60,' hrs =',ELAPSED_TIME/60/60/24,' dys'

9900 FORMAT (' Elapsed Time: ',F10.3,A6,F10.3,A6,F10.3,A6,F10.3,A4)

	CALL BIOUT(' ')
!MICROSOFT PROPRIETARY STUFF!

	CLOSE (5)
	CLOSE (7)
	CLOSE (8)
	CLOSE (9)

		DEALLOCATE(AA1,RZ1,AA2,RZ2,B12,C,S,H,GRAD,T1,	&
					T2,T3,CO,DS,DH )


!FORMAT STATEMENTS:
!-----------------------------------------------------------
292	FORMAT(I4,4D24.8)
900 FORMAT (' Elapsed Time: ',F10.3,A6,F10.3,A6,F10.3,A6,F10.3,A4)
1111  FORMAT(I3,6D12.4)
!-----------------------------------------------------------
END PROGRAM	cgauss