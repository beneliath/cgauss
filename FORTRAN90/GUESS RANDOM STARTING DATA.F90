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
SUBROUTINE Guess_Random_Starting_Data
!-----------------------------------------------------------
	USE          scalars_to_be_shared
	USE          vectors_to_be_shared
	USE program_settings_to_be_shared

	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	IMPLICIT          INTEGER    (I-N)
	
	INTRINSIC RANDOM_NUMBER
	INTRINSIC RANDOM_SEED
!-----------------------------------------------------------!	---:END OF FORMATTED COMMENTS:---

	CALL RANDOM_SEED()

   		CALL BIOUT('GENERATING RANDOM GAUSSIAN GEMINAL DATA: (printout follows...)')
		CALL BIOUT(' ')


	TEMPBD=BD/TWO

      DO I=1,M

	    CALL RANDOM_NUMBER(TEMP)
			AA1(I)=DBLE(TEMP)*1.0D-01		!alphas are randomized at order of magnitude   (-1)
        CALL RANDOM_NUMBER(TEMP)
			AA2(I)=DBLE(TEMP)*1.0D-01		!alphas are randomized at order of magnitude   (-1)
        CALL RANDOM_NUMBER(TEMP)
			B12(I)=DBLE(TEMP)*1.0D-03		!beta is randomized at order of magnitude      (-8)
		CALL RANDOM_NUMBER(TEMP)
			C  (I)=DBLE(TEMP)*1.0D-02		!lin coeff is randomized at order of magnitude (-2)

	   	RZ1(I)= TEMPBD
		RZ2(I)=-TEMPBD
			
      END DO

	CALL Cholesky_Deconditioning			!extract cholesky factors from alphas & betas
	CALL Simple_Wave_Function_Printout		!printout wavefunction w/ low significance

END SUBROUTINE Guess_Random_Starting_Data