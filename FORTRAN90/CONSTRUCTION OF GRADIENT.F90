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
SUBROUTINE Construct_Gradient
!-----------------------------------------------------------
	USE  scalars_to_be_shared
	USE  vectors_to_be_shared
	USE	matrices_to_be_shared

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	IMPLICIT          INTEGER     (I-N)
!-----------------------------------------------------------
		CALL Calculate_Total_Gradient

!		CALL BIOUT('DH-MATRIX:')
!			CALL TAB(DH,M*5,M,M*5,M)                  
!		CALL BIOUT(' ')

!		CALL BIOUT('DS-MATRIX:')
!			CALL TAB(DS,M*5,M,M*5,M)                  
!		CALL BIOUT(' ')

    SUMH=ZERO
	SUMS=ZERO
	DO J=1,M
		DO I=1,M
			SUMH=SUMH+C(I)*H(I,J)*C(J)
			SUMS=SUMS+C(I)*S(I,J)*C(J)
		END DO
	END DO

	ENG=SUMH/SUMS
	SCALE=ONE/SUMS

DEALLOCATE(GRAD)
ALLOCATE(GRADM(M,M*5),GRADA(M*5),GRADC(M),GRAD(M*6))

	DO I=1,M
		DO J=1,M*5
			GRADM(I,J)=DH(J,I)-ENG*DS(J,I)
			GRADM(I,J)=SCALE*GRADM(I,J)
		END DO
	END DO


	NN=5
	DO J=1,M
		DO K=1,NN
			IC=(J-1)*NN+K
			GRADA(IC)=ZERO
			DO I=1,M
				IF(I.NE.J) THEN
					GRADA(IC)=GRADA(IC)+TWO*C(I)*C(J)*GRADM(I,IC)
				ELSE
					GRADA(IC)=GRADA(IC)+C(I)*C(J)*GRADM(I,IC)
				END IF
			END DO
		END DO
	END DO

    DO I=1,M
		GRADC(I) = ZERO
	END DO
	
    DO J=1,M
		DO I=1,M
		  GRADC(I)=GRADC(I)+TWO*SCALE*( H(I,J)-ENG*S(I,J) )*C(J)
		END DO
	END DO

	K=1
	DO I=1,M*NN,NN
		DA1=GRADA(I)
		DA2=GRADA(I+1)
		DB=GRADA(I+2)
		GRADA(I)=DA1*(TWO*T1(K)+T3(K))+DA2*T3(K)-DB*T3(K)
		GRADA(I+1)=TWO*T2(K)*DA2
		GRADA(I+2)=DA1*T1(K)+(TWO*T3(K)+T1(K))*DA2-DB*T1(K)
	K=K+1
    END DO

      NX=M*NN
	DO I=1,NX
		GRAD(I) = GRADA(I)
	END DO
	
    DO I=1,M
		GRAD(I+NX) = GRADC(I)
	END DO

DEALLOCATE(GRADM,GRADA,GRADC)



END SUBROUTINE Construct_Gradient