!      PROJECT TITLE:  correlated gaussian geminals for molecular hydrogen
!MAIN EXECUTION FILE:  cgauss.F90
!
!   Program FILENAME:  cgauss.exe
!        Lead Author:  D. Gilmore
!    Other Author(s):  L. Adamowicz, D. Kinghorn
!
!        LAST EDITED:  23may96:wed/08:05
!

!NOTES FOR FUTURE EDITS:
!===================================================================================================
!	1.  unnecessary sharing of modules should be eliminated
!===================================================================================================
!
!

!-------------------------------------------------------&   !=======================================
SUBROUTINE	Print_Wave_Function_Stats
!-----------------------------------------------------------
	USE scalars_to_be_shared
	USE vectors_to_be_shared
	USE matrices_to_be_shared

	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	IMPLICIT          INTEGER    (I-N)
!-----------------------------------------------------------
!	Texp=ZERO
!	Vexp=ZERO
	HE  =ZERO
	SE  =ZERO
	Sexp=ZERO

	DO I=1,M
		DO J=1,M
			Sexp=Sexp + (CO(I)*S(I,J)*CO(J))
!			Texp=Texp + (CO(I)*Kinetic(I,J)*CO(J))
!			Vexp=Vexp + (CO(I)*Potential(I,J)*CO(J))
		END DO
	END DO
	
!	Texp=Texp/Sexp
!	Vexp=Vexp/Sexp

	CALL Calculate_Rayleigh_Quotient

    CALL BIOUT('==================================='	&
				//'=============')

	DO I=6,7
		WRITE(I,200)M,EMIN

		WRITE(I,300)R_ENER

!		WRITE(I,500)Texp

!		WRITE(I,600)Vexp

!		WRITE(I,400)-Vexp**TWO /(FOUR*Texp)

!		WRITE(I,700)(-TWO*Texp)/Vexp
	END DO
	
    CALL BIOUT('==================================='	&
				//'=============')

 200 FORMAT(' n = ',I4,'      E[(H-ES)c=0] =',F40.15)
 300 FORMAT('                c`Hc / c`Sc ='  ,F40.15)
!400 FORMAT('              -<V>^2 / 4<T> ='  ,F20.15)
!500 FORMAT('                        <T> ='  ,F20.15)
!600 FORMAT('                        <V> ='  ,F20.15)
!700 FORMAT('              -2 <T> /  <V> ='  ,F20.15)

END SUBROUTINE	Print_Wave_Function_Stats