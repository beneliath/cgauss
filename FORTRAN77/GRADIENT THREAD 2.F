SUBROUTINE Gradient_Thread_2(ARG)

USE MT

INTEGER(4) :: ARG

	CALL Gradient_Projection_2
	CALL Gradient_Projection_3

CALL EXITTHREAD(0)

END SUBROUTINE Gradient_Thread_2
