SUBROUTINE Gradient_Thread_1(ARG)

USE MT

INTEGER(4) :: ARG

AUTOMATIC

	CALL Gradient_Projection_0
	CALL Gradient_Projection_1

CALL EXITTHREAD(ARG)

END SUBROUTINE Gradient_Thread_1