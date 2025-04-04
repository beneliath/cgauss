SUBROUTINE Load_Restart_Point

	USE scalars_to_be_shared			  
	USE	vectors_to_be_shared			  
	USE	matrices_to_be_shared			  
	USE	strings_to_be_shared
	
	CHARACTER(LEN=256) :: rSTR			  

	CALL BIOUT(' STARTING FROM DESIGNATED RESTART POINT!')
	WRITE(6,*)' Point taken from File = ',RESTARTNAME
	WRITE(7,*)' Point taken from File = ',RESTARTNAME
	CALL BIOUT(' ')

	OPEN (5,FILE=RESTARTNAME)
	
		READ(5,*) rSTR
		write(6,*) rSTR
	
		DO WHILE (rSTR .NE. ':BEGIN:')
			READ(5,*) rSTR
		END DO

	READ(5,*) M
		WRITE(6,*)'NUMBER OF FUNCTIONS IN RESTART=',M
		WRITE(7,*)'NUMBER OF FUNCTIONS IN RESTART=',M
	CALL BIOUT(' ')

		DEALLOCATE(AA1,RZ1,AA2,RZ2,B12,C,S,H,GRAD,T1,	&
					T2,T3,CO,DS,DH )

		ALLOCATE(AA1(M),RZ1(M),AA2(M),RZ2(M),			&	!ALLOCATE NECESSARY MEMORY...
				B12(M),C(M),S(M,M),H(M,M),				&
				GRAD(M*6),T1(M),T2(M),T3(M),CO(M),		&
				DS(M*5,M),DH(M*5,M) )

	DO J=1,M
		READ(5,*) I,C(J),T1(J),T2(J),T3(J),RZ1(J),RZ2(J)
	END DO

   	CLOSE (5)

END SUBROUTINE Load_Restart_Point