DOUBLE PRECISION FUNCTION F0(Z)

	USE scalars_to_be_shared

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      IF(Z .LE. 16.3578D+00) THEN

      F0=DSQRT((((((C1 *Z		&
                  + C2)*Z		&
                  + C3)*Z		&
                  + C4)*Z		&
                  + C5)*Z		&
                  + ONE)/		&
             ((((((C6 *Z		&
                  + C7)*Z		&
                  + C8)*Z		&
                  + C9)*Z		&
                  + C10)*Z		&
                  + C11)*Z		&
                  + ONE))

		ELSE
		    R=TWO*Z
			R3=PI/TWO
			F0=DSQRT(R3/R)

	   END IF



END FUNCTION F0
