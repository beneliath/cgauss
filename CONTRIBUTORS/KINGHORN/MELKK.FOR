      SUBROUTINE  MELKK( N, NMAX, VLK, VLL, MASS, CHRG,P,
     &                   SKL, TKL, VKL, DSK, DSL, DTK, DTL, DVK, DVL,
     &                   LK, LL, AK, AL, INVAKL, MAK, MAL, IAAK, IAAL,
     &                   IALK, IALL, TMP1, TMP2, TMPV1, TMPV2,
     &                   PLL, PIA )
c
c --- Donald B. Kinghorn
c     Dept of Chemistry
c     Washington State University
c     April 3 1995
c
c     Last modified May 12 1995 DBK
c
	 IMPLICIT NONE
c     ..Scalar Arguments..
      INTEGER             N, NMAX
      DOUBLE PRECISION    SKL, TKL, VKL
c     ..Array Arguments..
      DOUBLE PRECISION    VLK(*), VLL(*), CHRG(*),
     &                    MASS(NMAX,*), P(NMAX,*),
     &                    DSK(*), DSL(*), DTK(*), DTL(*),
     &                    DVK(*), DVL(*),    
     &                    LK(NMAX,*), LL(NMAX,*),
     &                    AK(NMAX,*), AL(NMAX,*),
     &                    INVAKL(NMAX,*), MAK(NMAX,*),
     &                    MAL(NMAX,*),
     &                    IAAK(NMAX,*), IAAL(NMAX,*),
     &                    IALK(NMAX,*), IALL(NMAX,*),
     &                    TMP1(NMAX,*), TMP2(NMAX,*),
     &                    TMPV1(*), TMPV2(*),
     &                    PLL(NMAX,*), PIA(NMAX,*)      
c =====================================================================
c Purpose
c =======
c
c ****DIAGONAL ELEMENTS****
c ****THIS CODE IS MELKL WITH THE VECTORS FOR THE K AND L GRADIENT
c ****COMPONENTS ADDED TOGETER AT THE END. THE INPUT IS OF COURSE
c ****VLL=VLK. THE VECTORS DSL,DTL AND DVL ARE JUST TEMP VECTORS
c ****FOR STORING THE SYMMETRY TRANSFORMED PART OF THE FORMULAS 
c Correlated Gaussian integral and derivative formulas 
c This subroutine computes the kl(th) matrix element component
c for the overlap kinetic energy and potential energy and
c the N(N+1)/2 x 1 vectors of derivatives with respect to
c the elements in vechLk (VLK) and vechLl (VLL) 
c 
c Arguments
c =========
c
c N       (input) INTEGER
c         Size of matrices AK, AL, MASS, etc.  ( # of particles - 1 )
c
c NMAX    (input) INTEGER
c         Size of arrays AK, AL, MASS, etc.  ( NMAX >= N )     
c
c VLK     (input) DOUBLE PRECISION vector, dimension (N*(N+1)/2)
c         vech(LK) Independent variables for energy functional
c
c VLL     (input) DOUBLE PRECISION vector, dimension (N*(N+1)/2)
c         vech(LK)
c
c MASS    (input) DOUBLE PRECISION matrix, dimension (N,N)
c         Matrix of mass constants for the kinetic energy integrals
c         see MKMASS.FOR for details
c
c CHRG    (input) DOUBLE PRECISION vector, dimension (N*(N+1)/2)
c         Vector of charge products for the potential energy integrals
c         see MKCHRG.FOR for details
c
c P       (input) DOUBLE PRECISION matrix dimension (N,N)
c         Transformation matrix for the symmetry projection     
c
c SKL     (output) DOUBLE PRECISION
c         Overlap integral 
c
c TKL     (output) DOUBLE PRECISION
c         Kinetic energy integral
c
c VKL     (output) DOUBLE PRECISION
c         Potential energy integral
c
c DSK     (output) DOUBLE PRECISION vector, dimension (N*(N+1)/2)
c         Derivative of overlap w.r.t. vechLk
c    
c
c DSL     (output) DOUBLE PRECISION vector, dimension (N*(N+1)/2)
c         Derivative of overlap w.r.t. vechLl
c    
c
c DTK     (output) DOUBLE PRECISION vector, dimension (N*(N+1)/2)
c         Derivative of kinetic energy w.r.t. vechLk
c
c
c DTL     (output) DOUBLE PRECISION vector, dimension (N*(N+1)/2)
c         Derivative of kinetic energy w.r.t. vechLl
c    
c
c DVK     (output) DOUBLE PRECISION vector, dimension (N*(N+1)/2)
c         Derivative of potential energy w.r.t. vechLk
c    
c
c DVL     (output) DOUBLE PRECISION vector, dimension (N*(N+1)/2)
c         Derivative of potential energy w.r.t. vechLl
c              
c LK LL AK AL INVAKL MAK MAL IAAL IAAK IALK IALL TMP1 TMP2 PLL PIA 
c         are all NxN matrices used for
c         intermediate stages of computing the matrix elements
c
c TMPV1 TMPV2  are temporary vectors for storing Ln*Rn*vec[Akl^-1]
c         this term is used in the potential energy integral and derivative
c         also used for temporary storing various vech vectors
c
c
c =====================================================================
c
c     ..Parameters..
      DOUBLE PRECISION    ONE,    ZERO,   SIX,    TWO,    C2OSPI,
     &                    TWLV,   THREHV 
      PARAMETER         ( ONE = 1.0D+0, ZERO = 0.0D+0, SIX = 6.0D+0,
     &                    TWO = 2.0D+0, C2OSPI = 1.12837916709551D+0,
     &                    TWLV = 1.2D+1, THREHV = 1.5D+0 )
c     ..                                                                      
c     ..External Functions..

c     ..
c     ..External Subroutines..
      EXTERNAL            DTRMM, DPOTRF
c     ..
c     ..Intrinsic Functions..
      INTRINSIC           DSQRT, DABS
c     ..
c     ..Local Scalars..
      INTEGER             I, J, K, NN, INDX, INFO
c     NOT USED IN NEW DV CODE
c     INTEGER             L, C, R
      DOUBLE PRECISION    DETLK, DETLL, DETAKL, SUM
      
c     ..
c     ..Executable Statements..
c
c**********************************************************************
c
c     Build all of the intermediate matrices and products
c     LK, LL, AK, AL, INVAKL, MAK, MAL, IAAK, IAAL, IALK, IALL
c     Also, DETLK, DETLL, DETAKL for the overlap matrix element
c  
c**********************************************************************
c
c     NN = N(N+1)/2 is used for indixing over the derivative vech terms
c
      NN = N*(N+ONE)/TWO
c
c     initilaize LK and LL
c                                      

      DO J=1,N
          DO I=1,N
              LK(I,J)=ZERO
              LL(I,J)=ZERO
          END DO
      END DO
c
c     build the lower triangular matrices LK and LL from
c     the elements in VLK and VLL
c      
        
      DO J=1,N
          DO I=J,N
              INDX = I + (J-1) * (2*N-J)/2
              LK(I,J) = VLK( INDX )
              LL(I,J) = VLL( INDX )
          END DO
      END DO
c
c     Apply the transformation P to LL --- PLL = P'*LL
c
      DO J=1,N
          DO I=1,N
              PLL(I,J) = ZERO
              DO K=1,N
                  PLL(I,J) = PLL(I,J) + P(K,I) * LL(K,J)
              END DO
          END DO
      END DO
     
c
c     Build AK=LK*LL' and AL=P'*LL*LL'*P
c
      DO J=1,N
          DO I=1,N
              AK(I,J) = ZERO
              AL(I,J) = ZERO
              DO K=1,N
                  AK(I,J) = AK(I,J) + LK(I,K) * LK(J,K)
                  AL(I,J) = AL(I,J) + PLL(I,K) * PLL(J,K)
              END DO
          END DO
      END DO
 
c
c     compute AKL THE UPPER TRIANGLE OF AKL=AK+AL IS IN INVAKL
c     (we store the upper part of AKL for convience in the Cholesky
c     decomposition code below)
c
      DO J=1,N
          DO I=J,N
              INVAKL(J,I) = AK(I,J) + AL(I,J)
          END DO
      END DO     
      
c     compute det(LK) and det(LL)
      DETLK = ONE
      DETLL = ONE
      DO I=1,N
          DETLK = DETLK * LK(I,I)
          DETLL = DETLL * LL(I,I)
      END DO

c
c     compute the Cholesky decomposition of Akl (INVAKL)
c     This code is adapted from Numerical Recipes  choldc.for
c     det(AKL) is computed during the decomposition and stored in 
c     DETAKL
c
      DETAKL = ONE
      DO I=1,N
          DO J=I,N
              SUM = INVAKL(I,J)
              DO K=I-1,1,-1
                  SUM = SUM - INVAKL(I,K) * INVAKL(J,K)
              END DO
              IF(I.EQ.J)THEN
                  IF(SUM.LE.ZERO)PAUSE 'AKL IS NOT POSITIVE DEFINITE'
                  INVAKL(I,I) = DSQRT(SUM)
                  DETAKL = DETAKL * SUM
              ELSE
                  INVAKL(J,I) = SUM / INVAKL(I,I)
              END IF
        END DO
      END DO
c
c     INVAKL now contains the cholesky decomp,L, of AKL, where AKL=L*L'
c     DETAKL now contains det(AKL)
c
c     compute the inverse of AKL using the LAPACK subroutine
c
      CALL DPOTRI( 'L', N, INVAKL, NMAX, INFO )
c
c     the call above gives the lower triangle of INVAKL
c     to simplify the matrix products below we complete 
c     the symetric structure of INAKL.
c
      DO I=1,N-1
          DO J=I+1,N
              INVAKL(I,J) = INVAKL(J,I)
          END DO
      END DO
c
c     compute the matrix products MAK=MASS*AK, MAL=MASS*AL, 
c     IAAK=INVAKL*AK, IAAL=INVAKL*AL, 
c     IALK=INVAKL*LK, IALL=INVAKL*LL, PIA=P*INVAKL       
c
      DO J=1,N
          DO I=1,N
              MAK(I,J) = ZERO
              MAL(I,J) = ZERO
              IAAK(I,J) = ZERO
              IAAL(I,J) = ZERO
              IALK(I,J) = ZERO
              IALL(I,J) = ZERO
              PIA(I,J) = ZERO
              DO K=1,N
                  MAK(I,J) = MAK(I,J) + MASS(I,K)*AK(K,J)
                  MAL(I,J) = MAL(I,J) + MASS(I,K)*AL(K,J)
                  IAAK(I,J) = IAAK(I,J) + INVAKL(I,K)*AK(K,J)
                  IAAL(I,J) = IAAL(I,J) + INVAKL(I,K)*AL(K,J)
                  IALK(I,J) = IALK(I,J) + INVAKL(I,K)*LK(K,J)
                  IALL(I,J) = IALL(I,J) + INVAKL(I,K)*PLL(K,J)
                  PIA(I,J) = PIA(I,J) + P(I,K)*INVAKL(K,J)
              END DO
          END DO
      END DO    
c*********************************************************************
c     Now compute the integrals SKL TKL and VKL
c*********************************************************************
c
c     compute SKL = 2^(3*N/2) * ( abs(det(Lk)*det(Ll))/det(Akl) )^(3/2)
c
      SKL = DSQRT( TWO**(3*N) ) * 
     &        DSQRT( ( DABS( DETLK * DETLL ) / DETAKL )**3 )
      
c
c     compute TKL = SIX * SKL * trace( MASS * AK * INVAKL *AL )
c                 = SIX * SKL * trace( MAK*IAAL )
c
c
c     sum over all elements of Hadamard(MAK,IAAL') to get the trace
c
      TKL = ZERO
      DO J=1,N
          DO I=1,N
              TKL = TKL + MAK(I,J) * IAAL(J,I)
          END DO
      END DO
      TKL = SIX * SKL * TKL 
c  
c     Compute VKL = 2/sqrt(Pi) * SKL * vech(Q')'*(Ln*Rn * vec(Akl^-1))^[-1/2]
c     where vech(Q') is contained in CHRG and Rn*vec(Akl^-1) is 
c     evaluated explicitly and computed in a loop C2OSPI=2/sqrt(pi)
c     Put Rn*vec(Akl^-1) in DVK for use in the derivative   
c                 
      VKL = ZERO
      DO J=1,N
          DO I=J,N
              INDX = I + (J-1) * (2*N-J)/2
              IF(I .NE. J) THEN
                  DVK(INDX)=INVAKL(I,I)+INVAKL(J,J)-TWO*INVAKL(I,J)
                  VKL = VKL + CHRG( INDX ) *                            
     &            ONE/DSQRT( DVK( INDX ) )
              ELSE
                  DVK(INDX) = INVAKL(I,J)
                  VKL = VKL + CHRG( INDX ) * ONE/DSQRT( INVAKL(I,J) )
              END IF
          END DO
      END DO
      VKL = C2OSPI * SKL * VKL
c
c************************************************************************
c     Now compute the derivatives
c************************************************************************
c
c     store vech[(LK^-1)' -2*INVAKL*LK]     
c              = vech[1/diag[LK] - 2 * IALK]  in DSK
c     and vech[1/diag(LL) - 2 * PIA * PLL] in DSL
c     these vectors are needed in all of the derivatives
c     and are the main components of the SKL derivatives 
c    
c     first compute TMP1 = PIA * PLL
c
      DO J=1,N
          DO I=1,N
              TMP1(I,J) = ZERO
              DO K=1,N
                  TMP1(I,J) = TMP1(I,J) + PIA(I,K) * PLL(K,J)
              END DO
          END DO
      END DO
c
c
      DO J=1,N
          DO I=J,N
              INDX = I + (J-1) * (2*N-J)/2
              IF(I .NE. J) THEN
                  DSK(INDX) = -TWO * IALK(I,J)
                  DSL(INDX) = -TWO * TMP1(I,J)
              ELSE
                  DSK(INDX) = ONE/LK(I,I) - TWO * IALK(I,I)
                  DSL(INDX) = ONE/LL(I,I) - TWO * TMP1(I,I)
              END IF
          END DO
      END DO
c
c     Kinetic energy derivatives
c     DTK = 3/2*TKL*vech[1/diag[LK] - 2 * IALK] +
c           12*SKL*vech[IAAL*MAL*IALK]
c 
c     First put IAAL*MAL in TMP1 (and IAAK*MAK in TMP2)
      DO J=1,N
          DO I=1,N
              TMP1(I,J) = ZERO
              TMP2(I,J) = ZERO
              DO K=1,N
                  TMP1(I,J) = TMP1(I,J) + IAAL(I,K)*MAL(K,J)
                  TMP2(I,J) = TMP2(I,J) + IAAK(I,K)*MAK(K,J)
              END DO
          END DO
      END DO
c    
c     put P*IAAK*MAK in LL
c
      DO J=1,N
          DO I=1,N
              LL(I,J) = ZERO
              DO K=1,N
                  LL(I,J) = LL(I,J) + P(I,K) * TMP2(K,J)
              END DO
          END DO
      END DO              
c     
c     put vech[IAAL*MAL*IALK] in TMPV1 (vech[P*IAAK*MAK*IALL] in TMPV2)
      DO J=1,N
          DO I=J,N
              INDX = I + (J-1) * (2*N-J)/2
              TMPV1(INDX) = ZERO
              TMPV2(INDX) = ZERO
              DO K=1,N
                  TMPV1(INDX) = TMPV1(INDX) + TMP1(I,K)*IALK(K,J)
                  TMPV2(INDX) = TMPV2(INDX) + LL(I,K)*IALL(K,J)
              END DO
          END DO
      END DO    
c
c     Now build DTK and DTL
c     DTK = 3/2*TKL*vech[1/diag[LK] - 2 * IALK] +
c           12*SKL*vech[IAAL*MAL*IALK]     
c
      DO I=1,NN
          DTK(I) = THREHV * TKL * DSK(I) + TWLV * SKL * TMPV1(I)
          DTL(I) = THREHV * TKL * DSL(I) + TWLV * SKL * TMPV2(I)
      END DO                                               

**** TRY A NEW METHOD FOR DVK, DVL
** THE OLD CODE IS COMMENTED OUT WITH *'S 
C
C     BUILD THE MATRIX D (IN TMP1)
C     WHERE D(I,I)=SUM(K=1,N)[VECH[Q]*(VECH[B])^[-3/2]]
C     AND   D(I,J)=-VECH[Q]*VECH[B]
C     VECH[Q] IN CHRG, VECH[B] IN DVK
      DO J=1,N
          DO I=J,N
              TMP1(I,J) = ZERO
              IF (I .NE. J) THEN
                  INDX = I + (J-1)*(2*N-J)/2
                  TMP1(I,J) = - CHRG(INDX) * 
     &                ONE/DSQRT( DVK(INDX) * DVK(INDX) * DVK(INDX) )
              ELSE
                  DO K=1,(I-1)
                      INDX = I + (K-1)*(2*N-K)/2
                      TMP1(I,I) = TMP1(I,I) + CHRG(INDX) *
     &                 ONE/DSQRT( DVK(INDX) * DVK(INDX) * DVK(INDX) )
                  END DO
                  DO K=I,N
                      INDX = K + (I-1)*(2*N-I)/2
                      TMP1(I,I) = TMP1(I,I) + CHRG(INDX) *
     &                 ONE/DSQRT( DVK(INDX) * DVK(INDX) * DVK(INDX) )
                  END DO
              END IF
          END DO
      END DO
C     COMPLETE UPPER TRIANGLE OF D (TMP1)
      DO I=1,(N-1)
          DO J=(I+1),N
              TMP1(I,J) = TMP1(J,I)
          END DO
      END DO
C     DO THE (FULL) MATRIX MULTIPLY INVAKL*TMP1 IN TMP2
C         ******* AND P*INVAKL*TMP1 IN LL ********
      DO J=1,N
          DO I=1,N
              TMP2(I,J) = ZERO
              LL(I,J) = ZERO
              DO K=1,N
                  TMP2(I,J) = TMP2(I,J) + INVAKL(I,K)*TMP1(K,J)
                  LL(I,J) = LL(I,J) + PIA(I,K)*TMP1(K,J)
              END DO
          END DO
      END DO   
C     PUT VECH[INVAKL*D*INVAKL*LK] in TMPV1 
C         VECH[P*INVAKL*D*INVAKL*P'*LL in TMPV2
      DO J=1,N
          DO I=J,N
              INDX = I + (J-1) * (2*N-J)/2
              TMPV1(INDX) = ZERO
              TMPV2(INDX) = ZERO
              DO K=1,N
                  TMPV1(INDX) = TMPV1(INDX) + TMP2(I,K)*IALK(K,J)
                  TMPV2(INDX) = TMPV2(INDX) + LL(I,K)*IALL(K,J)
              END DO
          END DO
      END DO                                          
c
c     Now add (3/2)*VKL*vech[(LK^-1)' -2*INVAKL*LK] to the last result
c
      DO I=1,NN
          DVK(I) = THREHV * VKL * DSK(I) + C2OSPI * SKL * TMPV1(I)
          DVL(I) = THREHV * VKL * DSL(I) + C2OSPI * SKL * TMPV2(I)
      END DO  
          
*c
*c     Potential energy derivatives
*c     DVK = 3/2 VKL*vech[1/diag[LK] - 2 * IALK] + 2/sqrt(pi) * SKL
*c     (vechQ'[*](Ln*Rn*vec[INVAKL])^[3/2])'*Ln*Rn*kron[IAALk,INVAKL]Ln' 
*c
*c     First put (vechQ'[*](Ln*Rn*vec[INVAKL])^[3/2])' in DVK
*c     recall that DVK contains Ln*Rn*vec[INVAKL]     
*      DO J=1,N
*          DO I=J,N
*              INDX = I + (J-1) * (2*N-J)/2
*              DVK(INDX) = CHRG(INDX) * 
*     &            ONE/DSQRT( DVK(INDX) * DVK(INDX) * DVK(INDX) )
*          END DO
*      END DO
*c
*c     Now compute 2/sqrt(pi) * SKL * 
*c     (vechQ'[*](Ln*Rn*vec[INVAKL])^[3/2])'*Ln*Rn*kron[IAALk,INVAKL]Ln'
*c     and put it in TMPV1 (and TMPV2)
*c
*      DO K=1,N
*          DO L=K,N
*              
*              C = L + (K-1) * (2*N - K)/2
*              TMPV1(C) = ZERO
*              TMPV2(C) = ZERO
*              
*              DO I=1,N
*                  DO J=I,N
*                  
*                  R = J + (I-1) * (2*N - I)/2
*                  IF (I .NE. J) THEN
*                      TMPV1(C)=TMPV1(C)+DVK(R)*(IALK(I,K) * INVAKL(I,L)
*     &                                        + IALK(J,K) * INVAKL(J,L) 
*     &                                        - IALK(I,K) * INVAKL(J,L)
*     &                                        - IALK(J,K) * INVAKL(I,L))
*     
*                      TMPV2(C)=TMPV2(C)+DVK(R)*(IALL(I,K) * INVAKL(I,L)
*     &                                        + IALL(J,K) * INVAKL(J,L) 
*     &                                        - IALL(I,K) * INVAKL(J,L)
*     &                                        - IALL(J,K) * INVAKL(I,L))
*                  ELSE
*                      TMPV1(C)=TMPV1(C)+DVK(R)*(IALK(I,K) * INVAKL(I,L))
*                      TMPV2(C)=TMPV2(C)+DVK(R)*(IALL(I,K) * INVAKL(I,L))
*                  END IF
*                  
*                  END DO
*              END DO
*          END DO
*      END DO
*c
*c     Now add (3/2)*VKL*vech[(LK^-1)' -2*INVAKL*LK] to the last result
*c
*      DO I=1,NN
*          DVK(I) = THREHV * VKL * DSK(I) + C2OSPI * SKL * TMPV1(I)
*          DVL(I) = THREHV * VKL * DSL(I) + C2OSPI * SKL * TMPV2(I)
*      END DO
      
      
c
c     compute DSK = (3/2)*SKL*vech[(LK^-1)' -2*INVAKL*LK]     
c                 = 3/2 * SKL * vech[1/diag[LK] - 2 * IALK]
c
      DO I=1,NN
          DSK(I) = THREHV * SKL * DSK(I) 
          DSL(I) = THREHV * SKL * DSL(I)
      END DO

c ****THIS IS THE CODE THAT MAKES THIS FILE DIFFERENT FROM MELKL.FOR
c
c     Add DSL to DSK, DTL to DTK and DVL to DVK
c     
      DO I=1,NN
          DSK(I) = DSK(I) + DSL(I) 
          DTK(I) = DTK(I) + DTL(I)
          DVK(I) = DVK(I) + DVL(I)
      END DO
c
c     End of MELKK
c
      END

