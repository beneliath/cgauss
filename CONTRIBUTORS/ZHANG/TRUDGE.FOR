

      SUBROUTINE TRUDGE(K,FUNC)                                         
      IMPLICIT REAL*8(A-H,O-Z)
C                                                                       
C     PROGRAM MINIMIZES A FUNCTION OF N VARIABLES BY A MODFIED POWELL   
C     METHOD OF SEARCH ALONG CONJUGATE DIRECTIONS. COLUMNS OF V MATRIX  
C     ARE CURRENT SEARCH DIRECTIONS.COLUMNS 1 TO K ARE (APROXIMATELY)   
C     CONJUGATE.                                                        
C                                                                       
C     WHEN SUBROUTINE IS ENTERED WITH RESTART=.TRUE. THEN VALUES OF K   
C     AND ALL QUANTITIES IN COMMON T1 SHOULD HAVE BEEN SET PREVIOUSLY   
C     BY THE CALLING PROGRAM.                                           
C     WHEN TRUDGE IS CALLING WITH K NEGATIVE THEN ONLY REQUIRED VALUES  
C     ARE TIMLIM,TOLF,TOLR,NOISE,N,AND R.(NOISE NEEDED FOR SEARCH)      
C     IN THIS CASE INITIAL VALUES OF K,JSTART,RESTART,ALPH,V ARE IGNORED
C     IF TRUDGE IS CALLED WITH NON-NEGATIVE K AND RESTA-T=.FALSE. THEN  
C     VALUES OF TIMLIM,TOLF,TOLR,NOISE,N,R,AND V ARE REQUIRED.          
C     IN THIS CASE K IS RESET TO ZERO AFTER FIRST PASS.                 
C                                                                       
C     EXTERNAL FUNCTION FUNC(INOUT)                                     
C     INOUT=  2   FUNCTION READS MINIMAL STARTING DATA AND EVALUATES F(R
C     INOUT=  1   FUNCTION READS ALL RESTART DATA AND EVALUATES F(R)    
C     INOUT=  0   FUNCTION EVALUATES F(R)                               
C     INOUT= -1   FUNCTION PUNCHES DATA FOR FUTURE CALL FUNC(1)         
C     INOUT= -2   FUNCTION EVALUATES F(R) AND PRINTS OUTPUT.            
C                                                                       
      COMMON/T1/ N,NMAX,F,R(500),DR(500),CURVE,ALPHA,NOISE,EXACT,
     *           ALPH(500)
     1,RNOT(500), V(500,500), CURV(500),RESTART, 
     1 JSTART, FNOT, TIMLIM,      
     2 TOLF,TOLR 
      COMMON /PASSIT/ M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,EVEC,S,H,CO
      COMMON/EIGMAT/G(500,500),Q(500,500),SCRAP1(500),
     1 SCRAP2(500),
     *              SCRAP3(500)  
      DIMENSION P(500),U(500),EVEC(50)                                         
      LOGICAL EXACT, RESTART                                            
      REAL*8 NOISE
      EXTERNAL FUNC                                                     
      NMAX=500
	III=1
      
C      WRITE(6,*)TOLF,TOLR,NOISE
C      pause
                                                                
      IF(N.GT.NMAX) STOP 'N GT 180'
      FN=N                                                              
      INOUT=0                                                           
      IF(K.LT.0) GO TO 15                                               
      IF(.NOT.RESTART) CALL FUNC(INOUT)                                 
      DO 10 I=1,N                                                       
   10 IF(ALPH(I).LE.0.D+00) ALPH(I)=1.D-01                              
      IF(JSTART.LE.0.OR.JSTART.GT.N) JSTART=N                           
      J=JSTART                                                          
      RESTART=RESTART.AND.K.NE.N                                        
      GO TO 100                                                         
C                                                                       
C     IF K NEGATIVE SET INITIAL V MATRIX TO IDENTITY,STANDARD VALUES.   
C                                                                       
   15 K=N                                                               
      RESTART=.FALSE.                                                   
      CALL FUNC(INOUT)                                                  
      DO 30 I=1,N                                                       
      ALPH(I)=1.D-01                                                    
      DO 20 J=1,N                                                       
   20 V(I,J)=0.D+00                                                     
      CURV(I)=0.D+00                                                    
   30 V(I,I)=1.D+00                                                     
      J=N                                                               
C                                                                       
C     SEARCH ALONG J-TH DIRECTION FOR LOCAL MINIMUM                     
C                                                                       
  100 DO 105 I=1,N                                                      
  105 DR(I)=V(I,J)                                                      
      EXACT=J.LE.K                                                      
      ALPHA=ALPH(J)                                                     
      ALFSET=ALPHA                                                      
      CALL SEARCH(LEAVE,FUNC)                                          
C      open (unit=6,access='append',form='formatted')
C      WRITE(6,120) CURVE,ALPHA                                          
C      close (6)
  120 FORMAT(17H  NEW CURVATURE =,E13.5,10H   ALPHA =,F12.7)            
  121 FORMAT(17H  OLD CURVATURE =,E13.5,10H   ALFSET=,F12.7)            
  122 FORMAT(20X,30H ABNORMAL RETURN FROM SEARCH  )                     
      ALPHA=ABS(ALPHA)                                                 
      IF(ALPHA.LT.5.D-01*ALFSET) ALFSET=5.D-01*ALFSET                   
      IF(ALPHA.LT.25.D-02*ALFSET) ALFSET=5.D-1*ALFSET                   
      IF(ALPHA.GT.2.D+00*ALFSET) ALFSET=2.D+00*ALFSET                   
      IF(ALPHA.GT.4.D+00*ALFSET) ALFSET=2.D+00*ALFSET                   
      open (unit=6,access='append',form='formatted')
C      WRITE(6,121) CURV(J),ALFSET                                       
      close (6)
      CURV(J)=CURVE                                                     
      ALPH(J)=ALFSET                                                    
      IF(J.EQ.1) GO TO 175                                              
      J=J-1                                                             
      GO TO 100
C                                                                       
C     CALL FUNC(-1) CAUSES ALL RESTART DATA TO BE OUTPUT,               
C     THIS INCLUDES ALL DATA REQUIRED FOR FUTURE CALL FUNC(1).          
C                                                                       
      open (unit=6,access='append',form='formatted')
C  150 WRITE(6,155) K,TIMLIM,TOLF,TOLR,NOISE,J,FNOT                      
 150  close (6)
  155 FORMAT(I10,F10.0,3E10.1,I10,F20.10)                               
      INOUT=-2                                                          
      CALL FUNC(INOUT)                                                  
      DO 160 I=1,N                                                      
      open (unit=6,access='append',form='formatted')
C  160 WRITE(6,161)I,RNOT(I),CURV(I),ALPH(I)                             
  160 close (6)
  161 FORMAT(I5,3E20.10)                                                
      I1=1                                                              
  165 I2=I1+7                                                           
      IF(I2.GT.N) I2=N                                                  
      DO 166 L=1,N                                                      
      open (unit=6,access='append',form='formatted')
C  166 WRITE(6,168)(V(L,I),I=I1,I2)                                      
  166    close (6)
  168 FORMAT(8F10.7)                                                    
      I1=I2+1                                                           
      IF(I1.LE.N) GO TO 165                                             
      K=0                                                               
      RETURN
C                                                                       
C     END OF SEQUENCE OF SEARCH ALONG ALL DIRECTIONS IN V MATRIX.       
C                                                                       
  175 IF(RESTART) GO TO 185                                             
      DO 180 I=1,N                                                      
  180 RNOT(I)=R(I)                                                      
      RESTART=.TRUE.                                                    
      FNOT=F                                                            
      K=0                                                               
      J=N                                                               
      GO TO 100                                                         
C                                                                       
C     GENERATE NEW CONJUGATE DIRECTION IF NOT CONVERGED.                
C                                                                       
  185 DRDR=0.D+00                                                       
      DO 188 I=1,N                                                      
      DUM=R(I)-RNOT(I)                                                  
      DR(I)=DUM                                                         
  188 DRDR=DUM*DUM+DRDR                                                 
      ALFSET=SQRT(DRDR)                                                
      ALPHA=ALFSET                                                      
      IF(ALPHA.EQ.0.D+00) GO TO 192                                     
      DO 190 I=1,N                                                      
      U(I)=DR(I)/ALFSET                                                 
  190 DR(I)=U(I)                                                        
      EXACT=.TRUE.                                                      
      CALL SEARCH(LEAVE,FUNC)
      open (unit=6,access='append',form='formatted')
      IF(LEAVE.NE.0) WRITE(6,122)                                       
C      WRITE(6,120) CURVE,ALPHA                                          
      close (6)
C      ALPHA=ABS(ALPHA)                                                 
  192 DELF=FNOT-F                                                       
      FNOT=F                                                            
      DRDR=0.D+00                                                       
      DO 195 I=1,N                                                      
      DRDR=(R(I)-RNOT(I))**2+DRDR
  195 RNOT(I)=R(I)                                                      
      DELR=SQRT(DRDR/FN)                                               
      open (unit=6,access='append',form='formatted')
      I=1
      J=1
	WRITE(6,*) 'ITERATION INDEX: ',III
	WRITE(7,*) 'ITERATION INDEX: ',III
      CALL BIOUT('------------------------------------------------------
     1---------------------')
	III=III+1
      CALL BIOUT('     LinCoeff:.   AA1:        AA2:        B12:        
	1RZ1:       RZ2:')
      CALL BIOUT('FNC: =================================================
     1=====================')
 606  IF (J .LE. N/5) THEN
 1111 FORMAT(I3,6D12.4)
          WRITE(6,1111) J,EVEC(J),R(I),R(I+1),R(I+2),R(I+3),R(I+4)
          WRITE(7,1111) J,EVEC(J),R(I),R(I+1),R(I+2),R(I+3),R(I+4)
          I=I+5
          J=J+1
      GOTO 606
      END IF
      CALL BIOUT('------------------------------------------------------
     1---------------------')
      WRITE(6,210) F, DELF, DELR 
      WRITE(7,210) F, DELF, DELR 
      CALL BIOUT(' ')

     
     
     
                                            
      close (6)
  210 FORMAT(4H  E=,F20.15,10H  DELTA E=,E15.9,10H  DELTA R=,E15.9)     
      IF(DELF.GT.TOLF.OR.DELR.GT.TOLR) GO TO 220                        
      open (unit=6,access='append',form='formatted')
      WRITE(6,215)                                                      
      close (6)
  215 FORMAT(//,' CONVERGED')                                           
      GO TO 150                                                         
  220 IF(K.GT.0) GO TO 230                                              
      K=1                                                               
      K1=0                                                              
      DO 225 I=1,N                                                      
  225 R(I)=DR(I)                                                        
      DR(N)=CURVE                                                       
      P(1)=ALPHA                                                        
      GO TO 375                                                         
C                                                                       
C     CALCULATE K+1 Q VECTORS.Q(I) DOT V(J)=KRONIKER DELTA(I,J),J=1,K   
C                                                                       
  230 DO 250 J=1,K                                                      
      DOT=0.D+00                                                        
      DO 240 I=1,N                                                      
  240 DOT=DR(I)*V(I,J)+DOT                                              
      P(J)=DOT                                                          
      DO 250 I=1,N                                                      
  250 DR(I)=DR(I)-DOT*V(I,J)                                            
      DOT=0.D+00                                                        
      DO 260 I=1,N                                                      
  260 DOT=DR(I)*U(I)+DOT                                                
      K1=K                                                              
      K=K1+1                                                            
      RESTART=K.NE.N                                                    
      DO 270 I=1,N                                                      
      DUM=DR(I)/DOT                                                     
      Q(I,K)=DUM                                                        
      DO 270 J=1,K1                                                     
  270 Q(I,J)=V(I,J)-P(J)*DUM                                            
C                                                                       
C     G MATRIX=SUM L=1,K+1 OF Q(L)*CURV(L)*TRANSPOSE Q(L)               
C                                                                       
      DO 280 I=1,N                                                      
      DUM=CURVE*Q(I,K)                                                  
      DO 280 J=1,N                                                      
  280 G(I,J)=DUM*Q(J,K)                                                 
      DO 300 L=1,K1                                                     
      DO 300 I=1,N                                                      
      DUM=CURV(L)*Q(I,L)                                                
      DO 300 J=1,N                                                      
  300 G(I,J)=DUM*Q(J,L)+G(I,J)                                          
      CALL EIGN(N,DR)                                                   
C                                                                       
C     SET P(J) EQUAL TO CORRESPONDING ALPH VALUE.                       
C                                                                       
      DO 350 J=1,K                                                      
      JJ=N+1-J                                                          
      P(J)=ALPH(J)                                                      
      DOT=0.D+00                                                        
      DO 310 I=1,N                                                      
      R(I)=Q(I,JJ)                                                      
  310 DOT=R(I)*V(I,J)+DOT                                               
      DOT=ABS(DOT)                                                     
      IF(DOT.GT.7.D-01) GO TO 350                                       
      LMAX=J                                                            
      DOTMAX=DOT                                                        
      DO 340 L=1,K1                                                     
      IF(L.EQ.J) GO TO 340                                              
      DOT=0.D+00                                                        
      DO 320 I=1,N                                                      
  320 DOT=R(I)*V(I,L)+DOT                                               
      DOT=ABS(DOT)                                                     
      IF(DOT.LT.7.D-01) GO TO 330                                       
      P(J)=ALPH(L)                                                      
      GO TO 350                                                         
  330 IF(DOT.LE.DOTMAX) GO TO 340                                       
      DOTMAX=DOT                                                        
      LMAX=L                                                            
  340 CONTINUE                                                          
      P(J)=ALPH(LMAX)                                                   
      DOT=0.D+00                                                        
      DO 345 I=1,N                                                      
  345 DOT=R(I)*U(I)+DOT                                                 
      IF(ABS(DOT).GT.DOTMAX) P(J)=ALPHA                                
  350 CONTINUE                                                          
C                                                                       
C     SET NEV CONJUGATE VECTORS IN V EXCEPT KTH VECTOR IS SAVED IN R.   
C                                                                       
      DO 370 J=1,K1                                                     
      JJ=N+1-J                                                          
      DO 360 I=1,N                                                      
  360 V(I,J)=Q(I,JJ)                                                    
      ALPH(J)=P(J)                                                      
  370 CURV(J)=DR(JJ)                                                    
C                                                                       
C     ORTHOGONALIZE REMAINING SEARCH DIRECTIONS TO CONJUGATE DIRECTIONS.
C                                                                       
  375 IF(K.EQ.N) GO TO 450                                              
      DO 420 J=K,N                                                      
      DOT=0.D+00                                                        
      DO 380 I=1,N                                                      
  380 DOT=V(I,J)*R(I)+DOT                                               
      DO 390 I=1,N                                                      
  390 V(I,J)=V(I,J)-DOT*R(I)                                            
      IF(K.EQ.1) GO TO 420                                              
      DO 410 L=1,K1                                                     
      DOT=0.D+00                                                        
      DO 400 I=1,N                                                      
  400 DOT=V(I,J)*V(I,L)+DOT                                             
      DO 410 I=1,N                                                      
  410 V(I,J)=V(I,J)-DOT*V(I,L)                                          
  420 CONTINUE                                                          
      CMAX=0.D+00                                                       
      DO 425 L=K,N                                                      
  425 IF(CURV(L).GT.CMAX) CMAX=CURV(L)                                  
C                                                                       
C     CONSTRUCT G MATRIX FOR REMAINING SEARCH DIRECTIONS.               
C                                                                       
      CMAX=8.D+00*CMAX                                                  
      DO 440 I=1,N                                                      
      DO 440 J=1,I                                                      
      DUM=-CMAX*R(I)*R(J)                                               
      DO 430 L=K,N                                                      
  430 DUM=CURV(L)*V(I,L)*V(J,L)+DUM                                     
      IF(K1.LE.0) GO TO 435                                             
      DO 433 L=1,K1                                                     
  433 DUM=-CMAX*V(I,L)*V(J,L)+DUM                                       
  435 G(I,J)=DUM                                                        
  440 G(J,I)=DUM                                                        
C                                                                       
C     SET KTH CONJUGATE DIRECTIONS BACK INTO V MATRIX.                  
C                                                                       
  450 DO 460 I=1,N                                                      
      V(I,K)=R(I)                                                       
  460 R(I)=RNOT(I)                                                      
      ALPH(K)=P(K)                                                      
      CURV(K)=DR(N+1-K)                                                 
      IF(K.EQ.N) GO TO 480                                              
      CALL EIGN(N,DR)                                                   
C                                                                       
C     SET REMAINING SEARCH DIRECTIONS INTO V MATRIX.                    
C                                                                       
      K1=K+1                                                            
      DO 470 J=K1,N                                                     
      CURV(J)=DR(J)                                                     
      DO 470 I=1,N                                                      
  470 V(I,J)=Q(I,J)                                                     
  480 J=N                                                               
      GO TO 100
      END                                                               
