C***************************************************************************
      SUBROUTINE SEARCH(LEAVE,GAUSS)                                    
      IMPLICIT REAL*8 (A-H,O-Z)
C                                                                       
C                                                                       
C     SUBROUTINE LOCATES A LOCAL MINIMUM ALONG A GIVEN ONE-DIMENSIONAL  
C     SEARCH DIRECTION.                                                 
C     F(ALF) = FUNC(X+ALF*DX) FIND OPTIMUM ALF VALUE.                   
C     LOOKMAX = MAXIMUM NUMBER OF FUNCTION EVALUATIONS BEFORE SEARCH IS 
C               TERMINATED                                              
C     ENTER SUBROUTINE WITH FOLLOWING                                   
C     N     = DIMENSIONALITY OF X SPACE                                 
C     X     = INITIAL POSITION (IN N DIMENSIONAL PARAMETER SPACE)       
C     FUNC  = F(0) , I.E. FUNCTION VALUE AT INITIAL X.                  
C     DX    = SEARCH DIRECTIONS(IN N-DIMENSIONAL SPACE)                 
C     AMIN  = ESTIMATED INITIAL STEP SIZE                               
C     NOISE = ACCURACY OF FUNC VALUES. VARIATIONS SMALLER THAN NOISE    
C             ARE NOT CONSIDERED TO BE SIGNIFICANT.                     
C     EXACT = .FALSE.  SUPPRESSES FINAL CHECK OF QUADRATIC FIT.         
C                                                                       
C     EXIT FROM SUBROUTINE WITH FOLLOWING                               
C     X     = CALCULATED MINIMUM                                        
C     FUNC  = F(AMIN) ,I.E. FUNCTION VALUE AT FINAL X.                  
C     CURVE = CALCULATED CURVATURE ALON SEARCH DIRECTIONS               
C     LEAVE = 0 IF NORMAL EXIT                                          
C           = 2 IF EXCESSIV NUMBER OF FUNCTIONS EVALUATIONS             
C                                                                       
C     N,NMAX,FUNC,AND X ARE SHARED WITH SUBROUTINE GAUSS.               
C     ALL VARIABLES IN COMMON T1 ARE SHARED WITH PROGRAM CALLING SEARCH 
C                                                                       
      COMMON/T1/ N,NMAX, FFF,X(500), DX(500), CURVE,AMIN ,NOISE,EXACT    
      COMMON/MIN/ FMIN, LOOKS
      DIMENSION XNOT(500), F(9), SET(9)                                  
      EQUIVALENCE (F(1),F1),(F(3),F3),(F(4),F4),(F(5),F5),(F(6),F6),    
     1  (F(7),F7),(F(9),F9),(SET(1),SET1),                              
     2  (SET(2),SET2),(SET(3),SET3),(SET(4),SET4),(SET(5),SET5),        
     3  (SET(6),SET6),(SET(7),SET7),(SET(8),SET8),(SET(9),SET9)         
      LOGICAL BACK, SET, SET1,SET2,SET3,SET4,SET5,SET6,SET7,SET8,SET9   
     1,EXACT,HOLD                                                       
      LOGICAL AIF
      REAL*8 NOISE                                                        
      LOOKS=0                                                           
      LEAVE=0                                                           
      LOOKMAX=100                                                        
      DO 5 I=1,N                                                        
    5 XNOT(I)=X(I)                                                      
      ALF=AMIN                                                          
      AMIN=0.D+00                                                       
      FNOT = Fff                                                       
      FMIN=Fff                                                       
      FNEXT=1.0D+30                                                     
      IR=4                                                              
   10 LOOKS=LOOKS+1                                                     
      IF(LOOKS.LE.LOOKMAX) GO TO 20                                     
      LEAVE=2                                                           
      open (unit=6,access='append',form='formatted')
      WRITE(6,13) LOOKS                                                 
      close(6)
   13 FORMAT(/' EXCESSIVE NUMBER OF FUNCTION EVALUATIONS =',I3/)        
C                                                                       
C     RETURN FROM SUBROUTINE                                            
C                                                                       
   15 DO 16 I=1,N                                                       
   16 X(I)=XNOT(I)+AMIN*DX(I)                                           
      Fff=FMIN                                                         
      RETURN                                                            
C                                                                       
C     EVALUATE F(ALF)=FUNC(XNOT+ALF*DX)                                 
C                                                                       
   20 DO 25 I=1,N                                                       
   25 X(I)=XNOT(I)+ALF*DX(I)                                            
      INOUT=0                                                           
      CALL GAUSS(INOUT)                                                 
C                                                                       
C     INOUT = 7 INDICATES PATHOLOGICAL X VALUE FOR GAUSS CALCULATION.   
C                                                                       
      IF(INOUT.EQ.7) GO TO 305                                          
      HOLD=.FALSE.                                                      
      IF(Fff.GE.FMIN) GO TO 28                                         
      FNEXT=FMIN                                                        
      FMIN=Fff                                                         
      ANEXT=AMIN                                                        
      AMIN=ALF                                                          
      GO TO 29                                                          
   28 IF(Fff.GE.FNEXT) GO TO 29                                        
      FNEXT=Fff                                                        
      ANEXT=ALF                                                         
   27 FORMAT(' LOOKS =',I3,' ALF =',F10.5, ' FFF =',F15.9, ' AMIN =',  
     1 F10.5,' FMIN =',F12.6,' ANEXT =',F10.5,' FNEXT =',F12.6)         
   29 CONTINUE
C      open (unit=6,access='append',form='formatted')
C      WRITE(6,27) LOOKS, ALF, FFF, AMIN, FMIN, ANEXT, FNEXT            
C      close (6)
      GO TO (29,29,29,150,180,260,320) IR                               
C                                                                       
C     F1,F3,F5,F7,F9, ARE FUNC VALUES AT EVENTLY SPACED POINTS WITH     
C     SEPARATION DELTA. SET(I)=TRUE INDICATES THAT F(I) HAS BEEN        
C     EVALUATED. DELTA MAY INCREASE OR DECREASE BY FACTORS OF 2. 0.     
C     F VALUES SHIFT LEFT OR RIGHT SO THAT F5 IS ALWAYS LOWEST.         
C                                                                       
  150 DELTA=ABS(AMIN-ANEXT)                                            
      BACK=.FALSE.                                                      
      DO 155 I=1,9                                                      
  155 SET(I)=.FALSE.                                                    
      F5=FMIN                                                           
      SET5=.TRUE.                                                       
      K=3                                                               
      IF(AMIN.LT.ANEXT) K=7                                             
      F(K)=FNEXT                                                        
      SET(K)=.TRUE.                                                     
  160 IF(SET7) GO TO 230                                                
      IF(SET3) GO TO 200                                                
      IF(BACK) GO TO 250                                                
  170 ALF=AMIN+DELTA                                                    
      IR=5                                                              
      GO TO 10                                                          
  180 F7=FFF                                                           
      SET7=.TRUE.                                                       
      IF(ALF.NE.AMIN) GO TO 160                                         
C                                                                       
C     SHIFT LEFT                                                        
C                                                                       
      DO 190 I=1,7                                                      
      SET(I)=SET(I+2)                                                   
  190 F(I)=F(I+2)                                                       
      SET8=.FALSE.                                                      
      SET9=.FALSE.                                                      
      GO TO 160                                                         
  200 CONTINUE
      IF(SET1) GOTO 210
      IF(.NOT.SET1) GO TO 170
C                                                                       
C     EXTRAPOLATE USING F1,F3, AND F5                                   
C     ICREASE STEP SIZE IF F9 IS PREDICTED TO BE LOWER THAN F7.         
C     DECREASE STEP SIZE IF F4 IS PREDICTED TO BE LOWER THAN F5.        
C                                                                       
  210 BACK=4.D+00*F3.LT.(F1+3.D+00*F5)                                  
      IF(SET9.OR.HOLD) GO TO 215                                        
      IF(5.D+00*F3.GT.(2.D+00*F1+3.D+00*F5)) GO TO 220
  215 IF(SET4) GO TO 170                                                
      AIF=6.D+00*F3.GT.(F1+5.D+00*F5)                                   
      IF(AIF) GO TO 170
      IF(.NOT.AIF) GO TO 305
  220 DELTA=DELTA+DELTA                                                 
      F4=F3                                                             
      F3=F1                                                             
      F6=F7                                                             
      F7=F9                                                             
      SET4=SET3                                                         
      SET3=SET1                                                         
      SET6=SET7                                                         
      SET7=SET9                                                         
      SET1=.FALSE.                                                      
      SET2=.FALSE.                                                      
      SET8=.FALSE.                                                      
      SET9=.FALSE.                                                      
      GO TO 160                                                         
  230 IF(SET3) GO TO 300                                                
      IF(SET9) GO TO 240
      IF(.NOT.SET9) GO TO 250
C                                                                       
C     EXTRAPOLATE USING F5, F7, AND F9.                                 
C     INCREASE ISTEP SIZE IF F1 IS PREDICTED TO B LOWER THAN F3         
C     DECREASE STEP SIZE IF F6 IS PREDICTED TO BE LOWER THAN F5.        
C                                                                       
  240 BACK=4.D+00*F7.GT.(F9+3.D+00*F5)                                  
      IF(SET1.OR.HOLD) GO TO 245                                        
      IF(5.D+00*F7.GT.(2.D+00*F9+3.D+00*F5)) GO TO 220
  245 IF(SET6) GO TO 250                                                
      AIF=6.D+00*F7.GT.(F9+5.D+00*F5)                                   
      IF(AIF) GO TO 250
      IF(.NOT.AIF) GO TO 305
  250 ALF=AMIN-DELTA                                                    
      IR=6                                                              
      GO TO 10                                                          
  260 F3=Fff                                                           
      SET3=.TRUE.                                                       
      IF(ALF.NE.AMIN) GO TO 160                                         
      DO 270 I=1,7                                                      
      SET(10-I)=SET(8-I)                                                
  270 F(10-I)=F(8-I)                                                    
      SET1=.FALSE.                                                      
      SET2=.FALSE.                                                      
      GO TO 160                                                         
C                                                                       
C     INTERPOLATE USING F3, F5, AND F7.                                 
C     EITHER IMMEDIATELY ATTEMPT QUADRATIC FIT OR FIRST REDUCE STEP SIZE
C     IF AMIN=ZERO STEP SIZE  MAY BE REDUCED BY A FACTOR OF 2,4,8,OR 16 
C     IN OTHER CASES STEP SIZE IS REDUCED BY A FACTOR OF 2.             
C                                                                       
  300 DUM=F3                                                            
      IF(F7.GT.F3) DUM=F7                                               
      SPREAD=DUM-F5                                                     
      TOL=(FNOT-FMIN)*5.D-01                                            
      IF(TOL.LT.1.D+03*NOISE) TOL=1.D+03*NOISE                          
      IF(.NOT.EXACT) TOL=TOL+TOL                                        
      DUM=F3-F5-F5+F7                                                   
      CURVE=DUM/(2.D+00*DELTA*DELTA)                                    
      DUM=(F3-F7)/DUM                                                   
      BACK=DUM.LT.0.D+00                                                
      IF(SPREAD.LT.TOL) GO TO 310                                       
      IF(ABS(DUM).GT.75.D-02.OR.AMIN.NE.0.D+00) GO TO 305              
      SET1=SET4                                                         
      F1=F4                                                             
      SET9=SET6                                                         
      F9=F6                                                             
      SET2=.FALSE.                                                      
      SET3=.FALSE.                                                      
      SET4=.FALSE.                                                      
      SET6=.FALSE.                                                      
      SET7=.FALSE.                                                      
      SET8=.FALSE.                                                      
      DELTA=25.D-02*DELTA                                               
      IF(ABS(DUM).GT.375.D-03) GO TO 160                               
      SET1=.FALSE.                                                      
      SET9=.FALSE.                                                      
      DELTA=5.D-1*DELTA                                                 
      IF(ABS(DUM).GT.1875.D-04) GO TO 160                              
      DELTA=5.D-1*DELTA                                                 
      GO TO 160                                                         
C                                                                       
C     CUT STEP SIZE IN HALF                                             
C                                                                       
  305 DELTA=5.D-01*DELTA                                                
      HOLD=.TRUE.                                                       
      SET1=SET3                                                         
      SET3=SET4                                                         
      SET9=SET7                                                         
      SET7=SET6                                                         
      SET2=.FALSE.                                                      
      SET4=.FALSE.                                                      
      SET6=.FALSE.                                                      
      SET8=.FALSE.                                                      
      F1=F3                                                             
      F3=F4                                                             
      F9=F7                                                             
      F7=F6                                                             
      GO TO 160                                                         
C                                                                       
C     ATTEMPT QUADRATIC FIT                                             
C                                                                       
  310 ALF=AMIN+DUM*DELTA*5.D-1                                          
      FQUAD=F5-(F3-F7)*DUM*125.D-03                                     
      IR=7                                                              
      GO TO 10                                                          
  320 IF(.NOT.EXACT) GO TO 15                                           
      TOL=(FNOT-FMIN)*25.D-03                                           
      IF(TOL.LT.NOISE) TOL=NOISE                                        
      IF(ABS(Fff-FQUAD).LT.TOL) GO TO 15                              
      AIF=ALF.EQ.AMIN.OR.ALF.EQ.ANEXT
      IF(AIF) GO TO 150
      IF(.NOT.AIF) GO TO 305
      END                                                               
