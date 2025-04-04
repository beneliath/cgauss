



      SUBROUTINE EIGN(N,EIG)                                            
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 IPOSV,IVPOS,IORD
C                  ****  VERSION  3/20/73  ****                         
C        *****  DESCRIPTION                               *****         
C        *****   THIS IS A GIVENS-HOUSHOLDER ROUTINE      *****         
C        *****   RECEIVES A MATRIX-A-OF WHICH ONLY LOWER  *****         
C        *****   TRANGLE IS USED AND THIS IS DESTROYED    *****         
C        *****   EIGENVALUES ARE RETURNED IN ALGEBRAIC    *****         
C        *****   ASCENDING ORDER IN ARRAY-EIG-            *****         
C        *****   THE EIGENVECTORS ARE RETURNED STORED BY  *****         
C        *****   COLUMNS IN VEC.                          *****         
C        *****         PARAMETERS PASSED                  *****         
C        *****   RHO=TOLERANCE ON FINAL OFF-DIAG ELEMENTS *****         
C        *****   N  IS THE SIZE OF THE MATRIX             *****         
C        *****   OF THAT IN QCPE PROGRAM CINDO BY DOBOSH  *****         
C        *****  THIS ROUTINE IS SLIGHTLY MODIFIED VERSION *****         
C                                                                       
      DIMENSION W(500),EIG(500),GAMMA(500),BETA(500),P(500),
     *          Q(500),BETASQ(500) 
      COMMON /EIGMAT/A(500,500),VEC(500,500),IPOSV(500),
     *               IVPOS(500),IORD(500)   
      EQUIVALENCE (IPOSV,GAMMA),(IVPOS,BETA),(IORD,BETASQ),             
     *  (IVPOS,P),(IVPOS,Q)                                             
      RHO=1.0D-10                                                       
      RHOSQ=RHO*RHO                                                     
      IF(N.EQ.0) GO TO 560                                              
      SHIFT = 0.D+00                                                    
      N1=N-1                                                            
      N2=N-2                                                            
      GAMMA(1)=A(1,1)                                                   
      IF(N2.EQ.0) GO TO 270                                             
      IF(N2.LT.0)GO TO 280                                              
      DO 260 NR=1,N2                                                    
      B=A(NR+1,NR)                                                      
      S=0.D+00                                                          
      DO 130 I=NR,N2                                                    
  130 S=S+A(I+2,NR)**2                                                  
C        *****  PREPARE FOR POSSIBLE BYPASS OF TRANSFORMATION  *****    
      A(NR+1,NR)=0.D+00                                                 
      IF(S.LE.0) GO TO 250                                              
      S=S+B*B                                                           
      SGN=+1.D+00                                                       
      IF(B.GE.0) GO TO 160                                              
      SGN = -1.D+00                                                     
160   SQRTS=SQRT(S)                                                    
      D=SGN/(SQRTS+SQRTS)                                               
      TEMP=SQRT(5.D-01+B*D)                                            
      W(NR)=TEMP                                                        
      A(NR+1,NR)=TEMP                                                   
      D=D/TEMP                                                          
      B=-SGN*SQRTS                                                      
C        *****  D IS FACTOR OF PROPORTIONALITY.           *****         
C        *****  COMPUTE AND SAVE W VECTOR                               
      DO 170 I=NR,N2                                                    
      TEMP=D*A(I+2,NR)                                                  
      W(I+1)=TEMP                                                       
  170 A(I+2,NR)=TEMP                                                    
C        *****  PREMULTIPLY VECTOR-W-BY MATRIX -A- TO OBTAIN *****      
C        *****  -P- VECTOR.  SIMULTANEOUSLY ACCUMULATE DOT   *****      
C        *****  PRODUCT WP,(THE SCALAR K).                  *****       
      WTAW = 0.D+00                                                     
      DO 220 I=NR,N1                                                    
      SUM = 0.D+00                                                      
      DO 180 J=NR,I                                                     
  180 SUM=SUM+A(I+1,J+1)*W(J)                                           
      I1=I+1                                                            
      IF((N1-I1).LT.0) GO TO 210                                        
      DO 200 J=I1,N1                                                    
  200 SUM=SUM+A(J+1,I+1)*W(J)                                           
  210 P(I)=SUM                                                          
  220 WTAW=WTAW+SUM*W(I)                                                
C        *****  P VECTOR AND SCALAR K NOW STORED.NEXT     *****         
C        ***** COMPUTE -Q- VECTOR                         *****         
      DO 230 I=NR,N1                                                    
  230 Q(I)=P(I)-WTAW*W(I)                                               
C        *****  NOW FORM PAP MATRIX, REQUIRED PART        *****         
      DO 240 J=NR,N1                                                    
      QJ=Q(J)                                                           
      WJ=W(J)                                                           
      DO 240 I=J,N1                                                     
  240 A(I+1,J+1)=A(I+1,J+1)-2.D+00*(W(I)*QJ+WJ*Q(I))                    
  250 BETA(NR)=B                                                        
      BETASQ(NR)=B*B                                                    
  260 GAMMA(NR+1)=A(NR+1,NR+1)                                          
  270 B=A(N,N-1)                                                        
      BETA(N-1)=B                                                       
      BETASQ(N-1)=B*B                                                   
      GAMMA(N)=A(N,N)                                                   
  280 BETASQ(N)=0.D+00                                                  
C        *****  ADJOIN AN IDENTITY MATRIX TO BE           *****         
C        *****  POSTMULTIPLIED BY ROTATIONS               *****         
      DO 300 I=1,N                                                      
      DO 290 J=1,N                                                      
  290 VEC(I,J)=0.D+00                                                   
  300 VEC(I,I)=1.D+00                                                   
      M=N                                                               
      SUM=0.D+00                                                        
      NPAS=1                                                            
      GO TO 400                                                         
  310 SUM=SUM+SHIFT                                                     
      COSA=1.D+00                                                       
      G=GAMMA(1)-SHIFT                                                  
      PP=G                                                              
      PPBS=PP*PP+BETASQ(1)                                              
      PPBR=SQRT(PPBS)                                                  
      DO 370 J=1,M                                                      
      COSAP=COSA                                                        
      IF(PPBS.NE.0) GO TO 320                                           
      SINA = 0.D+00                                                     
      SINA2=0.D+00                                                      
      COSA=1.D+00                                                       
      GO TO 350                                                         
  320 SINA=BETA(J)/PPBR                                                 
      SINA2=BETASQ(J)/PPBS                                              
      COSA=PP/PPBR                                                      
C        *****  POSTMULTIPLY IDENTITY BY P-TRANSPOSE      *****         
      NT=J+NPAS                                                         
      IF(NT.LE.N) GO TO 330                                             
      NT = N                                                            
  330 DO 340 I=1,NT                                                     
      TEMP=COSA*VEC(I,J)+SINA*VEC(I,J+1)                                
      VEC(I,J+1)=-SINA*VEC(I,J)+COSA*VEC(I,J+1)                         
  340 VEC(I,J)=TEMP                                                     
 350  DIA=GAMMA(J+1)-SHIFT                                              
      U=SINA2*(G+DIA)                                                   
      GAMMA(J)=G+U                                                      
      G=DIA-U                                                           
      PP=DIA*COSA-SINA*COSAP*BETA(J)                                    
      IF(J.NE.M) GO TO 360                                              
      BETA(J)=SINA*PP                                                   
      BETASQ(J)=SINA2*PP*PP                                             
      GO TO 380                                                         
  360 PPBS=PP*PP+BETASQ(J+1)                                            
      PPBR=SQRT(PPBS)                                                  
      BETA(J)=SINA*PPBR                                                 
  370 BETASQ(J)=SINA2*PPBS                                              
 380  GAMMA(M+1)=G                                                      
C        *****  TEST FOR CONVERGENCE OF LAST DIAGONAL ELEMENT  *****    
      NPAS=NPAS+1                                                       
      IF(BETASQ(M).GT.RHOSQ) GO TO 410                                  
  390 EIG(M+1)=GAMMA(M+1)+SUM                                           
  400 BETA(M)=0.D+00                                                    
      BETASQ(M)=0.D+00                                                  
      M=M-1                                                             
      IF(M.EQ.0) GO TO 430                                              
      IF(BETASQ(M).LE.RHOSQ) GO TO 390                                  
C        *****  TAKE ROOT OF CORNER 2 BY 2 NEAREST TO     *****         
C        *****  LOEWR DIAGONAL IN VALUE AS ESTIMATE OF    *****         
C        *****  EIGENVALUE TO USE FOR SHIFT               *****         
  410 A2=GAMMA(M+1)                                                     
      R2=5.D-01*A2                                                      
      R1=5.D-01*GAMMA(M)                                                
      R12=R1+R2                                                         
      DIF=R1-R2                                                         
      TEMP=SQRT(DIF*DIF+BETASQ(M))                                     
      R1=R12+TEMP                                                       
      R2=R12-TEMP                                                       
      DIF=ABS(A2-R1)-ABS(A2-R2)                                       
      IF(DIF.LT.0) GO TO 420                                            
      SHIFT=R2                                                          
      GO TO 310                                                         
  420 SHIFT=R1                                                          
      GO TO 310                                                         
  430 EIG(1)=GAMMA(1)+SUM                                               
C        *****  INITIALIZE AUXILIARY TABLES REQUIRED FOR  *****         
C        *****  REARANGING THE VECTORS                    *****         
      DO 440 J=1,N                                                      
      IPOSV(J)=J                                                        
      IVPOS(J)=J                                                        
  440 IORD(J) = J                                                       
C        ***** USE A TRANSPOSITION SORT TO ORDER THE EIGENVALUES  ***** 
      M=N                                                               
      GO TO 470                                                         
  450 DO 460 J=1,M                                                      
      IF(EIG(J).LE.EIG(J+1)) GO TO 460                                  
      TEMP=EIG(J)                                                       
      EIG(J)=EIG(J+1)                                                   
      EIG(J+1)=TEMP                                                     
      ITEMP=IORD(J)                                                     
      IORD(J)=IORD(J+1)                                                 
      IORD(J+1)=ITEMP                                                   
  460 CONTINUE                                                          
  470 M=M-1                                                             
      IF(M.NE.0) GO TO 450                                              
      IF(N1.EQ.0) GO TO 500                                             
      DO 490 L=1,N1                                                     
      NV=IORD(L)                                                        
      NP=IPOSV(NV)                                                      
      IF(NP.EQ.L) GO TO 490                                             
      LV=IVPOS(L)                                                       
      IVPOS(NP)=LV                                                      
      IPOSV(LV)=NP                                                      
      DO 480 I=1,N                                                      
      TEMP=VEC(I,L)                                                     
      VEC(I,L)=VEC(I,NP)                                                
  480 VEC(I,NP) = TEMP                                                  
  490 CONTINUE                                                          
  500 CONTINUE                                                          
C        *****  BACK TRANSFORM THE VECTORS OF THE TRIPLE  *****         
C        *****  DIAGONAL MATRIX                                         
      DO 550 NRR=1,N                                                    
      K=N1                                                              
  510 K=K-1                                                             
      IF(K.LE.0) GO TO 540                                              
      SUM = 0.D+00                                                      
      DO 520 I=K,N1                                                     
  520 SUM=SUM+VEC(I+1,NRR)*A(I+1,K)                                     
      SUM=SUM+SUM                                                       
      DO 530 I=K,N1                                                     
  530 VEC(I+1,NRR)=VEC(I+1,NRR)-SUM*A(I+1,K)                            
      GO TO 510                                                         
  540 CONTINUE                                                          
  550 CONTINUE                                                          
  560 CONTINUE                                                          
      RETURN                                                            
      END    
