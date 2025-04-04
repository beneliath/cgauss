      SUBROUTINE GRAD(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,
     1                RZ2,B12,S,H,CO,G)	
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(MMAX=50)
      DIMENSION 
     1         AA1(MMAX),RX1(MMAX),RY1(MMAX),RZ1(MMAX),
     2         AA2(MMAX),RX2(MMAX),RY2(MMAX),RZ2(MMAX),
     3         B12(MMAX)
      DIMENSION S(MMAX,MMAX),H(MMAX,MMAX),CO(MMAX) 
C ---  FOR grdient
      DIMENSION GSA1(MMAX,MMAX),GSA2(MMAX,MMAX),GSB(MMAX,MMAX),
     1          GSR1(MMAX,MMAX),GSR2(MMAX,MMAX),DELTA(100,100),
     2          GHA1(MMAX,MMAX),GHA2(MMAX,MMAX),GHB(MMAX,MMAX),
     1          GHR1(MMAX,MMAX),GHR2(MMAX,MMAX),G(5*MMAX)
      COMMON /NORM/ SNORM(Mmax,Mmax)
      COMMON /VALUE/ EMIN
C
C  -- MMAX - MAX. LENGTH OF EXPANSION --
C
      DATA ONE,TWO,FOUR,SIX,EIGHT,HALF/1.,2.,4.,6.,8.,.5/               
C --- DATA ---
      DO 20 IL=1,100
      DO 20 IK=1,100
      DELTA(IL,IK)=0.0
20    CONTINUE
      DO 11 II=1,100
      DELTA(II,II)=1.0 
11    CONTINUE
C
C      DO 10 L=1,M
C      CALL GEMPRT(L,
C     1            AA1(L),RX1(L),RY1(L),RZ1(L),
C     2            AA2(L),RX2(L),RY2(L),RZ2(L),
C     3                                       B12(L))
C10    CONTINUE
      CALL GEOM(EX1,EY1,EZ1,EX2,EY2,EZ2)
C ----------------------------------------------------
      DO 100 L=1,M
      DO 200 K=1,M
C
C---- OVERLAP INEGRAL
C
      CALL OVERLP(RX1(L),RY1(L),RZ1(L),AA1(L),
     1            RX1(K),RY1(K),RZ1(K),AA1(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            RX2(K),RY2(K),RZ2(K),AA2(K),
     4            B12(L),B12(K),OV1,KLM)  
      CALL GROVER(L,K,RZ1(L),AA1(L),RZ1(K),AA1(K),                      
     1              RZ2(L),AA2(L),RZ2(K),AA2(K),B12(L),B12(K),
     2              DS1A1,DS1A2,DS1B,DS1R1,DS1R2) 
      CALL OVERLP(RX1(L),RY1(L),RZ1(L),AA1(L),
     1            RX2(K),RY2(K),RZ2(K),AA2(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            RX1(K),RY1(K),RZ1(K),AA1(K),
     4            B12(L),B12(K),OV2,KLM)  
      CALL GROVER(L,K,RZ1(L),AA1(L),RZ2(K),AA2(K),                      
     1              RZ2(L),AA2(L),RZ1(K),AA1(K),B12(L),B12(K),
     2              DS2A1,DS2A2,DS2B,DS2R1,DS2R2) 
      S(L,K)=OV1+OV2
      GSA1(L,K)=DS1A1+DS2A1
      GSA2(L,K)=DS1A2+DS2A2
      GSB(L,K) =DS1B +DS2B
      GSR1(L,K)=DS1R1+DS2R1
      GSR2(L,K)=DS1R2+DS2R2
C
C-----KINETIC ENERGY OF ELECTRON 1
C
      CALL GRKIN(OV1,L,K,
     1            RZ1(L),AA1(L),RZ1(K),AA1(K),                      
     2            RZ2(L),AA2(L),RZ2(K),AA2(K),
     3            B12(L),B12(K),DT1A1,DT1A2,DT1B,DT1R1,DT1R2)  
      CALL GRKIN(OV2,L,K, 
     1            RZ1(L),AA1(L),RZ2(K),AA2(K),                      
     2            RZ2(L),AA2(L),RZ1(K),AA1(K),
     3            B12(L),B12(K),DT2A1,DT2A2,DT2B,DT2R1,DT2R2)  
C
C----KINETIC ENERGY OF ELECTRON 2
C
      CALL GRKIN(OV1,L,K,
     1            RZ2(L),AA2(L),RZ2(K),AA2(K),                      
     2	          RZ1(L),AA1(L),RZ1(K),AA1(K),
     3            B12(L),B12(K),DT3A1,DT3A2,DT3B,DT3R1,DT3R2)  
      CALL GRKIN(OV2,L,K,
     1            RZ2(L),AA2(L),RZ1(K),AA1(K),                      
     2            RZ1(L),AA1(L),RZ2(K),AA2(K),
     3            B12(L),B12(K),DT4A1,DT4A2,DT4B,DT4R1,DT4R2) 
C
C---- N-E FOR ELECTRON 1
C
      CALL GRNUE(L,K,RZ1(L),AA1(L),
     1           RZ1(K),AA1(K),RZ2(L),AA2(L),
     2           RZ2(K),AA2(K),B12(L),B12(K),
     3           EZ1,DEN1A1,DEN1A2,DEN1B,DEN1R1,DEN1R2)                       
      CALL GRNUE(L,K,RZ1(L),AA1(L),
     1           RZ2(K),AA2(K),RZ2(L),AA2(L),
     2           RZ1(K),AA1(K),B12(L),B12(K),
     3           EZ1,DEN2A1,DEN2A2,DEN2B,DEN2R1,DEN2R2)                       
      CALL GRNUE(L,K,RZ1(L),AA1(L),
     1           RZ1(K),AA1(K),RZ2(L),AA2(L),
     2           RZ2(K),AA2(K),B12(L),B12(K),
     3           EZ2,DEN3A1,DEN3A2,DEN3B,DEN3R1,DEN3R2)
      CALL GRNUE(L,K,RZ1(L),AA1(L),
     1           RZ2(K),AA2(K),RZ2(L),AA2(L),
     2           RZ1(K),AA1(K),B12(L),B12(K),
     3           EZ2,DEN4A1,DEN4A2,DEN4B,DEN4R1,DEN4R2)                       
C
C---- N-E FOR ELECTRON 2
C
      CALL GRNUE(L,K,RZ2(L),AA2(L),
     1           RZ2(K),AA2(K),RZ1(L),AA1(L),
     2           RZ1(K),AA1(K),B12(L),B12(K),
     3           EZ1,DEN5A1,DEN5A2,DEN5B,DEN5R1,DEN5R2)                       
      CALL GRNUE(L,K,RZ2(L),AA2(L),
     1           RZ1(K),AA1(K),RZ1(L),AA1(L),
     2           RZ2(K),AA2(K),B12(L),B12(K),
     3           EZ1,DEN6A1,DEN6A2,DEN6B,DEN6R1,DEN6R2)                       
      CALL GRNUE(L,K,RZ2(L),AA2(L),
     1           RZ2(K),AA2(K),RZ1(L),AA1(L),
     2           RZ1(K),AA1(K),B12(L),B12(K),
     3           EZ2,DEN7A1,DEN7A2,DEN7B,DEN7R1,DEN7R2)                       
      CALL GRNUE(L,K,RZ2(L),AA2(L),
     1           RZ1(K),AA1(K),RZ1(L),AA1(L),
     2           RZ2(K),AA2(K),B12(L),B12(K),
     3           EZ2,DEN8A1,DEN8A2,DEN8B,DEN8R1,DEN8R2)                       
C 
      CALL GRELE(L,K,RZ1(L),AA1(L),
     1           RZ1(K),AA1(K),RZ2(L),AA2(L),
     2           RZ2(K),AA2(K),B12(L),B12(K),
     3           DER1A1,DER1A2,DER1B,DER1R1,DER1R2)  
      CALL GRELE(L,K,RZ1(L),AA1(L),
     1           RZ2(K),AA2(K),RZ2(L),AA2(L),
     2           RZ1(K),AA1(K),B12(L),B12(K),
     3           DER2A1,DER2A2,DER2B,DER2R1,DER2R2)  
      DTA1=(DT1A1+DT2A1+DT3A1+DT4A1)/TWO
      DTA2=(DT1A2+DT2A2+DT3A2+DT4A2)/TWO
      DTB= (DT1B+ DT2B+ DT3B+ DT4B)/TWO
      DTR1=(DT1R1+DT2R1+DT3R1+DT4R1)/TWO
      DTR2=(DT1R2+DT2R2+DT3R2+DT4R2)/TWO
      DENA1=DEN1A1+DEN2A1+DEN3A1+DEN4A1+DEN5A1+DEN6A1+DEN7A1+DEN8A1
      DENA2=DEN1A2+DEN2A2+DEN3A2+DEN4A2+DEN5A2+DEN6A2+DEN7A2+DEN8A2
      DENB =DEN1B +DEN2B +DEN3B +DEN4B +DEN5B +DEN6B +DEN7B +DEN8B
      DENR1=DEN1R1+DEN2R1+DEN3R1+DEN4R1+DEN5R1+DEN6R1+DEN7R1+DEN8R1
      DENR2=DEN1R2+DEN2R2+DEN3R2+DEN4R2+DEN5R2+DEN6R2+DEN7R2+DEN8R2
      DERA1=DER1A1+DER2A1
      DERA2=DER1A2+DER2A2
      DERB =DER1B +DER2B
      DERR1=DER1R1+DER2R1
      DERR2=DER1R2+DER2R2
      GHA1(L,K)=DERA1+DTA1-DENA1
      GHA2(L,K)=DERA2+DTA2-DENA2
      GHB(L,K) =DERB +DTB -DENB
      GHR1(L,K)=DERR1+DTR1-DENR1
      GHR2(L,K)=DERR2+DTR2-DENR2
 200  CONTINUE 
 100  CONTINUE 
C     CALL MATTAB(S,M)
      CALL TAB(S,M,M,MMAX,MMAX)
      CALL TAB(H,M,M,MMAX,MMAX)
      call ENERGY(M,H,S,EMIN,CO,KIM)
      DO 400 K=1,M
	GG1=0.0
 	GG2=0.0
	GG3=0.0
	GG4=0.0
	GG5=0.0
      DO 410 L=1,M
      CKL=CO(L)*CO(K)/(DSQRT(SNORM(K,K))*DSQRT(SNORM(L,L)) )	
      WW=2.0-DELTA(L,K)
	GG1=GG1+WW*CKL*(GHA1(L,K)-EMIN*GSA1(L,K))
	GG2=GG2+WW*CKL*(GHA2(L,K)-EMIN*GSA2(L,K))
	GG3=GG3+WW*CKL*(GHB(L,K) -EMIN*GSB(L,K) ) 
	GG4=GG4+WW*CKL*(GHR1(L,K)-EMIN*GSR1(L,K))
	GG5=GG5+WW*CKL*(GHR2(L,K)-EMIN*GSR2(L,K))
410	CONTINUE
	G(5*K -4) =GG1
        G(5*K -3) =GG2
        G(5*K -2) =GG3
        G(5*K -1) =GG4
	G(5*K ) =GG5
400	CONTINUE
C
      END