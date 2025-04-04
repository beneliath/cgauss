      PROGRAM OPT
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 NOISE
      COMMON /VALUE/ M,EMIN
      COMMON /GEMPARAM/ A(200),B(100)
      EXTERNAL F
      LOGICAL RESTART,EXACT
      COMMON /T1/ ND,NMAX,FFF,RR(180),DR(180),CURVE,
     1    ALPHA,NOISE,EXACT,ALPH(180),RNOT(180),V(180,180),
     2    CURV(180),RESTART,JSTART,FNOT,TIMLIM,
     3    TOLF,TOLR
      K=-1     
      TOLF=1.D-04
      TOLR=1.D-02
      NOISE=1.0D-06
      FFF=1.D+03
      RESTART=.FALSE.
      ND=2
      RR(1)=3
      RR(2)=4
      WRITE(6,*) RR(1)
      WRITE(6,*) RR(2)
      CALL TRUDGE(K,F)
      WRITE(6,*) RR(1)
      WRITE(6,*) RR(2)
      WRITE(6,*) ' OPT F = ',FFF
      RETURN
      END
