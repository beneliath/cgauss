      program aaaa
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/T1/ ND,NMAX,FFF,RR(20),DR(20),CURVE,ALPHA,NOISE,EXACT,
     *           ALPH(20),RNOT(20),V(20,20),CURV(20),RESTAR,JSTART,
     *           FNOT,TIMLIM,TOLF,TOLR
      EXTERNAL fu
      LOGICAL EXACT,RESTAR
      real*8 noise
C
      tolf = 1.d-04
      TOLR = 1.d-04
      NOISE = 1.d-08
      nd = 2
      nmax = 20
      KKKK=-1
      rr(1) = 2.d+00
      rr(2) = 2.d+00
      fff = 1.d+03
      CALL TRUDGE(KKKK,fu)
      STOP
      END
