      FUNCTION f1dim(x)
      INTEGER NMAX
      REAL*8 f1dim,func,x
      COMMON /f1com/ pcom,xicom,ncom
C      COMMON /PASSIT/ M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,S,H,CO
      PARAMETER (NMAX=50)
CU    USES func
      INTEGER j,ncom
      REAL*8 pcom(NMAX),xicom(NMAX),xt(NMAX)
      do 11 j=1,ncom
        xt(j)=pcom(j)+x*xicom(j)
11    continue
      CALL FUNK(XT,F1DIM)
c      f1dim=func(xt)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software $|2jaU-1j.
