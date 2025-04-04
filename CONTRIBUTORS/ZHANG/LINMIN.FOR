      SUBROUTINE linmin(p,xi,n,fret)
C      COMMON /PASSIT/ M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,S,H,CO
      INTEGER n,NMAX
      REAL*8 fret,p(n),xi(n),TOL
      PARAMETER (NMAX=50,TOL=1.0D-04)
CU    USES brent,f1dim,mnbrak
      INTEGER j,ncom
      REAL*8 ax,bx,fa,fb,fx,xmin,xx,pcom(NMAX),xicom(NMAX),brent
      COMMON /f1com/ pcom,xicom,ncom
      EXTERNAL f1dim
      ncom=n
      do 11 j=1,n
        pcom(j)=p(j)
        xicom(j)=xi(j)
11    continue
      ax=0.0D+00  
      xx=1.0D+00
      call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
      fret=brent(ax,xx,bx,f1dim,TOL,xmin)
      do 12 j=1,n
        xi(j)=xmin*xi(j)
        p(j)=p(j)+xi(j)
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software $|2jaU-1j.
