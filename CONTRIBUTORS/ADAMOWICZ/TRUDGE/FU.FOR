
      subroutine fu(i)
      implicit real*8(a-h,o-z)
c ***
      COMMON/T1/ N,NMAX,F, R(20),DR(20),CURVE,ALPHA,NOISE,EXACT,ALPH(20)A    239
     1,RNOT(20), V(20,20), CURV(20),RESTART, JSTART, FNOT, TIMLIM,      A    240
     2 TOLF,TOLR                                                        A    241
c ***
      i = 0
      x = r(1)
      y = r(2)
      f = (x-4.d+00)**2 + (y-10.d+00)**2 + x**2 * y**2
      write (6,*) ' x,y,f',x,y,f
      return
      end      
