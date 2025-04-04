SUBROUTINE dELECREP(OVERLAP,D1,D2,D3,D4,D5,				&
					AZ,A1,BZ,A2,			&
                    CZ,A3,DZ,A4,A5,A6,		&
					DER1,DER2,DER3,DER4,DER5)                       

	USE scalars_to_be_shared

    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	IMPLICIT INTEGER (I-N)


!*** FOUR CENTER ELECTRON REPULSION INTEGRAL ***
!-----------------------------------------------------

      ALK1=A1+A2                                                       
      ALK2=A3+A4                                                       
      BLK12=A5+A6
      ALKXX=ALK2+BLK12

	  QALK = ALK1+ALK2

      DETAB=ALK1*ALK2+QALK*BLK12
	
	  RRZ1=(A1*AZ+A2*BZ)/ALK1
      RRZ2=(A3*CZ+A4*DZ)/ALK2

	  DFT0=(RRZ1-RRZ2)

      DD=DFT0*DFT0


!-----------------------------------------------------

	Z=ALK1*ALK1*ALK2*ALK2*DD/(QALK*DETAB) 
	PHIF=((((((C6*Z+C7)*Z+C8)*Z+C9)*Z+C10)*Z+C11)*Z+ONE)
	THETAF=(((((C1*Z+C2)*Z+C3)*Z+C4)*Z+C5)*Z+ONE)
	OKL=DSQRT(DETAB)
	PKL=DSQRT(ONE/(PI*QALK))
	ANKL=OKL*PKL
	AMKL=TWO*OVERLAP*ANKL


IF (Z .LE. 16.3578D+00) THEN
	FF=DSQRT(THETAF/PHIF)

!DERDA1:
!======================================================
	  DDETAB=ALKXX

      t7 = PHIF*PHIF
      t11 = Z*Z*Z*Z*Z
      t12 = ALK1*ALK1
      t13 = ALK2*ALK2
      t14 = t12*t13
      t16 = ONE/ALK1
      t20 = A1*AZ+A2*BZ
      t30 = t20*t16-RRZ2
      t33 = ONE/QALK
      t34 = ONE/DETAB
      t38 = t30*t30
      t40 = QALK*QALK
      t51 = DETAB*DETAB

      t56 = t14*TWO*(AZ*t16-t20/t12)*t30*t33*t34-t14*t38*t34/t40+TWO*ALK1*t13*t38*t33*t34-	&
	  ALKXX*t14*t38/t51*t33

      t60 = Z*Z*Z*Z
      t61 = t60*t56
      t64 = Z*Z*Z
      t65 = t64*t56
      t68 = Z*Z
      t69 = t68*t56
      t72 = Z*t56
      t98 = PI*(ALK2+ALK1)

	  t100 = FF*TWO
	  t101 = THREE*C3
	  t102 = FOUR*C2
	  t103 = FIVE*C1
	  t104 = ONE/PHIF
	  t105 = THREE*C9
	  t106 = FOUR*C8
	  t107 = FIVE*C7
	  t108 = SIX*C6
	  t109 = THETAF/t7
	  t110 = AMKL*DSQRT(PHIF)*HALF/DSQRT(THETAF)
	  t111 = PKL*HALF/DSQRT(DETAB)
	  t112 = TWO*C4
	  t113 = TWO*C10
	  t114 = OKL*HALF/t98/DSQRT(t98)*PI

      DER1 = t110*(-t109*(t108*t11*t56+t107*t61+t106*t65+t105*t69+t113*t72+C11*t56)+	&
	  t104*(t103*t61+t102*t65+t101*t69+t112*t72+C5*t56))+t100*(OVERLAP*(-t114+t111*		&
	  DDETAB)+ANKL*D1)


!DERDA2:
!======================================================
      u16 = 1/ALK2
      u20 = A3*CZ+A4*DZ
      u30 = t20/ALK1-u20*u16
      u38 = u30**TWO
      u39 = u38*t34
      u48 = ALK1+BLK12
	  u50 = t34*t33

      u55 = t14*TWO*(-CZ*u16+u20/t13)*u30*u50-t14*u39/t40+TWO*t12*ALK2*u39*t33-t14*u48*	&
	  u38/t51*t33

      u60 = t60*u55
      u64 = t64*u55
      u68 = t68*u55
      u71 = Z*u55

      DER2 = t110*(-t109*(t108*t11*u55+t107*u60+t106*u64+t105*u68+t113*u71+C11*u55)+	&
	  t104*(t103*u60+t102*u64+t101*u68+t112*u71+C5*u55))+t100*(OVERLAP*(-t114+t111*		&
	  u48)+ANKL*D2)


!DERDA5:
!======================================================
	  v26 = (t20/ALK1-u20/ALK2)
      v27 = v26*v26
      v29 = ONE/t51
      v31 = t14*v27*v29
      v50 = t13*v27*v29

      DER3 = t110*(-t109*(-t108*t11*v31-v31*t107*t60-t106*t64*v31-t105*t68*v31-t113*Z*	&
	  v31-C11*t12*v50)+t104*(-t103*t60*v31-t102*t64*v31-t101*t68*v31-t112*Z*v31-C5*t12*	&
	  v50))+t100*(OVERLAP*t111*QALK+ANKL*D3)


!DERDAZ:
!======================================================
      w33 = ALK1*t13*v26*u50
      w30 = TWO*A1
	  w38 = t60*w30
      w44 = t64*w30
      w50 = t68*w30
      w55 = Z*A1
      w59 = A1*ALK1
      w63 = t13*v26*u50

	  w100 = t100*ANKL
	  w101 = C5*TWO
	  w102 = FOUR*C4
	  w103 = C11*TWO
	  w104 = FOUR*C10
	  w105 = t108*t11*TWO

      DER4 = t110*(-t109*(w105*A1*w33+t107*w38*w33+t106*w44*w33+t105*w50*w33+w104*w55*w33+	&
	  w103*w59*w63)+t104*(t103*w38*w33+t102*w44*w33+t101*w50*w33+w102*w55*w33+w101*w59*		&
	  w63))+w100*D4


!DERDBZ:
!======================================================
      x33 = t12*ALK2*v26*u50
      x30 = TWO*A3
	  x38 = t60*x30
      x44 = t64*x30
      x50 = t68*x30
      x55 = Z*A3
      x59 = A3*t12
      x63 = ALK2*v26*u50

      DER5 = t110*(-t109*(-w105*A3*x33-t107*x38*x33-t106*x44*x33-t105*x50*x33-w104*x55*x33-	&
	  w103*x59*x63)+t104*(-t103*x38*x33-t102*x44*x33-t101*x50*x33-w102*x55*x33-w101*x59*	&
	  x63))+w100*D5

!============================================================================================================
!============================================================================================================
!============================================================================================================

ELSE

		    R=TWO*Z
			R3=PI/TWO
			FF=DSQRT(R3/R)

 
!DERDA1:
!======================================================

      tt2 = Z*Z
      tt7 = ALK1*ALK1
      tt8 = ALK2*ALK2
      tt9 = tt7*tt8
      tt11 = ONE/ALK1
      tt15 = A1*AZ+A2*BZ
	  tt20 = A3*CZ+A4*DZ
      tt25 = tt15*tt11-tt20/ALK2
      tt27 = ONE/DETAB
      tt29 = ONE/QALK
      tt33 = tt25*tt25
      tt34 = tt33*tt27
      tt35 = QALK*QALK
      tt45 = DETAB*DETAB
      tt59 = PI*(ALK2+ALK1)

	  tt100 = AMKL*QUARTER*DSQRT(PI/tt2/Z)
	  tt101 = OKL*HALF/tt59/DSQRT(tt59)*PI

      DER1 = -tt100*(tt9*TWO*(AZ*tt11-tt15/tt7)*tt25*tt27*tt29-tt9*tt34/tt35+TWO*ALK1*tt8*	&
	  tt34*tt29-ALKXX*tt7*tt8*tt33/tt45*tt29)+t100*(OVERLAP*(-tt101+t111*DDETAB)+ANKL*D1)

!DERDA2:
!======================================================

      uu11 = ONE/ALK2
      uu29 = ONE/QALK
      uu35 = QALK*QALK
      uu43 = ALK1+BLK12
	  uu50 = PKL*HALF/DSQRT(DETAB)
	  uu51 = FF*TWO

      DER2 = -tt100*(tt9*TWO*(-CZ*uu11+tt20/tt8)*tt25*tt27*uu29-tt9*tt34/uu35+TWO*tt7*	&
	  ALK2*tt34*uu29-tt9*uu43*tt33/tt45*uu29)+uu51*(OVERLAP*(-tt101+uu50*uu43)+ANKL*D2)


!DERDA5:
!======================================================

      DER3 = tt100*tt7*tt8*tt33/tt45+uu51*(OVERLAP*uu50*QALK+ANKL*D3)


!DERDAZ:
!======================================================

	ww100 = uu51*ANKL
	ww101 = DFT0/DETAB/QALK
	ww102 = tt100*TWO
	ww103 = ww102*tt8*ww101

      DER4 = -ww103*A1*ALK1+ww100*D4


!DERDBZ:
!======================================================

      DER5 = ww103*A3*ALK2+ww100*D5


END IF

END SUBROUTINE dELECREP
