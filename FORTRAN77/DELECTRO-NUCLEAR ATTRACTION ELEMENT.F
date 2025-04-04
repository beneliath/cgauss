SUBROUTINE dENUCLEAR(OVERLAP,				&
			D1,D2,D3,D4,D5,					&
			AZ,A1,BZ,A2,		&
            CZ,A3,DZ,A4,A5,A6,	&
            EZ,SIGNUM,				&
			DEN1,DEN2,DEN3,DEN4,DEN5)                       

	USE scalars_to_be_shared

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	IMPLICIT INTEGER (I-N)


!*** FOUR CENTER NUCLEAR ATTRACTION INTEGRAL ***
!-----------------------------------------------------

      ALK1=A1+A2                                                       
      ALK2=A3+A4                                                       
      BLK12=A5+A6  
      ALKXX=ALK2+BLK12
      DETAB=ALK1*ALK2+(ALK1+ALK2)*BLK12


!--- ONE CENTER NUCLEAR INTEGRAL ---
!-----------------------------------------------------

      RZ1=(A1*AZ+A2*BZ)/ALK1
      RZ2=(A3*CZ+A4*DZ)/ALK2


!--- FOUR-CENTER ELECTRON REPULSION INTEGRAL ---
!-----------------------------------------------------

      RRZ=(ALK1*(ALK2+BLK12)*RZ1+ALK2*BLK12*RZ2)/DETAB
	   

      RE2=(RRZ-EZ)*(RRZ-EZ)

      DELTA2=(ALK2+BLK12)/(FOUR*DETAB)

      RE2=RE2/(FOUR*DELTA2)

!-----------------------------------------------------

	Z=RE2
	PHIF=((((((C6*Z+C7)*Z+C8)*Z+C9)*Z+C10)*Z+C11)*Z+ONE)
	THETAF=(((((C1*Z+C2)*Z+C3)*Z+C4)*Z+C5)*Z+ONE)

	OKL=DSQRT(DETAB)

	QKL=(RRZ-EZ)*(RRZ-EZ)

	RRRN=SIGNUM*(BD/TWO)+RRZ

	A1Z=A1*AZ
	A2Z=A2*BZ
	A3Z=A3*CZ
	A4Z=A4*DZ
	AZSUM=A1Z+A2Z+A3Z+A4Z
	AZSUM1=A1Z+A2Z
	AZSUM2=A3Z+A4Z
	FF=DSQRT(THETAF/PHIF)
	PKL=DSQRT(1.00D+00/(PI*(ALK2+BLK12)))
	ANKL=OKL*PKL
	AMKL=TWO*OVERLAP*ANKL


IF (Z .LE. 16.3578D+00) THEN
	FF=DSQRT(THETAF/PHIF)

!DENDA1:
!======================================================

      t7 = PHIF*PHIF
      t11 = Z*Z*Z*Z*Z
      t14 = Z*Z*Z*Z
      t17 = Z*Z*Z
      t20 = Z*Z
      t33 = DETAB*DETAB
	  t34 = (ALKXX*AZSUM1+BLK12*AZSUM2)/t33
	  t35 = DETAB*TWO

      t41 = QKL+t35*(ALKXX*AZ/DETAB-ALKXX*t34)*RRRN/ALKXX
	  t42 = FF*TWO
	  t43 = THREE*C3
	  t44 = FOUR*C2
	  t45 = FIVE*C1
	  t46 = ONE/PHIF
	  t47 = TWO*C10
	  t48 = THREE*C9
	  t49 = FOUR*C8
	  t50 = FIVE*C7
	  t51 = SIX*C6
	  t52 = THETAF/t7
	  t53 = AMKL*DSQRT(PHIF)*HALF/DSQRT(THETAF)
	  t54 = TWO*C4
	  t55 = (t51*t11+t50*t14+t49*t17+t48*t20+t47*Z+C11)
	  t56 = (t45*t14+t44*t17+t43*t20+t54*Z+C5)
	  t57 = t52*t55
	  t58 = PKL*HALF/DSQRT(DETAB)

      DEN1 = t53*(-t57*t41+t46*t56*t41)+t42*(OVERLAP*t58*ALKXX+ANKL*D1)

!DENDA2:
!======================================================

      u25 = ALK1+BLK12
      u27 = ONE/ALKXX
      u30 = ALKXX*ALKXX
	  u31 = DETAB*QKL/u30
	  u32 = RRRN*u27
	  u33 = QKL*u27

      u49 = u25*u33-u31+t35*((A1Z+A2Z+BLK12*CZ)/DETAB-u25*t34)*u32
      
	  u73 = PI*(ALK2+BLK12)
	  u74 = -OKL*HALF/u73/DSQRT(u73)*PI

      DEN2 = t53*(-t57*u49+t46*t56*u49)+t42*(OVERLAP*(u74+t58*u25)+ANKL*D2)

!DENDA5:
!======================================================

      v25 = ALK1+ALK2

      v47 = v25*u33-u31+t35*(AZSUM/DETAB-v25*t34)*u32

	  v48 = t54*Z
	  v49 = (t45*t14+t44*t17+t43*t20+v48+C5)

      DEN3 = t53*(-t57*v47+t46*v49*v47)+t42*(OVERLAP*(u74+t58*v25)+ANKL*D3)

!DENDAZ:
!======================================================

      w27 = TWO*A1*RRRN
	  w28 = t42*ANKL

      DEN4 = t53*(-t57*w27+t46*t56*w27)+w28*D4

!DENDBZ:
!======================================================

      x30 = A3*BLK12*RRRN/ALKXX

      DEN5 = t53*(-t57*TWO*x30+t46*v49*TWO*x30)+w28*D5

ELSE


!WRITE(6,*)'SPECIAL CASE ERROR FUNCTION HIT'


		    R=TWO*Z
			R3=PI/TWO
			FF=DSQRT(R3/R)


!DENDA1:
!======================================================
	DDETAB=ALKXX

	DZZ=(QKL+DETAB*(TWO*(ALKXX*AZ/DETAB-ALKXX*(ALKXX*AZSUM1+BLK12*AZSUM2)/(DETAB*DETAB))*RRRN)/ALKXX)

	DTHETAF=((FIVE*C1*Z**FOUR+FOUR*C2*Z**THREE+THREE*C3*Z**TWO+TWO*C4*Z+C5)*DZZ)

	DPHIF=((SIX*C6*Z**FIVE+FIVE*C7*Z**FOUR+FOUR*C8*Z**THREE+THREE*C9*Z**TWO+TWO*C10*Z+C11)*DZZ)

	DF0=-QUARTER*DSQRT(PI/(Z*Z*Z))*DZZ

	DOKL=((ONE/(TWO*DSQRT(DETAB)))*DDETAB)

	DNkl=(PKL*DOKL)

	DMKL=(TWO*(OVERLAP*DNKL+ANKL*D1))

	DEN1=AMKL*DF0 + FF*DMKL



!DENDA2:
!======================================================
	DDETAB=(ALK1+BLK12)

	DZZ=((DDETAB)*QKL/ALKXX-DETAB*QKL/(ALKXX*ALKXX)+DETAB*(TWO*((A1Z+A2Z+BLK12*CZ)/DETAB-(ALK1+BLK12)*(ALKXX*AZSUM1+BLK12*AZSUM2)/(DETAB*DETAB))*RRRN)/ALKXX)

	DTHETAF=((FIVE*C1*Z**FOUR+FOUR*C2*Z**THREE+THREE*C3*Z**TWO+TWO*C4*Z+C5)*DZZ)

	DPHIF=((SIX*C6*Z**FIVE+FIVE*C7*Z**FOUR+FOUR*C8*Z**THREE+THREE*C9*Z**TWO+TWO*C10*Z+C11)*DZZ)

	DF0=-QUARTER*DSQRT(PI/(Z*Z*Z))*DZZ

	DPkl=((-ONE/(TWO*((PI*ALK2+PI*BLK12)*DSQRT((PI*ALK2+PI*BLK12)))))*PI)

	DOKL=((ONE/(TWO*DSQRT(DETAB)))*DDETAB)

	DNkl=(OKL*DPKL+PKL*DOKL)

	DMKL=(TWO*(OVERLAP*DNKL+ANKL*D2))

	DEN2=AMKL*DF0 + FF*DMKL


!DENDA5:
!======================================================
	DDETAB=(ALK1+ALK2)

	DZZ=((DDETAB)*QKL/ALKXX-DETAB*QKL/(ALKXX*ALKXX)+DETAB*(TWO*(AZSUM/DETAB-(ALK1+ALK2)*(ALKXX*AZSUM1+BLK12*AZSUM2)/(DETAB*DETAB))*RRRN)/ALKXX)

	DTHETAF=((FIVE*C1*Z**FOUR+FOUR*C2*Z**THREE+THREE*C3*Z**TWO+TWO*C4*Z+C5)*DZZ)

	DPHIF=((SIX*C6*Z**FIVE+FIVE*C7*Z**FOUR+FOUR*C8*Z**THREE+THREE*C9*Z**TWO+TWO*C10*Z+C11)*DZZ)

	DF0=-QUARTER*DSQRT(PI/(Z*Z*Z))*DZZ

	DPkl=((-ONE/(TWO*((PI*ALK2+PI*BLK12)*DSQRT((PI*ALK2+PI*BLK12)))))*PI)

	DOKL=((ONE/(TWO*DSQRT(DETAB)))*DDETAB)

	DNkl=(OKL*DPKL+PKL*DOKL)

	DMKL=(TWO*(OVERLAP*DNKL+ANKL*D3))

	DEN3=AMKL*DF0+FF*DMKL


!DENDAZ:
!======================================================
	DZZ=(TWO*A1*RRRN)

	DTHETAF=((FIVE*C1*Z**FOUR+FOUR*C2*Z**THREE+THREE*C3*Z**TWO+TWO*C4*Z+C5)*DZZ)

	DPHIF=((SIX*C6*Z**FIVE+FIVE*C7*Z**FOUR+FOUR*C8*Z**THREE+THREE*C9*Z**TWO+TWO*C10*Z+C11)*DZZ)

	DF0=-QUARTER*DSQRT(PI/(Z*Z*Z))*DZZ

	DMKL=(TWO*(ANKL*D4))

	DEN4=AMKL*DF0+FF*DMKL


!DENDBZ:
!======================================================
	DZZ=(TWO*A3*BLK12*RRRN/ALKXX)

	DTHETAF=((FIVE*C1*Z**FOUR+FOUR*C2*Z**THREE+THREE*C3*Z**TWO+TWO*C4*Z+C5)*DZZ)

	DPHIF=((SIX*C6*Z**FIVE+FIVE*C7*Z**FOUR+FOUR*C8*Z**THREE+THREE*C9*Z**TWO+TWO*C10*Z+C11)*DZZ)

	DF0=-QUARTER*DSQRT(PI/(Z*Z*Z))*DZZ

	DMKL=(TWO*(ANKL*D5))

	DEN5=AMKL*DF0 + FF*DMKL

END IF

END SUBROUTINE dENUCLEAR
