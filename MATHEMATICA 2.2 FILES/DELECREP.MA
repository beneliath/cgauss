(*^

::[	frontEndVersion = "Microsoft Windows Mathematica Notebook Front End Version 2.2";
	microsoftWindowsStandardFontEncoding;
	fontset = title, "Arial", 24, L0, center, nohscroll, bold;
	fontset = subtitle, "Arial", 18, L0, center, nohscroll, bold;
	fontset = subsubtitle, "Arial", 14, L0, center, nohscroll, bold;
	fontset = section, "Arial", 14, L0, bold, grayBox;
	fontset = subsection, "Arial", 12, L0, bold, blackBox;
	fontset = subsubsection, "Arial", 10, L0, bold, whiteBox;
	fontset = text, "Arial", 12, L0;
	fontset = smalltext, "Arial", 10, L0;
	fontset = input, "Arial", 10, L0, B32768, grayDot;
	fontset = output, "Courier New", 9, L0, nowordwrap, R32768;
	fontset = message, "Courier New", 9, L0, nowordwrap, bold, italic, R65280;
	fontset = print, "Courier New", 8, L0, nowordwrap;
	fontset = info, "Courier New", 10, L0, nowordwrap;
	fontset = postscript, "Courier New", 8, L0, nowordwrap;
	fontset = name, "Arial", 10, L0, nohscroll, italic, B65280;
	fontset = header, "Times New Roman", 10, L0, right, nohscroll;
	fontset = footer, "Times New Roman", 10, L0, right, nohscroll;
	fontset = help, "Arial", 10, L0, nohscroll;
	fontset = clipboard, "Arial", 12, L0, nohscroll;
	fontset = completions, "Arial", 12, L0, nowordwrap, nohscroll;
	fontset = graphics, "Courier New", 10, L0, nowordwrap, nohscroll;
	fontset = special1, "Arial", 12, L0, nowordwrap, nohscroll;
	fontset = special2, "Arial", 12, L0, center, nowordwrap, nohscroll;
	fontset = special3, "Arial", 12, L0, right, nowordwrap, nohscroll;
	fontset = special4, "Arial", 12, L0, nowordwrap, nohscroll;
	fontset = special5, "Arial", 12, L0, nowordwrap, nohscroll;
	fontset = leftheader, "Arial", 12, L0, nowordwrap, nohscroll;
	fontset = leftfooter, "Arial", 12, L0, nowordwrap, nohscroll;
	fontset = reserved1, "Courier New", 10, L0, nowordwrap, nohscroll;]
:[font = subtitle; inactive; nohscroll; center; ]
Derivatives of cgauss electron repulsion integral for H2
:[font = subsubsection; inactive; ]
Extension to:  "Gaussian Correlation Functions:  Two-Electron Systems"
		Lester & Krauss, Journal of Chemical Physics, Vol. 41, #5, 01Sep64
:[font = smalltext; inactive; ]
Conventions used here are taken from Lester & Krauss...
		
		General WaveFunction = Exp[-alpha1(r1-RA)**2 -alpha2(r2-RB)**2 -beta12(r1-r2)**2]
		
		alpha1(L)=ALPHA1,  alpha2(L)=ALPHA2,  beta12(L)=ALPHA5
		alpha1(K)=ALPHA3,  alpha2(K)=ALPHA4,  beta12(K)=ALPHA6

		GAMMA1=ALPHA1 + ALPHA2,  GAMMA2=ALPHA3 + ALPHA4
		GAMMA3=ALPHA5 + ALPHA6

		RAB=RAl - RBl,  RCD=RAk - RBk,  RA=RAl, RB=RBl, RC=RAk, RD=RBk

		Psi(L) = Exp [ -(  ALPHA1*(r1l-RA)**2 + ALPHA2*(r2l-RB)**2 + ALPHA5(r1l-r2l)**2  ) ]
		Psi(K) = Exp [ -(  ALPHA3*(r1lk-RC)**2 + ALPHA4*(r2k-RD)**2 + ALPHA6(r1k-r2k)**2  ) ]

	Lester & Krauss show that for molecular hydrogen, <Psi(K) | Psi(L)> is equal to...

		(Pi**3 * K)/(  (GAMMA1*GAMMA2 + GAMMA3(GAMMA1+GAMMA2))**(3/2)  )
	
	where K is equal to...

		Exp [ - ( ((ALPHA1*ALPHA2*RAB^2)/(GAMMA1))  +  ((ALPHA3*ALPHA4*RCD^2)/(GAMMA2))
                      + ((GAMMA1*GAMMA2*GAMMA3*PQ^2)/(GAMMA1*GAMMA2 + GAMMA3(GAMMA1+GAMMA2)))
               	           ) ]

	where PQ = P-Q, and...

		P=(ALPHA1*RA + ALPHA2*RB)/GAMMA1
		Q=(ALPHA3*RC + ALPHA4*RD)/GAMMA2

:[font = input; ]
Ar=RAl;
Br=RBl;
Cr=RAk;
Dr=RBk;
AB=Ar-Br;
CD=Cr-Dr;
a1=al1;
a2=al2;
a3=ak1;
a4=ak2;
a5=bl;
a6=bk;
gamma1=a1+a2;
gamma2=a3+a4;
gamma3=a5+a6;
Pr=(a1*Ar + a2*Br)/(a1+a2);
Qr=(a3*Cr + a4*Dr)/(a3+a4);
PQ=Pr-Qr;
K=Exp [ - (    	((a1*a2*AB^2)/(gamma1))
	        +    ((a3*a4*CD^2)/(gamma2))
                      +    ((gamma1*gamma2*gamma3*PQ^2)/(gamma1*gamma2 + gamma3(gamma1+gamma2)))
               ) ];
EREP=(2*Pi^(5/2) * K)/((gamma1*gamma2 + gamma3*(gamma1+gamma2))*(gamma1+gamma2)^(1/2))*Erf[
((gamma1^2 * gamma2^2 *PQ^2)/((gamma1+gamma2)*(gamma1*(gamma2+gamma3) +
	gamma2*gamma3)))];
:[font = input; startGroup; ]
FortranForm[ReplaceAll[D[EREP,al1], {(al1*RAl + al2*RBl)/(al1+al2) -> P,(ak1*RAk + ak2*RBk)/(ak1+ak2) -> Q,
		ak1+ak2+al1+al2->GAMMA1 + GAMMA2,ak1+ak2+bk+bl->GAMMA3 + GAMMA2,
		al1+al2->GAMMA1,ak1+ak2->GAMMA2,RAk-RBk->RCD,RAl-RBl->RAB,RAl->RA,RBl->RB,
		RAk->RC,RBk->RD,bk+bl->GAMMA3,al1->ALPHA1,al2->ALPHA2,ak1->ALPHA3,
		ak2->ALPHA4,bl->ALPHA5,bk->ALPHA6
		}]]
:[font = output; inactive; formatted; output; endGroup; nowordwrap; ]
4*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
       GAMMA1**4*GAMMA2**4*(P - Q)**4/
        ((GAMMA1 + GAMMA2)**2*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))**2) - 
       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**2*
    (-(GAMMA1**2*GAMMA2**2*(GAMMA2 + GAMMA3)*(P - Q)**2/
         ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))**2)) - 
      GAMMA1**2*GAMMA2**2*(P - Q)**2/
       ((GAMMA1 + GAMMA2)**2*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))) + 
      2*GAMMA1*GAMMA2**2*(P - Q)**2/
       ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))) + 
      2*GAMMA1**2*GAMMA2**2*(P - Q)*(RA/GAMMA1 - (ALPHA1*RA + ALPHA2*RB)/GAMMA1**2)/
       ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))))/
   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*(GAMMA2 + GAMMA3)*Pi**(5/2)*
    Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
      ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))))/
   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
  E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**(5/2)*
    Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
      ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))))/
   ((GAMMA1 + GAMMA2)**(3/2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**(5/2)*
    (GAMMA1*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)*(P - Q)**2/
       (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
      GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
      ALPHA1*ALPHA2*RAB**2/GAMMA1**2 - ALPHA2*RAB**2/GAMMA1 - 
      2*GAMMA1*GAMMA2*GAMMA3*(P - Q)*(RA/GAMMA1 - (ALPHA1*RA + ALPHA2*RB)/GAMMA1**2)/
       (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
    Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
      ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))))/
   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))


;[o]
4*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -       GAMMA1**4*GAMMA2**4*(P - Q)**4/
     -        ((GAMMA1 + GAMMA2)**2*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))**2) - 
     -       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**2*
     -    (-(GAMMA1**2*GAMMA2**2*(GAMMA2 + GAMMA3)*(P - Q)**2/
     -         ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))**2)) - 
     -      GAMMA1**2*GAMMA2**2*(P - Q)**2/
     -       ((GAMMA1 + GAMMA2)**2*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))) + 
     -      2*GAMMA1*GAMMA2**2*(P - Q)**2/
     -       ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))) + 
     -      2*GAMMA1**2*GAMMA2**2*(P - Q)*(RA/GAMMA1 - (ALPHA1*RA + ALPHA2*RB)/GAMMA1**2)/
     -       ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))))/
     -   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*(GAMMA2 + GAMMA3)*
     -    Pi**(5/2)*Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
     -      ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))))/
     -   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -  E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**(5/2)*
     -    Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
     -      ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))))/
     -   ((GAMMA1 + GAMMA2)**(3/2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**(5/2)*
     -    (GAMMA1*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)*(P - Q)**2/
     -       (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -      GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
     -      ALPHA1*ALPHA2*RAB**2/GAMMA1**2 - ALPHA2*RAB**2/GAMMA1 - 
     -      2*GAMMA1*GAMMA2*GAMMA3*(P - Q)*(RA/GAMMA1 - (ALPHA1*RA + ALPHA2*RB)/GAMMA1**2)/
     -       (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -    Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
     -      ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))))/
     -   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))
:[font = input; startGroup; ]
FortranForm[ReplaceAll[D[EREP,al2], {(al1*RAl + al2*RBl)/(al1+al2) -> P,(ak1*RAk + ak2*RBk)/(ak1+ak2) -> Q,
		ak1+ak2+al1+al2->GAMMA1 + GAMMA2,ak1+ak2+bk+bl->GAMMA3 + GAMMA2,
		al1+al2->GAMMA1,ak1+ak2->GAMMA2,RAk-RBk->RCD,RAl-RBl->RAB,RAl->RA,RBl->RB,
		RAk->RC,RBk->RD,bk+bl->GAMMA3,al1->ALPHA1,al2->ALPHA2,ak1->ALPHA3,
		ak2->ALPHA4,bl->ALPHA5,bk->ALPHA6
		}]]
:[font = output; inactive; formatted; output; endGroup; nowordwrap; ]
4*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
       GAMMA1**4*GAMMA2**4*(P - Q)**4/
        ((GAMMA1 + GAMMA2)**2*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))**2) - 
       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**2*
    (-(GAMMA1**2*GAMMA2**2*(GAMMA2 + GAMMA3)*(P - Q)**2/
         ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))**2)) - 
      GAMMA1**2*GAMMA2**2*(P - Q)**2/
       ((GAMMA1 + GAMMA2)**2*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))) + 
      2*GAMMA1*GAMMA2**2*(P - Q)**2/
       ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))) + 
      2*GAMMA1**2*GAMMA2**2*(P - Q)*(RB/GAMMA1 - (ALPHA1*RA + ALPHA2*RB)/GAMMA1**2)/
       ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))))/
   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*(GAMMA2 + GAMMA3)*Pi**(5/2)*
    Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
      ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))))/
   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
  E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**(5/2)*
    Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
      ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))))/
   ((GAMMA1 + GAMMA2)**(3/2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**(5/2)*
    (GAMMA1*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)*(P - Q)**2/
       (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
      GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
      ALPHA1*ALPHA2*RAB**2/GAMMA1**2 - ALPHA1*RAB**2/GAMMA1 - 
      2*GAMMA1*GAMMA2*GAMMA3*(P - Q)*(RB/GAMMA1 - (ALPHA1*RA + ALPHA2*RB)/GAMMA1**2)/
       (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
    Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
      ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))))/
   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))


;[o]
4*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -       GAMMA1**4*GAMMA2**4*(P - Q)**4/
     -        ((GAMMA1 + GAMMA2)**2*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))**2) - 
     -       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**2*
     -    (-(GAMMA1**2*GAMMA2**2*(GAMMA2 + GAMMA3)*(P - Q)**2/
     -         ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))**2)) - 
     -      GAMMA1**2*GAMMA2**2*(P - Q)**2/
     -       ((GAMMA1 + GAMMA2)**2*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))) + 
     -      2*GAMMA1*GAMMA2**2*(P - Q)**2/
     -       ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))) + 
     -      2*GAMMA1**2*GAMMA2**2*(P - Q)*(RB/GAMMA1 - (ALPHA1*RA + ALPHA2*RB)/GAMMA1**2)/
     -       ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))))/
     -   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*(GAMMA2 + GAMMA3)*
     -    Pi**(5/2)*Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
     -      ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))))/
     -   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -  E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**(5/2)*
     -    Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
     -      ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))))/
     -   ((GAMMA1 + GAMMA2)**(3/2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**(5/2)*
     -    (GAMMA1*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)*(P - Q)**2/
     -       (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -      GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
     -      ALPHA1*ALPHA2*RAB**2/GAMMA1**2 - ALPHA1*RAB**2/GAMMA1 - 
     -      2*GAMMA1*GAMMA2*GAMMA3*(P - Q)*(RB/GAMMA1 - (ALPHA1*RA + ALPHA2*RB)/GAMMA1**2)/
     -       (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -    Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
     -      ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))))/
     -   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))
:[font = input; startGroup; ]
FortranForm[ReplaceAll[D[EREP,bl], {(al1*RAl + al2*RBl)/(al1+al2) -> P,(ak1*RAk + ak2*RBk)/(ak1+ak2) -> Q,
		ak1+ak2+al1+al2->GAMMA1 + GAMMA2,ak1+ak2+bk+bl->GAMMA3 + GAMMA2,
		al1+al2->GAMMA1,ak1+ak2->GAMMA2,RAk-RBk->RCD,RAl-RBl->RAB,RAl->RA,RBl->RB,
		RAk->RC,RBk->RD,bk+bl->GAMMA3,al1->ALPHA1,al2->ALPHA2,ak1->ALPHA3,
		ak2->ALPHA4,bl->ALPHA5,bk->ALPHA6
		}]]
:[font = output; inactive; formatted; output; endGroup; nowordwrap; ]
-4*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
       GAMMA1**4*GAMMA2**4*(P - Q)**4/
        ((GAMMA1 + GAMMA2)**2*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))**2) - 
       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*GAMMA1**2*GAMMA2**2*Pi**2*
    (P - Q)**2/(Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)*
     (GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))**2) - 
  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Sqrt(GAMMA1 + GAMMA2)*
    Pi**(5/2)*Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
      ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))))/
   (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**(5/2)*
    (GAMMA1*GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3*(P - Q)**2/
       (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
      GAMMA1*GAMMA2*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
    Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
      ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))))/
   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))


;[o]
-4*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -       GAMMA1**4*GAMMA2**4*(P - Q)**4/
     -        ((GAMMA1 + GAMMA2)**2*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))**2) - 
     -       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*GAMMA1**2*GAMMA2**2*
     -    Pi**2*(P - Q)**2/
     -   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)*
     -     (GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))**2) - 
     -  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Sqrt(GAMMA1 + GAMMA2)*
     -    Pi**(5/2)*Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
     -      ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))))/
     -   (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**(5/2)*
     -    (GAMMA1*GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3*(P - Q)**2/
     -       (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -      GAMMA1*GAMMA2*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -    Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
     -      ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))))/
     -   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))
:[font = input; startGroup; ]
FortranForm[ReplaceAll[D[EREP,RAl], {(al1*RAl + al2*RBl)/(al1+al2) -> P,(ak1*RAk + ak2*RBk)/(ak1+ak2) -> Q,
		ak1+ak2+al1+al2->GAMMA1 + GAMMA2,ak1+ak2+bk+bl->GAMMA3 + GAMMA2,
		al1+al2->GAMMA1,ak1+ak2->GAMMA2,RAk-RBk->RCD,RAl-RBl->RAB,RAl->RA,RBl->RB,
		RAk->RC,RBk->RD,bk+bl->GAMMA3,al1->ALPHA1,al2->ALPHA2,ak1->ALPHA3,
		ak2->ALPHA4,bl->ALPHA5,bk->ALPHA6
		}]]
:[font = output; inactive; formatted; output; endGroup; nowordwrap; ]
8*ALPHA1*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
       GAMMA1**4*GAMMA2**4*(P - Q)**4/
        ((GAMMA1 + GAMMA2)**2*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))**2) - 
       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*GAMMA1*GAMMA2**2*Pi**2*
    (P - Q)/((GAMMA1 + GAMMA2)**(3/2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)*
     (GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))) + 
  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**(5/2)*
    (-2*ALPHA1*GAMMA2*GAMMA3*(P - Q)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
      2*ALPHA1*ALPHA2*RAB/GAMMA1)*Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
      ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))))/
   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))


;[o]
8*ALPHA1*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -       GAMMA1**4*GAMMA2**4*(P - Q)**4/
     -        ((GAMMA1 + GAMMA2)**2*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))**2) - 
     -       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*GAMMA1*GAMMA2**2*Pi**2*
     -    (P - Q)/((GAMMA1 + GAMMA2)**(3/2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)*
     -     (GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))) + 
     -  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**(5/2)*
     -    (-2*ALPHA1*GAMMA2*GAMMA3*(P - Q)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -      2*ALPHA1*ALPHA2*RAB/GAMMA1)*
     -    Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
     -      ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))))/
     -   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))
:[font = input; startGroup; ]
FortranForm[ReplaceAll[D[EREP,RBl], {(al1*RAl + al2*RBl)/(al1+al2) -> P,(ak1*RAk + ak2*RBk)/(ak1+ak2) -> Q,
		ak1+ak2+al1+al2->GAMMA1 + GAMMA2,ak1+ak2+bk+bl->GAMMA3 + GAMMA2,
		al1+al2->GAMMA1,ak1+ak2->GAMMA2,RAk-RBk->RCD,RAl-RBl->RAB,RAl->RA,RBl->RB,
		RAk->RC,RBk->RD,bk+bl->GAMMA3,al1->ALPHA1,al2->ALPHA2,ak1->ALPHA3,
		ak2->ALPHA4,bl->ALPHA5,bk->ALPHA6
		}]]
:[font = output; inactive; formatted; output; endGroup; nowordwrap; ]
8*ALPHA2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
       GAMMA1**4*GAMMA2**4*(P - Q)**4/
        ((GAMMA1 + GAMMA2)**2*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))**2) - 
       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*GAMMA1*GAMMA2**2*Pi**2*
    (P - Q)/((GAMMA1 + GAMMA2)**(3/2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)*
     (GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))) + 
  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**(5/2)*
    (-2*ALPHA2*GAMMA2*GAMMA3*(P - Q)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
      2*ALPHA1*ALPHA2*RAB/GAMMA1)*Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
      ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))))/
   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))


;[o]
8*ALPHA2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -       GAMMA1**4*GAMMA2**4*(P - Q)**4/
     -        ((GAMMA1 + GAMMA2)**2*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))**2) - 
     -       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*GAMMA1*GAMMA2**2*Pi**2*
     -    (P - Q)/((GAMMA1 + GAMMA2)**(3/2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)*
     -     (GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))) + 
     -  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**(5/2)*
     -    (-2*ALPHA2*GAMMA2*GAMMA3*(P - Q)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
     -      2*ALPHA1*ALPHA2*RAB/GAMMA1)*
     -    Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
     -      ((GAMMA1 + GAMMA2)*(GAMMA2*GAMMA3 + GAMMA1*(GAMMA2 + GAMMA3))))/
     -   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))
:[font = input; startGroup; backColorRed = 65280; backColorGreen = 65280; backColorBlue = 65280; fontColorRed = 0; fontColorGreen = 0; fontColorBlue = 32768; plain; fontName = "Arial"; fontSize = 10; ]
FortranForm[ReplaceAll[D[EREP,ak1], {(al1*RAl + al2*RBl)/(al1+al2) -> P,(ak1*RAk + ak2*RBk)/(ak1+ak2) -> Q,
		ak1+ak2+al1+al2->GAMMA1 + GAMMA2,al1+al2+bk+bl->GAMMA3 + GAMMA1,
		al1+al2->GAMMA1,ak1+ak2->GAMMA2,RAk-RBk->RCD,RAl-RBl->RAB,RAl->RA,RBl->RB,
		RAk->RC,RBk->RD,bk+bl->GAMMA3,al1->ALPHA1,al2->ALPHA2,ak1->ALPHA3,
		ak2->ALPHA4,bl->ALPHA5,bk->ALPHA6
		}]]
:[font = output; inactive; formatted; output; endGroup; nowordwrap; ]
4*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
       GAMMA1**4*GAMMA2**4*(P - Q)**4/
        ((GAMMA1 + GAMMA2)**2*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)**2) - 
       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**2*
    (-(GAMMA1**2*GAMMA2**2*(GAMMA1 + GAMMA3)*(P - Q)**2/
         ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)**2)) - 
      GAMMA1**2*GAMMA2**2*(P - Q)**2/
       ((GAMMA1 + GAMMA2)**2*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)) + 
      2*GAMMA1**2*GAMMA2*(P - Q)**2/
       ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)) + 
      2*GAMMA1**2*GAMMA2**2*(P - Q)*(-(RC/GAMMA2) + (ALPHA3*RC + ALPHA4*RD)/GAMMA2**2)/
       ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)))/
   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*(GAMMA1 + GAMMA3)*Pi**(5/2)*
    Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
      ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)))/
   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
  E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**(5/2)*
    Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
      ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)))/
   ((GAMMA1 + GAMMA2)**(3/2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**(5/2)*
    (GAMMA1*GAMMA2*GAMMA3*(GAMMA1 + GAMMA3)*(P - Q)**2/
       (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
      GAMMA1*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
      ALPHA3*ALPHA4*RCD**2/GAMMA2**2 - ALPHA4*RCD**2/GAMMA2 - 
      2*GAMMA1*GAMMA2*GAMMA3*(P - Q)*(-(RC/GAMMA2) + (ALPHA3*RC + ALPHA4*RD)/GAMMA2**2)/
       (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
    Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
      ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)))/
   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))


;[o]
4*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -       GAMMA1**4*GAMMA2**4*(P - Q)**4/
     -        ((GAMMA1 + GAMMA2)**2*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)**2) - 
     -       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**2*
     -    (-(GAMMA1**2*GAMMA2**2*(GAMMA1 + GAMMA3)*(P - Q)**2/
     -         ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)**2)) - 
     -      GAMMA1**2*GAMMA2**2*(P - Q)**2/
     -       ((GAMMA1 + GAMMA2)**2*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)) + 
     -      2*GAMMA1**2*GAMMA2*(P - Q)**2/
     -       ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)) + 
     -      2*GAMMA1**2*GAMMA2**2*(P - Q)*(-(RC/GAMMA2) + (ALPHA3*RC + ALPHA4*RD)/GAMMA2**2)/
     -       ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)))/
     -   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*(GAMMA1 + GAMMA3)*
     -    Pi**(5/2)*Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
     -      ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)))/
     -   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -  E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**(5/2)*
     -    Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
     -      ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)))/
     -   ((GAMMA1 + GAMMA2)**(3/2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**(5/2)*
     -    (GAMMA1*GAMMA2*GAMMA3*(GAMMA1 + GAMMA3)*(P - Q)**2/
     -       (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -      GAMMA1*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
     -      ALPHA3*ALPHA4*RCD**2/GAMMA2**2 - ALPHA4*RCD**2/GAMMA2 - 
     -      2*GAMMA1*GAMMA2*GAMMA3*(P - Q)*(-(RC/GAMMA2) + (ALPHA3*RC + ALPHA4*RD)/GAMMA2**2)/
     -       (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -    Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
     -      ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)))/
     -   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))
:[font = input; startGroup; backColorRed = 65280; backColorGreen = 65280; backColorBlue = 65280; fontColorRed = 0; fontColorGreen = 0; fontColorBlue = 32768; plain; fontName = "Arial"; fontSize = 10; ]
FortranForm[ReplaceAll[D[EREP,ak2], {(al1*RAl + al2*RBl)/(al1+al2) -> P,(ak1*RAk + ak2*RBk)/(ak1+ak2) -> Q,
		ak1+ak2+al1+al2->GAMMA1 + GAMMA2,al1+al2+bk+bl->GAMMA3 + GAMMA1,
		al1+al2->GAMMA1,ak1+ak2->GAMMA2,RAk-RBk->RCD,RAl-RBl->RAB,RAl->RA,RBl->RB,
		RAk->RC,RBk->RD,bk+bl->GAMMA3,al1->ALPHA1,al2->ALPHA2,ak1->ALPHA3,
		ak2->ALPHA4,bl->ALPHA5,bk->ALPHA6
		}]]
:[font = output; inactive; formatted; output; endGroup; nowordwrap; ]
4*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
       GAMMA1**4*GAMMA2**4*(P - Q)**4/
        ((GAMMA1 + GAMMA2)**2*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)**2) - 
       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**2*
    (-(GAMMA1**2*GAMMA2**2*(GAMMA1 + GAMMA3)*(P - Q)**2/
         ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)**2)) - 
      GAMMA1**2*GAMMA2**2*(P - Q)**2/
       ((GAMMA1 + GAMMA2)**2*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)) + 
      2*GAMMA1**2*GAMMA2*(P - Q)**2/
       ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)) + 
      2*GAMMA1**2*GAMMA2**2*(P - Q)*(-(RD/GAMMA2) + (ALPHA3*RC + ALPHA4*RD)/GAMMA2**2)/
       ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)))/
   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*(GAMMA1 + GAMMA3)*Pi**(5/2)*
    Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
      ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)))/
   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
  E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**(5/2)*
    Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
      ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)))/
   ((GAMMA1 + GAMMA2)**(3/2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**(5/2)*
    (GAMMA1*GAMMA2*GAMMA3*(GAMMA1 + GAMMA3)*(P - Q)**2/
       (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
      GAMMA1*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
      ALPHA3*ALPHA4*RCD**2/GAMMA2**2 - ALPHA3*RCD**2/GAMMA2 - 
      2*GAMMA1*GAMMA2*GAMMA3*(P - Q)*(-(RD/GAMMA2) + (ALPHA3*RC + ALPHA4*RD)/GAMMA2**2)/
       (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
    Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
      ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)))/
   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))


;[o]
4*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -       GAMMA1**4*GAMMA2**4*(P - Q)**4/
     -        ((GAMMA1 + GAMMA2)**2*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)**2) - 
     -       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**2*
     -    (-(GAMMA1**2*GAMMA2**2*(GAMMA1 + GAMMA3)*(P - Q)**2/
     -         ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)**2)) - 
     -      GAMMA1**2*GAMMA2**2*(P - Q)**2/
     -       ((GAMMA1 + GAMMA2)**2*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)) + 
     -      2*GAMMA1**2*GAMMA2*(P - Q)**2/
     -       ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)) + 
     -      2*GAMMA1**2*GAMMA2**2*(P - Q)*(-(RD/GAMMA2) + (ALPHA3*RC + ALPHA4*RD)/GAMMA2**2)/
     -       ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)))/
     -   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*(GAMMA1 + GAMMA3)*
     -    Pi**(5/2)*Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
     -      ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)))/
     -   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -  E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**(5/2)*
     -    Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
     -      ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)))/
     -   ((GAMMA1 + GAMMA2)**(3/2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**(5/2)*
     -    (GAMMA1*GAMMA2*GAMMA3*(GAMMA1 + GAMMA3)*(P - Q)**2/
     -       (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -      GAMMA1*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
     -      ALPHA3*ALPHA4*RCD**2/GAMMA2**2 - ALPHA3*RCD**2/GAMMA2 - 
     -      2*GAMMA1*GAMMA2*GAMMA3*(P - Q)*(-(RD/GAMMA2) + (ALPHA3*RC + ALPHA4*RD)/GAMMA2**2)/
     -       (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -    Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
     -      ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)))/
     -   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))
:[font = input; startGroup; backColorRed = 65280; backColorGreen = 65280; backColorBlue = 65280; fontColorRed = 0; fontColorGreen = 0; fontColorBlue = 32768; plain; fontName = "Arial"; fontSize = 10; ]
FortranForm[ReplaceAll[D[EREP,bk], {(al1*RAl + al2*RBl)/(al1+al2) -> P,(ak1*RAk + ak2*RBk)/(ak1+ak2) -> Q,
		ak1+ak2+al1+al2->GAMMA1 + GAMMA2,al1+al2+bk+bl->GAMMA3 + GAMMA1,
		al1+al2->GAMMA1,ak1+ak2->GAMMA2,RAk-RBk->RCD,RAl-RBl->RAB,RAl->RA,RBl->RB,
		RAk->RC,RBk->RD,bk+bl->GAMMA3,al1->ALPHA1,al2->ALPHA2,ak1->ALPHA3,
		ak2->ALPHA4,bl->ALPHA5,bk->ALPHA6
		}]]
:[font = output; inactive; formatted; output; endGroup; nowordwrap; ]
-4*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
       GAMMA1**4*GAMMA2**4*(P - Q)**4/
        ((GAMMA1 + GAMMA2)**2*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)**2) - 
       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*GAMMA1**2*GAMMA2**2*Pi**2*
    (P - Q)**2/(Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)**2*
     (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Sqrt(GAMMA1 + GAMMA2)*
    Pi**(5/2)*Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
      ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)))/
   (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**(5/2)*
    (GAMMA1*GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3*(P - Q)**2/
       (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
      GAMMA1*GAMMA2*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
    Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
      ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)))/
   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))


;[o]
-4*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -       GAMMA1**4*GAMMA2**4*(P - Q)**4/
     -        ((GAMMA1 + GAMMA2)**2*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)**2) - 
     -       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*GAMMA1**2*GAMMA2**2*
     -    Pi**2*(P - Q)**2/
     -   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)**2*
     -     (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Sqrt(GAMMA1 + GAMMA2)*
     -    Pi**(5/2)*Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
     -      ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)))/
     -   (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**(5/2)*
     -    (GAMMA1*GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3*(P - Q)**2/
     -       (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -      GAMMA1*GAMMA2*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -    Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
     -      ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)))/
     -   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))
:[font = input; startGroup; backColorRed = 65280; backColorGreen = 65280; backColorBlue = 65280; fontColorRed = 0; fontColorGreen = 0; fontColorBlue = 32768; plain; fontName = "Arial"; fontSize = 10; ]
FortranForm[ReplaceAll[D[EREP,RAk], {(al1*RAl + al2*RBl)/(al1+al2) -> P,(ak1*RAk + ak2*RBk)/(ak1+ak2) -> Q,
		ak1+ak2+al1+al2->GAMMA1 + GAMMA2,al1+al2+bk+bl->GAMMA3 + GAMMA1,
		al1+al2->GAMMA1,ak1+ak2->GAMMA2,RAk-RBk->RCD,RAl-RBl->RAB,RAl->RA,RBl->RB,
		RAk->RC,RBk->RD,bk+bl->GAMMA3,al1->ALPHA1,al2->ALPHA2,ak1->ALPHA3,
		ak2->ALPHA4,bl->ALPHA5,bk->ALPHA6
		}]]
:[font = output; inactive; formatted; output; endGroup; nowordwrap; ]
-8*ALPHA3*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
       GAMMA1**4*GAMMA2**4*(P - Q)**4/
        ((GAMMA1 + GAMMA2)**2*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)**2) - 
       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*GAMMA1**2*GAMMA2*Pi**2*
    (P - Q)/((GAMMA1 + GAMMA2)**(3/2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)*
     (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**(5/2)*
    (2*ALPHA3*GAMMA1*GAMMA3*(P - Q)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
      2*ALPHA3*ALPHA4*RCD/GAMMA2)*Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
      ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)))/
   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))


;[o]
-8*ALPHA3*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -       GAMMA1**4*GAMMA2**4*(P - Q)**4/
     -        ((GAMMA1 + GAMMA2)**2*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)**2) - 
     -       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*GAMMA1**2*GAMMA2*Pi**2*
     -    (P - Q)/((GAMMA1 + GAMMA2)**(3/2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)*
     -     (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**(5/2)*
     -    (2*ALPHA3*GAMMA1*GAMMA3*(P - Q)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -      2*ALPHA3*ALPHA4*RCD/GAMMA2)*
     -    Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
     -      ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)))/
     -   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))
:[font = input; startGroup; backColorRed = 65280; backColorGreen = 65280; backColorBlue = 65280; fontColorRed = 0; fontColorGreen = 0; fontColorBlue = 32768; plain; fontName = "Arial"; fontSize = 10; ]
FortranForm[ReplaceAll[D[EREP,RBk], {(al1*RAl + al2*RBl)/(al1+al2) -> P,(ak1*RAk + ak2*RBk)/(ak1+ak2) -> Q,
		ak1+ak2+al1+al2->GAMMA1 + GAMMA2,al1+al2+bk+bl->GAMMA3 + GAMMA1,
		al1+al2->GAMMA1,ak1+ak2->GAMMA2,RAk-RBk->RCD,RAl-RBl->RAB,RAl->RA,RBl->RB,
		RAk->RC,RBk->RD,bk+bl->GAMMA3,al1->ALPHA1,al2->ALPHA2,ak1->ALPHA3,
		ak2->ALPHA4,bl->ALPHA5,bk->ALPHA6
		}]]
:[font = output; inactive; formatted; output; endGroup; nowordwrap; ]
-8*ALPHA4*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
       GAMMA1**4*GAMMA2**4*(P - Q)**4/
        ((GAMMA1 + GAMMA2)**2*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)**2) - 
       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*GAMMA1**2*GAMMA2*Pi**2*
    (P - Q)/((GAMMA1 + GAMMA2)**(3/2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)*
     (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**(5/2)*
    (2*ALPHA4*GAMMA1*GAMMA3*(P - Q)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
      2*ALPHA3*ALPHA4*RCD/GAMMA2)*Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
      ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)))/
   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))


;[o]
-8*ALPHA4*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -       GAMMA1**4*GAMMA2**4*(P - Q)**4/
     -        ((GAMMA1 + GAMMA2)**2*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)**2) - 
     -       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*GAMMA1**2*GAMMA2*Pi**2*
     -    (P - Q)/((GAMMA1 + GAMMA2)**(3/2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)*
     -     (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -  2*E**(-(GAMMA1*GAMMA2*GAMMA3*(P - Q)**2/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -       ALPHA1*ALPHA2*RAB**2/GAMMA1 - ALPHA3*ALPHA4*RCD**2/GAMMA2)*Pi**(5/2)*
     -    (2*ALPHA4*GAMMA1*GAMMA3*(P - Q)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
     -      2*ALPHA3*ALPHA4*RCD/GAMMA2)*
     -    Erf(GAMMA1**2*GAMMA2**2*(P - Q)**2/
     -      ((GAMMA1 + GAMMA2)*(GAMMA1*(bk + bl + GAMMA2) + GAMMA2*GAMMA3)))/
     -   (Sqrt(GAMMA1 + GAMMA2)*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))
^*)