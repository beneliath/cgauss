(*^

::[	frontEndVersion = "Microsoft Windows Mathematica Notebook Front End Version 2.2";
	microsoftWindowsStandardFontEncoding;
	fontset = title, "Arial", 24, L0, center, nohscroll, bold;
	fontset = subtitle, "Arial", 18, L0, center, nohscroll, bold;
	fontset = subsubtitle, "Arial", 14, L0, center, nohscroll, bold;
	fontset = section, "Arial", 14, L0, bold, grayBox;
	fontset = subsection, "Arial", 12, L0, bold, blackBox;
	fontset = subsubsection, "Arial", 10, L0, bold, whiteBox;
	fontset = text, "Times New Roman", 8, L0, italic;
	fontset = smalltext, "Arial", 10, L0;
	fontset = input, "Arial", 10, L0, B32768, grayDot;
	fontset = output, "Courier New", 9, L0, nohscroll, R32768;
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
:[font = title; inactive; nohscroll; backColorRed = 65280; backColorGreen = 65280; backColorBlue = 65280; fontColorRed = 0; fontColorGreen = 0; fontColorBlue = 0; bold; fontName = "Arial"; fontSize = 24; ]
TESTED AND APPROVED!
:[font = subtitle; inactive; nohscroll; center; backColorRed = 65280; backColorGreen = 65280; backColorBlue = 65280; fontColorRed = 0; fontColorGreen = 0; fontColorBlue = 0; bold; fontName = "Arial"; fontSize = 18; ]
Derivatives of cgauss overlap integral for H2
:[font = subsubsection; inactive; backColorRed = 65280; backColorGreen = 65280; backColorBlue = 65280; fontColorRed = 0; fontColorGreen = 0; fontColorBlue = 0; bold; fontName = "Arial"; fontSize = 10; ]
Extension to:  "Gaussian Correlation Functions:  Two-Electron Systems"
		Lester & Krauss, Journal of Chemical Physics, Vol. 41, #5, 01Sep64
:[font = text; inactive; backColorRed = 65280; backColorGreen = 65280; backColorBlue = 65280; fontColorRed = 0; fontColorGreen = 0; fontColorBlue = 0; italic; fontName = "Times New Roman"; fontSize = 8; ]
Conventions used here are taken from Lester & Krauss...
	The Overlap Integral here is:  <Psi(K) | Psi(L)> over the electron coordinates...
		
		General WaveFunction = Exp[-alpha1(r1-RA)**2 -alpha2(r2-RB)**2 -beta12(r1-r2)**2]
		
		alpha1(L)=ALPHA1,  alpha2(L)=ALPHA3,  beta12(L)=ALPHA5
		alpha1(K)=ALPHA2,  alpha2(K)=ALPHA4,  beta12(K)=ALPHA6

		GAMMA1=ALPHA1 + ALPHA2,  GAMMA2=ALPHA3 + ALPHA4
		GAMMA3=ALPHA5 + ALPHA6

		RA=RAl, RB=RBl, RC=RAk, RD=RBk

		Psi(L) = Exp [ -(  ALPHA1*(r1l-RA)**2 + ALPHA3*(r2l-RB)**2 + ALPHA5(r1l-r2l)**2  ) ]
		Psi(K) = Exp [ -(  ALPHA2*(r1lk-RC)**2 + ALPHA4*(r2k-RD)**2 + ALPHA6(r1k-r2k)**2  ) ]

	Lester & Krauss show that for molecular hydrogen, <Psi(K) | Psi(L)> is equal to...

		(Pi**3 * K)/(  (GAMMA1*GAMMA2 + GAMMA3(GAMMA1+GAMMA2))**(3/2)  )
	
	where K is equal to...

		Exp [ - ( ((ALPHA1*ALPHA2*RAB^2)/(GAMMA1))  +  ((ALPHA3*ALPHA4*RCD^2)/(GAMMA2))
                      		+ ((GAMMA1*GAMMA2*GAMMA3*PQ^2)/(GAMMA1*GAMMA2 + GAMMA3(GAMMA1+GAMMA2)))
               	          	              ) ]

	where RAB^2=(RAx-RBx)^2 + (RAy-RBy)^2 + (RAz-RBz)^2 ...

	and where PQ = P-Q, with...

		P=(ALPHA1*RA + ALPHA2*RB)/GAMMA1
		Q=(ALPHA3*RC + ALPHA4*RD)/GAMMA2

	where RA=RAx+RAy+RAz ... RA,RB,RC,RD Being Vectors.
:[font = input; ]
a1=.;a2=.;a3=.;a4=.;a5=.;a6=.;

gamma1=a1+a2;
gamma2=a3+a4;
gamma3=a5+a6;

RAx=.;RAy=.;RAz=.;
RBx=.;RBy=.;RBz=.;
RCx=.;RCy=.;RCz=.;
RDx=.;RDy=.;RDz=.;

Px=(a1*RAx + a2*RBx)/(gamma1);
Py=(a1*RAy + a2*RBy)/(gamma1);
Pz=(a1*RAz + a2*RBz)/(gamma1);
Qx=(a3*RCx + a4*RDx)/(gamma2);
Qy=(a3*RCy + a4*RDy)/(gamma2);
Qz=(a3*RCz + a4*RDz)/(gamma2);

detAB=(gamma1*gamma2 + gamma3(gamma1+gamma2));

ABSQUARED=((RAx-RBx)^2 + (RAy-RBy)^2 + (RAz-RBz)^2);
CDSQUARED=((RCx-RDx)^2 + (RCy-RDy)^2 + (RCz-RDz)^2);
PQSQUARED= ((Px-Qx)^2 + (Py-Qy)^2  + (Pz-Qz)^2);

Theta= - (   (a1*a2)/(gamma1)  ABSQUARED
           	+ (a3*a4)/(gamma2)  CDSQUARED
           	+ ( (gamma1*gamma2*gamma3*   PQSQUARED)  / detAB   ));

K=Exp [Theta];

OV=(Pi^3 * K)/(detAB)^(3/2);
:[font = input; startGroup; ]
D[OV,a1];
ReplaceRepeated[%,{(Pi^3 * K)/((gamma1*gamma2 + gamma3(gamma1+gamma2))^(3/2))->OVERLAP}];
ReplaceRepeated[%,{Theta->THETA}];
ReplaceRepeated[%,{Exp[THETA]->KFactor}];
ReplaceRepeated[%,{(gamma1*gamma2 + gamma3(gamma1+gamma2))->DETAB}];
ReplaceRepeated[%,{(Pi^3 * KFactor)->(OVERLAP * DETAB^(3/2))}];
ReplaceRepeated[%,{(a3*RCz + a4*RDz)/(gamma2)->QZ}];
ReplaceRepeated[%,{(a3*RCy + a4*RDy)/(gamma2)->QY}];
ReplaceRepeated[%,{(a3*RCx + a4*RDx)/(gamma2)->QX}];
ReplaceRepeated[%,{(a1*RAz + a2*RBz)/(gamma1)->PZ}];
ReplaceRepeated[%,{(a1*RAy + a2*RBy)/(gamma1)->PY}];
ReplaceRepeated[%,{(a1*RAx + a2*RBx)/(gamma1)->PX}];
ReplaceRepeated[%,{((Px-Qx)^2 + (Py-Qy)^2  + (Pz-Qz)^2)->PQsquared}];
ReplaceRepeated[%,{(PX-QX)->PQx,(PY-QY)->PQy,(PZ-QZ)->PQx}];
ReplaceRepeated[%,{((RAx-RBx)^2 + (RAy-RBy)^2 + (RAz-RBz)^2)->ABsquared}];
ReplaceRepeated[%,{((RCx-RDx)^2 + (RCy-RDy)^2 + (RCz-RDz)^2)->CDsquared}];
ReplaceRepeated[%,{a1+a2->GAMMA1}];
ReplaceRepeated[%,{a3+a4->GAMMA2}];
ReplaceRepeated[%,{a5+a6->GAMMA3}];
ReplaceRepeated[%,{ (Pi^3 * KFactor) / (DETAB^(3/2)) ->OVERLAP}];
ReplaceRepeated[%,{RAx->0,RBx->0,RCx->0,RDx->0,
	RAy->0,RBy->0,RCy->0,RDy->0}];
FortranForm[%]
:[font = output; inactive; formatted; output; endGroup; nohscroll; ]
-3*(GAMMA2 + GAMMA3)*OVERLAP/(2*DETAB) + 
  OVERLAP*(ABsquared*a1*a2/GAMMA1**2 - ABsquared*a2/GAMMA1 - 
     GAMMA2*GAMMA3*(2*PQx**2 + PQy**2)/DETAB + 
     GAMMA1*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)*(2*PQx**2 + PQy**2)/DETAB**2 - 
     2*GAMMA1*GAMMA2*GAMMA3*PQx*(RAz/GAMMA1 - (a1*RAz + a2*RBz)/GAMMA1**2)/DETAB)


;[o]
-3*(GAMMA2 + GAMMA3)*OVERLAP/(2*DETAB) + 
     -  OVERLAP*(ABsquared*a1*a2/GAMMA1**2 - ABsquared*a2/GAMMA1 - 
     -     GAMMA2*GAMMA3*(2*PQx**2 + PQy**2)/DETAB + 
     -     GAMMA1*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)*(2*PQx**2 + PQy**2)/DETAB**2 - 
     -     2*GAMMA1*GAMMA2*GAMMA3*PQx*(RAz/GAMMA1 - (a1*RAz + a2*RBz)/GAMMA1**2)/DETAB)
:[font = input; startGroup; backColorRed = 65280; backColorGreen = 65280; backColorBlue = 65280; fontColorRed = 0; fontColorGreen = 0; fontColorBlue = 32768; plain; fontName = "Arial"; fontSize = 10; ]
D[OV,a2];
ReplaceRepeated[%,{(Pi^3 * K)/((gamma1*gamma2 + gamma3(gamma1+gamma2))^(3/2))->OVERLAP}];
ReplaceRepeated[%,{Theta->THETA}];
ReplaceRepeated[%,{Exp[THETA]->KFactor}];
ReplaceRepeated[%,{(gamma1*gamma2 + gamma3(gamma1+gamma2))->DETAB}];
ReplaceRepeated[%,{(Pi^3 * KFactor)->(OVERLAP * DETAB^(3/2))}];
ReplaceRepeated[%,{(a3*RCz + a4*RDz)/(gamma2)->QZ}];
ReplaceRepeated[%,{(a3*RCy + a4*RDy)/(gamma2)->QY}];
ReplaceRepeated[%,{(a3*RCx + a4*RDx)/(gamma2)->QX}];
ReplaceRepeated[%,{(a1*RAz + a2*RBz)/(gamma1)->PZ}];
ReplaceRepeated[%,{(a1*RAy + a2*RBy)/(gamma1)->PY}];
ReplaceRepeated[%,{(a1*RAx + a2*RBx)/(gamma1)->PX}];
ReplaceRepeated[%,{((Px-Qx)^2 + (Py-Qy)^2  + (Pz-Qz)^2)->PQsquared}];
ReplaceRepeated[%,{(PX-QX)->PQx,(PY-QY)->PQy,(PZ-QZ)->PQx}];
ReplaceRepeated[%,{((RAx-RBx)^2 + (RAy-RBy)^2 + (RAz-RBz)^2)->ABsquared}];
ReplaceRepeated[%,{((RCx-RDx)^2 + (RCy-RDy)^2 + (RCz-RDz)^2)->CDsquared}];
ReplaceRepeated[%,{a1+a2->GAMMA1}];
ReplaceRepeated[%,{a3+a4->GAMMA2}];
ReplaceRepeated[%,{a5+a6->GAMMA3}];
ReplaceRepeated[%,{ (Pi^3 * KFactor) / (DETAB^(3/2)) ->OVERLAP}];
ReplaceRepeated[%,{RAx->0,RBx->0,RCx->0,RDx->0,
	RAy->0,RBy->0,RCy->0,RDy->0}];
FortranForm[%]
:[font = output; inactive; formatted; output; endGroup; nohscroll; ]
-3*(GAMMA2 + GAMMA3)*OVERLAP/(2*DETAB) + 
  OVERLAP*(ABsquared*a1*a2/GAMMA1**2 - ABsquared*a1/GAMMA1 - 
     GAMMA2*GAMMA3*(2*PQx**2 + PQy**2)/DETAB + 
     GAMMA1*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)*(2*PQx**2 + PQy**2)/DETAB**2 - 
     2*GAMMA1*GAMMA2*GAMMA3*PQx*(RBz/GAMMA1 - (a1*RAz + a2*RBz)/GAMMA1**2)/DETAB)


;[o]
-3*(GAMMA2 + GAMMA3)*OVERLAP/(2*DETAB) + 
     -  OVERLAP*(ABsquared*a1*a2/GAMMA1**2 - ABsquared*a1/GAMMA1 - 
     -     GAMMA2*GAMMA3*(2*PQx**2 + PQy**2)/DETAB + 
     -     GAMMA1*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)*(2*PQx**2 + PQy**2)/DETAB**2 - 
     -     2*GAMMA1*GAMMA2*GAMMA3*PQx*(RBz/GAMMA1 - (a1*RAz + a2*RBz)/GAMMA1**2)/DETAB)
:[font = input; startGroup; backColorRed = 65280; backColorGreen = 65280; backColorBlue = 65280; fontColorRed = 0; fontColorGreen = 0; fontColorBlue = 32768; plain; fontName = "Arial"; fontSize = 10; ]
D[OV,a5];
ReplaceRepeated[%,{(Pi^3 * K)/((gamma1*gamma2 + gamma3(gamma1+gamma2))^(3/2))->OVERLAP}];
ReplaceRepeated[%,{Theta->THETA}];
ReplaceRepeated[%,{Exp[THETA]->KFactor}];
ReplaceRepeated[%,{(gamma1*gamma2 + gamma3(gamma1+gamma2))->DETAB}];
ReplaceRepeated[%,{(Pi^3 * KFactor)->(OVERLAP * DETAB^(3/2))}];
ReplaceRepeated[%,{(a3*RCz + a4*RDz)/(gamma2)->QZ}];
ReplaceRepeated[%,{(a3*RCy + a4*RDy)/(gamma2)->QY}];
ReplaceRepeated[%,{(a3*RCx + a4*RDx)/(gamma2)->QX}];
ReplaceRepeated[%,{(a1*RAz + a2*RBz)/(gamma1)->PZ}];
ReplaceRepeated[%,{(a1*RAy + a2*RBy)/(gamma1)->PY}];
ReplaceRepeated[%,{(a1*RAx + a2*RBx)/(gamma1)->PX}];
ReplaceRepeated[%,{((Px-Qx)^2 + (Py-Qy)^2  + (Pz-Qz)^2)->PQsquared}];
ReplaceRepeated[%,{(PX-QX)->PQx,(PY-QY)->PQy,(PZ-QZ)->PQx}];
ReplaceRepeated[%,{((RAx-RBx)^2 + (RAy-RBy)^2 + (RAz-RBz)^2)->ABsquared}];
ReplaceRepeated[%,{((RCx-RDx)^2 + (RCy-RDy)^2 + (RCz-RDz)^2)->CDsquared}];
ReplaceRepeated[%,{a1+a2->GAMMA1}];
ReplaceRepeated[%,{a3+a4->GAMMA2}];
ReplaceRepeated[%,{a5+a6->GAMMA3}];
ReplaceRepeated[%,{ (Pi^3 * KFactor) / (DETAB^(3/2)) ->OVERLAP}];
ReplaceRepeated[%,{RAx->0,RBx->0,RCx->0,RDx->0,
	RAy->0,RBy->0,RCy->0,RDy->0}];
FortranForm[%]
:[font = output; inactive; formatted; output; endGroup; nohscroll; ]
-3*(GAMMA1 + GAMMA2)*OVERLAP/(2*DETAB) + 
  OVERLAP*(-(GAMMA1*GAMMA2*(2*PQx**2 + PQy**2)/DETAB) + 
     GAMMA1*GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3*(2*PQx**2 + PQy**2)/DETAB**2)


;[o]
-3*(GAMMA1 + GAMMA2)*OVERLAP/(2*DETAB) + 
     -  OVERLAP*(-(GAMMA1*GAMMA2*(2*PQx**2 + PQy**2)/DETAB) + 
     -     GAMMA1*GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3*(2*PQx**2 + PQy**2)/DETAB**2)
:[font = input; startGroup; backColorRed = 65280; backColorGreen = 65280; backColorBlue = 65280; fontColorRed = 0; fontColorGreen = 0; fontColorBlue = 32768; plain; fontName = "Arial"; fontSize = 10; ]
D[OV,RAz];
ReplaceRepeated[%,{(Pi^3 * K)/((gamma1*gamma2 + gamma3(gamma1+gamma2))^(3/2))->OVERLAP}];
ReplaceRepeated[%,{Theta->THETA}];
ReplaceRepeated[%,{Exp[THETA]->KFactor}];
ReplaceRepeated[%,{(gamma1*gamma2 + gamma3(gamma1+gamma2))->DETAB}];
ReplaceRepeated[%,{(Pi^3 * KFactor)->(OVERLAP * DETAB^(3/2))}];
ReplaceRepeated[%,{(a3*RCz + a4*RDz)/(gamma2)->QZ}];
ReplaceRepeated[%,{(a3*RCy + a4*RDy)/(gamma2)->QY}];
ReplaceRepeated[%,{(a3*RCx + a4*RDx)/(gamma2)->QX}];
ReplaceRepeated[%,{(a1*RAz + a2*RBz)/(gamma1)->PZ}];
ReplaceRepeated[%,{(a1*RAy + a2*RBy)/(gamma1)->PY}];
ReplaceRepeated[%,{(a1*RAx + a2*RBx)/(gamma1)->PX}];
ReplaceRepeated[%,{((Px-Qx)^2 + (Py-Qy)^2  + (Pz-Qz)^2)->PQsquared}];
ReplaceRepeated[%,{(PX-QX)->PQx,(PY-QY)->PQy,(PZ-QZ)->PQx}];
ReplaceRepeated[%,{((RAx-RBx)^2 + (RAy-RBy)^2 + (RAz-RBz)^2)->ABsquared}];
ReplaceRepeated[%,{((RCx-RDx)^2 + (RCy-RDy)^2 + (RCz-RDz)^2)->CDsquared}];
ReplaceRepeated[%,{a1+a2->GAMMA1}];
ReplaceRepeated[%,{a3+a4->GAMMA2}];
ReplaceRepeated[%,{a5+a6->GAMMA3}];
ReplaceRepeated[%,{ (Pi^3 * KFactor) / (DETAB^(3/2)) ->OVERLAP}];
ReplaceRepeated[%,{RAx->0,RBx->0,RCx->0,RDx->0,
	RAy->0,RBy->0,RCy->0,RDy->0}];
FortranForm[%]
:[font = output; inactive; formatted; output; endGroup; nohscroll; ]
OVERLAP*(-2*a1*GAMMA2*GAMMA3*PQx/DETAB - 2*a1*a2*(RAz - RBz)/GAMMA1)


;[o]
OVERLAP*(-2*a1*GAMMA2*GAMMA3*PQx/DETAB - 2*a1*a2*(RAz - RBz)/GAMMA1)
:[font = input; startGroup; backColorRed = 65280; backColorGreen = 65280; backColorBlue = 65280; fontColorRed = 0; fontColorGreen = 0; fontColorBlue = 32768; plain; fontName = "Arial"; fontSize = 10; ]
D[OV,RBz];
ReplaceRepeated[%,{(Pi^3 * K)/((gamma1*gamma2 + gamma3(gamma1+gamma2))^(3/2))->OVERLAP}];
ReplaceRepeated[%,{Theta->THETA}];
ReplaceRepeated[%,{Exp[THETA]->KFactor}];
ReplaceRepeated[%,{(gamma1*gamma2 + gamma3(gamma1+gamma2))->DETAB}];
ReplaceRepeated[%,{(Pi^3 * KFactor)->(OVERLAP * DETAB^(3/2))}];
ReplaceRepeated[%,{(a3*RCz + a4*RDz)/(gamma2)->QZ}];
ReplaceRepeated[%,{(a3*RCy + a4*RDy)/(gamma2)->QY}];
ReplaceRepeated[%,{(a3*RCx + a4*RDx)/(gamma2)->QX}];
ReplaceRepeated[%,{(a1*RAz + a2*RBz)/(gamma1)->PZ}];
ReplaceRepeated[%,{(a1*RAy + a2*RBy)/(gamma1)->PY}];
ReplaceRepeated[%,{(a1*RAx + a2*RBx)/(gamma1)->PX}];
ReplaceRepeated[%,{((Px-Qx)^2 + (Py-Qy)^2  + (Pz-Qz)^2)->PQsquared}];
ReplaceRepeated[%,{(PX-QX)->PQx,(PY-QY)->PQy,(PZ-QZ)->PQx}];
ReplaceRepeated[%,{((RAx-RBx)^2 + (RAy-RBy)^2 + (RAz-RBz)^2)->ABsquared}];
ReplaceRepeated[%,{((RCx-RDx)^2 + (RCy-RDy)^2 + (RCz-RDz)^2)->CDsquared}];
ReplaceRepeated[%,{a1+a2->GAMMA1}];
ReplaceRepeated[%,{a3+a4->GAMMA2}];
ReplaceRepeated[%,{a5+a6->GAMMA3}];
ReplaceRepeated[%,{ (Pi^3 * KFactor) / (DETAB^(3/2)) ->OVERLAP}];
ReplaceRepeated[%,{RAx->0,RBx->0,RCx->0,RDx->0,
	RAy->0,RBy->0,RCy->0,RDy->0}];
FortranForm[%]
:[font = output; inactive; formatted; output; endGroup; nohscroll; ]
OVERLAP*(-2*a2*GAMMA2*GAMMA3*PQx/DETAB + 2*a1*a2*(RAz - RBz)/GAMMA1)


;[o]
OVERLAP*(-2*a2*GAMMA2*GAMMA3*PQx/DETAB + 2*a1*a2*(RAz - RBz)/GAMMA1)
^*)