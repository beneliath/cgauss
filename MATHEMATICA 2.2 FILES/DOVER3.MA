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
:[font = subtitle; inactive; nohscroll; center; backColorRed = 65280; backColorGreen = 65280; backColorBlue = 65280; fontColorRed = 0; fontColorGreen = 0; fontColorBlue = 0; bold; fontName = "Arial"; fontSize = 18; ]
Example of Direct Derivative of Fortran Code... w/ Test.
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
ax=.;ay=.;az=.;bx=.;by=.;bz=.;cx=.;cy=.;cz=.;dx=.;dy=.;dz=.;
alk1=a1+a2;
alk2=a3+a4;                                                  
blk12=a5+a6;
detab=alk1*alk2+(alk1+alk2)*blk12;
ov=Pi^3 / (detab*Sqrt[detab]);

xk1=Exp[-((a1*a2)/alk1)*((ax-bx)^2+(ay-by)^2+(az-bz)^2)];
xk2=Exp[-((a3*a4)/alk2)*((cx-dx)^2+(cy-dy)^2+(cz-dz)^2)];
rx1=(a1*ax+a2*bx)/alk1;
rx2=(a3*cx+a4*dx)/alk2;
ry1=(a1*ay+a2*by)/alk1;
ry2=(a3*cy+a4*dy)/alk2;
rz1=(a1*az+a2*bz)/alk1;
rz2=(a3*cz+a4*dz)/alk2;
dd=(rx1-rx2)^2+(ry1-ry2)^2+(rz1-rz2)^2;
xk=xk1*xk2*Exp[-(alk1*alk2*blk12)*dd/detab];

ov=ov*xk;
:[font = input; startGroup; ]
Ga1=D[ov,a1];
ReplaceRepeated[%,{xk1*xk2*Exp[-(alk1*alk2*blk12)*dd/detab]->KFactor}];
ReplaceRepeated[%,{alk1*alk2+(alk1+alk2)*blk12 -> DETAB}];
ReplaceRepeated[%,{(rx1-rx2)^2+(ry1-ry2)^2+(rz1-rz2)^2 -> DD}];
ReplaceRepeated[%,{Exp[-((a3*a4)/alk2)*((cx-dx)^2+(cy-dy)^2+(cz-dz)^2)] -> XK2}];
ReplaceRepeated[%,{Exp[-((a1*a2)/alk1)*((ax-bx)^2+(ay-by)^2+(az-bz)^2)] -> XK1}];
ReplaceRepeated[%,{(a3*cz+a4*dz)/alk2 -> RZ2}];
ReplaceRepeated[%,{(a1*az+a2*bz)/alk1 -> RZ1}];
ReplaceRepeated[%,{(a3*cy+a4*dy)/alk2 -> RY2}];
ReplaceRepeated[%,{(a1*ay+a2*by)/alk1 -> RY1}];
ReplaceRepeated[%,{(a3*cx+a4*dx)/alk2 -> RX2}];
ReplaceRepeated[%,{(a1*ax+a2*bx)/alk1 -> RX1}];
ReplaceRepeated[%,{a1+a2 -> ALK1}];
ReplaceRepeated[%,{a3+a4 -> ALK2}];
ReplaceRepeated[%,{a5+a6 -> BLK12}];
ReplaceRepeated[%,{(KFactor Pi^3)/((DETAB)^(3/2)) -> OVERLAP}];

Ga1=FortranForm[%]
:[font = output; inactive; formatted; output; endGroup; nowordwrap; ]
-3*(ALK2 + BLK12)*KFactor*Pi**3/(2*DETAB**(5/2)) + 
  OVERLAP*(-(a2*((ax - bx)**2 + (ay - by)**2 + (az - bz)**2)/ALK1) + 
     a1*a2*((ax - bx)**2 + (ay - by)**2 + (az - bz)**2)/ALK1**2 + 
     ALK1*ALK2*BLK12*(ALK2 + BLK12)*DD/DETAB**2 - ALK2*BLK12*DD/DETAB - 
     ALK1*ALK2*BLK12*(2*(ax/ALK1 - (ax*a1 + a2*bx)/ALK1**2)*
          (RX1 - RX2) + 2*(ay/ALK1 - (ay*a1 + a2*by)/ALK1**2)*
          (RY1 - RY2) + 2*(az/ALK1 - (az*a1 + a2*bz)/ALK1**2)*(RZ1 - RZ2)
         )/DETAB)


;[o]
-3*(ALK2 + BLK12)*KFactor*Pi**3/(2*DETAB**(5/2)) + 
     -  OVERLAP*(-(a2*((ax - bx)**2 + (ay - by)**2 + (az - bz)**2)/
     -        ALK1) + a1*a2*((ax - bx)**2 + (ay - by)**2 + (az - bz)**2)/
     -      ALK1**2 + ALK1*ALK2*BLK12*(ALK2 + BLK12)*DD/DETAB**2 - 
     -     ALK2*BLK12*DD/DETAB - 
     -     ALK1*ALK2*BLK12*(2*(ax/ALK1 - (ax*a1 + a2*bx)/ALK1**2)*
     -          (RX1 - RX2) + 
     -         2*(ay/ALK1 - (ay*a1 + a2*by)/ALK1**2)*(RY1 - RY2) + 
     -         2*(az/ALK1 - (az*a1 + a2*bz)/ALK1**2)*(RZ1 - RZ2))/DETAB)
:[font = input; startGroup; ]
a1=.22;a2=.44;a3=.99;a4=.37;a5=.04;a6=.009;
ax=0;ay=0;az=0.7;
bx=0;by=0;bz=0.7;
cx=0;cy=0;cz=-0.7;
dx=0;dy=0;dz=-0.7;
POINT1=N[ov];pass1=a1;
a1=.2200001;
POINT2=N[ov];pass2=a1;
FDG=(POINT2 - POINT1)  /  (pass2 - pass1);
a1=pass1;
N[FDG,15]
:[font = output; inactive; formatted; output; endGroup; nowordwrap; ]
-60.86888401138184


;[o]
-60.8688840113818
^*)