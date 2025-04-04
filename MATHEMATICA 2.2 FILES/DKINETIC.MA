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
:[font = input; ]
      GAMA1=al1+al2;
      GAMA2=ak1+ak2;                                                      
      GAMA3=bl+bk;
      PX=(al1*AX+al2*BX)/GAMA1;                                           
      QX=(ak1*CX+ak2*DX)/GAMA2;                                            
      PY=(al1*AY+al2*BY)/GAMA1;                                            
      QY=(ak1*CY+ak2*DY)/GAMA2;                                            
      PZ=(al1*AZ+al2*BZ)/GAMA1;                                            
      QZ=(ak1*CZ+ak2*DZ)/GAMA2;                                          
      QCX=QX-CX;
      QCY=QY-CY;
      QCZ=QZ-CZ;                              
      PAX=PX-AX;
      PAY=PY-AY;
      PAZ=PZ-AZ;                              
      PBX=PX-BX;
      PBY=PY-BY;
      PBZ=PZ-BZ;                              
      PQX=PX-QX;
      PQY=PY-QY;
      PQZ=PZ-QZ;                             
      PQX2=PQX*PQX;
      PQY2=PQY*PQY;
      PQZ2=PQZ*PQZ;                   
      G1PG2=GAMA1+GAMA2;
      G1RG2=GAMA1*GAMA2;
      XMIAN=G1RG2+GAMA3*G1PG2;     
      DELT=G1RG2*GAMA3/XMIAN*4;       
      GA1DEL=DELT/GAMA1;
      OVXA11=PAX-PQX*GA1DEL/4;        
      OVYA11=PAY-PQY*GA1DEL/4;       
      OVZA11=PAZ-PQZ*GA1DEL/4;
      OVXB11=PBX-PQX*GA1DEL/4;
      OVYB11=PBY-PQY*GA1DEL/4;               
      OVZB11=PBZ-PQZ*GA1DEL/4;                 
      OVXA12=OVXA11*OVXA11;            
      OVYA12=OVYA11*OVYA11;            
      OVZA12=OVZA11*OVZA11;          
      POM=0.5/GAMA1-GA1DEL/8/GAMA1;                       
      OVXA12=(OVXA12+POM)*OVLAP;               
      OVYA12=(OVYA12+POM)*OVLAP;             
      OVZA12=(OVZA12+POM)*OVLAP;            
      OVXB12=OVXB11*OVXB11;       
      OVYB12=OVYB11*OVYB11;      
      OVZB12=OVZB11*OVZB11;    
      OVXB12=(OVXB12+POM)*OVLAP;          
      OVYB12=(OVYB12+POM)*OVLAP;         
      OVZB12=(OVZB12+POM)*OVLAP;       
      GA2DEL=DELT/GAMA2;
      OVXC21=QCX+PQX*GA2DEL/4;          
      OVYC21=QCY+PQY*GA2DEL/4;          
      OVZC21=QCZ+PQZ*GA2DEL/4;        
      OVXC22=OVXC21*OVXC21*OVLAP;        
      OVYC22=OVYC21*OVYC21*OVLAP;       
      OVZC22=OVZC21*OVZC21*OVLAP;     
      POM=0.5/GAMA2-GA2DEL/8/GAMA2;           
      POM=POM*OVLAP;
      OVXC22=OVXC22+POM;   
      OVYC22=OVYC22+POM;  
      OVZC22=OVZC22+POM;
      POM=GA2DEL/8/GAMA1;   
      XA1XC2=1-0.5*PQX2*DELT;     
      YA1YC2=1-0.5*PQY2*DELT;    
      ZA1ZC2=1-0.5*PQZ2*DELT;  
      XA1XC2=POM*XA1XC2;
      YA1YC2=POM*YA1YC2;
      ZA1ZC2=POM*ZA1ZC2;
      XA1XC2=PAX*QCX+PQX/4*(PAX*GA2DEL-QCX*GA1DEL)+XA1XC2;
      YA1YC2=PAY*QCY+PQY/4*(PAY*GA2DEL-QCY*GA1DEL)+YA1YC2;
      ZA1ZC2=PAZ*QCz+PQZ/4*(PAZ*GA2DEL-QCZ*GA1DEL)+ZA1ZC2;
      XA1XC2=XA1XC2*OVLAP;
      YA1YC2=YA1YC2*OVLAP;
      ZA1ZC2=ZA1ZC2*OVLAP;
      ACX=AX-CX;
      ACY=AY-CY;
      ACZ=AZ-CZ;                              
      OVA12=OVXA12+OVYA12+OVZA12;
      OVA11=2*(OVXA11*ACX+OVYA11*ACY+OVZA11*ACZ)*OVLAP;
      OVC22=OVXC22+OVYC22+OVZC22;
      OVC21=-2*(OVXC21*ACX+OVYC21*ACY+OVZC21*ACZ)*OVLAP;
      A1C2=-2*(XA1XC2+YA1YC2+ZA1ZC2);
      TKIN=4*bl*bk*(OVA12+OVA11+OVC22+OVC21+A1C2);            
      A6DG3=bk/GAMA3;
      A5DG3=bl/GAMA3;
      TKIN=TKIN-4*A1*A1*A6DG3*OVA12;              
      OVB12=OVXB12+OVYB12+OVZB12;            
      TKIN=TKIN-4*A2*A2*A5DG3*OVB12;            
      TKIN=TKIN+(6*(A1*A6DG3+A2*A5DG3)+
           4*bl*bk*(ACX*ACX+ACY*ACY+ACZ*ACZ))*OVLAP;
:[font = input; startGroup; ]
FortranForm[ReplaceAll[D[TKIN,al1], {(al1*RAl + al2*RBl)/(al1+al2) -> P,(ak1*RAk + ak2*RBk)/(ak1+ak2) -> Q,
		ak1+ak2+al1+al2->GAMMA1 + GAMMA2,ak1+ak2+bk+bl->GAMMA3 + GAMMA2,
		al1+al2->GAMMA1,ak1+ak2->GAMMA2,RAk-RBk->RCD,RAl-RBl->RAB,RAl->RA,RBl->RB,
		RAk->RC,RBk->RD,bk+bl->GAMMA3,al1->ALPHA1,al2->ALPHA2,ak1->ALPHA3,
		ak2->ALPHA4,bl->ALPHA5,bk->ALPHA6
		}]]
:[font = output; inactive; formatted; output; endGroup; nowordwrap; ]
-4*ALPHA6*A1**2*((-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
         GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*(-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + AX/GAMMA1 + 
            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
              GAMMA3*(GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            (-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + AX/GAMMA1)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
      (-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
         GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*(-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + AY/GAMMA1 + 
            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
              GAMMA3*(GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            (-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + AY/GAMMA1)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
      (-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
         GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*(-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + AZ/GAMMA1 + 
            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
              GAMMA3*(GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + AZ/GAMMA1)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP)/GAMMA3 - 
  4*ALPHA5*A2**2*((-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
         GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*(-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + AX/GAMMA1 + 
            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
              GAMMA3*(GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            (-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + AX/GAMMA1)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-BX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
      (-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
         GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*(-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + AY/GAMMA1 + 
            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
              GAMMA3*(GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            (-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + AY/GAMMA1)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-BY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
      (-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
         GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*(-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + AZ/GAMMA1 + 
            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
              GAMMA3*(GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + AZ/GAMMA1)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-BZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP)/GAMMA3 + 
  4*ALPHA5*ALPHA6*(2*(-(GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
              (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3*(GAMMA2 + GAMMA3)/
           (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
        (-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + AX/GAMMA1)*GAMMA1*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
        ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
      (-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2 + 
        GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     2*(-(GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3*
            (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
        (-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + AY/GAMMA1)*GAMMA1*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
        ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
      (-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2 + 
        GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     2*(-(GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3*
            (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
        (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + AZ/GAMMA1)*GAMMA1*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
        ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
      (-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2 + 
        GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     3*(GAMMA1*GAMMA3*(GAMMA2 + GAMMA3)/
         (2*GAMMA2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
        GAMMA3/(2*GAMMA2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP - 
     2*((AX - CX)*(-(GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
                 (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3*(GAMMA2 + GAMMA3)/
              (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
           (-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + AX/GAMMA1)*GAMMA1*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AY - CY)*(-(GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
                 (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3*(GAMMA2 + GAMMA3)/
              (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
           (-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + AY/GAMMA1)*GAMMA1*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AZ - CZ)*(-(GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
                 (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3*(GAMMA2 + GAMMA3)/
              (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
           (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + AZ/GAMMA1)*GAMMA1*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     2*((AX - CX)*(-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + AX/GAMMA1 + 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3*
             (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           (-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + AX/GAMMA1)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AY - CY)*(-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + AY/GAMMA1 + 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3*
             (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           (-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + AY/GAMMA1)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AZ - CZ)*(-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + AZ/GAMMA1 + 
           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3*
             (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + AZ/GAMMA1)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     (-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
        GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        2*(-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + AX/GAMMA1 + 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3*
             (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           (-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + AX/GAMMA1)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
         (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     (-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
        GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        2*(-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + AY/GAMMA1 + 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3*
             (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           (-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + AY/GAMMA1)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
         (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     (-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
        GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        2*(-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + AZ/GAMMA1 + 
           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3*
             (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + AZ/GAMMA1)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
         (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP - 
     2*(((-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + AX/GAMMA1)*
            (-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2) + 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*
             (-4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA1*GAMMA3*(GAMMA2 + GAMMA3)/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
               4*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
               4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
               4*(-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + AX/GAMMA1)*GAMMA1*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
           GAMMA3*(2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
                    (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
               4.*(-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + AX/GAMMA1)*GAMMA1*
                 ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
                 GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               2.*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*
                 GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
           GAMMA3*(GAMMA2 + GAMMA3)*(1 - 
               2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
                    (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
           (-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + AX/GAMMA1)*
             (4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA1*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4)*OVLAP + 
        ((-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + AY/GAMMA1)*
            (-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2) + 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*
             (-4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA1*GAMMA3*(GAMMA2 + GAMMA3)/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
               4*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
               4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
               4*(-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + AY/GAMMA1)*GAMMA1*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
           GAMMA3*(2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
                    (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
               4.*(-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + AY/GAMMA1)*GAMMA1*
                 ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
                 GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               2.*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*
                 GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
           GAMMA3*(GAMMA2 + GAMMA3)*(1 - 
               2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
                    (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
           (-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + AY/GAMMA1)*
             (4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA1*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4)*OVLAP + 
        OVLAP*(((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*
             (-4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA1*GAMMA3*(GAMMA2 + GAMMA3)/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
               4*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
               4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
               4*(-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + AZ/GAMMA1)*GAMMA1*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
           GAMMA3*(2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
                    (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
               4.*(-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + AZ/GAMMA1)*GAMMA1*
                 ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
                 GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               2.*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA2*
                 GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
           GAMMA3*(GAMMA2 + GAMMA3)*(1 - 
               2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*
                 GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
           (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + AZ/GAMMA1)*
             (4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA1*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
           (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + AZ/GAMMA1)*QCz)))


;[o]
-4*ALPHA6*A1**2*((-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -         GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*(-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + AX/GAMMA1 + 
     -            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -              GAMMA3*(GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            (-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + AX/GAMMA1)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -      (-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -         GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*(-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + AY/GAMMA1 + 
     -            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -              GAMMA3*(GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            (-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + AY/GAMMA1)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -      (-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -         GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*(-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + AZ/GAMMA1 + 
     -            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -              GAMMA3*(GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + AZ/GAMMA1)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP)/GAMMA3 - 
     -  4*ALPHA5*A2**2*((-0.5/GAMMA1**2 + 
     -         GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -         GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*(-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + AX/GAMMA1 + 
     -            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -              GAMMA3*(GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            (-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + AX/GAMMA1)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-BX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -      (-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -         GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*(-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + AY/GAMMA1 + 
     -            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -              GAMMA3*(GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            (-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + AY/GAMMA1)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-BY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -      (-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -         GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*(-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + AZ/GAMMA1 + 
     -            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -              GAMMA3*(GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + AZ/GAMMA1)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-BZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP)/GAMMA3 + 
     -  4*ALPHA5*ALPHA6*(2*(-(GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -              (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3*(GAMMA2 + GAMMA3)/
     -           (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -        (-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + AX/GAMMA1)*GAMMA1*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
     -        ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -      (-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2 + 
     -        GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     -     2*(-(GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3*
     -            (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -        (-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + AY/GAMMA1)*GAMMA1*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
     -        ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -      (-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2 + 
     -        GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     -     2*(-(GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3*
     -            (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -        (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + AZ/GAMMA1)*GAMMA1*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
     -        ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -      (-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2 + 
     -        GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     -     3*(GAMMA1*GAMMA3*(GAMMA2 + GAMMA3)/
     -         (2*GAMMA2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -        GAMMA3/(2*GAMMA2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP - 
     -     2*((AX - CX)*(-(GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -                 (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3*(GAMMA2 + GAMMA3)/
     -              (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -           (-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + AX/GAMMA1)*GAMMA1*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
     -           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AY - CY)*(-(GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -                 (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3*(GAMMA2 + GAMMA3)/
     -              (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -           (-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + AY/GAMMA1)*GAMMA1*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
     -           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AZ - CZ)*(-(GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -                 (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3*(GAMMA2 + GAMMA3)/
     -              (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -           (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + AZ/GAMMA1)*GAMMA1*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
     -           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -     2*((AX - CX)*(-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + AX/GAMMA1 + 
     -           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3*
     -             (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -           (-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + AX/GAMMA1)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AY - CY)*(-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + AY/GAMMA1 + 
     -           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3*
     -             (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -           (-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + AY/GAMMA1)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AZ - CZ)*(-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + AZ/GAMMA1 + 
     -           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3*
     -             (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -           (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + AZ/GAMMA1)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -     (-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
     -         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -        GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        2*(-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + AX/GAMMA1 + 
     -           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3*
     -             (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -           (-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + AX/GAMMA1)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -         (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -     (-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
     -         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -        GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        2*(-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + AY/GAMMA1 + 
     -           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3*
     -             (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -           (-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + AY/GAMMA1)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -         (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -     (-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
     -         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -        GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        2*(-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + AZ/GAMMA1 + 
     -           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3*
     -             (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -           (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + AZ/GAMMA1)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -         (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP - 
     -     2*(((-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + AX/GAMMA1)*
     -            (-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2) + 
     -           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*
     -             (-4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA1*GAMMA3*(GAMMA2 + GAMMA3)/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -               4*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -               4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
     -               4*(-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + AX/GAMMA1)*GAMMA1*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
     -           GAMMA3*(2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -                    (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -               4.*(-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + AX/GAMMA1)*GAMMA1*
     -                 ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -                 GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               2.*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*
     -                 GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -           GAMMA3*(GAMMA2 + GAMMA3)*
     -             (1 - 2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -                    (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -           (-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + AX/GAMMA1)*
     -             (4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA1*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4)*OVLAP + 
     -        ((-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + AY/GAMMA1)*
     -            (-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2) + 
     -           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*
     -             (-4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA1*GAMMA3*(GAMMA2 + GAMMA3)/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -               4*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -               4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
     -               4*(-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + AY/GAMMA1)*GAMMA1*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
     -           GAMMA3*(2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -                    (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -               4.*(-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + AY/GAMMA1)*GAMMA1*
     -                 ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -                 GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               2.*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*
     -                 GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -           GAMMA3*(GAMMA2 + GAMMA3)*
     -             (1 - 2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -                    (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -           (-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + AY/GAMMA1)*
     -             (4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA1*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4)*OVLAP + 
     -        OVLAP*(((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*
     -             (-4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA1*GAMMA3*(GAMMA2 + GAMMA3)/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -               4*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -               4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
     -               4*(-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + AZ/GAMMA1)*GAMMA1*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
     -           GAMMA3*(2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -                    (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -               4.*(-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + AZ/GAMMA1)*GAMMA1*
     -                 ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -                 GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               2.*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*
     -                 GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -           GAMMA3*(GAMMA2 + GAMMA3)*
     -             (1 - 2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -                    (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -           (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + AZ/GAMMA1)*
     -             (4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA1*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
     -           (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + AZ/GAMMA1)*QCz)))
:[font = input; startGroup; ]
FortranForm[ReplaceAll[D[TKIN,al2], {(al1*RAl + al2*RBl)/(al1+al2) -> P,(ak1*RAk + ak2*RBk)/(ak1+ak2) -> Q,
		ak1+ak2+al1+al2->GAMMA1 + GAMMA2,ak1+ak2+bk+bl->GAMMA3 + GAMMA2,
		al1+al2->GAMMA1,ak1+ak2->GAMMA2,RAk-RBk->RCD,RAl-RBl->RAB,RAl->RA,RBl->RB,
		RAk->RC,RBk->RD,bk+bl->GAMMA3,al1->ALPHA1,al2->ALPHA2,ak1->ALPHA3,
		ak2->ALPHA4,bl->ALPHA5,bk->ALPHA6
		}]]
:[font = output; inactive; formatted; output; endGroup; nowordwrap; ]
-4*ALPHA6*A1**2*((-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
         GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*(-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + BX/GAMMA1 + 
            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
              GAMMA3*(GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            (-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + BX/GAMMA1)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
      (-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
         GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*(-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + BY/GAMMA1 + 
            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
              GAMMA3*(GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            (-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + BY/GAMMA1)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
      (-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
         GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*(-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + BZ/GAMMA1 + 
            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
              GAMMA3*(GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + BZ/GAMMA1)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP)/GAMMA3 - 
  4*ALPHA5*A2**2*((-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
         GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*(-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + BX/GAMMA1 + 
            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
              GAMMA3*(GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            (-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + BX/GAMMA1)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-BX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
      (-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
         GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*(-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + BY/GAMMA1 + 
            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
              GAMMA3*(GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            (-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + BY/GAMMA1)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-BY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
      (-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
         GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*(-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + BZ/GAMMA1 + 
            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
              GAMMA3*(GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + BZ/GAMMA1)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-BZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP)/GAMMA3 + 
  4*ALPHA5*ALPHA6*(2*(-(GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
              (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3*(GAMMA2 + GAMMA3)/
           (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
        (-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + BX/GAMMA1)*GAMMA1*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
        ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
      (-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2 + 
        GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     2*(-(GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3*
            (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
        (-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + BY/GAMMA1)*GAMMA1*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
        ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
      (-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2 + 
        GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     2*(-(GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3*
            (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
        (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + BZ/GAMMA1)*GAMMA1*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
        ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
      (-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2 + 
        GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     3*(GAMMA1*GAMMA3*(GAMMA2 + GAMMA3)/
         (2*GAMMA2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
        GAMMA3/(2*GAMMA2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP - 
     2*((AX - CX)*(-(GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
                 (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3*(GAMMA2 + GAMMA3)/
              (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
           (-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + BX/GAMMA1)*GAMMA1*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AY - CY)*(-(GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
                 (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3*(GAMMA2 + GAMMA3)/
              (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
           (-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + BY/GAMMA1)*GAMMA1*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AZ - CZ)*(-(GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
                 (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3*(GAMMA2 + GAMMA3)/
              (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
           (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + BZ/GAMMA1)*GAMMA1*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     2*((AX - CX)*(-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + BX/GAMMA1 + 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3*
             (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           (-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + BX/GAMMA1)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AY - CY)*(-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + BY/GAMMA1 + 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3*
             (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           (-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + BY/GAMMA1)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AZ - CZ)*(-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + BZ/GAMMA1 + 
           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3*
             (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + BZ/GAMMA1)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     (-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
        GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        2*(-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + BX/GAMMA1 + 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3*
             (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           (-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + BX/GAMMA1)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
         (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     (-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
        GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        2*(-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + BY/GAMMA1 + 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3*
             (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           (-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + BY/GAMMA1)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
         (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     (-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
        GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        2*(-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + BZ/GAMMA1 + 
           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3*
             (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + BZ/GAMMA1)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
         (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP - 
     2*(((-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + BX/GAMMA1)*
            (-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2) + 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*
             (-4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA1*GAMMA3*(GAMMA2 + GAMMA3)/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
               4*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
               4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
               4*(-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + BX/GAMMA1)*GAMMA1*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
           GAMMA3*(2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
                    (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
               4.*(-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + BX/GAMMA1)*GAMMA1*
                 ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
                 GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               2.*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*
                 GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
           GAMMA3*(GAMMA2 + GAMMA3)*(1 - 
               2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
                    (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
           (-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + BX/GAMMA1)*
             (4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA1*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4)*OVLAP + 
        ((-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + BY/GAMMA1)*
            (-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2) + 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*
             (-4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA1*GAMMA3*(GAMMA2 + GAMMA3)/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
               4*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
               4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
               4*(-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + BY/GAMMA1)*GAMMA1*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
           GAMMA3*(2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
                    (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
               4.*(-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + BY/GAMMA1)*GAMMA1*
                 ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
                 GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               2.*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*
                 GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
           GAMMA3*(GAMMA2 + GAMMA3)*(1 - 
               2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
                    (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
           (-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + BY/GAMMA1)*
             (4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA1*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4)*OVLAP + 
        OVLAP*(((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*
             (-4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA1*GAMMA3*(GAMMA2 + GAMMA3)/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
               4*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
               4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
               4*(-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + BZ/GAMMA1)*GAMMA1*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
           GAMMA3*(2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
                    (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
               4.*(-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + BZ/GAMMA1)*GAMMA1*
                 ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
                 GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               2.*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA2*
                 GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
           GAMMA3*(GAMMA2 + GAMMA3)*(1 - 
               2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*
                 GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
           (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + BZ/GAMMA1)*
             (4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA1*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
           (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + BZ/GAMMA1)*QCz)))


;[o]
-4*ALPHA6*A1**2*((-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -         GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*(-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + BX/GAMMA1 + 
     -            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -              GAMMA3*(GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            (-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + BX/GAMMA1)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -      (-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -         GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*(-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + BY/GAMMA1 + 
     -            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -              GAMMA3*(GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            (-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + BY/GAMMA1)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -      (-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -         GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*(-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + BZ/GAMMA1 + 
     -            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -              GAMMA3*(GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + BZ/GAMMA1)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP)/GAMMA3 - 
     -  4*ALPHA5*A2**2*((-0.5/GAMMA1**2 + 
     -         GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -         GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*(-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + BX/GAMMA1 + 
     -            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -              GAMMA3*(GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            (-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + BX/GAMMA1)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-BX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -      (-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -         GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*(-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + BY/GAMMA1 + 
     -            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -              GAMMA3*(GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            (-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + BY/GAMMA1)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-BY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -      (-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -         GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*(-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + BZ/GAMMA1 + 
     -            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -              GAMMA3*(GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + BZ/GAMMA1)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-BZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP)/GAMMA3 + 
     -  4*ALPHA5*ALPHA6*(2*(-(GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -              (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3*(GAMMA2 + GAMMA3)/
     -           (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -        (-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + BX/GAMMA1)*GAMMA1*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
     -        ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -      (-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2 + 
     -        GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     -     2*(-(GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3*
     -            (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -        (-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + BY/GAMMA1)*GAMMA1*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
     -        ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -      (-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2 + 
     -        GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     -     2*(-(GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3*
     -            (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -        (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + BZ/GAMMA1)*GAMMA1*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
     -        ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -      (-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2 + 
     -        GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     -     3*(GAMMA1*GAMMA3*(GAMMA2 + GAMMA3)/
     -         (2*GAMMA2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -        GAMMA3/(2*GAMMA2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP - 
     -     2*((AX - CX)*(-(GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -                 (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3*(GAMMA2 + GAMMA3)/
     -              (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -           (-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + BX/GAMMA1)*GAMMA1*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
     -           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AY - CY)*(-(GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -                 (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3*(GAMMA2 + GAMMA3)/
     -              (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -           (-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + BY/GAMMA1)*GAMMA1*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
     -           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AZ - CZ)*(-(GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -                 (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3*(GAMMA2 + GAMMA3)/
     -              (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -           (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + BZ/GAMMA1)*GAMMA1*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
     -           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -     2*((AX - CX)*(-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + BX/GAMMA1 + 
     -           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3*
     -             (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -           (-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + BX/GAMMA1)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AY - CY)*(-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + BY/GAMMA1 + 
     -           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3*
     -             (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -           (-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + BY/GAMMA1)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AZ - CZ)*(-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + BZ/GAMMA1 + 
     -           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3*
     -             (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -           (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + BZ/GAMMA1)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -     (-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
     -         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -        GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        2*(-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + BX/GAMMA1 + 
     -           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3*
     -             (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -           (-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + BX/GAMMA1)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -         (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -     (-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
     -         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -        GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        2*(-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + BY/GAMMA1 + 
     -           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3*
     -             (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -           (-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + BY/GAMMA1)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -         (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -     (-0.5/GAMMA1**2 + GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
     -         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -        GAMMA2*GAMMA3/(2*GAMMA1**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        2*(-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + BZ/GAMMA1 + 
     -           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3*
     -             (GAMMA2 + GAMMA3)/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -           (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + BZ/GAMMA1)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -         (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP - 
     -     2*(((-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + BX/GAMMA1)*
     -            (-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2) + 
     -           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*
     -             (-4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA1*GAMMA3*(GAMMA2 + GAMMA3)/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -               4*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -               4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
     -               4*(-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + BX/GAMMA1)*GAMMA1*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
     -           GAMMA3*(2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -                    (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -               4.*(-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + BX/GAMMA1)*GAMMA1*
     -                 ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -                 GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               2.*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*
     -                 GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -           GAMMA3*(GAMMA2 + GAMMA3)*
     -             (1 - 2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -                    (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -           (-((ALPHA1*AX + ALPHA2*BX)/GAMMA1**2) + BX/GAMMA1)*
     -             (4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA1*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4)*OVLAP + 
     -        ((-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + BY/GAMMA1)*
     -            (-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2) + 
     -           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*
     -             (-4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA1*GAMMA3*(GAMMA2 + GAMMA3)/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -               4*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -               4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
     -               4*(-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + BY/GAMMA1)*GAMMA1*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
     -           GAMMA3*(2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -                    (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -               4.*(-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + BY/GAMMA1)*GAMMA1*
     -                 ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -                 GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               2.*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*
     -                 GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -           GAMMA3*(GAMMA2 + GAMMA3)*
     -             (1 - 2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -                    (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -           (-((ALPHA1*AY + ALPHA2*BY)/GAMMA1**2) + BY/GAMMA1)*
     -             (4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA1*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4)*OVLAP + 
     -        OVLAP*(((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*
     -             (-4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA1*GAMMA3*(GAMMA2 + GAMMA3)/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -               4*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -               4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) + 
     -               4*(-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + BZ/GAMMA1)*GAMMA1*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
     -           GAMMA3*(2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -                    (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA2*GAMMA3*(GAMMA2 + GAMMA3)/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -               4.*(-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + BZ/GAMMA1)*GAMMA1*
     -                 ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -                 GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               2.*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*
     -                 GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -           GAMMA3*(GAMMA2 + GAMMA3)*
     -             (1 - 2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -                    (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -           (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + BZ/GAMMA1)*
     -             (4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA1*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
     -           (-((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1**2) + BZ/GAMMA1)*QCz)))
:[font = input; startGroup; backColorRed = 65280; backColorGreen = 65280; backColorBlue = 65280; fontColorRed = 0; fontColorGreen = 0; fontColorBlue = 32768; plain; fontName = "Arial"; fontSize = 10; ]
FortranForm[ReplaceAll[D[TKIN,ak1], {(al1*RAl + al2*RBl)/(al1+al2) -> P,(ak1*RAk + ak2*RBk)/(ak1+ak2) -> Q,
		ak1+ak2+al1+al2->GAMMA1 + GAMMA2,ak1+ak2+bk+bl->GAMMA3 + GAMMA2,
		al1+al2->GAMMA1,ak1+ak2->GAMMA2,RAk-RBk->RCD,RAl-RBl->RAB,RAl->RA,RBl->RB,
		RAk->RC,RBk->RD,bk+bl->GAMMA3,al1->ALPHA1,al2->ALPHA2,ak1->ALPHA3,
		ak2->ALPHA4,bl->ALPHA5,bk->ALPHA6
		}]]
:[font = output; inactive; formatted; output; endGroup; nowordwrap; ]
-4*ALPHA6*A1**2*(((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
         GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*((bk + bl + GAMMA1)*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
                (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
            ((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - CX/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
      ((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
         GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*((bk + bl + GAMMA1)*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
                (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
            ((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - CY/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
      ((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
         GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*((bk + bl + GAMMA1)*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
                (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
            ((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - CZ/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP)/GAMMA3 - 
  4*ALPHA5*A2**2*(((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
         GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*((bk + bl + GAMMA1)*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
                (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
            ((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - CX/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-BX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
      ((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
         GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*((bk + bl + GAMMA1)*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
                (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
            ((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - CY/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-BY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
      ((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
         GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*((bk + bl + GAMMA1)*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
                (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
            ((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - CZ/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-BZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP)/GAMMA3 + 
  4*ALPHA5*ALPHA6*(2*(-((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2) + CX/GAMMA2 - 
        GAMMA1*(bk + bl + GAMMA1)*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
            (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
        GAMMA1*((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - CX/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
      (-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2 + 
        GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     2*(-((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2) + CY/GAMMA2 - 
        GAMMA1*(bk + bl + GAMMA1)*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
            (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
        GAMMA1*((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - CY/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
      (-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2 + 
        GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     2*(-((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2) + CZ/GAMMA2 - 
        GAMMA1*(bk + bl + GAMMA1)*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
            (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
        GAMMA1*((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - CZ/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
      (-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2 + 
        GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     3*(-0.5/GAMMA2**2 + GAMMA1*(bk + bl + GAMMA1)*GAMMA3/
         (2*GAMMA2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
        GAMMA1*GAMMA3/(2*GAMMA2**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP - 
     2*((AX - CX)*(-((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2) + CX/GAMMA2 - 
           GAMMA1*(bk + bl + GAMMA1)*
             ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
           GAMMA1*((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - CX/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AY - CY)*(-((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2) + CY/GAMMA2 - 
           GAMMA1*(bk + bl + GAMMA1)*
             ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
           GAMMA1*((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - CY/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AZ - CZ)*(-((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2) + CZ/GAMMA2 - 
           GAMMA1*(bk + bl + GAMMA1)*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
               (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
           GAMMA1*((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - CZ/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     2*((AX - CX)*((bk + bl + GAMMA1)*
             ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
           ((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - CX/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AY - CY)*((bk + bl + GAMMA1)*
             ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
           ((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - CY/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AZ - CZ)*((bk + bl + GAMMA1)*
             ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
           ((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - CZ/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     ((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
        GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        2*((bk + bl + GAMMA1)*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
               (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
           ((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - CX/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
         (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     ((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
        GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        2*((bk + bl + GAMMA1)*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
               (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
           ((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - CY/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
         (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     ((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
        GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        2*((bk + bl + GAMMA1)*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
               (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
           ((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - CZ/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
         (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP - 
     2*(((-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*
            (-((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2) + CX/GAMMA2) + 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*
             (-4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA1*(bk + bl + GAMMA1)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
               4*(bk + bl + GAMMA1)*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
               4*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4*(-((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2) + CX/GAMMA2)*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
           GAMMA3*(2.*GAMMA1*(bk + bl + GAMMA1)*
                 ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*
                 GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
               2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
                    (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4.*GAMMA1*((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - CX/GAMMA2)*
                 ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
                 GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
           (bk + bl + GAMMA1)*GAMMA3*(1 - 
               2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
                    (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
           ((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - CX/GAMMA2)*
             (4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA1*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4)*OVLAP + 
        ((-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*
            (-((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2) + CY/GAMMA2) + 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*
             (-4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA1*(bk + bl + GAMMA1)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
               4*(bk + bl + GAMMA1)*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
               4*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4*(-((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2) + CY/GAMMA2)*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
           GAMMA3*(2.*GAMMA1*(bk + bl + GAMMA1)*
                 ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*
                 GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
               2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
                    (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4.*GAMMA1*((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - CY/GAMMA2)*
                 ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
                 GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
           (bk + bl + GAMMA1)*GAMMA3*(1 - 
               2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
                    (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
           ((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - CY/GAMMA2)*
             (4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA1*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4)*OVLAP + 
        (((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*
             (-4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA1*(bk + bl + GAMMA1)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
               4*(bk + bl + GAMMA1)*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
               4*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4*(-((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2) + CZ/GAMMA2)*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
           GAMMA3*(2.*GAMMA1*(bk + bl + GAMMA1)*
                 ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*
                 GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
               2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
                    (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4.*GAMMA1*((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - CZ/GAMMA2)*
                 ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
                 GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
           (bk + bl + GAMMA1)*GAMMA3*(1 - 
               2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**
                  2*GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
           ((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - CZ/GAMMA2)*
             (4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA1*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4)*OVLAP))


;[o]
-4*ALPHA6*A1**2*(((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -         GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*((bk + bl + GAMMA1)*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -                (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -            ((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - CX/GAMMA2)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -      ((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -         GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*((bk + bl + GAMMA1)*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -                (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -            ((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - CY/GAMMA2)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -      ((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -         GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*((bk + bl + GAMMA1)*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -                (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -            ((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - CZ/GAMMA2)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP)/GAMMA3 - 
     -  4*ALPHA5*A2**2*(((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -         GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*((bk + bl + GAMMA1)*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -                (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -            ((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - CX/GAMMA2)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-BX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -      ((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -         GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*((bk + bl + GAMMA1)*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -                (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -            ((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - CY/GAMMA2)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-BY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -      ((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -         GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*((bk + bl + GAMMA1)*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -                (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -            ((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - CZ/GAMMA2)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-BZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP)/GAMMA3 + 
     -  4*ALPHA5*ALPHA6*(2*(-((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2) + CX/GAMMA2 - 
     -        GAMMA1*(bk + bl + GAMMA1)*
     -          ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -        GAMMA1*((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - CX/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -      (-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2 + 
     -        GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     -     2*(-((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2) + CY/GAMMA2 - 
     -        GAMMA1*(bk + bl + GAMMA1)*
     -          ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -        GAMMA1*((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - CY/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -      (-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2 + 
     -        GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     -     2*(-((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2) + CZ/GAMMA2 - 
     -        GAMMA1*(bk + bl + GAMMA1)*
     -          ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -        GAMMA1*((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - CZ/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -      (-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2 + 
     -        GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     -     3*(-0.5/GAMMA2**2 + GAMMA1*(bk + bl + GAMMA1)*GAMMA3/
     -         (2*GAMMA2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -        GAMMA1*GAMMA3/(2*GAMMA2**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP - 
     -     2*((AX - CX)*(-((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2) + CX/GAMMA2 - 
     -           GAMMA1*(bk + bl + GAMMA1)*
     -             ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -           GAMMA1*((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - CX/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AY - CY)*(-((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2) + CY/GAMMA2 - 
     -           GAMMA1*(bk + bl + GAMMA1)*
     -             ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -           GAMMA1*((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - CY/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AZ - CZ)*(-((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2) + CZ/GAMMA2 - 
     -           GAMMA1*(bk + bl + GAMMA1)*
     -             ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -           GAMMA1*((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - CZ/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -     2*((AX - CX)*((bk + bl + GAMMA1)*
     -             ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -             GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -           ((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - CX/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AY - CY)*((bk + bl + GAMMA1)*
     -             ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -             GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -           ((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - CY/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AZ - CZ)*((bk + bl + GAMMA1)*
     -             ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -             GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -           ((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - CZ/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -     ((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
     -         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -        GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        2*((bk + bl + GAMMA1)*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -               (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -           ((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - CX/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -         (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -     ((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
     -         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -        GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        2*((bk + bl + GAMMA1)*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -               (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -           ((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - CY/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -         (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -     ((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
     -         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -        GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        2*((bk + bl + GAMMA1)*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -               (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -           ((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - CZ/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -         (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP - 
     -     2*(((-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*
     -            (-((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2) + CX/GAMMA2) + 
     -           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*
     -             (-4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA1*(bk + bl + GAMMA1)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -               4*(bk + bl + GAMMA1)*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -               4*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4*(-((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2) + CX/GAMMA2)*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
     -           GAMMA3*(2.*GAMMA1*(bk + bl + GAMMA1)*
     -                 ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*
     -                 GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -               2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -                    (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4.*GAMMA1*((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - CX/GAMMA2)*
     -                 ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -                 GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -           (bk + bl + GAMMA1)*GAMMA3*
     -             (1 - 2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -                    (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -           ((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - CX/GAMMA2)*
     -             (4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA1*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4)*OVLAP + 
     -        ((-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*
     -            (-((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2) + CY/GAMMA2) + 
     -           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*
     -             (-4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA1*(bk + bl + GAMMA1)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -               4*(bk + bl + GAMMA1)*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -               4*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4*(-((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2) + CY/GAMMA2)*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
     -           GAMMA3*(2.*GAMMA1*(bk + bl + GAMMA1)*
     -                 ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*
     -                 GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -               2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -                    (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4.*GAMMA1*((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - CY/GAMMA2)*
     -                 ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -                 GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -           (bk + bl + GAMMA1)*GAMMA3*
     -             (1 - 2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -                    (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -           ((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - CY/GAMMA2)*
     -             (4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA1*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4)*OVLAP + 
     -        (((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*
     -             (-4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA1*(bk + bl + GAMMA1)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -               4*(bk + bl + GAMMA1)*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -               4*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4*(-((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2) + CZ/GAMMA2)*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
     -           GAMMA3*(2.*GAMMA1*(bk + bl + GAMMA1)*
     -                 ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*
     -                 GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -               2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -                    (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4.*GAMMA1*((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - CZ/GAMMA2)*
     -                 ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -                 GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -           (bk + bl + GAMMA1)*GAMMA3*
     -             (1 - 2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -                    (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -           ((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - CZ/GAMMA2)*
     -             (4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA1*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4)*OVLAP))
:[font = input; startGroup; backColorRed = 65280; backColorGreen = 65280; backColorBlue = 65280; fontColorRed = 0; fontColorGreen = 0; fontColorBlue = 32768; plain; fontName = "Arial"; fontSize = 10; ]
FortranForm[ReplaceAll[D[TKIN,ak2], {(al1*RAl + al2*RBl)/(al1+al2) -> P,(ak1*RAk + ak2*RBk)/(ak1+ak2) -> Q,
		ak1+ak2+al1+al2->GAMMA1 + GAMMA2,ak1+ak2+bk+bl->GAMMA3 + GAMMA2,
		al1+al2->GAMMA1,ak1+ak2->GAMMA2,RAk-RBk->RCD,RAl-RBl->RAB,RAl->RA,RBl->RB,
		RAk->RC,RBk->RD,bk+bl->GAMMA3,al1->ALPHA1,al2->ALPHA2,ak1->ALPHA3,
		ak2->ALPHA4,bl->ALPHA5,bk->ALPHA6
		}]]
:[font = output; inactive; formatted; output; endGroup; nowordwrap; ]
-4*ALPHA6*A1**2*(((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
         GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*((bk + bl + GAMMA1)*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
                (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
            ((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - DX/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
      ((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
         GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*((bk + bl + GAMMA1)*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
                (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
            ((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - DY/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
      ((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
         GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*((bk + bl + GAMMA1)*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
                (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
            ((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - DZ/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP)/GAMMA3 - 
  4*ALPHA5*A2**2*(((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
         GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*((bk + bl + GAMMA1)*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
                (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
            ((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - DX/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-BX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
      ((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
         GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*((bk + bl + GAMMA1)*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
                (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
            ((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - DY/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-BY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
      ((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
         GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*((bk + bl + GAMMA1)*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
                (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
            ((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - DZ/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-BZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP)/GAMMA3 + 
  4*ALPHA5*ALPHA6*(2*(-((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2) + DX/GAMMA2 - 
        GAMMA1*(bk + bl + GAMMA1)*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
            (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
        GAMMA1*((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - DX/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
      (-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2 + 
        GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     2*(-((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2) + DY/GAMMA2 - 
        GAMMA1*(bk + bl + GAMMA1)*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
            (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
        GAMMA1*((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - DY/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
      (-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2 + 
        GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     2*(-((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2) + DZ/GAMMA2 - 
        GAMMA1*(bk + bl + GAMMA1)*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
            (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
        GAMMA1*((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - DZ/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
      (-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2 + 
        GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     3*(-0.5/GAMMA2**2 + GAMMA1*(bk + bl + GAMMA1)*GAMMA3/
         (2*GAMMA2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
        GAMMA1*GAMMA3/(2*GAMMA2**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP - 
     2*((AX - CX)*(-((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2) + DX/GAMMA2 - 
           GAMMA1*(bk + bl + GAMMA1)*
             ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
           GAMMA1*((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - DX/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AY - CY)*(-((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2) + DY/GAMMA2 - 
           GAMMA1*(bk + bl + GAMMA1)*
             ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
           GAMMA1*((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - DY/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AZ - CZ)*(-((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2) + DZ/GAMMA2 - 
           GAMMA1*(bk + bl + GAMMA1)*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
               (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
           GAMMA1*((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - DZ/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     ((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
        GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        2*((bk + bl + GAMMA1)*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
               (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
           ((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - DX/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
         (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     ((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
        GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        2*((bk + bl + GAMMA1)*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
               (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
           ((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - DY/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
         (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     2*((AX - CX)*((bk + bl + GAMMA1)*
             ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
           ((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - DX/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AY - CY)*((bk + bl + GAMMA1)*
             ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
           ((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - DY/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AZ - CZ)*((bk + bl + GAMMA1)*
             ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
           ((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - DZ/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     ((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
        GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        2*((bk + bl + GAMMA1)*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
               (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
           ((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - DZ/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
         (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP - 
     2*(((-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*
            (-((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2) + DX/GAMMA2) + 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*
             (-4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA1*(bk + bl + GAMMA1)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
               4*(bk + bl + GAMMA1)*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
               4*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4*(-((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2) + DX/GAMMA2)*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
           GAMMA3*(2.*GAMMA1*(bk + bl + GAMMA1)*
                 ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*
                 GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
               2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
                    (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4.*GAMMA1*((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - DX/GAMMA2)*
                 ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
                 GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
           (bk + bl + GAMMA1)*GAMMA3*(1 - 
               2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
                    (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
           ((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - DX/GAMMA2)*
             (4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA1*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4)*OVLAP + 
        ((-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*
            (-((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2) + DY/GAMMA2) + 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*
             (-4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA1*(bk + bl + GAMMA1)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
               4*(bk + bl + GAMMA1)*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
               4*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4*(-((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2) + DY/GAMMA2)*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
           GAMMA3*(2.*GAMMA1*(bk + bl + GAMMA1)*
                 ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*
                 GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
               2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
                    (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4.*GAMMA1*((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - DY/GAMMA2)*
                 ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
                 GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
           (bk + bl + GAMMA1)*GAMMA3*(1 - 
               2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
                    (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
           ((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - DY/GAMMA2)*
             (4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA1*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4)*OVLAP + 
        (((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*
             (-4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA1*(bk + bl + GAMMA1)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
               4*(bk + bl + GAMMA1)*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
               4*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4*(-((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2) + DZ/GAMMA2)*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
           GAMMA3*(2.*GAMMA1*(bk + bl + GAMMA1)*
                 ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*
                 GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
               2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
                    (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4.*GAMMA1*((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - DZ/GAMMA2)*
                 ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
                 GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
           (bk + bl + GAMMA1)*GAMMA3*(1 - 
               2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**
                  2*GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
           ((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - DZ/GAMMA2)*
             (4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA1*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4)*OVLAP))


;[o]
-4*ALPHA6*A1**2*(((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -         GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*((bk + bl + GAMMA1)*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -                (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -            ((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - DX/GAMMA2)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -      ((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -         GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*((bk + bl + GAMMA1)*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -                (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -            ((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - DY/GAMMA2)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -      ((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -         GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*((bk + bl + GAMMA1)*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -                (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -            ((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - DZ/GAMMA2)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP)/GAMMA3 - 
     -  4*ALPHA5*A2**2*(((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -         GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*((bk + bl + GAMMA1)*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -                (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -            ((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - DX/GAMMA2)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-BX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -      ((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -         GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*((bk + bl + GAMMA1)*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -                (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -            ((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - DY/GAMMA2)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-BY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -      ((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -         GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*((bk + bl + GAMMA1)*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -                (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -            ((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - DZ/GAMMA2)*GAMMA2*GAMMA3/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-BZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP)/GAMMA3 + 
     -  4*ALPHA5*ALPHA6*(2*(-((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2) + DX/GAMMA2 - 
     -        GAMMA1*(bk + bl + GAMMA1)*
     -          ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -        GAMMA1*((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - DX/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -      (-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2 + 
     -        GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     -     2*(-((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2) + DY/GAMMA2 - 
     -        GAMMA1*(bk + bl + GAMMA1)*
     -          ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -        GAMMA1*((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - DY/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -      (-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2 + 
     -        GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     -     2*(-((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2) + DZ/GAMMA2 - 
     -        GAMMA1*(bk + bl + GAMMA1)*
     -          ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -        GAMMA1*((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - DZ/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -      (-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2 + 
     -        GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     -     3*(-0.5/GAMMA2**2 + GAMMA1*(bk + bl + GAMMA1)*GAMMA3/
     -         (2*GAMMA2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -        GAMMA1*GAMMA3/(2*GAMMA2**2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP - 
     -     2*((AX - CX)*(-((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2) + DX/GAMMA2 - 
     -           GAMMA1*(bk + bl + GAMMA1)*
     -             ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -           GAMMA1*((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - DX/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AY - CY)*(-((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2) + DY/GAMMA2 - 
     -           GAMMA1*(bk + bl + GAMMA1)*
     -             ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -           GAMMA1*((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - DY/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AZ - CZ)*(-((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2) + DZ/GAMMA2 - 
     -           GAMMA1*(bk + bl + GAMMA1)*
     -             ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -           GAMMA1*((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - DZ/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -     ((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
     -         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -        GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        2*((bk + bl + GAMMA1)*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -               (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -           ((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - DX/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -         (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -     ((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
     -         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -        GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        2*((bk + bl + GAMMA1)*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -               (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -           ((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - DY/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -         (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -     2*((AX - CX)*((bk + bl + GAMMA1)*
     -             ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -             GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -           ((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - DX/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AY - CY)*((bk + bl + GAMMA1)*
     -             ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -             GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -           ((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - DY/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AZ - CZ)*((bk + bl + GAMMA1)*
     -             ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -             GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -           ((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - DZ/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -     ((bk + bl + GAMMA1)*GAMMA2*GAMMA3/
     -         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -        GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        2*((bk + bl + GAMMA1)*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -               (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -           ((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - DZ/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -         (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP - 
     -     2*(((-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*
     -            (-((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2) + DX/GAMMA2) + 
     -           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*
     -             (-4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA1*(bk + bl + GAMMA1)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -               4*(bk + bl + GAMMA1)*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -               4*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4*(-((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2) + DX/GAMMA2)*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
     -           GAMMA3*(2.*GAMMA1*(bk + bl + GAMMA1)*
     -                 ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*
     -                 GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -               2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -                    (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4.*GAMMA1*((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - DX/GAMMA2)*
     -                 ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -                 GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -           (bk + bl + GAMMA1)*GAMMA3*
     -             (1 - 2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -                    (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -           ((ALPHA3*CX + ALPHA4*DX)/GAMMA2**2 - DX/GAMMA2)*
     -             (4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA1*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4)*OVLAP + 
     -        ((-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*
     -            (-((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2) + DY/GAMMA2) + 
     -           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*
     -             (-4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA1*(bk + bl + GAMMA1)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -               4*(bk + bl + GAMMA1)*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -               4*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4*(-((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2) + DY/GAMMA2)*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
     -           GAMMA3*(2.*GAMMA1*(bk + bl + GAMMA1)*
     -                 ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*
     -                 GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -               2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -                    (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4.*GAMMA1*((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - DY/GAMMA2)*
     -                 ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -                 GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -           (bk + bl + GAMMA1)*GAMMA3*
     -             (1 - 2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -                    (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -           ((ALPHA3*CY + ALPHA4*DY)/GAMMA2**2 - DY/GAMMA2)*
     -             (4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA1*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4)*OVLAP + 
     -        (((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*
     -             (-4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA1*(bk + bl + GAMMA1)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -               4*(bk + bl + GAMMA1)*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -               4*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4*(-((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2) + DZ/GAMMA2)*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
     -           GAMMA3*(2.*GAMMA1*(bk + bl + GAMMA1)*
     -                 ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*
     -                 GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -               2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -                    (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4.*GAMMA1*((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - DZ/GAMMA2)*
     -                 ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -                 GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) - 
     -           (bk + bl + GAMMA1)*GAMMA3*
     -             (1 - 2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -                    (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -           ((ALPHA3*CZ + ALPHA4*DZ)/GAMMA2**2 - DZ/GAMMA2)*
     -             (4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA1*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4)*OVLAP))
:[font = input; startGroup; backColorRed = 65280; backColorGreen = 65280; backColorBlue = 65280; fontColorRed = 0; fontColorGreen = 0; fontColorBlue = 32768; plain; fontName = "Arial"; fontSize = 10; ]
FortranForm[ReplaceAll[D[TKIN,bl], {(al1*RAl + al2*RBl)/(al1+al2) -> P,(ak1*RAk + ak2*RBk)/(ak1+ak2) -> Q,
		ak1+ak2+al1+al2->GAMMA1 + GAMMA2,ak1+ak2+bk+bl->GAMMA3 + GAMMA2,
		al1+al2->GAMMA1,ak1+ak2->GAMMA2,RAk-RBk->RCD,RAl-RBl->RAB,RAl->RA,RBl->RB,
		RAk->RC,RBk->RD,bk+bl->GAMMA3,al1->ALPHA1,al2->ALPHA2,ak1->ALPHA3,
		ak2->ALPHA4,bl->ALPHA5,bk->ALPHA6
		}]]
:[font = output; inactive; formatted; output; endGroup; nowordwrap; ]
(4*ALPHA6*((AX - CX)**2 + (AY - CY)**2 + (AZ - CZ)**2) + 
     6*(-(ALPHA6*A1/GAMMA3**2) - ALPHA5*A2/GAMMA3**2 + A2/GAMMA3))*OVLAP - 
  4*ALPHA6*A1**2*((GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
         GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*(((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
              (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
      (GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
         GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*(((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
              (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
      (GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
         GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*(((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
              (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP)/GAMMA3 + 
  4*ALPHA6*A1**2*((0.5/GAMMA1 - GAMMA2*GAMMA3/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
      (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
      (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP)/GAMMA3**2 - 
  4*ALPHA5*A2**2*((GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
         GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*(((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
              (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-BX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
      (GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
         GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*(((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
              (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-BY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
      (GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
         GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*(((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
              (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-BZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP)/GAMMA3 + 
  4*ALPHA5*A2**2*((0.5/GAMMA1 - GAMMA2*GAMMA3/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         (-BX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
      (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         (-BY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
      (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         (-BZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP)/GAMMA3**2 - 
  4*A2**2*((0.5/GAMMA1 - GAMMA2*GAMMA3/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         (-BX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
      (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         (-BY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
      (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         (-BZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP)/GAMMA3 + 
  4*ALPHA5*ALPHA6*(3*(GAMMA1*(GAMMA1 + GAMMA2)*GAMMA3/
         (2*GAMMA2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
        GAMMA1/(2*GAMMA2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     2*(-(GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*
            (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
        GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
      (-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2 + 
        GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     2*(-(GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*
            (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
        GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
      (-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2 + 
        GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     2*(-(GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*
            (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
        GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
      (-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2 + 
        GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP - 
     2*((AX - CX)*(-(GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
                 (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*(GAMMA1 + GAMMA2)*GAMMA3/
              (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
           GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AY - CY)*(-(GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
                 (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*(GAMMA1 + GAMMA2)*GAMMA3/
              (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
           GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AZ - CZ)*(-(GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
                 (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*(GAMMA1 + GAMMA2)*GAMMA3/
              (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
           GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     2*((AX - CX)*(((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*
             GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AY - CY)*(((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*
             GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AZ - CZ)*(((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
             (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     (GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
        GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        2*(((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
             (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
         (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     (GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
        GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        2*(((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
             (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
         (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     (GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
        GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        2*(((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
             (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
         (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP - 
     2*((GAMMA3*(2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
                    (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
               2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
                    (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA2/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*
             (-4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA1*(GAMMA1 + GAMMA2)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
               4*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
               4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA1/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 - 
           (GAMMA1 + GAMMA2)*GAMMA3*(1 - 
               2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
                    (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
           (1 - 2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
                   (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA2*GAMMA3/
               (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
        (GAMMA3*(2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
                    (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
               2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
                    (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA2/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*
             (-4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA1*(GAMMA1 + GAMMA2)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
               4*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
               4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA1/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 - 
           (GAMMA1 + GAMMA2)*GAMMA3*(1 - 
               2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
                    (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
           (1 - 2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
                   (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA2*GAMMA3/
               (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
        (GAMMA3*(2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
                    (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
               2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
                    (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA2/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*
             (-4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA1*(GAMMA1 + GAMMA2)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
               4*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
               4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA1/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 - 
           (GAMMA1 + GAMMA2)*GAMMA3*(1 - 
               2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
                    (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
           (1 - 2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
                   (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA2*GAMMA3/
               (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP)) + 
  4*ALPHA6*((-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2 + 
         GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
          (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2*OVLAP + 
     (-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2 + 
         GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
          (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2*OVLAP + 
     (-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2 + 
         GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
          (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2*OVLAP + 
     3*(0.5/GAMMA2 - GAMMA1*GAMMA3/(2*GAMMA2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*
      OVLAP - 2*((AX - CX)*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2 + 
           GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AY - CY)*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2 + 
           GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AZ - CZ)*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2 + 
           GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
     (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
     2*((AX - CX)*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AY - CY)*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AZ - CZ)*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP - 
     2*(((-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2) + 
           GAMMA3*(1 - 2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
                    (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*
             (4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA1*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4)*OVLAP + 
        ((-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2) + 
           GAMMA3*(1 - 2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
                    (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*
             (4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA1*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4)*OVLAP + 
        OVLAP*(GAMMA3*(1 - 2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
                    (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*
             (4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA1*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
           (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*QCz)))


;[o]
(4*ALPHA6*((AX - CX)**2 + (AY - CY)**2 + (AZ - CZ)**2) + 
     -     6*(-(ALPHA6*A1/GAMMA3**2) - ALPHA5*A2/GAMMA3**2 + A2/GAMMA3))*OVLAP - 
     -  4*ALPHA6*A1**2*((GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -         GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*(((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -              (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -      (GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -         GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*(((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -              (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -      (GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -         GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*(((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -              (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP)/GAMMA3 + 
     -  4*ALPHA6*A1**2*((0.5/GAMMA1 - GAMMA2*GAMMA3/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
     -      (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
     -      (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP)/GAMMA3**2 - 
     -  4*ALPHA5*A2**2*((GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -         GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*(((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -              (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-BX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -      (GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -         GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*(((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -              (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-BY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -      (GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -         GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*(((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -              (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-BZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP)/GAMMA3 + 
     -  4*ALPHA5*A2**2*((0.5/GAMMA1 - GAMMA2*GAMMA3/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         (-BX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
     -      (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         (-BY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
     -      (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         (-BZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP)/GAMMA3**2 - 
     -  4*A2**2*((0.5/GAMMA1 - GAMMA2*GAMMA3/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         (-BX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
     -      (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         (-BY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
     -      (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         (-BZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP)/GAMMA3 + 
     -  4*ALPHA5*ALPHA6*(3*(GAMMA1*(GAMMA1 + GAMMA2)*GAMMA3/
     -         (2*GAMMA2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -        GAMMA1/(2*GAMMA2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -     2*(-(GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*
     -            (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -        GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -      (-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2 + 
     -        GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     -     2*(-(GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*
     -            (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -        GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -      (-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2 + 
     -        GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     -     2*(-(GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*
     -            (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -        GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -      (-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2 + 
     -        GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP - 
     -     2*((AX - CX)*(-(GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -                 (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*(GAMMA1 + GAMMA2)*GAMMA3/
     -              (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -           GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AY - CY)*(-(GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -                 (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*(GAMMA1 + GAMMA2)*GAMMA3/
     -              (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -           GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AZ - CZ)*(-(GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -                 (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*(GAMMA1 + GAMMA2)*GAMMA3/
     -              (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -           GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -     2*((AX - CX)*(((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*
     -             GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2\
     -             - ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AY - CY)*(((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*
     -             GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2\
     -             - ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AZ - CZ)*(((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*
     -             GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2\
     -             - ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -     (GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
     -         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -        GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        2*(((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -             (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -         (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -     (GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
     -         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -        GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        2*(((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -             (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -         (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -     (GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
     -         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -        GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        2*(((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -             (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -         (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP - 
     -     2*((GAMMA3*(2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -                    (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -               2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -                    (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA2/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*
     -             (-4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA1*(GAMMA1 + GAMMA2)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -               4*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -               4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA1/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 - 
     -           (GAMMA1 + GAMMA2)*GAMMA3*
     -             (1 - 2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -                    (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -           (1 - 2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -                   (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA2*GAMMA3/
     -               (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -        (GAMMA3*(2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -                    (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -               2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -                    (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA2/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*
     -             (-4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA1*(GAMMA1 + GAMMA2)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -               4*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -               4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA1/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 - 
     -           (GAMMA1 + GAMMA2)*GAMMA3*
     -             (1 - 2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -                    (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -           (1 - 2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -                   (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA2*GAMMA3/
     -               (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -        (GAMMA3*(2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -                    (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -               2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -                    (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA2/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*
     -             (-4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA1*(GAMMA1 + GAMMA2)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -               4*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -               4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA1/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 - 
     -           (GAMMA1 + GAMMA2)*GAMMA3*
     -             (1 - 2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -                    (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -           (1 - 2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -                   (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA2*GAMMA3/
     -               (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP)) + 
     -  4*ALPHA6*((-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2 + 
     -         GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
     -          (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2*OVLAP + 
     -     (-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2 + 
     -         GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
     -          (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2*OVLAP + 
     -     (-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2 + 
     -         GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
     -          (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2*OVLAP + 
     -     3*(0.5/GAMMA2 - GAMMA1*GAMMA3/(2*GAMMA2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*
     -      OVLAP - 2*((AX - CX)*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2 + 
     -           GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AY - CY)*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2 + 
     -           GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AZ - CZ)*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2 + 
     -           GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -     (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
     -     (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
     -     2*((AX - CX)*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AY - CY)*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AZ - CZ)*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -     (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP - 
     -     2*(((-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2) + 
     -           GAMMA3*(1 - 2.*GAMMA1*
     -                 ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*
     -                 GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*
     -             (4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA1*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4)*OVLAP + 
     -        ((-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2) + 
     -           GAMMA3*(1 - 2.*GAMMA1*
     -                 ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*
     -                 GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*
     -             (4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA1*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4)*OVLAP + 
     -        OVLAP*(GAMMA3*(1 - 2.*GAMMA1*
     -                 ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*
     -                 GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*
     -             (4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA1*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
     -           (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*QCz)))
:[font = input; startGroup; backColorRed = 65280; backColorGreen = 65280; backColorBlue = 65280; fontColorRed = 0; fontColorGreen = 0; fontColorBlue = 32768; plain; fontName = "Arial"; fontSize = 10; ]
FortranForm[ReplaceAll[D[TKIN,bk], {(al1*RAl + al2*RBl)/(al1+al2) -> P,(ak1*RAk + ak2*RBk)/(ak1+ak2) -> Q,
		ak1+ak2+al1+al2->GAMMA1 + GAMMA2,ak1+ak2+bk+bl->GAMMA3 + GAMMA2,
		al1+al2->GAMMA1,ak1+ak2->GAMMA2,RAk-RBk->RCD,RAl-RBl->RAB,RAl->RA,RBl->RB,
		RAk->RC,RBk->RD,bk+bl->GAMMA3,al1->ALPHA1,al2->ALPHA2,ak1->ALPHA3,
		ak2->ALPHA4,bl->ALPHA5,bk->ALPHA6
		}]]
:[font = output; inactive; formatted; output; endGroup; nowordwrap; ]
(4*ALPHA5*((AX - CX)**2 + (AY - CY)**2 + (AZ - CZ)**2) + 
     6*(-(ALPHA6*A1/GAMMA3**2) - ALPHA5*A2/GAMMA3**2 + A1/GAMMA3))*OVLAP - 
  4*ALPHA6*A1**2*((GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
         GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*(((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
              (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
      (GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
         GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*(((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
              (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
      (GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
         GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*(((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
              (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP)/GAMMA3 + 
  4*ALPHA6*A1**2*((0.5/GAMMA1 - GAMMA2*GAMMA3/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
      (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
      (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP)/GAMMA3**2 - 
  4*A1**2*((0.5/GAMMA1 - GAMMA2*GAMMA3/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
      (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
      (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP)/GAMMA3 - 
  4*ALPHA5*A2**2*((GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
         GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*(((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
              (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-BX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
      (GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
         GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*(((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
              (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-BY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
      (GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
         GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         2*(((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
              (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
          (-BZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP)/GAMMA3 + 
  4*ALPHA5*A2**2*((0.5/GAMMA1 - GAMMA2*GAMMA3/
          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         (-BX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
      (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         (-BY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
      (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
         (-BZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP)/GAMMA3**2 + 
  4*ALPHA5*ALPHA6*(3*(GAMMA1*(GAMMA1 + GAMMA2)*GAMMA3/
         (2*GAMMA2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
        GAMMA1/(2*GAMMA2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     2*(-(GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*
            (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
        GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
      (-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2 + 
        GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     2*(-(GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*
            (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
        GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
      (-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2 + 
        GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     2*(-(GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*
            (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
        GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
      (-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2 + 
        GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP - 
     2*((AX - CX)*(-(GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
                 (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*(GAMMA1 + GAMMA2)*GAMMA3/
              (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
           GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AY - CY)*(-(GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
                 (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*(GAMMA1 + GAMMA2)*GAMMA3/
              (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
           GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AZ - CZ)*(-(GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
                 (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*(GAMMA1 + GAMMA2)*GAMMA3/
              (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
           GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     2*((AX - CX)*(((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*
             GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AY - CY)*(((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*
             GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AZ - CZ)*(((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
             (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     (GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
        GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        2*(((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
             (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
         (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     (GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
        GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        2*(((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
             (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
         (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     (GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
        GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        2*(((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
             (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
         (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP - 
     2*((GAMMA3*(2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
                    (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
               2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
                    (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA2/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*
             (-4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA1*(GAMMA1 + GAMMA2)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
               4*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
               4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA1/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 - 
           (GAMMA1 + GAMMA2)*GAMMA3*(1 - 
               2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
                    (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
           (1 - 2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
                   (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA2*GAMMA3/
               (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
        (GAMMA3*(2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
                    (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
               2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
                    (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA2/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*
             (-4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA1*(GAMMA1 + GAMMA2)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
               4*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
               4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA1/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 - 
           (GAMMA1 + GAMMA2)*GAMMA3*(1 - 
               2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
                    (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
           (1 - 2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
                   (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA2*GAMMA3/
               (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
        (GAMMA3*(2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
                    (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
               2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
                    (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA2/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*
             (-4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA1*(GAMMA1 + GAMMA2)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
               4*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
               4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA1/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 - 
           (GAMMA1 + GAMMA2)*GAMMA3*(1 - 
               2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
                    (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
           (1 - 2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
                   (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA2*GAMMA3/
               (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP)) + 
  4*ALPHA5*((-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2 + 
         GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
          (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2*OVLAP + 
     (-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2 + 
         GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
          (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2*OVLAP + 
     (-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2 + 
         GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
          (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2*OVLAP + 
     3*(0.5/GAMMA2 - GAMMA1*GAMMA3/(2*GAMMA2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*
      OVLAP - 2*((AX - CX)*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2 + 
           GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AY - CY)*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2 + 
           GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AZ - CZ)*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2 + 
           GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
     (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
     2*((AX - CX)*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AY - CY)*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (AZ - CZ)*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
        (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP - 
     2*(((-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2) + 
           GAMMA3*(1 - 2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
                    (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*
             (4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA1*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4)*OVLAP + 
        ((-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2) + 
           GAMMA3*(1 - 2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
                    (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*
             (4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA1*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4)*OVLAP + 
        OVLAP*(GAMMA3*(1 - 2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
                    (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*
             (4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA1*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
               4*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
           (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*QCz)))


;[o]
(4*ALPHA5*((AX - CX)**2 + (AY - CY)**2 + (AZ - CZ)**2) + 
     -     6*(-(ALPHA6*A1/GAMMA3**2) - ALPHA5*A2/GAMMA3**2 + A1/GAMMA3))*OVLAP - 
     -  4*ALPHA6*A1**2*((GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -         GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*(((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -              (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -      (GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -         GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*(((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -              (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -      (GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -         GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*(((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -              (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP)/GAMMA3 + 
     -  4*ALPHA6*A1**2*((0.5/GAMMA1 - GAMMA2*GAMMA3/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
     -      (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
     -      (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP)/GAMMA3**2 - 
     -  4*A1**2*((0.5/GAMMA1 - GAMMA2*GAMMA3/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
     -      (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
     -      (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP)/GAMMA3 - 
     -  4*ALPHA5*A2**2*((GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -         GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*(((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -              (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-BX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -      (GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -         GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*(((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -              (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-BY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -      (GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -         GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         2*(((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -              (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2/
     -             (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -          (-BZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP)/GAMMA3 + 
     -  4*ALPHA5*A2**2*((0.5/GAMMA1 - GAMMA2*GAMMA3/
     -          (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         (-BX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -            ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
     -      (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         (-BY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -            ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
     -      (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -         (-BZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -            ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -              GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP)/GAMMA3**2 + 
     -  4*ALPHA5*ALPHA6*(3*(GAMMA1*(GAMMA1 + GAMMA2)*GAMMA3/
     -         (2*GAMMA2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -        GAMMA1/(2*GAMMA2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -     2*(-(GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*
     -            (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -        GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -      (-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2 + 
     -        GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     -     2*(-(GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*
     -            (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -        GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -      (-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2 + 
     -        GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP + 
     -     2*(-(GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*
     -            (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -        GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -      (-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2 + 
     -        GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
     -         (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*OVLAP - 
     -     2*((AX - CX)*(-(GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -                 (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*(GAMMA1 + GAMMA2)*GAMMA3/
     -              (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -           GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AY - CY)*(-(GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -                 (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*(GAMMA1 + GAMMA2)*GAMMA3/
     -              (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -           GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AZ - CZ)*(-(GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -                 (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*(GAMMA1 + GAMMA2)*GAMMA3/
     -              (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -           GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -     2*((AX - CX)*(((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*
     -             GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2\
     -             - ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AY - CY)*(((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*
     -             GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2\
     -             - ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AZ - CZ)*(((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*
     -             GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2\
     -             - ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -     (GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
     -         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -        GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        2*(((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*
     -             (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -         (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -     (GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
     -         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -        GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        2*(((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*
     -             (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -         (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -     (GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
     -         (2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) - 
     -        GAMMA2/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        2*(((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*
     -             (GAMMA1 + GAMMA2)*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))*
     -         (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP - 
     -     2*((GAMMA3*(2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -                    (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -               2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -                    (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA2/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*
     -             (-4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA1*(GAMMA1 + GAMMA2)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -               4*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -               4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA1/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 - 
     -           (GAMMA1 + GAMMA2)*GAMMA3*
     -             (1 - 2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -                    (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -           (1 - 2.*GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -                   (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*GAMMA2*GAMMA3/
     -               (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -        (GAMMA3*(2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -                    (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -               2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -                    (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA2/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*
     -             (-4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA1*(GAMMA1 + GAMMA2)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -               4*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -               4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA1/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 - 
     -           (GAMMA1 + GAMMA2)*GAMMA3*
     -             (1 - 2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -                    (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -           (1 - 2.*GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -                   (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*GAMMA2*GAMMA3/
     -               (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -        (GAMMA3*(2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -                    (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 - 
     -               2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -                    (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA2/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*
     -             (-4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA1*(GAMMA1 + GAMMA2)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -               4*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*(GAMMA1 + GAMMA2)*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2 + 
     -               4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA1/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 - 
     -           (GAMMA1 + GAMMA2)*GAMMA3*
     -             (1 - 2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -                    (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)**2) + 
     -           (1 - 2.*GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -                   (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*GAMMA2*GAMMA3/
     -               (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP)) + 
     -  4*ALPHA5*((-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2 + 
     -         GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
     -          (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2*OVLAP + 
     -     (-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2 + 
     -         GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
     -          (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2*OVLAP + 
     -     (-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2 + 
     -         GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
     -          (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2*OVLAP + 
     -     3*(0.5/GAMMA2 - GAMMA1*GAMMA3/(2*GAMMA2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*
     -      OVLAP - 2*((AX - CX)*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2 + 
     -           GAMMA1*((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AY - CY)*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2 + 
     -           GAMMA1*((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AZ - CZ)*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2 + 
     -           GAMMA1*((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -     (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
     -     (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP + 
     -     2*((AX - CX)*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1 - 
     -           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AY - CY)*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1 - 
     -           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (AZ - CZ)*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)))*OVLAP + 
     -     (0.5/GAMMA1 - GAMMA2*GAMMA3/(2*GAMMA1*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -        (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - 
     -           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
     -            (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))**2)*OVLAP - 
     -     2*(((-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2) + 
     -           GAMMA3*(1 - 2.*GAMMA1*
     -                 ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)**2*
     -                 GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -           ((ALPHA1*AX + ALPHA2*BX)/GAMMA1 - (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*
     -             (4*(-AX + (ALPHA1*AX + ALPHA2*BX)/GAMMA1)*GAMMA1*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4*(-CX + (ALPHA3*CX + ALPHA4*DX)/GAMMA2)*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4)*OVLAP + 
     -        ((-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2) + 
     -           GAMMA3*(1 - 2.*GAMMA1*
     -                 ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)**2*
     -                 GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -           ((ALPHA1*AY + ALPHA2*BY)/GAMMA1 - (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*
     -             (4*(-AY + (ALPHA1*AY + ALPHA2*BY)/GAMMA1)*GAMMA1*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4*(-CY + (ALPHA3*CY + ALPHA4*DY)/GAMMA2)*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4)*OVLAP + 
     -        OVLAP*(GAMMA3*(1 - 2.*GAMMA1*
     -                 ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)**2*
     -                 GAMMA2*GAMMA3/(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/
     -            (2*(GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3)) + 
     -           ((ALPHA1*AZ + ALPHA2*BZ)/GAMMA1 - (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*
     -             (4*(-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*GAMMA1*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3) - 
     -               4*(-CZ + (ALPHA3*CZ + ALPHA4*DZ)/GAMMA2)*GAMMA2*GAMMA3/
     -                (GAMMA1*GAMMA2 + (GAMMA1 + GAMMA2)*GAMMA3))/4 + 
     -           (-AZ + (ALPHA1*AZ + ALPHA2*BZ)/GAMMA1)*QCz)))
^*)