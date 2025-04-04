SUBROUTINE DTKINET(OVLAP,D1,D2,D3,D4,D5,					&
				   AZ,A1,BZ,A2,					&
                   CZ,A3,DZ,A4,A5,A6,			&
				   DT1,DT2,DT3,DT4,DT5)       

	USE scalars_to_be_shared

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	IMPLICIT INTEGER (I-N)


!--- FOUR CENTER KINETIC ENERGY INTEGRAL ---                     
!-----------------------------------------------------

      GAMA1=A1+A2
      GAMA2=A3+A4                                                      
      GAMA3=A5+A6

      QQQ=A1*AZ+A2*BZ
	  YYY=A3*CZ+A4*DZ
	  PZ=QQQ/GAMA1                                            
      QZ=YYY/GAMA2                                           

      QCZ=QZ-CZ                                 

      PAZ=PZ-AZ                                 

      PBZ=PZ-BZ                                 

      PQZ=PZ-QZ

      PQZ2=PQZ*PQZ                      

      G1PG2=GAMA1+GAMA2                                                
      G1RG2=GAMA1*GAMA2                                                
      XMIAN=G1RG2+GAMA3*G1PG2                                          
      DELT=G1RG2*GAMA3/XMIAN*FOUR                                      
      GA1DEL=DELT/GAMA1                                                

      OVZA11=PAZ-PQZ*GA1DEL/FOUR                 

      OVZB11=PBZ-PQZ*GA1DEL/FOUR                                     

      OVZA12=OVZA11*OVZA11                                             

      POM=HALF/GAMA1-GA1DEL/EIGHT/GAMA1                                

      OVXA12=POM*OVLAP                                        
      OVYA12=POM*OVLAP                                        

      OVZA12=(OVZA12+POM)*OVLAP                                        

      OVZB12=OVZB11*OVZB11                                             

      OVXB12=POM*OVLAP                                         
      OVYB12=POM*OVLAP                                         

      OVZB12=(OVZB12+POM)*OVLAP                                        

      GA2DEL=DELT/GAMA2                                                

      OVZC21=QCZ+PQZ*GA2DEL/FOUR                                        

      OVZC22=OVZC21*OVZC21*OVLAP                                        

      POM=HALF/GAMA2-GA2DEL/EIGHT/GAMA2                                 
      POM=POM*OVLAP                                                     

      OVXC22=POM                                                 
      OVYC22=POM                                                 

      OVZC22=OVZC22+POM                                                 

      POM=GA2DEL/EIGHT/GAMA1                                            

      ZZZ=ONE-HALF*PQZ2*DELT                                         

      XA1XC2=POM
      YA1YC2=POM

      ZA1ZC2=POM*ZZZ                                                 

	  FPT4=PAZ*GA2DEL-QCZ*GA1DEL

	  FPT5=PAZ*QCZ

      ZA1ZC2=FPT5+PQZ/FOUR*FPT4+ZA1ZC2            

      XA1XC2=XA1XC2*OVLAP                                               
      YA1YC2=YA1YC2*OVLAP                                               
      ZA1ZC2=ZA1ZC2*OVLAP                                               

      ACZ=AZ-CZ                                 

      OVA12=OVXA12+OVYA12+OVZA12                                        

	  FTP6=OVZA11*ACZ

      OVA11=TWO*FTP6*OVLAP                

      OVC22=OVXC22+OVYC22+OVZC22                                       

	  FTP7=OVZC21*ACZ

      OVC21=-TWO*FTP7*OVLAP              

      A1C2=-TWO*(XA1XC2+YA1YC2+ZA1ZC2)                                 

	  FPT0=A5*A6

	  FPT2=FOUR*FPT0

	  FPT3=OVA12+OVA11+OVC22+OVC21+A1C2

      IF(GAMA3 .EQ. 0) THEN                                               
			A6DG3=HALF                                                       
			A5DG3=HALF                                                       
		ELSE
			A6DG3=A6/GAMA3                                                   
			A5DG3=A5/GAMA3
	  END IF

      OVB12=OVXB12+OVYB12+OVZB12                                       

	  FPT1=A1*A6DG3+A2*A5DG3

	  FTP10=SIX*FPT1

!-----------------------------------------------------
!DERIVATIVES:
!-----------------------------------------------------
!DTKDA1:
!======================================================

      t1 = ONE/GAMA1
      t2 = HALF*t1
	  t3 = EIGHTH*t1
      t5 = GA1DEL*t3
	  t6 = (t2-t5)
      t7 = t6*D1
      t8 = GAMA1*GAMA1
      t9 = ONE/t8
      t10 = HALF*t9
      t11 = GAMA2*GAMA3
      t12 = GAMA2+GAMA3
      t16 = XMIAN*XMIAN
      t17 = ONE/t16
      t19 = t11*t12*HALF*t1*t17
      t21 = ONE/XMIAN
      t23 = t11*HALF*t9*t21
	  t24 = -t10+t19+t23
      t25 = OVLAP*t24
      t26 = OVZA11*OVZA11
	  t27 = (t26+t2-t5)
      t28 = t27*D1
      t34 = AZ-QQQ*t1
      t36 = PQZ*FOUR
      t38 = GAMA3*t12
      t39 = t38*t17
      t45 = t34*t1-QUARTER*(-t36*GAMA2*t39+GA1DEL*t34*t1)
      t48 = OVLAP*(TWO*OVZA11*t45+t24)
      t55 = GAMA1*GAMA3
      t57 = ONE/GAMA2
      t66 = HALF*t57
      t60 = t55*t12*t66*t17
      t63 = GAMA3*t66*t21
	  t67 = GA2DEL*EIGHTH
      t68 = t67*t57
      t80 = QUARTER*(t36*t21*(-t55*t12*t21+GAMA3)+GA2DEL*t34*t1)
      t84 = OVZC21*OVZC21
      t103 = AZ*t1-QQQ*t9
      t115 = FOUR*GAMA3
      t127 = TWO*GAMA3
      t129 = G1RG2*PQZ2
	  t135 = GA2DEL*t3
	  t140 = TWO*t135
      ttt39 = QUARTER*PQZ
	  t150 = FPT5+ttt39*FPT4
	  t160 = t66-t68
	  t161 = t150+ZZZ*t135
	  t162 = (t84+t160)
	  t163 = TWO*t160

	  t164 = TWO*OVZC21
	  t166 = PAZ*t21
	  t167 = t129*t17
	  t168 = t129*t21
	  t169 = HALF*t21
	  t170 = HALF*t17
	  t171 = FOUR*t55
	  t172 = PAZ*t17
	  t173 = QCZ*t17
	  t175 = ONE-t127*t168
	  t177 = PQZ2*t21
	  t178 = t115*G1RG2*PQZ

      t165 = TWO*t7+TWO*t25+t28+t48+TWO*(t45*ACZ*OVLAP+D1*FTP6)+TWO*OVLAP*(t60-t63)+t163*	&
	  D1+OVLAP*(t164*t80+t60-t63)+t162*D1-TWO*(t80*ACZ*OVLAP+D1*FTP7)-TWO*(-OVLAP*t39+t140*	&
	  D1+OVLAP*(QCZ*t34*t1+QUARTER*(PQZ*(GA2DEL*t103-t171*t12*t172+FOUR*t11*t12*t173+t115*	&
	  t166)+FPT4*t34*t1)-t38*(t175)*t170+GAMA3*(t127*t12*t167-TWO*t11*t177-t178*t103*t21)*	&
	  t169)+t161*D1)

      t174 = A1*A1
      t176 = ONE/GAMA3
      t183 = OVZB11*OVZB11
      t188 = A2*A2
      t197 = ACZ*ACZ
	  t200 = t188*A5*t176
	  t201 = FPT2*t197
	  t202 = t183+t2-t5
	  t203 = FOUR*t174*A6

	  t1010 = t203*t176
	  t1030 = FOUR*t200
	  t1040 = TWO*OVZB11

      DT1 = t165*FPT2-EIGHT*A1*A6DG3*OVA12-(TWO*t7+TWO*t25+t28+t48)*t1010-(TWO*t25+TWO*t7+&
	  OVLAP*(t1040*t45+t24)+(t202)*D1)*t1030+(SIX*(FPT1)+t201)*D1+OVLAP*SIX*A6DG3


!DTKDA2:
!======================================================

      tt7 = t6*D2
      tt9 = GAMA1+GAMA3
      tt16 = t11*tt9*HALF*t1*t17
      tt19 = t1*t21
      tt20 = GAMA3*HALF*tt19
      tt22 = OVLAP*(tt16-tt20)
      tt25 = t27*D2
      tt29 = t11*tt9*t21
      tt34 = CZ*t57
      tt38 = GAMA2*GAMA2
      tt39 = ONE/tt38
      tt40 = YYY*tt39
      tt41 = -tt34+tt40
      tt44 = QUARTER*(t36*t21*(-tt29+GAMA3)+GA1DEL*tt41)
      tt47 = OVLAP*(-TWO*OVZA11*tt44+tt16-tt20)
      tt54 = HALF*tt39
      tt59 = t55*tt9*HALF*t57*t17
      tt62 = t55*HALF*tt39*t21
      tt72 = GAMA3*tt9
      tt78 = tt34-tt40+QUARTER*(-t36*GAMA1*tt72*t17+GA2DEL*tt41)
      tt100 = tt34-tt40
      tt107 = FOUR*GAMA2
	  tt160 = QCZ*t21
      
	  tt170 = EIGHTH*GA2DEL*t1

	  tt165 = TWO*tt7+TWO*tt22+tt25+tt47+TWO*(-tt44*ACZ*OVLAP+D2*FTP6)+TWO*OVLAP*(-tt54+	&
	  tt59+tt62)+t163*D2+OVLAP*(t164*tt78-tt54+tt59+tt62)+t162*D2-TWO*(tt78*ACZ*OVLAP+D2*	&
	  FTP7)-TWO*(TWO*(-OVLAP*GAMA3*tt9*t170+tt170*D2)+OVLAP*(PAZ*tt100+QUARTER*(PQZ*(-t171*	&
	  tt9*t172+tt107*GAMA3*tt9*t173-t115*tt160-tt107*GAMA3*tt100*t21)+FPT4*tt41)-tt72*t175*	&
	  t170+GAMA3*(t127*tt9*t167-TWO*t55*t177-t178*tt41*t21)*t169)+(t150+ZZZ*tt170)*D2)

      tt177 = tt29-GAMA3

	  tt200 = FOUR*t188*A5*t176
	  tt201 = FTP10+t201

      DT2 = tt165*FPT2-(TWO*tt7+TWO*tt22+tt25+tt47)*t1010-(OVLAP*tt19*tt177+TWO*	&
	  tt7+OVLAP*(-t1040*tt44+t2*t21*tt177)+t202*D2)*tt200+tt201*D2


!DTKDA5:
!======================================================

      ttt10 = t6*D3
	  ttt11 = t11*G1PG2
      ttt18 = ttt11*HALF*t1*t17
      ttt22 = GAMA2*HALF*t1*t21
      ttt24 = OVLAP*(ttt18-ttt22)
      ttt27 = t27*D3
      ttt31 = G1PG2*t21
      ttt34 = t21*(-t11*ttt31+GAMA2)
      ttt35 = t36*ttt34
      ttt38 = OVLAP*(-TWO*OVZA11*QUARTER*ttt35+ttt18-ttt22)
      ttt40 = ttt39*FOUR
      ttt41 = ACZ*OVLAP
	  ttt42 = t55*G1PG2
      ttt53 = ttt42*t66*t17
      ttt56 = GAMA1*t66*t21
      ttt68 = t21*(-t55*ttt31+GAMA1)
      ttt82 = GAMA3*G1PG2
      
	  ttt108 = t175*HALF
      
	  t3000 = TWO*t67*t1

	  ttt139 = TWO*ttt10+TWO*ttt24+ttt27+ttt38+TWO*(-ttt40*ttt34*ttt41+D3*FTP6)+TWO*OVLAP*	&
	  (ttt53-ttt56)+t163*D3+OVLAP*(t164*QUARTER*t36*ttt68+ttt53-ttt56)+t162*D3-TWO*(ttt40*	&
	  ttt68*ttt41+D3*FTP7)-TWO*(TWO*OVLAP*(-ttt82*t170+t169)+t3000*D3+OVLAP*(ttt39*FOUR*t21*&
	  (-ttt42*t166+ttt11*tt160+GAMA1*PAZ-GAMA2*QCZ)-ttt82*ttt108*t17+GAMA3*(t127*G1PG2*t167-&
	  TWO*t168)*t169+ttt108*t21)+t161*D3)

      ttt145 = GAMA3*GAMA3
      ttt146 = ONE/ttt145
      
	  DT3 = FOUR*A6*FPT3+ttt139*FPT2+t203*ttt146*OVA12-(TWO*ttt10+TWO*ttt24+ttt27+ttt38)*	&
	  t1010-FOUR*t176*(-t200+t188)*OVB12-(TWO*ttt24+TWO*ttt10+OVLAP*(-t1040*QUARTER*ttt35+	&
	  ttt18-ttt22)+t202*D3)*t1030+tt201*D3+OVLAP*(FOUR*t197*A6+SIX*(-A2*A5*ttt146-A1*A6*	&
	  ttt146+A2*t176))

				   
!DTKDAZ:
!======================================================

      tttt7 = t6*D4
      tttt10 = t27*D4
      tttt11 = OVLAP*TWO
      tttt12 = A1*t1
	  tttt13 = QUARTER*GA1DEL
      tttt14 = tttt13*tttt12
      tttt51 = tttt12-ONE
      tttt15 = tttt51-tttt14
      tttt17 = tttt11*OVZA11*tttt15
      tttt33 = QUARTER*GA2DEL

	  t1000 = OVLAP*EIGHT*ACZ*FPT0
	  t2000 = PQZ/t16

      DT4 = (TWO*tttt7+tttt10+tttt17+TWO*((tttt15*ACZ+OVZA11)*OVLAP+D4*FTP6)+t163*D4+tttt11*&
	  OVZC21*tttt33*tttt12+t162*D4-TWO*((tttt33*tttt12*ACZ+OVZC21)*OVLAP+D4*FTP7)-TWO*(t140*&
	  D4+OVLAP*(QCZ*tttt51+QUARTER*(PQZ*tttt51*GA2DEL+FPT4*tttt12)-TWO*A1*GAMA2*ttt145*		&
	  t2000)+t161*D4))*FPT2-(TWO*tttt7+tttt10+tttt17)*t1010-(TWO*tttt7+tttt11*OVZB11*		&
	  (tttt12-tttt14)+t202*D4)*tt200+tt201*D4+t1000


!DTKDBZ:
!======================================================

      tnt8 = TWO*t6*D5
      tnt11 = t27*D5
      tnt16 = A3*t57
      tnt17 = tttt13*tnt16
      tnt18 = tttt11*OVZA11*tnt17
      tnt50 = tnt16-ONE
      tnt35 = tnt50-QUARTER*GA2DEL*tnt16

      DT5 = FOUR*(tnt8+tnt11+tnt18+TWO*((tttt13*tnt16*ACZ-OVZA11)*OVLAP+D5*FTP6)+t163*D5+	&
	  tttt11*OVZC21*tnt35+t162*D5-TWO*((tnt35*ACZ-OVZC21)*OVLAP+D5*FTP7)-TWO*(t3000*D5+		&
	  OVLAP*(PAZ*tnt50+QUARTER*(-t36*tnt50*t11*t21-FPT4*tnt16)+TWO*A3*GAMA1*ttt145*t2000)+	&
	  t161*D5))*FPT0-(tnt8+tnt11+tnt18)*t1010-(tnt8+tttt11*OVZB11*tnt17+t202*D5)*tt200+		&
	  tt201*D5-t1000


END SUBROUTINE DTKINET