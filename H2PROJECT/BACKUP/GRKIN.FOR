C**************************************************
      SUBROUTINE grkin(OVLAP,L,K,R1LZ,AL1,R1KZ,Ak1,                      
     1              R2LZ,AL2,R2KZ,AK2,BL,BK,
     2              DTA1,DTA2,DTB,DTR1,DTR2) 
C**************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DELTA(100,100)
      DATA ONE,TWO,FOUR,SIX,EIGHT,HALF/1.,2.,4.,6.,8.,.5/               
      PI=3.1415926535898D+00

C ***                                                               ***
C ***    SUBROUTINE  CALCULATES FIRST DERIVATIVES OF KINETIC ENERGY ***
C ***                   
     
      DO 20 IL=1,100
      DO 20 IK=1,100
      DELTA(IL,IK)=0.0
20    CONTINUE
      DO 11 II=1,100
      DELTA(II,II)=1.0 
11    CONTINUE
      ALK1=AL1+Ak1                                                       
      ALK2=AL2+Ak2                                                       
      BLK12=BL+BK  
      DETAB=ALK1*ALK2+(ALK1+ALK2)*BLK12
	DETA1=(1.0+DELTA(L,K))*(ALK2 +BLK12)
	DETA2=(1.0+DELTA(L,K))*(ALK1 +BLK12)
	DETB=(1.0+DELTA(L,K))*(ALK1 +ALK2)
      RZ1=(AL1*R1LZ+Ak1*R1KZ)/ALK1
      RZ1A1= ( DELTA(L,K)*R1LZ+R1KZ )/ALK1
     1     -( AL1*R1LZ+Ak1*R1KZ )*(1.0+DELTA(L,K))/ALK1**2.0
      RZ1R1=(AL1*DELTA(L,K) + AK1 )/ALK1	
C
      RZ2=(AL2*R2LZ+Ak2*R2KZ)/ALK2
      RZ2A2=( DELTA(L,K)*R2LZ+R2KZ )/ALK2
     1     -( AL2*R2LZ+Ak2*R2KZ )*(1.0+DELTA(L,K))/ALK2**2.0
      RZ2R2=(AL2*DELTA(L,K) + AK2 )/ALK2	
      CALL grover(L,K,R1LZ,AL1,R1KZ,Ak1,                      
     1              R2LZ,AL2,R2KZ,AK2,BL,BK,
     2              DSA1,DSA2,DSB,DSR1,DSR2) 
C
      OVZA11=(RZ1-R1LZ)-(RZ1-RZ2)*ALK2*BLK12/DETAB                
      ZA11A1=RZ1A1 -( RZ1A1*ALK2*BLK12/DETAB 
     1  -(RZ1-RZ2)*ALK2*BLK12*DETA1/DETAB**2.0 )
      ZA11A2=- ( -RZ2A2*ALK2*BLK12/DETAB 
     1           +(RZ1-RZ2)*(1.0+DELTA(L,K))*BLK12/DETAB
     2           -(RZ1-RZ2)*ALK2*BLK12*DETA2/DETAB**2.0 )
      ZA11B =- ( (RZ1-RZ2)*ALK2*(1.0+DELTA(L,K))/DETAB
     1           -(RZ1-RZ2)*ALK2*BLK12*DETB /DETAB**2.0 )
      ZA11R1=(RZ1R1-DELTA(L,K) ) -RZ1R1*ALK2*BLK12/DETAB
      ZA11R2= RZ2R2*ALK2*BLK12/DETAB 
c
      OVZB11=(RZ1-R1KZ)-(RZ1-RZ2)*ALK2*BLK12/DETAB   
      ZB11A1=ZA11A1
      ZB11A2=ZA11A2
      ZB11B =ZA11B
      ZB11R1=(RZ1R1-1.0 ) -RZ1R1*ALK2*BLK12/DETAB
      ZB11R2=ZA11R2
C
      POM=0.5*(1.0/Alk1-ALK2*BLK12/(DETAB*ALK1) )                         
      POMA1=0.5* (
     1            -(1.0+DELTA(L,K))/ALK1**2
     2            +ALK2*BLK12/(DETAB*ALK1)**2.0 
     3            *(DETA1*ALK1+(1.0+DELTA(L,K))*DETAB )  )
      POMA2=0.5* ( -BLK12/ALK1
     1           *( (1.0+DELTA(L,K))*DETAB -ALK2*DETA2 )/DETAB**2.0  )
      POMB =0.5* ( -ALK2/ALK1
     1           *( (1.0+DELTA(L,K))*DETAB -BLK12*DETB )/DETAB**2.0  )
C
      OVZA12=(OVZA11*OVZA11+POM)*OVLAP                                        
      OVA12= (OVZA11**2.0+3.0*POM)*OVLAP
      DA12A1= (2.0*OVZA11*ZA11A1+3.0*POMA1)*OVLAP
     1      +(OVZA11**2.0+3.0*POM)*DSA1
      DA12A2= (2.0*OVZA11*ZA11A2+3.0*POMA2)*OVLAP
     1      +(OVZA11**2.0+3.0*POM)*DSA2
      DA12B = (2.0*OVZA11*ZA11B+3.0*POMB)*OVLAP
     1      +(OVZA11**2.0+3.0*POM)*DSB
      DA12R1= (2.0*OVZA11*ZA11R1)*OVLAP
     1      +(OVZA11**2.0+3.0*POM)*DSR1
      DA12R2= (2.0*OVZA11*ZA11R2)*OVLAP
     1      +(OVZA11**2.0+3.0*POM)*DSR2

C
      OVA11=2*(OVZA11*(R1LZ-R2LZ))*OVLAP                
      DA11A1=2*(ZA11A1*(R1LZ-R2LZ))*OVLAP 
     1      +2*(OVZA11*(R1LZ-R2LZ))*DSA1
      DA11A2=2*(ZA11A2*(R1LZ-R2LZ))*OVLAP 
     1      +2*(OVZA11*(R1LZ-R2LZ))*DSA2
      DA11B =2*(ZA11B *(R1LZ-R2LZ))*OVLAP 
     1      +2*(OVZA11*(R1LZ-R2LZ))*DSB
      DA11R1=2*(ZA11R1*(R1LZ-R2LZ)+ OVZA11*DELTA(L,K))*OVLAP 
     1      +2*(OVZA11*(R1LZ-R2LZ))*DSR1
      DA11R2=2*(ZA11R2*(R1LZ-R2LZ)- OVZA11*DELTA(L,K))*OVLAP 
     1      +2*(OVZA11*(R1LZ-R2LZ))*DSR2
C
      OVZB12=(OVZB11*OVZB11+POM)*OVLAP                                       
      OVB12= (OVZB11**2.0+3.0*POM)*OVLAP 
      DB12A1= (2.0*OVZB11*ZB11A1+3.0*POMA1)*OVLAP
     1      +(OVZB11**2.0+3.0*POM)*DSA1
      DB12A2= (2.0*OVZB11*ZB11A2+3.0*POMA2)*OVLAP
     1      +(OVZB11**2.0+3.0*POM)*DSA2
      DB12B = (2.0*OVZB11*ZB11B+3.0*POMB)*OVLAP
     1      +(OVZB11**2.0+3.0*POM)*DSB
      DB12R1= (2.0*OVZB11*ZB11R1)*OVLAP
     1      +(OVZB11**2.0+3.0*POM)*DSR1
      DB12R2= (2.0*OVZB11*ZB11R2)*OVLAP
     1      +(OVZB11**2.0+3.0*POM)*DSR2
C
      OVZC21=(RZ2-R2LZ)+(RZ1-RZ2)*ALK1*BLK12/DETAB
      ZC21A2=RZ2A2 +( -RZ2A2*ALK1*BLK12/DETAB 
     1  -(RZ1-RZ2)*ALK1*BLK12*DETA2/DETAB**2.0 )
      ZC21A1=+ ( RZ1A1*ALK1*BLK12/DETAB 
     1           +(RZ1-RZ2)*(1.0+DELTA(L,K))*BLK12/DETAB
     2           -(RZ1-RZ2)*ALK1*BLK12*DETA1/DETAB**2.0 )
      ZC21B =+ ( (RZ1-RZ2)*ALK1*(1.0+DELTA(L,K))/DETAB
     1           -(RZ1-RZ2)*ALK1*BLK12*DETB /DETAB**2.0 )
      ZC21R2=(RZ2R2-DELTA(L,K) ) -RZ2R2*ALK1*BLK12/DETAB
      ZC21R1= RZ1R1*ALK1*BLK12/DETAB 
C              
      OVC21=-2*(OVZC21*(R1LZ-R2LZ))*OVLAP              
      DC21A1=-2*(ZC21A1*(R1LZ-R2LZ))*OVLAP 
     1       -2*(OVZC21*(R1LZ-R2LZ))*DSA1
      DC21A2=-2*(ZC21A2*(R1LZ-R2LZ))*OVLAP 
     1       -2*(OVZC21*(R1LZ-R2LZ))*DSA2
      DC21B =-2*(ZC21B *(R1LZ-R2LZ))*OVLAP 
     1       -2*(OVZC21*(R1LZ-R2LZ))*DSB
      DC21R1=-2*(ZC21R1*(R1LZ-R2LZ)+ OVZC21*DELTA(L,K))*OVLAP 
     1       -2*(OVZC21*(R1LZ-R2LZ))*DSR1
      DC21R2=-2*(ZC21R2*(R1LZ-R2LZ)- OVZC21*DELTA(L,K))*OVLAP 
     1       -2*(OVZC21*(R1LZ-R2LZ))*DSR2
C
      POM=0.5*(1.0/ALK2-ALK1*BLK12/(DETAB*ALK2) )                       
      POMA2=0.5* (
     1            -(1.0+DELTA(L,K))/ALK2**2
     2            +ALK1*BLK12/(DETAB*ALK2)**2.0 
     3            *(DETA2*ALK2+(1.0+DELTA(L,K))*DETAB )  )
      POMA1=0.5* ( -BLK12/ALK2
     1           *( (1.0+DELTA(L,K))*DETAB -ALK1*DETA1 )/DETAB**2.0  )
      POMB =0.5* ( -ALK1/ALK2
     1           *( (1.0+DELTA(L,K))*DETAB -BLK12*DETB )/DETAB**2.0  )
C
      OVZC22=(OVZC21*OVZC21+POM)*OVLAP                       
      OVC22=(OVZC21**2.0+3.0*POM)*OVLAP
      DC22A1= (2.0*OVZC21*ZC21A1+3.0*POMA1)*OVLAP
     1      +(OVZC21**2.0+3.0*POM)*DSA1
      DC22A2= (2.0*OVZC21*ZC21A2+3.0*POMA2)*OVLAP
     1      +(OVZC21**2.0+3.0*POM)*DSA2
      DC22B = (2.0*OVZC21*ZC21B+3.0*POMB)*OVLAP
     1      +(OVZC21**2.0+3.0*POM)*DSB
      DC22R1= (2.0*OVZC21*ZC21R1)*OVLAP
     1      +(OVZC21**2.0+3.0*POM)*DSR1
      DC22R2= (2.0*OVZC21*ZC21R2)*OVLAP
     1      +(OVZC21**2.0+3.0*POM)*DSR2
C
      XX1=(0.5*BLK12/DETAB)
     1  *(1.0-2.0*(RZ1-RZ2)**2.0*ALK1*ALK2*BLK12/DETAB)
C
      YY1=0.5*BLK12/DETAB
      YY1A1=-0.5*BLK12*DETA1/DETAB**2
      YY1A2=-0.5*BLK12*DETA2/DETAB**2
      YY1B =0.5*( (1+DELTA(L,K))*DETAB -BLK12*DETB )/DETAB**2
C
      YY2=1.0-2.0*(RZ1-RZ2)**2.0*ALK1*ALK2*BLK12/DETAB
      YY2A1= -4.0*(RZ1-RZ2)*RZ1A1*ALK1*ALK2*BLK12/DETAB
     1       -2.0*(RZ1-RZ2)**2.0*(1.0+DELTA(L,K))*ALK2*BLK12/DETAB
     2       +2.0*(RZ1-RZ2)**2.0*ALK1*ALK2*BLK12*DETA1/DETAB**2.0
      YY2A2=  4.0*(RZ1-RZ2)*RZ2A2*ALK1*ALK2*BLK12/DETAB
     1       -2.0*(RZ1-RZ2)**2.0*(1.0+DELTA(L,K))*ALK1*BLK12/DETAB
     2       +2.0*(RZ1-RZ2)**2.0*ALK1*ALK2*BLK12*DETA2/DETAB**2.0
      YY2B = -2.0*(RZ1-RZ2)**2.0*(1.0+DELTA(L,K))*ALK1*ALK2/DETAB
     1       +2.0*(RZ1-RZ2)**2.0*ALK1*ALK2*BLK12*DETB/DETAB**2.0
      YY2R1= -4.0*(RZ1-RZ2)*RZ1R1*ALK1*ALK2*BLK12/DETAB
      YY2R2=  4.0*(RZ1-RZ2)*RZ2R2*ALK1*ALK2*BLK12/DETAB
C
      XX1A1=YY1A1*YY2+YY1*YY2A1
      XX1A2=YY1A2*YY2+YY1*YY2A2
      XX1B=YY1B*YY2+YY1*YY2B
      XX1R1=YY2R1*YY1
      XX1R2=YY2R2*YY1
C 
      XTEM=(RZ1-R1LZ)*ALK1*BLK12-(RZ2-R2LZ)*ALK2*BLK12
      XTEMA1=RZ1A1*ALK1*BLK12 +(RZ1-R1LZ)*(1.0+DELTA(L,K))*BLK12
      XTEMA2=-RZ2A2*ALK2*BLK12-(RZ2-R2LZ)*(1.0+DELTA(L,K))*BLK12
      XTEMB=(RZ1-R1LZ)*ALK1*(1.0+DELTA(L,K))
     1     -(RZ2-R2LZ)*ALK2*(1.0+DELTA(L,K))
      XTEMR1= (RZ1R1-DELTA(L,K))*ALK1*BLK12
      XTEMR2=-(RZ2R2-DELTA(L,K))*ALK2*BLK12
      XX2=(RZ1-R1LZ)*(RZ2-R2LZ)+(RZ1-RZ2)*XTEM/DETAB           
      XX2A1=RZ1A1*(RZ2-R2LZ)
     1     +RZ1A1*XTEM/DETAB+(RZ1-RZ2)*XTEMA1/DETAB
     2     -(RZ1-RZ2)*XTEM*DETA1/DETAB**2
      XX2A2=RZ2A2*(RZ1-R1LZ)
     1     -RZ2A2*XTEM/DETAB+(RZ1-RZ2)*XTEMA2/DETAB
     2     -(RZ1-RZ2)*XTEM*DETA2/DETAB**2
      XX2B=(RZ1-RZ2)*XTEMB/DETAB
     1     -(RZ1-RZ2)*XTEM*DETB/DETAB**2
      XX2R1=(RZ1R1-DELTA(L,K))*(RZ2-R2LZ) 
     1     +(RZ1R1*XTEM+(RZ1-RZ2)*XTEMR1)/DETAB
      XX2R2=(RZ2R2-DELTA(L,K))*(RZ1-R1LZ)
     1     +(-RZ2R2*XTEM+(RZ1-RZ2)*XTEMR2)/DETAB

      za1zc2=(XX1+XX2)*OVLAP
c     WRITE(6,*) 'JOHN,ZA1ZC2'
c     WRITE(6,*)  ZA1ZC2
	A1C2=-2*(XX1+XX2+BLK12/DETAB)*OVLAP                                 
      A1C2A1=-2*(
     1           (XX1A1+XX2A1-BLK12*DETA1/DETAB**2.0)*OVLAP
     2          +(XX1+XX2+BLK12/DETAB)*DSA1   )
      A1C2A2=-2*(
     1           (XX1A2+XX2A2-BLK12*DETA2/DETAB**2.0)*OVLAP
     2          +(XX1+XX2+BLK12/DETAB)*DSA2   )
      A1C2B =-2*(
     1     (XX1B+XX2B+(1.0+DELTA(L,K))/DETAB-BLK12*DETB/DETAB**2.0)*OVLAP
     2          +(XX1+XX2+BLK12/DETAB)*DSB   )
      A1C2R1=-2*(
     1           (XX1R1+XX2R1)*OVLAP
     2          +(XX1+XX2+BLK12/DETAB)*DSR1   )
      A1C2R2=-2*(
     1           (XX1R2+XX2R2)*OVLAP
     2          +(XX1+XX2+BLK12/DETAB)*DSR2   )
C
      TKIN=4*BL*BK*(OVA11 +OVA12+OVC22+OVC21+A1C2)                   
      DTA1=4*BL*BK*(DA11A1 +DA12A1+A1C2A1+DC22A1+DC21A1)                   
      DTA2=4*BL*BK*(DA11A2 +DA12A2+A1C2A2+DC22A2+DC21A2)                   
      DTR1=4*BL*BK*(DA11R1 +DA12R1+DC22R1+DC21R1+A1C2R1)                   
      DTR2=4*BL*BK*(DA11R2 +DA12R2+DC22R2+DC21R2+A1C2R2)                   
      DTB =4*(BL+DELTA(L,K)*BK)*(OVA11 +OVA12+OVC22+OVC21+A1C2)
     1    +4*BL*BK*(DA11B +DA12B+DC22B+DC21B+A1C2B)
      IF(BLK12) 10,12,10                                               
  12  A6DG3=0.5                                                     
      A5DG3=0.5                                                       
      GO TO 14                                                         
  10  A6DG3=BK/BLK12                                                   
      A5DG3=BL/BLK12                                                   
      A6DG3B=1/BLK12-BK*(1.0+DELTA(L,K))/BLK12**2.0
      A5DG3B=DELTA(L,K)/BLK12-BL*(1.0+DELTA(L,K))/BLK12**2.0
  14  CONTINUE                                                         
      TKIN=TKIN-4*AL1*AL1*A6DG3*OVA12                                 
      TKIN=TKIN-4*Ak1*Ak1*A5DG3*OVB12                                 
      DTA1=DTA1
     1     - ( 8*AL1*DELTA(L,K)*A6DG3*OVA12 +4*AL1*AL1*A6DG3*DA12A1 )
     2     - ( 8*AK1*A5DG3*OVB12 +4*AK1*AK1*A5DG3*DB12A1 )
      DTA2=DTA2
     1     - ( 4*AL1*AL1*A6DG3*DA12A2 )
     2     - ( 4*AK1*AK1*A5DG3*DB12A2 )
      DTR1=DTR1
     1     - ( 4*AL1*AL1*A6DG3*DA12R1 )
     2     - ( 4*AK1*AK1*A5DG3*DB12R1 )
      DTR2=DTR2
     1     - ( 4*AL1*AL1*A6DG3*DA12R2 )
     2     - ( 4*AK1*AK1*A5DG3*DB12R2 )
      DTB =DTB
     1     - ( 4*AL1*AL1*A6DG3B*OVA12 +4*AL1*AL1*A6DG3*DA12B )
     2     - ( 4*AK1*AK1*A5DG3B*OVB12 +4*AK1*AK1*A5DG3*DB12B )

      TKIN=TKIN+(6*(AL1*A6DG3+AK1*A5DG3)+                              
     *     4*BL*BK*(R1LZ-R2LZ)**2.0)*OVLAP                 
      DTA1=DTA1
     1   + (6*(AL1*A6DG3+AK1*A5DG3)+ 4*BL*BK*(R1LZ-R2LZ)**2.0)*DSA1
     2   +  6*(DELTA(L,K)*A6DG3+A5DG3) *OVLAP
      DTA2=DTA2
     1   + (6*(AL1*A6DG3+AK1*A5DG3)+ 4*BL*BK*(R1LZ-R2LZ)**2.0)*DSA2
      DTR1=DTR1
     1   + (6*(AL1*A6DG3+AK1*A5DG3)+ 4*BL*BK*(R1LZ-R2LZ)**2.0)*DSR1
     2   + (4*BL*BK*2*(R1LZ-R2LZ)*DELTA(L,K))*OVLAP
      DTR2=DTR2
     1   + (6*(AL1*A6DG3+AK1*A5DG3)+ 4*BL*BK*(R1LZ-R2LZ)**2.0)*DSR2
     2   + ( -4*BL*BK*2*(R1LZ-R2LZ)*DELTA(L,K))*OVLAP
      DTB =DTB
     1   + (6*(AL1*A6DG3+AK1*A5DG3)+ 4*BL*BK*(R1LZ-R2LZ)**2.0)*DSB
     2   + ( 6*(AL1*A6DG3B+AK1*A5DG3B) 
     3      +4*(BL+DELTA(L,K)*BK)*(R1LZ-R2LZ)**2.0)*OVLAP
	WRITE(7,*) 'GRKIN.F DTA1,DTA2,DTB,DTR1,DTR2'  
	WRITE(6,*) 'GRKIN.F DTA1,DTA2,DTB,DTR1,DTR2'
      WRITE(7,*)  DTA1,DTA2,DTB,DTR1,DTR2
      WRITE(6,*)  DTA1,DTA2,DTB,DTR1,DTR2
	RETURN
	END
