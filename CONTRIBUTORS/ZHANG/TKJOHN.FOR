

       SUBROUTINE TKjohn(OVLAP,R1LZ,AL1,R1KZ,AK1,                      
     *             R2LZ,AL2,R2KZ,AK2,BL,BK,TKIN)       
      IMPLICIT real*8 (a-h, o-z)               
C ***                                                               
C ***     FOUR CENTER KINETIC ENERGY INTEGRAL                       
C ***                                                               
C     COMMON/A2/ PI,PI2,PI3,PI32,PI52                                  
      DATA ONE,TWO,FOUR,SIX,EIGHT,HALF/1.,2.,4.,6.,8.,.5/              
      ALK1=AL1+AK1                                                      
      ALK2=AL2+AK2                                                      
      BLK12=BL+BK                                                      
      R1Z=(AL1*R1LZ+AK1*R1KZ)/ALK1                                            
      R2Z=(AL2*R2LZ+AK2*R2KZ)/ALK2                                     
      DETAB=ALK1*ALK2+BLK12*(ALK1+ALK2)                                   
      OVZA11=(R1Z-R1LZ)-(R1Z-R2Z)*ALK2*BLK12/DETAB                
      OVZB11=(R1Z-R1KZ)-(R1Z-R2Z)*ALK2*BLK12/DETAB   
      POM=0.5*(1.0/Alk1-ALK2*BLK12/(DETAB*ALK1) )                         
      OVZA12=(OVZA11*OVZA11+POM)*OVLAP                                        
      OVA12=OVZA12+2.0*POM*OVLAP
      OVA11=2*(OVZA11*(R1LZ-R2LZ))*OVLAP                
      OVZB12=(OVZB11*OVZB11+POM)*OVLAP                                       
      OVB12=OVZB12+2.0*POM*OVLAP 
      OVZC21=(R2Z-R2LZ)+(R1Z-R2Z)*ALK1*BLK12/DETAB             
      OVC21=-2*(OVZC21*(R1LZ-R2LZ))*OVLAP              
      POM=0.5*(1.0/ALK2-ALK1*BLK12/(DETAB*ALK2) )                       
      OVZC22=(OVZC21*OVZC21+POM)*OVLAP                       
      OVC22=OVZC22+2.0*POM*OVLAP
      XX1=(0.5*BLK12/DETAB)
     1  *(1.0-0.5*(R1Z-R2Z)**2.0*ALK1*ALK2*BLK12/DETAB*4.0)
      XX2=(R1Z-R1LZ)*(R2Z-R2LZ)+(R1Z-R2Z)*
     1  ((R1Z-R1LZ)*ALK1*BLK12-(R2Z-R2LZ)*ALK2*BLK12)/DETAB           
      za1zc2=(XX1+XX2)*OVLAP
c     WRITE(6,*) 'JOHN,ZA1ZC2'
c     WRITE(6,*)  ZA1ZC2
	A1C2=-2*(XX1+XX2+BLK12/DETAB)*OVLAP                                 

      TKIN=4*BL*BK*(OVA11 +OVA12+OVC22+OVC21+A1C2)                   
       write(6,*) "tk,john tkin=",tkin
      IF(BLK12) 10,12,10                                               
  12  A6DG3=0.5                                                     
      A5DG3=0.5                                                       
      GO TO 14                                                         
  10  A6DG3=BK/BLK12                                                   
      A5DG3=BL/BLK12                                                   
  14  CONTINUE                                                         
      TKIN=TKIN-4*AL1*AL1*A6DG3*OVA12                                 
      TKIN=TKIN-4*Ak1*Ak1*A5DG3*OVB12                                 
      TKIN=TKIN+(6*(AL1*A6DG3+AK1*A5DG3)+                              
     *     4*BL*BK*(R1LZ-R2LZ)**2.0)*OVLAP                 

      write(6,*) "tk,john tkin=",tkin
      RETURN  
900   FORMAT(1x,F10.5)                                                         
      END                                                                   
