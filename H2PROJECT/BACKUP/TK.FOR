      SUBROUTINE TK(OVLAP,AX,AY,AZ,A1,BX,BY,BZ,A2,                      
     *             CX,CY,CZ,A3,DX,DY,DZ,A4,A5,A6,TKIN)       
      IMPLICIT real*8 (a-h, o-z)               
C ***                                                               
C ***     FOUR CENTER KINETIC ENERGY INTEGRAL                       
C ***                                                               
C     COMMON/A2/ PI,PI2,PI3,PI32,PI52                                  
      DATA ONE,TWO,FOUR,SIX,EIGHT,HALF/1.,2.,4.,6.,8.,.5/              
C     write(6,900) OVLAP,AX,AY,AZ,A1,BX,BY,BZ,A2,                      
C    *             CX,CY,CZ,A3,DX,DY,DZ,A4,A5,A6       
      GAMA1=A1+A2                                                      
      GAMA2=A3+A4                                                      
      GAMA3=A5+A6                                                       
      PX=(A1*AX+A2*BX)/GAMA1                                           
      QX=(A3*CX+A4*DX)/GAMA2                                            
      PY=(A1*AY+A2*BY)/GAMA1                                            
      QY=(A3*CY+A4*DY)/GAMA2                                            
      PZ=(A1*AZ+A2*BZ)/GAMA1                                            
      QZ=(A3*CZ+A4*DZ)/GAMA2                                           
      QCX=QX-CX 
      QCY=QY-CY 
      QCZ=QZ-CZ                                 
      PAX=PX-AX 
      PAY=PY-AY 
      PAZ=PZ-AZ                                 
      PBX=PX-BX 
      PBY=PY-BY 
      PBZ=PZ-BZ                                 
      PQX=PX-QX 
      PQY=PY-QY 
      PQZ=PZ-QZ                                
      PQX2=PQX*PQX
      PQY2=PQY*PQY 
      PQZ2=PQZ*PQZ                      
      G1PG2=GAMA1+GAMA2                                                
      G1RG2=GAMA1*GAMA2                                                
      XMIAN=G1RG2+GAMA3*G1PG2                                          
      DELT=G1RG2*GAMA3/XMIAN*FOUR                                      
      GA1DEL=DELT/GAMA1                                                
      OVXA11=PAX-PQX*GA1DEL/FOUR                                       
      OVYA11=PAY-PQY*GA1DEL/FOUR                                       
      OVZA11=PAZ-PQZ*GA1DEL/FOUR                 
      OVXB11=PBX-PQX*GA1DEL/FOUR                 
      OVYB11=PBY-PQY*GA1DEL/FOUR                                 
      OVZB11=PBZ-PQZ*GA1DEL/FOUR                                     
      OVXA12=OVXA11*OVXA11                                            
      OVYA12=OVYA11*OVYA11                                             
      OVZA12=OVZA11*OVZA11                                             
      POM=HALF/GAMA1-GA1DEL/EIGHT/GAMA1                                
      OVXA12=(OVXA12+POM)*OVLAP                                        
      OVYA12=(OVYA12+POM)*OVLAP                                        
      OVZA12=(OVZA12+POM)*OVLAP                                        
      OVXB12=OVXB11*OVXB11                                             
      OVYB12=OVYB11*OVYB11                                             
      OVZB12=OVZB11*OVZB11                                             
      OVXB12=(OVXB12+POM)*OVLAP                                         
      OVYB12=(OVYB12+POM)*OVLAP                                         
      OVZB12=(OVZB12+POM)*OVLAP                                        
      GA2DEL=DELT/GAMA2                                                
      OVXC21=QCX+PQX*GA2DEL/FOUR                                       
      OVYC21=QCY+PQY*GA2DEL/FOUR                                        
      OVZC21=QCZ+PQZ*GA2DEL/FOUR                                        
      OVXC22=OVXC21*OVXC21*OVLAP                                        
      OVYC22=OVYC21*OVYC21*OVLAP                                        
      OVZC22=OVZC21*OVZC21*OVLAP                                        
      POM=HALF/GAMA2-GA2DEL/EIGHT/GAMA2                                 
      POM=POM*OVLAP                                                     
      OVXC22=OVXC22+POM                                                 
      OVYC22=OVYC22+POM                                                 
      OVZC22=OVZC22+POM                                                 
      POM=GA2DEL/EIGHT/GAMA1                                            
      XA1XC2=ONE-HALF*PQX2*DELT                                         
      YA1YC2=ONE-HALF*PQY2*DELT                                         
      ZA1ZC2=ONE-HALF*PQZ2*DELT                                         
      XA1XC2=POM*XA1XC2                                                 
      YA1YC2=POM*YA1YC2                                                 
      ZA1ZC2=POM*ZA1ZC2                                                 
      XA1XC2=PAX*QCX+PQX/FOUR*(PAX*GA2DEL-QCX*GA1DEL)+XA1XC2            
      YA1YC2=PAY*QCY+PQY/FOUR*(PAY*GA2DEL-QCY*GA1DEL)+YA1YC2            
      ZA1ZC2=PAZ*QCz+PQZ/FOUR*(PAZ*GA2DEL-QCZ*GA1DEL)+ZA1ZC2            
      XA1XC2=XA1XC2*OVLAP                                               
      YA1YC2=YA1YC2*OVLAP                                               
      ZA1ZC2=ZA1ZC2*OVLAP                                               
      ACX=AX-CX 
      ACY=AY-CY 
      ACZ=AZ-CZ                                 
      OVA12=OVXA12+OVYA12+OVZA12                                        
      OVA11=TWO*(OVXA11*ACX+OVYA11*ACY+OVZA11*ACZ)*OVLAP                
      OVC22=OVXC22+OVYC22+OVZC22                                       
      OVC21=-TWO*(OVXC21*ACX+OVYC21*ACY+OVZC21*ACZ)*OVLAP              
      A1C2=-TWO*(XA1XC2+YA1YC2+ZA1ZC2)                                 
      TKIN=FOUR*A5*A6*(OVA12+OVA11+OVC22+OVC21+A1C2)                   
      IF(GAMA3) 10,12,10                                               
  12  A6DG3=HALF                                                       
      A5DG3=HALF                                                       
      GO TO 14                                                         
  10  A6DG3=A6/GAMA3                                                   
      A5DG3=A5/GAMA3                                                   
  14  CONTINUE                                                         
      TKIN=TKIN-FOUR*A1*A1*A6DG3*OVA12                                 
      OVB12=OVXB12+OVYB12+OVZB12                                       
      TKIN=TKIN-FOUR*A2*A2*A5DG3*OVB12                                 
      TKIN=TKIN+(SIX*(A1*A6DG3+A2*A5DG3)+                              
     *     FOUR*A5*A6*(ACX*ACX+ACY*ACY+ACZ*ACZ))*OVLAP                 
      RETURN  
900   FORMAT(1x,F10.5)                                                         
      END                                                                   
