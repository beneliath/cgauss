C     Originallly Coded:  Z. Zhang
C           Modified By:  D. Gilmore ---> (24 May 95)
C     -----------------------------------------------------
C     L. ADAMOWICZ RESEARCH / UNIVERSITY OF ARIZONA / CHEM.
C     -----------------------------------------------------
      SUBROUTINE HAMILTON(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,
     1                    S,H,CO)
C     -----------------------------------------------------
C     M = Number of Gaussian Geminals used
C     AA1 = Exponent Alpha _sub(1)
C     RX1 = x-Component of Position for Gaussian of AA1
C     RY1 = y-Component of Position for Gaussian of AA1  
C     RZ1 = z-Component of Position for Gaussian of AA1
C     AA2 = Exponent Alpha _sub(2)
C     RX2 = x-Component of Position for Gaussian of AA2
C     RY2 = y-Component of Position for Gaussian of AA2  
C     RZ2 = z-Component of Position for Gaussian of AA2
C     B12 = Correlation Exponent
C     S = S-Matrix ?
C     H = H-Matrix ?  
C     CO = ?
C     -----------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /VALUE/ EMIN
      COMMON /SMALLONE/ RJUICE
      PARAMETER(MMAX=50)
      DIMENSION 
     1         AA1(MMAX),RX1(MMAX),RY1(MMAX),RZ1(MMAX),
     2         AA2(MMAX),RX2(MMAX),RY2(MMAX),RZ2(MMAX),
     3         B12(MMAX)
C	DIMENSION DETS(50),CC(50)
C
C     ...This was left behind by Zhang... (???)
C     DIMENSION S(MMAX*(MMAX+1)/2),H(MMAX*(MMAX+1)/2) 
C
      DIMENSION S(MMAX,MMAX),H(MMAX,MMAX),CO(MMAX) 
      DATA ONE,TWO,FOUR,SIX,EIGHT,HALF/1.,2.,4.,6.,8.,.5/
C         ...DATA, Sets the DATA "..."/ ? , Variables equal
C            to the DATA ? / "..." Values. 
C              
C     -----------------------------------------------------
C     --- PRINT OUT GEMINAL DATA: ---
C     -----------------------------------------------------     
C      DO L=1,M
C          CALL GEMPRT(L,
C     1            AA1(L),RX1(L),RY1(L),RZ1(L),
C     2            AA2(L),RX2(L),RY2(L),RZ2(L),
C     3            B12(L))
C             ...GEMPRT, Prints out the data supplied in
C                File = GUESS.F
C
C      END DO
C
      CALL GEOM(EX1,EY1,EZ1,EX2,EY2,EZ2)
C     -----------------------------------------------------	  

      DO L=1,M
          DO K=1,M
C             .............................................
C
C      WRITE(6,*) 'Loop Index: ',K,' in ',L
C      WRITE(7,*) 'Loop Index: ',K,' in ',L
C     -----------------------------------------------------
C     --- Calculate OVERLAP INTEGRAL ---
C     -----------------------------------------------------
C
C     Symmetrization=(1+P_[elec.1,elec.2])(1+P_[H1,H2])
C                   =(1 + P_[elec.1,elec.2] + P_[H1,H2] + 
C                     P_[elec.1,elec.2]*P_[H1,H2])
C
C     ...Ordinary Overlap Integral: 1
C
              CALL OVERLP(RX1(L),RY1(L),RZ1(L),AA1(L),
     1            RX1(K),RY1(K),RZ1(K),AA1(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            RX2(K),RY2(K),RZ2(K),AA2(K),
     4            B12(L),B12(K),OV1,KLM)
C
C                        KLM = ?
C
C     ...Electron Symmetrization Overlap Integral: P_[elec.1,elec.2]
C  
              CALL OVERLP(RX1(L),RY1(L),RZ1(L),AA1(L),
     1            RX2(K),RY2(K),RZ2(K),AA2(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            RX1(K),RY1(K),RZ1(K),AA1(K),
     4            B12(L),B12(K),OV2,KLM)
C
C     ...Nuclear Symmetrization Overlap Integral: P_[H1,H2]

              CALL OVERLP(RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
     4            B12(L),B12(K),OV3,KLM)
C
C     ...Electro/Nuclear Symmetrization Overlap Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL OVERLP(RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
     4            B12(L),B12(K),OV4,KLM)
C              
              S(L,K)=OV1+OV2+OV3+OV4
C              S(L,K)=OV1+OV2
C
C                 ...S(L,K) = The Lth,Kth Element of the
C                             Overlap Matrix.
C
C     -----------------------------------------------------
C     --- Calculate KINETIC ENERGY OF ELECTRON 1 ---
C     -----------------------------------------------------
              CALL TK(OV1,
     1            RX1(L),RY1(L),RZ1(L),AA1(L),
     2            RX1(K),RY1(K),RZ1(K),AA1(K),                      
     3            RX2(L),RY2(L),RZ2(L),AA2(L),
     4            RX2(K),RY2(K),RZ2(K),AA2(K),
     5            B12(L),B12(K),T1)  
C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C  
              CALL TK(OV2,
     1            RX1(L),RY1(L),RZ1(L),AA1(L),
     2            RX2(K),RY2(K),RZ2(K),AA2(K),                      
     3            RX2(L),RY2(L),RZ2(L),AA2(L),
     4            RX1(K),RY1(K),RZ1(K),AA1(K),
     5            B12(L),B12(K),T2)  
C
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
C
              CALL TK(OV3,
     1            RX1(L),RY1(L),RZ1(L),AA1(L),
     2            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
     3            RX2(L),RY2(L),RZ2(L),AA2(L),
     4            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
     5            B12(L),B12(K),T3)
C
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL TK(OV4,
     1            RX1(L),RY1(L),RZ1(L),AA1(L),
     2            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
     3            RX2(L),RY2(L),RZ2(L),AA2(L),
     4            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
     5            B12(L),B12(K),T4)
C          
C     -----------------------------------------------------
C     --- Calculate KINETIC ENERGY OF ELECTRON 2 ---
C     -----------------------------------------------------
              CALL TK(OV1,
     1            RX2(L),RY2(L),RZ2(L),AA2(L),
     2            RX2(K),RY2(K),RZ2(K),AA2(K),                      
     3            RX1(L),RY1(L),RZ1(L),AA1(L),
     4            RX1(K),RY1(K),RZ1(K),AA1(K),
     5            B12(L),B12(K),T5)  
C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C  
              CALL TK(OV2,
     1            RX2(L),RY2(L),RZ2(L),AA2(L),
     2            RX1(K),RY1(K),RZ1(K),AA1(K),                      
     3            RX1(L),RY1(L),RZ1(L),AA1(L),
     4            RX2(K),RY2(K),RZ2(K),AA2(K),
     5            B12(L),B12(K),T6) 
C
C     ...Nuclear Symmetrization Integral: P_[H1,H2]

              CALL TK(OV3,
     1            RX2(L),RY2(L),RZ2(L),AA2(L),
     2            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
     3            RX1(L),RY1(L),RZ1(L),AA1(L),
     4            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
     5            B12(L),B12(K),T7)
C
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL TK(OV4,
     1            RX2(L),RY2(L),RZ2(L),AA2(L),
     2            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
     3            RX1(L),RY1(L),RZ1(L),AA1(L),
     4            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
     5            B12(L),B12(K),T8)
C  
C     -----------------------------------------------------
C     --- Calculate Nuc./Elec. Attractn. FOR ELECTRON 1 ---
C     -----------------------------------------------------
              CALL ENUCL(RX1(L),RY1(L),RZ1(L),AA1(L),
     1           RX1(K),RY1(K),RZ1(K),AA1(K),                      
     2           RX2(L),RY2(L),RZ2(L),AA2(L),
     3           RX2(K),RY2(K),RZ2(K),AA2(K),
     4           B12(L),B12(K),
     5           EX1,EY1,EZ1,EN1,KLM)
C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C  
              CALL ENUCL(RX1(L),RY1(L),RZ1(L),AA1(L),
     1           RX2(K),RY2(K),RZ2(K),AA2(K),                      
     2           RX2(L),RY2(L),RZ2(L),AA2(L),
     3           RX1(K),RY1(K),RZ1(K),AA1(K),
     4           B12(L),B12(K),
     5           EX1,EY1,EZ1,EN2,KLM)                       
C
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
C
              CALL ENUCL(RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
     4            B12(L),B12(K),
     5            EX1,EY1,EZ1,EN3,KLM)
C
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL ENUCL(RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
     4            B12(L),B12(K),
     5            EX1,EY1,EZ1,EN4,KLM)
C  
C
C
              CALL ENUCL(RX1(L),RY1(L),RZ1(L),AA1(L),
     1           RX1(K),RY1(K),RZ1(K),AA1(K),                      
     2           RX2(L),RY2(L),RZ2(L),AA2(L),
     3           RX2(K),RY2(K),RZ2(K),AA2(K),
     4           B12(L),B12(K),
     5           EX2,EY2,EZ2,EN5,KLM)
C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C
              CALL ENUCL(RX1(L),RY1(L),RZ1(L),AA1(L),
     1           RX2(K),RY2(K),RZ2(K),AA2(K),                      
     2           RX2(L),RY2(L),RZ2(L),AA2(L),
     3           RX1(K),RY1(K),RZ1(K),AA1(K),
     4           B12(L),B12(K),
     5           EX2,EY2,EZ2,EN6,KLM)                       
C
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
C
              CALL ENUCL(RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
     4            B12(L),B12(K),
     5            EX2,EY2,EZ2,EN7,KLM)
C
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL ENUCL(RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
     4            B12(L),B12(K),
     5            EX2,EY2,EZ2,EN8,KLM)
C 
C     -----------------------------------------------------
C     --- Calculate Nuc./Elec. Attractn. FOR ELECTRON 2 ---
C     -----------------------------------------------------
              CALL ENUCL(RX2(L),RY2(L),RZ2(L),AA2(L),
     1           RX2(K),RY2(K),RZ2(K),AA2(K),                      
     2           RX1(L),RY1(L),RZ1(L),AA1(L),
     3           RX1(K),RY1(K),RZ1(K),AA1(K),
     4           B12(L),B12(K),
     5           EX1,EY1,EZ1,EN9,KLM)                       
C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C
              CALL ENUCL(RX2(L),RY2(L),RZ2(L),AA2(L),
     1           RX1(K),RY1(K),RZ1(K),AA1(K),                      
     2           RX1(L),RY1(L),RZ1(L),AA1(L),
     3           RX2(K),RY2(K),RZ2(K),AA2(K),
     4           B12(L),B12(K),
     5           EX1,EY1,EZ1,EN10,KLM)                       
C
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
C
              CALL ENUCL(RX2(L),RY2(L),RZ2(L),AA2(L),
     1            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
     2            RX1(L),RY1(L),RZ1(L),AA1(L),
     3            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
     4            B12(L),B12(K),
     5            EX1,EY1,EZ1,EN11,KLM)
C
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL ENUCL(RX2(L),RY2(L),RZ2(L),AA2(L),
     1            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
     2            RX1(L),RY1(L),RZ1(L),AA1(L),
     3            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
     4            B12(L),B12(K),
     5            EX1,EY1,EZ1,EN12,KLM)
C 
C
C 
              CALL ENUCL(RX2(L),RY2(L),RZ2(L),AA2(L),
     1           RX2(K),RY2(K),RZ2(K),AA2(K),                      
     2           RX1(L),RY1(L),RZ1(L),AA1(L),
     3           RX1(K),RY1(K),RZ1(K),AA1(K),
     4           B12(L),B12(K),
     5           EX2,EY2,EZ2,EN13,KLM)                       
C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C
              CALL ENUCL(RX2(L),RY2(L),RZ2(L),AA2(L),
     1           RX1(K),RY1(K),RZ1(K),AA1(K),                      
     2           RX1(L),RY1(L),RZ1(L),AA1(L),
     3           RX2(K),RY2(K),RZ2(K),AA2(K),
     4           B12(L),B12(K),
     5           EX2,EY2,EZ2,EN14,KLM)  
C
C     ...Nuclear Symmetrization Integral: P_[H1,H2]
C
              CALL ENUCL(RX2(L),RY2(L),RZ2(L),AA2(L),
     1            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
     2            RX1(L),RY1(L),RZ1(L),AA1(L),
     3            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
     4            B12(L),B12(K),
     5            EX2,EY2,EZ2,EN15,KLM)
C
C     ...Electro/Nuclear Symmetrization Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL ENUCL(RX2(L),RY2(L),RZ2(L),AA2(L),
     1            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
     2            RX1(L),RY1(L),RZ1(L),AA1(L),
     3            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
     4            B12(L),B12(K),
     5            EX2,EY2,EZ2,EN16,KLM)
C 
C                     
C     -----------------------------------------------------
C     --- Calculate the Elec./Elec. Repulsion ---
C     -----------------------------------------------------
              CALL ELREP(RX1(L),RY1(L),RZ1(L),AA1(L),
     1           RX1(K),RY1(K),RZ1(K),AA1(K),                      
     2           RX2(L),RY2(L),RZ2(L),AA2(L),
     3           RX2(K),RY2(K),RZ2(K),AA2(K),
     4           B12(L),B12(K),ER1,KLM)  
C
C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
C
              CALL ELREP(RX1(L),RY1(L),RZ1(L),AA1(L),
     1           RX2(K),RY2(K),RZ2(K),AA2(K),                      
     2           RX2(L),RY2(L),RZ2(L),AA2(L),
     3           RX1(K),RY1(K),RZ1(K),AA1(K),
     4           B12(L),B12(K),ER2,KLM)  
C
C     ...Nuclear Symmetrization Overlap Integral: P_[H1,H2]
C
              CALL ELREP(RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
     4            B12(L),B12(K),ER3,KLM)
C
C     ...Electro/Nuclear Symmetrization Overlap Integral:
C         P_[elec.1,elec.2]*P_[H1,H2]
C
              CALL ELREP(RX1(L),RY1(L),RZ1(L),AA1(L),
     1            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
     2            RX2(L),RY2(L),RZ2(L),AA2(L),
     3            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
     4            B12(L),B12(K),ER4,KLM)
C              

C
              T=(T1+T2+T3+T4+T5+T6+T7+T8)/TWO
C              T=(T1+T2+T5+T6)/TWO
C              WRITE(6,7707) L,K,EN1
C              WRITE(6,7707) L,K,EN2
C              WRITE(6,7707) L,K,EN3
C              WRITE(6,7707) L,K,EN4
C              WRITE(6,7707) L,K,EN5
C              WRITE(6,7707) L,K,EN6
C              WRITE(6,7707) L,K,EN7
C              WRITE(6,7707) L,K,EN8
C
C              WRITE(7,7707) L,K,EN1
C              WRITE(7,7707) L,K,EN2
C              WRITE(7,7707) L,K,EN3
C              WRITE(7,7707) L,K,EN4
C              WRITE(7,7707) L,K,EN5
C              WRITE(7,7707) L,K,EN6
C              WRITE(7,7707) L,K,EN7
C              WRITE(7,7707) L,K,EN8
C
C
C 7707 format(1x,'matrix element: l = ',i3,' k = ',i3,
C     1    '    en = 'f15.9)
C
C                 ...Total Kinetic Energy is equal to half
C                    of the sums of individual K.E. calcs.
C                    because of double counts in the nested
C                    loop.
C
              EN=EN1+EN2+EN3+EN4+EN5+EN6+EN7+EN8
     1            +EN9+EN10+EN11+EN12+EN13+EN14+EN15+EN16
C
C                 ...Total ElectroNuclear Attraction Energy
C                    is equal to the sum of the individual
C                    ElectroNuclear calcs.  
C

              ER=ER1+ER2+ER3+ER4
C
C                 ...Total Elec./Elec. Repulsion Energy is 
C                    equal to the sum of the individual
C                    Elec./Elec. calcs. 
C
              H(L,K)=ER+T-EN+(1.0D+00/1.4D+00)*S(L,K)
C                 ...Hamiltonian Matrix Element (L,K) equals
C                    Electron Repulsion + Kinetic Energy -
C                    ElectroNuclear Attraction.
C
C             .............................................
          END DO 
      END DO
C     ----------------------------------------------------- 
C      CALL TAB(S,M,M,MMAX,MMAX)                  
C         ...Prints out the Overlap Matrix / "S-Matrix"
C      CALL TAB(H,M,M,MMAX,MMAX)
C         ...Prints out the Hamiltonian Matrix / "H-Matrix"
      CALL ENERGY(M,H,S,EMIN,CO,KIM)
C         ...Calculates the Ground State Energy & Prints It
C     -----------------------------------------------------      
      JUICE=EMIN
      RETURN
      END
