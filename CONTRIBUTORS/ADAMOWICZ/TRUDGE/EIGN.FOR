      SUBROUTINE EIGN(N,EIG)                                            AA     1
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 IPOSV,IVPOS,IORD
C                  ****  VERSION  3/20/73  ****                         A      3
C        *****  DESCRIPTION                               *****         A      4
C        *****   THIS IS A GIVENS-HOUSHOLDER ROUTINE      *****         A      5
C        *****   RECEIVES A MATRIX-A-OF WHICH ONLY LOWER  *****         A      6
C        *****   TRANGLE IS USED AND THIS IS DESTROYED    *****         A      7
C        *****   EIGENVALUES ARE RETURNED IN ALGEBRAIC    *****         A      8
C        *****   ASCENDING ORDER IN ARRAY-EIG-            *****         A      9
C        *****   THE EIGENVECTORS ARE RETURNED STORED BY  *****         A     10
C        *****   COLUMNS IN VEC.                          *****         A     11
C        *****         PARAMETERS PASSED                  *****         A     12
C        *****   RHO=TOLERANCE ON FINAL OFF-DIAG ELEMENTS *****         A     13
C        *****   N  IS THE SIZE OF THE MATRIX             *****         A     14
C        *****   OF THAT IN QCPE PROGRAM CINDO BY DOBOSH  *****         A     15
C        *****  THIS ROUTINE IS SLIGHTLY MODIFIED VERSION *****         A     16
C                                                                       A     17
      DIMENSION W(20),EIG(20),GAMMA(20),BETA(20),P(20),Q(20),BETASQ(20)       AA     2
      COMMON /EIGMAT/A(20,20),VEC(20,20),IPOSV(20),IVPOS(20),IORD(20)   AA     3
      EQUIVALENCE (IPOSV,GAMMA),(IVPOS,BETA),(IORD,BETASQ),             AA     4
     *  (IVPOS,P),(IVPOS,Q)                                             AA     5
      RHO=1.0D-10                                                       AA     6
      RHOSQ=RHO*RHO                                                     A     23
      IF(N.EQ.0) GO TO 560                                              A     24
      SHIFT = 0.D+00                                                       A
      N1=N-1                                                            A     26
      N2=N-2                                                            A     27
      GAMMA(1)=A(1,1)                                                   A     28
      IF(N2.EQ.0) GO TO 270                                             A     29
      IF(N2.LT.0)GO TO 280                                              A     30
      DO 260 NR=1,N2                                                    A     31
      B=A(NR+1,NR)                                                      A     32
      S=0.D+00                                                              A
      DO 130 I=NR,N2                                                    A     34
  130 S=S+A(I+2,NR)**2                                                  A     35
C        *****  PREPARE FOR POSSIBLE BYPASS OF TRANSFORMATION  *****    A     36
      A(NR+1,NR)=0.D+00                                                     A
      IF(S.LE.0) GO TO 250                                              A     38
      S=S+B*B                                                           A     39
      SGN=+1.D+00                                                           A
      IF(B.GE.0) GO TO 160                                              A     41
      SGN = -1.D+00                                                        A
160   SQRTS=DSQRT(S)                                                     AA
      D=SGN/(SQRTS+SQRTS)                                               A     44
      TEMP=DSQRT(5.D-01+B*D)                                                 AA
      W(NR)=TEMP                                                        A     46
      A(NR+1,NR)=TEMP                                                   A     47
      D=D/TEMP                                                          A     48
      B=-SGN*SQRTS                                                      A     49
C        *****  D IS FACTOR OF PROPORTIONALITY.           *****         A     50
C        *****  COMPUTE AND SAVE W VECTOR                               A     51
      DO 170 I=NR,N2                                                    A     52
      TEMP=D*A(I+2,NR)                                                  A     53
      W(I+1)=TEMP                                                       A     54
  170 A(I+2,NR)=TEMP                                                    A     55
C        *****  PREMULTIPLY VECTOR-W-BY MATRIX -A- TO OBTAIN *****      A     56
C        *****  -P- VECTOR.  SIMULTANEOUSLY ACCUMULATE DOT   *****      A     57
C        *****  PRODUCT WP,(THE SCALAR K).                  *****       A     58
      WTAW = 0.D+00                                                        A
      DO 220 I=NR,N1                                                    A     60
      SUM = 0.D+00                                                         A
      DO 180 J=NR,I                                                     A     62
  180 SUM=SUM+A(I+1,J+1)*W(J)                                           A     63
      I1=I+1                                                            A     64
      IF((N1-I1).LT.0) GO TO 210                                        A     65
      DO 200 J=I1,N1                                                    A     66
  200 SUM=SUM+A(J+1,I+1)*W(J)                                           A     67
  210 P(I)=SUM                                                          A     68
  220 WTAW=WTAW+SUM*W(I)                                                A     69
C        *****  P VECTOR AND SCALAR K NOW STORED.NEXT     *****         A     70
C        ***** COMPUTE -Q- VECTOR                         *****         A     71
      DO 230 I=NR,N1                                                    A     72
  230 Q(I)=P(I)-WTAW*W(I)                                               A     73
C        *****  NOW FORM PAP MATRIX, REQUIRED PART        *****         A     74
      DO 240 J=NR,N1                                                    A     75
      QJ=Q(J)                                                           A     76
      WJ=W(J)                                                            A     7
      DO 240 I=J,N1                                                     A     78
  240 A(I+1,J+1)=A(I+1,J+1)-2.D+00*(W(I)*QJ+WJ*Q(I))                        A
  250 BETA(NR)=B                                                        A     80
      BETASQ(NR)=B*B                                                    A     81
  260 GAMMA(NR+1)=A(NR+1,NR+1)                                          A     82
  270 B=A(N,N-1)                                                        A     83
      BETA(N-1)=B                                                       A     84
      BETASQ(N-1)=B*B                                                   A     85
      GAMMA(N)=A(N,N)                                                   A     86
  280 BETASQ(N)=0.D+00                                                      A
C        *****  ADJOIN AN IDENTITY MATRIX TO BE           *****         A     88
C        *****  POSTMULTIPLIED BY ROTATIONS               *****         A     89
      DO 300 I=1,N                                                      A     90
      DO 290 J=1,N                                                      A     91
  290 VEC(I,J)=0.D+00                                                        A
  300 VEC(I,I)=1.D+00                                                        A
      M=N                                                                A     9
      SUM=0.D+00                                                            A
      NPAS=1                                                            A     96
      GO TO 400                                                         A     97
  310 SUM=SUM+SHIFT                                                     A     98
      COSA=1.D+00                                                           A
      G=GAMMA(1)-SHIFT                                                  A    100
      PP=G                                                               A    10
      PPBS=PP*PP+BETASQ(1)                                              A    102
      PPBR=DSQRT(PPBS)                                                   AA
      DO 370 J=1,M                                                      A    104
      COSAP=COSA                                                        A    105
      IF(PPBS.NE.0) GO TO 320                                           A    106
      SINA = 0.D+00                                                        A
      SINA2=0.D+00                                                          A
      COSA=1.D+00                                                           A
      GO TO 350                                                         A    110
  320 SINA=BETA(J)/PPBR                                                 A    111
      SINA2=BETASQ(J)/PPBS                                              A    112
      COSA=PP/PPBR                                                      A    113
C        *****  POSTMULTIPLY IDENTITY BY P-TRANSPOSE      *****         A    114
      NT=J+NPAS                                                         A    115
      IF(NT.LE.N) GO TO 330                                             A    116
      NT = N                                                            A    117
  330 DO 340 I=1,NT                                                     A    118
      TEMP=COSA*VEC(I,J)+SINA*VEC(I,J+1)                                A    119
      VEC(I,J+1)=-SINA*VEC(I,J)+COSA*VEC(I,J+1)                         A    120
  340 VEC(I,J)=TEMP                                                     A    121
 350  DIA=GAMMA(J+1)-SHIFT                                              A    122
      U=SINA2*(G+DIA)                                                   A    123
      GAMMA(J)=G+U                                                      A    124
      G=DIA-U                                                           A    125
      PP=DIA*COSA-SINA*COSAP*BETA(J)                                    A    126
      IF(J.NE.M) GO TO 360                                              A    127
      BETA(J)=SINA*PP                                                   A    128
      BETASQ(J)=SINA2*PP*PP                                             A    129
      GO TO 380                                                         A    130
  360 PPBS=PP*PP+BETASQ(J+1)                                            A    131
      PPBR=DSQRT(PPBS)                                                   AA    1
      BETA(J)=SINA*PPBR                                                 A    133
  370 BETASQ(J)=SINA2*PPBS                                              A    134
 380  GAMMA(M+1)=G                                                      A    135
C        *****  TEST FOR CONVERGENCE OF LAST DIAGONAL ELEMENT  *****    A    136
      NPAS=NPAS+1                                                       A    137
      IF(BETASQ(M).GT.RHOSQ) GO TO 410                                  A    138
  390 EIG(M+1)=GAMMA(M+1)+SUM                                           A    139
  400 BETA(M)=0.D+00                                                        A
      BETASQ(M)=0.D+00                                                      A
      M=M-1                                                             A    142
      IF(M.EQ.0) GO TO 430                                              A    143
      IF(BETASQ(M).LE.RHOSQ) GO TO 390                                  A    144
C        *****  TAKE ROOT OF CORNER 2 BY 2 NEAREST TO     *****         A    145
C        *****  LOEWR DIAGONAL IN VALUE AS ESTIMATE OF    *****         A    146
C        *****  EIGENVALUE TO USE FOR SHIFT               *****         A    147
  410 A2=GAMMA(M+1)                                                     A    148
      R2=5.D-01*A2                                                          A
      R1=5.D-01*GAMMA(M)                                                    A
      R12=R1+R2                                                         A    151
      DIF=R1-R2                                                         A    152
      TEMP=DSQRT(DIF*DIF+BETASQ(M))                                      AA    1
      R1=R12+TEMP                                                       A    154
      R2=R12-TEMP                                                       A    155
      DIF=DABS(A2-R1)-DABS(A2-R2)                                         AA
      IF(DIF.LT.0) GO TO 420                                            A    157
      SHIFT=R2                                                          A    158
      GO TO 310                                                         A    159
  420 SHIFT=R1                                                          A    160
      GO TO 310                                                         A    161
  430 EIG(1)=GAMMA(1)+SUM                                               A    162
C        *****  INITIALIZE AUXILIARY TABLES REQUIRED FOR  *****         A    163
C        *****  REARANGING THE VECTORS                    *****         A    164
      DO 440 J=1,N                                                      A    165
      IPOSV(J)=J                                                        A    166
      IVPOS(J)=J                                                        A    167
  440 IORD(J) = J                                                       A    168
C        ***** USE A TRANSPOSITION SORT TO ORDER THE EIGENVALUES  ***** A    169
      M=N                                                                A    17
      GO TO 470                                                         A    171
  450 DO 460 J=1,M                                                      A    172
      IF(EIG(J).LE.EIG(J+1)) GO TO 460                                  A    173
      TEMP=EIG(J)                                                       A    174
      EIG(J)=EIG(J+1)                                                   A    175
      EIG(J+1)=TEMP                                                     A    176
      ITEMP=IORD(J)                                                     A    177
      IORD(J)=IORD(J+1)                                                 A    178
      IORD(J+1)=ITEMP                                                   A    179
  460 CONTINUE                                                          A    180
  470 M=M-1                                                             A    181
      IF(M.NE.0) GO TO 450                                              A    182
      IF(N1.EQ.0) GO TO 500                                             A    183
      DO 490 L=1,N1                                                     A    184
      NV=IORD(L)                                                        A    185
      NP=IPOSV(NV)                                                      A    186
      IF(NP.EQ.L) GO TO 490                                             A    187
      LV=IVPOS(L)                                                       A    188
      IVPOS(NP)=LV                                                      A    189
      IPOSV(LV)=NP                                                      A    190
      DO 480 I=1,N                                                      A    191
      TEMP=VEC(I,L)                                                     A    192
      VEC(I,L)=VEC(I,NP)                                                A    193
  480 VEC(I,NP) = TEMP                                                  A    194
  490 CONTINUE                                                          A    195
  500 CONTINUE                                                          A    196
C        *****  BACK TRANSFORM THE VECTORS OF THE TRIPLE  *****         A    197
C        *****  DIAGONAL MATRIX                                         A    198
      DO 550 NRR=1,N                                                    A    199
      K=N1                                                              A    200
  510 K=K-1                                                             A    201
      IF(K.LE.0) GO TO 540                                              A    202
      SUM = 0.D+00                                                         A
      DO 520 I=K,N1                                                     A    204
  520 SUM=SUM+VEC(I+1,NRR)*A(I+1,K)                                     A    205
      SUM=SUM+SUM                                                       A    206
      DO 530 I=K,N1                                                     A    207
  530 VEC(I+1,NRR)=VEC(I+1,NRR)-SUM*A(I+1,K)                            A    208
      GO TO 510                                                         A    209
  540 CONTINUE                                                          A    210
  550 CONTINUE                                                          A    211
  560 CONTINUE                                                          A    212
      RETURN                                                            A    213
      END                                                               A    214
