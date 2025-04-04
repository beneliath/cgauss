      SUBROUTINE PLODDING(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,S,H)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(MMAX=50)
      COMMON /T1/ ND,NMAX,FFF,RR(20),DR(20),CURVE,ALPHA,NOISE,EXACT,
     *           ALPH(20),RNOT(20),V(20,20),CURV(20),RESTAR,JSTART,
     *           FNOT,TIMLIM,TOLF,TOLR
      COMMON /VALUE/ EMIN
      DIMENSION 
     1         AA1(MMAX),RX1(MMAX),RY1(MMAX),RZ1(MMAX),
     2         AA2(MMAX),RX2(MMAX),RY2(MMAX),RZ2(MMAX),
     3         B12(MMAX)
      EXTERNAL FU
      LOGICAL EXACT,RESTAR
      REAL*8 NOISE
C
      TOLF = 1.D-06
      TOLR = 1.D-06
      NOISE = 1.D-08
      ND = 6
      NMAX = 20
      KKKK=-1 
      RR(1) = AA1(1)
      RR(2) = B12(1)
      RR(3) = AA1(2)
      RR(4) = B12(2)
      RR(5) = AA1(3)
      RR(6) = B12(3)
C      FFF=1.0D+03
      FFF = EMIN
      CALL TRUDGE(KKKK,FU,M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,S,H)
      STOP
      END
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

      SUBROUTINE FU(I,M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,S,H)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(MMAX=50)
c ***
      COMMON/T1/ N,NMAX,F,R(500),DR(500),CURVE,ALPHA,NOISE,EXACT,ALPH(50
	10),RNOT(500), V(500,500),CURV(500),RESTART,JSTART,FNOT,TIMLIM,                   
     2 TOLF,TOLR                                                        
      COMMON /VALUE/ EMIN
      DIMENSION 
     1         AA1(MMAX),RX1(MMAX),RY1(MMAX),RZ1(MMAX),
     2         AA2(MMAX),RX2(MMAX),RY2(MMAX),RZ2(MMAX),
     3         B12(MMAX)
C ***

      I = 0
      AA1(1) = R(1)
      AA2(1) = R(1)
      B12(1) =R(2)
      AA1(2) = R(3)
      AA2(2) = R(3)
      B12(2) =R(4)
      AA1(3) = R(5)
      AA2(3) = R(5)
      B12(3) =R(6)
C      fu = 
C      f = (x-4.d+00)*(x-4.d+00) + (y-10.d+00)*(y-10.d+00)
      CALL HAMILTON(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,S,H,CO)
      F=EMIN
      WRITE (6,*) ' AA1,AA2,B12 (1) = ',AA1(1),AA2(1),B12(1)
      WRITE (6,*) ' AA1,AA2,B12 (2) = ',AA1(2),AA2(2),B12(2)
      WRITE (6,*) ' AA1,AA2,B12 (3) = ',AA1(3),AA2(3),B12(3)
      WRITE (6,*) '            EMIN = ',F,'/'
      WRITE (7,*) ' AA1,AA2,B12 (1) = ',AA1(1),AA2(1),B12(1)
      WRITE (7,*) ' AA1,AA2,B12 (2) = ',AA1(2),AA2(2),B12(2)
      WRITE (7,*) ' AA1,AA2,B12 (3) = ',AA1(3),AA2(3),B12(3)
      WRITE (7,*) '            EMIN = ',F,'/'

      RETURN
      END   
      SUBROUTINE SEARCH(LEAVE,GAUSS,M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,
     1    B12,S,H)                                  
      IMPLICIT REAL*8(A-H,O-Z)
C                                                                        A    50
C                                                                       A    507
C     SUBROUTINE LOCATES A LOCAL MINIMUM ALONG A GIVEN ONE-DIMENSIONAL  A    508
C     SEARCH DIRECTION.                                                 A    509
C     F(ALF) = FUNC(X+ALF*DX) FIND OPTIMUM ALF VALUE.                   A    510
C     LOOKMAX = MAXIMUM NUMBER OF FUNCTION EVALUATIONS BEFORE SEARCH IS A    511
C               TERMINATED                                              A    512
C     ENTER SUBROUTINE WITH FOLLOWING                                   A    513
C     N     = DIMENSIONALITY OF X SPACE                                 A    514
C     X     = INITIAL POSITION (IN N DIMENSIONAL PARAMETER SPACE)       A    515
C     FUNC  = F(0) , I.E. FUNCTION VALUE AT INITIAL X.                  A    516
C     DX    = SEARCH DIRECTIONS(IN N-DIMENSIONAL SPACE)                 A    517
C     AMIN  = ESTIMATED INITIAL STEP SIZE                               A    518
C     NOISE = ACCURACY OF FUNC VALUES. VARIATIONS SMALLER THAN NOISE    A    519
C             ARE NOT CONSIDERED TO BE SIGNIFICANT.                     A    520
C     EXACT = .FALSE.  SUPPRESSES FINAL CHECK OF QUADRATIC FIT.         A    521
C                                                                        A    52
C     EXIT FROM SUBROUTINE WITH FOLLOWING                               A    523
C     X     = CALCULATED MINIMUM                                        A    524
C     FUNC  = F(AMIN) ,I.E. FUNCTION VALUE AT FINAL X.                  A    525
C     CURVE = CALCULATED CURVATURE ALON SEARCH DIRECTIONS               A    526
C     LEAVE = 0 IF NORMAL EXIT                                          A    527
C           = 2 IF EXCESSIV NUMBER OF FUNCTIONS EVALUATIONS             A    528
C                                                                       A    529
C     N,NMAX,FUNC,AND X ARE SHARED WITH SUBROUTINE GAUSS.               A    530
C     ALL VARIABLES IN COMMON T1 ARE SHARED WITH PROGRAM CALLING SEARCH A    531
C                                                                       A    532
      COMMON/T1/ N,NMAX, FUNC,X(20), DX(20), CURVE,AMIN ,NOISE,EXACT    A    533
      DIMENSION XNOT(20), F(9), SET(9)                                  A    534
      EQUIVALENCE (F(1),F1),(F(3),F3),(F(4),F4),(F(5),F5),(F(6),F6),    A    535
     1  (F(7),F7),(F(9),F9),(SET(1),SET1),                              A    536
     2  (SET(2),SET2),(SET(3),SET3),(SET(4),SET4),(SET(5),SET5),        A    537
     3  (SET(6),SET6),(SET(7),SET7),(SET(8),SET8),(SET(9),SET9)         A    538
      LOGICAL BACK, SET, SET1,SET2,SET3,SET4,SET5,SET6,SET7,SET8,SET9   A    539
     1,EXACT,HOLD                                                       A    540
      LOGICAL AIF
      REAL*8 NOISE                                                      A    541
      LOOKS=0                                                           A    542
      LEAVE=0                                                           A    543
      LOOKMAX=100                                                        A    544
      DO 5 I=1,N                                                        A    545
    5 XNOT(I)=X(I)                                                      A    546
      ALF=AMIN                                                          A    547
      AMIN=0.D+00                                                          A
      FNOT = FUNC                                                       A    549
      FMIN=FUNC                                                         A    550
      FNEXT=1.0D+30                                                     A    551
      IR=4                                                              A    552
   10 LOOKS=LOOKS+1                                                     A    553
      IF(LOOKS.LE.LOOKMAX) GO TO 20                                     A    554
      LEAVE=2                                                           A    555
      WRITE(6,13) LOOKS                                                 A    556
   13 FORMAT(/' EXCESSIVE NUMBER OF FUNCTION EVALUATIONS =',I3/)        A    557
C                                                                       A    558
C     RETURN FROM SUBROUTINE                                            A    559
C                                                                       A    560
   15 DO 16 I=1,N                                                       A    561
   16 X(I)=XNOT(I)+AMIN*DX(I)                                           A    562
      FUNC=FMIN                                                         A    563
      RETURN                                                            A    564
C                                                                        A    56
C     EVALUATE F(ALF)=FUNC(XNOT+ALF*DX)                                 A    566
C                                                                        A    56
   20 DO 25 I=1,N                                                       A    568
   25 X(I)=XNOT(I)+ALF*DX(I)                                            A    569
      INOUT=0                                                           A    570
      CALL GAUSS(INOUT,M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,S,H)                                                 A    571
C                                                                        A    57
C     INOUT = 7 INDICATES PATHOLOGICAL X VALUE FOR GAUSS CALCULATION.   A    573
C                                                                       A    574
      IF(INOUT.EQ.7) GO TO 305                                          A    575
      HOLD=.FALSE.                                                      A    576
      IF(FUNC.GE.FMIN) GO TO 28                                         A    577
      FNEXT=FMIN                                                        A    578
      FMIN=FUNC                                                         A    579
      ANEXT=AMIN                                                        A    580
      AMIN=ALF                                                          A    581
      GO TO 29                                                          A    582
   28 IF(FUNC.GE.FNEXT) GO TO 29                                        A    583
      FNEXT=FUNC                                                        A    584
      ANEXT=ALF                                                          A    58
   27 FORMAT(' LOOKS =',I3,' ALF =',F10.5, ' FUNC =',F15.9, ' AMIN =',  A    586
     1 F10.5,' FMIN =',F12.6,' ANEXT =',F10.5,' FNEXT =',F12.6)         A    587
   29 WRITE(6,27) LOOKS, ALF, FUNC, AMIN, FMIN, ANEXT, FNEXT
      WRITE(7,27) LOOKS, ALF, FUNC, AMIN, FMIN, ANEXT, FNEXT
      GO TO (29,29,29,150,180,260,320) IR                               A    589
C                                                                       A    590
C     F1,F3,F5,F7,F9, ARE FUNC VALUES AT EVENTLY SPACED POINTS WITH     A    591
C     SEPARATION DELTA. SET(I)=TRUE INDICATES THAT F(I) HAS BEEN        A    592
C     EVALUATED. DELTA MAY INCREASE OR DECREASE BY FACTORS OF 2. 0.     A    593
C     F VALUES SHIFT LEFT OR RIGHT SO THAT F5 IS ALWAYS LOWEST.         A    594
C                                                                        A    59
  150 DELTA=DABS(AMIN-ANEXT)                                              A    5
      BACK=.FALSE.                                                      A    597
      DO 155 I=1,9                                                       A    59
  155 SET(I)=.FALSE.                                                     A    59
      F5=FMIN                                                           A    600
      SET5=.TRUE.                                                       A    601
      K=3                                                                A    60
      IF(AMIN.LT.ANEXT) K=7                                             A    603
      F(K)=FNEXT                                                        A    604
      SET(K)=.TRUE.                                                     A    605
  160 IF(SET7) GO TO 230                                                A    606
      IF(SET3) GO TO 200                                                A    607
      IF(BACK) GO TO 250                                                A    608
  170 ALF=AMIN+DELTA                                                    A    609
      IR=5                                                               A    61
      GO TO 10                                                          A    611
  180 F7=FUNC                                                           A    612
      SET7=.TRUE.                                                       A    613
      IF(ALF.NE.AMIN) GO TO 160                                         A    614
C                                                                       A    615
C     SHIFT LEFT                                                        A    616
C                                                                       A    617
      DO 190 I=1,7                                                      A    618
      SET(I)=SET(I+2)                                                   A    619
  190 F(I)=F(I+2)                                                       A    620
      SET8=.FALSE.                                                      A    621
      SET9=.FALSE.                                                      A    622
      GO TO 160                                                         A    623
  200 CONTINUE
      IF(SET1) GOTO 210
      IF(.NOT.SET1) GO TO 170
C                                                                        A    62
C     EXTRAPOLATE USING F1,F3, AND F5                                   A    626
C     ICREASE STEP SIZE IF F9 IS PREDICTED TO BE LOWER THAN F7.         A    627
C     DECREASE STEP SIZE IF F4 IS PREDICTED TO BE LOWER THAN F5.        A    628
C                                                                       A    629
  210 BACK=4.D+00*F3.LT.(F1+3.D+00*F5)                                        A
      IF(SET9.OR.HOLD) GO TO 215                                        A    631
      IF(5.D+00*F3.GT.(2.D+00*F1+3.D+00*F5)) GO TO 220
  215 IF(SET4) GO TO 170                                                A    633
      AIF=6.D+00*F3.GT.(F1+5.D+00*F5)                                     634
      IF(AIF) GO TO 170
      IF(.NOT.AIF) GO TO 305
  220 DELTA=DELTA+DELTA                                                 A    635
      F4=F3                                                             A    636
      F3=F1                                                             A    637
      F6=F7                                                             A    638
      F7=F9                                                             A    639
      SET4=SET3                                                         A    640
      SET3=SET1                                                         A    641
      SET6=SET7                                                         A    642
      SET7=SET9                                                         A    643
      SET1=.FALSE.                                                      A    644
      SET2=.FALSE.                                                      A    645
      SET8=.FALSE.                                                      A    646
      SET9=.FALSE.                                                      A    647
      GO TO 160                                                         A    648
  230 IF(SET3) GO TO 300                                                A    649
      IF(SET9) GO TO 240
      IF(.NOT.SET9) GO TO 250
C                                                                       A    651
C     EXTRAPOLATE USING F5, F7, AND F9.                                 A    652
C     INCREASE ISTEP SIZE IF F1 IS PREDICTED TO B LOWER THAN F3         A    653
C     DECREASE STEP SIZE IF F6 IS PREDICTED TO BE LOWER THAN F5.        A    654
C                                                                        A    65
  240 BACK=4.D+00*F7.GT.(F9+3.D+00*F5)                                        A
      IF(SET1.OR.HOLD) GO TO 245                                        A    657
      IF(5.D+00*F7.GT.(2.D+00*F9+3.D+00*F5)) GO TO 220
  245 IF(SET6) GO TO 250                                                A    659
      AIF=6.D+00*F7.GT.(F9+5.D+00*F5)                                    660
      IF(AIF) GO TO 250
      IF(.NOT.AIF) GO TO 305
  250 ALF=AMIN-DELTA                                                    A    661
      IR=6                                                              A    662
      GO TO 10                                                          A    663
  260 F3=FUNC                                                           A    664
      SET3=.TRUE.                                                       A    665
      IF(ALF.NE.AMIN) GO TO 160                                         A    666
      DO 270 I=1,7                                                      A    667
      SET(10-I)=SET(8-I)                                                A    668
  270 F(10-I)=F(8-I)                                                    A    669
      SET1=.FALSE.                                                      A    670
      SET2=.FALSE.                                                      A    671
      GO TO 160                                                         A    672
C                                                                        A    67
C     INTERPOLATE USING F3, F5, AND F7.                                 A    674
C     EITHER IMMEDIATELY ATTEMPT QUADRATIC FIT OR FIRST REDUCE STEP SIZEA    675
C     IF AMIN=ZERO STEP SIZE  MAY BE REDUCED BY A FACTOR OF 2,4,8,OR 16 A    676
C     IN OTHER CASES STEP SIZE IS REDUCED BY A FACTOR OF 2.             A    677
C                                                                       A    678
  300 DUM=F3                                                            A    679
      IF(F7.GT.F3) DUM=F7                                               A    680
      SPREAD=DUM-F5                                                     A    681
      TOL=(FNOT-FMIN)*5.D-01                                               A
      IF(TOL.LT.1.D+03*NOISE) TOL=1.D+03*NOISE                          A    683
      IF(.NOT.EXACT) TOL=TOL+TOL                                        A    684
      DUM=F3-F5-F5+F7                                                   A    685
      CURVE=DUM/(2.D+00*DELTA*DELTA)                                       A
      DUM=(F3-F7)/DUM                                                   A    687
      BACK=DUM.LT.0.D+00                                                   A
      IF(SPREAD.LT.TOL) GO TO 310                                       A    689
      IF(DABS(DUM).GT.75.D-02.OR.AMIN.NE.0.D+00) GO TO 305                     A
      SET1=SET4                                                          A    69
      F1=F4                                                             A    692
      SET9=SET6                                                         A    693
      F9=F6                                                             A    694
      SET2=.FALSE.                                                       A    69
      SET3=.FALSE.                                                      A    696
      SET4=.FALSE.                                                      A    697
      SET6=.FALSE.                                                      A    698
      SET7=.FALSE.                                                      A    699
      SET8=.FALSE.                                                      A    700
      DELTA=25.D-02*DELTA                                                  A
      IF(DABS(DUM).GT.375.D-03) GO TO 160                                   A
      SET1=.FALSE.                                                       A    70
      SET9=.FALSE.                                                      A    704
      DELTA=5.D-1*DELTA                                                   A    7
      IF(DABS(DUM).GT.1875.D-04) GO TO 160                                  A
      DELTA=5.D-1*DELTA                                                   A    7
      GO TO 160                                                         A    708
C                                                                       A    709
C     CUT STEP SIZE IN HALF                                             A    710
C                                                                        A    71
  305 DELTA=5.D-01*DELTA                                                   A
      HOLD=.TRUE.                                                       A    713
      SET1=SET3                                                          A    71
      SET3=SET4                                                         A    715
      SET9=SET7                                                          A    71
      SET7=SET6                                                         A    717
      SET2=.FALSE.                                                      A    718
      SET4=.FALSE.                                                      A    719
      SET6=.FALSE.                                                      A    720
      SET8=.FALSE.                                                      A    721
      F1=F3                                                             A    722
      F3=F4                                                             A    723
      F9=F7                                                             A    724
      F7=F6                                                             A    725
      GO TO 160                                                         A    726
C                                                                       A    727
C     ATTEMPT QUADRATIC FIT                                             A    728
C                                                                        A    72
  310 ALF=AMIN+DUM*DELTA*5.D-1                                            A    7
      FQUAD=F5-(F3-F7)*DUM*125.D-03                                        A
      IR=7                                                              A    732
      GO TO 10                                                          A    733
  320 IF(.NOT.EXACT) GO TO 15                                           A    734
      TOL=(FNOT-FMIN)*25.D-03                                             A    7
      IF(TOL.LT.NOISE) TOL=NOISE                                        A    736
      IF(DABS(FUNC-FQUAD).LT.TOL) GO TO 15                               A    73
      AIF=ALF.EQ.AMIN.OR.ALF.EQ.ANEXT
      IF(AIF) GO TO 150
      IF(.NOT.AIF) GO TO 305
      END                                                                A    73


      SUBROUTINE TRUDGE(K,FUNC,M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,S,
     1    H)                                        
      IMPLICIT REAL*8(A-H,O-Z)
C                                                                       A    216
C     PROGRAM MINIMIZES A FUNCTION OF N VARIABLES BY A MODFIED POWELL   A    217
C     METHOD OF SEARCH ALONG CONJUGATE DIRECTIONS. COLUMNS OF V MATRIX  A    218
C     ARE CURRENT SEARCH DIRECTIONS.COLUMNS 1 TO K ARE (APROXIMATELY)   A    219
C     CONJUGATE.                                                        A    220
C                                                                       A    221
C     WHEN SUBROUTINE IS ENTERED WITH RESTART=.TRUE. THEN VALUES OF K   A    222
C     AND ALL QUANTITIES IN COMMON T1 SHOULD HAVE BEEN SET PREVIOUSLY   A    223
C     BY THE CALLING PROGRAM.                                           A    224
C     WHEN TRUDGE IS CALLING WITH K NEGATIVE THEN ONLY REQUIRED VALUES  
C     ARE TIMLIM,TOLF,TOLR,NOISE,N,AND R.(NOISE NEEDED FOR SEARCH)      
C     IN THIS CASE INITIAL VALUES OF K,JSTART,RESTART,ALPH,V ARE IGNORED
C     IF TRUDGE IS CALLED WITH NON-NEGATIVE K AND RESTA-T=.FALSE. THEN  
C     VALUES OF TIMLIM,TOLF,TOLR,NOISE,N,R,AND V ARE REQUIRED.          
C     IN THIS CASE K IS RESET TO ZERO AFTER FIRST PASS.                 
C                                                                       
C     EXTERNAL FUNCTION FUNC(INOUT)                                     
C     INOUT=  2   FUNCTION READS MINIMAL STARTING DATA AND EVALUATES F(R
C     INOUT=  1   FUNCTION READS ALL RESTART DATA AND EVALUATES F(R)    
C     INOUT=  0   FUNCTION EVALUATES F(R)                               
C     INOUT= -1   FUNCTION PUNCHES DATA FOR FUTURE CALL FUNC(1)         
C     INOUT= -2   FUNCTION EVALUATES F(R) AND PRINTS OUTPUT.            
C                                                                       
      COMMON/T1/ N,NMAX,F, R(20),DR(20),CURVE,ALPHA,NOISE,EXACT,ALPH(20)
     1,RNOT(20), V(20,20), CURV(20),RESTART, JSTART, FNOT, TIMLIM,      
     2 TOLF,TOLR                                                       
      COMMON/EIGMAT/G(20,20),Q(20,20),SCRAP1(20),SCRAP2(20),SCRAP3(20)  
      DIMENSION P(20),U(20)                                             
      LOGICAL EXACT, RESTART                                           
      REAL*8 NOISE
      EXTERNAL FUNC                                                    
      NMAX=20                                                         
      IF(N.GT.NMAX) STOP
      FN=N                                                             
      INOUT=0                                                          
      IF(K.LT.0) GO TO 15                                             
      IF(.NOT.RESTART) CALL FUNC(INOUT,M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,
     1    RZ2,B12,S,H)                                
      DO 10 I=1,N                                                       
   10 IF(ALPH(I).LE.0.D+00) ALPH(I)=1.D-01                                   
      IF(JSTART.LE.0.OR.JSTART.GT.N) JSTART=N                           
      J=JSTART                                                          
      RESTART=RESTART.AND.K.NE.N                                        A    256
      GO TO 100                                                         A    257
C                                                                       A    258
C     IF K NEGATIVE SET INITIAL V MATRIX TO IDENTITY,STANDARD VALUES.   A    259
C                                                                       A    260
   15 K=N                                                               A    261
      RESTART=.FALSE.                                                   A    262
      CALL FUNC(INOUT,M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,S,H)                                                  A    263
      DO 30 I=1,N                                                       A    264
      ALPH(I)=1.D-01                                                      A    2
      DO 20 J=1,N                                                       A    266
   20 V(I,J)=0.D+00                                                        A
      CURV(I)=0.D+00                                                       A
   30 V(I,I)=1.D+00                                                        A
      J=N                                                               A    270
C                                                                        A    27
C     SEARCH ALONG J-TH DIRECTION FOR LOCAL MINIMUM                     A    272
C                                                                        A    27
  100 DO 105 I=1,N                                                      A    274
  105 DR(I)=V(I,J)                                                      A    275
      EXACT=J.LE.K                                                      A    276
      ALPHA=ALPH(J)                                                     A    277
      ALFSET=ALPHA                                                      A    278
      CALL SEARCH(LEAVE,FUNC,M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,S,H)                                          
C      WRITE(6,120) CURVE,ALPHA                                          A    280
  120 FORMAT(17H  NEW CURVATURE =,E13.5,10H   ALPHA =,F12.7)            A    281
  121 FORMAT(17H  OLD CURVATURE =,E13.5,10H   ALFSET=,F12.7)            A    282
  122 FORMAT(20X,30H ABNORMAL RETURN FROM SEARCH  )                     A    283
      ALPHA=DABS(ALPHA)                                                  A    28
      IF(ALPHA.LT.5.D-01*ALFSET) ALFSET=5.D-01*ALFSET                         A
      IF(ALPHA.LT.25.D-02*ALFSET) ALFSET=5.D-1*ALFSET                         A
      IF(ALPHA.GT.2.D+00*ALFSET) ALFSET=2.D+00*ALFSET                         A
      IF(ALPHA.GT.4.D+00*ALFSET) ALFSET=2.D+00*ALFSET                         A
      WRITE(6,121) CURV(J),ALFSET                                       A    289
      CURV(J)=CURVE                                                     A    290
      ALPH(J)=ALFSET                                                    A    291
      IF(J.EQ.1) GO TO 175                                              A    292
      J=J-1                                                             A    293
      GO TO 100
C                                                                       A    296
C     CALL FUNC(-1) CAUSES ALL RESTART DATA TO BE OUTPUT,               A    297
C     THIS INCLUDES ALL DATA REQUIRED FOR FUTURE CALL FUNC(1).          A    298
C                                                                       A    299
  150 WRITE(6,155) K,TIMLIM,TOLF,TOLR,NOISE,J,FNOT                      A    300
  155 FORMAT(I10,F10.0,3E10.1,I10,F20.10)                               A    301
      INOUT=-2                                                           A    30
      CALL FUNC(INOUT,M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,S,H)                                                  A    303
      DO 160 I=1,N                                                      A    306
  160 WRITE(6,161)I,RNOT(I),CURV(I),ALPH(I)                             A    307
  161 FORMAT(I5,3E20.10)                                                A    308
      I1=1                                                              A    309
  165 I2=I1+7                                                            A    31
      IF(I2.GT.N) I2=N                                                  A    311
      DO 166 L=1,N                                                      A    312
  166 WRITE(6,168)(V(L,I),I=I1,I2)                                      A    313
  167 FORMAT(5X,8F15.10)                                                A    314
  168 FORMAT(8F10.7)                                                    A    315
      I1=I2+1                                                           A    316
      IF(I1.LE.N) GO TO 165                                             A    317
      K=0                                                               A    318
      RETURN
C                                                                       A    321
C     END OF SEQUENCE OF SEARCH ALONG ALL DIRECTIONS IN V MATRIX.       A    322
C                                                                       A    323
  175 IF(RESTART) GO TO 185                                             A    324
      DO 180 I=1,N                                                      A    325
  180 RNOT(I)=R(I)                                                      A    326
      RESTART=.TRUE.                                                    A    327
      FNOT=F                                                            A    328
      K=0                                                               A    329
      J=N                                                               A    330
      GO TO 100                                                         A    331
C                                                                        A    33
C     GENERATE NEW CONJUGATE DIRECTION IF NOT CONVERGED.                A    333
C                                                                        A    33
  185 DRDR=0.D+00                                                          A
      DO 188 I=1,N                                                      A    336
      DUM=R(I)-RNOT(I)                                                  A    337
      DR(I)=DUM                                                         A    338
  188 DRDR=DUM*DUM+DRDR                                                 A    339
      ALFSET=DSQRT(DRDR)                                                  A    3
      ALPHA=ALFSET                                                      A    341
      IF(ALPHA.EQ.0.D+00) GO TO 192                                        A
      DO 190 I=1,N                                                      A    343
      U(I)=DR(I)/ALFSET                                                 A    344
  190 DR(I)=U(I)                                                        A    345
      EXACT=.TRUE.                                                      A    346
      CALL SEARCH(LEAVE,FUNC,M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,S,H)
      IF(LEAVE.NE.0) WRITE(6,122)                                       A    348
C      WRITE(6,120) CURVE,ALPHA                                          A    349
      ALPHA=DABS(ALPHA)                                                  A    35
  192 DELF=FNOT-F                                                       A    351
      FNOT=F                                                            A    352
      DRDR=0.D+00                                                          A
      DO 195 I=1,N                                                      A    354
      DRDR=(R(I)-RNOT(I))**2+DRDR
  195 RNOT(I)=R(I)                                                      A    356
      DELR=DSQRT(DRDR/FN)                                                A    35
      WRITE(6,210) F, DELF, DELR                                        A    358
  210 FORMAT(4H  F=,F20.10,10H  DELTA F=,E15.5,10H  DELTA R=,E15.5)     A    359
      IF(DELF.GT.TOLF.OR.DELR.GT.TOLR) GO TO 220                        A    360
      WRITE(6,215)                                                      A    361
  215 FORMAT(//,' CONVERGED')                                           A    362
      GO TO 150                                                         A    363
  220 IF(K.GT.0) GO TO 230                                              A    364
      K=1                                                               A    365
      K1=0                                                              A    366
      DO 225 I=1,N                                                      A    367
  225 R(I)=DR(I)                                                        A    368
      DR(N)=CURVE                                                       A    369
      P(1)=ALPHA                                                        A    370
      GO TO 375                                                         A    371
C                                                                       A    372
C     CALCULATE K+1 Q VECTORS.Q(I) DOT V(J)=KRONIKER DELTA(I,J),J=1,K   A    373
C                                                                       A    374
  230 DO 250 J=1,K                                                      A    375
      DOT=0.D+00                                                           A
      DO 240 I=1,N                                                      A    377
  240 DOT=DR(I)*V(I,J)+DOT                                              A    378
      P(J)=DOT                                                          A    379
      DO 250 I=1,N                                                      A    380
  250 DR(I)=DR(I)-DOT*V(I,J)                                            A    381
      DOT=0.D+00                                                           A
      DO 260 I=1,N                                                      A    383
  260 DOT=DR(I)*U(I)+DOT                                                A    384
      K1=K                                                              A    385
      K=K1+1                                                            A    386
      RESTART=K.NE.N                                                    A    387
      DO 270 I=1,N                                                      A    388
      DUM=DR(I)/DOT                                                     A    389
      Q(I,K)=DUM                                                        A    390
      DO 270 J=1,K1                                                     A    391
  270 Q(I,J)=V(I,J)-P(J)*DUM                                            A    392
C                                                                       A    393
C     G MATRIX=SUM L=1,K+1 OF Q(L)*CURV(L)*TRANSPOSE Q(L)               A    394
C                                                                        A    39
      DO 280 I=1,N                                                      A    396
      DUM=CURVE*Q(I,K)                                                  A    397
      DO 280 J=1,N                                                      A    398
  280 G(I,J)=DUM*Q(J,K)                                                 A    399
      DO 300 L=1,K1                                                     A    400
      DO 300 I=1,N                                                      A    401
      DUM=CURV(L)*Q(I,L)                                                A    402
      DO 300 J=1,N                                                      A    403
  300 G(I,J)=DUM*Q(J,L)+G(I,J)                                          A    404
      CALL EIGN(N,DR)                                                   A    405
C                                                                       A    406
C     SET P(J) EQUAL TO CORRESPONDING ALPH VALUE.                       A    407
C                                                                       A    408
      DO 350 J=1,K                                                      A    409
      JJ=N+1-J                                                          A    410
      P(J)=ALPH(J)                                                      A    411
      DOT=0.D+00                                                           A
      DO 310 I=1,N                                                      A    413
      R(I)=Q(I,JJ)                                                      A    414
  310 DOT=R(I)*V(I,J)+DOT                                               A    415
      DOT=DABS(DOT)                                                      A    41
      IF(DOT.GT.7.D-01) GO TO 350                                          A
      LMAX=J                                                             A    41
      DOTMAX=DOT                                                        A    419
      DO 340 L=1,K1                                                     A    420
      IF(L.EQ.J) GO TO 340                                              A    421
      DOT=0.D+00                                                           A
      DO 320 I=1,N                                                      A    423
  320 DOT=R(I)*V(I,L)+DOT                                               A    424
      DOT=DABS(DOT)                                                      A    42
      IF(DOT.LT.7.D-01) GO TO 330                                          A
      P(J)=ALPH(L)                                                      A    427
      GO TO 350                                                         A    428
  330 IF(DOT.LE.DOTMAX) GO TO 340                                       A    429
      DOTMAX=DOT                                                        A    430
      LMAX=L                                                            A    431
  340 CONTINUE                                                          A    432
      P(J)=ALPH(LMAX)                                                   A    433
      DOT=0.D+00                                                           A
      DO 345 I=1,N                                                      A    435
  345 DOT=R(I)*U(I)+DOT                                                 A    436
      IF(DABS(DOT).GT.DOTMAX) P(J)=ALPHA                                 A    43
  350 CONTINUE                                                          A    438
C                                                                        A    43
C     SET NEV CONJUGATE VECTORS IN V EXCEPT KTH VECTOR IS SAVED IN R.   A    440
C                                                                       A    441
      DO 370 J=1,K1                                                     A    442
      JJ=N+1-J                                                          A    443
      DO 360 I=1,N                                                      A    444
  360 V(I,J)=Q(I,JJ)                                                    A    445
      ALPH(J)=P(J)                                                       A    44
  370 CURV(J)=DR(JJ)                                                    A    447
C                                                                        A    44
C     ORTHOGONALIZE REMAINING SEARCH DIRECTIONS TO CONJUGATE DIRECTIONS.A    449
C                                                                       A    450
  375 IF(K.EQ.N) GO TO 450                                              A    451
      DO 420 J=K,N                                                      A    452
      DOT=0.D+00                                                           A
      DO 380 I=1,N                                                      A    454
  380 DOT=V(I,J)*R(I)+DOT                                               A    455
      DO 390 I=1,N                                                      A    456
  390 V(I,J)=V(I,J)-DOT*R(I)                                            A    457
      IF(K.EQ.1) GO TO 420                                              A    458
      DO 410 L=1,K1                                                     A    459
      DOT=0.D+00                                                           A
      DO 400 I=1,N                                                      A    461
  400 DOT=V(I,J)*V(I,L)+DOT                                             A    462
      DO 410 I=1,N                                                      A    463
  410 V(I,J)=V(I,J)-DOT*V(I,L)                                          A    464
  420 CONTINUE                                                          A    465
      CMAX=0.D+00                                                          A
      DO 425 L=K,N                                                      A    467
  425 IF(CURV(L).GT.CMAX) CMAX=CURV(L)                                  A    468
C                                                                       A    469
C     CONSTRUCT G MATRIX FOR REMAINING SEARCH DIRECTIONS.               A    470
C                                                                        A    47
      CMAX=8.D+00*CMAX                                                     A
      DO 440 I=1,N                                                      A    473
      DO 440 J=1,I                                                      A    474
      DUM=-CMAX*R(I)*R(J)                                               A    475
      DO 430 L=K,N                                                      A    476
  430 DUM=CURV(L)*V(I,L)*V(J,L)+DUM                                     A    477
      IF(K1.LE.0) GO TO 435                                             A    478
      DO 433 L=1,K1                                                     A    479
  433 DUM=-CMAX*V(I,L)*V(J,L)+DUM                                       A    480
  435 G(I,J)=DUM                                                        A    481
  440 G(J,I)=DUM                                                        A    482
C                                                                       A    483
C     SET KTH CONJUGATE DIRECTIONS BACK INTO V MATRIX.                  A    484
C                                                                       A    485
  450 DO 460 I=1,N                                                      A    486
      V(I,K)=R(I)                                                       A    487
  460 R(I)=RNOT(I)                                                      A    488
      ALPH(K)=P(K)                                                       A    48
      CURV(K)=DR(N+1-K)                                                 A    490
      IF(K.EQ.N) GO TO 480                                              A    491
      CALL EIGN(N,DR)                                                   A    492
C                                                                       A    493
C     SET REMAINING SEARCH DIRECTIONS INTO V MATRIX.                    A    494
C                                                                        A    49
      K1=K+1                                                            A    496
      DO 470 J=K1,N                                                     A    497
      CURV(J)=DR(J)                                                     A    498
      DO 470 I=1,N                                                      A    499
  470 V(I,J)=Q(I,J)                                                     A    500
  480 J=N                                                               A    501
      GO TO 100
      END                                                                A    50
