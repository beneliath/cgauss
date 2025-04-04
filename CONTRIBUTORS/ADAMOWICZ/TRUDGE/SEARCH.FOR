      SUBROUTINE SEARCH(LEAVE,GAUSS)                                    A    505
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
      LOOKMAX=30                                                        A    544
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
      CALL GAUSS(INOUT)                                                 A    571
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
   29 WRITE(6,27) LOOKS, ALF, FUNC, AMIN, FMIN, ANEXT, FNEXT            A    588
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
