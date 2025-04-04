

      SUBROUTINE TRUDGE(K,FUNC)                                         A    215
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
C     WHEN TRUDGE IS CALLING WITH K NEGATIVE THEN ONLY REQUIRED VALUES  A    225
C     ARE TIMLIM,TOLF,TOLR,NOISE,N,AND R.(NOISE NEEDED FOR SEARCH)      A    226
C     IN THIS CASE INITIAL VALUES OF K,JSTART,RESTART,ALPH,V ARE IGNOREDA    227
C     IF TRUDGE IS CALLED WITH NON-NEGATIVE K AND RESTA-T=.FALSE. THEN  A    228
C     VALUES OF TIMLIM,TOLF,TOLR,NOISE,N,R,AND V ARE REQUIRED.          A    229
C     IN THIS CASE K IS RESET TO ZERO AFTER FIRST PASS.                 A    230
C                                                                       A    231
C     EXTERNAL FUNCTION FUNC(INOUT)                                     A    232
C     INOUT=  2   FUNCTION READS MINIMAL STARTING DATA AND EVALUATES F(RA    233
C     INOUT=  1   FUNCTION READS ALL RESTART DATA AND EVALUATES F(R)    A    234
C     INOUT=  0   FUNCTION EVALUATES F(R)                               A    235
C     INOUT= -1   FUNCTION PUNCHES DATA FOR FUTURE CALL FUNC(1)         A    236
C     INOUT= -2   FUNCTION EVALUATES F(R) AND PRINTS OUTPUT.            A    237
C                                                                        A    23
      COMMON/T1/ N,NMAX,F, R(20),DR(20),CURVE,ALPHA,NOISE,EXACT,ALPH(20)A    239
     1,RNOT(20), V(20,20), CURV(20),RESTART, JSTART, FNOT, TIMLIM,      A    240
     2 TOLF,TOLR                                                        A    241
      COMMON/EIGMAT/G(20,20),Q(20,20),SCRAP1(20),SCRAP2(20),SCRAP3(20)  A    242
      DIMENSION P(20),U(20)                                             A    243
      LOGICAL EXACT, RESTART                                            A    244
      REAL*8 NOISE
      EXTERNAL FUNC                                                     A    245
      NMAX=20                                                           A    246
      IF(N.GT.NMAX) STOP
      FN=N                                                              A    248
      INOUT=0                                                           A    249
      IF(K.LT.0) GO TO 15                                               A    250
      IF(.NOT.RESTART) CALL FUNC(INOUT)                                 A    251
      DO 10 I=1,N                                                       A    252
   10 IF(ALPH(I).LE.0.D+00) ALPH(I)=1.D-01                                   A
      IF(JSTART.LE.0.OR.JSTART.GT.N) JSTART=N                           A    254
      J=JSTART                                                          A    255
      RESTART=RESTART.AND.K.NE.N                                        A    256
      GO TO 100                                                         A    257
C                                                                       A    258
C     IF K NEGATIVE SET INITIAL V MATRIX TO IDENTITY,STANDARD VALUES.   A    259
C                                                                       A    260
   15 K=N                                                               A    261
      RESTART=.FALSE.                                                   A    262
      CALL FUNC(INOUT)                                                  A    263
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
      CALL SEARCH(LEAVE,FUNC)                                           A    279
      WRITE(6,120) CURVE,ALPHA                                          A    280
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
      CALL FUNC(INOUT)                                                  A    303
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
      CALL SEARCH(LEAVE,FUNC)
      IF(LEAVE.NE.0) WRITE(6,122)                                       A    348
      WRITE(6,120) CURVE,ALPHA                                          A    349
      ALPHA=DABS(ALPHA)                                                  A    35
  192 DELF=FNOT-F                                                       A    351
      FNOT=F                                                            A    352
      DRDR=0.D+00                                                          A
      DO 195 I=1,N                                                      A    354
      DRDR=(R(I)-RNOT(I))**2+DRDR
  195 RNOT(I)=R(I)                                                      A    356
      DELR=DSQRT(DRDR/FN)                                                A    35
      WRITE(6,210) F, DELF, DELR
	
   
     
     
                                             A    358
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
