	TITLE	C:\users\dgilmore\Graduate.School\Adamowicz.Research\H2ANL\ENERGY.FOR
	.386P
include listing.inc
if @Version gt 510
.model FLAT
else
_DATA	SEGMENT DWORD USE32 PUBLIC 'DATA'
_DATA	ENDS
_CONST	SEGMENT DWORD USE32 PUBLIC 'DATA'
_CONST	ENDS
_BSS	SEGMENT DWORD USE32 PUBLIC 'BSS'
_BSS	ENDS
_TEXT	SEGMENT PARA USE32 PUBLIC 'CODE'
_TEXT	ENDS
$$SYMBOLS	SEGMENT BYTE USE32 'DEBSYM'
$$SYMBOLS	ENDS
$$TYPES	SEGMENT BYTE USE32 'DEBTYP'
$$TYPES	ENDS
FLAT	GROUP _DATA
	ASSUME	CS: FLAT, DS: FLAT, SS: FLAT
endif
EXTRN	__FFljj:NEAR
EXTRN	_BLANKM@12:NEAR
EXTRN	__Deallocate_i:NEAR
EXTRN	_make_info:NEAR
EXTRN	_crash:NEAR
EXTRN	__Allocate_i:NEAR
EXTRN	_VECH@16:NEAR
EXTRN	_YACOBI@24:NEAR
EXTRN	_DEVEC@16:NEAR
EXTRN	_TRANSPOSE@20:NEAR
EXTRN	_MATRIXM@28:NEAR
_DATA	SEGMENT
_C:\users\dgilmore\Graduate.School\Adamowicz.Research\H2ANL\ENERGY.FOR DB 
	DB	'C:\users\dgilmore\Graduate.School\Adamowicz.Research\H2ANL\EN'
	DB	'ERGY.FOR', 00H
_DATA	ENDS
PUBLIC	_ENERGY@24
EXTRN	___FHlink@4:NEAR
EXTRN	___FHunlink@0:NEAR
EXTRN	___FHovfl@0:NEAR
EXTRN	__fltused:NEAR
EXTRN	__FIsqrt:NEAR
EXTRN	___FHranger@0:NEAR
_BSS	SEGMENT
_SS	DQ	010000H DUP (?)
_SHALF	DQ	010000H DUP (?)
_SNHALF	DQ	010000H DUP (?)
_HH	DQ	010000H DUP (?)
_XT	DQ	010000H DUP (?)
_XY	DQ	010000H DUP (?)
_W	DQ	010000H DUP (?)
_IX	DD	01H DUP (?)
	ALIGN	8

_CC	DQ	01H DUP (?)
_C	DQ	010000H DUP (?)
_SNORM	DQ	00H DUP (?)
_MMAX	DD	01H DUP (?)
_MNMAX	DD	01H DUP (?)
_NPMAX	DD	01H DUP (?)
_NNP1D2	DD	01H DUP (?)
_AA	DQ	0101d0H DUP (?)
_TP	DQ	0101d0H DUP (?)
_BIG	DQ	0100H DUP (?)
_JB	DD	0100H DUP (?)
_UT	DQ	010000H DUP (?)
_U	DQ	010000H DUP (?)
_TEM	DQ	010000H DUP (?)
_BSS	ENDS
_DATA	SEGMENT
	ORG $+2
_7SNORM	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H
_.SG42_	DB	00H
	ORG $+3
_.SG43_	DB	00H
	ORG $+7
$T98374	DB	'_ENERGY@24', 00H
_DATA	ENDS
_CONST	SEGMENT
$T98385	DQ	00000000000000000r		; 0
$T98406	DQ	03ff0000000000000r		; 1
_CONST	ENDS
_TEXT	SEGMENT
; File C:\users\dgilmore\Graduate.School\Adamowicz.Research\H2ANL\ENERGY.FOR
$T98306 = -140
$T98307 = -140
$T98308 = -144
$T98373 = -136
$T98309 = -140
$T98310 = -144
$T98311 = -140
$T98312 = -144
$T98313 = -148
$T98318 = -140
$T98319 = -144
$T98320 = -148
$T98321 = -140
$T98322 = -144
$T98323 = -148
$T98324 = -140
$T98326 = -140
$T98327 = -144
$T98328 = -148
$T98329 = -140
_M$ = 8
$T98330 = -144
_7H$ = -44
$T98331 = -148
_H$ = 12
$T98332 = -140
_7S$ = -88
$T98333 = -144
_S$ = 16
$T98334 = -140
_EMIN$ = 20
$T98335 = -144
_7CO$ = -120
$T98336 = -148
_CO$ = 24
$T98337 = -140
_KLM$ = 28
$T98338 = -140
$T98339 = -140
$T98340 = -144
$T98341 = -140
$T98342 = -144
$T98343 = -148
$T98344 = -140
$T98345 = -144
$T98346 = -148
$T98348 = -140
$T98349 = -144
$T98350 = -148
_ENERGY@24 PROC NEAR

; 6    :       SUBROUTINE ENERGY(M,H,S,EMIN,CO,KLM)

	sub	esp, 156				; 0000009cH
	lea	eax, DWORD PTR $T98373[esp+156]
	push	ebx
	push	esi
	push	edi
	push	ebp
	push	eax
	call	___FHlink@4
	mov	DWORD PTR $T98373[esp+176], OFFSET FLAT:_C:\users\dgilmore\Graduate.School\Adamowicz.Research\H2ANL\ENERGY.FOR
	mov	DWORD PTR $T98373[esp+180], OFFSET FLAT:$T98374
	mov	ecx, DWORD PTR _M$[esp+168]
	lea	eax, DWORD PTR _7CO$[esp+172]
	mov	DWORD PTR $T98373[esp+184], 7

; 7    :       IMPLICIT REAL*8 (A-H,O-Z)

	mov	edx, DWORD PTR [ecx]
	push	edx
	push	1
	push	1
	push	eax
	call	_make_info
	mov	ecx, DWORD PTR _M$[esp+184]
	add	esp, 16					; 00000010H
	mov	esi, DWORD PTR [ecx]
	push	esi
	mov	edi, esi
	push	1
	lea	eax, DWORD PTR _7S$[esp+180]
	push	esi
	push	1
	push	2
	push	eax
	call	_make_info
	lea	eax, DWORD PTR _7H$[esp+196]
	add	esp, 24					; 00000018H
	push	edi
	push	1
	push	esi
	push	1
	push	2
	mov	esi, 37					; 00000025H
	push	eax
	call	_make_info
	mov	eax, DWORD PTR _M$[esp+192]
	mov	DWORD PTR $T98373[esp+208], esi
	add	esp, 24					; 00000018H
	mov	eax, DWORD PTR [eax]

; 8    : 	INTEGER M
; 9    : C 
; 10   : C     N=number of particles
; 11   : C     M= number of functions used in the expansion
; 12   : C     IPER=number of permutations
; 13   : C     NP = N*(N-1)/2 = number of interacting pairs
; 14   : C     MN = M*N = number of orbital exponents
; 15   : C     MNP = M*N*(N-1)/2 = number of correlation exponents
; 16   : C
; 17   :       PARAMETER ( NMAX = 5 )
; 18   : C      PARAMETER ( MMAX = 3 )
; 19   : C      PARAMETER ( NPMAX = NMAX*(NMAX-1)/2 )
; 20   : C      PARAMETER ( MNMAX = MMAX*NMAX )
; 21   : C 
; 22   :       REAL*8, DIMENSION(M) ::	CO
; 23   : 	REAL*8, DIMENSION(M,M) :: S,H
; 24   : 	REAL*8 SNORM[ALLOCATABLE](:,:)
; 25   : 
; 26   : C      DIMENSION  SNORM(MMAX,MMAX)
; 27   :       DIMENSION AA(66000),BIG(256),JB(256),TEMP(256),TP(66000)
; 28   :       DIMENSION UT(256,256),U(256,256),TEM(256,256),SS(256,256)
; 29   :       DIMENSION TEM1(256,256),TEM2(256,256),TEM3(256,256)
; 30   :       DIMENSION SHALF(256,256),SNHALF(256,256),CHECKS(256,256)
; 31   :       DIMENSION HH(256,256),XY(256,256),XT(256,256),W(256,256)
; 32   :       DIMENSION CT(256,256),CHECKH(256,256),C(256,256),AB(256)
; 33   : C
; 34   : C      COMMON /MOLPARAM/ N, M, NP 
; 35   : C      COMMON /FFF/ S, H
; 36   : c      COMMON /NORM/ SNORM
; 37   : 	ALLOCATE(SNORM(M,M))

	push	eax
	push	1
	push	eax
	push	1
	push	2
	push	OFFSET FLAT:_7SNORM
	call	_make_info
	add	esp, 24					; 00000018H
	mov	ecx, DWORD PTR _SNORM
	test	ecx, ecx
	je	SHORT $L40
	push	OFFSET FLAT:_.SG43_
	push	OFFSET FLAT:_.SG42_
	call	_crash
	add	esp, 8
$L40:
	mov	DWORD PTR $T98373[esp+184], 37		; 00000025H
	push	0
	push	8
	push	OFFSET FLAT:_7SNORM
	call	__Allocate_i
	mov	ecx, DWORD PTR _M$[esp+180]
	mov	DWORD PTR _SNORM, eax
	mov	DWORD PTR $T98373[esp+196], 38		; 00000026H
	add	esp, 12					; 0000000cH
	mov	edx, DWORD PTR [ecx]

; 38   : 	MMAX=M

	mov	DWORD PTR _MMAX, edx
	mov	DWORD PTR $T98373[esp+184], 39		; 00000027H
	mov	esi, DWORD PTR [ecx]

; 39   : 	MNMAX=M*NMAX

	mov	edi, esi
	imul	edi, 5
	jno	SHORT $L98376
	call	___FHovfl@0
$L98376:
	mov	ecx, DWORD PTR _KLM$[esp+168]
	mov	DWORD PTR _MNMAX, edi
	mov	DWORD PTR $T98373[esp+184], 40		; 00000028H

; 40   : 	NPMAX=NMAX*(NMAX-1)/2

	mov	DWORD PTR _NPMAX, 10			; 0000000aH
	mov	DWORD PTR $T98373[esp+184], 42		; 0000002aH

; 41   : C
; 42   :       IF(KLM.EQ.5) GOTO 1414

	cmp	DWORD PTR [ecx], 5
	je	$L50
	mov	DWORD PTR $T98373[esp+184], 92		; 0000005cH

; 43   : C
; 44   : C     PRINT-OUT THE S & H MATRICES (UNNORMALIZED):
; 45   : C     Uncomment all code between the (===) lines to execute...
; 46   : C     ===============================================================
; 47   : C      CALL BIOUT(' ')
; 48   : C      CALL BIOUT(' Overlap (S) matrix: ')
; 49   : C      CALL BIOUT(' ------------------- ')
; 50   : C          CALL TAB(S,M,M,MMAX,MMAX)
; 51   : C      CALL BIOUT(' ')
; 52   : C      CALL BIOUT(' Hamiltonian (H) matrix: ')
; 53   : C      CALL BIOUT(' ----------------------- ')
; 54   : C          CALL TAB(H,M,M,MMAX,MMAX)
; 55   : C      CALL BIOUT(' ')
; 56   : C     ===============================================================
; 57   : C      STOP
; 58   : C
; 59   : C     ---------------------------------------------------------------
; 60   : C     NORMALIZATION OF S & H MATRICES:
; 61   : C     ---------------------------------------------------------------
; 62   : C     Copy the contents of matrix S into matrix SNORM...
; 63   : c      CALL MATRIXCOPY(M,S,MMAX,SNORM,MMAX)
; 64   : C
; 65   : C     Normalization of the S & H matrix elements...
; 66   : c      DO 59 I=1,M
; 67   : c      DO 59 J=1,M
; 68   : c      S(I,J) = S(I,J)/( DSQRT(SNORM(I,I))*DSQRT(SNORM(J,J)) )
; 69   : c      H(I,J) = H(I,J)/( DSQRT(SNORM(I,I))*DSQRT(SNORM(J,J)) )
; 70   : C      S(I,J) = S(I,J)/( (SNORM(I,I)**(1/2))*(SNORM(J,J)**(1/2)) )
; 71   : C      H(I,J) = H(I,J)/( (SNORM(I,I)**(1/2))*(SNORM(J,J)**(1/2)) )
; 72   : c 59   CONTINUE
; 73   : C
; 74   : C
; 75   : C     PRINT-OUTS OF NORMALIZED S & H MATRICES:
; 76   : C     Uncomment all code between the (===) lines to execute...
; 77   : C     ===============================================================
; 78   : C      CALL BIOUT(' ')
; 79   : C      CALL BIOUT(' Normalized Overlap (S) Matrix: ')
; 80   : C      CALL BIOUT(' ------------------------------ ')
; 81   : C         CALL TAB(S,M,M,MMAX,MMAX)
; 82   : C      CALL BIOUT(' ')
; 83   : C      CALL BIOUT(' Normalized Hamiltonian (H) Matrix: ')
; 84   : C      CALL BIOUT(' ---------------------------------- ')
; 85   : C         CALL TAB(H,M,M,MMAX,MMAX)
; 86   : C      CALL BIOUT(' ')
; 87   : C     ===============================================================
; 88   : C      STOP
; 89   : C
; 90   : C
; 91   : C     NNP1D2=number of lower triangular elements for YACOBI
; 92   :       NNP1D2=M*(M+1)/2

	mov	edi, esi
	inc	edi
	jno	SHORT $L98378
	call	___FHovfl@0
$L98378:
	imul	esi, edi
	jno	SHORT $L98379
	call	___FHovfl@0
$L98379:
	mov	eax, esi
	push	OFFSET FLAT:_AA
	mov	esi, DWORD PTR _S$[esp+172]
	push	OFFSET FLAT:_MMAX
	mov	ecx, DWORD PTR _M$[esp+176]
	push	esi
	cdq
	sub	eax, edx
	push	ecx
	mov	DWORD PTR $T98373[esp+200], 98		; 00000062H
	sar	eax, 1
	mov	DWORD PTR _NNP1D2, eax

; 93   : C
; 94   : C     ---------------------------------------------------------------
; 95   : C     DIAGNOLIZE THE S MATRIX
; 96   : C     ---------------------------------------------------------------
; 97   : C     Define the one dimensional array AA(II) as the S matrix for YACOBI...
; 98   :       CALL VECH (M,S,MMAX,AA)

	call	_VECH@16
	mov	ecx, DWORD PTR _M$[esp+168]
	push	OFFSET FLAT:_JB
	mov	DWORD PTR $T98373[esp+188], 102		; 00000066H

; 99   : C
; 100  : C
; 101  : C     Call YACOBI to get the TRANSPOSE of the ROTATION MATRIX...
; 102  :       CALL YACOBI(AA,TP,M,NNP1D2,BIG,JB)

	push	OFFSET FLAT:_BIG
	push	OFFSET FLAT:_NNP1D2
	push	ecx
	mov	edi, 256				; 00000100H
	push	OFFSET FLAT:_TP
	push	OFFSET FLAT:_AA
	call	_YACOBI@24
	lea	ecx, DWORD PTR $T98306[esp+172]
	mov	eax, DWORD PTR _M$[esp+168]
	mov	DWORD PTR $T98306[esp+172], edi
	push	ecx
	mov	DWORD PTR $T98373[esp+188], 105		; 00000069H

; 103  : C
; 104  : C     Construct the transpose matrix UT FROM the vector TP...
; 105  :       CALL DEVEC (M,TP,UT,256)

	push	OFFSET FLAT:_UT
	push	OFFSET FLAT:_TP
	push	eax
	call	_DEVEC@16
	lea	ecx, DWORD PTR $T98307[esp+172]
	lea	edx, DWORD PTR $T98308[esp+172]
	mov	DWORD PTR $T98307[esp+172], edi
	push	ecx
	mov	eax, DWORD PTR _M$[esp+172]
	push	OFFSET FLAT:_U
	mov	DWORD PTR $T98308[esp+180], edi
	push	edx
	mov	DWORD PTR $T98373[esp+196], 120		; 00000078H

; 106  : C
; 107  : C
; 108  : C     PRINT-OUT THE TRANSPOSE ROTATION MATRIX:
; 109  : C     Uncomment all code between the (===) lines to execute...
; 110  : C     ===============================================================
; 111  : C      CALL BIOUT(' ')
; 112  : C      CALL BIOUT(' TRANSPOSE matrix ')
; 113  : C      CALL BIOUT(' ---------------------- ')
; 114  : C         CALL TAB(UT,M,M,100,100)
; 115  : C     ===============================================================
; 116  : C      STOP
; 117  : C
; 118  : C
; 119  : C     Construct the regular rotation matrix from its transpose (UT)...
; 120  :       CALL TRANSPOSE(M,UT,256,U,256)

	push	OFFSET FLAT:_UT
	push	eax
	call	_TRANSPOSE@20
	lea	ecx, DWORD PTR $T98309[esp+172]
	mov	DWORD PTR $T98309[esp+172], edi
	lea	edx, DWORD PTR $T98310[esp+172]
	push	ecx
	mov	DWORD PTR $T98310[esp+176], edi
	push	OFFSET FLAT:_TEM
	mov	eax, DWORD PTR _M$[esp+176]
	push	OFFSET FLAT:_MMAX
	mov	DWORD PTR $T98373[esp+196], 137		; 00000089H

; 121  : C
; 122  : C
; 123  : C     SHOW THAT THE ROTATION IS A UNITARY TRANSFORMATION:
; 124  : C     Uncomment all code between the (===) lines to execute...
; 125  : C     ===============================================================
; 126  : C      CALL BIOUT(' ')
; 127  : C      CALL BIOUT(' Unitary Transformation?: ')
; 128  : C      CALL BIOUT(' ------------------------ ' )
; 129  : C          CALL MATRIXM (M,U,100,UT,100,TEM3,100)
; 130  : C          CALL TAB(TEM3,M,M,100,100)
; 131  : C      CALL BIOUT(' ')
; 132  : C     ===============================================================
; 133  : C      STOP
; 134  : C
; 135  : C
; 136  : C     Perform the operations to diagonalize S
; 137  :       CALL MATRIXM(M,UT,256,S,MMAX,TEM,256)

	push	esi
	push	edx
	push	OFFSET FLAT:_UT
	push	eax
	call	_MATRIXM@28
	lea	ecx, DWORD PTR $T98311[esp+172]
	mov	DWORD PTR $T98311[esp+172], edi
	lea	edx, DWORD PTR $T98312[esp+172]
	push	ecx
	mov	DWORD PTR $T98312[esp+176], edi
	push	OFFSET FLAT:_SS
	mov	DWORD PTR $T98313[esp+180], edi
	push	edx
	lea	edi, DWORD PTR $T98313[esp+184]
	push	OFFSET FLAT:_U
	mov	eax, DWORD PTR _M$[esp+184]
	push	edi
	mov	DWORD PTR $T98373[esp+204], 138		; 0000008aH

; 138  :       CALL MATRIXM(M,TEM,256,U,256,SS,256)

	push	OFFSET FLAT:_TEM
	push	eax
	mov	esi, 1
	call	_MATRIXM@28
	mov	ecx, 160				; 000000a0H
	mov	eax, DWORD PTR _M$[esp+168]
	mov	DWORD PTR $T98373[esp+184], ecx

; 139  : 
; 140  : C
; 141  : C
; 142  : C     PRINT-OUT THE DIAGONALIZED S MATRIX:
; 143  : C     Uncomment all code between the (===) lines to execute...
; 144  : C     ===============================================================
; 145  : C      CALL BIOUT(' ')
; 146  : C      CALL BIOUT(' Diagonalized (S) matrix: ')
; 147  : C      CALL BIOUT(' ----------------------- ' )
; 148  : C          CALL TAB(SS,M,M,100,100)
; 149  : C      CALL BIOUT(' ')
; 150  : C     ===============================================================
; 151  : C      STOP
; 152  : C
; 153  : C
; 154  : C     ---------------------------------------------------------------
; 155  : C     DIAGONALIZE THE H MATRIX:  PART I
; 156  : C     ---------------------------------------------------------------
; 157  : C     Create S**(1/2) & S**(-1/2) matrices...
; 158  : C
; 159  : C
; 160  :       DO 3130 I=1,M

	mov	edx, DWORD PTR [eax]
	mov	DWORD PTR -140+[esp+172], edx
	test	edx, edx
	jle	$L66
$L98314:
	mov	ecx, 161				; 000000a1H
	mov	eax, DWORD PTR _M$[esp+168]
	mov	ebx, 1
	mov	DWORD PTR $T98373[esp+184], ecx

; 161  :           DO 3130 J=1,M

	mov	edi, DWORD PTR [eax]
	test	edi, edi
	jle	SHORT $L67
$L98315:
	mov	DWORD PTR $T98373[esp+184], 162		; 000000a2H

; 162  :               SHALF(I,J)=0.0D+00

	cmp	esi, 256				; 00000100H
	jg	SHORT $L98380
	cmp	esi, 1
	jge	SHORT $L98381
$L98380:
	call	___FHranger@0
$L98381:
	mov	ebp, esi
	cmp	ebx, 256				; 00000100H
	jg	SHORT $L98383
	cmp	ebx, 1
	jge	SHORT $L98384
$L98383:
	call	___FHranger@0
$L98384:
	mov	ecx, ebx
	xor	eax, eax
	shl	ecx, 8
	dec	edi
	add	ecx, ebp
	inc	ebx
	mov	DWORD PTR $T98373[esp+184], 161		; 000000a1H
	mov	DWORD PTR _SHALF[ecx*8-2056], eax
	test	edi, edi
	mov	DWORD PTR _SHALF[ecx*8-2052], eax

; 161  :           DO 3130 J=1,M

	jne	SHORT $L98315
	mov	DWORD PTR $T98373[esp+184], 160		; 000000a0H
$L67:

; 139  : 
; 140  : C
; 141  : C
; 142  : C     PRINT-OUT THE DIAGONALIZED S MATRIX:
; 143  : C     Uncomment all code between the (===) lines to execute...
; 144  : C     ===============================================================
; 145  : C      CALL BIOUT(' ')
; 146  : C      CALL BIOUT(' Diagonalized (S) matrix: ')
; 147  : C      CALL BIOUT(' ----------------------- ' )
; 148  : C          CALL TAB(SS,M,M,100,100)
; 149  : C      CALL BIOUT(' ')
; 150  : C     ===============================================================
; 151  : C      STOP
; 152  : C
; 153  : C
; 154  : C     ---------------------------------------------------------------
; 155  : C     DIAGONALIZE THE H MATRIX:  PART I
; 156  : C     ---------------------------------------------------------------
; 157  : C     Create S**(1/2) & S**(-1/2) matrices...
; 158  : C
; 159  : C
; 160  :       DO 3130 I=1,M

	mov	eax, DWORD PTR -140+[esp+172]
	inc	esi
	mov	DWORD PTR $T98373[esp+184], 160		; 000000a0H
	dec	eax
	mov	DWORD PTR -140+[esp+172], eax
	test	eax, eax
	jne	$L98314
	mov	DWORD PTR $T98373[esp+184], 163		; 000000a3H
$L66:

; 163  : 3130          CONTINUE

	mov	ecx, 164				; 000000a4H
	mov	eax, DWORD PTR _M$[esp+168]
	mov	edi, 1
	mov	DWORD PTR $T98373[esp+184], ecx

; 164  :           DO 3140  II=1,M

	mov	esi, DWORD PTR [eax]
	test	esi, esi
	jle	$L73
$L98316:
	mov	DWORD PTR $T98373[esp+184], 167		; 000000a7H

; 165  : C      WRITE(6,*) 'SS (',II,',',II,') = ',SS(II,II)
; 166  : C      WRITE(7,*) 'SS (',II,',',II,') = ',SS(II,II)
; 167  :               SHALF(II,II)=DSQRT(SS(II,II))

	cmp	edi, 256				; 00000100H
	jg	SHORT $L98386
	cmp	edi, 1
	jge	SHORT $L98387
$L98386:
	call	___FHranger@0
$L98387:
	mov	DWORD PTR -140+[esp+172], edi
	cmp	edi, 256				; 00000100H
	jg	SHORT $L98388
	cmp	edi, 1
	jge	SHORT $L98389
$L98388:
	call	___FHranger@0
$L98389:
	cmp	edi, 256				; 00000100H
	mov	ebx, edi
	mov	ebp, edi
	jg	SHORT $L98390
	cmp	ebp, 1
	jge	SHORT $L98391
$L98390:
	call	___FHranger@0
$L98391:
	shl	ebp, 8
	cmp	edi, 256				; 00000100H
	jg	SHORT $L98392
	cmp	edi, 1
	jge	SHORT $L98393
$L98392:
	call	___FHranger@0
$L98393:
	shl	ebx, 8
	add	ebp, edi
	dec	esi
	inc	edi
	fld	QWORD PTR _SS[ebp*8-2056]
	call	__FIsqrt
	mov	DWORD PTR $T98373[esp+184], 164		; 000000a4H
	add	ebx, DWORD PTR -140+[esp+172]
	test	esi, esi
	fstp	QWORD PTR _SHALF[ebx*8-2056]

; 164  :           DO 3140  II=1,M

	jne	$L98316
	mov	DWORD PTR $T98373[esp+184], 169		; 000000a9H
$L73:

; 168  : C              SHALF(II,II)=(SS(II,II)**(1/2))
; 169  : 3140  CONTINUE

	mov	esi, 256				; 00000100H
	lea	ecx, DWORD PTR $T98318[esp+172]
	mov	DWORD PTR $T98318[esp+172], esi
	push	ecx
	lea	eax, DWORD PTR $T98319[esp+176]
	push	OFFSET FLAT:_TEM
	lea	edx, DWORD PTR $T98320[esp+180]
	push	eax
	mov	ecx, DWORD PTR _M$[esp+180]
	push	OFFSET FLAT:_SHALF
	mov	DWORD PTR $T98319[esp+188], esi
	push	edx
	mov	DWORD PTR $T98320[esp+192], esi
	push	OFFSET FLAT:_U
	mov	DWORD PTR $T98373[esp+208], 171		; 000000abH

; 170  : C
; 171  :       CALL MATRIXM(M,U,256,SHALF,256,TEM,256)

	push	ecx
	call	_MATRIXM@28
	lea	ecx, DWORD PTR $T98321[esp+172]
	mov	DWORD PTR $T98321[esp+172], esi
	lea	edx, DWORD PTR $T98322[esp+172]
	push	ecx
	lea	eax, DWORD PTR $T98323[esp+176]
	push	OFFSET FLAT:_SHALF
	mov	DWORD PTR $T98322[esp+180], esi
	push	edx
	mov	ecx, DWORD PTR _M$[esp+180]
	push	OFFSET FLAT:_UT
	mov	DWORD PTR $T98323[esp+188], esi
	push	eax
	mov	DWORD PTR $T98373[esp+204], 173		; 000000adH

; 172  : C
; 173  :       CALL MATRIXM(M,TEM,256,UT,256,SHALF,256)

	push	OFFSET FLAT:_TEM
	push	ecx
	call	_MATRIXM@28
	mov	ecx, DWORD PTR _M$[esp+168]
	mov	DWORD PTR $T98324[esp+172], esi
	mov	DWORD PTR $T98373[esp+184], 189		; 000000bdH

; 174  : C
; 175  : C
; 176  : C     PRINT-OUT THE S**(1/2) MATRIX:
; 177  : C     Uncomment all code between the (===) lines to execute...
; 178  : C     ===============================================================
; 179  : C      CALL BIOUT(' ')
; 180  : C      CALL BIOUT(' S**(1/2) matrix: ')
; 181  : C      CALL BIOUT(' ---------------- ')
; 182  : C          CALL TAB(SHALF,M,M,100,100)
; 183  : C      CALL BIOUT(' ')
; 184  : C     ===============================================================
; 185  : C      STOP
; 186  : C
; 187  : C
; 188  : C     Blank the values in matrix SNHALF -->0.0D+00
; 189  :       CALL BLANKM(M,SNHALF,256)

	lea	esi, DWORD PTR $T98324[esp+172]
	push	esi
	push	OFFSET FLAT:_SNHALF
	push	ecx
	call	_BLANKM@12
	mov	ecx, 193				; 000000c1H
	mov	eax, DWORD PTR _M$[esp+168]
	mov	esi, 1
	mov	DWORD PTR $T98373[esp+184], ecx

; 190  : C
; 191  : C
; 192  : C     Construct the S**(-1/2) matrix...
; 193  :       DO 3160 II=1,M

	mov	edx, DWORD PTR [eax]
	mov	DWORD PTR -140+[esp+172], edx
	test	edx, edx
	jle	$L77
$L98325:
	mov	DWORD PTR $T98373[esp+184], 194		; 000000c2H

; 194  :           SNHALF(II,II)=DSQRT(SS(II,II))

	cmp	esi, 256				; 00000100H
	jg	SHORT $L98394
	cmp	esi, 1
	jge	SHORT $L98395
$L98394:
	call	___FHranger@0
$L98395:
	mov	ebx, esi
	cmp	esi, 256				; 00000100H
	jg	SHORT $L98396
	cmp	esi, 1
	jge	SHORT $L98397
$L98396:
	call	___FHranger@0
$L98397:
	cmp	esi, 256				; 00000100H
	mov	edi, esi
	mov	ebp, esi
	jg	SHORT $L98398
	cmp	ebp, 1
	jge	SHORT $L98399
$L98398:
	call	___FHranger@0
$L98399:
	shl	ebp, 8
	cmp	esi, 256				; 00000100H
	jg	SHORT $L98400
	cmp	esi, 1
	jge	SHORT $L98401
$L98400:
	call	___FHranger@0
$L98401:
	shl	edi, 8
	add	ebp, esi
	add	edi, ebx
	mov	ebx, esi
	fld	QWORD PTR _SS[ebp*8-2056]
	call	__FIsqrt
	fstp	QWORD PTR _SNHALF[edi*8-2056]
	mov	DWORD PTR $T98373[esp+184], 196		; 000000c4H

; 195  : C          SNHALF(II,II)=(SS(II,II)**(1/2))
; 196  :           SNHALF(II,II)=1.0D+00/SNHALF(II,II)

	cmp	ebx, 256				; 00000100H
	jg	SHORT $L98402
	cmp	ebx, 1
	jge	SHORT $L98403
$L98402:
	call	___FHranger@0
$L98403:
	cmp	esi, 256				; 00000100H
	mov	edi, esi
	jg	SHORT $L98404
	cmp	edi, 1
	jge	SHORT $L98405
$L98404:
	call	___FHranger@0
$L98405:
	fld	QWORD PTR $T98406
	cmp	esi, 256				; 00000100H
	mov	ebp, esi
	jg	SHORT $L98407
	cmp	ebp, 1
	jge	SHORT $L98408
$L98407:
	call	___FHranger@0
$L98408:
	shl	ebp, 8
	cmp	esi, 256				; 00000100H
	jg	SHORT $L98409
	cmp	esi, 1
	jge	SHORT $L98410
$L98409:
	call	___FHranger@0
$L98410:
	shl	edi, 8
	add	ebp, esi
	add	edi, ebx
	inc	esi
	fdiv	QWORD PTR _SNHALF[ebp*8-2056]
	dec	DWORD PTR -140+[esp+172]
	mov	DWORD PTR $T98373[esp+184], 193		; 000000c1H
	cmp	DWORD PTR -140+[esp+172], 0
	fstp	QWORD PTR _SNHALF[edi*8-2056]

; 190  : C
; 191  : C
; 192  : C     Construct the S**(-1/2) matrix...
; 193  :       DO 3160 II=1,M

	jne	$L98325
	mov	DWORD PTR $T98373[esp+184], 197		; 000000c5H
$L77:

; 197  : 3160  CONTINUE

	mov	esi, 256				; 00000100H
	lea	ecx, DWORD PTR $T98326[esp+172]
	mov	DWORD PTR $T98326[esp+172], esi
	push	ecx
	lea	eax, DWORD PTR $T98327[esp+176]
	push	OFFSET FLAT:_TEM
	lea	edx, DWORD PTR $T98328[esp+180]
	push	eax
	mov	ecx, DWORD PTR _M$[esp+180]
	push	OFFSET FLAT:_SNHALF
	mov	DWORD PTR $T98327[esp+188], esi
	push	edx
	mov	DWORD PTR $T98328[esp+192], esi
	push	OFFSET FLAT:_U
	mov	DWORD PTR $T98373[esp+208], 199		; 000000c7H

; 198  : C
; 199  :       CALL MATRIXM(M,U,256,SNHALF,256,TEM,256)

	push	ecx
	call	_MATRIXM@28
	lea	ecx, DWORD PTR $T98329[esp+172]
	mov	DWORD PTR $T98329[esp+172], esi
	lea	edx, DWORD PTR $T98330[esp+172]
	push	ecx
	lea	eax, DWORD PTR $T98331[esp+176]
	push	OFFSET FLAT:_SNHALF
	mov	DWORD PTR $T98330[esp+180], esi
	push	edx
	mov	ecx, DWORD PTR _M$[esp+180]
	push	OFFSET FLAT:_UT
	mov	DWORD PTR $T98331[esp+188], esi
	push	eax
	mov	DWORD PTR $T98373[esp+204], 202		; 000000caH

; 200  : C
; 201  : C
; 202  :       CALL MATRIXM(M,TEM,256,UT,256,SNHALF,256)

	push	OFFSET FLAT:_TEM
	push	ecx
	call	_MATRIXM@28
	lea	ecx, DWORD PTR $T98332[esp+172]
	mov	DWORD PTR $T98332[esp+172], esi
	mov	edx, DWORD PTR _H$[esp+168]
	push	ecx
	lea	eax, DWORD PTR $T98333[esp+176]
	push	OFFSET FLAT:_TEM
	mov	ecx, DWORD PTR _M$[esp+176]
	push	OFFSET FLAT:_MMAX
	mov	DWORD PTR $T98333[esp+184], esi
	push	edx
	mov	DWORD PTR $T98373[esp+200], 232		; 000000e8H

; 203  : C
; 204  : C
; 205  : C     PRINT-OUT THE S**(-1/2) MATRIX:
; 206  : C     Uncomment all code between the (===) lines to execute...
; 207  : C     ===============================================================
; 208  : C      CALL BIOUT(' ')
; 209  : C      CALL BIOUT(' S**(-1/2) matrix: ')
; 210  : C      CALL BIOUT(' ---------------- ')
; 211  : C          CALL TAB(SNHALF,M,M,100,100)
; 212  : C      CALL BIOUT(' ')
; 213  : C     ===============================================================
; 214  : C      STOP
; 215  : C
; 216  : C
; 217  : C     CHECK FOR IDENTITY MATRIX:
; 218  : C     Uncomment all code between the (===) lines to execute...
; 219  : C     ===============================================================
; 220  : C      CALL MATRIXM(M,SNHALF,100,SHALF,100,CHECKS,100)
; 221  : C      CALL BIOUT(' ')
; 222  : C      CALL BIOUT(' S**(-1/2) * S**(1/2) = I CHECK:')
; 223  : C      CALL BIOUT(' -------------------------------')
; 224  : C      CALL TAB(CHECKS,M,M,100,100)
; 225  : C     ===============================================================
; 226  : C      STOP 
; 227  : C
; 228  : C
; 229  : C     CREATE A NEW MATRIX HH=S**(-1/2) * H * S**(-1/2)
; 230  : C
; 231  : C
; 232  :       CALL MATRIXM(M,SNHALF,256,H,MMAX,TEM,256)

	push	eax
	push	OFFSET FLAT:_SNHALF
	push	ecx
	call	_MATRIXM@28
	lea	ecx, DWORD PTR $T98334[esp+172]
	mov	DWORD PTR $T98334[esp+172], esi
	lea	edx, DWORD PTR $T98335[esp+172]
	push	ecx
	lea	eax, DWORD PTR $T98336[esp+176]
	push	OFFSET FLAT:_HH
	mov	DWORD PTR $T98335[esp+180], esi
	push	edx
	mov	ecx, DWORD PTR _M$[esp+180]
	push	OFFSET FLAT:_SNHALF
	mov	DWORD PTR $T98336[esp+188], esi
	push	eax
	mov	DWORD PTR $T98373[esp+204], 235		; 000000ebH

; 233  : C
; 234  : C
; 235  :       CALL MATRIXM(M,TEM,256,SNHALF,256,HH,256)

	push	OFFSET FLAT:_TEM
	push	ecx
	call	_MATRIXM@28
	lea	ecx, DWORD PTR $T98337[esp+172]
	push	OFFSET FLAT:_AA
	mov	edx, DWORD PTR _M$[esp+172]
	push	ecx
	mov	DWORD PTR $T98337[esp+180], esi
	push	OFFSET FLAT:_HH
	mov	DWORD PTR $T98373[esp+196], 254		; 000000feH

; 236  : C
; 237  : C
; 238  : C     PRINT-OUT THE HH MATRIX:
; 239  : C     Uncomment all code between the (===) lines to execute...
; 240  : C     ===============================================================
; 241  : C      CALL BIOUT(' ')
; 242  : C      CALL BIOUT(' HH matrix: ')
; 243  : C      CALL BIOUT(' ---------- ')
; 244  : C          CALL TAB(HH,M,M,100,100)
; 245  : C      CALL BIOUT(' ')
; 246  : C     ===============================================================
; 247  : C      STOP
; 248  : C
; 249  : C
; 250  : C     DIAGNOLIZE HH MATRIX TO GET EIGENVALUE W AND EIGENVECTOR C-----
; 251  : C
; 252  : C
; 253  : C     Define the one dimensional array AA(II) as the HH matrix for YABOCI...
; 254  :       CALL VECH (M,HH,256,AA)

	push	edx
	call	_VECH@16
	mov	ecx, DWORD PTR _M$[esp+168]
	push	OFFSET FLAT:_JB
	mov	DWORD PTR $T98373[esp+188], 258		; 00000102H

; 255  : C
; 256  : C
; 257  : C     Call YACOBI to get the TRANSPOSE of the ROTATION MATRIX...
; 258  :       CALL YACOBI(AA,TP,M,NNP1D2,BIG,JB)

	push	OFFSET FLAT:_BIG
	push	OFFSET FLAT:_NNP1D2
	push	ecx
	push	OFFSET FLAT:_TP
	push	OFFSET FLAT:_AA
	call	_YACOBI@24
	lea	ecx, DWORD PTR $T98338[esp+172]
	mov	edx, DWORD PTR _M$[esp+168]
	mov	DWORD PTR $T98338[esp+172], esi
	push	ecx
	mov	DWORD PTR $T98373[esp+188], 262		; 00000106H

; 259  : C
; 260  : C
; 261  : C     Construct the transpose matrix XT FROM the vector TP...
; 262  :       CALL DEVEC (M,TP,XT,256)

	push	OFFSET FLAT:_XT
	push	OFFSET FLAT:_TP
	push	edx
	call	_DEVEC@16
	lea	ecx, DWORD PTR $T98339[esp+172]
	lea	edx, DWORD PTR $T98340[esp+172]
	mov	DWORD PTR $T98339[esp+172], esi
	push	ecx
	mov	DWORD PTR $T98340[esp+176], esi
	push	OFFSET FLAT:_XY
	mov	DWORD PTR $T98373[esp+192], 277		; 00000115H

; 263  : C
; 264  : C
; 265  : C     PRINT-OUT THE TRANSPOSE ROTATION MATRIX:
; 266  : C     Uncomment all code between the (===) lines to execute...
; 267  : C     ===============================================================
; 268  : C      CALL BIOUT(' ')
; 269  : C      CALL BIOUT(' TRANSPOSE matrix ')
; 270  : C      CALL BIOUT(' ---------------------- ')
; 271  : C         CALL TAB(XT,M,M,100,100)
; 272  : C     ===============================================================
; 273  : C      STOP
; 274  : C
; 275  : C
; 276  : C     Construct the regular rotation matrix from its transpose (XT)...
; 277  :       CALL TRANSPOSE(M,XT,256,XY,256)

	push	edx
	push	OFFSET FLAT:_XT
	mov	eax, DWORD PTR _M$[esp+184]
	push	eax
	call	_TRANSPOSE@20
	lea	ecx, DWORD PTR $T98341[esp+172]
	mov	DWORD PTR $T98341[esp+172], esi
	lea	edx, DWORD PTR $T98342[esp+172]
	push	ecx
	lea	eax, DWORD PTR $T98343[esp+176]
	push	OFFSET FLAT:_TEM
	mov	DWORD PTR $T98342[esp+180], esi
	push	edx
	mov	ecx, DWORD PTR _M$[esp+180]
	push	OFFSET FLAT:_HH
	mov	DWORD PTR $T98343[esp+188], esi
	push	eax
	mov	DWORD PTR $T98373[esp+204], 309		; 00000135H

; 278  : C
; 279  : C
; 280  : C     SHOW THAT THE ROTATION IS A UNITARY TRANSFORMATION:
; 281  : C     Uncomment all code between the (===) lines to execute...
; 282  : C     ===============================================================
; 283  : C      CALL BIOUT(' ')
; 284  : C      CALL BIOUT(' Unitary Transformation?: ')
; 285  : C      CALL BIOUT(' ------------------------ ' )
; 286  : C          CALL MATRIXM (M,XY,100,XT,100,TEM3,100)
; 287  : C          CALL TAB(TEM3,M,M,100,100)
; 288  : C      CALL BIOUT(' ')
; 289  : C     ===============================================================
; 290  : C      STOP
; 291  : C
; 292  : C
; 293  : C     CHECK FOR IDENTITY MATRIX:
; 294  : C     Uncomment all code between the (===) lines to execute...
; 295  : C     ===============================================================
; 296  : C      CALL MATRIXM(M,SNHALF,100,XY,100,TEM,100)
; 297  : C      CALL TRANSPOSE(M,TEM,100,TEM1,100)
; 298  : C      CALL MATRIXM(M,TEM1,100,S,MMAX,TEM2,100)
; 299  : C      CALL MATRIXM(M,TEM2,100,TEM,100,TEM3,100)
; 300  : C      CALL BIOUT(' ')
; 301  : C      CALL BIOUT(' (UT * ST**(-1/2)) * S * (S**(-1/2) * U) = I CHECK:')
; 302  : C      CALL BIOUT(' -------------------------------------------------')
; 303  : C      CALL TAB(TEM3,M,M,100,100)
; 304  : C     ===============================================================
; 305  : C      STOP 
; 306  : C
; 307  : C
; 308  : C     Perform the operations to diagonalize HH
; 309  :       CALL MATRIXM(M,XT,256,HH,256,TEM,256)

	push	OFFSET FLAT:_XT
	push	ecx
	call	_MATRIXM@28
	lea	ecx, DWORD PTR $T98344[esp+172]
	mov	DWORD PTR $T98344[esp+172], esi
	lea	edx, DWORD PTR $T98345[esp+172]
	push	ecx
	mov	DWORD PTR $T98345[esp+176], esi
	push	OFFSET FLAT:_W
	mov	DWORD PTR $T98346[esp+180], esi
	push	edx
	lea	esi, DWORD PTR $T98346[esp+184]
	push	OFFSET FLAT:_XY
	mov	eax, DWORD PTR _M$[esp+184]
	push	esi
	mov	DWORD PTR $T98373[esp+204], 310		; 00000136H

; 310  :       CALL MATRIXM(M,TEM,256,XY,256,W,256)

	push	OFFSET FLAT:_TEM
	push	eax
	call	_MATRIXM@28
	mov	esi, DWORD PTR _EMIN$[esp+168]
	mov	ecx, DWORD PTR _W+4
	mov	DWORD PTR $T98373[esp+184], 326		; 00000146H

; 311  : C
; 312  : C
; 313  : C     PRINT-OUT THE DIAGONALIZED W MATRIX:
; 314  : C     Uncomment all code between the (===) lines to execute...
; 315  : C     ===============================================================
; 316  : C      CALL BIOUT(' ')
; 317  : C      CALL BIOUT(' Diagonalized (W) matrix: ')
; 318  : C      CALL BIOUT(' ----------------------- ' )
; 319  : C          CALL TAB(W,M,M,100,100)
; 320  : C      CALL BIOUT(' ')
; 321  : C     ===============================================================
; 322  : C      STOP
; 323  : C
; 324  : C
; 325  : C     Pull out the Ground State Energy from the Eigenvalues...
; 326  :       EMIN=W(1,1)

	mov	edi, 1
	mov	DWORD PTR [esi+4], ecx
	mov	eax, DWORD PTR _W
	mov	ecx, 329				; 00000149H
	mov	DWORD PTR [esi], eax
	mov	DWORD PTR _IX, edi
	mov	DWORD PTR $T98373[esp+184], 327		; 00000147H

; 327  :       IX=1

	mov	eax, DWORD PTR _M$[esp+168]
	mov	DWORD PTR $T98373[esp+184], ecx

; 328  : C
; 329  :       DO 3211 I=1,M

	mov	ebx, DWORD PTR [eax]
	test	ebx, ebx
	jle	$L84
$L98347:
	cmp	edi, 256				; 00000100H
	mov	ebp, edi
	mov	DWORD PTR $T98373[esp+184], 331		; 0000014bH

; 330  :           CC=0.0D+00
; 331  :           CC=W(I,I)

	jg	SHORT $L98411
	cmp	ebp, 1
	jge	SHORT $L98412
$L98411:
	call	___FHranger@0
$L98412:
	shl	ebp, 8
	cmp	edi, 256				; 00000100H
	jg	SHORT $L98413
	cmp	edi, 1
	jge	SHORT $L98414
$L98413:
	call	___FHranger@0
$L98414:
	fld	QWORD PTR [esi]
	add	ebp, edi
	mov	DWORD PTR $T98373[esp+184], 332		; 0000014cH
	mov	eax, DWORD PTR _W[ebp*8-2052]
	mov	DWORD PTR _CC+4, eax
	mov	eax, DWORD PTR _W[ebp*8-2056]
	mov	DWORD PTR _CC, eax

; 332  :           IF(CC .LT. EMIN) THEN

	fcomp	QWORD PTR _CC
	fwait
	fnstsw	ax
	test	ah, 65					; 00000041H
	jne	SHORT $L85
	mov	ecx, DWORD PTR _CC+4
	mov	edx, DWORD PTR _CC
	mov	DWORD PTR [esi+4], ecx
	mov	DWORD PTR _IX, edi
	mov	DWORD PTR $T98373[esp+184], 333		; 0000014dH

; 333  :               EMIN=CC

	mov	DWORD PTR $T98373[esp+184], 334		; 0000014eH
	mov	DWORD PTR $T98373[esp+184], 329		; 00000149H
	mov	DWORD PTR [esi], edx
$L85:

; 328  : C
; 329  :       DO 3211 I=1,M

	dec	ebx
	inc	edi
	mov	DWORD PTR $T98373[esp+184], 329		; 00000149H
	test	ebx, ebx
	jne	$L98347
	mov	DWORD PTR $T98373[esp+184], 336		; 00000150H
$L84:

; 334  :               IX=I
; 335  :           ENDIF
; 336  : 3211  CONTINUE

	mov	eax, 256				; 00000100H
	lea	ecx, DWORD PTR $T98348[esp+172]
	mov	DWORD PTR $T98348[esp+172], eax
	push	ecx
	lea	edx, DWORD PTR $T98349[esp+176]
	push	OFFSET FLAT:_C
	mov	DWORD PTR $T98349[esp+180], eax
	push	edx
	mov	ecx, DWORD PTR _M$[esp+180]
	push	OFFSET FLAT:_XY
	mov	DWORD PTR $T98373[esp+200], 364		; 0000016cH

; 337  : C
; 338  : C
; 339  : C     Print out the lowest eigenvalue = Ground State Energy
; 340  : C      CALL BIOUT(' ')
; 341  : C      WRITE(7,1856) EMIN
; 342  : C      WRITE(6,1856) EMIN
; 343  : C      RJUICE=EMIN
; 344  : C
; 345  : C
; 346  : C      STOP 
; 347  : C
; 348  : 
; 349  : 
; 350  : 
; 351  : 
; 352  : 
; 353  : 
; 354  : 
; 355  : 
; 356  : 
; 357  : 
; 358  : 
; 359  : 
; 360  : 
; 361  : 
; 362  : 
; 363  : C ***   FIND EIGENVECTOR C
; 364  :       CALL MATRIXM(M,SNHALF,256,XY,256,C,256)

	mov	DWORD PTR $T98350[esp+188], eax
	lea	eax, DWORD PTR $T98350[esp+188]
	push	eax
	push	OFFSET FLAT:_SNHALF
	push	ecx
	call	_MATRIXM@28
	mov	ecx, 365				; 0000016dH
	mov	edx, DWORD PTR _M$[esp+168]
	mov	esi, 1
	mov	DWORD PTR $T98373[esp+184], ecx

; 365  :       DO 600  I=1,M

	mov	eax, DWORD PTR [edx]
	mov	DWORD PTR -140+[esp+172], eax
	test	eax, eax
	jle	$L90
$L98351:
	mov	ecx, 366				; 0000016eH
	mov	eax, DWORD PTR _M$[esp+168]
	mov	edi, 1
	mov	DWORD PTR $T98373[esp+184], ecx

; 366  :       DO 600  J=1,M

	mov	eax, DWORD PTR [eax]
	mov	DWORD PTR -144+[esp+172], eax
	test	eax, eax
	jle	$L91
$L98352:
	cmp	esi, 256				; 00000100H
	mov	ebp, esi
	mov	DWORD PTR $T98373[esp+184], 367		; 0000016fH

; 367  :       C(I,J)=0.0D+00

	jg	SHORT $L98415
	cmp	ebp, 1
	jge	SHORT $L98416
$L98415:
	call	___FHranger@0
$L98416:
	cmp	edi, 256				; 00000100H
	mov	ebx, edi
	jg	SHORT $L98417
	cmp	ebx, 1
	jge	SHORT $L98418
$L98417:
	call	___FHranger@0
$L98418:
	shl	ebx, 8
	xor	ecx, ecx
	add	ebx, ebp
	mov	edx, 368				; 00000170H
	mov	eax, DWORD PTR _M$[esp+168]
	mov	DWORD PTR $T98373[esp+184], edx
	mov	DWORD PTR _C[ebx*8-2056], ecx
	mov	DWORD PTR _C[ebx*8-2052], ecx

; 368  :       DO 650  K=1,M

	mov	ebx, 1
	mov	eax, DWORD PTR [eax]
	mov	DWORD PTR -156+[esp+172], eax
	cmp	eax, ecx
	jle	$L93
$L98353:
	mov	DWORD PTR $T98373[esp+184], 369		; 00000171H

; 369  :       C(I,J)=C(I,J)+SNHALF(I,K)*XY(K,J)

	cmp	esi, 256				; 00000100H
	jg	SHORT $L98419
	cmp	esi, 1
	jge	SHORT $L98420
$L98419:
	call	___FHranger@0
$L98420:
	mov	DWORD PTR -152+[esp+172], esi
	cmp	edi, 256				; 00000100H
	jg	SHORT $L98421
	cmp	edi, 1
	jge	SHORT $L98422
$L98421:
	call	___FHranger@0
$L98422:
	cmp	ebx, 256				; 00000100H
	mov	ebp, ebx
	mov	DWORD PTR -148+[esp+172], edi
	jg	SHORT $L98423
	cmp	ebp, 1
	jge	SHORT $L98424
$L98423:
	call	___FHranger@0
$L98424:
	shl	ebp, 8
	cmp	esi, 256				; 00000100H
	jg	SHORT $L98425
	cmp	esi, 1
	jge	SHORT $L98426
$L98425:
	call	___FHranger@0
$L98426:
	add	ebp, esi
	fld	QWORD PTR _SNHALF[ebp*8-2056]
	cmp	edi, 256				; 00000100H
	mov	ebp, edi
	jg	SHORT $L98427
	cmp	ebp, 1
	jge	SHORT $L98428
$L98427:
	call	___FHranger@0
$L98428:
	shl	ebp, 8
	cmp	ebx, 256				; 00000100H
	jg	SHORT $L98429
	cmp	ebx, 1
	jge	SHORT $L98430
$L98429:
	call	___FHranger@0
$L98430:
	add	ebp, ebx
	fmul	QWORD PTR _XY[ebp*8-2056]
	cmp	edi, 256				; 00000100H
	mov	ebp, edi
	jg	SHORT $L98431
	cmp	ebp, 1
	jge	SHORT $L98432
$L98431:
	call	___FHranger@0
$L98432:
	shl	ebp, 8
	cmp	esi, 256				; 00000100H
	jg	SHORT $L98433
	cmp	esi, 1
	jge	SHORT $L98434
$L98433:
	call	___FHranger@0
$L98434:
	add	ebp, esi
	mov	eax, DWORD PTR -148+[esp+172]
	shl	eax, 8
	inc	ebx
	fadd	QWORD PTR _C[ebp*8-2056]
	add	eax, DWORD PTR -152+[esp+172]
	dec	DWORD PTR -156+[esp+172]
	mov	DWORD PTR $T98373[esp+184], 368		; 00000170H
	fstp	QWORD PTR _C[eax*8-2056]

; 368  :       DO 650  K=1,M

	cmp	DWORD PTR -156+[esp+172], 0
	jne	$L98353
	mov	DWORD PTR $T98373[esp+184], 366		; 0000016eH
$L93:

; 366  :       DO 600  J=1,M

	mov	eax, DWORD PTR -144+[esp+172]
	inc	edi
	mov	DWORD PTR $T98373[esp+184], 366		; 0000016eH
	dec	eax
	mov	DWORD PTR -144+[esp+172], eax
	test	eax, eax
	jne	$L98352
	mov	DWORD PTR $T98373[esp+184], 365		; 0000016dH
$L91:

; 365  :       DO 600  I=1,M

	mov	eax, DWORD PTR -140+[esp+172]
	inc	esi
	mov	DWORD PTR $T98373[esp+184], 365		; 0000016dH
	dec	eax
	mov	DWORD PTR -140+[esp+172], eax
	test	eax, eax
	jne	$L98351
	mov	DWORD PTR $T98373[esp+184], 371		; 00000173H
$L90:

; 370  :  650   CONTINUE
; 371  :  600   CONTINUE

	mov	eax, DWORD PTR _M$[esp+168]
	mov	esi, 1
	mov	DWORD PTR $T98373[esp+184], 371		; 00000173H
	mov	ecx, DWORD PTR [eax]
	mov	eax, 380				; 0000017cH
	mov	edi, ecx
	mov	DWORD PTR -140+[esp+172], ecx
	mov	DWORD PTR $T98373[esp+184], eax
	test	edi, edi

; 372  : C      open(unit=6,access='append',form='formatted')
; 373  : C      CALL BIOUT('CHECK THE EIGENVECTOR C')
; 374  : C      DO 220  I=1,M
; 375  : C      DO 220  J=1,M
; 376  : C      WRITE(6,*) I,J,C(I,J)
; 377  : C      WRITE(7,*) I,J,C(I,J)
; 378  : C 220   CONTINUE
; 379  : C     CLOSE(6)
; 380  :       DO 226  I=1,M

	jle	SHORT $L50
$L98354:
	cmp	esi, DWORD PTR -140+[esp+172]
	mov	ebx, esi
	mov	DWORD PTR $T98373[esp+184], 381		; 0000017dH

; 381  :       CO(I) = C(I,IX)

	jg	SHORT $L98435
	cmp	ebx, 1
	jge	SHORT $L98436
$L98435:
	call	___FHranger@0
$L98436:
	mov	ebp, DWORD PTR _IX
	cmp	ebp, 256				; 00000100H
	jg	SHORT $L98437
	cmp	ebp, 1
	jge	SHORT $L98438
$L98437:
	call	___FHranger@0
$L98438:
	shl	ebp, 8
	cmp	esi, 256				; 00000100H
	jg	SHORT $L98439
	cmp	esi, 1
	jge	SHORT $L98440
$L98439:
	call	___FHranger@0
$L98440:
	add	ebp, esi
	mov	ecx, DWORD PTR _CO$[esp+168]
	dec	edi
	inc	esi
	mov	eax, DWORD PTR _C[ebp*8-2052]
	test	edi, edi
	mov	DWORD PTR $T98373[esp+184], 380		; 0000017cH
	mov	DWORD PTR [ecx+ebx*8-4], eax
	mov	eax, DWORD PTR _C[ebp*8-2056]
	mov	DWORD PTR [ecx+ebx*8-8], eax

; 372  : C      open(unit=6,access='append',form='formatted')
; 373  : C      CALL BIOUT('CHECK THE EIGENVECTOR C')
; 374  : C      DO 220  I=1,M
; 375  : C      DO 220  J=1,M
; 376  : C      WRITE(6,*) I,J,C(I,J)
; 377  : C      WRITE(7,*) I,J,C(I,J)
; 378  : C 220   CONTINUE
; 379  : C     CLOSE(6)
; 380  :       DO 226  I=1,M

	jne	SHORT $L98354
	mov	DWORD PTR $T98373[esp+184], 382		; 0000017eH
$L50:

; 382  :  226   CONTINUE

	push	0
	mov	ecx, DWORD PTR _SNORM
	mov	DWORD PTR $T98373[esp+188], 436		; 000001b4H

; 383  : C      DO 230 I=1,M
; 384  : C      DO 230 J=1,M
; 385  : C      CT(I,J)=C(J,I)
; 386  : C 230   CONTINUE
; 387  : C      CALL MATRIXM(M,CT,100,S,100,TEM,100)
; 388  : C      DO 700  I=1,M
; 389  : C      DO 700  J=1,M
; 390  : C      TEM(I,J)=0
; 391  : C      DO 750  K=1,M
; 392  : C      TEM(I,J)=TEM(I,J)+CT(I,K)*S(K,J)
; 393  : C 750   CONTINUE
; 394  : C 700   CONTINUE
; 395  : C      CALL MATRIXM(M,TEM,100,C,100,CHECKS,100)
; 396  : C      DO 100  I=1,M
; 397  : C      DO 100  J=1,M
; 398  : C      CHECKS(I,J)=0
; 399  : C      DO 150  K=1,M
; 400  : C      CHECKS(I,J)=CHECKS(I,J)+TEM(I,K)*C(K,J)
; 401  : C 150   CONTINUE
; 402  : C 100   CONTINUE
; 403  : C      CALL TAB(CHECKS,M,M,100,100)
; 404  : C     STOP
; 405  : C       open(unit=6,access='append',form='formatted')
; 406  : C     WRITE(6,*) 'CHECK S'
; 407  : C     DO 240 I=1,M
; 408  : C     DO 240 J=1,M
; 409  : C     WRITE(6,*) CHECKS(I,J)
; 410  : C240  CONTINUE
; 411  : C     close(6)
; 412  : C      DO 90 L=1,M
; 413  : C      DO 90 K=1,M
; 414  : C      SSS = SSS + CO(L)*CO(K)*S(L,K)
; 415  : C90    CONTINUE 
; 416  : C *** 
; 417  : C      open (unit=6,access='append',form='formatted')
; 418  : C      WRITE(6,*)' NORMALIZATION = ', SSS
; 419  : C      WRITE(7,*)' NORMALIZATION = ', SSS
; 420  : C      CALL TAB (C,M,M,100,100)
; 421  : C      CALL BIOUT('******************************************')
; 422  : C      CALL BIOUT('          EXPANSION COEFFICIENTS')
; 423  : C      CALL BIOUT('==========================================')
; 424  : C      DO 33 J = 1, M
; 425  : 	
; 426  : C      WRITE(6,698) J, CO(J)
; 427  : C      WRITE(7,698) J, CO(J)
; 428  : C698   FORMAT(10X, 'CO(', I2,') =' D17.10)
; 429  : C33    CONTINUE 
; 430  : C      GO TO 1515
; 431  :  1414  CONTINUE
; 432  : C      EMIN = 10.0D+00
; 433  : C      open (unit=6,access='append',form='formatted')
; 434  : C      write (7,*) ' emin = 10.00 '
; 435  : C      close (6)
; 436  : 1515  DEALLOCATE(SNORM)

	push	8
	push	OFFSET FLAT:_7SNORM
	push	ecx
	call	__Deallocate_i
	mov	DWORD PTR $T98373[esp+200], 445		; 000001bdH
	mov	DWORD PTR _SNORM, 0
	add	esp, 16					; 00000010H

; 437  : 	RETURN
; 438  : C
; 439  : C     FORMATS:
; 440  : C     --------------------------------------------
; 441  :  11   FORMAT(/,/)
; 442  : 1855  FORMAT (/)
; 443  : 1856  FORMAT (' GroundState Energy (CALCULATED) = ',F18.12,' au')
; 444  : C     --------------------------------------------
; 445  :       END

	call	___FHunlink@0
	pop	ebp
	pop	edi
	pop	esi
	pop	ebx
	add	esp, 156				; 0000009cH
	ret	24					; 00000018H
	npad	7
_ENERGY@24 ENDP
_TEXT	ENDS
END
