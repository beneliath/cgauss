	TITLE	C:\users\dgilmore\Graduate.School\Adamowicz.Research\H2ANL\ENUCL.FOR
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
PUBLIC	_DOVN
PUBLIC	_DEN
EXTRN	__FFljj:NEAR
EXTRN	_f90_write:NEAR
EXTRN	_lwrite_ch:NEAR
EXTRN	_lwrite_d:NEAR
EXTRN	_f90_endwrite:NEAR
EXTRN	_F0@4:NEAR
_DATA	SEGMENT
_DOVN	DB	028H DUP (?)
_DEN	DB	050H DUP (?)
	ORG $-120
_C:\users\dgilmore\Graduate.School\Adamowicz.Research\H2ANL\ENUCL.FOR DB 
	DB	'C:\users\dgilmore\Graduate.School\Adamowicz.Research\H2ANL\EN'
	DB	'UCL.FOR', 00H
_DATA	ENDS
PUBLIC	_ENUCL@100
EXTRN	__FIsqrt:NEAR
EXTRN	__FIfexp:NEAR
EXTRN	___FHlink@4:NEAR
EXTRN	___FHunlink@0:NEAR
EXTRN	__fltused:NEAR
_BSS	SEGMENT
_DPKL	DQ	01H DUP (?)
_ALK1	DQ	01H DUP (?)
_ALK2	DQ	01H DUP (?)
_BLK12	DQ	01H DUP (?)
_ALKXX	DQ	01H DUP (?)
_DETAB	DQ	01H DUP (?)
_RZ1	DQ	01H DUP (?)
_RZ2	DQ	01H DUP (?)
_RRZ	DQ	01H DUP (?)
_RE2	DQ	01H DUP (?)
_DELTA2	DQ	01H DUP (?)
_Z	DQ	01H DUP (?)
_PHIF	DQ	01H DUP (?)
_THETAF	DQ	01H DUP (?)
_OKL	DQ	01H DUP (?)
_QKL	DQ	01H DUP (?)
_RRR	DQ	01H DUP (?)
_RRRN	DQ	01H DUP (?)
_A1Z	DQ	01H DUP (?)
_A2Z	DQ	01H DUP (?)
_A3Z	DQ	01H DUP (?)
_A4Z	DQ	01H DUP (?)
_AZSUM	DQ	01H DUP (?)
_AZSUM1	DQ	01H DUP (?)
_AZSUM2	DQ	01H DUP (?)
_DDETAB	DQ	01H DUP (?)
_DZZ	DQ	01H DUP (?)
_DTHETAF DQ	01H DUP (?)
_DPHIF	DQ	01H DUP (?)
_DF0	DQ	01H DUP (?)
_FF	DQ	01H DUP (?)
_PKL	DQ	01H DUP (?)
_ANKL	DQ	01H DUP (?)
_AMKL	DQ	01H DUP (?)
_DOKL	DQ	01H DUP (?)
_DNKL	DQ	01H DUP (?)
_DMKL	DQ	01H DUP (?)
_BSS	ENDS
_DATA	SEGMENT
	ORG $+3
_.SG89_	DB	'RE2='
	ORG $+4
$T98309	DB	'_ENUCL@100', 00H
_DATA	ENDS
_CONST	SEGMENT
$T98311	DQ	00000000000000000r		; 0
$T98312	DQ	0400921fb54442d28r		; 3.14159
$T98313	DQ	04000000000000000r		; 2
$T98314	DQ	04010000000000000r		; 4
$T98315	DQ	03f12e1a07588fc71r		; 7.20267e-005
$T98316	DQ	03f4f4f8e903a27b1r		; 0.000955529
$T98317	DQ	03f84c5ef42558089r		; 0.0101432
$T98318	DQ	03fb2e7fbe90c4bfdr		; 0.0738523
$T98319	DQ	03fd5a92bbb78443dr		; 0.33845
$T98320	DQ	03fec287351b9e49br		; 0.879938
$T98321	DQ	03ff0000000000000r		; 1
$T98322	DQ	03f0d9b282dc36e12r		; 5.64691e-005
$T98323	DQ	03f48da32e37f14e9r		; 0.000758433
$T98324	DQ	03f7f88563ebdb857r		; 0.00769838
$T98325	DQ	03fb01c78cb5ca40ar		; 0.0629344
$T98326	DQ	03fcb4c7959a1f6a4r		; 0.213271
$T98327	DQ	03fe66ae7d566cf42r		; 0.70055
$T98328	DQ	03f3280f91c9a24cbr		; 0.000282346
$T98329	DQ	03f68da32e37f14e9r		; 0.00303373
$T98330	DQ	03f97a640af0e4a41r		; 0.0230951
$T98331	DQ	03fc01c78cb5ca40ar		; 0.125869
$T98332	DQ	03f3c5270b04d7aaar		; 0.00043216
$T98333	DQ	03f7391b91a2458cfr		; 0.00477764
$T98334	DQ	03fa4c5ef42558089r		; 0.0405726
$T98335	DQ	03fcc5bf9dd9271fcr		; 0.221557
$T98336	DQ	03fe5a92bbb78443dr		; 0.676901
$T98337	DQ	0bff0000000000000r		; -1
_CONST	ENDS
_TEXT	SEGMENT
; File C:\users\dgilmore\Graduate.School\Adamowicz.Research\H2ANL\ENUCL.FOR
$T98308 = -16
_DFLAG$ = 8
_OVERLAP$ = 12
_AZ$ = 24
_A1$ = 28
_BZ$ = 40
_A2$ = 44
_CZ$ = 56
_A3$ = 60
_DZ$ = 72
_A4$ = 76
_A5$ = 80
_A6$ = 84
_EZ$ = 96
_EN$ = 100
_KLM$ = 104
_ENUCL@100 PROC NEAR

; 6    :       SUBROUTINE ENUCL(DFLAG,OVERLAP,AX,AY,AZ,A1,BX,BY,BZ,A2,                      

	sub	esp, 192				; 000000c0H
	lea	eax, DWORD PTR $T98308[esp+192]
	push	ebx
	push	esi
	push	edi
	push	ebp
	push	eax
	call	___FHlink@4
	mov	ebx, DWORD PTR _A2$[esp+204]
	mov	esi, DWORD PTR _A1$[esp+204]
	mov	DWORD PTR $T98308[esp+212], OFFSET FLAT:_C:\users\dgilmore\Graduate.School\Adamowicz.Research\H2ANL\ENUCL.FOR
	fld	QWORD PTR [ebx]
	fadd	QWORD PTR [esi]
	fstp	QWORD PTR _ALK1
	mov	ecx, DWORD PTR _A4$[esp+204]
	mov	edi, DWORD PTR _A3$[esp+204]
	mov	DWORD PTR $T98308[esp+216], OFFSET FLAT:$T98309
	fld	QWORD PTR [ecx]
	fadd	QWORD PTR [edi]
	fstp	QWORD PTR _ALK2
	mov	ecx, DWORD PTR _A6$[esp+204]
	mov	edx, DWORD PTR _A5$[esp+204]
	mov	DWORD PTR $T98308[esp+220], 46		; 0000002eH

; 7    :      1                 CX,CY,CZ,A3,DX,DY,DZ,A4,A5,A6,
; 8    :      2                 EX,EY,EZ,EN,KLM)                       
; 9    : C     -----------------------------------------------------
; 10   : C     ENUCL INTGL_1 = <psi(L)@A|(-Z_C/|r1 - R_C|)|psi{k)@B>
; 11   : C     -----------------------------------------------------
; 12   :       IMPLICIT REAL*8 (A-H,O-Z)
; 13   : 	INTEGER DFLAG
; 14   : C     -----------------------------------------------------
; 15   : C     *** FOUR CENTER NUCLEAR ATTRACTION INTEGRAL ***
; 16   : C     -----------------------------------------------------
; 17   : 	PARAMETER(C1=0.564691197D-04)
; 18   : 	PARAMETER(C2=0.758433197D-03)
; 19   : 	PARAMETER(C3=0.769838037D-02)
; 20   : 	PARAMETER(C4=0.629344460D-01)
; 21   : 	PARAMETER(C5=0.213271302D+00)
; 22   : 	PARAMETER(C6=0.720266520D-04)
; 23   : 	PARAMETER(C7=0.955528842D-03)
; 24   : 	PARAMETER(C8=0.101431553D-01)
; 25   : 	PARAMETER(C9=0.738522953D-01)
; 26   : 	PARAMETER(C10=0.338450368D+00)
; 27   : 	PARAMETER(C11=0.879937801D+00)
; 28   : 	PARAMETER(ONE=1.00D+00)
; 29   : 	PARAMETER(TWO=2.00D+00)
; 30   : 	PARAMETER(FOUR=4.00D+00)
; 31   : 	PARAMETER(SIX=6.00D+00)
; 32   : 	PARAMETER(EIGHT=8.00D+00)
; 33   : 	PARAMETER(HALF=0.50D+00)
; 34   : 	PARAMETER(PI=3.1415926535898D+00)
; 35   : 	PARAMETER(RX1=0.)
; 36   : 	PARAMETER(RX2=0.)
; 37   : 	PARAMETER(RY1=0.)
; 38   : 	PARAMETER(RY2=0.)
; 39   : 	PARAMETER(RRX=0.)
; 40   : 	PARAMETER(RRY=0.)
; 41   : C
; 42   : 	COMMON /DEN/ DENDA1,DENDA2,DENDA5,DENDAZ,DENDBZ,
; 43   : 	1			 DENDA1P,DENDA2P,DENDA5P,DENDAZP,DENDBZP
; 44   : 	COMMON /DOVn/ DOVnDA1,DOVnDA2,DOVnDA5,DOVnDAZ,DOVnDBZ
; 45   : C     -----------------------------------------------------
; 46   :       ALK1=A1+A2                                                       

	fld	QWORD PTR [ecx]
	fadd	QWORD PTR [edx]
	fstp	QWORD PTR _BLK12
	fld	QWORD PTR _BLK12
	fadd	QWORD PTR _ALK2
	fstp	QWORD PTR -32+[esp+208]
	mov	ecx, DWORD PTR -32+[esp+212]
	mov	edx, DWORD PTR -32+[esp+208]
	fld	QWORD PTR _ALK2
	fadd	QWORD PTR _ALK1
	fmul	QWORD PTR _BLK12
	fld	QWORD PTR _ALK2
	fmul	QWORD PTR _ALK1
	faddp	ST(1), ST(0)
	fstp	QWORD PTR _DETAB
	mov	DWORD PTR $T98308[esp+220], 47		; 0000002fH

; 47   :       ALK2=A3+A4                                                       

	fld	QWORD PTR -32+[esp+208]
	fcomp	QWORD PTR $T98311
	fwait
	fnstsw	ax
	mov	DWORD PTR $T98308[esp+220], 48		; 00000030H

; 48   :       BLK12=A5+A6  

	mov	DWORD PTR $T98308[esp+220], 49		; 00000031H

; 49   :       ALKXX=ALK2+BLK12

	mov	DWORD PTR $T98308[esp+220], 50		; 00000032H
	mov	DWORD PTR $T98308[esp+220], 51		; 00000033H
	test	ah, 65					; 00000041H
	mov	DWORD PTR _ALKXX+4, ecx
	mov	DWORD PTR _ALKXX, edx

; 50   :       DETAB=ALK1*ALK2+(ALK1+ALK2)*BLK12
; 51   :       IF(ALKXX.LE.0.0D+00) GOTO 300

	jne	$L80
	fld	QWORD PTR _DETAB
	fcomp	QWORD PTR $T98311
	fwait
	fnstsw	ax
	mov	DWORD PTR $T98308[esp+220], 52		; 00000034H

; 52   :       IF(DETAB.LE.0.0D+00) GOTO 300                       

	test	ah, 65					; 00000041H
	jne	$L80
	fld	QWORD PTR [esi]
	fld	QWORD PTR [ebx]
	mov	ecx, DWORD PTR _AZ$[esp+204]
	mov	edx, DWORD PTR _BZ$[esp+204]
	mov	DWORD PTR $T98308[esp+220], 56		; 00000038H

; 53   : C     -----------------------------------------------------
; 54   : C     --- ONE CENTER NUCLEAR INTEGRAL ---
; 55   : C     -----------------------------------------------------
; 56   :       RZ1=(A1*AZ+A2*BZ)/ALK1

	fmul	QWORD PTR [edx]
	fxch	ST(1)
	fmul	QWORD PTR [ecx]
	mov	ecx, DWORD PTR _A4$[esp+204]
	mov	edx, DWORD PTR _DZ$[esp+204]
	mov	DWORD PTR $T98308[esp+220], 57		; 00000039H
	faddp	ST(1), ST(0)
	fld	QWORD PTR [ecx]
	fmul	QWORD PTR [edx]
	fxch	ST(1)
	fdiv	QWORD PTR _ALK1
	mov	eax, DWORD PTR _CZ$[esp+204]
	mov	ecx, DWORD PTR -32+[esp+212]
	mov	DWORD PTR -24+[esp+212], ecx
	mov	edx, DWORD PTR -32+[esp+208]
	mov	DWORD PTR -24+[esp+208], edx
	mov	DWORD PTR $T98308[esp+220], 62		; 0000003eH
	fstp	QWORD PTR _RZ1

; 57   :       RZ2=(A3*CZ+A4*DZ)/ALK2

	fld	QWORD PTR [edi]
	fmul	QWORD PTR [eax]
	fld	QWORD PTR -24+[esp+208]
	fmul	QWORD PTR $T98312
	fxch	ST(1)
	faddp	ST(2), ST(0)
	fdivr	QWORD PTR _DETAB
	fxch	ST(1)
	fdiv	QWORD PTR _ALK2
	fstp	QWORD PTR _RZ2

; 58   : C     -----------------------------------------------------
; 59   : C     --- FOUR-CENTER ELECTRON REPULSION INTEGRAL ---
; 60   : C     -----------------------------------------------------
; 61   : C
; 62   : 	EN=TWO*OVERLAP*DSQRT(DETAB/(PI*(ALK2+BLK12)))

	call	__FIsqrt
	fld	QWORD PTR _BLK12
	fmul	QWORD PTR _ALK2
	mov	ebp, DWORD PTR _OVERLAP$[esp+204]
	mov	ecx, DWORD PTR _EN$[esp+204]
	fld	QWORD PTR _ALK1
	fmul	QWORD PTR -24+[esp+208]
	fxch	ST(1)
	fmul	QWORD PTR _RZ2
	fxch	ST(1)
	fmul	QWORD PTR _RZ1
	fxch	ST(2)
	fmul	QWORD PTR [ebp]
	mov	DWORD PTR $T98308[esp+220], 64		; 00000040H
	mov	DWORD PTR $T98308[esp+220], 66		; 00000042H
	fmul	QWORD PTR $T98313
	fstp	QWORD PTR [ecx]

; 63   : C
; 64   :       RRZ=(ALK1*(ALK2+BLK12)*RZ1+ALK2*BLK12*RZ2)/DETAB 

	faddp	ST(1), ST(0)
	mov	ecx, DWORD PTR _EZ$[esp+204]
	fdiv	QWORD PTR _DETAB
	fstp	QWORD PTR _RRZ

; 65   : C
; 66   :       RE2=(RRZ-EZ)**2

	fld	QWORD PTR _RRZ
	fsub	QWORD PTR [ecx]
	fld	QWORD PTR $T98313
	call	__FIfexp
	fstp	QWORD PTR _RE2
	mov	DWORD PTR $T98308[esp+220], 67		; 00000043H

; 67   : 	WRITE(6,*)'RE2=',RE2

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	4
	push	OFFSET FLAT:_.SG89_
	call	_lwrite_ch
	add	esp, 8
	mov	ecx, DWORD PTR _RE2+4
	mov	edx, DWORD PTR _RE2
	npad	1
	push	ecx
	push	edx
	call	_lwrite_d
	add	esp, 8
	call	_f90_endwrite
	fld	QWORD PTR _BLK12
	fadd	QWORD PTR _ALK2
	fld	QWORD PTR _DETAB
	fmul	QWORD PTR $T98314
	mov	DWORD PTR $T98308[esp+220], 70		; 00000046H

; 68   : 
; 69   : C
; 70   :       DELTA2=(ALK2+BLK12)/(FOUR*DETAB)

	mov	DWORD PTR $T98308[esp+220], 71		; 00000047H
	fdivp	ST(1), ST(0)
	mov	DWORD PTR $T98308[esp+220], 76		; 0000004cH
	push	OFFSET FLAT:_RE2
	fstp	QWORD PTR _DELTA2

; 71   :       RE2=RE2/(FOUR*DELTA2)

	fld	QWORD PTR _DELTA2
	fmul	QWORD PTR $T98314
	fdivr	QWORD PTR _RE2
	fstp	QWORD PTR _RE2

; 72   : C
; 73   : 
; 74   : 
; 75   : 
; 76   :       EN=EN*F0(RE2) 

	call	_F0@4
	mov	ecx, DWORD PTR _EN$[esp+204]
	mov	edx, DWORD PTR _DFLAG$[esp+204]
	mov	DWORD PTR $T98308[esp+220], 78		; 0000004eH
	fmul	QWORD PTR [ecx]
	cmp	DWORD PTR [edx], 1
	fstp	QWORD PTR [ecx]

; 77   : C     -----------------------------------------------------
; 78   : 	IF (DFLAG .EQ. 1) THEN

	jne	$L80
	fld	QWORD PTR _RE2
	fmul	QWORD PTR $T98315
	mov	ecx, DWORD PTR _RE2+4
	mov	edx, DWORD PTR _RE2
	mov	DWORD PTR $T98308[esp+220], 80		; 00000050H

; 79   : C
; 80   : 	Z=RE2

	fadd	QWORD PTR $T98316
	mov	DWORD PTR $T98308[esp+220], 81		; 00000051H
	mov	DWORD PTR _Z+4, ecx
	fmul	QWORD PTR _RE2
	mov	DWORD PTR _Z, edx

; 81   : 	PHIF=((((((C6*Z+C7)*Z+C8)*Z+C9)*Z+C10)*Z+C11)*Z+1.00D+00)

	fadd	QWORD PTR $T98317
	fmul	QWORD PTR _RE2
	fadd	QWORD PTR $T98318
	fmul	QWORD PTR _RE2
	mov	DWORD PTR $T98308[esp+220], 82		; 00000052H
	fld	QWORD PTR _RE2
	fmul	QWORD PTR $T98322
	fxch	ST(1)
	fadd	QWORD PTR $T98319
	fxch	ST(1)
	fadd	QWORD PTR $T98323
	fxch	ST(1)
	fmul	QWORD PTR _RE2
	fxch	ST(1)
	fmul	QWORD PTR _RE2
	fxch	ST(1)
	fadd	QWORD PTR $T98320
	fxch	ST(1)
	fadd	QWORD PTR $T98324
	fxch	ST(1)
	fmul	QWORD PTR _RE2
	fxch	ST(1)
	fmul	QWORD PTR _RE2
	fxch	ST(1)
	fadd	QWORD PTR $T98321
	fxch	ST(1)
	fadd	QWORD PTR $T98325
	fld	QWORD PTR _DETAB
	fxch	ST(2)
	fstp	QWORD PTR _PHIF

; 82   : 	THETAF=(((((C1*Z+C2)*Z+C3)*Z+C4)*Z+C5)*Z+1.00D+00)

	fmul	QWORD PTR _RE2
	fadd	QWORD PTR $T98326
	fmul	QWORD PTR _RE2
	fadd	QWORD PTR $T98321
	fstp	QWORD PTR _THETAF
	call	__FIsqrt
	fstp	QWORD PTR -72+[esp+208]
	fld	QWORD PTR _RRZ
	mov	ecx, DWORD PTR -72+[esp+212]
	mov	edx, DWORD PTR -72+[esp+208]
	fld	QWORD PTR $T98313
	mov	DWORD PTR $T98308[esp+220], 84		; 00000054H

; 83   : C
; 84   : 	OKL=DSQRT(DETAB)

	mov	DWORD PTR $T98308[esp+220], 86		; 00000056H
	mov	DWORD PTR _OKL+4, ecx
	mov	ecx, DWORD PTR _EZ$[esp+204]
	fxch	ST(1)
	fsub	QWORD PTR [ecx]
	mov	DWORD PTR _OKL, edx

; 85   : C
; 86   : 	QKL=(RRZ-EZ)**2.00D+00

	fxch	ST(1)
	call	__FIfexp
	fstp	QWORD PTR _QKL
	fld	QWORD PTR _RRZ
	fadd	QWORD PTR $T98327
	fld	QWORD PTR _RRZ
	fsub	QWORD PTR $T98327
	fld	QWORD PTR [esi]
	fxch	ST(2)
	fstp	QWORD PTR _RRR
	fstp	QWORD PTR _RRRN
	mov	DWORD PTR $T98308[esp+220], 88		; 00000058H

; 87   : C
; 88   : 	RRR=0.70055D+00+RRZ

	fld	QWORD PTR [ebx]
	mov	ecx, DWORD PTR _AZ$[esp+204]
	mov	edx, DWORD PTR _BZ$[esp+204]
	mov	DWORD PTR $T98308[esp+220], 89		; 00000059H

; 89   : 	RRRN=-0.70055D+00+RRZ

	fmul	QWORD PTR [edx]
	fld	QWORD PTR [edi]
	fxch	ST(2)
	fmul	QWORD PTR [ecx]
	fxch	ST(1)
	fstp	QWORD PTR _A2Z
	fstp	QWORD PTR _A1Z
	mov	DWORD PTR $T98308[esp+220], 90		; 0000005aH

; 90   : 	A1Z=A1*AZ

	mov	ecx, DWORD PTR _CZ$[esp+204]
	mov	edx, DWORD PTR _A4$[esp+204]
	mov	DWORD PTR $T98308[esp+220], 91		; 0000005bH

; 91   : 	A2Z=A2*BZ

	fmul	QWORD PTR [ecx]
	fld	QWORD PTR [edx]
	mov	DWORD PTR $T98308[esp+220], 92		; 0000005cH

; 92   : 	A3Z=A3*CZ

	fxch	ST(1)
	fstp	QWORD PTR _A3Z
	mov	DWORD PTR $T98308[esp+220], 93		; 0000005dH

; 93   : 	A4Z=A4*DZ

	mov	DWORD PTR $T98308[esp+220], 94		; 0000005eH
	mov	ecx, DWORD PTR _DZ$[esp+204]
	fmul	QWORD PTR [ecx]
	fstp	QWORD PTR _A4Z
	fld	QWORD PTR _A4Z
	fadd	QWORD PTR _A3Z
	fld	QWORD PTR _A2Z
	fld	QWORD PTR _A2Z
	fadd	QWORD PTR _A1Z
	fxch	ST(2)
	fstp	QWORD PTR -24+[esp+208]

; 94   : 	AZSUM=A1Z+A2Z+A3Z+A4Z

	fadd	QWORD PTR -24+[esp+208]
	fxch	ST(1)
	fstp	QWORD PTR _AZSUM1
	fadd	QWORD PTR _A1Z
	mov	ecx, DWORD PTR -24+[esp+212]
	mov	edx, DWORD PTR -24+[esp+208]
	fld	QWORD PTR _ALKXX
	fmul	QWORD PTR _AZSUM1
	fxch	ST(1)
	fstp	QWORD PTR _AZSUM
	fld	QWORD PTR _BLK12
	fmul	QWORD PTR -24+[esp+208]
	fld	QWORD PTR _DETAB
	fxch	ST(1)
	faddp	ST(2), ST(0)
	fld	QWORD PTR $T98313
	mov	DWORD PTR $T98308[esp+220], 95		; 0000005fH

; 95   : 	AZSUM1=A1Z+A2Z

	fxch	ST(2)
	fstp	QWORD PTR -80+[esp+208]
	mov	DWORD PTR $T98308[esp+220], 96		; 00000060H

; 96   : 	AZSUM2=A3Z+A4Z

	mov	DWORD PTR _AZSUM2+4, ecx
	mov	DWORD PTR _AZSUM2, edx
	fxch	ST(1)
	call	__FIfexp
	fstp	QWORD PTR -88+[esp+208]
	mov	ecx, DWORD PTR _AZ$[esp+204]
	fld	QWORD PTR [ecx]
	fmul	QWORD PTR _ALKXX
	fld	QWORD PTR _ALKXX
	fmul	QWORD PTR -80+[esp+208]
	fxch	ST(1)
	fdiv	QWORD PTR _DETAB
	mov	DWORD PTR $T98308[esp+220], 102		; 00000066H
	fxch	ST(1)
	fdiv	QWORD PTR -88+[esp+208]
	fsubp	ST(1), ST(0)
	fstp	QWORD PTR -40+[esp+208]

; 97   : C
; 98   : C	DENDA1:
; 99   : C	======================================================
; 100  : 	DDETAB=ALKXX
; 101  : C
; 102  : 	DZZ=QKL+DETAB*(2.00D+00*(ALKXX*AZ/DETAB-ALKXX*(ALKXX*AZSUM1

	fld	QWORD PTR -40+[esp+208]
	fmul	QWORD PTR _RRRN
	fld	QWORD PTR _RE2
	fmul	QWORD PTR $T98328
	fxch	ST(1)
	fmul	QWORD PTR $T98313
	fxch	ST(1)
	fadd	QWORD PTR $T98329
	fxch	ST(1)
	fmul	QWORD PTR _DETAB
	fxch	ST(1)
	fmul	QWORD PTR _RE2
	fxch	ST(1)
	fdiv	QWORD PTR _ALKXX
	fadd	QWORD PTR _QKL
	fxch	ST(1)
	fadd	QWORD PTR $T98330
	mov	DWORD PTR $T98308[esp+220], 105		; 00000069H
	fmul	QWORD PTR _RE2
	fxch	ST(1)
	fstp	QWORD PTR _DZZ
	fadd	QWORD PTR $T98331
	fld	QWORD PTR _DZZ
	fld	QWORD PTR _RE2
	fmul	QWORD PTR $T98332
	fxch	ST(2)
	fmul	QWORD PTR _RE2
	fxch	ST(2)
	fadd	QWORD PTR $T98333
	fxch	ST(2)
	fadd	QWORD PTR $T98326
	fxch	ST(2)
	fmul	QWORD PTR _RE2
	fadd	QWORD PTR $T98334
	fxch	ST(2)
	fstp	QWORD PTR -160+[esp+208]

; 103  : 	-+BLK12*AZSUM2)/DETAB**2.00D+00)*RRRN)/ALKXX
; 104  : C
; 105  : 	DTHETAF=(5.00D+00*C1*Z**4.00D+00+4.00D+00*C2*Z**3.00D+00+3.00D+00*

	fmul	QWORD PTR -160+[esp+208]
	fxch	ST(1)
	fmul	QWORD PTR _RE2
	fadd	QWORD PTR $T98335
	fxch	ST(1)
	fstp	QWORD PTR _DTHETAF
	fmul	QWORD PTR _RE2
	fld	QWORD PTR _DZZ
	fxch	ST(1)
	fadd	QWORD PTR $T98336
	fld	QWORD PTR _THETAF
	fxch	ST(1)
	fmul	QWORD PTR _RE2
	mov	DWORD PTR $T98308[esp+220], 108		; 0000006cH
	fadd	QWORD PTR $T98320
	fstp	QWORD PTR -168+[esp+208]

; 106  : 	-C3*Z**2.00D+00+2.00D+00*C4*Z+C5)*DZZ
; 107  : C
; 108  : 	DPHIF=(6.00D+00*C6*Z**5.00D+00+5.00D+00*C7*Z**4.00D+00+4.00D+00*C8

	fxch	ST(1)
	fmul	QWORD PTR -168+[esp+208]
	fstp	QWORD PTR _DPHIF
	call	__FIsqrt
	fmul	QWORD PTR $T98313
	fld	QWORD PTR _PHIF
	call	__FIsqrt
	fdivrp	ST(1), ST(0)
	fstp	QWORD PTR -192+[esp+208]
	fld	QWORD PTR $T98321
	fdiv	QWORD PTR _PHIF
	fstp	QWORD PTR -176+[esp+208]
	fld	QWORD PTR _PHIF
	fld	QWORD PTR $T98313
	call	__FIfexp
	fdivr	QWORD PTR _THETAF
	mov	DWORD PTR $T98308[esp+220], 111		; 0000006fH
	fstp	QWORD PTR -184+[esp+208]

; 109  : 	-*Z**3.00D+00+3.00D+00*C9*Z**2.00D+00+2.00D+00*C10*Z+C11)*DZZ
; 110  : C
; 111  : 	DF0=(DSQRT(PHIF)/(2.00D+00*DSQRT(THETAF)))*(-(THETAF/PHIF**2.00D+0

	fld	QWORD PTR _DTHETAF
	fmul	QWORD PTR -176+[esp+208]
	fld	QWORD PTR _DPHIF
	fmul	QWORD PTR -184+[esp+208]
	fld	QWORD PTR _THETAF
	fdiv	QWORD PTR _PHIF
	mov	DWORD PTR $T98308[esp+220], 114		; 00000072H
	fxch	ST(1)
	fsubp	ST(2), ST(0)
	fxch	ST(1)
	fmul	QWORD PTR -192+[esp+208]
	fstp	QWORD PTR _DF0

; 112  : 	-0)*DPHIF+(1.00D+00/PHIF)*DTHETAF)
; 113  : C
; 114  : 	FF=DSQRT(THETAF/PHIF)

	call	__FIsqrt
	fstp	QWORD PTR _FF
	fld	QWORD PTR _BLK12
	fadd	QWORD PTR _ALK2
	mov	DWORD PTR $T98308[esp+220], 117		; 00000075H

; 115  : C
; 116  : C
; 117  : 	PKL=DSQRT(1.00D+00/(PI*(ALK2+BLK12)))

	fmul	QWORD PTR $T98312
	fdivr	QWORD PTR $T98321
	call	__FIsqrt
	fstp	QWORD PTR _PKL
	fld	QWORD PTR _PKL
	fmul	QWORD PTR -72+[esp+208]
	fld	QWORD PTR [ebp]
	fld	QWORD PTR -72+[esp+208]
	fmul	QWORD PTR $T98313
	fxch	ST(2)
	fstp	QWORD PTR _ANKL
	fmul	QWORD PTR _ANKL
	fxch	ST(1)
	fdivr	QWORD PTR $T98321
	mov	DWORD PTR $T98308[esp+220], 119		; 00000077H

; 118  : C
; 119  : 	ANKL=OKL*PKL

	mov	DWORD PTR $T98308[esp+220], 120		; 00000078H

; 120  : 	AMKL=2.00D+00*OVERLAP*ANKL

	fxch	ST(1)
	fmul	QWORD PTR $T98313
	fxch	ST(1)
	fstp	QWORD PTR -96+[esp+208]
	fstp	QWORD PTR _AMKL
	fld	QWORD PTR _ALKXX
	fmul	QWORD PTR -96+[esp+208]
	fld	QWORD PTR _PKL
	fld	QWORD PTR [ebp]
	fxch	ST(2)
	fstp	QWORD PTR -24+[esp+208]
	fmul	QWORD PTR -24+[esp+208]
	fld	QWORD PTR _DOVN
	fmul	QWORD PTR _ANKL
	fxch	ST(1)
	fstp	QWORD PTR -24+[esp+208]
	fld	QWORD PTR _FF
	fxch	ST(2)
	fmul	QWORD PTR -24+[esp+208]
	faddp	ST(1), ST(0)
	fmul	QWORD PTR $T98313
	fstp	QWORD PTR -24+[esp+208]
	fmul	QWORD PTR -24+[esp+208]
	fstp	QWORD PTR -48+[esp+208]
	mov	DWORD PTR $T98308[esp+220], 127		; 0000007fH

; 121  : 	DOKL=(1.00D+00/(2.00D+00*DSQRT(DETAB)))*DDETAB
; 122  : C
; 123  : 	DNkl=(PKL*DOKL)
; 124  : C
; 125  : 	DMKL=2.00D+00*(OVERLAP*DNKL + ANKL*DOVnDA1)
; 126  : C
; 127  : 	DENDA1=AMKL*DF0 + FF*DMKL

	fld	QWORD PTR _DF0
	fmul	QWORD PTR _AMKL
	fld	QWORD PTR _BLK12
	fadd	QWORD PTR _ALK1
	fxch	ST(1)
	fadd	QWORD PTR -48+[esp+208]
	mov	ecx, DWORD PTR _CZ$[esp+204]
	fxch	ST(1)
	fstp	QWORD PTR -32+[esp+208]
	fstp	QWORD PTR _DEN
	fld	QWORD PTR [ecx]
	fmul	QWORD PTR _BLK12
	fld	QWORD PTR -32+[esp+208]
	fmul	QWORD PTR -80+[esp+208]
	fxch	ST(1)
	fadd	QWORD PTR _A2Z
	fld	QWORD PTR _ALKXX
	fxch	ST(2)
	fdiv	QWORD PTR -88+[esp+208]
	fxch	ST(1)
	fadd	QWORD PTR _A1Z
	fld	QWORD PTR $T98313
	fxch	ST(1)
	fdiv	QWORD PTR _DETAB
	fsubrp	ST(2), ST(0)
	fxch	ST(1)
	fstp	QWORD PTR -56+[esp+208]
	call	__FIfexp
	fld	QWORD PTR _DETAB
	fmul	QWORD PTR _QKL
	fld	QWORD PTR -32+[esp+208]
	fmul	QWORD PTR _QKL
	fxch	ST(1)
	fdivrp	ST(2), ST(0)
	mov	DWORD PTR $T98308[esp+220], 133		; 00000085H
	fdiv	QWORD PTR _ALKXX
	fstp	QWORD PTR -64+[esp+208]
	fstp	QWORD PTR -144+[esp+208]

; 128  : C
; 129  : C	DENDA2:
; 130  : C	======================================================
; 131  : 	DDETAB=ALK1+BLK12
; 132  : C
; 133  : 	DZZ=(DDETAB)*QKL/ALKXX-DETAB*QKL/ALKXX**2.00D+00+DETAB*(

	fld	QWORD PTR -56+[esp+208]
	fmul	QWORD PTR _RRRN
	mov	DWORD PTR $T98308[esp+220], 137		; 00000089H
	fmul	QWORD PTR $T98313
	fmul	QWORD PTR _DETAB
	fdiv	QWORD PTR _ALKXX
	mov	DWORD PTR $T98308[esp+220], 140		; 0000008cH
	fsub	QWORD PTR -144+[esp+208]
	fadd	QWORD PTR -64+[esp+208]
	fstp	QWORD PTR _DZZ

; 134  : 	-2.00D+00*((A1Z+A2Z+BLK12*CZ)/DETAB-(ALK1+BLK12)*(ALKX
; 135  :      -X*AZSUM1+BLK12*AZSUM2)/DETAB**2.00D+00)*RRRN)/ALKXX
; 136  : C
; 137  : 	DTHETAF=(5.00D+00*C1*Z**4.00D+00+4.00D+00*C2*Z**3.00D+00+3.00D+00*

	fld	QWORD PTR _DZZ
	fmul	QWORD PTR -160+[esp+208]
	fld	QWORD PTR _DZZ
	fmul	QWORD PTR -168+[esp+208]
	fxch	ST(1)
	fstp	QWORD PTR _DTHETAF

; 138  : 	-C3*Z**2.00D+00+2.00D+00*C4*Z+C5)*DZZ
; 139  : C
; 140  : 	DPHIF=(6.00D+00*C6*Z**5.00D+00+5.00D+00*C7*Z**4.00D+00+4.00D+00*C8

	fstp	QWORD PTR _DPHIF
	fld	QWORD PTR _DTHETAF
	fmul	QWORD PTR -176+[esp+208]
	fld	QWORD PTR _DPHIF
	fmul	QWORD PTR -184+[esp+208]
	fld	QWORD PTR _BLK12
	fmul	QWORD PTR $T98312
	fxch	ST(1)
	fsubp	ST(2), ST(0)
	fld	QWORD PTR _ALK2
	fmul	QWORD PTR $T98312
	fxch	ST(2)
	fmul	QWORD PTR -192+[esp+208]
	mov	DWORD PTR $T98308[esp+220], 143		; 0000008fH

; 141  : 	-*Z**3.00D+00+3.00D+00*C9*Z**2.00D+00+2.00D+00*C10*Z+C11)*DZZ
; 142  : C
; 143  : 	DF0=(DSQRT(PHIF)/(2.00D+00*DSQRT(THETAF)))*(-(THETAF/PHIF**2.00D+0

	fstp	QWORD PTR _DF0
	faddp	ST(1), ST(0)
	fstp	QWORD PTR -24+[esp+208]
	fld	QWORD PTR -24+[esp+208]
	call	__FIsqrt
	fmul	QWORD PTR -24+[esp+208]
	fld	QWORD PTR -96+[esp+208]
	fmul	QWORD PTR -32+[esp+208]
	fld	QWORD PTR -72+[esp+208]
	fxch	ST(2)
	fmul	QWORD PTR $T98313
	fxch	ST(1)
	fstp	QWORD PTR -24+[esp+208]
	fdivr	QWORD PTR $T98337
	fmul	QWORD PTR $T98312
	fld	QWORD PTR _PKL
	fmul	QWORD PTR -24+[esp+208]
	fxch	ST(1)
	fstp	QWORD PTR -112+[esp+208]
	fld	QWORD PTR [ebp]
	fxch	ST(2)
	fmul	QWORD PTR -112+[esp+208]
	fld	QWORD PTR _DOVN+8
	fmul	QWORD PTR _ANKL
	fxch	ST(1)
	fstp	QWORD PTR -32+[esp+208]
	fxch	ST(1)
	fadd	QWORD PTR -32+[esp+208]
	fld	QWORD PTR _FF
	fxch	ST(1)
	fstp	QWORD PTR -24+[esp+208]
	fxch	ST(2)
	fmul	QWORD PTR -24+[esp+208]
	mov	DWORD PTR $T98308[esp+220], 155		; 0000009bH
	faddp	ST(1), ST(0)
	fld	QWORD PTR _DF0
	fmul	QWORD PTR _AMKL
	fxch	ST(1)
	fmul	QWORD PTR $T98313
	fld	QWORD PTR _ALK2
	fadd	QWORD PTR _ALK1
	fxch	ST(1)
	fstp	QWORD PTR -24+[esp+208]
	fstp	QWORD PTR -152+[esp+208]
	fld	QWORD PTR -152+[esp+208]
	fmul	QWORD PTR -80+[esp+208]
	fxch	ST(2)
	fmul	QWORD PTR -24+[esp+208]
	fxch	ST(2)
	fdiv	QWORD PTR -88+[esp+208]
	fxch	ST(2)
	fstp	QWORD PTR -72+[esp+208]

; 144  : 	-0)*DPHIF+(1.00D+00/PHIF)*DTHETAF)
; 145  : C
; 146  : 	DPkl=(-1.00D+00/(2.00D+00*((PI*ALK2+PI*BLK12)*DSQRT((PI*ALK2+PI*BL
; 147  : 	-K12)))))*PI
; 148  : C
; 149  : 	DOKL=(1.00D+00/(2.00D+00*DSQRT(DETAB)))*DDETAB
; 150  : C
; 151  : 	DNkl=(OKL*DPKL + PKL*DOKL)
; 152  : C
; 153  : 	DMKL=2.00D+00*(OVERLAP*DNKL + ANKL*DOVnDA2)
; 154  : C
; 155  : 	DENDA2=AMKL*DF0 + FF*DMKL

	fadd	QWORD PTR -72+[esp+208]
	fld	QWORD PTR _AZSUM
	fdiv	QWORD PTR _DETAB
	fxch	ST(1)
	fstp	QWORD PTR _DEN+8
	fsubrp	ST(1), ST(0)
	fld	QWORD PTR -152+[esp+208]
	fmul	QWORD PTR _QKL
	fxch	ST(1)
	fstp	QWORD PTR -80+[esp+208]
	fdiv	QWORD PTR _ALKXX
	fstp	QWORD PTR -88+[esp+208]
	mov	DWORD PTR $T98308[esp+220], 161		; 000000a1H

; 156  : C
; 157  : C	DENDA5:
; 158  : C	======================================================
; 159  : 	DDETAB=ALK1+ALK2
; 160  : C
; 161  : 	DZZ=(DDETAB)*QKL/ALKXX-DETAB*QKL/ALKXX**2.00D+00+DETAB*(2.00D+00*

	fld	QWORD PTR -80+[esp+208]
	fmul	QWORD PTR _RRRN
	mov	DWORD PTR $T98308[esp+220], 165		; 000000a5H
	fmul	QWORD PTR $T98313
	fmul	QWORD PTR _DETAB
	fdiv	QWORD PTR _ALKXX
	mov	DWORD PTR $T98308[esp+220], 168		; 000000a8H
	fsub	QWORD PTR -144+[esp+208]
	fadd	QWORD PTR -88+[esp+208]
	fstp	QWORD PTR _DZZ

; 162  : 	-(AZSUM/DETAB-(ALK1+ALK2)*(ALKXX*AZSUM1
; 163  :      -+BLK12*AZSUM2)/DETAB**2.00D+00)*RRRN)/ALKXX
; 164  : C
; 165  : 	DTHETAF=(5.00D+00*C1*Z**4.00D+00+4.00D+00*C2*Z**3.00D+00+3.00D+00*

	fld	QWORD PTR _DZZ
	fmul	QWORD PTR -160+[esp+208]
	fld	QWORD PTR _DZZ
	fmul	QWORD PTR -168+[esp+208]
	fxch	ST(1)
	fstp	QWORD PTR _DTHETAF

; 166  : 	-C3*Z**2.00D+00+2.00D+00*C4*Z+C5)*DZZ
; 167  : C
; 168  : 	DPHIF=(6.00D+00*C6*Z**5.00D+00+5.00D+00*C7*Z**4.00D+00+4.00D+00*C8

	fstp	QWORD PTR _DPHIF
	fld	QWORD PTR _DTHETAF
	fmul	QWORD PTR -176+[esp+208]
	fld	QWORD PTR _DPHIF
	fmul	QWORD PTR -184+[esp+208]
	fld	QWORD PTR -96+[esp+208]
	fmul	QWORD PTR -152+[esp+208]
	fxch	ST(1)
	fsubp	ST(2), ST(0)
	fld	QWORD PTR _PKL
	fxch	ST(1)
	fstp	QWORD PTR -120+[esp+208]
	fmul	QWORD PTR -120+[esp+208]
	fld	QWORD PTR [ebp]
	fxch	ST(2)
	fmul	QWORD PTR -192+[esp+208]
	fxch	ST(1)
	fadd	QWORD PTR -32+[esp+208]
	fld	QWORD PTR _DOVN+16
	fmul	QWORD PTR _ANKL
	fxch	ST(2)
	fstp	QWORD PTR _DF0
	fstp	QWORD PTR -128+[esp+208]
	fld	QWORD PTR _FF
	fxch	ST(2)
	fmul	QWORD PTR -128+[esp+208]
	mov	DWORD PTR $T98308[esp+220], 171		; 000000abH

; 169  : 	-*Z**3.00D+00+3.00D+00*C9*Z**2.00D+00+2.00D+00*C10*Z+C11)*DZZ
; 170  : C
; 171  : 	DF0=(DSQRT(PHIF)/(2.00D+00*DSQRT(THETAF)))*(-(THETAF/PHIF**2.00D+0

	fld	QWORD PTR _DF0
	fmul	QWORD PTR _AMKL
	fxch	ST(2)
	faddp	ST(1), ST(0)
	mov	DWORD PTR $T98308[esp+220], 183		; 000000b7H
	mov	DWORD PTR $T98308[esp+220], 187		; 000000bbH
	fmul	QWORD PTR $T98313
	fld	QWORD PTR [esi]
	fmul	QWORD PTR _RRRN
	fxch	ST(1)
	fstp	QWORD PTR -24+[esp+208]
	fmul	QWORD PTR $T98313
	fxch	ST(2)
	fmul	QWORD PTR -24+[esp+208]
	mov	DWORD PTR $T98308[esp+220], 189		; 000000bdH
	fxch	ST(2)
	fstp	QWORD PTR _DZZ
	fld	QWORD PTR _DZZ
	fmul	QWORD PTR -160+[esp+208]
	fxch	ST(2)
	fstp	QWORD PTR -32+[esp+208]

; 172  : 	-0)*DPHIF+(1.00D+00/PHIF)*DTHETAF)
; 173  : C
; 174  : 	DPkl=(-1.00D+00/(2.00D+00*((PI*ALK2+PI*BLK12)*DSQRT((PI*ALK2+PI*BL
; 175  : 	-K12)))))*PI
; 176  : C
; 177  : 	DOKL=(1.00D+00/(2.00D+00*DSQRT(DETAB)))*DDETAB
; 178  : C
; 179  : 	DNkl=(OKL*DPKL + PKL*DOKL)
; 180  : C
; 181  : 	DMKL=2.00D+00*(OVERLAP*DNKL + ANKL*DOVnDA5)
; 182  : C
; 183  : 	DENDA5=AMKL*DF0 + FF*DMKL

	fadd	QWORD PTR -32+[esp+208]
	fld	QWORD PTR _DZZ
	fmul	QWORD PTR -168+[esp+208]
	fxch	ST(2)
	fstp	QWORD PTR _DTHETAF
	fstp	QWORD PTR _DEN+16

; 184  : C
; 185  : C	DENDAZ:
; 186  : C	======================================================
; 187  : 	DZZ=2.00D+00*A1*RRRN
; 188  : C
; 189  : 	DTHETAF=(5.00D+00*C1*Z**4.00D+00+4.00D+00*C2*Z**3.00D+00+3.00D+00*

	fstp	QWORD PTR _DPHIF
	fld	QWORD PTR _DTHETAF
	fmul	QWORD PTR -176+[esp+208]
	fld	QWORD PTR _DPHIF
	fmul	QWORD PTR -184+[esp+208]
	fld	QWORD PTR _DOVN+24
	fmul	QWORD PTR _ANKL
	fxch	ST(1)
	fsubp	ST(2), ST(0)
	fld	QWORD PTR _FF
	fxch	ST(1)
	fmul	QWORD PTR $T98313
	fxch	ST(2)
	fmul	QWORD PTR -192+[esp+208]
	mov	DWORD PTR $T98308[esp+220], 192		; 000000c0H

; 190  : 	-C3*Z**2.00D+00+2.00D+00*C4*Z+C5)*DZZ
; 191  : C
; 192  : 	DPHIF=(6.00D+00*C6*Z**5.00D+00+5.00D+00*C7*Z**4.00D+00+4.00D+00*C8

	fxch	ST(2)
	fstp	QWORD PTR -24+[esp+208]
	fmul	QWORD PTR -24+[esp+208]
	fxch	ST(1)
	fstp	QWORD PTR _DF0
	fstp	QWORD PTR -24+[esp+208]
	fld	QWORD PTR _DF0
	fmul	QWORD PTR _AMKL
	fld	QWORD PTR [edi]
	fmul	QWORD PTR _BLK12
	fxch	ST(1)
	fadd	QWORD PTR -24+[esp+208]
	mov	DWORD PTR $T98308[esp+220], 195		; 000000c3H

; 193  : 	-*Z**3.00D+00+3.00D+00*C9*Z**2.00D+00+2.00D+00*C10*Z+C11)*DZZ
; 194  : C
; 195  : 	DF0=(DSQRT(PHIF)/(2.00D+00*DSQRT(THETAF)))*(-(THETAF/PHIF**2.00D+0

	fxch	ST(1)
	fstp	QWORD PTR -96+[esp+208]
	fstp	QWORD PTR _DEN+24
	fld	QWORD PTR -96+[esp+208]
	fmul	QWORD PTR _RRRN
	mov	DWORD PTR $T98308[esp+220], 200		; 000000c8H

; 196  : 	-0)*DPHIF+(1.00D+00/PHIF)*DTHETAF)
; 197  : C
; 198  : 	DMKL=2.00D+00*(ANKL*DOVnDAZ)
; 199  : C
; 200  : 	DENDAZ=AMKL*DF0 + FF*DMKL

	mov	DWORD PTR $T98308[esp+220], 204		; 000000ccH

; 201  : C
; 202  : C	DENDBZ:
; 203  : C	======================================================
; 204  : 	DZZ=2.00D+00*A3*BLK12*RRRN/ALKXX

	fmul	QWORD PTR $T98313
	mov	DWORD PTR $T98308[esp+220], 206		; 000000ceH
	mov	DWORD PTR $T98308[esp+220], 209		; 000000d1H
	fdiv	QWORD PTR _ALKXX
	fstp	QWORD PTR _DZZ

; 205  : C
; 206  : 	DTHETAF=(5.00D+00*C1*Z**4.00D+00+4.00D+00*C2*Z**3.00D+00+3.00D+00*

	fld	QWORD PTR _DZZ
	fmul	QWORD PTR -160+[esp+208]
	fld	QWORD PTR _DZZ
	fmul	QWORD PTR -168+[esp+208]
	fxch	ST(1)
	fstp	QWORD PTR _DTHETAF

; 207  : 	-C3*Z**2.00D+00+2.00D+00*C4*Z+C5)*DZZ
; 208  : C
; 209  : 	DPHIF=(6.00D+00*C6*Z**5.00D+00+5.00D+00*C7*Z**4.00D+00+4.00D+00*C8

	fstp	QWORD PTR _DPHIF
	fld	QWORD PTR _DTHETAF
	fmul	QWORD PTR -176+[esp+208]
	fld	QWORD PTR _DPHIF
	fmul	QWORD PTR -184+[esp+208]
	fld	QWORD PTR _DOVN+32
	fmul	QWORD PTR _ANKL
	fxch	ST(1)
	fsubp	ST(2), ST(0)
	fld	QWORD PTR _FF
	fxch	ST(1)
	fmul	QWORD PTR $T98313
	fxch	ST(2)
	fmul	QWORD PTR -192+[esp+208]
	mov	DWORD PTR $T98308[esp+220], 212		; 000000d4H

; 210  : 	-*Z**3.00D+00+3.00D+00*C9*Z**2.00D+00+2.00D+00*C10*Z+C11)*DZZ
; 211  : C
; 212  : 	DF0=(DSQRT(PHIF)/(2.00D+00*DSQRT(THETAF)))*(-(THETAF/PHIF**2.00D+0

	fxch	ST(2)
	fstp	QWORD PTR -136+[esp+208]
	fmul	QWORD PTR -136+[esp+208]
	fxch	ST(1)
	fstp	QWORD PTR _DF0
	fstp	QWORD PTR -104+[esp+208]
	mov	DWORD PTR $T98308[esp+220], 217		; 000000d9H

; 213  : 	-0)*DPHIF+(1.00D+00/PHIF)*DTHETAF)
; 214  : C
; 215  : 	DMKL=2.00D+00*(ANKL*DOVnDBZ)
; 216  : C
; 217  : 	DENDBZ=AMKL*DF0 + FF*DMKL

	fld	QWORD PTR _DF0
	fmul	QWORD PTR _AMKL
	fld	QWORD PTR -40+[esp+208]
	fmul	QWORD PTR _RRR
	fxch	ST(1)
	fadd	QWORD PTR -104+[esp+208]
	fxch	ST(1)
	fmul	QWORD PTR $T98313
	mov	DWORD PTR $T98308[esp+220], 224		; 000000e0H
	mov	DWORD PTR $T98308[esp+220], 227		; 000000e3H
	fmul	QWORD PTR _DETAB
	fxch	ST(1)
	fstp	QWORD PTR _DEN+32

; 218  : C
; 219  : C
; 220  : C	DENDA1P:
; 221  : C	======================================================
; 222  : 	DDETAB=ALKXX
; 223  : C
; 224  : 	DZZ=QKL+DETAB*(2.00D+00*(ALKXX*AZ/DETAB-ALKXX*(ALKXX*AZSUM1

	fdiv	QWORD PTR _ALKXX
	mov	DWORD PTR $T98308[esp+220], 230		; 000000e6H
	fadd	QWORD PTR _QKL
	fstp	QWORD PTR _DZZ

; 225  : 	-+BLK12*AZSUM2)/DETAB**2.00D+00)*RRR)/ALKXX
; 226  : C
; 227  : 	DTHETAF=(5.00D+00*C1*Z**4.00D+00+4.00D+00*C2*Z**3.00D+00+3.00D+00*

	fld	QWORD PTR _DZZ
	fmul	QWORD PTR -160+[esp+208]
	fld	QWORD PTR _DZZ
	fmul	QWORD PTR -168+[esp+208]
	fxch	ST(1)
	fstp	QWORD PTR _DTHETAF

; 228  : 	-C3*Z**2.00D+00+2.00D+00*C4*Z+C5)*DZZ
; 229  : C
; 230  : 	DPHIF=(6.00D+00*C6*Z**5.00D+00+5.00D+00*C7*Z**4.00D+00+4.00D+00*C8

	fstp	QWORD PTR _DPHIF
	fld	QWORD PTR _DTHETAF
	fmul	QWORD PTR -176+[esp+208]
	fld	QWORD PTR _DPHIF
	fmul	QWORD PTR -184+[esp+208]
	mov	DWORD PTR $T98308[esp+220], 233		; 000000e9H

; 231  : 	-*Z**3.00D+00+3.00D+00*C9*Z**2.00D+00+2.00D+00*C10*Z+C11)*DZZ
; 232  : C
; 233  : 	DF0=(DSQRT(PHIF)/(2.00D+00*DSQRT(THETAF)))*(-(THETAF/PHIF**2.00D+0

	mov	DWORD PTR $T98308[esp+220], 242		; 000000f2H
	fsubp	ST(1), ST(0)
	mov	DWORD PTR $T98308[esp+220], 248		; 000000f8H
	fmul	QWORD PTR -192+[esp+208]
	fstp	QWORD PTR _DF0

; 234  : 	-0)*DPHIF+(1.00D+00/PHIF)*DTHETAF)
; 235  : C
; 236  : 	DOKL=(1.00D+00/(2.00D+00*DSQRT(DETAB)))*DDETAB
; 237  : C
; 238  : 	DNkl=(PKL*DOKL)
; 239  : C
; 240  : 	DMKL=2.00D+00*(OVERLAP*DNKL + ANKL*DOVnDA1)
; 241  : C
; 242  : 	DENDA1P=AMKL*DF0 + FF*DMKL

	fld	QWORD PTR _DF0
	fmul	QWORD PTR _AMKL
	fld	QWORD PTR -56+[esp+208]
	fmul	QWORD PTR _RRR
	fxch	ST(1)
	fadd	QWORD PTR -48+[esp+208]
	fxch	ST(1)
	fmul	QWORD PTR $T98313
	mov	DWORD PTR $T98308[esp+220], 253		; 000000fdH
	fmul	QWORD PTR _DETAB
	fxch	ST(1)
	fstp	QWORD PTR _DEN+40

; 243  : C
; 244  : C	DENDA2P:
; 245  : C	======================================================
; 246  : 	DDETAB=ALK1+BLK12
; 247  : C
; 248  : 	DZZ=(DDETAB)*QKL/ALKXX-DETAB*QKL/ALKXX**2.00D+00+DETAB*(2.00D+00*

	fdiv	QWORD PTR _ALKXX
	mov	DWORD PTR $T98308[esp+220], 256		; 00000100H
	fsub	QWORD PTR -144+[esp+208]
	fadd	QWORD PTR -64+[esp+208]
	fstp	QWORD PTR _DZZ

; 249  : 	-((A1Z+A2Z+BLK12*CZ)/
; 250  : 	-DETAB-(ALK1+BLK12)*(ALKXX*AZSUM1+BLK12*AZSUM2)/
; 251  :      -DETAB**2.00D+00)*RRR)/ALKXX
; 252  : C
; 253  : 	DTHETAF=(5.00D+00*C1*Z**4.00D+00+4.00D+00*C2*Z**3.00D+00+3.00D+00*

	fld	QWORD PTR _DZZ
	fmul	QWORD PTR -160+[esp+208]
	fld	QWORD PTR _DZZ
	fmul	QWORD PTR -168+[esp+208]
	fxch	ST(1)
	fstp	QWORD PTR _DTHETAF

; 254  : 	-C3*Z**2.00D+00+2.00D+00*C4*Z+C5)*DZZ
; 255  : C
; 256  : 	DPHIF=(6.00D+00*C6*Z**5.00D+00+5.00D+00*C7*Z**4.00D+00+4.00D+00*C8

	fstp	QWORD PTR _DPHIF
	fld	QWORD PTR _DTHETAF
	fmul	QWORD PTR -176+[esp+208]
	fld	QWORD PTR _DPHIF
	fmul	QWORD PTR -184+[esp+208]
	mov	DWORD PTR $T98308[esp+220], 259		; 00000103H

; 257  : 	-*Z**3.00D+00+3.00D+00*C9*Z**2.00D+00+2.00D+00*C10*Z+C11)*DZZ
; 258  : C
; 259  : 	DF0=(DSQRT(PHIF)/(2.00D+00*DSQRT(THETAF)))*(-(THETAF/PHIF**2.00D+0

	mov	DWORD PTR $T98308[esp+220], 271		; 0000010fH
	fsubp	ST(1), ST(0)
	mov	DWORD PTR $T98308[esp+220], 275		; 00000113H
	mov	ecx, DWORD PTR -152+[esp+212]
	fmul	QWORD PTR -192+[esp+208]
	mov	edx, DWORD PTR -152+[esp+208]
	mov	DWORD PTR _DDETAB+4, ecx
	mov	DWORD PTR $T98308[esp+220], 277		; 00000115H
	mov	DWORD PTR _DDETAB, edx
	fstp	QWORD PTR _DF0

; 260  : 	-0)*DPHIF+(1.00D+00/PHIF)*DTHETAF)
; 261  : C
; 262  : 	DPkl=(-1.00D+00/(2.00D+00*((PI*ALK2+PI*BLK12)*DSQRT((PI*ALK2+PI*BL
; 263  : 	-K12)))))*PI
; 264  : C
; 265  : 	DOKL=(1.00D+00/(2.00D+00*DSQRT(DETAB)))*DDETAB
; 266  : C
; 267  : 	DNkl=(OKL*DPKL + PKL*DOKL)
; 268  : C
; 269  : 	DMKL=2.00D+00*(OVERLAP*DNKL + ANKL*DOVnDA2)
; 270  : C
; 271  : 	DENDA2P=AMKL*DF0 + FF*DMKL

	fld	QWORD PTR _DF0
	fmul	QWORD PTR _AMKL
	fld	QWORD PTR -80+[esp+208]
	fmul	QWORD PTR _RRR
	fxch	ST(1)
	fadd	QWORD PTR -72+[esp+208]
	fxch	ST(1)
	fmul	QWORD PTR $T98313
	mov	DWORD PTR $T98308[esp+220], 282		; 0000011aH
	fmul	QWORD PTR _DETAB
	fxch	ST(1)
	fstp	QWORD PTR _DEN+48

; 272  : C
; 273  : C	DENDA5P:
; 274  : C	======================================================
; 275  : 	DDETAB=ALK1+ALK2
; 276  : C
; 277  : 	DZZ=(DDETAB)*QKL/ALKXX-DETAB*QKL/ALKXX**2.00D+00+DETAB*(2.00D+00*

	fdiv	QWORD PTR _ALKXX
	mov	DWORD PTR $T98308[esp+220], 285		; 0000011dH
	fsub	QWORD PTR -144+[esp+208]
	fadd	QWORD PTR -88+[esp+208]
	fstp	QWORD PTR _DZZ

; 278  : 	-(AZSUM/DETAB-(ALK1+ALK2)*
; 279  :      -(ALKXX*AZSUM1+BLK12*AZSUM2)/DETAB**2.00D+00)*RRR)/
; 280  :      -ALKXX
; 281  : C
; 282  : 	DTHETAF=(5.00D+00*C1*Z**4.00D+00+4.00D+00*C2*Z**3.00D+00+3.00D+00*

	fld	QWORD PTR _DZZ
	fmul	QWORD PTR -160+[esp+208]
	fld	QWORD PTR _DZZ
	fmul	QWORD PTR -168+[esp+208]
	fxch	ST(1)
	fstp	QWORD PTR _DTHETAF

; 283  : 	-C3*Z**2.00D+00+2.00D+00*C4*Z+C5)*DZZ
; 284  : C
; 285  : 	DPHIF=(6.00D+00*C6*Z**5.00D+00+5.00D+00*C7*Z**4.00D+00+4.00D+00*C8

	fstp	QWORD PTR _DPHIF
	fld	QWORD PTR _DTHETAF
	fmul	QWORD PTR -176+[esp+208]
	fld	QWORD PTR _DPHIF
	fmul	QWORD PTR -184+[esp+208]
	mov	DWORD PTR $T98308[esp+220], 288		; 00000120H

; 286  : 	-*Z**3.00D+00+3.00D+00*C9*Z**2.00D+00+2.00D+00*C10*Z+C11)*DZZ
; 287  : C
; 288  : 	DF0=(DSQRT(PHIF)/(2.00D+00*DSQRT(THETAF)))*(-(THETAF/PHIF**2.00D+0

	mov	ecx, DWORD PTR -112+[esp+212]
	mov	edx, DWORD PTR -112+[esp+208]
	fsubp	ST(1), ST(0)
	mov	DWORD PTR $T98308[esp+220], 291		; 00000123H
	mov	DWORD PTR $T98308[esp+220], 294		; 00000126H
	fmul	QWORD PTR -192+[esp+208]
	mov	DWORD PTR _DPKL+4, ecx
	mov	ecx, DWORD PTR -120+[esp+212]
	mov	DWORD PTR _DPKL, edx
	fstp	QWORD PTR _DF0

; 289  : 	-0)*DPHIF+(1.00D+00/PHIF)*DTHETAF)
; 290  : C
; 291  : 	DPkl=(-1.00D+00/(2.00D+00*((PI*ALK2+PI*BLK12)*DSQRT((PI*ALK2+PI*BL
; 292  : 	-K12)))))*PI
; 293  : C
; 294  : 	DOKL=(1.00D+00/(2.00D+00*DSQRT(DETAB)))*DDETAB

	mov	DWORD PTR _DOKL+4, ecx
	mov	edx, DWORD PTR -120+[esp+208]
	mov	ecx, DWORD PTR -128+[esp+212]
	fld	QWORD PTR _DF0
	fmul	QWORD PTR _AMKL
	mov	DWORD PTR _DOKL, edx
	mov	DWORD PTR _DNKL+4, ecx
	mov	DWORD PTR $T98308[esp+220], 296		; 00000128H

; 295  : C
; 296  : 	DNkl=(OKL*DPKL + PKL*DOKL)

	fadd	QWORD PTR -32+[esp+208]
	fld	QWORD PTR [esi]
	fmul	QWORD PTR _RRR
	fxch	ST(1)
	fstp	QWORD PTR _DEN+56
	fmul	QWORD PTR $T98313
	mov	DWORD PTR $T98308[esp+220], 300		; 0000012cH
	mov	DWORD PTR $T98308[esp+220], 304		; 00000130H
	mov	DWORD PTR $T98308[esp+220], 306		; 00000132H
	fstp	QWORD PTR _DZZ
	fld	QWORD PTR _DZZ
	fmul	QWORD PTR -160+[esp+208]
	fld	QWORD PTR _DZZ
	fmul	QWORD PTR -168+[esp+208]
	fxch	ST(1)
	fstp	QWORD PTR _DTHETAF
	fstp	QWORD PTR _DPHIF
	fld	QWORD PTR _DTHETAF
	fmul	QWORD PTR -176+[esp+208]
	fld	QWORD PTR _DPHIF
	fmul	QWORD PTR -184+[esp+208]
	mov	DWORD PTR $T98308[esp+220], 309		; 00000135H
	mov	DWORD PTR $T98308[esp+220], 312		; 00000138H
	fsubp	ST(1), ST(0)
	mov	DWORD PTR $T98308[esp+220], 317		; 0000013dH
	mov	edx, DWORD PTR -128+[esp+208]
	fmul	QWORD PTR -192+[esp+208]
	mov	DWORD PTR $T98308[esp+220], 321		; 00000141H
	mov	DWORD PTR _DNKL, edx

; 297  : C
; 298  : 	DMKL=2.00D+00*(OVERLAP*DNKL + ANKL*DOVnDA5)
; 299  : C
; 300  : 	DENDA5P=AMKL*DF0 + FF*DMKL
; 301  : C
; 302  : C	DENDAZP:
; 303  : C	======================================================
; 304  : 	DZZ=2.00D+00*A1*RRR
; 305  : C
; 306  : 	DTHETAF=(5.00D+00*C1*Z**4.00D+00+4.00D+00*C2*Z**3.00D+00+3.00D+00*
; 307  : 	-C3*Z**2.00D+00+2.00D+00*C4*Z+C5)*DZZ
; 308  : C
; 309  : 	DPHIF=(6.00D+00*C6*Z**5.00D+00+5.00D+00*C7*Z**4.00D+00+4.00D+00*C8
; 310  : 	-*Z**3.00D+00+3.00D+00*C9*Z**2.00D+00+2.00D+00*C10*Z+C11)*DZZ
; 311  : C
; 312  : 	DF0=(DSQRT(PHIF)/(2.00D+00*DSQRT(THETAF)))*(-(THETAF/PHIF**2.00D+0

	fstp	QWORD PTR _DF0

; 313  : 	-0)*DPHIF+(1.00D+00/PHIF)*DTHETAF)
; 314  : C
; 315  : 	DMKL=2.00D+00*(ANKL*DOVnDAZ)
; 316  : C
; 317  : 	DENDAZP=AMKL*DF0 + FF*DMKL

	fld	QWORD PTR _DF0
	fmul	QWORD PTR _AMKL
	fld	QWORD PTR -96+[esp+208]
	fmul	QWORD PTR _RRR
	fxch	ST(1)
	fadd	QWORD PTR -24+[esp+208]
	fxch	ST(1)
	fmul	QWORD PTR $T98313
	mov	DWORD PTR $T98308[esp+220], 323		; 00000143H
	mov	DWORD PTR $T98308[esp+220], 326		; 00000146H
	fdiv	QWORD PTR _ALKXX
	fstp	QWORD PTR _DZZ
	fstp	QWORD PTR _DEN+64

; 318  : C
; 319  : C	DENDBZP:
; 320  : C	======================================================
; 321  : 	DZZ=2.00D+00*A3*BLK12*RRR/ALKXX
; 322  : C
; 323  : 	DTHETAF=(5.00D+00*C1*Z**4.00D+00+4.00D+00*C2*Z**3.00D+00+3.00D+00*

	fld	QWORD PTR _DZZ
	fmul	QWORD PTR -160+[esp+208]
	fld	QWORD PTR _DZZ
	fmul	QWORD PTR -168+[esp+208]
	fxch	ST(1)
	fstp	QWORD PTR _DTHETAF

; 324  : 	-C3*Z**2.00D+00+2.00D+00*C4*Z+C5)*DZZ
; 325  : C
; 326  : 	DPHIF=(6.00D+00*C6*Z**5.00D+00+5.00D+00*C7*Z**4.00D+00+4.00D+00*C8

	fstp	QWORD PTR _DPHIF
	fld	QWORD PTR _DTHETAF
	fmul	QWORD PTR -176+[esp+208]
	fld	QWORD PTR _DPHIF
	fmul	QWORD PTR -184+[esp+208]
	mov	DWORD PTR $T98308[esp+220], 329		; 00000149H

; 327  : 	-*Z**3.00D+00+3.00D+00*C9*Z**2.00D+00+2.00D+00*C10*Z+C11)*DZZ
; 328  : C
; 329  : 	DF0=(DSQRT(PHIF)/(2.00D+00*DSQRT(THETAF)))*(-(THETAF/PHIF**2.00D+0

	mov	ecx, DWORD PTR -136+[esp+212]
	mov	edx, DWORD PTR -136+[esp+208]
	fsubp	ST(1), ST(0)
	mov	DWORD PTR $T98308[esp+220], 332		; 0000014cH
	mov	DWORD PTR $T98308[esp+220], 334		; 0000014eH
	fmul	QWORD PTR -192+[esp+208]
	mov	DWORD PTR _DMKL+4, ecx
	mov	DWORD PTR $T98308[esp+220], 337		; 00000151H
	mov	DWORD PTR _DMKL, edx
	fstp	QWORD PTR _DF0

; 330  : 	-0)*DPHIF+(1.00D+00/PHIF)*DTHETAF)
; 331  : C
; 332  : 	DMKL=2.00D+00*(ANKL*DOVnDBZ)
; 333  : C
; 334  : 	DENDBZP=AMKL*DF0 + FF*DMKL

	fld	QWORD PTR _DF0
	fmul	QWORD PTR _AMKL
	fadd	QWORD PTR -104+[esp+208]
	fstp	QWORD PTR _DEN+72
$L80:

; 335  : C
; 336  : C
; 337  : 	END IF

	mov	DWORD PTR $T98308[esp+220], 339		; 00000153H

; 338  : C
; 339  :  300	KLM=5

	mov	ecx, DWORD PTR _KLM$[esp+204]
	mov	DWORD PTR $T98308[esp+220], 343		; 00000157H
	mov	DWORD PTR [ecx], 5

; 340  : C
; 341  : C
; 342  : 1500  RETURN                                                            
; 343  :       END 

	call	___FHunlink@0
	pop	ebp
	pop	edi
	pop	esi
	pop	ebx
	add	esp, 192				; 000000c0H
	ret	100					; 00000064H
	npad	6
	npad	7
_ENUCL@100 ENDP
_TEXT	ENDS
END
