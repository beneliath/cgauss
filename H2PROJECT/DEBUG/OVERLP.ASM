	TITLE	C:\users\dgilmore\Graduate.School\Adamowicz.Research\H2ANL\OVERLP.FOR
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
PUBLIC	_DOV
EXTRN	__FFljj:NEAR
_DATA	SEGMENT
_DOV	DB	028H DUP (?)
	ORG $-40
_C:\users\dgilmore\Graduate.School\Adamowicz.Research\H2ANL\OVERLP.FOR DB 
	DB	'C:\users\dgilmore\Graduate.School\Adamowicz.Research\H2ANL\OV'
	DB	'ERLP.FOR', 00H
_DATA	ENDS
PUBLIC	_OVERLP@88
EXTRN	__FIsqrt:NEAR
EXTRN	__FIexp:NEAR
EXTRN	__FIfexp:NEAR
EXTRN	___FHlink@4:NEAR
EXTRN	___FHunlink@0:NEAR
EXTRN	__fltused:NEAR
_BSS	SEGMENT
_RZ1	DQ	01H DUP (?)
_RZ2	DQ	01H DUP (?)
_DD	DQ	01H DUP (?)
_CKL	DQ	01H DUP (?)
_XK	DQ	01H DUP (?)
_ALK1	DQ	01H DUP (?)
_ALK2	DQ	01H DUP (?)
_BLK12	DQ	01H DUP (?)
_DETSL	DQ	01H DUP (?)
_CCL	DQ	01H DUP (?)
_DETSK	DQ	01H DUP (?)
_CCK	DQ	01H DUP (?)
_DETAB	DQ	01H DUP (?)
_ZERO	DQ	01H DUP (?)
_AB2	DQ	01H DUP (?)
_CD2	DQ	01H DUP (?)
_RX1	DQ	01H DUP (?)
_RX2	DQ	01H DUP (?)
_RY1	DQ	01H DUP (?)
_RY2	DQ	01H DUP (?)
_BSS	ENDS
_DATA	SEGMENT
	ORG $+2
$T98311	DB	'_OVERLP@88', 00H
_DATA	ENDS
_CONST	SEGMENT
$T98313	DQ	04000000000000000r		; 2
$T98314	DQ	04020000000000000r		; 8
$T98315	DQ	04008000000000000r		; 3
$T98316	DQ	04010000000000000r		; 4
$T98317	DQ	03fe8000000000000r		; 0.75
$T98318	DQ	03ff8000000000000r		; 1.5
$T98319	DQ	03fd0000000000000r		; 0.25
$T98320	DQ	04004000000000000r		; 2.5
$T98321	DQ	0c020000000000000r		; -8
$T98322	DQ	0c030000000000000r		; -16
$T98323	DQ	0c000000000000000r		; -2
_CONST	ENDS
_TEXT	SEGMENT
; File C:\users\dgilmore\Graduate.School\Adamowicz.Research\H2ANL\OVERLP.FOR
$T98310 = -16
_DFLAG$ = 12
_AX$ = 16
_AY$ = 20
_AZ$ = 24
_A1$ = 28
_BX$ = 32
_BY$ = 36
_BZ$ = 40
_A2$ = 44
_CX$ = 48
_CY$ = 52
_CZ$ = 56
_A3$ = 60
_DX$ = 64
_DY$ = 68
_DZ$ = 72
_A4$ = 76
_A5$ = 80
_A6$ = 84
_OV$ = 88
_KLM$ = 92
_OVERLP@88 PROC NEAR

; 3    :       SUBROUTINE OVERLP(M,DFLAG,AX,AY,AZ,A1,BX,BY,BZ,A2,                      

	sub	esp, 256				; 00000100H
	lea	eax, DWORD PTR $T98310[esp+256]
	push	ebx
	push	esi
	push	edi
	push	ebp
	push	eax
	call	___FHlink@4
	mov	DWORD PTR $T98310[esp+276], OFFSET FLAT:_C:\users\dgilmore\Graduate.School\Adamowicz.Research\H2ANL\OVERLP.FOR
	mov	DWORD PTR $T98310[esp+280], OFFSET FLAT:$T98311
	mov	edi, DWORD PTR _A1$[esp+268]
	mov	esi, DWORD PTR _A3$[esp+268]
	mov	DWORD PTR $T98310[esp+284], 16		; 00000010H

; 4    :      1                  CX,CY,CZ,A3,DX,DY,DZ,A4,
; 5    :      2                  A5,A6,OV,KLM) 
; 6    : C     ----------------------------------------------------- 
; 7    :       IMPLICIT REAL*8 (A-H,O-Z)
; 8    : 	REAL*8 DETSL,DETSK,CCL,CCK
; 9    : 	INTEGER M,DFLAG
; 10   : 	COMMON /DOV/ DOVDA1,DOVDA2,DOVDA5,DOVDAZ,DOVDBZ
; 11   : 
; 12   : C     -----------------------------------------------------
; 13   :       DATA ONE,TWO,THREE,FOUR,FIVE,SIX,EIGHT/1.,2.,3.,4.,5.,6.,8./
; 14   : 	DATA HALF/0.5/
; 15   : C     -----------------------------------------------------
; 16   : 		ALK1=2.0D+00*A1

	fld	QWORD PTR [edi]
	fmul	QWORD PTR $T98313
	fld	QWORD PTR [esi]
	fmul	QWORD PTR $T98313
	fxch	ST(1)
	fstp	QWORD PTR _ALK1
	fstp	QWORD PTR _ALK2
	mov	DWORD PTR $T98310[esp+284], 17		; 00000011H

; 17   : 		ALK2=2.0D+00*A3                                                       

	mov	DWORD PTR $T98310[esp+284], 18		; 00000012H

; 18   : 		BLK12=2.0D+00*A5  

	mov	DWORD PTR $T98310[esp+284], 19		; 00000013H
	mov	ebp, DWORD PTR _A5$[esp+268]
	fld	QWORD PTR [ebp]
	fmul	QWORD PTR $T98313
	fld	QWORD PTR _ALK2
	fadd	QWORD PTR _ALK1
	fld	QWORD PTR _ALK2
	fmul	QWORD PTR _ALK1
	fxch	ST(2)
	fstp	QWORD PTR _BLK12

; 19   : 		DETSL=(ALK1*ALK2)+((ALK1+ALK2)*BLK12)

	fmul	QWORD PTR _BLK12
	mov	ecx, DWORD PTR _AZ$[esp+268]
	mov	edx, DWORD PTR _CZ$[esp+268]
	faddp	ST(1), ST(0)
	fld	QWORD PTR [ecx]
	fsub	QWORD PTR [edx]
	fxch	ST(1)
	fstp	QWORD PTR _DETSL
	fstp	QWORD PTR -128+[esp+272]
	fld	QWORD PTR -128+[esp+272]
	fld	QWORD PTR $T98313
	call	__FIfexp
	fstp	QWORD PTR -248+[esp+272]
	fld	QWORD PTR [ebp]
	fmul	QWORD PTR [esi]
	mov	DWORD PTR $T98310[esp+284], 20		; 00000014H

; 20   : 	    CCL=8*A1*A3*A5*((AZ-CZ)**2)/DETSL

	mov	DWORD PTR $T98310[esp+284], 22		; 00000016H
	fmul	QWORD PTR [edi]
	mov	ecx, DWORD PTR _A2$[esp+268]
	fmul	QWORD PTR -248+[esp+272]
	fld	QWORD PTR [ecx]
	fmul	QWORD PTR $T98313
	fxch	ST(1)
	fmul	QWORD PTR $T98314
	mov	DWORD PTR $T98310[esp+284], 23		; 00000017H
	mov	edx, DWORD PTR _A4$[esp+268]
	fdiv	QWORD PTR _DETSL
	mov	DWORD PTR $T98310[esp+284], 24		; 00000018H
	mov	ebx, DWORD PTR _A6$[esp+268]
	fstp	QWORD PTR _CCL

; 21   : C
; 22   : 		ALK1=2.0D+00*A2

	fstp	QWORD PTR _ALK1

; 23   : 		ALK2=2.0D+00*A4                                                       

	fld	QWORD PTR [edx]
	fmul	QWORD PTR $T98313
	fld	QWORD PTR [ebx]
	fmul	QWORD PTR $T98313
	fxch	ST(1)
	fstp	QWORD PTR _ALK2

; 24   : 		BLK12=2.0D+00*A6 

	fstp	QWORD PTR _BLK12
	fld	QWORD PTR _ALK2
	fadd	QWORD PTR _ALK1
	fld	QWORD PTR _ALK2
	fmul	QWORD PTR _ALK1
	fxch	ST(1)
	fmul	QWORD PTR _BLK12
	mov	ecx, DWORD PTR _BZ$[esp+268]
	mov	edx, DWORD PTR _DZ$[esp+268]
	mov	DWORD PTR $T98310[esp+284], 25		; 00000019H

; 25   : 		DETSK=(ALK1*ALK2)+((ALK1+ALK2)*BLK12)

	faddp	ST(1), ST(0)
	fld	QWORD PTR [ecx]
	fsub	QWORD PTR [edx]
	fld	QWORD PTR $T98313
	fxch	ST(2)
	fstp	QWORD PTR _DETSK
	mov	DWORD PTR $T98310[esp+284], 26		; 0000001aH

; 26   : 	    CCK=8*A2*A4*A6*((BZ-DZ)**2)/DETSK

	fxch	ST(1)
	call	__FIfexp
	fmul	QWORD PTR [ebx]
	mov	ecx, DWORD PTR _A4$[esp+268]
	mov	edx, DWORD PTR _A2$[esp+268]
	mov	DWORD PTR $T98310[esp+284], 28		; 0000001cH
	fmul	QWORD PTR [ecx]
	fld	QWORD PTR [edx]
	fadd	QWORD PTR [edi]
	fxch	ST(1)
	fmul	QWORD PTR [edx]
	mov	DWORD PTR $T98310[esp+284], 29		; 0000001dH
	fmul	QWORD PTR $T98314
	fxch	ST(1)
	fstp	QWORD PTR _ALK1
	fdiv	QWORD PTR _DETSK
	mov	DWORD PTR $T98310[esp+284], 30		; 0000001eH
	fstp	QWORD PTR _CCK

; 27   : C     -----------------------------------------------------
; 28   :       ALK1=A1+A2
; 29   :       ALK2=A3+A4                                                       

	fld	QWORD PTR [ecx]
	fadd	QWORD PTR [esi]
	fld	QWORD PTR [ebx]
	fadd	QWORD PTR [ebp]
	fxch	ST(1)
	fstp	QWORD PTR _ALK2

; 30   :       BLK12=A5+A6  

	fstp	QWORD PTR _BLK12
	fld	QWORD PTR _ALK2
	fadd	QWORD PTR _ALK1
	fld	QWORD PTR _ALK2
	fmul	QWORD PTR _ALK1
	fxch	ST(1)
	fstp	QWORD PTR -192+[esp+272]
	fstp	QWORD PTR -184+[esp+272]
	mov	DWORD PTR $T98310[esp+284], 31		; 0000001fH

; 31   :       DETAB=ALK1*ALK2+(ALK1+ALK2)*BLK12

	fld	QWORD PTR -192+[esp+272]
	fmul	QWORD PTR _BLK12
	fadd	QWORD PTR -184+[esp+272]
	fcom	QWORD PTR _ZERO
	fstp	QWORD PTR _DETAB
	fwait
	fnstsw	ax
	mov	DWORD PTR $T98310[esp+284], 33		; 00000021H

; 32   : C
; 33   :       IF(DETAB.LE.ZERO) GOTO 300

	test	ah, 65					; 00000041H
	je	SHORT $L98308

; 98   : 	-/ALK2+2.0D+00*ALK1*A3*BLK12*(RZ1-RZ2)/DETAB)
; 99   : C     -----------------------------------------------------
; 100  :       GOTO 1500
; 101  : C     -----------------------------------------------------
; 102  : 300   KLM=5

	mov	ecx, DWORD PTR _KLM$[esp+268]
	mov	DWORD PTR $T98310[esp+284], 102		; 00000066H
	mov	DWORD PTR [ecx], 5
	jmp	$L79
$L98308:
	fld	QWORD PTR _DETSK
	fld	QWORD PTR $T98315
	mov	DWORD PTR $T98310[esp+284], 35		; 00000023H

; 34   : C     -----------------------------------------------------
; 35   : 	OV=(DABS(DSQRT(DABS(DSQRT(DETSL**3.0D+00)))) *

	call	__FIfexp
	call	__FIsqrt
	fabs
	call	__FIsqrt
	fabs
	fld	QWORD PTR _DETSL
	fld	QWORD PTR $T98315
	call	__FIfexp
	call	__FIsqrt
	fabs
	call	__FIsqrt
	fabs
	fmulp	ST(1), ST(0)
	fld	QWORD PTR _DETAB
	call	__FIsqrt
	fmul	QWORD PTR _DETAB
	mov	ebx, DWORD PTR _OV$[esp+268]
	mov	eax, DWORD PTR _AZ$[esp+268]
	mov	ecx, DWORD PTR _BZ$[esp+268]
	mov	edx, DWORD PTR _AY$[esp+268]
	fdivp	ST(1), ST(0)
	fstp	QWORD PTR [ebx]
	fld	QWORD PTR [eax]
	fsub	QWORD PTR [ecx]
	fld	QWORD PTR [edx]
	mov	eax, DWORD PTR _BY$[esp+268]
	fsub	QWORD PTR [eax]
	fxch	ST(1)
	fstp	QWORD PTR -136+[esp+272]
	fld	QWORD PTR $T98313
	call	__FIfexp
	mov	ecx, DWORD PTR _AX$[esp+268]
	mov	edx, DWORD PTR _BX$[esp+268]
	fld	QWORD PTR [ecx]
	fsub	QWORD PTR [edx]
	fld	QWORD PTR $T98313
	call	__FIfexp
	faddp	ST(1), ST(0)
	fld	QWORD PTR -136+[esp+272]
	fld	QWORD PTR $T98313
	call	__FIfexp
	faddp	ST(1), ST(0)
	mov	DWORD PTR $T98310[esp+284], 39		; 00000027H
	mov	eax, DWORD PTR _CZ$[esp+268]
	fstp	QWORD PTR -224+[esp+272]

; 36   : 	-	DABS(DSQRT(DABS(DSQRT(DETSK**3.0D+00))))) /
; 37   : 	-	(DETAB*DSQRT(DETAB))
; 38   : C     -----------------------------------------------------
; 39   : 	AB2=(AX-BX)**2 + (AY-BY)**2 + (AZ-BZ)**2

	mov	ecx, DWORD PTR -224+[esp+276]
	mov	edx, DWORD PTR -224+[esp+272]
	mov	DWORD PTR _AB2+4, ecx
	mov	ecx, DWORD PTR _DZ$[esp+268]
	mov	DWORD PTR _AB2, edx
	mov	edx, DWORD PTR _CY$[esp+268]
	fld	QWORD PTR [eax]
	fsub	QWORD PTR [ecx]
	mov	eax, DWORD PTR _DY$[esp+268]
	fld	QWORD PTR [edx]
	fsub	QWORD PTR [eax]
	fxch	ST(1)
	fstp	QWORD PTR -144+[esp+272]
	fld	QWORD PTR $T98313
	call	__FIfexp
	mov	ecx, DWORD PTR _CX$[esp+268]
	mov	edx, DWORD PTR _DX$[esp+268]
	fld	QWORD PTR [ecx]
	fsub	QWORD PTR [edx]
	fld	QWORD PTR $T98313
	call	__FIfexp
	faddp	ST(1), ST(0)
	fld	QWORD PTR -144+[esp+272]
	fld	QWORD PTR $T98313
	call	__FIfexp
	faddp	ST(1), ST(0)
	mov	DWORD PTR $T98310[esp+284], 40		; 00000028H
	mov	eax, DWORD PTR _AX$[esp+268]
	fstp	QWORD PTR -232+[esp+272]

; 40   : 	CD2=(CX-DX)**2 + (CY-DY)**2 + (CZ-DZ)**2

	mov	ecx, DWORD PTR -232+[esp+276]
	mov	edx, DWORD PTR -232+[esp+272]
	mov	DWORD PTR _CD2+4, ecx
	mov	ecx, DWORD PTR _A2$[esp+268]
	mov	DWORD PTR _CD2, edx
	mov	edx, DWORD PTR _BX$[esp+268]
	mov	DWORD PTR $T98310[esp+284], 41		; 00000029H
	fld	QWORD PTR [eax]
	fmul	QWORD PTR [edi]
	fld	QWORD PTR [ecx]
	fmul	QWORD PTR [edx]
	mov	ecx, DWORD PTR _A4$[esp+268]
	mov	edx, DWORD PTR _DX$[esp+268]
	faddp	ST(1), ST(0)
	fstp	QWORD PTR -32+[esp+272]

; 41   :       RX1=(A1*AX+A2*BX)/ALK1

	mov	eax, DWORD PTR _CX$[esp+268]
	fld	QWORD PTR -32+[esp+272]
	fdiv	QWORD PTR _ALK1
	fstp	QWORD PTR _RX1
	fld	QWORD PTR [ecx]
	fmul	QWORD PTR [edx]
	fld	QWORD PTR [esi]
	fmul	QWORD PTR [eax]
	faddp	ST(1), ST(0)
	fstp	QWORD PTR -96+[esp+272]
	mov	DWORD PTR $T98310[esp+284], 42		; 0000002aH

; 42   :       RX2=(A3*CX+A4*DX)/ALK2

	fld	QWORD PTR -96+[esp+272]
	fdiv	QWORD PTR _ALK2
	mov	edx, DWORD PTR _A2$[esp+268]
	mov	eax, DWORD PTR _BY$[esp+268]
	mov	ecx, DWORD PTR _AY$[esp+268]
	fstp	QWORD PTR _RX2
	fld	QWORD PTR [ecx]
	fmul	QWORD PTR [edi]
	fld	QWORD PTR [edx]
	fmul	QWORD PTR [eax]
	mov	DWORD PTR $T98310[esp+284], 43		; 0000002bH
	mov	ecx, DWORD PTR _A4$[esp+268]
	mov	edx, DWORD PTR _DY$[esp+268]
	faddp	ST(1), ST(0)
	fstp	QWORD PTR -24+[esp+272]

; 43   :       RY1=(A1*AY+A2*BY)/ALK1

	fld	QWORD PTR -24+[esp+272]
	fdiv	QWORD PTR _ALK1
	mov	eax, DWORD PTR _CY$[esp+268]
	fstp	QWORD PTR _RY1
	fld	QWORD PTR [ecx]
	fmul	QWORD PTR [edx]
	fld	QWORD PTR [esi]
	fmul	QWORD PTR [eax]
	mov	DWORD PTR $T98310[esp+284], 44		; 0000002cH
	mov	ecx, DWORD PTR _AZ$[esp+268]
	faddp	ST(1), ST(0)
	mov	edx, DWORD PTR _A2$[esp+268]
	fstp	QWORD PTR -88+[esp+272]

; 44   :       RY2=(A3*CY+A4*DY)/ALK2

	fld	QWORD PTR -88+[esp+272]
	fdiv	QWORD PTR _ALK2
	mov	eax, DWORD PTR _BZ$[esp+268]
	fstp	QWORD PTR _RY2
	fld	QWORD PTR [edi]
	fmul	QWORD PTR [ecx]
	fld	QWORD PTR [edx]
	fmul	QWORD PTR [eax]
	mov	DWORD PTR $T98310[esp+284], 45		; 0000002dH
	mov	ecx, DWORD PTR _A4$[esp+268]
	mov	edx, DWORD PTR _DZ$[esp+268]
	faddp	ST(1), ST(0)
	fstp	QWORD PTR -40+[esp+272]

; 45   :       RZ1=(A1*AZ+A2*BZ)/ALK1

	fld	QWORD PTR -40+[esp+272]
	fdiv	QWORD PTR _ALK1
	mov	eax, DWORD PTR _CZ$[esp+268]
	fstp	QWORD PTR _RZ1
	fld	QWORD PTR [ecx]
	fmul	QWORD PTR [edx]
	fld	QWORD PTR [esi]
	fmul	QWORD PTR [eax]
	mov	DWORD PTR $T98310[esp+284], 46		; 0000002eH
	faddp	ST(1), ST(0)
	fstp	QWORD PTR -80+[esp+272]

; 46   :       RZ2=(A3*CZ+A4*DZ)/ALK2

	fld	QWORD PTR -80+[esp+272]
	fdiv	QWORD PTR _ALK2
	fstp	QWORD PTR _RZ2
	fld	QWORD PTR _RZ1
	fsub	QWORD PTR _RZ2
	fld	QWORD PTR _RY1
	fsub	QWORD PTR _RY2
	fld	QWORD PTR _RX1
	fsub	QWORD PTR _RX2
	fxch	ST(2)
	fstp	QWORD PTR -240+[esp+272]
	fstp	QWORD PTR -152+[esp+272]
	fstp	QWORD PTR -160+[esp+272]
	fld	QWORD PTR -240+[esp+272]
	fld	QWORD PTR $T98313
	mov	DWORD PTR $T98310[esp+284], 48		; 00000030H

; 47   : C
; 48   :       DD=((RX1-RX2)**2)+((RY1-RY2)**2)+((RZ1-RZ2)**2)

	call	__FIfexp
	fld	QWORD PTR -152+[esp+272]
	fld	QWORD PTR $T98313
	call	__FIfexp
	faddp	ST(1), ST(0)
	fld	QWORD PTR -160+[esp+272]
	fld	QWORD PTR $T98313
	call	__FIfexp
	faddp	ST(1), ST(0)
	mov	ecx, DWORD PTR _A4$[esp+268]
	mov	edx, DWORD PTR _A2$[esp+268]
	fld	QWORD PTR [ecx]
	fmul	QWORD PTR [esi]
	fxch	ST(1)
	fstp	QWORD PTR _DD
	fstp	QWORD PTR -104+[esp+272]
	fld	QWORD PTR [edx]
	fmul	QWORD PTR [edi]
	mov	DWORD PTR $T98310[esp+284], 50		; 00000032H
	fld	QWORD PTR -184+[esp+272]
	fmul	QWORD PTR _BLK12
	fxch	ST(1)
	fstp	QWORD PTR -48+[esp+272]

; 49   : C
; 50   : 	CKL=(((A1*A2)/ALK1)*((AX-BX)**2+(AY-BY)**2+(AZ-BZ)**2)+((A3*A4

	fmul	QWORD PTR _DD
	fld	QWORD PTR -104+[esp+272]
	fdiv	QWORD PTR _ALK2
	fmul	QWORD PTR -232+[esp+272]
	fxch	ST(1)
	fdiv	QWORD PTR _DETAB
	mov	DWORD PTR $T98310[esp+284], 54		; 00000036H
	faddp	ST(1), ST(0)
	fld	QWORD PTR -48+[esp+272]
	fdiv	QWORD PTR _ALK1
	fmul	QWORD PTR -224+[esp+272]
	fld	QWORD PTR _CCK
	fadd	QWORD PTR _CCL
	fxch	ST(1)
	faddp	ST(2), ST(0)
	fdiv	QWORD PTR $T98313
	fxch	ST(1)
	fstp	QWORD PTR _CKL

; 51   : 	-)/ALK2)*((CX-DX)**2+(CY-DY)**2+(CZ-DZ)**2)+((ALK1*ALK2*BLK12)*DD/D
; 52   :      -ETAB))
; 53   : C
; 54   : 	XK=EXP(-CKL+((CCK+CCL)/2.0D+00))

	fsub	QWORD PTR _CKL
	call	__FIexp
	fstp	QWORD PTR _XK
	mov	DWORD PTR $T98310[esp+284], 56		; 00000038H

; 55   : C     -----------------------------------------------------
; 56   :       OV=OV*XK

	fld	QWORD PTR [ebx]
	fmul	QWORD PTR _XK
	mov	DWORD PTR $T98310[esp+284], 60		; 0000003cH
	mov	ecx, DWORD PTR _DFLAG$[esp+268]
	fstp	QWORD PTR [ebx]

; 57   : C     -----------------------------------------------------
; 58   : C	 DERIVATIVES OF THE OVERLAP:
; 59   : C     -----------------------------------------------------
; 60   : 	IF (DFLAG .EQ. 1) THEN

	cmp	DWORD PTR [ecx], 1
	jne	$L79
	fld	QWORD PTR [ebp]
	fadd	QWORD PTR [esi]
	fld	QWORD PTR _DETSL
	fld	QWORD PTR $T98313
	fxch	ST(2)
	fmul	QWORD PTR $T98316
	fstp	QWORD PTR -64+[esp+272]
	fxch	ST(1)
	call	__FIfexp
	fstp	QWORD PTR -168+[esp+272]
	fld	QWORD PTR [ebp]
	fmul	QWORD PTR -248+[esp+272]
	fld	QWORD PTR _ALK1
	fld	QWORD PTR $T98313
	fxch	ST(2)
	fstp	QWORD PTR -72+[esp+272]
	fxch	ST(1)
	call	__FIfexp
	fstp	QWORD PTR -216+[esp+272]
	fld	QWORD PTR _ALK2
	fadd	QWORD PTR _BLK12
	fld	QWORD PTR _DETAB
	fld	QWORD PTR $T98313
	fxch	ST(2)
	fstp	QWORD PTR -56+[esp+272]
	fxch	ST(1)
	call	__FIfexp
	fstp	QWORD PTR -176+[esp+272]
	fld	QWORD PTR _DETSK
	fld	QWORD PTR $T98317
	call	__FIfexp
	fstp	QWORD PTR -256+[esp+272]
	fld	QWORD PTR _DETAB
	fld	QWORD PTR $T98318
	call	__FIfexp
	fld	QWORD PTR _DETSL
	fld	QWORD PTR $T98319
	call	__FIfexp
	fmulp	ST(1), ST(0)
	fstp	QWORD PTR -120+[esp+272]
	fld	QWORD PTR -120+[esp+272]
	fmul	QWORD PTR $T98316
	fld	QWORD PTR _DETSL
	fld	QWORD PTR $T98317
	fxch	ST(2)
	fstp	QWORD PTR -112+[esp+272]
	fxch	ST(1)
	call	__FIfexp
	fstp	QWORD PTR -200+[esp+272]
	fld	QWORD PTR _DETAB
	fld	QWORD PTR $T98320
	call	__FIfexp
	fmul	QWORD PTR $T98313
	mov	DWORD PTR $T98310[esp+284], 62		; 0000003eH
	fld	QWORD PTR [ebp]
	fmul	QWORD PTR [esi]
	fxch	ST(1)
	fstp	QWORD PTR -208+[esp+272]

; 61   : C
; 62   :       DOVDA1=OV*(-(AB2*A2/ALK1)+AB2*A1*A2/ALK1**2.0D+00+ALK1*ALK2*BLK12*

	fmul	QWORD PTR [edi]
	fld	QWORD PTR [esi]
	fmul	QWORD PTR -72+[esp+272]
	fxch	ST(1)
	fmul	QWORD PTR -64+[esp+272]
	fxch	ST(1)
	fmul	QWORD PTR $T98314
	fxch	ST(1)
	fmul	QWORD PTR -248+[esp+272]
	fxch	ST(1)
	fdiv	QWORD PTR _DETSL
	mov	ecx, DWORD PTR _AY$[esp+268]
	fxch	ST(1)
	fmul	QWORD PTR $T98321
	fld	QWORD PTR [ecx]
	fdiv	QWORD PTR _ALK1
	fxch	ST(1)
	fdiv	QWORD PTR -168+[esp+272]
	mov	edx, DWORD PTR _AX$[esp+268]
	faddp	ST(2), ST(0)
	fld	QWORD PTR -24+[esp+272]
	fdiv	QWORD PTR -216+[esp+272]
	fsubp	ST(1), ST(0)
	fxch	ST(1)
	fdiv	QWORD PTR $T98313
	fxch	ST(1)
	fmul	QWORD PTR -152+[esp+272]
	fld	QWORD PTR [edx]
	fdiv	QWORD PTR _ALK1
	mov	eax, DWORD PTR _AZ$[esp+268]
	fld	QWORD PTR -32+[esp+272]
	fdiv	QWORD PTR -216+[esp+272]
	fsubp	ST(1), ST(0)
	fld	QWORD PTR [eax]
	fdiv	QWORD PTR _ALK1
	fxch	ST(1)
	fmul	QWORD PTR -160+[esp+272]
	fld	QWORD PTR -40+[esp+272]
	fdiv	QWORD PTR -216+[esp+272]
	fxch	ST(1)
	faddp	ST(3), ST(0)
	fsubp	ST(1), ST(0)
	fmul	QWORD PTR -240+[esp+272]
	faddp	ST(1), ST(0)
	fld	QWORD PTR _ALK2
	fmul	QWORD PTR _ALK1
	fxch	ST(1)
	fmul	QWORD PTR $T98313
	fxch	ST(1)
	fmul	QWORD PTR -56+[esp+272]
	fxch	ST(1)
	fmul	QWORD PTR _ALK2
	fxch	ST(1)
	fmul	QWORD PTR _BLK12
	fxch	ST(1)
	fmul	QWORD PTR _ALK1
	fxch	ST(1)
	fmul	QWORD PTR _DD
	fxch	ST(1)
	fmul	QWORD PTR _BLK12
	fxch	ST(1)
	fdiv	QWORD PTR -176+[esp+272]
	fxch	ST(1)
	fdiv	QWORD PTR _DETAB
	fsubp	ST(2), ST(0)
	mov	ecx, DWORD PTR _A2$[esp+268]
	faddp	ST(1), ST(0)
	fld	QWORD PTR _ALK2
	fmul	QWORD PTR _BLK12
	fld	QWORD PTR [ecx]
	fmul	QWORD PTR -224+[esp+272]
	fxch	ST(1)
	fmul	QWORD PTR _DD
	fxch	ST(1)
	fdiv	QWORD PTR _ALK1
	fxch	ST(1)
	fdiv	QWORD PTR _DETAB
	fsubp	ST(2), ST(0)
	fld	QWORD PTR -48+[esp+272]
	fmul	QWORD PTR -224+[esp+272]
	fxch	ST(1)
	fsubp	ST(2), ST(0)
	fdiv	QWORD PTR -216+[esp+272]
	faddp	ST(1), ST(0)
	fld	QWORD PTR -56+[esp+272]
	fmul	QWORD PTR -256+[esp+272]
	fxch	ST(1)
	fmul	QWORD PTR [ebx]
	fxch	ST(1)
	fmul	QWORD PTR -200+[esp+272]
	fld	QWORD PTR -64+[esp+272]
	fmul	QWORD PTR -256+[esp+272]
	fxch	ST(1)
	fmul	QWORD PTR _XK
	fxch	ST(1)
	fmul	QWORD PTR _XK
	fxch	ST(1)
	fmul	QWORD PTR $T98315
	fxch	ST(1)
	fmul	QWORD PTR $T98315
	fxch	ST(1)
	fdiv	QWORD PTR -208+[esp+272]
	fsubp	ST(2), ST(0)
	fdiv	QWORD PTR -112+[esp+272]
	faddp	ST(1), ST(0)
	fld	QWORD PTR [ebp]
	fadd	QWORD PTR [edi]
	fxch	ST(1)
	fstp	QWORD PTR _DOV
	fmul	QWORD PTR $T98316
	fld	QWORD PTR _ALK2
	fld	QWORD PTR $T98313
	fxch	ST(2)
	fstp	QWORD PTR -32+[esp+272]
	fxch	ST(1)
	call	__FIfexp
	fstp	QWORD PTR -40+[esp+272]
	fld	QWORD PTR _ALK1
	fadd	QWORD PTR _BLK12
	fstp	QWORD PTR -24+[esp+272]
	mov	DWORD PTR $T98310[esp+284], 74		; 0000004aH

; 63   : 	-(ALK2+BLK12)*DD/DETAB**2.0D+00-ALK2*BLK12*DD/DETAB+(-8.0D+00*A1*A3
; 64   :      -*A5*(4.0D+00*A3+4.0D+00*A5)*(AZ-CZ)**2.0D+00/DETSL**2.0D+00+8.0D
; 65   :      -+00*A3*A5*(AZ-CZ)**2.0D+00/DETSL)/2.0D+00-ALK1*ALK2*BLK12*(2.0D+
; 66   :      -00*(AX/ALK1-(AX*A1+A2*BX)/ALK1**2.0D+00)*(RX1-RX2)+2.0D+00*(AY/ALK
; 67   :      -1-(AY*A1+A2*BY)/ALK1**2.0D+00)*(RY1-RY2)+2.0D+00*(AZ/ALK1-(AZ*A1+A
; 68   :      -2*BZ)/ALK1**2.0D+00)*(RZ1-RZ2))/DETAB)+3.0D+00*(4.0D+00*A3+4.0D+00
; 69   :      -*A5)*DETSK**(3.0D+00/4.0D+00)*XK/(4.0D+00*DETAB**(3.0D+00/2.0D+0
; 70   :      -0)*DETSL**(1.0D+00/4.0D+00))-3.0D+00*(ALK2+BLK12)*DETSK**(3.0D
; 71   :      -+00/4.0D+00)*DETSL**(3.0D+00/4.0D+00)*XK/(2.0D+00*DETAB**(5.0D+0
; 72   :      -0/2.0D+00))
; 73   : C
; 74   : 	DOVDA2=OV*(-(A4*CD2/ALK2)+A3*A4*CD2/ALK2**2.0D+00+ALK1*ALK2*BLK12*

	fld	QWORD PTR [ebp]
	fmul	QWORD PTR [esi]
	fld	QWORD PTR [edi]
	fmul	QWORD PTR -72+[esp+272]
	fxch	ST(1)
	fmul	QWORD PTR [edi]
	fxch	ST(1)
	fmul	QWORD PTR $T98314
	fxch	ST(1)
	fmul	QWORD PTR -248+[esp+272]
	fxch	ST(1)
	fdiv	QWORD PTR _DETSL
	fxch	ST(1)
	fmul	QWORD PTR -32+[esp+272]
	fld	QWORD PTR -80+[esp+272]
	fdiv	QWORD PTR -40+[esp+272]
	mov	ecx, DWORD PTR _CZ$[esp+268]
	fxch	ST(1)
	fmul	QWORD PTR $T98321
	fld	QWORD PTR [ecx]
	fdiv	QWORD PTR _ALK2
	fxch	ST(1)
	fdiv	QWORD PTR -168+[esp+272]
	mov	edx, DWORD PTR _CY$[esp+268]
	faddp	ST(3), ST(0)
	fsubp	ST(1), ST(0)
	fld	QWORD PTR -88+[esp+272]
	fdiv	QWORD PTR -40+[esp+272]
	fxch	ST(1)
	fmul	QWORD PTR -240+[esp+272]
	fxch	ST(2)
	fdiv	QWORD PTR $T98313
	fld	QWORD PTR [edx]
	fdiv	QWORD PTR _ALK2
	mov	eax, DWORD PTR _CX$[esp+268]
	fsubp	ST(2), ST(0)
	fld	QWORD PTR -96+[esp+272]
	fdiv	QWORD PTR -40+[esp+272]
	fxch	ST(2)
	fmul	QWORD PTR -152+[esp+272]
	fld	QWORD PTR [eax]
	fdiv	QWORD PTR _ALK2
	fxch	ST(1)
	faddp	ST(4), ST(0)
	fsubp	ST(2), ST(0)
	fxch	ST(1)
	fmul	QWORD PTR -160+[esp+272]
	faddp	ST(2), ST(0)
	fld	QWORD PTR _ALK2
	fmul	QWORD PTR _ALK1
	fxch	ST(2)
	fmul	QWORD PTR $T98313
	fxch	ST(2)
	fmul	QWORD PTR -24+[esp+272]
	fxch	ST(2)
	fmul	QWORD PTR _ALK2
	fxch	ST(2)
	fmul	QWORD PTR _BLK12
	fxch	ST(2)
	fmul	QWORD PTR _ALK1
	fld	QWORD PTR _ALK1
	fmul	QWORD PTR _BLK12
	fxch	ST(3)
	fmul	QWORD PTR _DD
	fxch	ST(3)
	fmul	QWORD PTR _DD
	fxch	ST(3)
	fdiv	QWORD PTR -176+[esp+272]
	mov	ecx, DWORD PTR _A4$[esp+268]
	fxch	ST(1)
	fmul	QWORD PTR _BLK12
	fxch	ST(3)
	fdiv	QWORD PTR _DETAB
	fxch	ST(3)
	fdiv	QWORD PTR _DETAB
	fsubp	ST(2), ST(0)
	fld	QWORD PTR [ecx]
	fmul	QWORD PTR -232+[esp+272]
	fxch	ST(1)
	faddp	ST(2), ST(0)
	fld	QWORD PTR -104+[esp+272]
	fmul	QWORD PTR -232+[esp+272]
	fxch	ST(3)
	fsubp	ST(2), ST(0)
	fdiv	QWORD PTR _ALK2
	fsubp	ST(1), ST(0)
	fxch	ST(1)
	fdiv	QWORD PTR -40+[esp+272]
	faddp	ST(1), ST(0)
	fld	QWORD PTR -24+[esp+272]
	fmul	QWORD PTR -256+[esp+272]
	fld	QWORD PTR -32+[esp+272]
	fmul	QWORD PTR -256+[esp+272]
	fxch	ST(1)
	fmul	QWORD PTR -200+[esp+272]
	fxch	ST(1)
	fmul	QWORD PTR _XK
	fxch	ST(1)
	fmul	QWORD PTR _XK
	fxch	ST(1)
	fmul	QWORD PTR $T98315
	fxch	ST(1)
	fmul	QWORD PTR $T98315
	fxch	ST(1)
	fdiv	QWORD PTR -112+[esp+272]
	fxch	ST(2)
	fmul	QWORD PTR [ebx]
	fxch	ST(1)
	fdiv	QWORD PTR -208+[esp+272]
	fsubp	ST(1), ST(0)
	fld	QWORD PTR [esi]
	fadd	QWORD PTR [edi]
	fxch	ST(2)
	faddp	ST(1), ST(0)
	fld	QWORD PTR [ebp]
	fmul	QWORD PTR [esi]
	fxch	ST(1)
	fstp	QWORD PTR _DOV+8
	fmul	QWORD PTR [edi]
	fxch	ST(1)
	fmul	QWORD PTR $T98313
	fld	QWORD PTR [esi]
	fmul	QWORD PTR [edi]
	fxch	ST(1)
	fstp	QWORD PTR -24+[esp+272]
	fmul	QWORD PTR -248+[esp+272]
	fxch	ST(1)
	fmul	QWORD PTR -24+[esp+272]
	fxch	ST(1)
	fmul	QWORD PTR $T98314
	fxch	ST(1)
	fmul	QWORD PTR -248+[esp+272]
	fxch	ST(1)
	fdiv	QWORD PTR _DETSL
	mov	DWORD PTR $T98310[esp+284], 86		; 00000056H

; 75   : 	-(ALK1+BLK12)*DD/DETAB**2.0D+00-ALK1*BLK12*DD/DETAB+(-8.0D+00*A1*A3
; 76   :      -*A5*(4.0D+00*A1+4.0D+00*A5)*(AZ-CZ)**2.0D+00/DETSL**2.0D+00+8.0D
; 77   :      -+00*A1*A5*(AZ-CZ)**2.0D+00/DETSL)/2.0D+00-ALK1*ALK2*BLK12*(2.0D+
; 78   :      -00*(-(CX/ALK2)+(A3*CX+A4*DX)/ALK2**2.0D+00)*(RX1-RX2)+2.0D+00*(-(C
; 79   :      -Y/ALK2)+(A3*CY+A4*DY)/ALK2**2.0D+00)*(RY1-RY2)+2.0D+00*(-(CZ/ALK2)
; 80   :      -+(A3*CZ+A4*DZ)/ALK2**2.0D+00)*(RZ1-RZ2))/DETAB)+3.0D+00*(4.0D+00*A
; 81   :      -1+4.0D+00*A5)*DETSK**(3.0D+00/4.0D+00)*XK/(4.0D+00*DETAB**(3.0D+
; 82   :      -00/2.0D+00)*DETSL**(1.0D+00/4.0D+00))-3.0D+00*(ALK1+BLK12)*DETSK
; 83   :      -**(3.0D+00/4.0D+00)*DETSL**(3.0D+00/4.0D+00)*XK/(2.0D+00*DETAB
; 84   :      -**(5.0D+00/2.0D+00))
; 85   : C
; 86   : 	DOVDA5=(ALK1*ALK2*(ALK1+ALK2)*BLK12*DD/DETAB**2.0D+00-ALK1*ALK2*DD

	fxch	ST(1)
	fmul	QWORD PTR $T98322
	fld	QWORD PTR _ALK2
	fmul	QWORD PTR _ALK1
	fxch	ST(1)
	fdiv	QWORD PTR -168+[esp+272]
	faddp	ST(2), ST(0)
	fmul	QWORD PTR -192+[esp+272]
	fld	QWORD PTR -184+[esp+272]
	fmul	QWORD PTR _DD
	fxch	ST(2)
	fdiv	QWORD PTR $T98313
	fxch	ST(1)
	fmul	QWORD PTR _BLK12
	fxch	ST(2)
	fdiv	QWORD PTR _DETAB
	fxch	ST(2)
	fmul	QWORD PTR _DD
	fld	QWORD PTR -192+[esp+272]
	fmul	QWORD PTR -256+[esp+272]
	fxch	ST(1)
	fdiv	QWORD PTR -176+[esp+272]
	faddp	ST(2), ST(0)
	fmul	QWORD PTR -200+[esp+272]
	fld	QWORD PTR -24+[esp+272]
	fmul	QWORD PTR -256+[esp+272]
	fxch	ST(3)
	fsubp	ST(2), ST(0)
	fmul	QWORD PTR _XK
	fld	QWORD PTR -120+[esp+272]
	fmul	QWORD PTR $T98313
	fxch	ST(3)
	fmul	QWORD PTR _XK
	fxch	ST(1)
	fmul	QWORD PTR $T98315
	fxch	ST(1)
	fmul	QWORD PTR $T98315
	fxch	ST(1)
	fdiv	QWORD PTR -208+[esp+272]
	fxch	ST(2)
	fmul	QWORD PTR [ebx]
	fld	QWORD PTR [ebp]
	fmul	QWORD PTR [esi]
	fxch	ST(3)
	fsubp	ST(1), ST(0)
	fxch	ST(3)
	fdivp	ST(1), ST(0)
	faddp	ST(2), ST(0)
	fmul	QWORD PTR [edi]
	fmul	QWORD PTR -128+[esp+272]
	fxch	ST(1)
	fstp	QWORD PTR _DOV+16
	fstp	QWORD PTR -24+[esp+272]
	mov	DWORD PTR $T98310[esp+284], 94		; 0000005eH

; 87   : 	-/DETAB+(-16.0D+00*A1*A3*(2.0D+00*A1+2.0D+00*A3)*A5*(AZ-CZ)**2.0D+0
; 88   :      -0/DETSL**2.0D+00+8.0D+00*A1*A3*(AZ-CZ)**2.0D+00/DETSL)/2.0D+00
; 89   :      -)*OV+3.0D+00*(2.0D+00*A1+2.0D+00*A3)*DETSK**(3.0D+00/4.0D+00)*XK
; 90   :      -/(2.0D+00*DETAB**(3.0D+00/2.0D+00)*DETSL**(1.0D+00/4.0D+00))-3.0
; 91   :      -D+00*(ALK1+ALK2)*DETSK**(3.0D+00/4.0D+00)*DETSL**(3.0D+00/4.0D
; 92   :      -+00)*XK/(2.0D+00*DETAB**(5.0D+00/2.0D+00))
; 93   : C
; 94   : 	DOVDAZ=OV*(-2.0D+00*A1*A2*(AZ-BZ)/ALK1+8.0D+00*A1*A3*A5*(AZ-CZ)/DE

	fld	QWORD PTR [edi]
	fmul	QWORD PTR _ALK2
	mov	ecx, DWORD PTR _A2$[esp+268]
	fmul	QWORD PTR -240+[esp+272]
	fld	QWORD PTR [ecx]
	fmul	QWORD PTR [edi]
	fxch	ST(1)
	fmul	QWORD PTR _BLK12
	fld	QWORD PTR -24+[esp+272]
	fmul	QWORD PTR $T98314
	fxch	ST(2)
	fmul	QWORD PTR -136+[esp+272]
	fxch	ST(2)
	fdiv	QWORD PTR _DETSL
	fxch	ST(1)
	fmul	QWORD PTR $T98323
	fxch	ST(2)
	fmul	QWORD PTR $T98313
	fxch	ST(2)
	fdiv	QWORD PTR _DETAB
	mov	DWORD PTR $T98310[esp+284], 97		; 00000061H
	fxch	ST(2)
	fdiv	QWORD PTR _ALK1
	fsubp	ST(2), ST(0)
	fld	QWORD PTR [esi]
	fmul	QWORD PTR _ALK1
	fxch	ST(1)
	faddp	ST(2), ST(0)
	mov	edx, DWORD PTR _A4$[esp+268]
	fmul	QWORD PTR -240+[esp+272]
	fld	QWORD PTR [edx]
	fmul	QWORD PTR [esi]
	fxch	ST(1)
	fmul	QWORD PTR _BLK12
	fxch	ST(1)
	fmul	QWORD PTR -144+[esp+272]
	fxch	ST(1)
	fmul	QWORD PTR $T98313
	fld	QWORD PTR -24+[esp+272]
	fmul	QWORD PTR $T98314
	fxch	ST(1)
	fdiv	QWORD PTR _DETAB
	fxch	ST(3)
	fmul	QWORD PTR [ebx]
	fxch	ST(1)
	fdiv	QWORD PTR _DETSL
	fxch	ST(2)
	fmul	QWORD PTR $T98313
	fxch	ST(1)
	fstp	QWORD PTR _DOV+24

; 95   : 	-TSL-2.0D+00*ALK2*A1*BLK12*(RZ1-RZ2)/DETAB)
; 96   : C
; 97   : 	DOVDBZ=OV*(-8.0D+00*A1*A3*A5*(AZ-CZ)/DETSL-2.0D+00*A3*A4*(CZ-DZ)

	fdiv	QWORD PTR _ALK2
	mov	DWORD PTR $T98310[esp+284], 110		; 0000006eH
	fsubp	ST(2), ST(0)
	fsubp	ST(1), ST(0)
	fmul	QWORD PTR [ebx]
	fstp	QWORD PTR _DOV+32
$L79:

; 103  : C     -----------------------------------------------------
; 104  : 1500  END IF
; 105  : 
; 106  : 
; 107  : c	WRITE(6,*) OV
; 108  : 
; 109  : 	RETURN                                                            
; 110  :       END

	mov	DWORD PTR $T98310[esp+284], 110		; 0000006eH
	call	___FHunlink@0
	pop	ebp
	pop	edi
	pop	esi
	pop	ebx
	add	esp, 256				; 00000100H
	ret	88					; 00000058H
	npad	4
_OVERLP@88 ENDP
_TEXT	ENDS
END
