	TITLE	C:\users\dgilmore\Graduate.School\Adamowicz.Research\H2ANL\F0.FOR
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
FLAT	GROUP _DATA
	ASSUME	CS: FLAT, DS: FLAT, SS: FLAT
endif
EXTRN	__FFljj:NEAR
_DATA	SEGMENT
_C:\users\dgilmore\Graduate.School\Adamowicz.Research\H2ANL\F0.FOR DB 'C:'
	DB	'\users\dgilmore\Graduate.School\Adamowicz.Research\H2ANL\F0.F'
	DB	'OR', 00H
_DATA	ENDS
PUBLIC	_F0@4
EXTRN	__FIsqrt:NEAR
EXTRN	__fltused:NEAR
_BSS	SEGMENT
_R	DQ	01H DUP (?)
_R2	DQ	01H DUP (?)
_PI	DQ	01H DUP (?)
_R3	DQ	01H DUP (?)
_BSS	ENDS
_CONST	SEGMENT
$T98309	DQ	040305b98c7e28241r		; 16.3578
$T98310	DQ	00000000000000000r		; 0
$T98311	DQ	04000000000000000r		; 2
$T98312	DQ	0400921fb54442d28r		; 3.14159
$T98313	DQ	03ff921fb54442d28r		; 1.5708
$T98314	DQ	03f0d9b282dc36e12r		; 5.64691e-005
$T98315	DQ	03f48da32e37f14e9r		; 0.000758433
$T98316	DQ	03f7f88563ebdb857r		; 0.00769838
$T98317	DQ	03fb01c78cb5ca40ar		; 0.0629344
$T98318	DQ	03fcb4c7959a1f6a4r		; 0.213271
$T98319	DQ	03ff0000000000000r		; 1
$T98320	DQ	03f12e1a07588fc71r		; 7.20267e-005
$T98321	DQ	03f4f4f8e903a27b1r		; 0.000955529
$T98322	DQ	03f84c5ef42558089r		; 0.0101432
$T98323	DQ	03fb2e7fbe90c4bfdr		; 0.0738523
$T98324	DQ	03fd5a92bbb78443dr		; 0.33845
$T98325	DQ	03fec287351b9e49br		; 0.879938
_CONST	ENDS
_TEXT	SEGMENT
; File C:\users\dgilmore\Graduate.School\Adamowicz.Research\H2ANL\F0.FOR
_F0$ = -8
_Z$ = 8
_F0@4	PROC NEAR

; 2    :       FUNCTION F0(Z)

	mov	ecx, DWORD PTR _Z$[esp-4]
	sub	esp, 8

; 3    :       IMPLICIT REAL*8(A-H,O-Z)
; 4    : C ***
; 5    :       IF(Z-16.3578D+00) 1,1,2

	fld	QWORD PTR [ecx]
	fsub	QWORD PTR $T98309
	fcomp	QWORD PTR $T98310
	fwait
	fnstsw	ax
	fld	QWORD PTR [ecx]
	test	ah, 65					; 00000041H

; 6    :     1 CONTINUE
; 7    :       F0=DSQRT((((((0.564691197D-04 *Z
; 8    :      *             +0.758433197D-03)*Z
; 9    :      *             +0.769838037D-02)*Z
; 10   :      *             +0.629344460D-01)*Z
; 11   :      *             +0.213271302D-00)*Z
; 12   :      *             +1.0D+00)/
; 13   :      *        ((((((0.720266520D-04 *Z
; 14   :      *             +0.955528842D-03)*Z
; 15   :      *             +0.101431553D-01)*Z
; 16   :      *             +0.738522953D-01)*Z
; 17   :      *             +0.338450368D+00)*Z
; 18   :      *             +0.879937801D+00)*Z
; 19   :      *             +1.0D+00))
; 20   :       GO TO 3
; 21   :     2 CONTINUE
; 22   :       R=2.D+00*Z

	jne	SHORT $L98308
	fmul	QWORD PTR $T98311
	fstp	QWORD PTR _R

; 23   :       R2=DSQRT(R)

	fld	QWORD PTR _R
	call	__FIsqrt
	fstp	QWORD PTR _R2

; 24   :       PI=3.1415926535898D+00

	fld	QWORD PTR $T98313
	mov	DWORD PTR _PI, 1413754152		; 54442d28H
	mov	DWORD PTR _PI+4, 1074340347		; 400921fbH

; 25   :       R3=DSQRT(PI/2.D+00)

	call	__FIsqrt
	fstp	QWORD PTR _R3

; 26   :       F0=R3/R2

	fld	QWORD PTR _R3
	fdiv	QWORD PTR _R2
	fstp	QWORD PTR _F0$[esp+8]
	fld	QWORD PTR _F0$[esp+8]

; 27   :     3 CONTINUE
; 28   :       RETURN
; 29   :       END

	add	esp, 8
	ret	4
$L98308:
	fmul	QWORD PTR $T98314
	fadd	QWORD PTR $T98315
	fmul	QWORD PTR [ecx]
	fadd	QWORD PTR $T98316
	fmul	QWORD PTR [ecx]
	fadd	QWORD PTR $T98317
	fmul	QWORD PTR [ecx]
	fadd	QWORD PTR $T98318
	fmul	QWORD PTR [ecx]
	fadd	QWORD PTR $T98319
	fld	QWORD PTR [ecx]
	fmul	QWORD PTR $T98320
	fadd	QWORD PTR $T98321
	fmul	QWORD PTR [ecx]
	fadd	QWORD PTR $T98322
	fmul	QWORD PTR [ecx]
	fadd	QWORD PTR $T98323
	fmul	QWORD PTR [ecx]
	fadd	QWORD PTR $T98324
	fmul	QWORD PTR [ecx]
	fadd	QWORD PTR $T98325
	fmul	QWORD PTR [ecx]
	fadd	QWORD PTR $T98319
	fdivp	ST(1), ST(0)
	call	__FIsqrt
	fstp	QWORD PTR _F0$[esp+8]
	fld	QWORD PTR _F0$[esp+8]
	add	esp, 8
	ret	4
	npad	7
	npad	7
_F0@4	ENDP
_TEXT	ENDS
END
