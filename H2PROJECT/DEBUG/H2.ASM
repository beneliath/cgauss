	TITLE	C:\users\dgilmore\Graduate.School\Adamowicz.Research\H2ANL\H2.FOR
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
PUBLIC	_FINALGRAD
EXTRN	__FFsSTOP:NEAR
EXTRN	_f90_open:NEAR
EXTRN	_SEED@4:NEAR
EXTRN	_f90_write:NEAR
EXTRN	_lwrite_ch:NEAR
EXTRN	_f90_endwrite:NEAR
EXTRN	__FFljj:NEAR
EXTRN	_GUESS@48:NEAR
EXTRN	_GETDAT@12:NEAR
EXTRN	_f90_read:NEAR
EXTRN	_lread_ch:NEAR
EXTRN	_f90_close:NEAR
EXTRN	_f90_endread:NEAR
EXTRN	_GETTIM@16:NEAR
EXTRN	_f90_rewind:NEAR
EXTRN	_HAMILTONP@92:NEAR
EXTRN	_RANDOM@4:NEAR
EXTRN	_write_i2:NEAR
EXTRN	__Deallocate_i:NEAR
EXTRN	_lread_i:NEAR
EXTRN	_BIOUT@8:NEAR
EXTRN	_lread_d:NEAR
EXTRN	_write_i:NEAR
EXTRN	_write_d:NEAR
EXTRN	_make_info:NEAR
EXTRN	_NGRADOPT@48:NEAR
EXTRN	_lwrite_d:NEAR
EXTRN	__Allocate_i:NEAR
EXTRN	_HAMILTON@80:NEAR
EXTRN	_GRADOPT@48:NEAR
EXTRN	_lwrite_i:NEAR
EXTRN	__FIstrcpy@24:NEAR
_DATA	SEGMENT
_FINALGRAD DQ	01H DUP (?)
	ORG $-8
_C:\users\dgilmore\Graduate.School\Adamowicz.Research\H2ANL\H2.FOR DB 'C:'
	DB	'\users\dgilmore\Graduate.School\Adamowicz.Research\H2ANL\H2.F'
	DB	'OR', 00H
	ORG $+2
_D:\MSDEV\INCLUDE/msflib.fd DB 'D:\MSDEV\INCLUDE/msflib.fd', 00H
_DATA	ENDS
PUBLIC	_main
EXTRN	__FIclne@16:NEAR
EXTRN	__FIcleq@16:NEAR
EXTRN	__FIlog:NEAR
EXTRN	__FIsqrt:NEAR
EXTRN	__ftol:NEAR
EXTRN	__fltused:NEAR
EXTRN	__FIfexp:NEAR
_BSS	SEGMENT
_XRX1	DQ	01H DUP (?)
_XRY1	DQ	01H DUP (?)
_RESTARTOUTPUT DB 042H DUP (?)
	ALIGN	8

_XRX2	DQ	01H DUP (?)
_XRY2	DQ	01H DUP (?)
_OPTIMIZEDWAVE DB 042H DUP (?)
	ALIGN	8

_BD	DQ	01H DUP (?)
_ISYMM	DD	01H DUP (?)
	ALIGN	8

_FMAX	DQ	01H DUP (?)
_RESTART DD	01H DUP (?)
	ALIGN	8

_FMIN	DQ	01H DUP (?)
_STGROW	DD	01H DUP (?)
_ITEMP	DD	01H DUP (?)
_ZMAX	DQ	01H DUP (?)
_PREOPT	DD	01H DUP (?)
	ALIGN	8

_TENRGY	DQ	01H DUP (?)
_BL0	DQ	01H DUP (?)
_ISTEND	DD	01H DUP (?)
_INDEX	DD	01H DUP (?)
_BL1	DQ	01H DUP (?)
_STGROWPER DQ	01H DUP (?)
_BLURBEGIN DQ	01H DUP (?)
_BLUREND DQ	01H DUP (?)
_BETAMAX DQ	01H DUP (?)
_BETAMIN DQ	01H DUP (?)
_TMPYEAR DW	01H DUP (?)
	ALIGN	4

_IVIEW	DD	01H DUP (?)
_TMPMONTH DW	01H DUP (?)
	ALIGN	8

_ZNUC	DQ	01H DUP (?)
_TMPDAY	DW	01H DUP (?)
	ALIGN	4

_RSTR	DB	042H DUP (?)
	ALIGN	8

_XB12	DQ	01H DUP (?)
_A	DQ	01H DUP (?)
_B	DQ	01H DUP (?)
_TMPHOUR DW	01H DUP (?)
	ALIGN	4

_TMPMINUTE DW	01H DUP (?)
	ALIGN	4

_TEMP	DD	01H DUP (?)
_TMPSECOND DW	01H DUP (?)
	ALIGN	8

_TEMP1	DQ	01H DUP (?)
_TMPHUND DW	01H DUP (?)
	ALIGN	8

_TA1	DQ	01H DUP (?)
_KA	DD	01H DUP (?)
	ALIGN	8

_TA2	DQ	01H DUP (?)
_TB12	DQ	01H DUP (?)
_TRZ1	DQ	01H DUP (?)
_TRZ2	DQ	01H DUP (?)
_TEMP2	DQ	01H DUP (?)
_TSTORE1T DQ	00H DUP (?)
_WAVEOPT DD	01H DUP (?)
	ALIGN	8

_TSTORE2T DQ	00H DUP (?)
_TSTORE3T DQ	00H DUP (?)
_TOF1	DQ	01H DUP (?)
_RZ1STORET DQ	00H DUP (?)
_XAA1	DQ	01H DUP (?)
_METHODNAME DB	042H DUP (?)
	ALIGN	4

_MMAX	DD	01H DUP (?)
_RZ2STORET DQ	00H DUP (?)
_EVECSTORET DQ	00H DUP (?)
_RX1	DQ	00H DUP (?)
_TOF2	DQ	01H DUP (?)
_THENISYMM DQ	01H DUP (?)
_EBEST	DQ	01H DUP (?)
_XAA2	DQ	01H DUP (?)
_RY1	DQ	00H DUP (?)
_ELSEISYMM DQ	01H DUP (?)
_XRZ1	DQ	01H DUP (?)
_DFLAG	DD	01H DUP (?)
	ALIGN	8

_RZ1	DQ	00H DUP (?)
_RELOADFLAG DD	01H DUP (?)
	ALIGN	8

_XRZ2	DQ	01H DUP (?)
_GENNAME DB	042H DUP (?)
	ALIGN	8

_AA2	DQ	00H DUP (?)
_RX2	DQ	00H DUP (?)
_ENDBASIS DD	01H DUP (?)
	ALIGN	8

_THENDFLAG DQ	01H DUP (?)
_RY2	DQ	00H DUP (?)
_TSTORE1 DQ	00H DUP (?)
_RZ2	DQ	00H DUP (?)
_TSTORE2 DQ	00H DUP (?)
_EMIN	DQ	01H DUP (?)
_B12	DQ	00H DUP (?)
_TSTORE3 DQ	00H DUP (?)
_EVEC	DQ	00H DUP (?)
_RZ1STORE DQ	00H DUP (?)
_S	DQ	00H DUP (?)
_RZ2STORE DQ	00H DUP (?)
_XC	DQ	01H DUP (?)
_H	DQ	00H DUP (?)
_EVECSTORE DQ	00H DUP (?)
_CO	DQ	00H DUP (?)
_GRAD	DQ	00H DUP (?)
_T1	DQ	00H DUP (?)
_T2	DQ	00H DUP (?)
_IMODE	DD	01H DUP (?)
_STARTBASIS DD	01H DUP (?)
_DEFLAG	DQ	01H DUP (?)
_T3	DQ	00H DUP (?)
_METHODOPT DD	01H DUP (?)
	ALIGN	8

_ENER	DQ	01H DUP (?)
_LOOP	DD	01H DUP (?)
_METHODGEN DD	01H DUP (?)
_ENRGY	DQ	00H DUP (?)
_M	DD	01H DUP (?)
	ALIGN	8

_AA1	DQ	00H DUP (?)
_ERROR	DD	01H DUP (?)
_FILENAME DB	042H DUP (?)
	ALIGN	4

_RESTARTNAME DB	042H DUP (?)
_BSS	ENDS
_DATA	SEGMENT
	ORG $+1
_7TSTORE3T DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H
_7TSTORE2 DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H
	ORG $+4
_7T3	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H
_7RY2	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H
_7RX1	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H
_7RZ2STORET DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H
_7RZ1STORET DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H
_7TSTORE2T DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H
_7TSTORE1T DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H
_7EVECSTORET DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H
_7RZ2STORE DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H
	ORG $+4
_7RZ1STORE DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H
	ORG $+4
_7TSTORE3 DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H
	ORG $+4
_7TSTORE1 DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H
	ORG $+4
_7EVECSTORE DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H
	ORG $+4
_7T2	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H
_7T1	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H
_7GRAD	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H
_7ENRGY	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H
_7CO	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H
_7S	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H
	ORG $+4
_7EVEC	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H
_7B12	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H
_7RZ2	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H
_7RX2	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H
_7AA2	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H
_7RZ1	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H
_7RY1	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H
_7AA1	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H
_7H	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H
_.SG62_	DB	'H2_PROJECT.OUT'
	ORG $+2
_.SG64_	DB	'RESTART61.DAT'
	ORG $+3
_.SG66_	DB	'RESTART_H2_PROJECT.DAT'
	ORG $+2
_.SG68_	DB	'OPTIMAL_WAVEFUNCTION.DAT'
_.SG83_	DB	':BEGIN:'
	ORG $+1
_.SG87_	DB	':END:'
	ORG $+3
_.SG91_	DB	'NUMBER_OF_FUNCTIONS='
_.SG95_	DB	'FUNCTION_GENERATION='
_.SG98_	DB	'random'
	ORG $+2
_.SG101_ DB	'SYMMETRY='
	ORG $+3
_.SG104_ DB	'ON'
	ORG $+2
_.SG109_ DB	'OPTIMIZATION_METHOD='
_.SG112_ DB	'none'
_.SG115_ DB	'imsl:duming'
	ORG $+1
_.SG118_ DB	'imsl:duminf'
	ORG $+1
_.SG121_ DB	'OPTIMIZED_WAVEFUNCTION_OUTPUT_FILE='
	ORG $+1
_.SG124_ DB	'OUTPUT_FILE='
_.SG127_ DB	'RESTART_FILE='
	ORG $+3
_.SG130_ DB	'START_FROM_RESTART='
	ORG $+1
_.SG133_ DB	'true'
_.SG136_ DB	'false'
	ORG $+3
_.SG139_ DB	'RESTART_FROM='
	ORG $+3
_.SG142_ DB	'PRE-OPTIMIZE_BEFORE_STOCHASTIC_GROWTH='
	ORG $+2
_.SG145_ DB	'true'
_.SG148_ DB	'false'
	ORG $+3
_.SG151_ DB	'STOCHASTIC_END='
	ORG $+1
_.SG154_ DB	'STOCHASTIC_PERCENTAGE_GROW='
	ORG $+1
_.SG157_ DB	'true'
_.SG160_ DB	'false'
	ORG $+3
_.SG163_ DB	'STOCHASTIC_GROW_PERCENTAGE='
	ORG $+1
_.SG193_ DB	'!NOT ENOUGH MEMORY SPACE AVAILABLE ON THIS BOX! ...ABORT'
	DB	'ING!'
_.SG195_ DB	'UNKNOWN'
	ORG $+1
_.SG196_ DB	'SEQUENTIAL'
	ORG $+2
_.SG199_ DB	'Correlated Gaussian Geminal Output:  H_2'
_.SG201_ DB	'==================================================='
	ORG $+1
_.SG202_ DB	'Adamowicz Research /            [Univ. of AZ, CHEM]'
	ORG $+1
_.SG203_ DB	'Programming:[D. Gilmore, L. Adamowicz, D. Kinghorn]'
	ORG $+1
_.SG204_ DB	' Current date and time stamp...'
	ORG $+1
_.SG214_ DB	'(''     '',I2,''/'',I2.2,''/'',I4.4)'
	ORG $+2
_.SG216_ DB	'(''     '',I2,'':'',I2.2,'':'',I2.2,'':'',I2.2,'' '',A,'''
	DB	'm'')'
	ORG $+3
_.SG217_ DB	'OUTPUT STREAM BEGINS...'
	ORG $+1
_.SG218_ DB	'Correlated Gaussian Geminal Output:  H_2'
_.SG219_ DB	'==================================================='
	ORG $+1
_.SG220_ DB	'Adamowicz Research /            [Univ. of AZ, CHEM]'
	ORG $+1
_.SG221_ DB	'Programming:[D. Gilmore, L. Adamowicz, D. Kinghorn]'
	ORG $+1
_.SG222_ DB	' PROGRAM BEGINS:  TIME STAMP...'
	ORG $+1
_.SG223_ DB	'(''     '',I2,''/'',I2.2,''/'',I4.4)'
	ORG $+2
_.SG224_ DB	'(''     '',I2,'':'',I2.2,'':'',I2.2,'':'',I2.2,'' '',A,'''
	DB	'm'')'
	ORG $+3
_.SG225_ DB	'OUTPUT STREAM BEGINS...'
	ORG $+1
_.SG227_ DB	' '
	ORG $+3
_.SG231_ DB	'IMSL:DUMINF'
	ORG $+1
_.SG234_ DB	'IMSL:DUMING'
	ORG $+1
_.SG238_ DB	'READ FROM FILE'
	ORG $+2
_.SG241_ DB	'RANDOM'
	ORG $+2
_.SG244_ DB	'RESTART'
	ORG $+1
_.SG245_ DB	'NUMBER OF FUNCTIONS='
_.SG247_ DB	'FUNCTION GENERATION='
_.SG248_ DB	'OPTIMIZATION METHOD='
_.SG249_ DB	'    OUTPUT FILENAME='
_.SG250_ DB	' NEW RESTART OUTPUT='
_.SG253_ DB	'       RESTART FROM='
_.SG254_ DB	'NUMBER OF FUNCTIONS='
_.SG255_ DB	'FUNCTION GENERATION='
_.SG256_ DB	'OPTIMIZATION METHOD='
_.SG257_ DB	'    OUTPUT FILENAME='
_.SG258_ DB	' NEW RESTART OUTPUT='
_.SG261_ DB	'       RESTART FROM='
_.SG262_ DB	' '
	ORG $+3
_.SG268_ DB	'STARTING FROM DESIGNATED RESTART POINT!'
	ORG $+1
_.SG269_ DB	' Point taken from File = '
	ORG $+3
_.SG270_ DB	' Point taken from File = '
	ORG $+3
_.SG271_ DB	' '
	ORG $+3
_.SG275_ DB	':BEGIN:'
	ORG $+1
_.SG277_ DB	'NUMBER OF FUNCTIONS IN RESTART='
	ORG $+1
_.SG278_ DB	'NUMBER OF FUNCTIONS IN RESTART='
	ORG $+1
_.SG279_ DB	' '
	ORG $+3
_.SG284_ DB	'!NOT ENOUGH MEMORY SPACE AVAILABLE ON THIS BOX! ...ABORT'
	DB	'ING!'
_.SG289_ DB	'UNKNOWN'
	ORG $+1
_.SG290_ DB	'SEQUENTIAL'
	ORG $+2
_.SG299_ DB	' '
	ORG $+3
_.SG306_ DB	'ENERGY='
	ORG $+1
_.SG315_ DB	'STARTING WITH '
	ORG $+2
_.SG316_ DB	' FUNCTIONS...'
	ORG $+3
_.SG317_ DB	'STARTING ENERGY='
_.SG412_ DB	' '
	ORG $+3
_.SG413_ DB	'BEST DIAGNOSTIC VALUE IN '
	ORG $+3
_.SG414_ DB	' GUESSES='
	ORG $+3
_.SG415_ DB	'***************************'
	ORG $+1
_.SG416_ DB	'ADDING THIS FUNCTION TO Psi'
	ORG $+1
_.SG417_ DB	'***************************'
	ORG $+1
_.SG434_ DB	'NOW WE HAVE'
	ORG $+1
_.SG435_ DB	' FUNCTIONS...'
	ORG $+3
_.SG436_ DB	'EXACT ENERGY='
	ORG $+3
_.SG437_ DB	'********************************'
_.SG438_ DB	'PERFORMING GRADIENT OPTIMIZATION'
_.SG439_ DB	'********************************'
_.SG444_ DB	'PRE-OPTIMIZATION ENERGY --->'
_.SG451_ DB	'STOCHASTIC GROWTH PERCENTAGE='
	ORG $+3
_.SG452_ DB	'STOCHASTIC GROWTH PERCENTAGE='
	ORG $+3
_.SG453_ DB	' '
	ORG $+3
_.SG512_ DB	' '
	ORG $+3
_.SG513_ DB	'BEST DIAGNOSTIC VALUE IN '
	ORG $+3
_.SG514_ DB	' GUESSES='
	ORG $+3
_.SG515_ DB	'***************************'
	ORG $+1
_.SG516_ DB	'ADDING THIS FUNCTION TO Psi'
	ORG $+1
_.SG517_ DB	'***************************'
	ORG $+1
_.SG527_ DB	'NOW WE HAVE'
	ORG $+1
_.SG528_ DB	' FUNCTIONS...'
	ORG $+3
_.SG529_ DB	'EXACT ENERGY='
	ORG $+3
_.SG530_ DB	'********************************'
_.SG531_ DB	'PERFORMING GRADIENT OPTIMIZATION'
_.SG532_ DB	'********************************'
_.SG537_ DB	'PRE-OPTIMIZATION ENERGY --->'
_.SG542_ DB	'(''     '',I2,''/'',I2.2,''/'',I4.4)'
	ORG $+2
_.SG543_ DB	'(''     '',I2,'':'',I2.2,'':'',I2.2,'':'',I2.2,'' '',A,'''
	DB	'm'')'
	ORG $+3
_.SG544_ DB	'(''     '',I2,''/'',I2.2,''/'',I4.4)'
	ORG $+2
_.SG545_ DB	'(''     '',I2,'':'',I2.2,'':'',I2.2,'':'',I2.2,'' '',A,'''
	DB	'm'')'
	ORG $+3
_.SG546_ DB	' '
	ORG $+3
_.SG572_ DB	'!NOT ENOUGH MEMORY SPACE AVAILABLE ON THIS BOX! ...ABORT'
	DB	'ING!'
_.SG580_ DB	'HIT RELOAD!'
	ORG $+1
_.SG591_ DB	'LOWEST ENERGY = '
_.SG592_ DB	'UNKNOWN'
	ORG $+1
_.SG593_ DB	'SEQUENTIAL'
	ORG $+2
_.SG595_ DB	':BEGIN:'
	ORG $+1
_.SG598_ DB	':END:'
	ORG $+3
_.SG599_ DB	'ENERGY='
	ORG $+1
_.SG602_ DB	' '
	ORG $+3
_.SG603_ DB	' Square-Magnitude of the GRADIENT is'
_.SG606_ DB	' Square-Magnitude of the GRADIENT is'
_.SG607_ DB	' '
	ORG $+3
_.SG608_ DB	'     LinCoeff:.   AA1:        AA2:        B12:        RZ'
	DB	'1:       RZ2:'
	ORG $+3
_.SG609_ DB	'FNC: ==================================================='
	DB	'==================='
	ORG $+1
_.SG613_ DB	'(I3,6D12.4)'
	ORG $+1
_.SG616_ DB	'(I3,6D12.4)'
	ORG $+1
_.SG617_ DB	'--------------------------------------------------------'
	DB	'-------------------'
	ORG $+1
_.SG618_ DB	' '
	ORG $+3
_.SG620_ DB	'           OPTIMIZED ENERGY = '
	ORG $+2
_.SG621_ DB	'           OPTIMIZED ENERGY = '
	ORG $+2
_.SG622_ DB	'           BEST LIT. RESULT = '
	ORG $+2
_.SG623_ DB	'           BEST LIT. RESULT = '
	ORG $+2
_.SG624_ DB	'           ------------------------------------------'
	ORG $+3
_.SG625_ DB	'           ------------------------------------------'
	ORG $+3
_.SG626_ DB	'                 DIFFERENCE = '
	ORG $+2
_.SG627_ DB	'                 DIFFERENCE = '
	ORG $+2
_.SG628_ DB	' '
	ORG $+3
_.SG629_ DB	' PROGRAM FINISHES:  TIME STAMP...'
	ORG $+3
_.SG630_ DB	' PROGRAM FINISHES:  TIME STAMP...'
	ORG $+3
_.SG631_ DB	'(''     '',I2,''/'',I2.2,''/'',I4.4)'
	ORG $+2
_.SG632_ DB	'(''     '',I2,'':'',I2.2,'':'',I2.2,'':'',I2.2,'' '',A,'''
	DB	'm'')'
	ORG $+3
_.SG633_ DB	'(''     '',I2,''/'',I2.2,''/'',I4.4)'
	ORG $+2
_.SG634_ DB	'(''     '',I2,'':'',I2.2,'':'',I2.2,'':'',I2.2,'' '',A,'''
	DB	'm'')'
_DATA	ENDS
_CONST	SEGMENT
$T98426	DQ	04049000000000000r		; 50
$T98427	DQ	03ff0000000000000r		; 1
$T98428	DQ	00000000000000000r		; 0
$T98429	DQ	03fdf7bf4273b68c3r		; 0.491941
$T98430	DQ	03ff4058c16db8c5ar		; 1.25135
$T98431	DQ	03fe2f79b34c97dbbr		; 0.592725
$T98432	DQ	03fe7aa996ea39b5ar		; 0.739575
$T98433	DQ	03ffa97501546cfaer		; 1.66194
$T98434	DQ	03fd9b248b671f79fr		; 0.401507
$T98435	DQ	03fbc1a354f66e206r		; 0.109775
$T98436	DQ	03fa525e3b3804478r		; 0.0413047
$T98437	DQ	03fcbac3086ccbb2fr		; 0.216192
$T98438	DQ	03fe28cc9b0e7093ar		; 0.579686
$T98439	DQ	0bfe42045c9ddece2r		; -0.62894
$T98440	DQ	0bfb710bde63e3feer		; -0.0900992
$T98441	DQ	04000000000000000r		; 2
$T98442	DQ	03ff66ae7d566cf42r		; 1.4011
$T98443	DQ	04024000000000000r		; 10
$T98444	DQ	0c01c000000000000r		; -7
$T98445	DQ	03ffb333333333333r		; 1.7
$T98446	DQ	03fe6666666666666r		; 0.7
$T98447	DQ	03ff4cccccccccccdr		; 1.3
$T98448	DQ	04014000000000000r		; 5
$T98449	DQ	0c014000000000000r		; -5
$T98450	DQ	03fe66ae7d566cf42r		; 0.70055
$T98451	DQ	04031000000000000r		; 17
$T98452	DQ	03f1a36e2eb1c432dr		; 0.0001
$T98454	DQ	03fe0000000000000r		; 0.5
$T98455	DQ	0bff0000000000000r		; -1
$T98456	DQ	0c000000000000000r		; -2
$T98457	DQ	04059000000000000r		; 100
$T98458	DQ	0bff2caa73b7eed57r		; -1.17448
_CONST	ENDS
_TEXT	SEGMENT
; File C:\users\dgilmore\Graduate.School\Adamowicz.Research\H2ANL\H2.FOR
_.T518_$ = -60
$T98313 = -8
_.T521_$ = -68
_.T524_$ = -64
_K$ = -76
_.T281_$ = -40
_J$ = -16
_I$ = -80
_.T418_$ = -48
_.T167_$ = -44
_.T428_$ = -56
_.T558_$ = -28
_.T431_$ = -52
_.T566_$ = -32
_.T575_$ = -72
_main	PROC NEAR

; 1    : C     L. ADAMOWICZ RESEARCH / UNIVERSITY OF ARIZONA /  CHEM

	sub	esp, 80					; 00000050H
	push	ebx
	push	esi
	push	edi
	mov	esi, 1
	push	ebp
	mov	DWORD PTR _IMODE, esi
; File D:\MSDEV\INCLUDE/msflib.fd
; File C:\users\dgilmore\Graduate.School\Adamowicz.Research\H2ANL\H2.FOR

; 2    : C     -----------------------------------------------------
; 3    : C     Coded:									   D. Gilmore
; 4    : C
; 5    : C     --- PROGRAM CORE: ---
; 6    : C     -----------------------------------------------------   
; 7    :       IMPLICIT REAL*8 (A-H,O-Z)
; 8    : 	INTEGER ERROR,METHODGEN,METHODOPT,DFLAG,ISYMM,RESTART,WAVEOPT,
; 9    : 	-		STARTBASIS,ENDBASIS,LOOP,RELOADFLAG,INDEX,PREOPT,STGROW
; 10   :       INTEGER*2 TMPDAY,TMPMONTH,TMPYEAR,TMPHOUR,TMPMINUTE,TMPSECOND,
; 11   : 	-			TMPHUND
; 12   : 	CHARACTER	rSTR*66,FILENAME*66,METHODNAME*66,GENNAME*66,
; 13   : 	-			RESTARTNAME*66,RESTARTOUTPUT*66,OPTIMIZEDWAVE*66
; 14   : C     -----------------------------------------------------
; 15   :       REAL*8  AA1[ALLOCATABLE](:),RX1[ALLOCATABLE](:),
; 16   : 	-		RY1[ALLOCATABLE](:),RZ1[ALLOCATABLE](:),
; 17   :      -		AA2[ALLOCATABLE](:),RX2[ALLOCATABLE](:),
; 18   :      -		RY2[ALLOCATABLE](:),RZ2[ALLOCATABLE](:),
; 19   :      -		B12[ALLOCATABLE](:),EVEC[ALLOCATABLE](:)
; 20   : C
; 21   :       REAL*8	S [ALLOCATABLE](:,:),H[ALLOCATABLE](:,:),
; 22   : 	-		CO[ALLOCATABLE](:),ENRGY[ALLOCATABLE](:)
; 23   : 	REAL*8	GRAD[ALLOCATABLE](:),T1[ALLOCATABLE](:),
; 24   : 	-		T2  [ALLOCATABLE](:),T3[ALLOCATABLE](:)
; 25   : 	REAL*8  EVECSTORE[ALLOCATABLE](:,:),TSTORE1[ALLOCATABLE](:,:),
; 26   : 	-		TSTORE2[ALLOCATABLE](:,:),TSTORE3[ALLOCATABLE](:,:),
; 27   : 	-		RZ1STORE[ALLOCATABLE](:,:),RZ2STORE[ALLOCATABLE](:,:)
; 28   : 	REAL*8  EVECSTORET[ALLOCATABLE](:),TSTORE1T[ALLOCATABLE](:),
; 29   : 	-		TSTORE2T[ALLOCATABLE](:),TSTORE3T[ALLOCATABLE](:),
; 30   : 	-		RZ1STORET[ALLOCATABLE](:),RZ2STORET[ALLOCATABLE](:)
; 31   : 	REAL*4 TEMP
; 32   : 	REAL*8 Q(2,2),Y(2,2),INDX(2),F(2)
; 33   : C     -----------------------------------------------------
; 34   : 	COMMON /FINALGRAD/ GG
; 35   : C     -----------------------------------------------------
; 36   :       INCLUDE 'FLIB.FD'
; 37   : C     -----------------------------------------------------
; 38   : C     =====================================================
; 39   : C	DEFAULT VALUES:
; 40   : C     =====================================================
; 41   : 	IMODE=1
; 42   : 	METHODOPT=1

	push	14					; 0000000eH
	mov	DWORD PTR _METHODOPT, esi

; 43   : 	METHODGEN=1

	push	esi
	mov	DWORD PTR _METHODGEN, esi

; 44   :       M=2

	mov	DWORD PTR _M, 2

; 45   : 	FILENAME='H2_PROJECT.OUT'

	push	OFFSET FLAT:_.SG62_
	push	66					; 00000042H
	push	esi
	push	OFFSET FLAT:_FILENAME
	call	__FIstrcpy@24

; 46   : 	RESTARTNAME='RESTART61.DAT'

	push	13					; 0000000dH
	push	esi
	push	OFFSET FLAT:_.SG64_
	push	66					; 00000042H
	push	esi
	push	OFFSET FLAT:_RESTARTNAME
	call	__FIstrcpy@24

; 47   : 	RESTARTOUTPUT='RESTART_H2_PROJECT.DAT'

	push	22					; 00000016H
	push	esi
	push	OFFSET FLAT:_.SG66_
	push	66					; 00000042H
	push	esi
	push	OFFSET FLAT:_RESTARTOUTPUT
	call	__FIstrcpy@24

; 48   : 	OPTIMIZEDWAVE='OPTIMAL_WAVEFUNCTION.DAT'

	push	24					; 00000018H
	push	esi
	push	OFFSET FLAT:_.SG68_
	push	66					; 00000042H
	push	esi
	push	OFFSET FLAT:_OPTIMIZEDWAVE
	call	__FIstrcpy@24

; 49   : 	ISYMM=1

	xor	eax, eax
	mov	DWORD PTR _ISYMM, esi

; 50   : 	RESTART=0

	mov	DWORD PTR _RESTART, eax
	mov	DWORD PTR _STGROW, eax

; 51   : 	STGROW=0
; 52   : 	PREOPT=1

	mov	ecx, DWORD PTR _IMODE
	mov	DWORD PTR _PREOPT, esi

; 53   : 	ISTEND=256

	mov	DWORD PTR _STGROWPER, eax
	cmp	ecx, eax
	mov	DWORD PTR _ISTEND, 256			; 00000100H

; 54   : 	STGROWPER=50

	mov	DWORD PTR _STGROWPER+4, 1078525952	; 40490000H

; 55   : C     =====================================================
; 56   : 	IF (IMODE) THEN

	je	$L75
$L77:

; 57   :  100	READ(*,*) rSTR

	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	-33000					; ffff7f18H
	call	_f90_read
	add	esp, 28					; 0000001cH
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	_lread_ch
	add	esp, 8
	call	_f90_endread

; 58   : 	IF (rSTR .NE. ':BEGIN:') GOTO 100

	push	7
	push	OFFSET FLAT:_.SG83_
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	__FIclne@16
	test	al, al
	jne	SHORT $L77
	mov	esi, 1
	xor	ebx, ebx
	mov	edi, 1072693248				; 3ff00000H
	mov	ebp, 2
$L85:

; 59   :  200  READ(*,*) rSTR

	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	-33000					; ffff7f18H
	call	_f90_read
	add	esp, 28					; 0000001cH
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	_lread_ch
	add	esp, 8
	call	_f90_endread

; 60   : 	IF (rSTR .EQ. ':END:') GOTO 300

	push	5
	push	OFFSET FLAT:_.SG87_
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	__FIcleq@16
	test	al, al
	jne	$L75

; 61   : 		IF (rSTR .EQ. 'NUMBER_OF_FUNCTIONS=') READ(*,*) M

	push	20					; 00000014H
	push	OFFSET FLAT:_.SG91_
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	__FIcleq@16
	test	al, al
	je	SHORT $L90
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	-33000					; ffff7f18H
	call	_f90_read
	add	esp, 28					; 0000001cH
	push	OFFSET FLAT:_M
	call	_lread_i
	add	esp, 4
	call	_f90_endread
$L90:

; 62   : C
; 63   : 		IF (rSTR .EQ. 'FUNCTION_GENERATION=') THEN

	push	20					; 00000014H
	push	OFFSET FLAT:_.SG95_
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	__FIcleq@16
	test	al, al
	je	SHORT $L94

; 64   : 			READ(*,*) rSTR

	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	-33000					; ffff7f18H
	call	_f90_read
	add	esp, 28					; 0000001cH
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	_lread_ch
	add	esp, 8
	call	_f90_endread

; 65   : 			IF (rSTR .EQ. 'random') METHODGEN=1

	push	6
	push	OFFSET FLAT:_.SG98_
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	__FIcleq@16
	test	al, al
	je	SHORT $L94
	mov	DWORD PTR _METHODGEN, esi
$L94:

; 66   : 		END IF
; 67   : C
; 68   : 		IF (rSTR .EQ. 'SYMMETRY=') THEN

	push	9
	push	OFFSET FLAT:_.SG101_
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	__FIcleq@16
	test	al, al
	je	SHORT $L100

; 69   : 			READ(*,*) rSTR

	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	-33000					; ffff7f18H
	call	_f90_read
	add	esp, 28					; 0000001cH
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	_lread_ch
	add	esp, 8
	call	_f90_endread

; 70   : 			IF (rSTR .EQ. 'ON') THEN ISYMM=1

	push	ebp
	push	OFFSET FLAT:_.SG104_
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	__FIcleq@16
	test	al, al
	je	SHORT $L103
	mov	DWORD PTR _THENISYMM, ebx
	mov	DWORD PTR _THENISYMM+4, edi
$L103:

; 71   : 				ELSE ISYMM=0

	mov	DWORD PTR _ELSEISYMM, ebx
	mov	DWORD PTR _ELSEISYMM+4, ebx
$L100:

; 72   : 		END IF
; 73   : C
; 74   : 		IF (rSTR .EQ. 'OPTIMIZATION_METHOD=') THEN

	push	20					; 00000014H
	push	OFFSET FLAT:_.SG109_
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	__FIcleq@16
	test	al, al
	je	SHORT $L108

; 75   : 			READ(*,*) rSTR

	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	-33000					; ffff7f18H
	call	_f90_read
	add	esp, 28					; 0000001cH
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	_lread_ch
	add	esp, 8
	call	_f90_endread

; 76   : 			IF (rSTR .EQ. 'none') METHODOPT=0

	push	4
	push	OFFSET FLAT:_.SG112_
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	__FIcleq@16
	test	al, al
	je	SHORT $L111
	mov	DWORD PTR _METHODOPT, ebx
$L111:

; 77   : 			IF (rSTR .EQ. 'imsl:duming') METHODOPT=1

	push	11					; 0000000bH
	push	OFFSET FLAT:_.SG115_
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	__FIcleq@16
	test	al, al
	je	SHORT $L114
	mov	DWORD PTR _METHODOPT, esi
$L114:

; 78   : 			IF (rSTR .EQ. 'imsl:duminf') METHODOPT=2

	push	11					; 0000000bH
	push	OFFSET FLAT:_.SG118_
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	__FIcleq@16
	test	al, al
	je	SHORT $L108
	mov	DWORD PTR _METHODOPT, ebp
$L108:

; 79   : 		END IF
; 80   : C
; 81   : 		IF (rSTR .EQ. 'OPTIMIZED_WAVEFUNCTION_OUTPUT_FILE=') READ

	push	35					; 00000023H
	push	OFFSET FLAT:_.SG121_
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	__FIcleq@16
	test	al, al
	je	SHORT $L120
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	-33000					; ffff7f18H
	call	_f90_read
	add	esp, 28					; 0000001cH
	push	66					; 00000042H
	push	OFFSET FLAT:_OPTIMIZEDWAVE
	call	_lread_ch
	add	esp, 8
	call	_f90_endread
$L120:

; 82   : 	-			(*,*) OPTIMIZEDWAVE
; 83   : C
; 84   : 		IF (rSTR .EQ. 'OUTPUT_FILE=') READ(*,*) FILENAME

	push	12					; 0000000cH
	push	OFFSET FLAT:_.SG124_
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	__FIcleq@16
	test	al, al
	je	SHORT $L123
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	-33000					; ffff7f18H
	call	_f90_read
	add	esp, 28					; 0000001cH
	push	66					; 00000042H
	push	OFFSET FLAT:_FILENAME
	call	_lread_ch
	add	esp, 8
	call	_f90_endread
$L123:

; 85   : C
; 86   : 		IF (rSTR .EQ. 'RESTART_FILE=') READ(*,*) RESTARTOUTPUT

	push	13					; 0000000dH
	push	OFFSET FLAT:_.SG127_
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	__FIcleq@16
	test	al, al
	je	SHORT $L126
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	-33000					; ffff7f18H
	call	_f90_read
	add	esp, 28					; 0000001cH
	push	66					; 00000042H
	push	OFFSET FLAT:_RESTARTOUTPUT
	call	_lread_ch
	add	esp, 8
	call	_f90_endread
$L126:

; 87   : C
; 88   : 		IF (rSTR .EQ. 'START_FROM_RESTART=') THEN

	push	19					; 00000013H
	push	OFFSET FLAT:_.SG130_
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	__FIcleq@16
	test	al, al
	je	SHORT $L129

; 89   : 			READ(*,*) rSTR

	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	-33000					; ffff7f18H
	call	_f90_read
	add	esp, 28					; 0000001cH
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	_lread_ch
	add	esp, 8
	call	_f90_endread

; 90   : 			IF (rSTR .EQ. 'true') RESTART=1

	push	4
	push	OFFSET FLAT:_.SG133_
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	__FIcleq@16
	test	al, al
	je	SHORT $L132
	mov	DWORD PTR _RESTART, esi
$L132:

; 91   : 			IF (rSTR .EQ. 'false') RESTART=0

	push	5
	push	OFFSET FLAT:_.SG136_
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	__FIcleq@16
	test	al, al
	je	SHORT $L129
	mov	DWORD PTR _RESTART, ebx
$L129:

; 92   : 		END IF
; 93   : C
; 94   : 		IF (rSTR .EQ. 'RESTART_FROM=') READ(*,*) RESTARTNAME

	push	13					; 0000000dH
	push	OFFSET FLAT:_.SG139_
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	__FIcleq@16
	test	al, al
	je	SHORT $L138
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	-33000					; ffff7f18H
	call	_f90_read
	add	esp, 28					; 0000001cH
	push	66					; 00000042H
	push	OFFSET FLAT:_RESTARTNAME
	call	_lread_ch
	add	esp, 8
	call	_f90_endread
$L138:

; 95   : C
; 96   : 		IF (rSTR .EQ. 'PRE-OPTIMIZE_BEFORE_STOCHASTIC_GROWTH=') THEN

	push	38					; 00000026H
	push	OFFSET FLAT:_.SG142_
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	__FIcleq@16
	test	al, al
	je	SHORT $L141

; 97   : 			READ(*,*) rSTR

	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	-33000					; ffff7f18H
	call	_f90_read
	add	esp, 28					; 0000001cH
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	_lread_ch
	add	esp, 8
	call	_f90_endread

; 98   : 			IF (rSTR .EQ. 'true') PREOPT=1

	push	4
	push	OFFSET FLAT:_.SG145_
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	__FIcleq@16
	test	al, al
	je	SHORT $L144
	mov	DWORD PTR _PREOPT, esi
$L144:

; 99   : 			IF (rSTR .EQ. 'false') PREOPT=0

	push	5
	push	OFFSET FLAT:_.SG148_
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	__FIcleq@16
	test	al, al
	je	SHORT $L141
	mov	DWORD PTR _PREOPT, ebx
$L141:

; 100  : 		END IF
; 101  : C
; 102  : 		IF (rSTR .EQ. 'STOCHASTIC_END=') READ(*,*) ISTEND

	push	15					; 0000000fH
	push	OFFSET FLAT:_.SG151_
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	__FIcleq@16
	test	al, al
	je	SHORT $L150
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	-33000					; ffff7f18H
	call	_f90_read
	add	esp, 28					; 0000001cH
	push	OFFSET FLAT:_ISTEND
	call	_lread_i
	add	esp, 4
	call	_f90_endread
$L150:

; 103  : C
; 104  : 		IF (rSTR .EQ. 'STOCHASTIC_PERCENTAGE_GROW=') THEN

	push	27					; 0000001bH
	push	OFFSET FLAT:_.SG154_
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	__FIcleq@16
	test	al, al
	je	SHORT $L153

; 105  : 			READ (*,*) rSTR

	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	-33000					; ffff7f18H
	call	_f90_read
	add	esp, 28					; 0000001cH
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	_lread_ch
	add	esp, 8
	call	_f90_endread

; 106  : 			IF (rSTR .EQ. 'true') STGROW=1

	push	4
	push	OFFSET FLAT:_.SG157_
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	__FIcleq@16
	test	al, al
	je	SHORT $L156
	mov	DWORD PTR _STGROW, esi
$L156:

; 107  : 			IF (rSTR .EQ. 'false') STGROW=0

	push	5
	push	OFFSET FLAT:_.SG160_
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	__FIcleq@16
	test	al, al
	je	SHORT $L153
	mov	DWORD PTR _STGROW, ebx
$L153:

; 108  : 		END IF
; 109  : C
; 110  : 		IF (rSTR .EQ. 'STOCHASTIC_GROW_PERCENTAGE=') 

	push	27					; 0000001bH
	push	OFFSET FLAT:_.SG163_
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	__FIcleq@16
	test	al, al
	je	$L85

; 111  : 	-			READ(*,*) STGROWPER

	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	-33000					; ffff7f18H
	call	_f90_read
	add	esp, 28					; 0000001cH
	push	OFFSET FLAT:_STGROWPER
	call	_lread_d
	add	esp, 4
	call	_f90_endread

; 112  : C
; 113  : 	GOTO 200

	jmp	$L85
$L75:

; 114  : C     =====================================================
; 115  :  300	END IF
; 116  :  	MMAX=M

	mov	DWORD PTR _.T167_$[esp+96], 0
	mov	eax, DWORD PTR _M
	push	eax
	mov	DWORD PTR _MMAX, eax

; 117  : C
; 118  :       ALLOCATE(	AA1(MMAX),RX1(MMAX),RY1(MMAX),RZ1(MMAX),

	push	1
	push	1
	push	OFFSET FLAT:_7RX1
	call	_make_info
	lea	eax, DWORD PTR _.T167_$[esp+112]
	add	esp, 16					; 00000010H
	push	eax
	push	8
	push	OFFSET FLAT:_7RX1
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _RX1, eax
	mov	eax, DWORD PTR _MMAX
	push	eax
	push	1
	push	1
	push	OFFSET FLAT:_7RY1
	call	_make_info
	lea	eax, DWORD PTR _.T167_$[esp+112]
	add	esp, 16					; 00000010H
	push	eax
	push	8
	push	OFFSET FLAT:_7RY1
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _RY1, eax
	mov	eax, DWORD PTR _MMAX
	push	eax
	push	1
	push	1
	push	OFFSET FLAT:_7RZ1
	call	_make_info
	lea	eax, DWORD PTR _.T167_$[esp+112]
	add	esp, 16					; 00000010H
	push	eax
	push	8
	push	OFFSET FLAT:_7RZ1
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _RZ1, eax
	mov	eax, DWORD PTR _MMAX
	push	eax
	push	1
	push	1
	push	OFFSET FLAT:_7AA2
	call	_make_info
	lea	eax, DWORD PTR _.T167_$[esp+112]
	add	esp, 16					; 00000010H
	push	eax
	push	8
	push	OFFSET FLAT:_7AA2
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _AA2, eax
	mov	eax, DWORD PTR _MMAX
	push	eax
	push	1
	push	1
	push	OFFSET FLAT:_7RX2
	call	_make_info
	lea	eax, DWORD PTR _.T167_$[esp+112]
	add	esp, 16					; 00000010H
	push	eax
	push	8
	push	OFFSET FLAT:_7RX2
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _RX2, eax
	mov	eax, DWORD PTR _MMAX
	push	eax
	push	1
	push	1
	push	OFFSET FLAT:_7RY2
	call	_make_info
	lea	eax, DWORD PTR _.T167_$[esp+112]
	add	esp, 16					; 00000010H
	push	eax
	push	8
	push	OFFSET FLAT:_7RY2
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _RY2, eax
	mov	eax, DWORD PTR _MMAX
	push	eax
	push	1
	push	1
	push	OFFSET FLAT:_7RZ2
	call	_make_info
	lea	eax, DWORD PTR _.T167_$[esp+112]
	add	esp, 16					; 00000010H
	push	eax
	push	8
	push	OFFSET FLAT:_7RZ2
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _RZ2, eax
	mov	eax, DWORD PTR _MMAX
	push	eax
	push	1
	push	1
	push	OFFSET FLAT:_7B12
	call	_make_info
	lea	eax, DWORD PTR _.T167_$[esp+112]
	add	esp, 16					; 00000010H
	push	eax
	push	8
	push	OFFSET FLAT:_7B12
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _B12, eax
	mov	eax, DWORD PTR _MMAX
	push	eax
	push	1
	push	1
	push	OFFSET FLAT:_7EVEC
	call	_make_info
	lea	eax, DWORD PTR _.T167_$[esp+112]
	add	esp, 16					; 00000010H
	push	eax
	push	8
	push	OFFSET FLAT:_7EVEC
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _EVEC, eax
	mov	eax, DWORD PTR _MMAX
	push	eax
	push	1
	push	eax
	push	1
	push	2
	push	OFFSET FLAT:_7S
	call	_make_info
	lea	eax, DWORD PTR _.T167_$[esp+120]
	add	esp, 24					; 00000018H
	push	eax
	push	8
	push	OFFSET FLAT:_7S
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _S, eax
	mov	eax, DWORD PTR _MMAX
	push	eax
	push	1
	push	eax
	push	1
	push	2
	push	OFFSET FLAT:_7H
	call	_make_info
	lea	eax, DWORD PTR _.T167_$[esp+120]
	add	esp, 24					; 00000018H
	push	eax
	push	8
	push	OFFSET FLAT:_7H
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _H, eax
	mov	eax, DWORD PTR _MMAX
	push	eax
	push	1
	push	1
	push	OFFSET FLAT:_7CO
	call	_make_info
	lea	eax, DWORD PTR _.T167_$[esp+112]
	add	esp, 16					; 00000010H
	push	eax
	push	8
	push	OFFSET FLAT:_7CO
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _CO, eax
	mov	eax, DWORD PTR _MMAX
	lea	eax, DWORD PTR [eax+eax*2]
	add	eax, eax
	push	eax
	push	1
	push	1
	push	OFFSET FLAT:_7GRAD
	call	_make_info
	lea	ecx, DWORD PTR _.T167_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7GRAD
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _MMAX
	mov	DWORD PTR _GRAD, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7T1
	call	_make_info
	lea	ecx, DWORD PTR _.T167_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7T1
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _MMAX
	mov	DWORD PTR _T1, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7T2
	call	_make_info
	lea	ecx, DWORD PTR _.T167_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7T2
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _MMAX
	mov	DWORD PTR _T2, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7T3
	call	_make_info
	lea	ecx, DWORD PTR _.T167_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7T3
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _LOOP
	mov	DWORD PTR _T3, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7ENRGY
	call	_make_info
	lea	ecx, DWORD PTR _.T167_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7ENRGY
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _MMAX
	mov	DWORD PTR _ENRGY, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7AA1
	call	_make_info
	lea	ecx, DWORD PTR _.T167_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7AA1
	call	__Allocate_i
	mov	ecx, DWORD PTR _.T167_$[esp+108]
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _AA1, eax
	mov	DWORD PTR _ERROR, ecx

; 119  : 	-			AA2(MMAX),RX2(MMAX),RY2(MMAX),RZ2(MMAX),
; 120  : 	-			B12(MMAX),EVEC(MMAX),S(MMAX,MMAX),H(MMAX,MMAX),
; 121  :      -			CO(MMAX),GRAD(MMAX*6),T1(MMAX),T2(MMAX),
; 122  :      -			T3(MMAX),ENRGY(LOOP),STAT=ERROR)
; 123  : C
; 124  : 	IF (ERROR .NE. 0) STOP '!NOT ENOUGH MEMORY SPACE AVAILABLE ON THIS

	test	ecx, ecx
	je	SHORT $L190
	push	OFFSET FLAT:_.SG193_
	push	60					; 0000003cH
	call	__FFsSTOP
	add	esp, 8
$L190:

; 125  : 	- BOX! ...ABORTING!'
; 126  : C
; 127  : 	OPEN (7, FILE = FILENAME, ACCESS = 'SEQUENTIAL',

	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	10					; 0000000aH
	push	OFFSET FLAT:_.SG196_
	push	7
	push	OFFSET FLAT:_.SG195_
	push	66					; 00000042H
	push	OFFSET FLAT:_FILENAME
	push	7
	call	_f90_open
	add	esp, 124				; 0000007cH

; 128  :      -		STATUS = 'UNKNOWN')
; 129  : C     -----------------------------------------------------
; 130  :       WRITE(6,*)'Correlated Gaussian Geminal Output:  H_2'

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	40					; 00000028H
	push	OFFSET FLAT:_.SG199_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 131  :       WRITE(6,*)'==================================================='

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	51					; 00000033H
	push	OFFSET FLAT:_.SG201_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 132  :       WRITE(6,*)'Adamowicz Research /            [Univ. of AZ, CHEM]'

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	51					; 00000033H
	push	OFFSET FLAT:_.SG202_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 133  :       WRITE(6,*)'Programming:[D. Gilmore, L. Adamowicz, D. Kinghorn]' 

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	51					; 00000033H
	push	OFFSET FLAT:_.SG203_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 134  : C      
; 135  :       WRITE (6,*) ' Current date and time stamp...'

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	31					; 0000001fH
	push	OFFSET FLAT:_.SG204_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 136  :       CALL GETDAT(tmpyear, tmpmonth, tmpday)

	push	OFFSET FLAT:_TMPDAY
	push	OFFSET FLAT:_TMPMONTH
	push	OFFSET FLAT:_TMPYEAR
	call	_GETDAT@12

; 137  :       CALL GETTIM(tmphour, tmpminute, tmpsecond, tmphund)  

	push	OFFSET FLAT:_TMPHUND
	push	OFFSET FLAT:_TMPSECOND
	push	OFFSET FLAT:_TMPMINUTE
	push	OFFSET FLAT:_TMPHOUR
	call	_GETTIM@16

; 138  :       WRITE (6, 9001) tmpmonth, tmpday, tmpyear

	push	0
	push	0
	push	0
	push	-1
	push	OFFSET FLAT:_.SG214_
	push	6
	call	_f90_write
	mov	ax, WORD PTR _TMPMONTH
	add	esp, 24					; 00000018H
	push	eax
	call	_write_i2
	mov	ax, WORD PTR _TMPDAY
	add	esp, 4
	push	eax
	call	_write_i2
	mov	ax, WORD PTR _TMPYEAR
	add	esp, 4
	push	eax
	call	_write_i2
	add	esp, 4
	call	_f90_endwrite

; 139  :       WRITE (6, 9002) tmphour, tmpminute, tmpsecond, tmphund

	push	0
	push	0
	push	0
	push	-1
	push	OFFSET FLAT:_.SG216_
	push	6
	call	_f90_write
	mov	ax, WORD PTR _TMPHOUR
	add	esp, 24					; 00000018H
	push	eax
	call	_write_i2
	mov	ax, WORD PTR _TMPMINUTE
	add	esp, 4
	push	eax
	call	_write_i2
	mov	ax, WORD PTR _TMPSECOND
	add	esp, 4
	push	eax
	call	_write_i2
	mov	ax, WORD PTR _TMPHUND
	add	esp, 4
	push	eax
	call	_write_i2
	add	esp, 4
	call	_f90_endwrite

; 140  :       WRITE(6,*)'OUTPUT STREAM BEGINS...'

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	23					; 00000017H
	push	OFFSET FLAT:_.SG217_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 141  :       WRITE(7,*)'Correlated Gaussian Geminal Output:  H_2'

	push	0
	push	0
	push	0
	push	0
	push	0
	push	7
	call	_f90_write
	add	esp, 24					; 00000018H
	push	40					; 00000028H
	push	OFFSET FLAT:_.SG218_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 142  :       WRITE(7,*)'==================================================='

	push	0
	push	0
	push	0
	push	0
	push	0
	push	7
	call	_f90_write
	add	esp, 24					; 00000018H
	push	51					; 00000033H
	push	OFFSET FLAT:_.SG219_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 143  :       WRITE(7,*)'Adamowicz Research /            [Univ. of AZ, CHEM]'

	push	0
	push	0
	push	0
	push	0
	push	0
	push	7
	call	_f90_write
	add	esp, 24					; 00000018H
	push	51					; 00000033H
	push	OFFSET FLAT:_.SG220_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 144  :       WRITE(7,*)'Programming:[D. Gilmore, L. Adamowicz, D. Kinghorn]' 

	push	0
	push	0
	push	0
	push	0
	push	0
	push	7
	call	_f90_write
	add	esp, 24					; 00000018H
	push	51					; 00000033H
	push	OFFSET FLAT:_.SG221_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 145  : C      
; 146  :       WRITE (7,*) ' PROGRAM BEGINS:  TIME STAMP...'

	push	0
	push	0
	push	0
	push	0
	push	0
	push	7
	call	_f90_write
	add	esp, 24					; 00000018H
	push	31					; 0000001fH
	push	OFFSET FLAT:_.SG222_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 147  :       CALL GETDAT(tmpyear, tmpmonth, tmpday)

	push	OFFSET FLAT:_TMPDAY
	push	OFFSET FLAT:_TMPMONTH
	push	OFFSET FLAT:_TMPYEAR
	call	_GETDAT@12

; 148  :       CALL GETTIM(tmphour, tmpminute, tmpsecond, tmphund)  

	push	OFFSET FLAT:_TMPHUND
	push	OFFSET FLAT:_TMPSECOND
	push	OFFSET FLAT:_TMPMINUTE
	push	OFFSET FLAT:_TMPHOUR
	call	_GETTIM@16

; 149  :       WRITE (7, 9003) tmpmonth, tmpday, tmpyear

	push	0
	push	0
	push	0
	push	-1
	push	OFFSET FLAT:_.SG223_
	push	7
	call	_f90_write
	mov	ax, WORD PTR _TMPMONTH
	add	esp, 24					; 00000018H
	push	eax
	call	_write_i2
	mov	ax, WORD PTR _TMPDAY
	add	esp, 4
	push	eax
	call	_write_i2
	mov	ax, WORD PTR _TMPYEAR
	add	esp, 4
	push	eax
	call	_write_i2
	add	esp, 4
	call	_f90_endwrite

; 150  :       WRITE (7, 9004) tmphour, tmpminute, tmpsecond, tmphund

	push	0
	push	0
	push	0
	push	-1
	push	OFFSET FLAT:_.SG224_
	push	7
	call	_f90_write
	mov	ax, WORD PTR _TMPHOUR
	add	esp, 24					; 00000018H
	push	eax
	call	_write_i2
	mov	ax, WORD PTR _TMPMINUTE
	add	esp, 4
	push	eax
	call	_write_i2
	mov	ax, WORD PTR _TMPSECOND
	add	esp, 4
	push	eax
	call	_write_i2
	mov	ax, WORD PTR _TMPHUND
	add	esp, 4
	push	eax
	call	_write_i2
	add	esp, 4
	call	_f90_endwrite

; 151  :       WRITE(7,*)'OUTPUT STREAM BEGINS...'

	push	0
	push	0
	push	0
	push	0
	push	0
	push	7
	call	_f90_write
	add	esp, 24					; 00000018H
	push	23					; 00000017H
	push	OFFSET FLAT:_.SG225_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 152  : 	CALL BIOUT(' ')

	push	1
	push	OFFSET FLAT:_.SG227_
	call	_BIOUT@8

; 153  : C     -----------------------------------------------------
; 154  : 	IF(METHODOPT.EQ.0) METHODNAME='IMSL:DUMINF'

	mov	eax, DWORD PTR _METHODOPT
	test	eax, eax
	jne	SHORT $L228
	push	11					; 0000000bH
	push	1
	push	OFFSET FLAT:_.SG231_
	push	66					; 00000042H
	push	1
	push	OFFSET FLAT:_METHODNAME
	call	__FIstrcpy@24
$L228:

; 155  : 	IF(METHODOPT.EQ.1) METHODNAME='IMSL:DUMING'

	mov	eax, DWORD PTR _METHODOPT
	cmp	eax, 1
	jne	SHORT $L232
	push	11					; 0000000bH
	push	1
	push	OFFSET FLAT:_.SG234_
	push	66					; 00000042H
	push	1
	push	OFFSET FLAT:_METHODNAME
	call	__FIstrcpy@24
$L232:

; 156  : 	IF(METHODGEN.EQ.0) GENNAME='READ FROM FILE'

	mov	eax, DWORD PTR _METHODGEN
	test	eax, eax
	jne	SHORT $L235
	push	14					; 0000000eH
	push	1
	push	OFFSET FLAT:_.SG238_
	push	66					; 00000042H
	push	1
	push	OFFSET FLAT:_GENNAME
	call	__FIstrcpy@24
$L235:

; 157  : 	IF(METHODGEN.EQ.1) GENNAME='RANDOM'

	mov	eax, DWORD PTR _METHODGEN
	cmp	eax, 1
	jne	SHORT $L239
	push	6
	push	1
	push	OFFSET FLAT:_.SG241_
	push	66					; 00000042H
	push	1
	push	OFFSET FLAT:_GENNAME
	call	__FIstrcpy@24
$L239:

; 158  : 	IF(RESTART) THEN

	mov	eax, DWORD PTR _RESTART
	test	eax, eax
	je	SHORT $L242

; 159  : 		M=0

	mov	DWORD PTR _M, 0

; 160  : 		GENNAME="RESTART"

	push	7
	push	1
	push	OFFSET FLAT:_.SG244_
	push	66					; 00000042H
	push	1
	push	OFFSET FLAT:_GENNAME
	call	__FIstrcpy@24
$L242:

; 161  : 	END IF
; 162  : 		WRITE(6,*)'NUMBER OF FUNCTIONS=',M

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	20					; 00000014H
	push	OFFSET FLAT:_.SG245_
	call	_lwrite_ch
	add	esp, 8
	mov	eax, DWORD PTR _M
	push	eax
	call	_lwrite_i
	add	esp, 4
	call	_f90_endwrite

; 163  : 		WRITE(6,*)'FUNCTION GENERATION=',GENNAME

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	20					; 00000014H
	push	OFFSET FLAT:_.SG247_
	call	_lwrite_ch
	add	esp, 8
	push	66					; 00000042H
	push	OFFSET FLAT:_GENNAME
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 164  : 		WRITE(6,*)'OPTIMIZATION METHOD=',METHODNAME

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	20					; 00000014H
	push	OFFSET FLAT:_.SG248_
	call	_lwrite_ch
	add	esp, 8
	push	66					; 00000042H
	push	OFFSET FLAT:_METHODNAME
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 165  : 		WRITE(6,*)'    OUTPUT FILENAME=',FILENAME

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	20					; 00000014H
	push	OFFSET FLAT:_.SG249_
	call	_lwrite_ch
	add	esp, 8
	push	66					; 00000042H
	push	OFFSET FLAT:_FILENAME
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 166  : 	    WRITE(6,*)' NEW RESTART OUTPUT=',RESTARTOUTPUT

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	20					; 00000014H
	push	OFFSET FLAT:_.SG250_
	call	_lwrite_ch
	add	esp, 8
	push	66					; 00000042H
	push	OFFSET FLAT:_RESTARTOUTPUT
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 167  : 	IF(RESTART) WRITE(6,*)'       RESTART FROM=',RESTARTNAME

	mov	eax, DWORD PTR _RESTART
	test	eax, eax
	je	SHORT $L251
	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	20					; 00000014H
	push	OFFSET FLAT:_.SG253_
	call	_lwrite_ch
	add	esp, 8
	push	66					; 00000042H
	push	OFFSET FLAT:_RESTARTNAME
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite
$L251:

; 168  : C
; 169  : 		WRITE(7,*)'NUMBER OF FUNCTIONS=',M

	push	0
	push	0
	push	0
	push	0
	push	0
	push	7
	call	_f90_write
	add	esp, 24					; 00000018H
	push	20					; 00000014H
	push	OFFSET FLAT:_.SG254_
	call	_lwrite_ch
	add	esp, 8
	mov	eax, DWORD PTR _M
	push	eax
	call	_lwrite_i
	add	esp, 4
	call	_f90_endwrite

; 170  : 		WRITE(7,*)'FUNCTION GENERATION=',GENNAME

	push	0
	push	0
	push	0
	push	0
	push	0
	push	7
	call	_f90_write
	add	esp, 24					; 00000018H
	push	20					; 00000014H
	push	OFFSET FLAT:_.SG255_
	call	_lwrite_ch
	add	esp, 8
	push	66					; 00000042H
	push	OFFSET FLAT:_GENNAME
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 171  : 		WRITE(7,*)'OPTIMIZATION METHOD=',METHODNAME

	push	0
	push	0
	push	0
	push	0
	push	0
	push	7
	call	_f90_write
	add	esp, 24					; 00000018H
	push	20					; 00000014H
	push	OFFSET FLAT:_.SG256_
	call	_lwrite_ch
	add	esp, 8
	push	66					; 00000042H
	push	OFFSET FLAT:_METHODNAME
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 172  : 		WRITE(7,*)'    OUTPUT FILENAME=',FILENAME

	push	0
	push	0
	push	0
	push	0
	push	0
	push	7
	call	_f90_write
	add	esp, 24					; 00000018H
	push	20					; 00000014H
	push	OFFSET FLAT:_.SG257_
	call	_lwrite_ch
	add	esp, 8
	push	66					; 00000042H
	push	OFFSET FLAT:_FILENAME
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 173  : 	    WRITE(7,*)' NEW RESTART OUTPUT=',RESTARTOUTPUT

	push	0
	push	0
	push	0
	push	0
	push	0
	push	7
	call	_f90_write
	add	esp, 24					; 00000018H
	push	20					; 00000014H
	push	OFFSET FLAT:_.SG258_
	call	_lwrite_ch
	add	esp, 8
	push	66					; 00000042H
	push	OFFSET FLAT:_RESTARTOUTPUT
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 174  : 	IF(RESTART) WRITE(7,*)'       RESTART FROM=',RESTARTNAME

	mov	eax, DWORD PTR _RESTART
	test	eax, eax
	je	SHORT $L259
	push	0
	push	0
	push	0
	push	0
	push	0
	push	7
	call	_f90_write
	add	esp, 24					; 00000018H
	push	20					; 00000014H
	push	OFFSET FLAT:_.SG261_
	call	_lwrite_ch
	add	esp, 8
	push	66					; 00000042H
	push	OFFSET FLAT:_RESTARTNAME
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite
$L259:

; 175  : 	CALL BIOUT(' ')

	push	1
	push	OFFSET FLAT:_.SG262_
	call	_BIOUT@8

; 176  : C     -----------------------------------------------------
; 177  : C     LOAD DATA:
; 178  : C     -----------------------------------------------------
; 179  : 	IF (RESTART.EQ.1) GOTO 212

	mov	eax, DWORD PTR _RESTART
	cmp	eax, 1
	jne	$L98342

; 181  :       GOTO 313
; 182  :  212	CALL BIOUT('STARTING FROM DESIGNATED RESTART POINT!')

	push	39					; 00000027H
	push	OFFSET FLAT:_.SG268_
	call	_BIOUT@8

; 183  : 	WRITE(6,*)' Point taken from File = ',RESTARTNAME

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	25					; 00000019H
	push	OFFSET FLAT:_.SG269_
	call	_lwrite_ch
	add	esp, 8
	push	66					; 00000042H
	push	OFFSET FLAT:_RESTARTNAME
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 184  : 	WRITE(7,*)' Point taken from File = ',RESTARTNAME

	push	0
	push	0
	push	0
	push	0
	push	0
	push	7
	call	_f90_write
	add	esp, 24					; 00000018H
	push	25					; 00000019H
	push	OFFSET FLAT:_.SG270_
	call	_lwrite_ch
	add	esp, 8
	push	66					; 00000042H
	push	OFFSET FLAT:_RESTARTNAME
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 185  : 	CALL BIOUT(' ')

	push	1
	push	OFFSET FLAT:_.SG271_
	call	_BIOUT@8

; 186  : 	CLOSE (8)

	push	0
	push	0
	push	0
	push	8
	call	_f90_close
	add	esp, 16					; 00000010H

; 187  : 	OPEN (5,FILE=RESTARTNAME)

	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	66					; 00000042H
	push	OFFSET FLAT:_RESTARTNAME
	push	5
	call	_f90_open
	add	esp, 124				; 0000007cH
$L273:

; 188  :  414		READ(5,*)rSTR

	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	5
	call	_f90_read
	add	esp, 28					; 0000001cH
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	_lread_ch
	add	esp, 8
	call	_f90_endread

; 189  : 		IF (rSTR .NE. ':BEGIN:') GOTO 414

	push	7
	push	OFFSET FLAT:_.SG275_
	push	66					; 00000042H
	push	OFFSET FLAT:_RSTR
	call	__FIclne@16
	test	al, al
	jne	SHORT $L273

; 190  : 	READ(5,*)M

	xor	esi, esi
	push	esi
	push	esi
	push	esi
	push	esi
	push	esi
	push	esi
	push	5
	call	_f90_read
	add	esp, 28					; 0000001cH
	push	OFFSET FLAT:_M
	call	_lread_i
	add	esp, 4
	call	_f90_endread

; 191  : 		WRITE(6,*)'NUMBER OF FUNCTIONS IN RESTART=',M

	push	esi
	push	esi
	push	esi
	push	esi
	push	esi
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	31					; 0000001fH
	push	OFFSET FLAT:_.SG277_
	call	_lwrite_ch
	add	esp, 8
	mov	eax, DWORD PTR _M
	push	eax
	call	_lwrite_i
	add	esp, 4
	call	_f90_endwrite

; 192  : 		WRITE(7,*)'NUMBER OF FUNCTIONS IN RESTART=',M

	push	esi
	push	esi
	push	esi
	push	esi
	push	esi
	push	7
	call	_f90_write
	add	esp, 24					; 00000018H
	push	31					; 0000001fH
	push	OFFSET FLAT:_.SG278_
	call	_lwrite_ch
	add	esp, 8
	mov	eax, DWORD PTR _M
	push	eax
	call	_lwrite_i
	add	esp, 4
	call	_f90_endwrite

; 193  : 	CALL BIOUT(' ')

	push	1
	push	OFFSET FLAT:_.SG279_
	call	_BIOUT@8

; 194  :   	MMAX=M

	mov	eax, DWORD PTR _M
	push	esi
	push	8
	mov	ecx, DWORD PTR _RX1
	push	OFFSET FLAT:_7RX1
	mov	DWORD PTR _MMAX, eax

; 195  : 	DEALLOCATE (AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,EVEC,S,H,CO,

	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _RY1
	mov	DWORD PTR _RX1, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7RY1
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _RZ1
	mov	DWORD PTR _RY1, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7RZ1
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _AA2
	mov	DWORD PTR _RZ1, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7AA2
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _RX2
	mov	DWORD PTR _AA2, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7RX2
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _RY2
	mov	DWORD PTR _RX2, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7RY2
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	DWORD PTR _RY2, esi
	push	esi
	mov	ecx, DWORD PTR _RZ2
	push	8
	push	OFFSET FLAT:_7RZ2
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _B12
	mov	DWORD PTR _RZ2, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7B12
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _EVEC
	mov	DWORD PTR _B12, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7EVEC
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _S
	mov	DWORD PTR _EVEC, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7S
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _H
	mov	DWORD PTR _S, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7H
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _CO
	mov	DWORD PTR _H, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7CO
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _GRAD
	mov	DWORD PTR _CO, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7GRAD
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _T1
	mov	DWORD PTR _GRAD, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7T1
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _T2
	mov	DWORD PTR _T1, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7T2
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _T3
	mov	DWORD PTR _T2, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7T3
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _AA1
	mov	DWORD PTR _T3, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7AA1
	push	ecx
	call	__Deallocate_i
	mov	DWORD PTR _.T281_$[esp+112], esi
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _MMAX
	mov	DWORD PTR _AA1, esi

; 196  : 	-			GRAD,T1,T2,T3)
; 197  :       ALLOCATE(	AA1(MMAX),RX1(MMAX),RY1(MMAX),RZ1(MMAX),

	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7RX1
	call	_make_info
	lea	ecx, DWORD PTR _.T281_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7RX1
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _MMAX
	mov	DWORD PTR _RX1, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7RY1
	call	_make_info
	lea	ecx, DWORD PTR _.T281_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7RY1
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _MMAX
	mov	DWORD PTR _RY1, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7RZ1
	call	_make_info
	lea	ecx, DWORD PTR _.T281_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7RZ1
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _MMAX
	mov	DWORD PTR _RZ1, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7AA2
	call	_make_info
	lea	ecx, DWORD PTR _.T281_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7AA2
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _MMAX
	mov	DWORD PTR _AA2, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7RX2
	call	_make_info
	lea	ecx, DWORD PTR _.T281_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7RX2
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _MMAX
	mov	DWORD PTR _RX2, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7RY2
	call	_make_info
	lea	ecx, DWORD PTR _.T281_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7RY2
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _MMAX
	mov	DWORD PTR _RY2, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7RZ2
	call	_make_info
	lea	ecx, DWORD PTR _.T281_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7RZ2
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _MMAX
	mov	DWORD PTR _RZ2, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7B12
	call	_make_info
	lea	ecx, DWORD PTR _.T281_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7B12
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _B12, eax
	mov	ecx, DWORD PTR _MMAX
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7EVEC
	call	_make_info
	lea	ecx, DWORD PTR _.T281_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7EVEC
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _EVEC, eax
	mov	eax, DWORD PTR _MMAX
	push	eax
	push	1
	push	eax
	push	1
	push	2
	push	OFFSET FLAT:_7S
	call	_make_info
	lea	ecx, DWORD PTR _.T281_$[esp+120]
	add	esp, 24					; 00000018H
	push	ecx
	push	8
	push	OFFSET FLAT:_7S
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _S, eax
	mov	eax, DWORD PTR _MMAX
	push	eax
	push	1
	push	eax
	push	1
	push	2
	push	OFFSET FLAT:_7H
	call	_make_info
	lea	ecx, DWORD PTR _.T281_$[esp+120]
	add	esp, 24					; 00000018H
	push	ecx
	push	8
	push	OFFSET FLAT:_7H
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _MMAX
	mov	DWORD PTR _H, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7CO
	call	_make_info
	lea	ecx, DWORD PTR _.T281_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7CO
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _MMAX
	mov	DWORD PTR _CO, eax
	lea	eax, DWORD PTR [ecx+ecx*2]
	add	eax, eax
	push	eax
	push	1
	push	1
	push	OFFSET FLAT:_7GRAD
	call	_make_info
	lea	ecx, DWORD PTR _.T281_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7GRAD
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _MMAX
	mov	DWORD PTR _GRAD, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7T1
	call	_make_info
	lea	ecx, DWORD PTR _.T281_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7T1
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _MMAX
	mov	DWORD PTR _T1, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7T2
	call	_make_info
	lea	ecx, DWORD PTR _.T281_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7T2
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _MMAX
	mov	DWORD PTR _T2, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7T3
	call	_make_info
	lea	ecx, DWORD PTR _.T281_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7T3
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _MMAX
	mov	DWORD PTR _T3, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7AA1
	call	_make_info
	lea	ecx, DWORD PTR _.T281_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7AA1
	call	__Allocate_i
	mov	ecx, DWORD PTR _.T281_$[esp+108]
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _AA1, eax
	mov	DWORD PTR _ERROR, ecx

; 198  : 	-			AA2(MMAX),RX2(MMAX),RY2(MMAX),RZ2(MMAX),
; 199  : 	-			B12(MMAX),EVEC(MMAX),S(MMAX,MMAX),H(MMAX,MMAX),
; 200  :      -			CO(MMAX),GRAD(MMAX*6),T1(MMAX),T2(MMAX),
; 201  :      -			T3(MMAX),STAT=ERROR)
; 202  : C
; 203  : 	IF (ERROR .NE. 0) STOP '!NOT ENOUGH MEMORY SPACE AVAILABLE ON THIS

	test	ecx, ecx
	je	SHORT $L282
	push	OFFSET FLAT:_.SG284_
	push	60					; 0000003cH
	call	__FFsSTOP
	add	esp, 8
$L282:

; 204  : 	- BOX! ...ABORTING!'
; 205  : 	DO J=1,M

	mov	esi, DWORD PTR _M
	mov	edi, 1
	test	esi, esi
	jle	$L286
	xor	ebx, ebx
$L98308:

; 206  : 		READ(5,*) I,EVEC(J),T1(J),T2(J),T3(J),RZ1(J),RZ2(J)

	push	ebx
	dec	esi
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	ebx
	push	5
	call	_f90_read
	lea	eax, DWORD PTR _I$[esp+124]
	add	esp, 28					; 0000001cH
	push	eax
	call	_lread_i
	add	esp, 4
	mov	eax, DWORD PTR _7EVEC+8
	add	eax, edi
	mov	ecx, DWORD PTR _EVEC
	lea	eax, DWORD PTR [ecx+eax*8]
	push	eax
	call	_lread_d
	add	esp, 4
	mov	ecx, DWORD PTR _7T1+8
	add	ecx, edi
	mov	eax, DWORD PTR _T1
	lea	edx, DWORD PTR [eax+ecx*8]
	push	edx
	call	_lread_d
	add	esp, 4
	mov	eax, DWORD PTR _7T2+8
	add	eax, edi
	mov	ecx, DWORD PTR _T2
	lea	eax, DWORD PTR [ecx+eax*8]
	push	eax
	call	_lread_d
	add	esp, 4
	mov	ecx, DWORD PTR _7T3+8
	add	ecx, edi
	mov	eax, DWORD PTR _T3
	lea	edx, DWORD PTR [eax+ecx*8]
	push	edx
	call	_lread_d
	add	esp, 4
	mov	eax, DWORD PTR _7RZ1+8
	add	eax, edi
	mov	ecx, DWORD PTR _RZ1
	lea	eax, DWORD PTR [ecx+eax*8]
	push	eax
	call	_lread_d
	add	esp, 4
	mov	eax, DWORD PTR _7RZ2+8
	add	eax, edi
	mov	ecx, DWORD PTR _RZ2
	lea	eax, DWORD PTR [ecx+eax*8]
	push	eax
	call	_lread_d
	add	esp, 4
	call	_f90_endread

; 207  : 	RX1(J)=0.00D+00

	mov	ecx, DWORD PTR _7RX1+8
	mov	eax, DWORD PTR _RX1
	add	ecx, edi
	mov	DWORD PTR [eax+ecx*8], ebx
	mov	DWORD PTR [eax+ecx*8+4], ebx

; 208  : 	RX2(J)=0.00D+00

	mov	ecx, DWORD PTR _7RX2+8
	add	ecx, edi
	mov	eax, DWORD PTR _RX2
	mov	DWORD PTR [eax+ecx*8], ebx
	mov	DWORD PTR [eax+ecx*8+4], ebx

; 209  : 	RY1(J)=0.00D+00

	mov	eax, DWORD PTR _7RY1+8
	add	eax, edi
	mov	ecx, DWORD PTR _RY1
	mov	DWORD PTR [ecx+eax*8], ebx
	mov	DWORD PTR [ecx+eax*8+4], ebx

; 210  : 	RY1(J)=0.00D+00

	mov	eax, DWORD PTR _7RY1+8
	add	eax, edi
	inc	edi
	test	esi, esi
	mov	DWORD PTR [ecx+eax*8], ebx
	mov	DWORD PTR [ecx+eax*8+4], ebx

; 204  : 	- BOX! ...ABORTING!'
; 205  : 	DO J=1,M

	jne	$L98308
$L286:

; 211  : 	END DO
; 212  : 	
; 213  :      	CLOSE (5)

	push	0
	push	0
	push	0
	push	5
	call	_f90_close
	add	esp, 16					; 00000010H

; 214  : 	OPEN (8, FILE = RESTARTOUTPUT, ACCESS = 'SEQUENTIAL',

	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	10					; 0000000aH
	push	OFFSET FLAT:_.SG290_
	push	7
	push	OFFSET FLAT:_.SG289_
	push	66					; 00000042H
	push	OFFSET FLAT:_RESTARTOUTPUT
	push	8
	call	_f90_open
	add	esp, 124				; 0000007cH
	mov	eax, DWORD PTR _RESTART

; 215  :      -		STATUS = 'UNKNOWN')
; 216  : C     -----------------------------------------------------
; 217  : C
; 218  : C     HAMILTONIAN SUBROUTINE:  Calcs. overlap & energy...
; 219  : C     -----------------------------------------------------
; 220  : 
; 221  : 
; 222  : 	IF (RESTART .EQ. 1) GOTO 323

	cmp	eax, 1
	je	$L98421
	jmp	SHORT $L267
$L98342:

; 180  : 	CALL GUESS(METHODGEN,M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,EVEC)

	mov	eax, DWORD PTR _EVEC
	mov	ecx, DWORD PTR _B12
	push	eax
	mov	edx, DWORD PTR _RZ2
	push	ecx
	mov	eax, DWORD PTR _RY2
	push	edx
	mov	ecx, DWORD PTR _RX2
	push	eax
	mov	edx, DWORD PTR _AA2
	push	ecx
	mov	eax, DWORD PTR _RZ1
	push	edx
	mov	ecx, DWORD PTR _RY1
	push	eax
	mov	edx, DWORD PTR _RX1
	push	ecx
	mov	eax, DWORD PTR _AA1
	push	edx
	push	eax
	push	OFFSET FLAT:_M
	push	OFFSET FLAT:_METHODGEN
	call	_GUESS@48
$L267:

; 223  : C	WRITE(6,*)'Ts ALTERED'
; 224  :  313	DO I=1,M

	mov	edx, DWORD PTR _M
	mov	esi, DWORD PTR _T3
	mov	DWORD PTR _I$[esp+96], 1
	test	edx, edx
	jle	$L98420
$L98309:

; 225  : 			T1(I)=AA1(I)	

	mov	edi, DWORD PTR _7AA1+8
	mov	eax, DWORD PTR _I$[esp+96]
	add	edi, eax
	mov	ebx, DWORD PTR _7T1+8
	mov	eax, DWORD PTR _AA1
	mov	ecx, DWORD PTR _I$[esp+96]
	add	ebx, ecx
	dec	edx
	mov	ebp, DWORD PTR [eax+edi*8+4]
	mov	ecx, DWORD PTR _T1
	mov	eax, DWORD PTR [eax+edi*8]
	mov	edi, DWORD PTR _T2
	mov	DWORD PTR [ecx+ebx*8+4], ebp
	mov	DWORD PTR [ecx+ebx*8], eax

; 226  : 			T2(I)=AA2(I)	

	mov	ecx, DWORD PTR _7AA2+8
	mov	eax, DWORD PTR _I$[esp+96]
	mov	ebx, DWORD PTR _7T2+8
	add	ecx, eax
	mov	ebp, DWORD PTR _I$[esp+96]
	mov	eax, DWORD PTR _AA2
	add	ebx, ebp
	mov	ebp, DWORD PTR [eax+ecx*8+4]
	mov	eax, DWORD PTR [eax+ecx*8]
	mov	ecx, DWORD PTR _B12
	mov	DWORD PTR [edi+ebx*8+4], ebp
	mov	DWORD PTR [edi+ebx*8], eax

; 227  : 			T3(I)=B12(I)	

	mov	edi, DWORD PTR _7B12+8
	mov	eax, DWORD PTR _I$[esp+96]
	mov	ebx, DWORD PTR _7T3+8
	add	edi, eax
	add	ebx, eax
	mov	ebp, DWORD PTR [ecx+edi*8+4]
	mov	ecx, DWORD PTR [ecx+edi*8]
	mov	DWORD PTR [esi+ebx*8+4], ebp
	mov	DWORD PTR [esi+ebx*8], ecx

; 223  : C	WRITE(6,*)'Ts ALTERED'
; 224  :  313	DO I=1,M

	mov	eax, DWORD PTR _I$[esp+96]
	inc	eax
	test	edx, edx
	mov	DWORD PTR _I$[esp+96], eax
	jne	$L98309
$L98420:

; 228  : 	END DO
; 229  : C
; 230  : 	EVEC(1) = 4.919405348897607E-001

	mov	ecx, DWORD PTR _7EVEC+8
	mov	eax, DWORD PTR _EVEC
	mov	DWORD PTR [eax+ecx*8+8], 658204867	; 273b68c3H
	mov	DWORD PTR [eax+ecx*8+12], 1071610868	; 3fdf7bf4H

; 231  : 	EVEC(2) = 1.251354302680531D+00

	mov	ecx, DWORD PTR _7EVEC+8
	mov	DWORD PTR [eax+ecx*8+16], 383487066	; 16db8c5aH
	mov	DWORD PTR [eax+ecx*8+20], 1072956812	; 3ff4058cH

; 232  : 
; 233  : 	T1(1) = 5.927253752706102E-001

	mov	eax, DWORD PTR _7T1+8
	mov	ecx, DWORD PTR _T1
	mov	DWORD PTR [ecx+eax*8+8], 885620155	; 34c97dbbH
	mov	DWORD PTR [ecx+eax*8+12], 1071839131	; 3fe2f79bH

; 234  : 	T1(2) = 7.395751152915679E-001

	mov	eax, DWORD PTR _7T1+8
	mov	DWORD PTR [ecx+eax*8+16], 1856215898	; 6ea39b5aH
	mov	DWORD PTR [ecx+eax*8+20], 1072147097	; 3fe7aa99H

; 235  : 
; 236  : 	T2(1) = 1.661941607581848D+00

	mov	ecx, DWORD PTR _7T2+8
	mov	eax, DWORD PTR _T2
	mov	DWORD PTR [eax+ecx*8+8], 356962222	; 1546cfaeH
	mov	DWORD PTR [eax+ecx*8+12], 1073387344	; 3ffa9750H

; 237  : 	T2(2) = 4.015065938655144E-001

	mov	ecx, DWORD PTR _7T2+8
	mov	DWORD PTR [eax+ecx*8+16], -1234045025	; b671f79fH
	mov	DWORD PTR [eax+ecx*8+20], 1071231560	; 3fd9b248H

; 238  : 
; 239  : 	T3(1) = 1.097749060489904E-001

	mov	eax, DWORD PTR _7T3+8
	mov	ecx, DWORD PTR _RZ1
	mov	DWORD PTR [esi+eax*8+8], 1332142598	; 4f66e206H
	mov	DWORD PTR [esi+eax*8+12], 1069292085	; 3fbc1a35H

; 240  : 	T3(2) = 4.130469862151015E-002

	mov	eax, DWORD PTR _7T3+8
	mov	DWORD PTR [esi+eax*8+16], -1283439496	; b3804478H
	mov	DWORD PTR [esi+eax*8+20], 1067787747	; 3fa525e3H

; 241  : 
; 242  : 	RZ1(1) = 2.161923082544148E-001

	mov	eax, DWORD PTR _7RZ1+8
	mov	DWORD PTR [ecx+eax*8+8], -2033403089	; 86ccbb2fH
	mov	DWORD PTR [ecx+eax*8+12], 1070312496	; 3fcbac30H

; 243  : 	RZ1(2) = 5.796860175253535E-001

	mov	eax, DWORD PTR _7RZ1+8
	mov	DWORD PTR [ecx+eax*8+16], -1327036102	; b0e7093aH
	mov	DWORD PTR [ecx+eax*8+20], 1071811785	; 3fe28cc9H

; 244  : 
; 245  : 	RZ2(1) = -6.289395277703103E-001

	mov	ecx, DWORD PTR _7RZ2+8
	mov	eax, DWORD PTR _RZ2
	mov	DWORD PTR [eax+ecx*8+8], -908202782	; c9ddece2H
	mov	DWORD PTR [eax+ecx*8+12], -1075568571	; bfe42045H

; 246  : 	RZ2(2) = -9.009920951052639E-002

	mov	ecx, DWORD PTR _7RZ2+8
	xor	edx, edx
	mov	DWORD PTR [eax+ecx*8+16], -432128018	; e63e3feeH
	mov	DWORD PTR [eax+ecx*8+20], -1078521667	; bfb710bdH

; 247  : 
; 248  : 	RX1(1) = 0.00D+00

	mov	ecx, DWORD PTR _7RX1+8
	mov	eax, DWORD PTR _RX1
	mov	DWORD PTR [eax+ecx*8+8], edx
	mov	DWORD PTR [eax+ecx*8+12], edx

; 249  : 	RX1(2) = 0.00D+00

	mov	ecx, DWORD PTR _7RX1+8
	mov	DWORD PTR [eax+ecx*8+16], edx
	mov	DWORD PTR [eax+ecx*8+20], edx

; 250  : 
; 251  : 	RX2(1) = 0.00D+00

	mov	eax, DWORD PTR _7RX2+8
	mov	ecx, DWORD PTR _RX2
	mov	DWORD PTR [ecx+eax*8+8], edx
	mov	DWORD PTR [ecx+eax*8+12], edx

; 252  : 	RX2(2) = 0.00D+00

	mov	eax, DWORD PTR _7RX2+8
	mov	DWORD PTR [ecx+eax*8+16], edx
	mov	DWORD PTR [ecx+eax*8+20], edx
	jmp	SHORT $L293
$L98421:
	mov	esi, DWORD PTR _T3
$L293:

; 253  : C
; 254  :  323  DO I=1,M

	mov	DWORD PTR _I$[esp+96], 1
	mov	edi, DWORD PTR _M
	test	edi, edi
	jle	$L296
$L98310:
	mov	eax, DWORD PTR _7T1+8
	mov	ecx, DWORD PTR _I$[esp+96]
	add	eax, ecx
	mov	edx, DWORD PTR _T1
	mov	ecx, DWORD PTR _7T3+8
	mov	ebx, DWORD PTR _I$[esp+96]
	lea	eax, DWORD PTR [edx+eax*8]
	add	ecx, ebx

; 255  : 	    AA1(I)=T1(I)**2 + T1(I)*T3(I)

	fld	QWORD PTR [esi+ecx*8]
	fadd	QWORD PTR [eax]
	mov	ecx, DWORD PTR _AA1
	fmul	QWORD PTR [eax]
	mov	eax, DWORD PTR _7AA1+8
	add	eax, ebx
	fstp	QWORD PTR [ecx+eax*8]
	mov	eax, DWORD PTR _7T3+8
	mov	ebx, DWORD PTR _I$[esp+96]
	add	eax, ebx
	mov	edx, DWORD PTR _7T1+8
	add	edx, ebx
	mov	ecx, DWORD PTR _T1
	lea	eax, DWORD PTR [esi+eax*8]

; 256  :           AA2(I)=T2(I)**2 + T3(I)**2 + T1(I)*T3(I)

	fld	QWORD PTR [ecx+edx*8]
	fadd	QWORD PTR [eax]
	mov	ecx, DWORD PTR _7T2+8
	add	ecx, ebx
	fmul	QWORD PTR [eax]
	mov	eax, DWORD PTR _T2
	fld	QWORD PTR [eax+ecx*8]
	fld	QWORD PTR $T98441
	call	__FIfexp
	faddp	ST(1), ST(0)
	mov	ecx, DWORD PTR _7AA2+8
	mov	ebx, DWORD PTR _I$[esp+96]
	add	ecx, ebx
	mov	eax, DWORD PTR _AA2
	push	0
	push	0
	fstp	QWORD PTR [eax+ecx*8]

; 257  :           B12(I)=-T1(I)*T3(I)

	mov	eax, DWORD PTR _7T3+8
	mov	ebx, DWORD PTR _I$[esp+104]
	add	eax, ebx
	mov	ecx, DWORD PTR _T1
	push	0
	push	0
	fld	QWORD PTR [esi+eax*8]
	mov	eax, DWORD PTR _7T1+8
	push	0
	add	eax, ebx
	push	6
	fmul	QWORD PTR [ecx+eax*8]
	mov	eax, DWORD PTR _7B12+8
	mov	ecx, DWORD PTR _B12
	add	eax, ebx
	fchs
	fstp	QWORD PTR [ecx+eax*8]

; 258  : 
; 259  : 	write(6,*)I,AA1(I),AA2(I),B12(I),RZ1(I),RZ2(I)

	call	_f90_write
	mov	ebx, DWORD PTR _I$[esp+120]
	add	esp, 24					; 00000018H
	push	ebx
	call	_lwrite_i
	mov	ebx, DWORD PTR _I$[esp+100]
	add	esp, 4
	mov	eax, DWORD PTR _7AA1+8
	mov	ecx, DWORD PTR _AA1
	add	eax, ebx
	mov	edx, DWORD PTR [ecx+eax*8+4]
	mov	eax, DWORD PTR [ecx+eax*8]
	push	edx
	push	eax
	call	_lwrite_d
	mov	ebx, DWORD PTR _I$[esp+104]
	add	esp, 8
	mov	ecx, DWORD PTR _7AA2+8
	mov	eax, DWORD PTR _AA2
	add	ecx, ebx
	mov	edx, DWORD PTR [eax+ecx*8+4]
	mov	ecx, DWORD PTR [eax+ecx*8]
	push	edx
	push	ecx
	call	_lwrite_d
	mov	ebx, DWORD PTR _I$[esp+104]
	add	esp, 8
	mov	ecx, DWORD PTR _7B12+8
	mov	eax, DWORD PTR _B12
	add	ecx, ebx
	mov	edx, DWORD PTR [eax+ecx*8+4]
	mov	ecx, DWORD PTR [eax+ecx*8]
	push	edx
	push	ecx
	call	_lwrite_d
	mov	ebx, DWORD PTR _I$[esp+104]
	add	esp, 8
	mov	edx, DWORD PTR _7RZ1+8
	mov	eax, DWORD PTR _RZ1
	add	edx, ebx
	mov	ecx, DWORD PTR [eax+edx*8+4]
	mov	edx, DWORD PTR [eax+edx*8]
	push	ecx
	push	edx
	call	_lwrite_d
	mov	ebx, DWORD PTR _I$[esp+104]
	add	esp, 8
	mov	eax, DWORD PTR _7RZ2+8
	mov	ecx, DWORD PTR _RZ2
	add	eax, ebx
	dec	edi
	mov	edx, DWORD PTR [ecx+eax*8+4]
	mov	eax, DWORD PTR [ecx+eax*8]
	push	edx
	push	eax
	call	_lwrite_d
	add	esp, 8
	call	_f90_endwrite

; 253  : C
; 254  :  323  DO I=1,M

	mov	ebx, DWORD PTR _I$[esp+96]
	mov	esi, DWORD PTR _T3
	inc	ebx
	test	edi, edi
	mov	DWORD PTR _I$[esp+96], ebx
	jne	$L98310
$L296:

; 260  : C
; 261  : C		WRITE(6,292)I,AA1(I),AA2(I),RZ1(I),RZ2(I)
; 262  : 	END DO
; 263  : C
; 264  : C	STOP
; 265  : C
; 266  : 	CALL BIOUT(' ')

	push	1
	push	OFFSET FLAT:_.SG299_
	call	_BIOUT@8

; 267  : C	 CALL BIOUT('BEGINNING TIMING ROUTINE...Please Stand By...')
; 268  : C	 CALL BIOUT('---------------------------------------------')
; 269  : C      CALL GETTIM(tmphour, tmpminute, tmpsecond, tmphund)  
; 270  : C	 CALL BIOUT(':BEGIN:')
; 271  : C      WRITE (7, 9004) tmphour, tmpminute, tmpsecond, tmphund
; 272  : C      WRITE (6, 9004) tmphour, tmpminute, tmpsecond, tmphund
; 273  : 	DFLAG=0

	xor	eax, eax
	mov	ecx, DWORD PTR _PREOPT
	mov	DWORD PTR _DFLAG, eax
	cmp	ecx, 1

; 274  : 	IF (PREOPT .EQ. 1) THEN DFLAG=1

	jne	SHORT $L301
	mov	DWORD PTR _THENDFLAG, eax
	mov	DWORD PTR _THENDFLAG+4, 1072693248	; 3ff00000H
$L301:

; 275  : C	 CALL BIOUT(':FUNCTION CALCULATION IN PROGRESS:')
; 276  : 	CALL HAMILTON(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,EVEC,S,H,

	push	OFFSET FLAT:_DFLAG
	mov	eax, DWORD PTR _T3
	push	eax
	mov	ecx, DWORD PTR _T2
	push	ecx
	mov	edx, DWORD PTR _T1
	push	edx
	mov	eax, DWORD PTR _GRAD
	push	eax
	mov	ecx, DWORD PTR _CO
	push	OFFSET FLAT:_EMIN
	mov	edx, DWORD PTR _H
	push	ecx
	mov	eax, DWORD PTR _S
	push	edx
	mov	ecx, DWORD PTR _EVEC
	push	eax
	mov	edx, DWORD PTR _B12
	push	ecx
	mov	eax, DWORD PTR _RZ2
	push	edx
	mov	ecx, DWORD PTR _RY2
	push	eax
	mov	edx, DWORD PTR _RX2
	push	ecx
	mov	eax, DWORD PTR _AA2
	push	edx
	mov	ecx, DWORD PTR _RZ1
	push	eax
	mov	edx, DWORD PTR _RY1
	push	ecx
	mov	eax, DWORD PTR _RX1
	push	edx
	mov	ecx, DWORD PTR _AA1
	push	eax
	push	ecx
	push	OFFSET FLAT:_M
	call	_HAMILTON@80

; 277  : 	-				CO,EMIN,GRAD,T1,T2,T3,DFLAG)
; 278  : 
; 279  : 	WRITE(6,*)'ENERGY=',EMIN

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	7
	push	OFFSET FLAT:_.SG306_
	call	_lwrite_ch
	add	esp, 8
	mov	ecx, DWORD PTR _EMIN+4
	mov	edx, DWORD PTR _EMIN
	npad	1
	push	ecx
	push	edx
	call	_lwrite_d
	add	esp, 8
	call	_f90_endwrite

; 280  : 	STOP

	push	0
	push	0
	call	__FFsSTOP
	add	esp, 8
	mov	ecx, DWORD PTR _PREOPT

; 281  : C
; 282  : 	IF (PREOPT .EQ. 1) THEN

	cmp	ecx, 1
	jne	$L307

; 283  : 			CALL GRADOPT(M,T1,RX1,RY1,RZ1,T2,RX2,RY2,RZ2,T3,EVEC,GRAD)

	mov	eax, DWORD PTR _GRAD
	mov	ecx, DWORD PTR _EVEC
	push	eax
	mov	edx, DWORD PTR _T3
	push	ecx
	mov	eax, DWORD PTR _RZ2
	push	edx
	mov	ecx, DWORD PTR _RY2
	push	eax
	mov	edx, DWORD PTR _RX2
	push	ecx
	mov	eax, DWORD PTR _T2
	push	edx
	mov	ecx, DWORD PTR _RZ1
	push	eax
	mov	edx, DWORD PTR _RY1
	push	ecx
	mov	eax, DWORD PTR _RX1
	push	edx
	mov	ecx, DWORD PTR _T1
	push	eax
	mov	ebp, 1
	push	ecx
	push	OFFSET FLAT:_M
	call	_GRADOPT@48

; 284  : 	DO NNN=1,M

	mov	esi, DWORD PTR _M
	mov	edi, DWORD PTR _T1
	mov	ebx, DWORD PTR _T3
	test	esi, esi
	jle	$L98422
	jmp	SHORT $L98311
$L98466:
	fstp	QWORD PTR [ecx+eax*8]
$L98311:
	mov	eax, DWORD PTR _7T1+8
	mov	ecx, DWORD PTR _7T3+8
	add	eax, ebp
	add	ecx, ebp
	dec	esi
	fld	QWORD PTR [ebx+ecx*8]
	fadd	QWORD PTR [edi+eax*8]
	lea	eax, DWORD PTR [edi+eax*8]
	mov	ecx, DWORD PTR _AA1

; 285  : 		AA1(NNN)=T1(NNN)**2 +T1(NNN)*T3(NNN)

	fmul	QWORD PTR [eax]
	mov	eax, DWORD PTR _7AA1+8
	add	eax, ebp
	fstp	QWORD PTR [ecx+eax*8]
	mov	eax, DWORD PTR _7T3+8
	mov	ecx, DWORD PTR _7T1+8
	add	eax, ebp
	add	ecx, ebp
	fld	QWORD PTR [edi+ecx*8]
	fadd	QWORD PTR [ebx+eax*8]
	lea	eax, DWORD PTR [ebx+eax*8]
	mov	ecx, DWORD PTR _7T2+8

; 286  : 		AA2(NNN)=T2(NNN)**2 + T3(NNN)**2 +T1(NNN)*T3(NNN)

	add	ecx, ebp
	fmul	QWORD PTR [eax]
	mov	eax, DWORD PTR _T2
	fld	QWORD PTR [eax+ecx*8]
	fld	QWORD PTR $T98441
	call	__FIfexp
	faddp	ST(1), ST(0)
	mov	ecx, DWORD PTR _7AA2+8
	mov	eax, DWORD PTR _AA2
	add	ecx, ebp
	fstp	QWORD PTR [eax+ecx*8]

; 287  : 		B12(NNN)=-T1(NNN)*T3(NNN)

	mov	eax, DWORD PTR _7T1+8
	mov	ecx, DWORD PTR _B12
	add	eax, ebp
	fld	QWORD PTR [edi+eax*8]
	mov	eax, DWORD PTR _7T3+8
	add	eax, ebp
	fmul	QWORD PTR [ebx+eax*8]
	mov	eax, DWORD PTR _7B12+8
	add	eax, ebp
	inc	ebp
	fchs
	test	esi, esi

; 284  : 	DO NNN=1,M

	jne	$L98466
	fstp	QWORD PTR [ecx+eax*8]
$L98422:

; 288  : 	END DO
; 289  :       DEFLAG=0

	xor	eax, eax
	push	OFFSET FLAT:_DFLAG
	push	ebx
	mov	ecx, DWORD PTR _T2
	push	ecx
	mov	edx, DWORD PTR _GRAD
	push	edi
	mov	DWORD PTR _DEFLAG, eax
	push	edx
	mov	ecx, DWORD PTR _H
	push	OFFSET FLAT:_EMIN
	mov	edx, DWORD PTR _S
	mov	DWORD PTR _DEFLAG+4, eax

; 290  : 	CALL HAMILTON(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,EVEC,S,H,

	mov	eax, DWORD PTR _CO
	push	eax
	push	ecx
	push	edx
	mov	eax, DWORD PTR _EVEC
	push	eax
	mov	ecx, DWORD PTR _B12
	push	ecx
	mov	edx, DWORD PTR _RZ2
	push	edx
	mov	eax, DWORD PTR _RY2
	push	eax
	mov	ecx, DWORD PTR _RX2
	push	ecx
	mov	edx, DWORD PTR _AA2
	push	edx
	mov	eax, DWORD PTR _RZ1
	push	eax
	mov	ecx, DWORD PTR _RY1
	push	ecx
	mov	edx, DWORD PTR _RX1
	push	edx
	mov	eax, DWORD PTR _AA1
	push	eax
	push	OFFSET FLAT:_M
	call	_HAMILTON@80
$L307:

; 291  : 	-				CO,EMIN,GRAD,T1,T2,T3,DFLAG)
; 292  : 	END IF
; 293  : C
; 294  : 	ENER=EMIN

	mov	eax, DWORD PTR _EMIN+4
	mov	ecx, DWORD PTR _EMIN
	push	0
	mov	DWORD PTR _ENER+4, eax
	push	0
	mov	DWORD PTR _ENER, ecx

; 295  : 	WRITE(6,*)'STARTING WITH ',M,' FUNCTIONS...'

	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	14					; 0000000eH
	push	OFFSET FLAT:_.SG315_
	call	_lwrite_ch
	add	esp, 8
	mov	ecx, DWORD PTR _M
	push	ecx
	call	_lwrite_i
	add	esp, 4
	push	13					; 0000000dH
	push	OFFSET FLAT:_.SG316_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 296  : 	WRITE(6,*)'STARTING ENERGY=',ENER

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	16					; 00000010H
	push	OFFSET FLAT:_.SG317_
	call	_lwrite_ch
	add	esp, 8
	mov	ecx, DWORD PTR _ENER+4
	mov	edx, DWORD PTR _ENER
	mov	esi, 1
	push	ecx
	push	edx
	call	_lwrite_d
	add	esp, 8
	call	_f90_endwrite

; 297  : C
; 298  : 	DO J=1,M

	mov	edx, DWORD PTR _M
	test	edx, edx
	jle	SHORT $L318
$L98312:

; 299  : 		EVEC(J)=CO(J)

	mov	ecx, DWORD PTR _7CO+8
	mov	ebp, DWORD PTR _CO
	add	ecx, esi
	mov	edi, DWORD PTR _7EVEC+8
	add	edi, esi
	mov	ebx, DWORD PTR _EVEC
	mov	eax, DWORD PTR [ebp+ecx*8+4]
	dec	edx
	mov	ecx, DWORD PTR [ebp+ecx*8]
	inc	esi
	mov	DWORD PTR [ebx+edi*8+4], eax
	test	edx, edx
	mov	DWORD PTR [ebx+edi*8], ecx

; 297  : C
; 298  : 	DO J=1,M

	jne	SHORT $L98312
$L318:

; 300  : 	END DO
; 301  : C
; 302  : 	XRX1=0.00D+00

	mov	DWORD PTR _XRX1, 0
	mov	DWORD PTR _XRX1+4, 0

; 303  : 	XRY1=0.00D+00

	mov	DWORD PTR _XRY1, 0
	mov	DWORD PTR _XRY1+4, 0

; 304  : 	XRX2=0.00D+00

	mov	DWORD PTR _XRX2, 0
	mov	DWORD PTR _XRX2+4, 0

; 305  : 	XRY2=0.00D+00

	mov	DWORD PTR _XRY2, 0
	mov	DWORD PTR _XRY2+4, 0

; 306  : C
; 307  : C      CALL GETTIM(tmphour, tmpminute, tmpsecond, tmphund)  
; 308  : C      WRITE (7, 9004) tmphour, tmpminute, tmpsecond, tmphund
; 309  : C      WRITE (6, 9004) tmphour, tmpminute, tmpsecond, tmphund
; 310  : C
; 311  : C	 DFLAG=1
; 312  : C	 CALL BIOUT(':GRADIENT CALCULATION IN PROGRESS:')
; 313  : C	 CALL HAMILTON(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,EVEC,S,H,
; 314  : C	-				CO,EMIN,GRAD,T1,T2,T3,DFLAG)
; 315  : C      CALL GETTIM(tmphour, tmpminute, tmpsecond, tmphund)  
; 316  : C      WRITE (7, 9004) tmphour, tmpminute, tmpsecond, tmphund
; 317  : C      WRITE (6, 9004) tmphour, tmpminute, tmpsecond, tmphund
; 318  : CC
; 319  : C	CALL BIOUT(':END:')
; 320  : C	CALL BIOUT(' ')
; 321  : C     -----------------------------------------------------
; 322  : C      WRITE(6,*) '           UNOPTIMIZED ENERGY = ',EMIN
; 323  : C      WRITE(7,*) '           UNOPTIMIZED ENERGY = ',EMIN
; 324  : C      CALL BIOUT(' ')
; 325  : C      WRITE (7,*) ' UNOPTIMIZED ENERGY ROUTINE FINISHES:  TIME STAMP...'
; 326  : C      WRITE (6,*) ' UNOPTIMIZED ENERGY ROUTINE FINISHES:  TIME STAMP...'
; 327  : C
; 328  : C	CALL BIOUT(' ')
; 329  : C      CALL BIOUT('     LinCoeff:.   AA1:        AA2:        B12:        
; 330  : C	-RZ1:       RZ2:')
; 331  : C      CALL BIOUT('FNC: =================================================
; 332  : C     -=====================')
; 333  : C      DO J=1,M
; 334  : C          WRITE(6,1111) J,EVEC(J),AA1(J),AA2(J),B12(J),RZ1(J),RZ2(J)
; 335  : C          WRITE(7,1111) J,EVEC(J),AA1(J),AA2(J),B12(J),RZ1(J),RZ2(J)
; 336  : C	END DO
; 337  : C
; 338  : C
; 339  : 	CALL SEED(RND$TIMESEED)

	mov	DWORD PTR $T98313[esp+96], -1
	lea	eax, DWORD PTR $T98313[esp+96]
	push	eax
	call	_SEED@4

; 340  : C
; 341  : 	BD=1.4011D+00

	mov	ecx, -714682558				; d566cf42H
	mov	edx, -1071906816			; c01c0000H
	mov	DWORD PTR _FMAX, 0
	mov	DWORD PTR _FMAX+4, 1076101120		; 40240000H
	mov	DWORD PTR _FMIN, 0
	mov	DWORD PTR _BD, ecx
	mov	DWORD PTR _FMIN+4, edx
	mov	DWORD PTR _BD+4, 1073113831		; 3ff66ae7H

; 342  :       FMAX=10.00D+00
; 343  : 	FMIN=-7.00D+00
; 344  : 	ZMAX=1.70D+00

	mov	DWORD PTR _ZMAX, 858993459		; 33333333H
	mov	DWORD PTR _ZMAX+4, 1073427251		; 3ffb3333H

; 345  : 	BL0=0.7D+00

	mov	DWORD PTR _BL0, 1717986918		; 66666666H
	mov	DWORD PTR _BL0+4, 1072064102		; 3fe66666H

; 346  : 	BL1=1.3D+00

	mov	DWORD PTR _BL1, -858993459		; cccccccdH
	mov	DWORD PTR _BL1+4, 1073007820		; 3ff4ccccH

; 347  : 	BLURBEGIN=5.

	mov	DWORD PTR _BLURBEGIN, 0
	mov	DWORD PTR _BLUREND, 0
	mov	eax, 1075052544				; 40140000H
	mov	DWORD PTR _BLUREND+4, edx
	mov	DWORD PTR _BLURBEGIN+4, eax
	mov	DWORD PTR _ZNUC, ecx

; 348  : 	BLUREND=-7.
; 349  : 	BETAMAX=5.

	mov	DWORD PTR _BETAMAX, 0
	mov	DWORD PTR _BETAMIN, 0
	mov	DWORD PTR _BETAMIN+4, -1072431104	; c0140000H
	mov	DWORD PTR _IVIEW, 1000			; 000003e8H
	push	OFFSET FLAT:_TEMP
	mov	DWORD PTR _BETAMAX+4, eax

; 350  : 	BETAMIN=-5.
; 351  : C
; 352  : C	IVIEW=4000
; 353  : 	IVIEW=1000
; 354  : 	ZNUC=BD/2.00D+00

	mov	DWORD PTR _ZNUC+4, 1072065255		; 3fe66ae7H

; 355  : C
; 356  : 	A=FMAX-FMIN

	mov	DWORD PTR _A, 0
	mov	DWORD PTR _A+4, 1076953088		; 40310000H

; 357  : 	B=FMIN

	mov	DWORD PTR _B, 0
	mov	DWORD PTR _B+4, edx

; 358  : C
; 359  : 	CALL RANDOM(TEMP)

	call	_RANDOM@4

; 360  : 	TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fld	QWORD PTR $T98452
	fxch	ST(1)
	fstp	QWORD PTR _TEMP1
	call	__FIlog
	fstp	QWORD PTR -8+[esp+96]

; 361  : C
; 362  : 	K=-DLOG(0.0001)/(ZMAX-ZNUC)

	fld	QWORD PTR _ZMAX
	fsub	QWORD PTR _ZNUC
	fdivr	QWORD PTR -8+[esp+96]
	fchs
	call	__ftol
	fld	QWORD PTR -8+[esp+96]
	fdiv	QWORD PTR _ZNUC
	mov	DWORD PTR _K$[esp+96], eax

; 363  : 	KA=-DLOG(0.0001)/(ZNUC)

	fchs
	call	__ftol
	mov	ecx, DWORD PTR _STGROW
	mov	DWORD PTR _KA, eax

; 364  : C
; 365  : C
; 366  : C
; 367  : 	IF (STGROW .EQ. 0) THEN

	test	ecx, ecx
	jne	$L344

; 368  : C
; 369  : 	DO JJJ=M,ISTEND

	mov	eax, DWORD PTR _ISTEND
	mov	ecx, DWORD PTR _M
	sub	eax, ecx
	inc	eax
	mov	DWORD PTR -36+[esp+96], eax
	test	eax, eax
	jle	$L344
$L98315:

; 370  : 	DO I=1,IVIEW

	mov	DWORD PTR _I$[esp+96], 1
	mov	esi, DWORD PTR _IVIEW
	test	esi, esi
	jle	$L349
	xor	edi, edi
	mov	ebx, 1072693248				; 3ff00000H
	mov	ebp, -1074790400			; bff00000H
$L98316:

; 371  : 	IF (TEMP1.GE..5) THEN

	fld	QWORD PTR _TEMP1
	fcomp	QWORD PTR $T98454
	fwait
	fnstsw	ax
	test	ah, 1

; 372  : C
; 373  : 		CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	jne	$L352
	call	_RANDOM@4

; 374  : 		TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1
$L353:

; 375  :  4141		CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	call	_RANDOM@4

; 376  : 		TEMP2=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP2

; 377  : 	IF (TEMP2.LT..5) GOTO 4141

	fld	QWORD PTR _TEMP2
	fcomp	QWORD PTR $T98454
	fwait
	fnstsw	ax
	test	ah, 1
	jne	SHORT $L353

; 378  : 	TOF1=(FMAX-FMIN)*TEMP1+FMIN

	fld	QWORD PTR _FMAX
	fsub	QWORD PTR _FMIN
	fld	QWORD PTR $T98443
	fxch	ST(1)
	fmul	QWORD PTR _TEMP1
	fadd	QWORD PTR _FMIN
	fstp	QWORD PTR _TOF1

; 379  : 	XAA1=TEMP2*10**TOF1

	fld	QWORD PTR _TOF1
	call	__FIfexp
	fmul	QWORD PTR _TEMP2
	push	OFFSET FLAT:_TEMP
	fstp	QWORD PTR _XAA1

; 380  : C
; 381  : 		CALL RANDOM(TEMP)

	call	_RANDOM@4

; 382  : 		TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1
$L359:

; 383  :  5151		CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	call	_RANDOM@4

; 384  : 		TEMP2=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP2

; 385  : 	IF (TEMP2.LT..5) GOTO 5151

	fld	QWORD PTR _TEMP2
	fcomp	QWORD PTR $T98454
	fwait
	fnstsw	ax
	test	ah, 1
	jne	SHORT $L359
	fld	QWORD PTR _FMAX
	fsub	QWORD PTR _FMIN
	fstp	QWORD PTR -16+[esp+96]

; 386  : 	TOF2=(FMAX-FMIN)*TEMP1+FMIN

	fld	QWORD PTR -16+[esp+96]
	fmul	QWORD PTR _TEMP1
	fadd	QWORD PTR _FMIN
	fstp	QWORD PTR _TOF2

; 387  : 	XAA2=TEMP2*10**TOF2

	fld	QWORD PTR $T98443
	fld	QWORD PTR _TOF2
	call	__FIfexp
	fmul	QWORD PTR _TEMP2
	fstp	QWORD PTR _XAA2
	fld	QWORD PTR _FMIN
	fabs
	fstp	QWORD PTR -8+[esp+96]
	fild	DWORD PTR _K$[esp+96]
	fstp	QWORD PTR -24+[esp+96]

; 388  : C
; 389  : C	XXX1=(TOF1+DABS(FMIN)+1.+D-(DABS(FMIN)+1.))/C
; 390  : C	XXX2=(TOF2+DABS(FMIN)+1.+D-(DABS(FMIN)+1.))/C
; 391  : C	XRZ1=-(-SIGMAP *DLOG(XXX1))+ZNUC
; 392  : 	XRZ1=-DLOG((TOF1+DABS(FMIN))/(FMAX-FMIN))/K + ZNUC

	fld	QWORD PTR -8+[esp+96]
	fadd	QWORD PTR _TOF1
	fdiv	QWORD PTR -16+[esp+96]
	call	__FIlog
	fdiv	QWORD PTR -24+[esp+96]
	fsubr	QWORD PTR _ZNUC
	fstp	QWORD PTR _XRZ1

; 393  : C	XRZ2=-(-SIGMAP *DLOG(XXX2))+ZNUC
; 394  : 	XRZ2=-DLOG((TOF2+DABS(FMIN))/(FMAX-FMIN))/K + ZNUC

	fld	QWORD PTR _TOF2
	fadd	QWORD PTR -8+[esp+96]
	fdiv	QWORD PTR -16+[esp+96]
	call	__FIlog
	fdiv	QWORD PTR -24+[esp+96]
	push	OFFSET FLAT:_TEMP
	fsubr	QWORD PTR _ZNUC
	fstp	QWORD PTR _XRZ2

; 395  : 	CALL RANDOM(TEMP)

	call	_RANDOM@4

; 396  : 	TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1

; 397  : 	IF (TEMP1.GE..5) XRZ1=-1.*XRZ1

	fld	QWORD PTR _TEMP1
	fcomp	QWORD PTR $T98454
	fwait
	fnstsw	ax
	test	ah, 1
	jne	SHORT $L366
	fld	QWORD PTR _XRZ1
	fchs
	fstp	QWORD PTR _XRZ1
$L366:

; 398  : 	CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	call	_RANDOM@4

; 399  : 	TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1

; 400  : 	IF (TEMP1.GE..5) XRZ2=-1*XRZ2

	fld	QWORD PTR _TEMP1
	fcomp	QWORD PTR $T98454
	fwait
	fnstsw	ax
	test	ah, 1
	jne	SHORT $L368
	fld	QWORD PTR _XRZ2
	fchs
	fstp	QWORD PTR _XRZ2
$L368:

; 401  : 	IF (BL0.GT.0.) THEN

	fld	QWORD PTR _BL0
	fcomp	QWORD PTR $T98428
	fwait
	fnstsw	ax
	test	ah, 65					; 00000041H
	jne	$L351

; 402  : 		CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	call	_RANDOM@4

; 403  : 	TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1

; 404  : 		XC=1

	fld	QWORD PTR _TEMP1
	fcomp	QWORD PTR $T98454
	fwait
	fnstsw	ax
	test	ah, 1
	mov	DWORD PTR _XC, edi
	mov	DWORD PTR _XC+4, ebx

; 405  : 		IF (TEMP1.GE..5) XC=-1.

	jne	SHORT $L373
	mov	DWORD PTR _XC, edi
	mov	DWORD PTR _XC+4, ebp
$L373:

; 406  : 	CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	call	_RANDOM@4

; 407  : 		TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1

; 408  : 		IF (TOF1 .LE. BLURBEGIN) THEN

	fld	QWORD PTR _BLURBEGIN
	fcomp	QWORD PTR _TOF1
	fwait
	fnstsw	ax
	test	ah, 1
	jne	SHORT $L375

; 409  : 		XRZ1=XRZ1+XC*DSQRT(-2.*BL0**2 * DLOG(TEMP1))* 

	fld	QWORD PTR _BLURBEGIN
	fsub	QWORD PTR _TOF1
	fld	QWORD PTR _TEMP1
	call	__FIlog
	fld	QWORD PTR _BL0
	fld	QWORD PTR $T98441
	call	__FIfexp
	fmulp	ST(1), ST(0)
	fmul	QWORD PTR $T98456
	call	__FIsqrt
	fmulp	ST(1), ST(0)
	fld	QWORD PTR _BLURBEGIN
	fsub	QWORD PTR _BLUREND
	fxch	ST(1)
	fmul	QWORD PTR _XC
	fdivrp	ST(1), ST(0)
	fadd	QWORD PTR _XRZ1
	fstp	QWORD PTR _XRZ1
$L375:

; 410  :  	-		(BLURBEGIN-TOF1)/(BLURBEGIN-BLUREND)
; 411  : 		END IF	
; 412  : C
; 413  : 		CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	call	_RANDOM@4

; 414  : 	TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1

; 415  : 		XC=1

	fld	QWORD PTR _TEMP1
	fcomp	QWORD PTR $T98454
	fwait
	fnstsw	ax
	test	ah, 1
	mov	DWORD PTR _XC, edi
	mov	DWORD PTR _XC+4, ebx

; 416  : 		IF (TEMP1.GE..5) XC=-1.

	jne	SHORT $L377
	mov	DWORD PTR _XC, edi
	mov	DWORD PTR _XC+4, ebp
$L377:

; 417  : 	CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	call	_RANDOM@4

; 418  : 		TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1

; 419  : 		IF (TOF2 .LE. BLURBEGIN) THEN

	fld	QWORD PTR _BLURBEGIN
	fcomp	QWORD PTR _TOF2
	fwait
	fnstsw	ax
	test	ah, 1
	jne	$L351

; 420  : 		XRZ2=XRZ2+XC*DSQRT(-2.*BL0**2 * DLOG(TEMP1))* 

	fld	QWORD PTR _BLURBEGIN
	fsub	QWORD PTR _TOF2
	fld	QWORD PTR _TEMP1
	call	__FIlog
	fld	QWORD PTR _BL0
	jmp	$L98461
$L352:

; 421  :   	-		(BLURBEGIN-TOF2)/(BLURBEGIN-BLUREND)
; 422  : 		END IF
; 423  : C
; 424  : 	END IF
; 425  : C
; 426  : C
; 427  : 	ELSE

	call	_RANDOM@4

; 428  : 		CALL RANDOM(TEMP)
; 429  : 		TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1
$L381:

; 430  :  2121		CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	call	_RANDOM@4

; 431  : 		TEMP2=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP2

; 432  : 	IF (TEMP2.LT..5) GOTO 2121

	fld	QWORD PTR _TEMP2
	fcomp	QWORD PTR $T98454
	fwait
	fnstsw	ax
	test	ah, 1
	jne	SHORT $L381

; 433  : 	TOF1=(FMAX-FMIN)*TEMP1+FMIN

	fld	QWORD PTR _FMAX
	fsub	QWORD PTR _FMIN
	fld	QWORD PTR $T98443
	fxch	ST(1)
	fmul	QWORD PTR _TEMP1
	fadd	QWORD PTR _FMIN
	fstp	QWORD PTR _TOF1

; 434  : 	XAA1=TEMP2*10**TOF1

	fld	QWORD PTR _TOF1
	call	__FIfexp
	fmul	QWORD PTR _TEMP2
	push	OFFSET FLAT:_TEMP
	fstp	QWORD PTR _XAA1

; 435  : C
; 436  : 		CALL RANDOM(TEMP)

	call	_RANDOM@4

; 437  : 		TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1
$L384:

; 438  :  3131		CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	call	_RANDOM@4

; 439  : 		TEMP2=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP2

; 440  : 	IF (TEMP2.LT..5) GOTO 3131

	fld	QWORD PTR _TEMP2
	fcomp	QWORD PTR $T98454
	fwait
	fnstsw	ax
	test	ah, 1
	jne	SHORT $L384
	fld	QWORD PTR _FMAX
	fsub	QWORD PTR _FMIN
	fstp	QWORD PTR -24+[esp+96]

; 441  : 	TOF2=(FMAX-FMIN)*TEMP1+FMIN

	fld	QWORD PTR -24+[esp+96]
	fmul	QWORD PTR _TEMP1
	fadd	QWORD PTR _FMIN
	fstp	QWORD PTR _TOF2

; 442  : 	XAA2=TEMP2*10**TOF2

	fld	QWORD PTR $T98443
	fld	QWORD PTR _TOF2
	call	__FIfexp
	fmul	QWORD PTR _TEMP2
	fstp	QWORD PTR _XAA2
	fld	QWORD PTR _FMIN
	fabs
	fstp	QWORD PTR -8+[esp+96]
	fild	DWORD PTR _KA
	fstp	QWORD PTR -16+[esp+96]

; 443  : C
; 444  : C	XXX1=(TOF1+DABS(FMIN)+1.+B-(DABS(FMIN)+1.))/A
; 445  : C	XXX2=(TOF2+DABS(FMIN)+1.+B-(DABS(FMIN)+1.))/A
; 446  : C	XRZ1=(-SIGMA *DLOG(XXX1))+ZNUC
; 447  : C	XRZ2=(-SIGMA *DLOG(XXX2))+ZNUC
; 448  : 	XRZ1=-DLOG((TOF1+DABS(FMIN))/(FMAX-FMIN))/KA - ZNUC

	fld	QWORD PTR -8+[esp+96]
	fadd	QWORD PTR _TOF1
	fdiv	QWORD PTR -24+[esp+96]
	call	__FIlog
	fdiv	QWORD PTR -16+[esp+96]
	fadd	QWORD PTR _ZNUC
	fchs
	fstp	QWORD PTR _XRZ1

; 449  : 	XRZ2=-DLOG((TOF2+DABS(FMIN))/(FMAX-FMIN))/KA - ZNUC

	fld	QWORD PTR _TOF2
	fadd	QWORD PTR -8+[esp+96]
	fdiv	QWORD PTR -24+[esp+96]
	call	__FIlog
	fdiv	QWORD PTR -16+[esp+96]
	push	OFFSET FLAT:_TEMP
	fadd	QWORD PTR _ZNUC
	fchs
	fstp	QWORD PTR _XRZ2

; 450  : 	CALL RANDOM(TEMP)

	call	_RANDOM@4

; 451  : 	TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1

; 452  : 	IF (TEMP1.GE..5) XRZ1=-1.*XRZ1

	fld	QWORD PTR _TEMP1
	fcomp	QWORD PTR $T98454
	fwait
	fnstsw	ax
	test	ah, 1
	jne	SHORT $L387
	fld	QWORD PTR _XRZ1
	fchs
	fstp	QWORD PTR _XRZ1
$L387:

; 453  : 	CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	call	_RANDOM@4

; 454  : 	TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1

; 455  : 	IF (TEMP1.GE..5) XRZ2=-1*XRZ2

	fld	QWORD PTR _TEMP1
	fcomp	QWORD PTR $T98454
	fwait
	fnstsw	ax
	test	ah, 1
	jne	SHORT $L389
	fld	QWORD PTR _XRZ2
	fchs
	fstp	QWORD PTR _XRZ2
$L389:

; 456  : C
; 457  : 	IF (BL1.GT.0.) THEN

	fld	QWORD PTR _BL1
	fcomp	QWORD PTR $T98428
	fwait
	fnstsw	ax
	test	ah, 65					; 00000041H
	jne	$L351

; 458  : 		CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	call	_RANDOM@4

; 459  : 	TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1

; 460  : 		XC=1

	fld	QWORD PTR _TEMP1
	fcomp	QWORD PTR $T98454
	fwait
	fnstsw	ax
	test	ah, 1
	mov	DWORD PTR _XC, edi
	mov	DWORD PTR _XC+4, ebx

; 461  : 		IF (TEMP1.GE..5) XC=-1.

	jne	SHORT $L393
	mov	DWORD PTR _XC, edi
	mov	DWORD PTR _XC+4, ebp
$L393:

; 462  : 	CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	call	_RANDOM@4

; 463  : 		TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1

; 464  : 		IF (TOF1 .LE. BLURBEGIN) THEN

	fld	QWORD PTR _BLURBEGIN
	fcomp	QWORD PTR _TOF1
	fwait
	fnstsw	ax
	test	ah, 1
	jne	SHORT $L395

; 465  : 		XRZ1=XRZ1+XC*DSQRT(-2.*BL1**2 * DLOG(TEMP1))* 

	fld	QWORD PTR _BLURBEGIN
	fsub	QWORD PTR _TOF1
	fld	QWORD PTR _TEMP1
	call	__FIlog
	fld	QWORD PTR _BL1
	fld	QWORD PTR $T98441
	call	__FIfexp
	fmulp	ST(1), ST(0)
	fmul	QWORD PTR $T98456
	call	__FIsqrt
	fmulp	ST(1), ST(0)
	fld	QWORD PTR _BLURBEGIN
	fsub	QWORD PTR _BLUREND
	fxch	ST(1)
	fmul	QWORD PTR _XC
	fdivrp	ST(1), ST(0)
	fadd	QWORD PTR _XRZ1
	fstp	QWORD PTR _XRZ1
$L395:

; 466  :   	-		(BLURBEGIN-TOF1)/(BLURBEGIN-BLUREND)
; 467  : 		END IF
; 468  : 	
; 469  : C
; 470  : 		CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	call	_RANDOM@4

; 471  : 	TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1

; 472  : 		XC=1

	fld	QWORD PTR _TEMP1
	fcomp	QWORD PTR $T98454
	fwait
	fnstsw	ax
	test	ah, 1
	mov	DWORD PTR _XC, edi
	mov	DWORD PTR _XC+4, ebx

; 473  : 		IF (TEMP1.GE..5) XC=-1.

	jne	SHORT $L397
	mov	DWORD PTR _XC, edi
	mov	DWORD PTR _XC+4, ebp
$L397:

; 474  : 	CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	call	_RANDOM@4

; 475  : 		TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1

; 476  : 		IF (TOF2 .LE. BLURBEGIN) THEN

	fld	QWORD PTR _BLURBEGIN
	fcomp	QWORD PTR _TOF2
	fwait
	fnstsw	ax
	test	ah, 1
	jne	SHORT $L351

; 477  : 		XRZ2=XRZ2+XC*DSQRT(-2.*BL1**2 * DLOG(TEMP1))* 

	fld	QWORD PTR _BLURBEGIN
	fsub	QWORD PTR _TOF2
	fld	QWORD PTR _TEMP1
	call	__FIlog
	fld	QWORD PTR _BL1
$L98461:
	fld	QWORD PTR $T98441
	call	__FIfexp
	fmulp	ST(1), ST(0)
	fmul	QWORD PTR $T98456
	call	__FIsqrt
	fmulp	ST(1), ST(0)
	fmul	QWORD PTR _XC
	fld	QWORD PTR _BLURBEGIN
	fsub	QWORD PTR _BLUREND
	fdivp	ST(1), ST(0)
	fadd	QWORD PTR _XRZ2
	fstp	QWORD PTR _XRZ2
$L351:

; 478  :   	-		(BLURBEGIN-TOF2)/(BLURBEGIN-BLUREND)
; 479  : 
; 480  : 		END IF
; 481  : C
; 482  : 	END IF
; 483  : 
; 484  : 
; 485  : C
; 486  : 	END IF
; 487  : 	CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	call	_RANDOM@4

; 488  : 	TEMP1=TEMP

	fld	DWORD PTR _TEMP
	push	OFFSET FLAT:_TEMP
	fstp	QWORD PTR _TEMP1

; 489  : 	CALL RANDOM(TEMP)

	call	_RANDOM@4

; 490  : 	TEMP2=TEMP

	fld	DWORD PTR _TEMP
	fld1
	fadd	QWORD PTR _BETAMIN
	fxch	ST(1)
	fstp	QWORD PTR _TEMP2
	fstp	QWORD PTR -8+[esp+96]

; 491  : 	XB12=DABS(TEMP1*10.**(((BETAMAX+1.)-(BETAMIN+1))*

	fld	QWORD PTR $T98443
	fld1
	fadd	QWORD PTR _BETAMAX
	fsub	QWORD PTR -8+[esp+96]
	fmul	QWORD PTR _TEMP2
	fadd	QWORD PTR -8+[esp+96]
	call	__FIfexp
	fmul	QWORD PTR _TEMP1
	push	OFFSET FLAT:_XB12
	push	OFFSET FLAT:_XRZ2
	push	OFFSET FLAT:_XRY2
	push	OFFSET FLAT:_XRX2
	fabs
	push	OFFSET FLAT:_XAA2
	push	OFFSET FLAT:_XRZ1
	fstp	QWORD PTR _XB12

; 492  : 	-TEMP2+(BETAMIN+1)))
; 493  : C	WRITE(6,*)'BETA=',XB12
; 494  : C	WRITE(6,292)I,XAA1,XAA2,XRZ1,XRZ2
; 495  : 	CALL HAMILTONp(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,

	push	OFFSET FLAT:_XRY1
	mov	eax, DWORD PTR _H
	push	OFFSET FLAT:_XRX1
	mov	ecx, DWORD PTR _S
	push	OFFSET FLAT:_XAA1
	mov	edx, DWORD PTR _EVEC
	push	OFFSET FLAT:_EMIN
	push	eax
	push	ecx
	mov	eax, DWORD PTR _B12
	push	edx
	mov	ecx, DWORD PTR _RZ2
	push	eax
	mov	edx, DWORD PTR _RY2
	push	ecx
	mov	eax, DWORD PTR _RX2
	push	edx
	mov	ecx, DWORD PTR _AA2
	push	eax
	mov	edx, DWORD PTR _RZ1
	push	ecx
	mov	eax, DWORD PTR _RY1
	push	edx
	mov	ecx, DWORD PTR _RX1
	push	eax
	mov	edx, DWORD PTR _AA1
	push	ecx
	push	edx
	push	OFFSET FLAT:_M
	call	_HAMILTONP@92

; 496  :      1                    EVEC,S,H,EMIN,XAA1,XRX1,XRY1,XRZ1,
; 497  :      2					XAA2,XRX2,XRY2,XRZ2,XB12)
; 498  : 	IF (I .EQ. 1) ENER=EMIN

	mov	ecx, DWORD PTR _I$[esp+96]
	cmp	ecx, 1
	jne	SHORT $L403
	mov	eax, DWORD PTR _EMIN+4
	mov	ecx, DWORD PTR _EMIN
	mov	DWORD PTR _ENER+4, eax
	mov	DWORD PTR _ENER, ecx
$L403:

; 499  : 	IF (EMIN.LT.ENER) THEN 

	fld	QWORD PTR _ENER
	fcomp	QWORD PTR _EMIN
	fwait
	fnstsw	ax
	test	ah, 65					; 00000041H
	jne	$L350

; 500  : 		ENER=EMIN

	mov	eax, DWORD PTR _EMIN+4
	mov	ecx, DWORD PTR _EMIN
	mov	edx, DWORD PTR _XAA1+4
	mov	DWORD PTR _ENER+4, eax
	mov	eax, DWORD PTR _XAA1
	mov	DWORD PTR _ENER, ecx

; 501  : 		TA1=XAA1

	mov	ecx, DWORD PTR _XAA2+4
	mov	DWORD PTR _TA1+4, edx
	mov	edx, DWORD PTR _XAA2
	mov	DWORD PTR _TA1, eax

; 502  : 		TA2=XAA2

	mov	eax, DWORD PTR _XB12+4
	mov	DWORD PTR _TA2+4, ecx
	mov	ecx, DWORD PTR _XB12
	mov	DWORD PTR _TA2, edx

; 503  : 		TB12=XB12

	mov	edx, DWORD PTR _XRZ1+4
	mov	DWORD PTR _TB12+4, eax
	mov	eax, DWORD PTR _XRZ1
	mov	DWORD PTR _TB12, ecx

; 504  : 		TRZ1=XRZ1

	mov	ecx, DWORD PTR _XRZ2+4
	mov	DWORD PTR _TRZ1+4, edx
	mov	edx, DWORD PTR _XRZ2
	mov	DWORD PTR _TRZ1, eax

; 505  : 		TRZ2=XRZ2

	mov	DWORD PTR _TRZ2+4, ecx
	mov	DWORD PTR _TRZ2, edx
$L350:

; 370  : 	DO I=1,IVIEW

	dec	esi
	mov	eax, DWORD PTR _I$[esp+96]
	inc	eax
	test	esi, esi
	mov	DWORD PTR _I$[esp+96], eax
	jne	$L98316
$L349:

; 506  : 	END IF
; 507  :       END DO
; 508  : C	STOP
; 509  : C
; 510  :  292	FORMAT(I4,4D24.8)
; 511  : 	CALL BIOUT(' ')

	push	1
	push	OFFSET FLAT:_.SG412_
	call	_BIOUT@8

; 512  : 	WRITE(6,*)'BEST DIAGNOSTIC VALUE IN ',IVIEW,' GUESSES=',ENER

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	25					; 00000019H
	push	OFFSET FLAT:_.SG413_
	call	_lwrite_ch
	add	esp, 8
	mov	eax, DWORD PTR _IVIEW
	push	eax
	call	_lwrite_i
	add	esp, 4
	push	9
	push	OFFSET FLAT:_.SG414_
	call	_lwrite_ch
	add	esp, 8
	mov	eax, DWORD PTR _ENER+4
	mov	ecx, DWORD PTR _ENER
	npad	1
	push	eax
	push	ecx
	call	_lwrite_d
	add	esp, 8
	call	_f90_endwrite

; 513  : C	WRITE(6,*)TA1,TRZ1,TA2,TRZ2,TB12
; 514  : 	WRITE(6,*)'***************************'

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	27					; 0000001bH
	push	OFFSET FLAT:_.SG415_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 515  : 	WRITE(6,*)'ADDING THIS FUNCTION TO Psi'

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	27					; 0000001bH
	push	OFFSET FLAT:_.SG416_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 516  : 	WRITE(6,*)'***************************'

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	27					; 0000001bH
	push	OFFSET FLAT:_.SG417_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 517  : 	ALLOCATE(EVECSTORET(M),TSTORE1T(M),

	mov	DWORD PTR _.T418_$[esp+96], 0
	mov	ecx, DWORD PTR _M
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7TSTORE1T
	call	_make_info
	lea	ecx, DWORD PTR _.T418_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7TSTORE1T
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _TSTORE1T, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7TSTORE2T
	call	_make_info
	lea	ecx, DWORD PTR _.T418_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7TSTORE2T
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _TSTORE2T, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7TSTORE3T
	call	_make_info
	lea	ecx, DWORD PTR _.T418_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7TSTORE3T
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _TSTORE3T, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7RZ1STORET
	call	_make_info
	lea	ecx, DWORD PTR _.T418_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7RZ1STORET
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _RZ1STORET, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7RZ2STORET
	call	_make_info
	lea	ecx, DWORD PTR _.T418_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7RZ2STORET
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _RZ2STORET, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7EVECSTORET
	call	_make_info
	lea	ecx, DWORD PTR _.T418_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7EVECSTORET
	call	__Allocate_i
	mov	ecx, DWORD PTR _.T418_$[esp+108]
	add	esp, 12					; 0000000cH
	mov	edx, DWORD PTR _M
	mov	DWORD PTR _EVECSTORET, eax
	mov	DWORD PTR _ERROR, ecx
	test	edx, edx

; 518  :      -			TSTORE2T(M),TSTORE3T(M),
; 519  :      -			RZ1STORET(M),RZ2STORET(M),
; 520  :      -			STAT=ERROR)
; 521  : 	DO III=1,M

	mov	ecx, 1
	jle	$L426
$L98318:

; 522  : 		EVECSTORET(III)=EVEC(III)

	mov	eax, DWORD PTR _7EVEC+8
	mov	ebp, DWORD PTR _EVEC
	add	eax, ecx
	mov	edi, DWORD PTR _7EVECSTORET+8
	add	edi, ecx
	mov	ebx, DWORD PTR _EVECSTORET
	mov	esi, DWORD PTR [ebp+eax*8+4]
	dec	edx
	mov	eax, DWORD PTR [ebp+eax*8]
	mov	DWORD PTR [ebx+edi*8+4], esi
	mov	ebp, DWORD PTR _AA1
	mov	DWORD PTR [ebx+edi*8], eax

; 523  : 		TSTORE1T(III)=AA1(III)

	mov	eax, DWORD PTR _7AA1+8
	mov	ebx, DWORD PTR _7TSTORE1T+8
	add	eax, ecx
	add	ebx, ecx
	mov	edi, DWORD PTR _TSTORE1T
	npad	1
	mov	esi, DWORD PTR [ebp+eax*8+4]
	mov	eax, DWORD PTR [ebp+eax*8]
	mov	DWORD PTR [edi+ebx*8+4], esi
	mov	DWORD PTR [edi+ebx*8], eax

; 524  : 		TSTORE2T(III)=AA2(III)

	mov	eax, DWORD PTR _AA2
	mov	ebp, DWORD PTR _7AA2+8
	mov	ebx, DWORD PTR _7TSTORE2T+8
	add	ebp, ecx
	add	ebx, ecx
	mov	edi, DWORD PTR _TSTORE2T
	npad	1
	mov	esi, DWORD PTR [eax+ebp*8+4]
	mov	eax, DWORD PTR [eax+ebp*8]
	mov	DWORD PTR [edi+ebx*8+4], esi
	mov	esi, DWORD PTR _B12
	mov	DWORD PTR [edi+ebx*8], eax

; 525  : 		TSTORE3T(III)=B12(III)

	mov	eax, DWORD PTR _7B12+8
	add	eax, ecx
	mov	ebx, DWORD PTR _7TSTORE3T+8
	add	ebx, ecx
	mov	edi, DWORD PTR _TSTORE3T
	mov	ebp, DWORD PTR [esi+eax*8+4]
	mov	eax, DWORD PTR [esi+eax*8]
	mov	DWORD PTR [edi+ebx*8+4], ebp
	mov	DWORD PTR [edi+ebx*8], eax

; 526  : 		RZ1STORET(III)=RZ1(III)

	mov	eax, DWORD PTR _RZ1
	mov	esi, DWORD PTR _7RZ1+8
	mov	ebx, DWORD PTR _7RZ1STORET+8
	add	esi, ecx
	add	ebx, ecx
	mov	edi, DWORD PTR _RZ1STORET
	npad	1
	mov	ebp, DWORD PTR [eax+esi*8+4]
	mov	eax, DWORD PTR [eax+esi*8]
	mov	DWORD PTR [edi+ebx*8+4], ebp
	mov	DWORD PTR [edi+ebx*8], eax

; 527  : 		RZ2STORET(III)=RZ2(III)

	mov	eax, DWORD PTR _RZ2
	mov	esi, DWORD PTR _7RZ2+8
	mov	ebx, DWORD PTR _7RZ2STORET+8
	add	esi, ecx
	add	ebx, ecx
	mov	edi, DWORD PTR _RZ2STORET
	inc	ecx
	mov	ebp, DWORD PTR [eax+esi*8+4]
	test	edx, edx
	mov	eax, DWORD PTR [eax+esi*8]
	mov	DWORD PTR [edi+ebx*8+4], ebp
	mov	DWORD PTR [edi+ebx*8], eax

; 518  :      -			TSTORE2T(M),TSTORE3T(M),
; 519  :      -			RZ1STORET(M),RZ2STORET(M),
; 520  :      -			STAT=ERROR)
; 521  : 	DO III=1,M

	jne	$L98318
$L426:

; 528  : 	END DO
; 529  : 	M=M+1

	mov	eax, DWORD PTR _M
	push	0
	inc	eax
	push	8
	push	OFFSET FLAT:_7RX1
	mov	ecx, DWORD PTR _RX1
	push	ecx
	mov	DWORD PTR _M, eax

; 530  : 	DEALLOCATE (AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,EVEC,S,H,CO,

	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _RY1
	mov	DWORD PTR _RX1, 0
	push	0
	push	8
	push	OFFSET FLAT:_7RY1
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _RZ1
	mov	DWORD PTR _RY1, 0
	push	0
	push	8
	push	OFFSET FLAT:_7RZ1
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _AA2
	mov	DWORD PTR _RZ1, 0
	push	0
	push	8
	push	OFFSET FLAT:_7AA2
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _RX2
	mov	DWORD PTR _AA2, 0
	push	0
	push	8
	push	OFFSET FLAT:_7RX2
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _RY2
	mov	DWORD PTR _RX2, 0
	push	0
	push	8
	push	OFFSET FLAT:_7RY2
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _RZ2
	mov	DWORD PTR _RY2, 0
	push	0
	push	8
	push	OFFSET FLAT:_7RZ2
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _B12
	mov	DWORD PTR _RZ2, 0
	push	0
	push	8
	push	OFFSET FLAT:_7B12
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _EVEC
	mov	DWORD PTR _B12, 0
	push	0
	push	8
	push	OFFSET FLAT:_7EVEC
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _S
	mov	DWORD PTR _EVEC, 0
	push	0
	push	8
	push	OFFSET FLAT:_7S
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _H
	mov	DWORD PTR _S, 0
	push	0
	push	8
	push	OFFSET FLAT:_7H
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _CO
	mov	DWORD PTR _H, 0
	push	0
	push	8
	push	OFFSET FLAT:_7CO
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	DWORD PTR _CO, 0
	push	0
	push	8
	push	OFFSET FLAT:_7GRAD
	mov	ecx, DWORD PTR _GRAD
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _T1
	mov	DWORD PTR _GRAD, 0
	push	0
	push	8
	push	OFFSET FLAT:_7T1
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _T2
	mov	DWORD PTR _T1, 0
	push	0
	push	8
	push	OFFSET FLAT:_7T2
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _T3
	mov	DWORD PTR _T2, 0
	push	0
	push	8
	push	OFFSET FLAT:_7T3
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _AA1
	mov	DWORD PTR _T3, 0
	push	0
	push	8
	push	OFFSET FLAT:_7AA1
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _AA1, 0

; 531  : 	-			GRAD,T1,T2,T3)
; 532  :       ALLOCATE(	AA1(M),RX1(M),RY1(M),RZ1(M),

	mov	DWORD PTR _.T428_$[esp+96], 0
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7RX1
	call	_make_info
	lea	ecx, DWORD PTR _.T428_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7RX1
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _RX1, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7RY1
	call	_make_info
	lea	ecx, DWORD PTR _.T428_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7RY1
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _RY1, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7RZ1
	call	_make_info
	lea	ecx, DWORD PTR _.T428_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7RZ1
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _RZ1, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7AA2
	call	_make_info
	lea	ecx, DWORD PTR _.T428_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7AA2
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _AA2, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7RX2
	call	_make_info
	lea	ecx, DWORD PTR _.T428_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7RX2
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _RX2, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7RY2
	call	_make_info
	lea	ecx, DWORD PTR _.T428_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7RY2
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _RY2, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7RZ2
	call	_make_info
	lea	ecx, DWORD PTR _.T428_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7RZ2
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _RZ2, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7B12
	call	_make_info
	lea	ecx, DWORD PTR _.T428_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7B12
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _B12, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7EVEC
	call	_make_info
	lea	ecx, DWORD PTR _.T428_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7EVEC
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _EVEC, eax
	npad	1
	push	ecx
	push	1
	push	ecx
	push	1
	push	2
	push	OFFSET FLAT:_7S
	call	_make_info
	lea	ecx, DWORD PTR _.T428_$[esp+120]
	add	esp, 24					; 00000018H
	push	ecx
	push	8
	push	OFFSET FLAT:_7S
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _S, eax
	npad	1
	push	ecx
	push	1
	push	ecx
	push	1
	push	2
	push	OFFSET FLAT:_7H
	call	_make_info
	lea	ecx, DWORD PTR _.T428_$[esp+120]
	add	esp, 24					; 00000018H
	push	ecx
	push	8
	push	OFFSET FLAT:_7H
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _H, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7CO
	call	_make_info
	lea	ecx, DWORD PTR _.T428_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7CO
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _CO, eax
	lea	eax, DWORD PTR [ecx+ecx*2]
	add	eax, eax
	push	eax
	push	1
	push	1
	push	OFFSET FLAT:_7GRAD
	call	_make_info
	lea	ecx, DWORD PTR _.T428_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7GRAD
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _GRAD, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7T1
	call	_make_info
	lea	ecx, DWORD PTR _.T428_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7T1
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _T1, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7T2
	call	_make_info
	lea	ecx, DWORD PTR _.T428_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7T2
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _T2, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7T3
	call	_make_info
	lea	ecx, DWORD PTR _.T428_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7T3
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _T3, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7AA1
	call	_make_info
	lea	ecx, DWORD PTR _.T428_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7AA1
	call	__Allocate_i
	mov	ecx, DWORD PTR _.T428_$[esp+108]
	add	esp, 12					; 0000000cH
	mov	edx, 1
	mov	DWORD PTR _AA1, eax
	mov	DWORD PTR _ERROR, ecx

; 533  : 	-			AA2(M),RX2(M),RY2(M),RZ2(M),
; 534  : 	-			B12(M),EVEC(M),S(M,M),H(M,M),
; 535  :      -			CO(M),GRAD(M*6),T1(M),T2(M),
; 536  :      -			T3(M),STAT=ERROR)
; 537  : 
; 538  : 	DO III=1,M-1

	mov	ecx, DWORD PTR _M
	dec	ecx
	test	ecx, ecx
	jle	$L429
$L98319:

; 539  : 		EVEC(III)=EVECSTORET(III)

	mov	esi, DWORD PTR _7EVECSTORET+8
	mov	ebp, DWORD PTR _EVECSTORET
	add	esi, edx
	mov	edi, DWORD PTR _7EVEC+8
	add	edi, edx
	mov	ebx, DWORD PTR _EVEC
	mov	eax, DWORD PTR [ebp+esi*8+4]
	dec	ecx
	mov	ebp, DWORD PTR [ebp+esi*8]
	mov	DWORD PTR [ebx+edi*8+4], eax
	mov	esi, DWORD PTR _TSTORE1T
	mov	DWORD PTR [ebx+edi*8], ebp

; 540  : 	AA1(III)=TSTORE1T(III)

	mov	ebp, DWORD PTR _7TSTORE1T+8
	mov	edi, DWORD PTR _7AA1+8
	add	ebp, edx
	add	edi, edx
	mov	ebx, DWORD PTR _AA1
	npad	1
	mov	eax, DWORD PTR [esi+ebp*8+4]
	mov	ebp, DWORD PTR [esi+ebp*8]
	mov	DWORD PTR [ebx+edi*8+4], eax
	mov	DWORD PTR [ebx+edi*8], ebp

; 541  : 	AA2(III)=TSTORE2T(III)

	mov	ebp, DWORD PTR _TSTORE2T
	mov	esi, DWORD PTR _7TSTORE2T+8
	mov	edi, DWORD PTR _7AA2+8
	add	esi, edx
	add	edi, edx
	mov	ebx, DWORD PTR _AA2
	npad	1
	mov	eax, DWORD PTR [ebp+esi*8+4]
	mov	ebp, DWORD PTR [ebp+esi*8]
	mov	DWORD PTR [ebx+edi*8+4], eax
	mov	DWORD PTR [ebx+edi*8], ebp

; 542  : 	B12(III)=TSTORE3T(III)

	mov	ebp, DWORD PTR _TSTORE3T
	mov	esi, DWORD PTR _7TSTORE3T+8
	mov	ebx, DWORD PTR _7B12+8
	add	esi, edx
	add	ebx, edx
	mov	edi, DWORD PTR _B12
	npad	1
	mov	eax, DWORD PTR [ebp+esi*8+4]
	mov	ebp, DWORD PTR [ebp+esi*8]
	mov	DWORD PTR [edi+ebx*8+4], eax
	mov	eax, DWORD PTR _RZ1STORET
	mov	DWORD PTR [edi+ebx*8], ebp

; 543  : 	RZ1(III)=RZ1STORET(III)

	mov	ebp, DWORD PTR _7RZ1STORET+8
	add	ebp, edx
	mov	ebx, DWORD PTR _7RZ1+8
	add	ebx, edx
	mov	edi, DWORD PTR _RZ1
	mov	esi, DWORD PTR [eax+ebp*8+4]
	mov	eax, DWORD PTR [eax+ebp*8]
	mov	DWORD PTR [edi+ebx*8+4], esi
	mov	DWORD PTR [edi+ebx*8], eax

; 544  : 	RZ2(III)=RZ2STORET(III)

	mov	eax, DWORD PTR _RZ2STORET
	mov	ebp, DWORD PTR _7RZ2STORET+8
	mov	edi, DWORD PTR _7RZ2+8
	add	ebp, edx
	add	edi, edx
	mov	ebx, DWORD PTR _RZ2
	npad	1
	mov	esi, DWORD PTR [eax+ebp*8+4]
	mov	eax, DWORD PTR [eax+ebp*8]
	mov	DWORD PTR [ebx+edi*8+4], esi
	mov	DWORD PTR [ebx+edi*8], eax

; 545  : 	RX1(III)=0.

	mov	eax, DWORD PTR _RX1
	mov	esi, DWORD PTR _7RX1+8
	add	esi, edx
	mov	DWORD PTR [eax+esi*8], 0
	mov	DWORD PTR [eax+esi*8+4], 0

; 546  : 	RY1(III)=0.

	mov	esi, DWORD PTR _7RY1+8
	mov	eax, DWORD PTR _RY1
	add	esi, edx
	mov	DWORD PTR [eax+esi*8], 0
	mov	DWORD PTR [eax+esi*8+4], 0

; 547  : 	RY2(III)=0.

	mov	esi, DWORD PTR _7RY2+8
	mov	eax, DWORD PTR _RY2
	add	esi, edx
	mov	DWORD PTR [eax+esi*8], 0
	mov	DWORD PTR [eax+esi*8+4], 0

; 548  : 	RX2(III)=0.

	mov	esi, DWORD PTR _7RX2+8
	mov	eax, DWORD PTR _RX2
	add	esi, edx
	inc	edx
	test	ecx, ecx
	mov	DWORD PTR [eax+esi*8], 0
	mov	DWORD PTR [eax+esi*8+4], 0

; 533  : 	-			AA2(M),RX2(M),RY2(M),RZ2(M),
; 534  : 	-			B12(M),EVEC(M),S(M,M),H(M,M),
; 535  :      -			CO(M),GRAD(M*6),T1(M),T2(M),
; 536  :      -			T3(M),STAT=ERROR)
; 537  : 
; 538  : 	DO III=1,M-1

	jne	$L98319
$L429:

; 549  : 	END DO
; 550  : 	AA1(M)=TA1

	mov	ecx, DWORD PTR _M
	mov	eax, DWORD PTR _7AA1+8
	add	ecx, eax
	mov	edx, DWORD PTR _AA1
	mov	ebx, DWORD PTR _TA1+4
	mov	eax, DWORD PTR _TA1
	mov	DWORD PTR [edx+ecx*8+4], ebx
	mov	DWORD PTR [edx+ecx*8], eax

; 551  : 	AA2(M)=TA2

	mov	ecx, DWORD PTR _M
	mov	ebp, DWORD PTR _7AA2+8
	mov	edx, DWORD PTR _AA2
	add	ecx, ebp
	mov	eax, DWORD PTR _TA2+4
	mov	ebx, DWORD PTR _TA2
	npad	1
	mov	DWORD PTR [edx+ecx*8+4], eax
	mov	eax, DWORD PTR _TB12+4
	mov	DWORD PTR [edx+ecx*8], ebx

; 552  : 	B12(M)=TB12

	mov	edx, DWORD PTR _7B12+8
	mov	ecx, DWORD PTR _M
	mov	ebx, DWORD PTR _TB12
	add	edx, ecx
	mov	ecx, DWORD PTR _B12
	mov	DWORD PTR [ecx+edx*8+4], eax
	mov	eax, DWORD PTR _TRZ1+4
	mov	DWORD PTR [ecx+edx*8], ebx

; 553  : 	RZ1(M)=TRZ1

	mov	edx, DWORD PTR _M
	mov	esi, DWORD PTR _7RZ1+8
	mov	ecx, DWORD PTR _RZ1
	add	edx, esi
	mov	ebx, DWORD PTR _TRZ1
	mov	DWORD PTR [ecx+edx*8+4], eax
	mov	eax, DWORD PTR _TRZ2+4
	mov	DWORD PTR [ecx+edx*8], ebx

; 554  : 	RZ2(M)=TRZ2

	mov	edx, DWORD PTR _7RZ2+8
	mov	edi, DWORD PTR _M
	mov	ecx, DWORD PTR _RZ2
	add	edx, edi
	mov	ebx, DWORD PTR _TRZ2
	mov	DWORD PTR [ecx+edx*8+4], eax
	mov	eax, DWORD PTR _RX1
	mov	DWORD PTR [ecx+edx*8], ebx

; 555  : 	RX1(M)=0.

	mov	ecx, DWORD PTR _7RX1+8
	add	ecx, DWORD PTR _M
	mov	DWORD PTR [eax+ecx*8], 0
	mov	DWORD PTR [eax+ecx*8+4], 0

; 556  : 	RY1(M)=0.

	mov	eax, DWORD PTR _M
	mov	ecx, DWORD PTR _7RY1+8
	add	eax, ecx
	mov	ecx, DWORD PTR _RY1
	mov	DWORD PTR [ecx+eax*8], 0
	mov	DWORD PTR [ecx+eax*8+4], 0

; 557  : 	RY2(M)=0.

	mov	eax, DWORD PTR _7RY2+8
	mov	edi, DWORD PTR _M
	add	eax, edi
	mov	ecx, DWORD PTR _RY2
	mov	DWORD PTR [ecx+eax*8], 0
	mov	DWORD PTR [ecx+eax*8+4], 0

; 558  : 	RX2(M)=0.

	mov	eax, DWORD PTR _7RX2+8
	mov	edi, DWORD PTR _M
	add	eax, edi
	mov	ecx, DWORD PTR _RX2
	mov	DWORD PTR [ecx+eax*8], 0
	mov	edx, DWORD PTR _TSTORE1T
	mov	DWORD PTR [ecx+eax*8+4], 0

; 559  : 
; 560  : 
; 561  : 	DEALLOCATE(EVECSTORET,TSTORE1T,TSTORE2T,TSTORE3T,

	mov	DWORD PTR _.T431_$[esp+96], 0
	lea	eax, DWORD PTR _.T431_$[esp+96]
	push	eax
	push	8
	push	OFFSET FLAT:_7TSTORE1T
	push	edx
	call	__Deallocate_i
	lea	ecx, DWORD PTR _.T431_$[esp+112]
	add	esp, 16					; 00000010H
	mov	DWORD PTR _TSTORE1T, 0
	push	ecx
	mov	edx, DWORD PTR _TSTORE2T
	push	8
	push	OFFSET FLAT:_7TSTORE2T
	push	edx
	call	__Deallocate_i
	lea	ecx, DWORD PTR _.T431_$[esp+112]
	add	esp, 16					; 00000010H
	mov	DWORD PTR _TSTORE2T, 0
	push	ecx
	mov	edx, DWORD PTR _TSTORE3T
	push	8
	push	OFFSET FLAT:_7TSTORE3T
	push	edx
	call	__Deallocate_i
	lea	ecx, DWORD PTR _.T431_$[esp+112]
	add	esp, 16					; 00000010H
	mov	DWORD PTR _TSTORE3T, 0
	push	ecx
	mov	edx, DWORD PTR _RZ1STORET
	push	8
	push	OFFSET FLAT:_7RZ1STORET
	push	edx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	DWORD PTR _RZ1STORET, 0
	lea	ecx, DWORD PTR _.T431_$[esp+96]
	mov	edx, DWORD PTR _RZ2STORET
	push	ecx
	push	8
	push	OFFSET FLAT:_7RZ2STORET
	push	edx
	call	__Deallocate_i
	lea	ecx, DWORD PTR _.T431_$[esp+112]
	add	esp, 16					; 00000010H
	mov	DWORD PTR _RZ2STORET, 0
	push	ecx
	mov	edx, DWORD PTR _EVECSTORET
	push	8
	push	OFFSET FLAT:_7EVECSTORET
	push	edx
	call	__Deallocate_i
	mov	ecx, DWORD PTR _.T431_$[esp+112]
	add	esp, 16					; 00000010H
	mov	edx, DWORD PTR _T3
	mov	DWORD PTR _ERROR, ecx
	push	OFFSET FLAT:_DFLAG
	mov	eax, DWORD PTR _T2
	push	edx
	mov	ecx, DWORD PTR _T1
	push	eax
	mov	edx, DWORD PTR _GRAD
	push	ecx
	mov	eax, DWORD PTR _CO
	push	edx
	mov	ecx, DWORD PTR _H
	push	OFFSET FLAT:_EMIN
	mov	edx, DWORD PTR _S
	mov	DWORD PTR _EVECSTORET, 0

; 562  :      -			RZ1STORET,RZ2STORET,STAT=ERROR)
; 563  : C
; 564  : 	DFLAG=0

	mov	DWORD PTR _DFLAG, 0

; 565  : 	CALL HAMILTON(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,EVEC,S,H,

	push	eax
	push	ecx
	push	edx
	mov	eax, DWORD PTR _EVEC
	push	eax
	mov	ecx, DWORD PTR _B12
	push	ecx
	mov	edx, DWORD PTR _RZ2
	push	edx
	mov	eax, DWORD PTR _RY2
	push	eax
	mov	ecx, DWORD PTR _RX2
	push	ecx
	mov	edx, DWORD PTR _AA2
	push	edx
	mov	eax, DWORD PTR _RZ1
	push	eax
	mov	ecx, DWORD PTR _RY1
	push	ecx
	mov	edx, DWORD PTR _RX1
	push	edx
	mov	eax, DWORD PTR _AA1
	push	eax
	push	OFFSET FLAT:_M
	call	_HAMILTON@80

; 566  : 	-				CO,EMIN,GRAD,T1,T2,T3,DFLAG)
; 567  : 	DO NNN=1,M

	mov	esi, DWORD PTR _M
	mov	edx, 1
	test	esi, esi
	jle	SHORT $L432
$L98320:

; 568  : 		EVEC(NNN)=CO(NNN)

	mov	ecx, DWORD PTR _7CO+8
	mov	eax, DWORD PTR _CO
	add	ecx, edx
	mov	ebx, DWORD PTR _7EVEC+8
	add	ebx, edx
	mov	edi, DWORD PTR _EVEC
	mov	ebp, DWORD PTR [eax+ecx*8+4]
	dec	esi
	mov	eax, DWORD PTR [eax+ecx*8]
	inc	edx
	mov	DWORD PTR [edi+ebx*8+4], ebp
	test	esi, esi
	mov	DWORD PTR [edi+ebx*8], eax

; 566  : 	-				CO,EMIN,GRAD,T1,T2,T3,DFLAG)
; 567  : 	DO NNN=1,M

	jne	SHORT $L98320
$L432:

; 569  : 	END DO
; 570  : 	WRITE(6,*)'NOW WE HAVE',M,' FUNCTIONS...'

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	11					; 0000000bH
	push	OFFSET FLAT:_.SG434_
	call	_lwrite_ch
	add	esp, 8
	mov	eax, DWORD PTR _M
	push	eax
	call	_lwrite_i
	add	esp, 4
	push	13					; 0000000dH
	push	OFFSET FLAT:_.SG435_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 571  : 	WRITE(6,*)'EXACT ENERGY=',EMIN

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	13					; 0000000dH
	push	OFFSET FLAT:_.SG436_
	call	_lwrite_ch
	add	esp, 8
	mov	eax, DWORD PTR _EMIN+4
	mov	ecx, DWORD PTR _EMIN
	npad	1
	push	eax
	push	ecx
	call	_lwrite_d
	add	esp, 8
	call	_f90_endwrite

; 572  : C	IDOIT=0
; 573  : C	IF (IDOIT.EQ.1)	THEN
; 574  : 	WRITE(6,*)'********************************'

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	32					; 00000020H
	push	OFFSET FLAT:_.SG437_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 575  : 	WRITE(6,*)'PERFORMING GRADIENT OPTIMIZATION'

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	32					; 00000020H
	push	OFFSET FLAT:_.SG438_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 576  : 	WRITE(6,*)'********************************'

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	32					; 00000020H
	push	OFFSET FLAT:_.SG439_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 577  : 	DO NNN=1,M

	mov	esi, DWORD PTR _M
	mov	edx, 1
	test	esi, esi
	jle	SHORT $L440
$L98321:

; 578  : 		EVEC(NNN)=CO(NNN)

	mov	eax, DWORD PTR _7CO+8
	mov	ecx, DWORD PTR _CO
	add	eax, edx
	mov	edi, DWORD PTR _7EVEC+8
	add	edi, edx
	mov	ebx, DWORD PTR _EVEC
	mov	ebp, DWORD PTR [ecx+eax*8+4]
	dec	esi
	mov	eax, DWORD PTR [ecx+eax*8]
	inc	edx
	mov	DWORD PTR [ebx+edi*8+4], ebp
	test	esi, esi
	mov	DWORD PTR [ebx+edi*8], eax

; 577  : 	DO NNN=1,M

	jne	SHORT $L98321
$L440:

; 579  : 	END DO
; 580  : 	DO NNN=1,M

	mov	edi, DWORD PTR _M
	mov	esi, 1
	test	edi, edi
	jle	$L442
	jmp	SHORT $L98322
$L98467:
	fstp	QWORD PTR [eax+ecx*8]
$L98322:

; 581  : 		T1(NNN)=DSQRT(AA1(NNN)+B12(NNN))

	mov	ecx, DWORD PTR _7B12+8
	mov	eax, DWORD PTR _B12
	add	ecx, esi
	dec	edi
	fld	QWORD PTR [eax+ecx*8]
	mov	ecx, DWORD PTR _7AA1+8
	mov	eax, DWORD PTR _AA1
	add	ecx, esi
	fadd	QWORD PTR [eax+ecx*8]
	call	__FIsqrt
	mov	ecx, DWORD PTR _7T1+8
	mov	eax, DWORD PTR _T1
	add	ecx, esi
	mov	edx, DWORD PTR _B12
	fstp	QWORD PTR [eax+ecx*8]
	mov	ebx, DWORD PTR _7B12+8
	add	ebx, esi
	fld	QWORD PTR [edx+ebx*8]
	fld	QWORD PTR $T98441
	lea	ebx, DWORD PTR [edx+ebx*8]

; 582  : 		T2(NNN)=DSQRT(AA2(NNN)+B12(NNN)- ((B12(NNN)**2)/(AA1(NNN)

	call	__FIfexp
	mov	eax, DWORD PTR _7AA1+8
	mov	ecx, DWORD PTR _AA1
	add	eax, esi
	fld	QWORD PTR [ecx+eax*8]
	fadd	QWORD PTR [ebx]
	mov	ecx, DWORD PTR _7AA2+8
	mov	eax, DWORD PTR _AA2
	add	ecx, esi
	fdivp	ST(1), ST(0)
	fsubr	QWORD PTR [eax+ecx*8]
	fadd	QWORD PTR [ebx]
	call	__FIsqrt
	mov	ecx, DWORD PTR _7T2+8
	mov	eax, DWORD PTR _T2
	add	ecx, esi
	mov	edx, DWORD PTR _B12
	fstp	QWORD PTR [eax+ecx*8]
	mov	eax, DWORD PTR _7B12+8
	mov	ecx, DWORD PTR _AA1
	add	eax, esi
	fld	QWORD PTR [edx+eax*8]
	lea	eax, DWORD PTR [edx+eax*8]

; 583  : 	-		+B12(NNN))))
; 584  : 		T3(NNN)=-B12(NNN)/DSQRT(AA1(NNN)+B12(NNN))

	mov	edx, DWORD PTR _7AA1+8
	add	edx, esi
	fld	QWORD PTR [ecx+edx*8]
	fadd	QWORD PTR [eax]
	call	__FIsqrt
	fdivp	ST(1), ST(0)
	mov	eax, DWORD PTR _T3
	mov	ecx, DWORD PTR _7T3+8
	add	ecx, esi
	inc	esi
	test	edi, edi
	fchs

; 579  : 	END DO
; 580  : 	DO NNN=1,M

	jne	$L98467
	fstp	QWORD PTR [eax+ecx*8]
$L442:

; 585  : C	WRITE(6,*)'A1 ->',AA1(NNN),T1(NNN)**2 +T1(NNN)*T3(NNN)
; 586  : C	WRITE(6,*)'A2 ->',AA2(NNN),T2(NNN)**2 +T3(NNN)**2 +T1(NNN)*T3(NNN)
; 587  : C	WRITE(6,*)'B12->',B12(NNN),-T1(NNN)*T3(NNN)
; 588  : C	WRITE(6,*)'-------------------------'
; 589  : 	END DO
; 590  : 	DEFLAG=1

	push	OFFSET FLAT:_DFLAG
	mov	eax, DWORD PTR _T3
	push	eax
	mov	ecx, DWORD PTR _T2
	push	ecx
	mov	edx, DWORD PTR _T1
	push	edx
	mov	eax, DWORD PTR _GRAD
	push	eax
	mov	ecx, DWORD PTR _CO
	push	OFFSET FLAT:_EMIN
	mov	edx, DWORD PTR _H
	push	ecx
	mov	eax, DWORD PTR _S
	push	edx
	mov	ecx, DWORD PTR _EVEC
	push	eax
	mov	edx, DWORD PTR _B12
	push	ecx
	mov	eax, DWORD PTR _RZ2
	push	edx
	mov	ecx, DWORD PTR _RY2
	push	eax
	mov	edx, DWORD PTR _RX2
	push	ecx
	mov	eax, DWORD PTR _AA2
	push	edx
	mov	ecx, DWORD PTR _RZ1
	push	eax
	mov	edx, DWORD PTR _RY1
	push	ecx
	mov	eax, DWORD PTR _RX1
	push	edx
	mov	ecx, DWORD PTR _AA1
	mov	DWORD PTR _DEFLAG, 0
	mov	DWORD PTR _DEFLAG+4, 1072693248		; 3ff00000H

; 591  : 	CALL HAMILTON(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,EVEC,S,H,

	push	eax
	push	ecx
	push	OFFSET FLAT:_M
	call	_HAMILTON@80

; 592  : 	-				CO,EMIN,GRAD,T1,T2,T3,DFLAG)
; 593  : 	WRITE(6,*)'PRE-OPTIMIZATION ENERGY --->',EMIN

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	28					; 0000001cH
	push	OFFSET FLAT:_.SG444_
	call	_lwrite_ch
	add	esp, 8
	mov	ecx, DWORD PTR _EMIN+4
	mov	edx, DWORD PTR _EMIN
	npad	1
	push	ecx
	push	edx
	call	_lwrite_d
	add	esp, 8
	call	_f90_endwrite

; 594  : 	CALL GRADOPT(M,T1,RX1,RY1,RZ1,T2,RX2,RY2,RZ2,T3,EVEC,GRAD)

	mov	ecx, DWORD PTR _GRAD
	mov	edx, DWORD PTR _EVEC
	push	ecx
	mov	eax, DWORD PTR _T3
	push	edx
	mov	ecx, DWORD PTR _RZ2
	push	eax
	mov	edx, DWORD PTR _RY2
	push	ecx
	mov	eax, DWORD PTR _RX2
	push	edx
	mov	ecx, DWORD PTR _T2
	push	eax
	mov	edx, DWORD PTR _RZ1
	push	ecx
	mov	eax, DWORD PTR _RY1
	push	edx
	mov	ecx, DWORD PTR _RX1
	push	eax
	mov	edx, DWORD PTR _T1
	push	ecx
	mov	esi, 1
	push	edx
	push	OFFSET FLAT:_M
	call	_GRADOPT@48

; 595  : 	DO NNN=1,M

	mov	edi, DWORD PTR _M
	mov	ebx, DWORD PTR _T1
	mov	ebp, DWORD PTR _T3
	test	edi, edi
	jle	$L445
	jmp	SHORT $L98323
$L98468:
	fstp	QWORD PTR [ecx+eax*8]
$L98323:
	mov	eax, DWORD PTR _7T1+8
	mov	ecx, DWORD PTR _7T3+8
	add	eax, esi
	add	ecx, esi
	dec	edi
	fld	QWORD PTR [ebp+ecx*8]
	fadd	QWORD PTR [ebx+eax*8]
	lea	eax, DWORD PTR [ebx+eax*8]
	mov	ecx, DWORD PTR _AA1

; 596  : 		AA1(NNN)=T1(NNN)**2 +T1(NNN)*T3(NNN)

	fmul	QWORD PTR [eax]
	mov	eax, DWORD PTR _7AA1+8
	add	eax, esi
	fstp	QWORD PTR [ecx+eax*8]
	mov	eax, DWORD PTR _7T3+8
	mov	ecx, DWORD PTR _7T1+8
	add	eax, esi
	add	ecx, esi
	fld	QWORD PTR [ebx+ecx*8]
	fadd	QWORD PTR [ebp+eax*8]
	lea	eax, DWORD PTR [ebp+eax*8]
	mov	ecx, DWORD PTR _7T2+8

; 597  : 		AA2(NNN)=T2(NNN)**2 + T3(NNN)**2 +T1(NNN)*T3(NNN)

	add	ecx, esi
	fmul	QWORD PTR [eax]
	mov	eax, DWORD PTR _T2
	fld	QWORD PTR [eax+ecx*8]
	fld	QWORD PTR $T98441
	call	__FIfexp
	faddp	ST(1), ST(0)
	mov	eax, DWORD PTR _7AA2+8
	mov	ecx, DWORD PTR _AA2
	add	eax, esi
	fstp	QWORD PTR [ecx+eax*8]

; 598  : 		B12(NNN)=-T1(NNN)*T3(NNN)

	mov	eax, DWORD PTR _7T1+8
	mov	ecx, DWORD PTR _B12
	add	eax, esi
	fld	QWORD PTR [ebx+eax*8]
	mov	eax, DWORD PTR _7T3+8
	add	eax, esi
	fmul	QWORD PTR [ebp+eax*8]
	mov	eax, DWORD PTR _7B12+8
	add	eax, esi
	inc	esi
	fchs
	test	edi, edi

; 595  : 	DO NNN=1,M

	jne	$L98468
	fstp	QWORD PTR [ecx+eax*8]
$L445:

; 599  : 	END DO
; 600  : C
; 601  : 	CALL HAMILTON(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,EVEC,S,H,

	push	OFFSET FLAT:_DFLAG
	mov	eax, DWORD PTR _T2
	push	ebp
	mov	ecx, DWORD PTR _GRAD
	push	eax
	mov	edx, DWORD PTR _CO
	push	ebx
	mov	eax, DWORD PTR _H
	push	ecx
	push	OFFSET FLAT:_EMIN
	push	edx
	mov	ecx, DWORD PTR _S
	push	eax
	mov	edx, DWORD PTR _EVEC
	push	ecx
	mov	eax, DWORD PTR _B12
	push	edx
	mov	ecx, DWORD PTR _RZ2
	push	eax
	mov	edx, DWORD PTR _RY2
	push	ecx
	mov	eax, DWORD PTR _RX2
	push	edx
	mov	ecx, DWORD PTR _AA2
	push	eax
	mov	edx, DWORD PTR _RZ1
	push	ecx
	mov	eax, DWORD PTR _RY1
	push	edx
	mov	ecx, DWORD PTR _RX1
	push	eax
	mov	edx, DWORD PTR _AA1
	push	ecx
	mov	esi, 1
	push	edx
	push	OFFSET FLAT:_M
	call	_HAMILTON@80

; 602  : 	-				CO,EMIN,GRAD,T1,T2,T3,DFLAG)
; 603  : 	DO J=1,M

	mov	edx, DWORD PTR _M
	test	edx, edx
	jle	SHORT $L348
$L98324:

; 604  : 		EVEC(J)=CO(J)

	mov	eax, DWORD PTR _7CO+8
	mov	ebp, DWORD PTR _CO
	add	eax, esi
	mov	edi, DWORD PTR _7EVEC+8
	add	edi, esi
	mov	ebx, DWORD PTR _EVEC
	mov	ecx, DWORD PTR [ebp+eax*8+4]
	dec	edx
	mov	eax, DWORD PTR [ebp+eax*8]
	inc	esi
	mov	DWORD PTR [ebx+edi*8+4], ecx
	test	edx, edx
	mov	DWORD PTR [ebx+edi*8], eax

; 602  : 	-				CO,EMIN,GRAD,T1,T2,T3,DFLAG)
; 603  : 	DO J=1,M

	jne	SHORT $L98324
$L348:

; 368  : C
; 369  : 	DO JJJ=M,ISTEND

	mov	eax, DWORD PTR -36+[esp+96]
	dec	eax
	mov	DWORD PTR -36+[esp+96], eax
	jne	$L98315
$L344:

; 605  : 	END DO
; 606  : C
; 607  : C	END IF
; 608  : 	END DO
; 609  : C
; 610  : C
; 611  : 	END IF
; 612  : C
; 613  : C
; 614  : C
; 615  : 	IF (STGROW .EQ. 1) THEN

	mov	eax, DWORD PTR _STGROW
	cmp	eax, 1
	jne	$L449

; 616  : 	WRITE(6,*)'STOCHASTIC GROWTH PERCENTAGE=',STGROWPER/100.0D+00

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	29					; 0000001dH
	push	OFFSET FLAT:_.SG451_
	call	_lwrite_ch
	fld	QWORD PTR _STGROWPER
	fdiv	QWORD PTR $T98457
	add	esp, 8
	sub	esp, 8
	fstp	QWORD PTR [esp]
	call	_lwrite_d
	add	esp, 8
	call	_f90_endwrite

; 617  : 	WRITE(7,*)'STOCHASTIC GROWTH PERCENTAGE=',STGROWPER/100.0D+00

	push	0
	push	0
	push	0
	push	0
	push	0
	push	7
	call	_f90_write
	add	esp, 24					; 00000018H
	push	29					; 0000001dH
	push	OFFSET FLAT:_.SG452_
	call	_lwrite_ch
	fld	QWORD PTR _STGROWPER
	fdiv	QWORD PTR $T98457
	add	esp, 8
	sub	esp, 8
	fstp	QWORD PTR [esp]
	call	_lwrite_d
	add	esp, 8
	call	_f90_endwrite

; 618  : 	CALL BIOUT(' ')

	push	1
	push	OFFSET FLAT:_.SG453_
	call	_BIOUT@8

; 619  : C
; 620  : 	DO WHILE (M .LE. ISTEND)

	mov	eax, DWORD PTR _ISTEND
	mov	ecx, DWORD PTR _M
	cmp	ecx, eax
	jg	$L449
$L455:

; 621  : 		ITEMP=INT((STGROWPER/100.0D+00)*M)

	fild	DWORD PTR _M
	fld	QWORD PTR _STGROWPER
	fdiv	QWORD PTR $T98457
	fmulp	ST(1), ST(0)
	call	__ftol
	mov	DWORD PTR _ITEMP, eax
	test	eax, eax

; 622  : 		IF (ITEMP .LE. 0) ITEMP=1

	jg	SHORT $L457
	mov	DWORD PTR _ITEMP, 1
$L457:

; 623  : 		IF (ITEMP+M .GT. ISTEND) ITEMP=ISTEND-M

	mov	eax, DWORD PTR _ITEMP
	mov	ecx, DWORD PTR _M
	add	eax, ecx
	mov	edx, DWORD PTR _ISTEND
	cmp	eax, edx
	jle	SHORT $L459
	mov	eax, edx
	sub	eax, ecx
	mov	DWORD PTR _ITEMP, eax
$L459:

; 624  : 	DO KKK=1,ITEMP

	mov	eax, DWORD PTR _ITEMP
	mov	DWORD PTR -36+[esp+96], eax
	test	eax, eax
	jle	$L462
$L98325:

; 625  : C
; 626  : C
; 627  : 	DO I=1,IVIEW

	mov	DWORD PTR _I$[esp+96], 1
	mov	ebp, DWORD PTR _IVIEW
	test	ebp, ebp
	jle	$L464
	xor	esi, esi
	mov	edi, 1072693248				; 3ff00000H
	mov	ebx, -1074790400			; bff00000H
$L98326:

; 628  : 	IF (TEMP1.GE..5) THEN

	fld	QWORD PTR _TEMP1
	fcomp	QWORD PTR $T98454
	fwait
	fnstsw	ax
	test	ah, 1

; 629  : C
; 630  : 		CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	jne	$L467
	call	_RANDOM@4

; 631  : 		TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1
$L468:

; 632  :  8141		CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	call	_RANDOM@4

; 633  : 		TEMP2=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP2

; 634  : 	IF (TEMP2.LT..5) GOTO 8141

	fld	QWORD PTR _TEMP2
	fcomp	QWORD PTR $T98454
	fwait
	fnstsw	ax
	test	ah, 1
	jne	SHORT $L468

; 635  : 	TOF1=(FMAX-FMIN)*TEMP1+FMIN

	fld	QWORD PTR _FMAX
	fsub	QWORD PTR _FMIN
	fld	QWORD PTR $T98443
	fxch	ST(1)
	fmul	QWORD PTR _TEMP1
	fadd	QWORD PTR _FMIN
	fstp	QWORD PTR _TOF1

; 636  : 	XAA1=TEMP2*10**TOF1

	fld	QWORD PTR _TOF1
	call	__FIfexp
	fmul	QWORD PTR _TEMP2
	push	OFFSET FLAT:_TEMP
	fstp	QWORD PTR _XAA1

; 637  : C
; 638  : 		CALL RANDOM(TEMP)

	call	_RANDOM@4

; 639  : 		TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1
$L471:

; 640  :  8151		CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	call	_RANDOM@4

; 641  : 		TEMP2=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP2

; 642  : 	IF (TEMP2.LT..5) GOTO 8151

	fld	QWORD PTR _TEMP2
	fcomp	QWORD PTR $T98454
	fwait
	fnstsw	ax
	test	ah, 1
	jne	SHORT $L471
	fld	QWORD PTR _FMAX
	fsub	QWORD PTR _FMIN
	fstp	QWORD PTR -16+[esp+96]

; 643  : 	TOF2=(FMAX-FMIN)*TEMP1+FMIN

	fld	QWORD PTR -16+[esp+96]
	fmul	QWORD PTR _TEMP1
	fadd	QWORD PTR _FMIN
	fstp	QWORD PTR _TOF2

; 644  : 	XAA2=TEMP2*10**TOF2

	fld	QWORD PTR $T98443
	fld	QWORD PTR _TOF2
	call	__FIfexp
	fmul	QWORD PTR _TEMP2
	fstp	QWORD PTR _XAA2
	fld	QWORD PTR _FMIN
	fabs
	fstp	QWORD PTR -8+[esp+96]
	fild	DWORD PTR _K$[esp+96]
	fstp	QWORD PTR -24+[esp+96]

; 645  : C
; 646  : C	XXX1=(TOF1+DABS(FMIN)+1.+D-(DABS(FMIN)+1.))/C
; 647  : C	XXX2=(TOF2+DABS(FMIN)+1.+D-(DABS(FMIN)+1.))/C
; 648  : C	XRZ1=-(-SIGMAP *DLOG(XXX1))+ZNUC
; 649  : 	XRZ1=-DLOG((TOF1+DABS(FMIN))/(FMAX-FMIN))/K + ZNUC

	fld	QWORD PTR -8+[esp+96]
	fadd	QWORD PTR _TOF1
	fdiv	QWORD PTR -16+[esp+96]
	call	__FIlog
	fdiv	QWORD PTR -24+[esp+96]
	fsubr	QWORD PTR _ZNUC
	fstp	QWORD PTR _XRZ1

; 650  : C	XRZ2=-(-SIGMAP *DLOG(XXX2))+ZNUC
; 651  : 	XRZ2=-DLOG((TOF2+DABS(FMIN))/(FMAX-FMIN))/K + ZNUC

	fld	QWORD PTR _TOF2
	fadd	QWORD PTR -8+[esp+96]
	fdiv	QWORD PTR -16+[esp+96]
	call	__FIlog
	fdiv	QWORD PTR -24+[esp+96]
	push	OFFSET FLAT:_TEMP
	fsubr	QWORD PTR _ZNUC
	fstp	QWORD PTR _XRZ2

; 652  : 	CALL RANDOM(TEMP)

	call	_RANDOM@4

; 653  : 	TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1

; 654  : 	IF (TEMP1.GE..5) XRZ1=-1.*XRZ1

	fld	QWORD PTR _TEMP1
	fcomp	QWORD PTR $T98454
	fwait
	fnstsw	ax
	test	ah, 1
	jne	SHORT $L474
	fld	QWORD PTR _XRZ1
	fchs
	fstp	QWORD PTR _XRZ1
$L474:

; 655  : 	CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	call	_RANDOM@4

; 656  : 	TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1

; 657  : 	IF (TEMP1.GE..5) XRZ2=-1*XRZ2

	fld	QWORD PTR _TEMP1
	fcomp	QWORD PTR $T98454
	fwait
	fnstsw	ax
	test	ah, 1
	jne	SHORT $L476
	fld	QWORD PTR _XRZ2
	fchs
	fstp	QWORD PTR _XRZ2
$L476:

; 658  : 	IF (BL0.GT.0.) THEN

	fld	QWORD PTR _BL0
	fcomp	QWORD PTR $T98428
	fwait
	fnstsw	ax
	test	ah, 65					; 00000041H
	jne	$L466

; 659  : 		CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	call	_RANDOM@4

; 660  : 	TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1

; 661  : 		XC=1

	fld	QWORD PTR _TEMP1
	fcomp	QWORD PTR $T98454
	fwait
	fnstsw	ax
	test	ah, 1
	mov	DWORD PTR _XC, esi
	mov	DWORD PTR _XC+4, edi

; 662  : 		IF (TEMP1.GE..5) XC=-1.

	jne	SHORT $L480
	mov	DWORD PTR _XC, esi
	mov	DWORD PTR _XC+4, ebx
$L480:

; 663  : 	CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	call	_RANDOM@4

; 664  : 		TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1

; 665  : 		IF (TOF1 .LE. BLURBEGIN) THEN

	fld	QWORD PTR _BLURBEGIN
	fcomp	QWORD PTR _TOF1
	fwait
	fnstsw	ax
	test	ah, 1
	jne	SHORT $L482

; 666  : 		XRZ1=XRZ1+XC*DSQRT(-2.*BL0**2 * DLOG(TEMP1))* 

	fld	QWORD PTR _BLURBEGIN
	fsub	QWORD PTR _TOF1
	fld	QWORD PTR _TEMP1
	call	__FIlog
	fld	QWORD PTR _BL0
	fld	QWORD PTR $T98441
	call	__FIfexp
	fmulp	ST(1), ST(0)
	fmul	QWORD PTR $T98456
	call	__FIsqrt
	fmulp	ST(1), ST(0)
	fld	QWORD PTR _BLURBEGIN
	fsub	QWORD PTR _BLUREND
	fxch	ST(1)
	fmul	QWORD PTR _XC
	fdivrp	ST(1), ST(0)
	fadd	QWORD PTR _XRZ1
	fstp	QWORD PTR _XRZ1
$L482:

; 667  :  	-		(BLURBEGIN-TOF1)/(BLURBEGIN-BLUREND)
; 668  : 		END IF	
; 669  : C
; 670  : 		CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	call	_RANDOM@4

; 671  : 	TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1

; 672  : 		XC=1

	fld	QWORD PTR _TEMP1
	fcomp	QWORD PTR $T98454
	fwait
	fnstsw	ax
	test	ah, 1
	mov	DWORD PTR _XC, esi
	mov	DWORD PTR _XC+4, edi

; 673  : 		IF (TEMP1.GE..5) XC=-1.

	jne	SHORT $L484
	mov	DWORD PTR _XC, esi
	mov	DWORD PTR _XC+4, ebx
$L484:

; 674  : 	CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	call	_RANDOM@4

; 675  : 		TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1

; 676  : 		IF (TOF2 .LE. BLURBEGIN) THEN

	fld	QWORD PTR _BLURBEGIN
	fcomp	QWORD PTR _TOF2
	fwait
	fnstsw	ax
	test	ah, 1
	jne	$L466

; 677  : 		XRZ2=XRZ2+XC*DSQRT(-2.*BL0**2 * DLOG(TEMP1))* 

	fld	QWORD PTR _BLURBEGIN
	fsub	QWORD PTR _TOF2
	fld	QWORD PTR _TEMP1
	call	__FIlog
	fld	QWORD PTR _BL0
	jmp	$L98463
$L467:

; 678  :   	-		(BLURBEGIN-TOF2)/(BLURBEGIN-BLUREND)
; 679  : 		END IF
; 680  : C
; 681  : 	END IF
; 682  : C
; 683  : C
; 684  : 	ELSE

	call	_RANDOM@4

; 685  : 		CALL RANDOM(TEMP)
; 686  : 		TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1
$L488:

; 687  :  9121		CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	call	_RANDOM@4

; 688  : 		TEMP2=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP2

; 689  : 	IF (TEMP2.LT..5) GOTO 9121

	fld	QWORD PTR _TEMP2
	fcomp	QWORD PTR $T98454
	fwait
	fnstsw	ax
	test	ah, 1
	jne	SHORT $L488

; 690  : 	TOF1=(FMAX-FMIN)*TEMP1+FMIN

	fld	QWORD PTR _FMAX
	fsub	QWORD PTR _FMIN
	fld	QWORD PTR $T98443
	fxch	ST(1)
	fmul	QWORD PTR _TEMP1
	fadd	QWORD PTR _FMIN
	fstp	QWORD PTR _TOF1

; 691  : 	XAA1=TEMP2*10**TOF1

	fld	QWORD PTR _TOF1
	call	__FIfexp
	fmul	QWORD PTR _TEMP2
	push	OFFSET FLAT:_TEMP
	fstp	QWORD PTR _XAA1

; 692  : C
; 693  : 		CALL RANDOM(TEMP)

	call	_RANDOM@4

; 694  : 		TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1
$L491:

; 695  :  9131		CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	call	_RANDOM@4

; 696  : 		TEMP2=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP2

; 697  : 	IF (TEMP2.LT..5) GOTO 9131

	fld	QWORD PTR _TEMP2
	fcomp	QWORD PTR $T98454
	fwait
	fnstsw	ax
	test	ah, 1
	jne	SHORT $L491
	fld	QWORD PTR _FMAX
	fsub	QWORD PTR _FMIN
	fstp	QWORD PTR -24+[esp+96]

; 698  : 	TOF2=(FMAX-FMIN)*TEMP1+FMIN

	fld	QWORD PTR -24+[esp+96]
	fmul	QWORD PTR _TEMP1
	fadd	QWORD PTR _FMIN
	fstp	QWORD PTR _TOF2

; 699  : 	XAA2=TEMP2*10**TOF2

	fld	QWORD PTR $T98443
	fld	QWORD PTR _TOF2
	call	__FIfexp
	fmul	QWORD PTR _TEMP2
	fstp	QWORD PTR _XAA2
	fld	QWORD PTR _FMIN
	fabs
	fstp	QWORD PTR -8+[esp+96]
	fild	DWORD PTR _KA
	fstp	QWORD PTR -16+[esp+96]

; 700  : C
; 701  : C	XXX1=(TOF1+DABS(FMIN)+1.+B-(DABS(FMIN)+1.))/A
; 702  : C	XXX2=(TOF2+DABS(FMIN)+1.+B-(DABS(FMIN)+1.))/A
; 703  : C	XRZ1=(-SIGMA *DLOG(XXX1))+ZNUC
; 704  : C	XRZ2=(-SIGMA *DLOG(XXX2))+ZNUC
; 705  : 	XRZ1=-DLOG((TOF1+DABS(FMIN))/(FMAX-FMIN))/KA - ZNUC

	fld	QWORD PTR -8+[esp+96]
	fadd	QWORD PTR _TOF1
	fdiv	QWORD PTR -24+[esp+96]
	call	__FIlog
	fdiv	QWORD PTR -16+[esp+96]
	fadd	QWORD PTR _ZNUC
	fchs
	fstp	QWORD PTR _XRZ1

; 706  : 	XRZ2=-DLOG((TOF2+DABS(FMIN))/(FMAX-FMIN))/KA - ZNUC

	fld	QWORD PTR _TOF2
	fadd	QWORD PTR -8+[esp+96]
	fdiv	QWORD PTR -24+[esp+96]
	call	__FIlog
	fdiv	QWORD PTR -16+[esp+96]
	push	OFFSET FLAT:_TEMP
	fadd	QWORD PTR _ZNUC
	fchs
	fstp	QWORD PTR _XRZ2

; 707  : 	CALL RANDOM(TEMP)

	call	_RANDOM@4

; 708  : 	TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1

; 709  : 	IF (TEMP1.GE..5) XRZ1=-1.*XRZ1

	fld	QWORD PTR _TEMP1
	fcomp	QWORD PTR $T98454
	fwait
	fnstsw	ax
	test	ah, 1
	jne	SHORT $L494
	fld	QWORD PTR _XRZ1
	fchs
	fstp	QWORD PTR _XRZ1
$L494:

; 710  : 	CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	call	_RANDOM@4

; 711  : 	TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1

; 712  : 	IF (TEMP1.GE..5) XRZ2=-1*XRZ2

	fld	QWORD PTR _TEMP1
	fcomp	QWORD PTR $T98454
	fwait
	fnstsw	ax
	test	ah, 1
	jne	SHORT $L496
	fld	QWORD PTR _XRZ2
	fchs
	fstp	QWORD PTR _XRZ2
$L496:

; 713  : C
; 714  : 	IF (BL1.GT.0.) THEN

	fld	QWORD PTR _BL1
	fcomp	QWORD PTR $T98428
	fwait
	fnstsw	ax
	test	ah, 65					; 00000041H
	jne	$L466

; 715  : 		CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	call	_RANDOM@4

; 716  : 	TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1

; 717  : 		XC=1

	fld	QWORD PTR _TEMP1
	fcomp	QWORD PTR $T98454
	fwait
	fnstsw	ax
	test	ah, 1
	mov	DWORD PTR _XC, esi
	mov	DWORD PTR _XC+4, edi

; 718  : 		IF (TEMP1.GE..5) XC=-1.

	jne	SHORT $L500
	mov	DWORD PTR _XC, esi
	mov	DWORD PTR _XC+4, ebx
$L500:

; 719  : 	CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	call	_RANDOM@4

; 720  : 		TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1

; 721  : 		IF (TOF1 .LE. BLURBEGIN) THEN

	fld	QWORD PTR _BLURBEGIN
	fcomp	QWORD PTR _TOF1
	fwait
	fnstsw	ax
	test	ah, 1
	jne	SHORT $L502

; 722  : 		XRZ1=XRZ1+XC*DSQRT(-2.*BL1**2 * DLOG(TEMP1))* 

	fld	QWORD PTR _BLURBEGIN
	fsub	QWORD PTR _TOF1
	fld	QWORD PTR _TEMP1
	call	__FIlog
	fld	QWORD PTR _BL1
	fld	QWORD PTR $T98441
	call	__FIfexp
	fmulp	ST(1), ST(0)
	fmul	QWORD PTR $T98456
	call	__FIsqrt
	fmulp	ST(1), ST(0)
	fld	QWORD PTR _BLURBEGIN
	fsub	QWORD PTR _BLUREND
	fxch	ST(1)
	fmul	QWORD PTR _XC
	fdivrp	ST(1), ST(0)
	fadd	QWORD PTR _XRZ1
	fstp	QWORD PTR _XRZ1
$L502:

; 723  :   	-		(BLURBEGIN-TOF1)/(BLURBEGIN-BLUREND)
; 724  : 		END IF
; 725  : 	
; 726  : C
; 727  : 		CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	call	_RANDOM@4

; 728  : 	TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1

; 729  : 		XC=1

	fld	QWORD PTR _TEMP1
	fcomp	QWORD PTR $T98454
	fwait
	fnstsw	ax
	test	ah, 1
	mov	DWORD PTR _XC, esi
	mov	DWORD PTR _XC+4, edi

; 730  : 		IF (TEMP1.GE..5) XC=-1.

	jne	SHORT $L504
	mov	DWORD PTR _XC, esi
	mov	DWORD PTR _XC+4, ebx
$L504:

; 731  : 	CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	call	_RANDOM@4

; 732  : 		TEMP1=TEMP

	fld	DWORD PTR _TEMP
	fstp	QWORD PTR _TEMP1

; 733  : 		IF (TOF2 .LE. BLURBEGIN) THEN

	fld	QWORD PTR _BLURBEGIN
	fcomp	QWORD PTR _TOF2
	fwait
	fnstsw	ax
	test	ah, 1
	jne	SHORT $L466

; 734  : 		XRZ2=XRZ2+XC*DSQRT(-2.*BL1**2 * DLOG(TEMP1))* 

	fld	QWORD PTR _BLURBEGIN
	fsub	QWORD PTR _TOF2
	fld	QWORD PTR _TEMP1
	call	__FIlog
	fld	QWORD PTR _BL1
$L98463:
	fld	QWORD PTR $T98441
	call	__FIfexp
	fmulp	ST(1), ST(0)
	fmul	QWORD PTR $T98456
	call	__FIsqrt
	fmulp	ST(1), ST(0)
	fmul	QWORD PTR _XC
	fld	QWORD PTR _BLURBEGIN
	fsub	QWORD PTR _BLUREND
	fdivp	ST(1), ST(0)
	fadd	QWORD PTR _XRZ2
	fstp	QWORD PTR _XRZ2
$L466:

; 735  :   	-		(BLURBEGIN-TOF2)/(BLURBEGIN-BLUREND)
; 736  : 
; 737  : 		END IF
; 738  : C
; 739  : 	END IF
; 740  : 
; 741  : 
; 742  : C
; 743  : 	END IF
; 744  : 	CALL RANDOM(TEMP)

	push	OFFSET FLAT:_TEMP
	call	_RANDOM@4

; 745  : 	TEMP1=TEMP

	fld	DWORD PTR _TEMP
	push	OFFSET FLAT:_TEMP
	fstp	QWORD PTR _TEMP1

; 746  : 	CALL RANDOM(TEMP)

	call	_RANDOM@4

; 747  : 	TEMP2=TEMP

	fld	DWORD PTR _TEMP
	fld1
	fadd	QWORD PTR _BETAMIN
	fxch	ST(1)
	fstp	QWORD PTR _TEMP2
	fstp	QWORD PTR -8+[esp+96]

; 748  : 	XB12=DABS(TEMP1*10.**(((BETAMAX+1.)-(BETAMIN+1))*

	fld	QWORD PTR $T98443
	fld1
	fadd	QWORD PTR _BETAMAX
	fsub	QWORD PTR -8+[esp+96]
	fmul	QWORD PTR _TEMP2
	fadd	QWORD PTR -8+[esp+96]
	call	__FIfexp
	fmul	QWORD PTR _TEMP1
	push	OFFSET FLAT:_XB12
	push	OFFSET FLAT:_XRZ2
	push	OFFSET FLAT:_XRY2
	push	OFFSET FLAT:_XRX2
	fabs
	push	OFFSET FLAT:_XAA2
	push	OFFSET FLAT:_XRZ1
	fstp	QWORD PTR _XB12

; 749  : 	-TEMP2+(BETAMIN+1)))
; 750  : C	WRITE(6,*)'BETA=',XB12
; 751  : C	WRITE(6,292)I,XAA1,XAA2,XRZ1,XRZ2
; 752  : 	CALL HAMILTONp(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,

	push	OFFSET FLAT:_XRY1
	mov	eax, DWORD PTR _H
	push	OFFSET FLAT:_XRX1
	mov	ecx, DWORD PTR _S
	push	OFFSET FLAT:_XAA1
	mov	edx, DWORD PTR _EVEC
	push	OFFSET FLAT:_EMIN
	push	eax
	push	ecx
	mov	eax, DWORD PTR _B12
	push	edx
	mov	ecx, DWORD PTR _RZ2
	push	eax
	mov	edx, DWORD PTR _RY2
	push	ecx
	mov	eax, DWORD PTR _RX2
	push	edx
	mov	ecx, DWORD PTR _AA2
	push	eax
	mov	edx, DWORD PTR _RZ1
	push	ecx
	mov	eax, DWORD PTR _RY1
	push	edx
	mov	ecx, DWORD PTR _RX1
	push	eax
	mov	edx, DWORD PTR _AA1
	push	ecx
	push	edx
	push	OFFSET FLAT:_M
	call	_HAMILTONP@92

; 753  :      1                    EVEC,S,H,EMIN,XAA1,XRX1,XRY1,XRZ1,
; 754  :      2					XAA2,XRX2,XRY2,XRZ2,XB12)
; 755  : 	IF (I .EQ. 1) ENER=EMIN

	mov	ecx, DWORD PTR _I$[esp+96]
	cmp	ecx, 1
	jne	SHORT $L508
	mov	eax, DWORD PTR _EMIN+4
	mov	ecx, DWORD PTR _EMIN
	mov	DWORD PTR _ENER+4, eax
	mov	DWORD PTR _ENER, ecx
$L508:

; 756  : 	IF (EMIN.LT.ENER) THEN 

	fld	QWORD PTR _ENER
	fcomp	QWORD PTR _EMIN
	fwait
	fnstsw	ax
	test	ah, 65					; 00000041H
	jne	$L465

; 757  : 		ENER=EMIN

	mov	eax, DWORD PTR _EMIN+4
	mov	ecx, DWORD PTR _EMIN
	mov	edx, DWORD PTR _XAA1+4
	mov	DWORD PTR _ENER+4, eax
	mov	eax, DWORD PTR _XAA1
	mov	DWORD PTR _ENER, ecx

; 758  : 		TA1=XAA1

	mov	ecx, DWORD PTR _XAA2+4
	mov	DWORD PTR _TA1+4, edx
	mov	edx, DWORD PTR _XAA2
	mov	DWORD PTR _TA1, eax

; 759  : 		TA2=XAA2

	mov	eax, DWORD PTR _XB12+4
	mov	DWORD PTR _TA2+4, ecx
	mov	ecx, DWORD PTR _XB12
	mov	DWORD PTR _TA2, edx

; 760  : 		TB12=XB12

	mov	edx, DWORD PTR _XRZ1+4
	mov	DWORD PTR _TB12+4, eax
	mov	eax, DWORD PTR _XRZ1
	mov	DWORD PTR _TB12, ecx

; 761  : 		TRZ1=XRZ1

	mov	ecx, DWORD PTR _XRZ2+4
	mov	DWORD PTR _TRZ1+4, edx
	mov	edx, DWORD PTR _XRZ2
	mov	DWORD PTR _TRZ1, eax

; 762  : 		TRZ2=XRZ2

	mov	DWORD PTR _TRZ2+4, ecx
	mov	DWORD PTR _TRZ2, edx
$L465:

; 625  : C
; 626  : C
; 627  : 	DO I=1,IVIEW

	dec	ebp
	mov	eax, DWORD PTR _I$[esp+96]
	inc	eax
	test	ebp, ebp
	mov	DWORD PTR _I$[esp+96], eax
	jne	$L98326
$L464:

; 763  : 	END IF
; 764  :       END DO
; 765  : 	CALL BIOUT(' ')

	push	1
	push	OFFSET FLAT:_.SG512_
	call	_BIOUT@8

; 766  : 	WRITE(6,*)'BEST DIAGNOSTIC VALUE IN ',IVIEW,' GUESSES=',ENER

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	25					; 00000019H
	push	OFFSET FLAT:_.SG513_
	call	_lwrite_ch
	add	esp, 8
	mov	eax, DWORD PTR _IVIEW
	push	eax
	call	_lwrite_i
	add	esp, 4
	push	9
	push	OFFSET FLAT:_.SG514_
	call	_lwrite_ch
	add	esp, 8
	mov	eax, DWORD PTR _ENER+4
	mov	ecx, DWORD PTR _ENER
	npad	1
	push	eax
	push	ecx
	call	_lwrite_d
	add	esp, 8
	call	_f90_endwrite

; 767  : C	WRITE(6,*)TA1,TRZ1,TA2,TRZ2,TB12
; 768  : 	WRITE(6,*)'***************************'

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	27					; 0000001bH
	push	OFFSET FLAT:_.SG515_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 769  : 	WRITE(6,*)'ADDING THIS FUNCTION TO Psi'

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	27					; 0000001bH
	push	OFFSET FLAT:_.SG516_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 770  : 	WRITE(6,*)'***************************'

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	27					; 0000001bH
	push	OFFSET FLAT:_.SG517_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 771  : 	ALLOCATE(EVECSTORET(M),TSTORE1T(M),

	mov	DWORD PTR _.T518_$[esp+96], 0
	mov	ecx, DWORD PTR _M
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7TSTORE1T
	call	_make_info
	lea	ecx, DWORD PTR _.T518_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7TSTORE1T
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _TSTORE1T, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7TSTORE2T
	call	_make_info
	lea	ecx, DWORD PTR _.T518_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7TSTORE2T
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _TSTORE2T, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7TSTORE3T
	call	_make_info
	lea	ecx, DWORD PTR _.T518_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7TSTORE3T
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _TSTORE3T, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7RZ1STORET
	call	_make_info
	lea	ecx, DWORD PTR _.T518_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7RZ1STORET
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _RZ1STORET, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7RZ2STORET
	call	_make_info
	lea	ecx, DWORD PTR _.T518_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7RZ2STORET
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _RZ2STORET, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7EVECSTORET
	call	_make_info
	lea	ecx, DWORD PTR _.T518_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7EVECSTORET
	call	__Allocate_i
	mov	ecx, DWORD PTR _.T518_$[esp+108]
	add	esp, 12					; 0000000cH
	mov	edx, 1
	mov	DWORD PTR _EVECSTORET, eax
	mov	DWORD PTR _ERROR, ecx

; 772  :      -			TSTORE2T(M),TSTORE3T(M),
; 773  :      -			RZ1STORET(M),RZ2STORET(M),
; 774  :      -			STAT=ERROR)
; 775  : 	DO III=1,M

	mov	ecx, DWORD PTR _M
	test	ecx, ecx
	jle	$L519
$L98327:

; 776  : 		EVECSTORET(III)=EVEC(III)

	mov	ebp, DWORD PTR _7EVEC+8
	mov	eax, DWORD PTR _EVEC
	add	ebp, edx
	mov	ebx, DWORD PTR _7EVECSTORET+8
	add	ebx, edx
	mov	edi, DWORD PTR _EVECSTORET
	mov	esi, DWORD PTR [eax+ebp*8+4]
	dec	ecx
	mov	eax, DWORD PTR [eax+ebp*8]
	mov	DWORD PTR [edi+ebx*8+4], esi
	mov	DWORD PTR [edi+ebx*8], eax

; 777  : 		TSTORE1T(III)=AA1(III)

	mov	eax, DWORD PTR _AA1
	mov	ebp, DWORD PTR _7AA1+8
	mov	ebx, DWORD PTR _7TSTORE1T+8
	add	ebp, edx
	add	ebx, edx
	mov	edi, DWORD PTR _TSTORE1T
	npad	1
	mov	esi, DWORD PTR [eax+ebp*8+4]
	mov	eax, DWORD PTR [eax+ebp*8]
	mov	DWORD PTR [edi+ebx*8+4], esi
	mov	esi, DWORD PTR _AA2
	mov	DWORD PTR [edi+ebx*8], eax

; 778  : 		TSTORE2T(III)=AA2(III)

	mov	eax, DWORD PTR _7AA2+8
	add	eax, edx
	mov	edi, DWORD PTR _7TSTORE2T+8
	add	edi, edx
	mov	ebx, DWORD PTR _TSTORE2T
	mov	ebp, DWORD PTR [esi+eax*8+4]
	mov	eax, DWORD PTR [esi+eax*8]
	mov	esi, DWORD PTR _B12
	mov	DWORD PTR [ebx+edi*8+4], ebp
	mov	DWORD PTR [ebx+edi*8], eax

; 779  : 		TSTORE3T(III)=B12(III)

	mov	eax, DWORD PTR _7B12+8
	add	eax, edx
	mov	ebx, DWORD PTR _7TSTORE3T+8
	add	ebx, edx
	mov	edi, DWORD PTR _TSTORE3T
	mov	ebp, DWORD PTR [esi+eax*8+4]
	mov	eax, DWORD PTR [esi+eax*8]
	mov	DWORD PTR [edi+ebx*8+4], ebp
	mov	DWORD PTR [edi+ebx*8], eax

; 780  : 		RZ1STORET(III)=RZ1(III)

	mov	eax, DWORD PTR _RZ1
	mov	esi, DWORD PTR _7RZ1+8
	mov	ebx, DWORD PTR _7RZ1STORET+8
	add	esi, edx
	add	ebx, edx
	mov	edi, DWORD PTR _RZ1STORET
	npad	1
	mov	ebp, DWORD PTR [eax+esi*8+4]
	mov	eax, DWORD PTR [eax+esi*8]
	mov	DWORD PTR [edi+ebx*8+4], ebp
	mov	DWORD PTR [edi+ebx*8], eax

; 781  : 		RZ2STORET(III)=RZ2(III)

	mov	eax, DWORD PTR _RZ2
	mov	ebp, DWORD PTR _7RZ2+8
	mov	ebx, DWORD PTR _7RZ2STORET+8
	add	ebp, edx
	add	ebx, edx
	mov	edi, DWORD PTR _RZ2STORET
	inc	edx
	mov	esi, DWORD PTR [eax+ebp*8+4]
	test	ecx, ecx
	mov	eax, DWORD PTR [eax+ebp*8]
	mov	DWORD PTR [edi+ebx*8+4], esi
	mov	DWORD PTR [edi+ebx*8], eax

; 772  :      -			TSTORE2T(M),TSTORE3T(M),
; 773  :      -			RZ1STORET(M),RZ2STORET(M),
; 774  :      -			STAT=ERROR)
; 775  : 	DO III=1,M

	jne	$L98327
$L519:

; 782  : 	END DO
; 783  : 	M=M+1

	mov	eax, DWORD PTR _M
	push	0
	inc	eax
	push	8
	push	OFFSET FLAT:_7RX1
	mov	ecx, DWORD PTR _RX1
	push	ecx
	mov	DWORD PTR _M, eax

; 784  : 	DEALLOCATE (AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,EVEC,S,H,CO,

	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _RY1
	mov	DWORD PTR _RX1, 0
	push	0
	push	8
	push	OFFSET FLAT:_7RY1
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _RZ1
	mov	DWORD PTR _RY1, 0
	push	0
	push	8
	push	OFFSET FLAT:_7RZ1
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _AA2
	mov	DWORD PTR _RZ1, 0
	push	0
	push	8
	push	OFFSET FLAT:_7AA2
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _RX2
	mov	DWORD PTR _AA2, 0
	push	0
	push	8
	push	OFFSET FLAT:_7RX2
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _RY2
	mov	DWORD PTR _RX2, 0
	push	0
	push	8
	push	OFFSET FLAT:_7RY2
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _RZ2
	mov	DWORD PTR _RY2, 0
	push	0
	push	8
	push	OFFSET FLAT:_7RZ2
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _B12
	mov	DWORD PTR _RZ2, 0
	push	0
	push	8
	push	OFFSET FLAT:_7B12
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _EVEC
	mov	DWORD PTR _B12, 0
	push	0
	push	8
	push	OFFSET FLAT:_7EVEC
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _S
	mov	DWORD PTR _EVEC, 0
	push	0
	push	8
	push	OFFSET FLAT:_7S
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _H
	mov	DWORD PTR _S, 0
	push	0
	push	8
	push	OFFSET FLAT:_7H
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _CO
	mov	DWORD PTR _H, 0
	push	0
	push	8
	push	OFFSET FLAT:_7CO
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	DWORD PTR _CO, 0
	push	0
	push	8
	push	OFFSET FLAT:_7GRAD
	mov	ecx, DWORD PTR _GRAD
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _T1
	mov	DWORD PTR _GRAD, 0
	push	0
	push	8
	push	OFFSET FLAT:_7T1
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _T2
	mov	DWORD PTR _T1, 0
	push	0
	push	8
	push	OFFSET FLAT:_7T2
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _T3
	mov	DWORD PTR _T2, 0
	push	0
	push	8
	push	OFFSET FLAT:_7T3
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _AA1
	mov	DWORD PTR _T3, 0
	push	0
	push	8
	push	OFFSET FLAT:_7AA1
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _AA1, 0

; 785  : 	-			GRAD,T1,T2,T3)
; 786  :       ALLOCATE(	AA1(M),RX1(M),RY1(M),RZ1(M),

	mov	DWORD PTR _.T521_$[esp+96], 0
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7RX1
	call	_make_info
	lea	ecx, DWORD PTR _.T521_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7RX1
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _RX1, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7RY1
	call	_make_info
	lea	ecx, DWORD PTR _.T521_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7RY1
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _RY1, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7RZ1
	call	_make_info
	lea	ecx, DWORD PTR _.T521_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7RZ1
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _RZ1, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7AA2
	call	_make_info
	lea	ecx, DWORD PTR _.T521_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7AA2
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _AA2, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7RX2
	call	_make_info
	lea	ecx, DWORD PTR _.T521_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7RX2
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _RX2, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7RY2
	call	_make_info
	lea	ecx, DWORD PTR _.T521_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7RY2
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _RY2, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7RZ2
	call	_make_info
	lea	ecx, DWORD PTR _.T521_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7RZ2
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _RZ2, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7B12
	call	_make_info
	lea	ecx, DWORD PTR _.T521_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7B12
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _B12, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7EVEC
	call	_make_info
	lea	ecx, DWORD PTR _.T521_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7EVEC
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _EVEC, eax
	npad	1
	push	ecx
	push	1
	push	ecx
	push	1
	push	2
	push	OFFSET FLAT:_7S
	call	_make_info
	lea	ecx, DWORD PTR _.T521_$[esp+120]
	add	esp, 24					; 00000018H
	push	ecx
	push	8
	push	OFFSET FLAT:_7S
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _S, eax
	npad	1
	push	ecx
	push	1
	push	ecx
	push	1
	push	2
	push	OFFSET FLAT:_7H
	call	_make_info
	lea	ecx, DWORD PTR _.T521_$[esp+120]
	add	esp, 24					; 00000018H
	push	ecx
	push	8
	push	OFFSET FLAT:_7H
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _H, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7CO
	call	_make_info
	lea	ecx, DWORD PTR _.T521_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7CO
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _CO, eax
	lea	eax, DWORD PTR [ecx+ecx*2]
	add	eax, eax
	push	eax
	push	1
	push	1
	push	OFFSET FLAT:_7GRAD
	call	_make_info
	lea	ecx, DWORD PTR _.T521_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7GRAD
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _GRAD, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7T1
	call	_make_info
	lea	ecx, DWORD PTR _.T521_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7T1
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _T1, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7T2
	call	_make_info
	lea	ecx, DWORD PTR _.T521_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7T2
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _T2, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7T3
	call	_make_info
	lea	ecx, DWORD PTR _.T521_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7T3
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _M
	mov	DWORD PTR _T3, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7AA1
	call	_make_info
	lea	ecx, DWORD PTR _.T521_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7AA1
	call	__Allocate_i
	mov	ecx, DWORD PTR _.T521_$[esp+108]
	add	esp, 12					; 0000000cH
	mov	edx, DWORD PTR _M
	mov	DWORD PTR _AA1, eax
	dec	edx
	mov	DWORD PTR _ERROR, ecx

; 787  : 	-			AA2(M),RX2(M),RY2(M),RZ2(M),
; 788  : 	-			B12(M),EVEC(M),S(M,M),H(M,M),
; 789  :      -			CO(M),GRAD(M*6),T1(M),T2(M),
; 790  :      -			T3(M),STAT=ERROR)
; 791  : 
; 792  : 	DO III=1,M-1

	mov	ecx, 1
	test	edx, edx
	jle	$L522
$L98328:

; 793  : 		EVEC(III)=EVECSTORET(III)

	mov	ebp, DWORD PTR _7EVECSTORET+8
	mov	eax, DWORD PTR _EVECSTORET
	add	ebp, ecx
	mov	edi, DWORD PTR _7EVEC+8
	add	edi, ecx
	mov	ebx, DWORD PTR _EVEC
	mov	esi, DWORD PTR [eax+ebp*8+4]
	dec	edx
	mov	eax, DWORD PTR [eax+ebp*8]
	mov	DWORD PTR [ebx+edi*8+4], esi
	mov	DWORD PTR [ebx+edi*8], eax

; 794  : 	AA1(III)=TSTORE1T(III)

	mov	eax, DWORD PTR _TSTORE1T
	mov	esi, DWORD PTR _7TSTORE1T+8
	mov	ebx, DWORD PTR _7AA1+8
	add	esi, ecx
	add	ebx, ecx
	mov	edi, DWORD PTR _AA1
	npad	1
	mov	ebp, DWORD PTR [eax+esi*8+4]
	mov	eax, DWORD PTR [eax+esi*8]
	mov	DWORD PTR [edi+ebx*8+4], ebp
	mov	DWORD PTR [edi+ebx*8], eax

; 795  : 	AA2(III)=TSTORE2T(III)

	mov	eax, DWORD PTR _TSTORE2T
	mov	esi, DWORD PTR _7TSTORE2T+8
	mov	ebx, DWORD PTR _7AA2+8
	add	esi, ecx
	add	ebx, ecx
	mov	edi, DWORD PTR _AA2
	npad	1
	mov	ebp, DWORD PTR [eax+esi*8+4]
	mov	eax, DWORD PTR [eax+esi*8]
	mov	DWORD PTR [edi+ebx*8+4], ebp
	mov	DWORD PTR [edi+ebx*8], eax

; 796  : 	B12(III)=TSTORE3T(III)

	mov	eax, DWORD PTR _TSTORE3T
	mov	esi, DWORD PTR _7TSTORE3T+8
	mov	ebx, DWORD PTR _7B12+8
	add	esi, ecx
	add	ebx, ecx
	mov	edi, DWORD PTR _B12
	npad	1
	mov	ebp, DWORD PTR [eax+esi*8+4]
	mov	eax, DWORD PTR [eax+esi*8]
	mov	DWORD PTR [edi+ebx*8+4], ebp
	mov	DWORD PTR [edi+ebx*8], eax

; 797  : 	RZ1(III)=RZ1STORET(III)

	mov	eax, DWORD PTR _RZ1STORET
	mov	esi, DWORD PTR _7RZ1STORET+8
	mov	edi, DWORD PTR _7RZ1+8
	add	esi, ecx
	add	edi, ecx
	mov	ebx, DWORD PTR _RZ1
	npad	1
	mov	ebp, DWORD PTR [eax+esi*8+4]
	mov	eax, DWORD PTR [eax+esi*8]
	mov	DWORD PTR [ebx+edi*8+4], ebp
	mov	ebp, DWORD PTR _RZ2STORET
	mov	DWORD PTR [ebx+edi*8], eax

; 798  : 	RZ2(III)=RZ2STORET(III)

	mov	edi, DWORD PTR _7RZ2+8
	mov	esi, DWORD PTR _7RZ2STORET+8
	add	edi, ecx
	add	esi, ecx
	mov	ebx, DWORD PTR _RZ2
	mov	eax, DWORD PTR [ebp+esi*8+4]
	mov	ebp, DWORD PTR [ebp+esi*8]
	mov	DWORD PTR [ebx+edi*8+4], eax
	mov	eax, DWORD PTR _RX1
	mov	DWORD PTR [ebx+edi*8], ebp

; 799  : 	RX1(III)=0.

	mov	esi, DWORD PTR _7RX1+8
	add	esi, ecx
	mov	DWORD PTR [eax+esi*8], 0
	mov	DWORD PTR [eax+esi*8+4], 0

; 800  : 	RY1(III)=0.

	mov	esi, DWORD PTR _7RY1+8
	mov	eax, DWORD PTR _RY1
	add	esi, ecx
	mov	DWORD PTR [eax+esi*8], 0
	mov	DWORD PTR [eax+esi*8+4], 0

; 801  : 	RY2(III)=0.

	mov	esi, DWORD PTR _7RY2+8
	mov	eax, DWORD PTR _RY2
	add	esi, ecx
	mov	DWORD PTR [eax+esi*8], 0
	mov	DWORD PTR [eax+esi*8+4], 0

; 802  : 	RX2(III)=0.

	mov	esi, DWORD PTR _7RX2+8
	mov	eax, DWORD PTR _RX2
	add	esi, ecx
	inc	ecx
	test	edx, edx
	mov	DWORD PTR [eax+esi*8], 0
	mov	DWORD PTR [eax+esi*8+4], 0

; 787  : 	-			AA2(M),RX2(M),RY2(M),RZ2(M),
; 788  : 	-			B12(M),EVEC(M),S(M,M),H(M,M),
; 789  :      -			CO(M),GRAD(M*6),T1(M),T2(M),
; 790  :      -			T3(M),STAT=ERROR)
; 791  : 
; 792  : 	DO III=1,M-1

	jne	$L98328
$L522:

; 803  : 	END DO
; 804  : 	AA1(M)=TA1

	mov	ecx, DWORD PTR _M
	mov	eax, DWORD PTR _7AA1+8
	add	ecx, eax
	mov	edx, DWORD PTR _AA1
	mov	ebx, DWORD PTR _TA1+4
	mov	eax, DWORD PTR _TA1
	mov	DWORD PTR [edx+ecx*8+4], ebx
	mov	DWORD PTR [edx+ecx*8], eax

; 805  : 	AA2(M)=TA2

	mov	ecx, DWORD PTR _M
	mov	ebp, DWORD PTR _7AA2+8
	mov	edx, DWORD PTR _AA2
	add	ecx, ebp
	mov	eax, DWORD PTR _TA2+4
	mov	ebx, DWORD PTR _TA2
	npad	1
	mov	DWORD PTR [edx+ecx*8+4], eax
	mov	eax, DWORD PTR _TB12+4
	mov	DWORD PTR [edx+ecx*8], ebx

; 806  : 	B12(M)=TB12

	mov	edx, DWORD PTR _7B12+8
	mov	ecx, DWORD PTR _M
	mov	ebx, DWORD PTR _TB12
	add	edx, ecx
	mov	ecx, DWORD PTR _B12
	mov	DWORD PTR [ecx+edx*8+4], eax
	mov	eax, DWORD PTR _TRZ1+4
	mov	DWORD PTR [ecx+edx*8], ebx

; 807  : 	RZ1(M)=TRZ1

	mov	edx, DWORD PTR _M
	mov	esi, DWORD PTR _7RZ1+8
	mov	ecx, DWORD PTR _RZ1
	add	edx, esi
	mov	ebx, DWORD PTR _TRZ1
	mov	DWORD PTR [ecx+edx*8+4], eax
	mov	eax, DWORD PTR _TRZ2+4
	mov	DWORD PTR [ecx+edx*8], ebx

; 808  : 	RZ2(M)=TRZ2

	mov	ecx, DWORD PTR _7RZ2+8
	mov	edx, DWORD PTR _M
	mov	ebx, DWORD PTR _TRZ2
	add	ecx, edx
	mov	edx, DWORD PTR _RZ2
	mov	DWORD PTR [edx+ecx*8+4], eax
	mov	DWORD PTR [edx+ecx*8], ebx

; 809  : 	RX1(M)=0.

	mov	ecx, DWORD PTR _RX1
	mov	eax, DWORD PTR _7RX1+8
	mov	edi, DWORD PTR _M
	add	eax, edi
	mov	DWORD PTR [ecx+eax*8], 0
	mov	DWORD PTR [ecx+eax*8+4], 0

; 810  : 	RY1(M)=0.

	mov	eax, DWORD PTR _M
	mov	ecx, DWORD PTR _7RY1+8
	add	eax, ecx
	mov	ecx, DWORD PTR _RY1
	mov	DWORD PTR [ecx+eax*8], 0
	mov	DWORD PTR [ecx+eax*8+4], 0

; 811  : 	RY2(M)=0.

	mov	eax, DWORD PTR _7RY2+8
	mov	edi, DWORD PTR _M
	add	eax, edi
	mov	ecx, DWORD PTR _RY2
	mov	DWORD PTR [ecx+eax*8], 0
	mov	DWORD PTR [ecx+eax*8+4], 0

; 812  : 	RX2(M)=0.

	mov	eax, DWORD PTR _7RX2+8
	mov	edi, DWORD PTR _M
	add	eax, edi
	mov	ecx, DWORD PTR _RX2
	mov	DWORD PTR [ecx+eax*8], 0
	mov	edx, DWORD PTR _TSTORE1T
	mov	DWORD PTR [ecx+eax*8+4], 0

; 813  : 
; 814  : 
; 815  : 	DEALLOCATE(EVECSTORET,TSTORE1T,TSTORE2T,TSTORE3T,

	mov	DWORD PTR _.T524_$[esp+96], 0
	lea	eax, DWORD PTR _.T524_$[esp+96]
	push	eax
	push	8
	push	OFFSET FLAT:_7TSTORE1T
	push	edx
	call	__Deallocate_i
	lea	ecx, DWORD PTR _.T524_$[esp+112]
	add	esp, 16					; 00000010H
	mov	DWORD PTR _TSTORE1T, 0
	push	ecx
	mov	edx, DWORD PTR _TSTORE2T
	push	8
	push	OFFSET FLAT:_7TSTORE2T
	push	edx
	call	__Deallocate_i
	lea	ecx, DWORD PTR _.T524_$[esp+112]
	add	esp, 16					; 00000010H
	mov	DWORD PTR _TSTORE2T, 0
	push	ecx
	mov	edx, DWORD PTR _TSTORE3T
	push	8
	push	OFFSET FLAT:_7TSTORE3T
	push	edx
	call	__Deallocate_i
	lea	ecx, DWORD PTR _.T524_$[esp+112]
	add	esp, 16					; 00000010H
	mov	DWORD PTR _TSTORE3T, 0
	push	ecx
	mov	edx, DWORD PTR _RZ1STORET
	push	8
	push	OFFSET FLAT:_7RZ1STORET
	push	edx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	DWORD PTR _RZ1STORET, 0
	lea	ecx, DWORD PTR _.T524_$[esp+96]
	mov	edx, DWORD PTR _RZ2STORET
	push	ecx
	push	8
	push	OFFSET FLAT:_7RZ2STORET
	push	edx
	call	__Deallocate_i
	lea	ecx, DWORD PTR _.T524_$[esp+112]
	add	esp, 16					; 00000010H
	mov	DWORD PTR _RZ2STORET, 0
	push	ecx
	mov	edx, DWORD PTR _EVECSTORET
	push	8
	push	OFFSET FLAT:_7EVECSTORET
	push	edx
	call	__Deallocate_i
	mov	ecx, DWORD PTR _.T524_$[esp+112]
	add	esp, 16					; 00000010H
	mov	edx, DWORD PTR _T3
	mov	DWORD PTR _ERROR, ecx
	push	OFFSET FLAT:_DFLAG
	mov	eax, DWORD PTR _T2
	push	edx
	mov	ecx, DWORD PTR _T1
	push	eax
	mov	edx, DWORD PTR _GRAD
	push	ecx
	mov	eax, DWORD PTR _CO
	push	edx
	mov	ecx, DWORD PTR _H
	push	OFFSET FLAT:_EMIN
	mov	edx, DWORD PTR _S
	mov	DWORD PTR _EVECSTORET, 0

; 816  :      -			RZ1STORET,RZ2STORET,STAT=ERROR)
; 817  : C
; 818  : 	DFLAG=0

	mov	DWORD PTR _DFLAG, 0

; 819  : 	CALL HAMILTON(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,EVEC,S,H,

	push	eax
	push	ecx
	push	edx
	mov	eax, DWORD PTR _EVEC
	push	eax
	mov	ecx, DWORD PTR _B12
	push	ecx
	mov	edx, DWORD PTR _RZ2
	push	edx
	mov	eax, DWORD PTR _RY2
	push	eax
	mov	ecx, DWORD PTR _RX2
	push	ecx
	mov	edx, DWORD PTR _AA2
	push	edx
	mov	eax, DWORD PTR _RZ1
	push	eax
	mov	ecx, DWORD PTR _RY1
	push	ecx
	mov	edx, DWORD PTR _RX1
	push	edx
	mov	eax, DWORD PTR _AA1
	push	eax
	push	OFFSET FLAT:_M
	call	_HAMILTON@80

; 820  : 	-				CO,EMIN,GRAD,T1,T2,T3,DFLAG)
; 821  : 	DO NNN=1,M

	mov	esi, DWORD PTR _M
	mov	edx, 1
	test	esi, esi
	jle	SHORT $L525
$L98329:

; 822  : 		EVEC(NNN)=CO(NNN)

	mov	ebp, DWORD PTR _7CO+8
	mov	eax, DWORD PTR _CO
	add	ebp, edx
	mov	edi, DWORD PTR _7EVEC+8
	add	edi, edx
	mov	ebx, DWORD PTR _EVEC
	mov	ecx, DWORD PTR [eax+ebp*8+4]
	dec	esi
	mov	eax, DWORD PTR [eax+ebp*8]
	inc	edx
	mov	DWORD PTR [ebx+edi*8+4], ecx
	test	esi, esi
	mov	DWORD PTR [ebx+edi*8], eax

; 820  : 	-				CO,EMIN,GRAD,T1,T2,T3,DFLAG)
; 821  : 	DO NNN=1,M

	jne	SHORT $L98329
$L525:

; 823  : 	END DO
; 824  : 	WRITE(6,*)'NOW WE HAVE',M,' FUNCTIONS...'

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	11					; 0000000bH
	push	OFFSET FLAT:_.SG527_
	call	_lwrite_ch
	add	esp, 8
	mov	eax, DWORD PTR _M
	push	eax
	call	_lwrite_i
	add	esp, 4
	push	13					; 0000000dH
	push	OFFSET FLAT:_.SG528_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 825  : 	WRITE(6,*)'EXACT ENERGY=',EMIN

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	13					; 0000000dH
	push	OFFSET FLAT:_.SG529_
	call	_lwrite_ch
	add	esp, 8
	mov	eax, DWORD PTR _EMIN+4
	mov	ecx, DWORD PTR _EMIN
	npad	1
	push	eax
	push	ecx
	call	_lwrite_d
	add	esp, 8
	call	_f90_endwrite

; 624  : 	DO KKK=1,ITEMP

	dec	DWORD PTR -36+[esp+96]
	jne	$L98325
$L462:

; 826  : 	END DO
; 827  : 	WRITE(6,*)'********************************'

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	32					; 00000020H
	push	OFFSET FLAT:_.SG530_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 828  : 	WRITE(6,*)'PERFORMING GRADIENT OPTIMIZATION'

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	32					; 00000020H
	push	OFFSET FLAT:_.SG531_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 829  : 	WRITE(6,*)'********************************'

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	32					; 00000020H
	push	OFFSET FLAT:_.SG532_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 830  : 	DO NNN=1,M

	mov	esi, DWORD PTR _M
	mov	edx, 1
	test	esi, esi
	jle	SHORT $L533
$L98330:

; 831  : 		EVEC(NNN)=CO(NNN)

	mov	ecx, DWORD PTR _7CO+8
	mov	ebp, DWORD PTR _CO
	add	ecx, edx
	mov	edi, DWORD PTR _7EVEC+8
	add	edi, edx
	mov	ebx, DWORD PTR _EVEC
	mov	eax, DWORD PTR [ebp+ecx*8+4]
	dec	esi
	mov	ecx, DWORD PTR [ebp+ecx*8]
	inc	edx
	mov	DWORD PTR [ebx+edi*8+4], eax
	test	esi, esi
	mov	DWORD PTR [ebx+edi*8], ecx

; 830  : 	DO NNN=1,M

	jne	SHORT $L98330
$L533:

; 832  : 	END DO
; 833  : 	DO NNN=1,M

	mov	esi, DWORD PTR _M
	mov	edi, 1
	test	esi, esi
	jle	$L535
	jmp	SHORT $L98331
$L98469:
	fstp	QWORD PTR [eax+ecx*8]
$L98331:

; 834  : 		T1(NNN)=DSQRT(AA1(NNN)+B12(NNN))

	mov	eax, DWORD PTR _7B12+8
	mov	ecx, DWORD PTR _B12
	add	eax, edi
	dec	esi
	fld	QWORD PTR [ecx+eax*8]
	mov	eax, DWORD PTR _7AA1+8
	mov	ecx, DWORD PTR _AA1
	add	eax, edi
	fadd	QWORD PTR [ecx+eax*8]
	call	__FIsqrt
	mov	ecx, DWORD PTR _7T1+8
	mov	eax, DWORD PTR _T1
	add	ecx, edi
	mov	edx, DWORD PTR _B12
	fstp	QWORD PTR [eax+ecx*8]
	mov	ebx, DWORD PTR _7B12+8
	add	ebx, edi
	fld	QWORD PTR [edx+ebx*8]
	fld	QWORD PTR $T98441
	lea	ebx, DWORD PTR [edx+ebx*8]

; 835  : 		T2(NNN)=DSQRT(AA2(NNN)+B12(NNN)- ((B12(NNN)**2)/(AA1(NNN)

	call	__FIfexp
	mov	eax, DWORD PTR _7AA1+8
	mov	ecx, DWORD PTR _AA1
	add	eax, edi
	fld	QWORD PTR [ecx+eax*8]
	fadd	QWORD PTR [ebx]
	mov	ecx, DWORD PTR _7AA2+8
	mov	eax, DWORD PTR _AA2
	add	ecx, edi
	fdivp	ST(1), ST(0)
	fsubr	QWORD PTR [eax+ecx*8]
	fadd	QWORD PTR [ebx]
	call	__FIsqrt
	mov	ecx, DWORD PTR _7T2+8
	mov	eax, DWORD PTR _T2
	add	ecx, edi
	mov	edx, DWORD PTR _B12
	fstp	QWORD PTR [eax+ecx*8]
	mov	eax, DWORD PTR _7B12+8
	mov	ecx, DWORD PTR _7AA1+8
	add	eax, edi
	add	ecx, edi
	fld	QWORD PTR [edx+eax*8]
	lea	eax, DWORD PTR [edx+eax*8]

; 836  : 	-		+B12(NNN))))
; 837  : 		T3(NNN)=-B12(NNN)/DSQRT(AA1(NNN)+B12(NNN))

	mov	edx, DWORD PTR _AA1
	fld	QWORD PTR [edx+ecx*8]
	fadd	QWORD PTR [eax]
	call	__FIsqrt
	fdivp	ST(1), ST(0)
	mov	eax, DWORD PTR _T3
	mov	ecx, DWORD PTR _7T3+8
	add	ecx, edi
	inc	edi
	test	esi, esi
	fchs

; 832  : 	END DO
; 833  : 	DO NNN=1,M

	jne	$L98469
	fstp	QWORD PTR [eax+ecx*8]
$L535:

; 838  : C	WRITE(6,*)'A1 ->',AA1(NNN),T1(NNN)**2 +T1(NNN)*T3(NNN)
; 839  : C	WRITE(6,*)'A2 ->',AA2(NNN),T2(NNN)**2 +T3(NNN)**2 +T1(NNN)*T3(NNN)
; 840  : C	WRITE(6,*)'B12->',B12(NNN),-T1(NNN)*T3(NNN)
; 841  : C	WRITE(6,*)'-------------------------'
; 842  : 	END DO
; 843  : 	DEFLAG=1

	push	OFFSET FLAT:_DFLAG
	mov	eax, DWORD PTR _T3
	push	eax
	mov	ecx, DWORD PTR _T2
	push	ecx
	mov	edx, DWORD PTR _T1
	push	edx
	mov	eax, DWORD PTR _GRAD
	push	eax
	mov	ecx, DWORD PTR _CO
	push	OFFSET FLAT:_EMIN
	mov	edx, DWORD PTR _H
	push	ecx
	mov	eax, DWORD PTR _S
	push	edx
	mov	ecx, DWORD PTR _EVEC
	push	eax
	mov	edx, DWORD PTR _B12
	push	ecx
	mov	eax, DWORD PTR _RZ2
	push	edx
	mov	ecx, DWORD PTR _RY2
	push	eax
	mov	edx, DWORD PTR _RX2
	push	ecx
	mov	eax, DWORD PTR _AA2
	push	edx
	mov	ecx, DWORD PTR _RZ1
	push	eax
	mov	edx, DWORD PTR _RY1
	push	ecx
	mov	eax, DWORD PTR _RX1
	push	edx
	mov	ecx, DWORD PTR _AA1
	mov	DWORD PTR _DEFLAG, 0
	mov	DWORD PTR _DEFLAG+4, 1072693248		; 3ff00000H

; 844  : 	CALL HAMILTON(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,EVEC,S,H,

	push	eax
	push	ecx
	push	OFFSET FLAT:_M
	call	_HAMILTON@80

; 845  : 	-				CO,EMIN,GRAD,T1,T2,T3,DFLAG)
; 846  : 	WRITE(6,*)'PRE-OPTIMIZATION ENERGY --->',EMIN

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	28					; 0000001cH
	push	OFFSET FLAT:_.SG537_
	call	_lwrite_ch
	add	esp, 8
	mov	ecx, DWORD PTR _EMIN+4
	mov	edx, DWORD PTR _EMIN
	npad	1
	push	ecx
	push	edx
	call	_lwrite_d
	add	esp, 8
	call	_f90_endwrite

; 847  : 	CALL GRADOPT(M,T1,RX1,RY1,RZ1,T2,RX2,RY2,RZ2,T3,EVEC,GRAD)

	mov	ecx, DWORD PTR _GRAD
	mov	edx, DWORD PTR _EVEC
	push	ecx
	mov	eax, DWORD PTR _T3
	push	edx
	mov	ecx, DWORD PTR _RZ2
	push	eax
	mov	edx, DWORD PTR _RY2
	push	ecx
	mov	eax, DWORD PTR _RX2
	push	edx
	mov	ecx, DWORD PTR _T2
	push	eax
	mov	edx, DWORD PTR _RZ1
	push	ecx
	mov	eax, DWORD PTR _RY1
	push	edx
	mov	ecx, DWORD PTR _RX1
	push	eax
	mov	edx, DWORD PTR _T1
	push	ecx
	mov	ebp, 1
	push	edx
	push	OFFSET FLAT:_M
	call	_GRADOPT@48

; 848  : 	DO NNN=1,M

	mov	esi, DWORD PTR _M
	mov	ebx, DWORD PTR _T1
	mov	edi, DWORD PTR _T3
	test	esi, esi
	jle	$L538
	jmp	SHORT $L98332
$L98470:
	fstp	QWORD PTR [eax+ecx*8]
$L98332:
	mov	eax, DWORD PTR _7T1+8
	mov	ecx, DWORD PTR _7T3+8
	add	eax, ebp
	add	ecx, ebp
	dec	esi
	fld	QWORD PTR [edi+ecx*8]
	fadd	QWORD PTR [ebx+eax*8]
	lea	eax, DWORD PTR [ebx+eax*8]
	mov	ecx, DWORD PTR _AA1

; 849  : 		AA1(NNN)=T1(NNN)**2 +T1(NNN)*T3(NNN)

	fmul	QWORD PTR [eax]
	mov	eax, DWORD PTR _7AA1+8
	add	eax, ebp
	fstp	QWORD PTR [ecx+eax*8]
	mov	eax, DWORD PTR _7T3+8
	mov	ecx, DWORD PTR _7T1+8
	add	eax, ebp
	add	ecx, ebp
	fld	QWORD PTR [ebx+ecx*8]
	fadd	QWORD PTR [edi+eax*8]
	lea	eax, DWORD PTR [edi+eax*8]
	mov	ecx, DWORD PTR _7T2+8

; 850  : 		AA2(NNN)=T2(NNN)**2 + T3(NNN)**2 +T1(NNN)*T3(NNN)

	add	ecx, ebp
	fmul	QWORD PTR [eax]
	mov	eax, DWORD PTR _T2
	fld	QWORD PTR [eax+ecx*8]
	fld	QWORD PTR $T98441
	call	__FIfexp
	faddp	ST(1), ST(0)
	mov	ecx, DWORD PTR _7AA2+8
	mov	eax, DWORD PTR _AA2
	add	ecx, ebp
	fstp	QWORD PTR [eax+ecx*8]

; 851  : 		B12(NNN)=-T1(NNN)*T3(NNN)

	mov	eax, DWORD PTR _7T1+8
	mov	ecx, DWORD PTR _7B12+8
	add	eax, ebp
	add	ecx, ebp
	fld	QWORD PTR [ebx+eax*8]
	mov	eax, DWORD PTR _7T3+8
	add	eax, ebp
	inc	ebp
	test	esi, esi
	fmul	QWORD PTR [edi+eax*8]
	mov	eax, DWORD PTR _B12
	fchs

; 848  : 	DO NNN=1,M

	jne	$L98470
	fstp	QWORD PTR [eax+ecx*8]
$L538:

; 852  : 	END DO
; 853  : C
; 854  : 	CALL HAMILTON(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,EVEC,S,H,

	push	OFFSET FLAT:_DFLAG
	mov	eax, DWORD PTR _T2
	push	edi
	mov	ecx, DWORD PTR _GRAD
	push	eax
	mov	edx, DWORD PTR _CO
	push	ebx
	mov	eax, DWORD PTR _H
	push	ecx
	push	OFFSET FLAT:_EMIN
	push	edx
	mov	ecx, DWORD PTR _S
	push	eax
	mov	edx, DWORD PTR _EVEC
	push	ecx
	mov	eax, DWORD PTR _B12
	push	edx
	mov	ecx, DWORD PTR _RZ2
	push	eax
	mov	edx, DWORD PTR _RY2
	push	ecx
	mov	eax, DWORD PTR _RX2
	push	edx
	mov	ecx, DWORD PTR _AA2
	push	eax
	mov	edx, DWORD PTR _RZ1
	push	ecx
	mov	eax, DWORD PTR _RY1
	push	edx
	mov	ecx, DWORD PTR _RX1
	push	eax
	mov	edx, DWORD PTR _AA1
	push	ecx
	mov	esi, 1
	push	edx
	push	OFFSET FLAT:_M
	call	_HAMILTON@80

; 855  : 	-				CO,EMIN,GRAD,T1,T2,T3,DFLAG)
; 856  : 	DO J=1,M

	mov	edx, DWORD PTR _M
	test	edx, edx
	jle	SHORT $L98387
$L98333:

; 857  : 		EVEC(J)=CO(J)

	mov	ebp, DWORD PTR _7CO+8
	mov	eax, DWORD PTR _CO
	add	ebp, esi
	mov	edi, DWORD PTR _7EVEC+8
	add	edi, esi
	mov	ebx, DWORD PTR _EVEC
	mov	ecx, DWORD PTR [eax+ebp*8+4]
	dec	edx
	mov	eax, DWORD PTR [eax+ebp*8]
	inc	esi
	mov	DWORD PTR [ebx+edi*8+4], ecx
	test	edx, edx
	mov	DWORD PTR [ebx+edi*8], eax

; 855  : 	-				CO,EMIN,GRAD,T1,T2,T3,DFLAG)
; 856  : 	DO J=1,M

	jne	SHORT $L98333
$L98387:
	mov	eax, DWORD PTR _ISTEND
	mov	ecx, DWORD PTR _M
	cmp	ecx, eax
	jle	$L455
$L449:

; 858  : 	END DO
; 859  : C
; 860  : 	END DO
; 861  : C
; 862  : C
; 863  : 	END IF
; 864  : 
; 865  : C
; 866  : C	
; 867  : 	CLOSE(7)

	push	0
	push	0
	push	0
	push	7
	call	_f90_close
	add	esp, 16					; 00000010H

; 868  : 	CLOSE(8)

	push	0
	push	0
	push	0
	push	8
	call	_f90_close
	add	esp, 16					; 00000010H

; 869  : 	CLOSE(9)

	push	0
	push	0
	push	0
	push	9
	call	_f90_close
	add	esp, 16					; 00000010H

; 870  : 	STOP

	push	0
	push	0
	call	__FFsSTOP
	add	esp, 8

; 871  : C
; 872  :       CALL GETDAT(tmpyear, tmpmonth, tmpday)

	push	OFFSET FLAT:_TMPDAY
	push	OFFSET FLAT:_TMPMONTH
	push	OFFSET FLAT:_TMPYEAR
	call	_GETDAT@12

; 873  :       CALL GETTIM(tmphour, tmpminute, tmpsecond, tmphund)  

	push	OFFSET FLAT:_TMPHUND
	push	OFFSET FLAT:_TMPSECOND
	push	OFFSET FLAT:_TMPMINUTE
	push	OFFSET FLAT:_TMPHOUR
	call	_GETTIM@16

; 874  : C
; 875  :       WRITE (6, 9003) tmpmonth, tmpday, tmpyear

	push	0
	push	0
	push	0
	push	-1
	push	OFFSET FLAT:_.SG542_
	push	6
	call	_f90_write
	mov	ax, WORD PTR _TMPMONTH
	add	esp, 24					; 00000018H
	push	eax
	call	_write_i2
	mov	ax, WORD PTR _TMPDAY
	add	esp, 4
	push	eax
	call	_write_i2
	mov	ax, WORD PTR _TMPYEAR
	add	esp, 4
	push	eax
	call	_write_i2
	add	esp, 4
	call	_f90_endwrite

; 876  :       WRITE (6, 9004) tmphour, tmpminute, tmpsecond, tmphund

	push	0
	push	0
	push	0
	push	-1
	push	OFFSET FLAT:_.SG543_
	push	6
	call	_f90_write
	mov	ax, WORD PTR _TMPHOUR
	add	esp, 24					; 00000018H
	push	eax
	call	_write_i2
	mov	ax, WORD PTR _TMPMINUTE
	add	esp, 4
	push	eax
	call	_write_i2
	mov	ax, WORD PTR _TMPSECOND
	add	esp, 4
	push	eax
	call	_write_i2
	mov	ax, WORD PTR _TMPHUND
	add	esp, 4
	push	eax
	call	_write_i2
	add	esp, 4
	call	_f90_endwrite

; 877  : C
; 878  :       WRITE (7, 9003) tmpmonth, tmpday, tmpyear

	push	0
	push	0
	push	0
	push	-1
	push	OFFSET FLAT:_.SG544_
	push	7
	call	_f90_write
	mov	ax, WORD PTR _TMPMONTH
	add	esp, 24					; 00000018H
	push	eax
	call	_write_i2
	mov	ax, WORD PTR _TMPDAY
	add	esp, 4
	push	eax
	call	_write_i2
	mov	ax, WORD PTR _TMPYEAR
	add	esp, 4
	push	eax
	call	_write_i2
	add	esp, 4
	call	_f90_endwrite

; 879  :       WRITE (7, 9004) tmphour, tmpminute, tmpsecond, tmphund

	push	0
	push	0
	push	0
	push	-1
	push	OFFSET FLAT:_.SG545_
	push	7
	call	_f90_write
	mov	ax, WORD PTR _TMPHOUR
	add	esp, 24					; 00000018H
	push	eax
	call	_write_i2
	mov	ax, WORD PTR _TMPMINUTE
	add	esp, 4
	push	eax
	call	_write_i2
	mov	ax, WORD PTR _TMPSECOND
	add	esp, 4
	push	eax
	call	_write_i2
	mov	ax, WORD PTR _TMPHUND
	add	esp, 4
	push	eax
	call	_write_i2
	add	esp, 4
	call	_f90_endwrite

; 880  : C
; 881  :       CALL BIOUT (' ')

	push	1
	push	OFFSET FLAT:_.SG546_
	call	_BIOUT@8

; 882  : C     -----------------------------------------------------
; 883  : C	OPTIMIZATION ROUTINE BEGINS...
; 884  : C     -----------------------------------------------------
; 885  : 	IF (WAVEOPT.EQ.0) THEN

	mov	eax, DWORD PTR _WAVEOPT
	test	eax, eax
	jne	$L98465

; 886  : 	IF (METHODOPT.EQ.0) THEN

	mov	eax, DWORD PTR _METHODOPT
	test	eax, eax
	jne	SHORT $L98464

; 887  : 	CALL NGRADOPT(M,T1,RX1,RY1,RZ1,T2,RX2,RY2,RZ2,T3,EVEC,GRAD)

	mov	eax, DWORD PTR _GRAD
	mov	ecx, DWORD PTR _EVEC
	push	eax
	mov	edi, DWORD PTR _T3
	push	ecx
	mov	eax, DWORD PTR _RZ2
	push	edi
	mov	ecx, DWORD PTR _RY2
	push	eax
	mov	edx, DWORD PTR _RX2
	push	ecx
	mov	eax, DWORD PTR _T2
	push	edx
	mov	ecx, DWORD PTR _RZ1
	push	eax
	mov	edx, DWORD PTR _RY1
	push	ecx
	mov	eax, DWORD PTR _RX1
	push	edx
	mov	esi, DWORD PTR _T1
	push	eax
	push	esi
	push	OFFSET FLAT:_M
	call	_NGRADOPT@48
$L98464:
	mov	esi, DWORD PTR _T1
	mov	edi, DWORD PTR _T3

; 888  : 	END IF
; 889  : 	IF (METHODOPT.EQ.1) THEN

	mov	eax, DWORD PTR _METHODOPT
	cmp	eax, 1
	jne	SHORT $L547

; 890  : 	CALL GRADOPT(M,T1,RX1,RY1,RZ1,T2,RX2,RY2,RZ2,T3,EVEC,GRAD)

	mov	eax, DWORD PTR _GRAD
	mov	ecx, DWORD PTR _EVEC
	push	eax
	mov	edx, DWORD PTR _RZ2
	push	ecx
	mov	eax, DWORD PTR _RY2
	push	edi
	mov	ecx, DWORD PTR _RX2
	push	edx
	push	eax
	push	ecx
	mov	edx, DWORD PTR _T2
	push	edx
	mov	eax, DWORD PTR _RZ1
	push	eax
	mov	ecx, DWORD PTR _RY1
	push	ecx
	mov	edx, DWORD PTR _RX1
	push	edx
	push	esi
	push	OFFSET FLAT:_M
	call	_GRADOPT@48
$L98465:
	mov	esi, DWORD PTR _T1
	mov	edi, DWORD PTR _T3
$L547:

; 891  : 	END IF
; 892  : 	END IF
; 893  : C
; 894  : 	IF (WAVEOPT.EQ.1) THEN

	mov	eax, DWORD PTR _WAVEOPT
	cmp	eax, 1
	jne	$L555

; 895  : 	RELOADFLAG=0

	mov	eax, DWORD PTR _LOOP
	mov	ecx, DWORD PTR _ENDBASIS
	mov	DWORD PTR _RELOADFLAG, 0

; 896  : 	ALLOCATE(EVECSTORE(ENDBASIS,LOOP),TSTORE1(ENDBASIS,LOOP),

	mov	DWORD PTR _.T558_$[esp+96], 0
	push	eax
	push	1
	push	ecx
	push	1
	push	2
	push	OFFSET FLAT:_7TSTORE1
	call	_make_info
	lea	ecx, DWORD PTR _.T558_$[esp+120]
	add	esp, 24					; 00000018H
	push	ecx
	push	8
	push	OFFSET FLAT:_7TSTORE1
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _LOOP
	mov	edx, DWORD PTR _ENDBASIS
	mov	DWORD PTR _TSTORE1, eax
	push	ecx
	push	1
	push	edx
	push	1
	push	2
	push	OFFSET FLAT:_7TSTORE2
	call	_make_info
	lea	ecx, DWORD PTR _.T558_$[esp+120]
	add	esp, 24					; 00000018H
	push	ecx
	push	8
	push	OFFSET FLAT:_7TSTORE2
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _LOOP
	mov	edx, DWORD PTR _ENDBASIS
	mov	DWORD PTR _TSTORE2, eax
	push	ecx
	push	1
	push	edx
	push	1
	push	2
	push	OFFSET FLAT:_7TSTORE3
	call	_make_info
	lea	ecx, DWORD PTR _.T558_$[esp+120]
	add	esp, 24					; 00000018H
	push	ecx
	push	8
	push	OFFSET FLAT:_7TSTORE3
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _LOOP
	mov	edx, DWORD PTR _ENDBASIS
	mov	DWORD PTR _TSTORE3, eax
	push	ecx
	push	1
	push	edx
	push	1
	push	2
	push	OFFSET FLAT:_7RZ1STORE
	call	_make_info
	lea	ecx, DWORD PTR _.T558_$[esp+120]
	add	esp, 24					; 00000018H
	push	ecx
	push	8
	push	OFFSET FLAT:_7RZ1STORE
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _LOOP
	mov	edx, DWORD PTR _ENDBASIS
	mov	DWORD PTR _RZ1STORE, eax
	push	ecx
	push	1
	push	edx
	push	1
	push	2
	push	OFFSET FLAT:_7RZ2STORE
	call	_make_info
	lea	ecx, DWORD PTR _.T558_$[esp+120]
	add	esp, 24					; 00000018H
	push	ecx
	push	8
	push	OFFSET FLAT:_7RZ2STORE
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _LOOP
	mov	edx, DWORD PTR _ENDBASIS
	mov	DWORD PTR _RZ2STORE, eax
	push	ecx
	push	1
	push	edx
	push	1
	push	2
	push	OFFSET FLAT:_7EVECSTORE
	call	_make_info
	lea	ecx, DWORD PTR _.T558_$[esp+120]
	add	esp, 24					; 00000018H
	push	ecx
	push	8
	push	OFFSET FLAT:_7EVECSTORE
	call	__Allocate_i
	mov	ecx, DWORD PTR _.T558_$[esp+108]
	mov	edx, DWORD PTR _ENDBASIS
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _EVECSTORE, eax
	mov	DWORD PTR _ERROR, ecx

; 897  :      -			TSTORE2(ENDBASIS,LOOP),TSTORE3(ENDBASIS,LOOP),
; 898  :      -			RZ1STORE(ENDBASIS,LOOP),RZ2STORE(ENDBASIS,LOOP),
; 899  :      -			STAT=ERROR)
; 900  : 	ALLOCATE(EVECSTORET(ENDBASIS),TSTORE1T(ENDBASIS),

	mov	DWORD PTR _.T566_$[esp+96], 0
	push	edx
	push	1
	push	1
	push	OFFSET FLAT:_7TSTORE1T
	call	_make_info
	lea	ecx, DWORD PTR _.T566_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7TSTORE1T
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _ENDBASIS
	mov	DWORD PTR _TSTORE1T, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7TSTORE2T
	call	_make_info
	lea	ecx, DWORD PTR _.T566_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7TSTORE2T
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _ENDBASIS
	mov	DWORD PTR _TSTORE2T, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7TSTORE3T
	call	_make_info
	lea	ecx, DWORD PTR _.T566_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7TSTORE3T
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _ENDBASIS
	mov	DWORD PTR _TSTORE3T, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7RZ1STORET
	call	_make_info
	lea	ecx, DWORD PTR _.T566_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7RZ1STORET
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _ENDBASIS
	mov	DWORD PTR _RZ1STORET, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7RZ2STORET
	call	_make_info
	lea	ecx, DWORD PTR _.T566_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7RZ2STORET
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _ENDBASIS
	mov	DWORD PTR _RZ2STORET, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7EVECSTORET
	call	_make_info
	lea	ecx, DWORD PTR _.T566_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7EVECSTORET
	call	__Allocate_i
	mov	ecx, DWORD PTR _.T566_$[esp+108]
	add	esp, 12					; 0000000cH
	mov	edx, DWORD PTR _STARTBASIS
	mov	DWORD PTR _EVECSTORET, eax
	mov	eax, DWORD PTR _ENDBASIS
	mov	DWORD PTR _ERROR, ecx

; 901  :      -			TSTORE2T(ENDBASIS),TSTORE3T(ENDBASIS),
; 902  :      -			RZ1STORET(ENDBASIS),RZ2STORET(ENDBASIS),
; 903  :      -			STAT=ERROR)
; 904  : 
; 905  : 		DO J=STARTBASIS,ENDBASIS

	sub	eax, edx
	mov	esi, DWORD PTR _T1
	inc	eax
	mov	edi, DWORD PTR _T3
	mov	DWORD PTR -8+[esp+96], eax
	mov	DWORD PTR _J$[esp+96], edx
	test	eax, eax
	jle	$L555
$L98334:

; 906  : 		M=J

	mov	eax, DWORD PTR _J$[esp+96]
	mov	ecx, DWORD PTR _ERROR
	mov	DWORD PTR _M, eax
	mov	DWORD PTR _MMAX, eax

; 907  : 		MMAX=M
; 908  : C
; 909  : 	IF (ERROR .NE. 0) STOP '!NOT ENOUGH MEMORY SPACE AVAILABLE ON THIS

	test	ecx, ecx
	je	SHORT $L570
	push	OFFSET FLAT:_.SG572_
	push	60					; 0000003cH
	call	__FFsSTOP
	add	esp, 8
$L570:

; 910  : 	- BOX! ...ABORTING!'
; 911  : C
; 912  : 	DO K=1,LOOP

	mov	DWORD PTR _K$[esp+96], 1
	mov	eax, DWORD PTR _LOOP
	mov	DWORD PTR -24+[esp+96], eax
	test	eax, eax
	jle	$L573
$L98335:

; 913  : C
; 914  : C
; 915  : 	DEALLOCATE (AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,EVEC,S,H,CO,

	push	0
	mov	eax, DWORD PTR _RX1
	push	8
	push	OFFSET FLAT:_7RX1
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _RY1
	mov	DWORD PTR _RX1, 0
	push	0
	push	8
	push	OFFSET FLAT:_7RY1
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _RZ1
	mov	DWORD PTR _RY1, 0
	push	0
	push	8
	push	OFFSET FLAT:_7RZ1
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _AA2
	mov	DWORD PTR _RZ1, 0
	push	0
	push	8
	push	OFFSET FLAT:_7AA2
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _RX2
	mov	DWORD PTR _AA2, 0
	push	0
	push	8
	push	OFFSET FLAT:_7RX2
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _RY2
	mov	DWORD PTR _RX2, 0
	push	0
	push	8
	push	OFFSET FLAT:_7RY2
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _RZ2
	mov	DWORD PTR _RY2, 0
	push	0
	push	8
	push	OFFSET FLAT:_7RZ2
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _B12
	mov	DWORD PTR _RZ2, 0
	push	0
	push	8
	push	OFFSET FLAT:_7B12
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _EVEC
	mov	DWORD PTR _B12, 0
	push	0
	push	8
	push	OFFSET FLAT:_7EVEC
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _S
	mov	DWORD PTR _EVEC, 0
	push	0
	push	8
	push	OFFSET FLAT:_7S
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _H
	mov	DWORD PTR _S, 0
	push	0
	push	8
	push	OFFSET FLAT:_7H
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _CO
	mov	DWORD PTR _H, 0
	push	0
	push	8
	push	OFFSET FLAT:_7CO
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _GRAD
	mov	DWORD PTR _CO, 0
	push	0
	push	8
	push	OFFSET FLAT:_7GRAD
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _T1
	mov	DWORD PTR _GRAD, 0
	push	0
	push	8
	push	OFFSET FLAT:_7T1
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _T2
	mov	DWORD PTR _T1, 0
	push	0
	push	8
	push	OFFSET FLAT:_7T2
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _T3
	mov	DWORD PTR _T2, 0
	push	0
	push	8
	push	OFFSET FLAT:_7T3
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _AA1
	mov	DWORD PTR _T3, 0
	push	0
	push	8
	push	OFFSET FLAT:_7AA1
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _MMAX
	mov	DWORD PTR _AA1, 0

; 916  : 	-			GRAD,T1,T2,T3)
; 917  :       ALLOCATE(	AA1(MMAX),RX1(MMAX),RY1(MMAX),RZ1(MMAX),

	mov	DWORD PTR _.T575_$[esp+96], 0
	push	eax
	push	1
	push	1
	push	OFFSET FLAT:_7RX1
	call	_make_info
	lea	eax, DWORD PTR _.T575_$[esp+112]
	add	esp, 16					; 00000010H
	push	eax
	push	8
	push	OFFSET FLAT:_7RX1
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _RX1, eax
	mov	eax, DWORD PTR _MMAX
	push	eax
	push	1
	push	1
	push	OFFSET FLAT:_7RY1
	call	_make_info
	lea	eax, DWORD PTR _.T575_$[esp+112]
	add	esp, 16					; 00000010H
	push	eax
	push	8
	push	OFFSET FLAT:_7RY1
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _RY1, eax
	mov	eax, DWORD PTR _MMAX
	push	eax
	push	1
	push	1
	push	OFFSET FLAT:_7RZ1
	call	_make_info
	lea	eax, DWORD PTR _.T575_$[esp+112]
	add	esp, 16					; 00000010H
	push	eax
	push	8
	push	OFFSET FLAT:_7RZ1
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _RZ1, eax
	mov	eax, DWORD PTR _MMAX
	push	eax
	push	1
	push	1
	push	OFFSET FLAT:_7AA2
	call	_make_info
	lea	eax, DWORD PTR _.T575_$[esp+112]
	add	esp, 16					; 00000010H
	push	eax
	push	8
	push	OFFSET FLAT:_7AA2
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _AA2, eax
	mov	eax, DWORD PTR _MMAX
	push	eax
	push	1
	push	1
	push	OFFSET FLAT:_7RX2
	call	_make_info
	lea	eax, DWORD PTR _.T575_$[esp+112]
	add	esp, 16					; 00000010H
	push	eax
	push	8
	push	OFFSET FLAT:_7RX2
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _RX2, eax
	mov	eax, DWORD PTR _MMAX
	push	eax
	push	1
	push	1
	push	OFFSET FLAT:_7RY2
	call	_make_info
	lea	eax, DWORD PTR _.T575_$[esp+112]
	add	esp, 16					; 00000010H
	push	eax
	push	8
	push	OFFSET FLAT:_7RY2
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _RY2, eax
	mov	eax, DWORD PTR _MMAX
	push	eax
	push	1
	push	1
	push	OFFSET FLAT:_7RZ2
	call	_make_info
	lea	eax, DWORD PTR _.T575_$[esp+112]
	add	esp, 16					; 00000010H
	push	eax
	push	8
	push	OFFSET FLAT:_7RZ2
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _RZ2, eax
	mov	eax, DWORD PTR _MMAX
	push	eax
	push	1
	push	1
	push	OFFSET FLAT:_7B12
	call	_make_info
	lea	eax, DWORD PTR _.T575_$[esp+112]
	add	esp, 16					; 00000010H
	push	eax
	push	8
	push	OFFSET FLAT:_7B12
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _B12, eax
	mov	eax, DWORD PTR _MMAX
	push	eax
	push	1
	push	1
	push	OFFSET FLAT:_7EVEC
	call	_make_info
	lea	eax, DWORD PTR _.T575_$[esp+112]
	add	esp, 16					; 00000010H
	push	eax
	push	8
	push	OFFSET FLAT:_7EVEC
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _EVEC, eax
	mov	eax, DWORD PTR _MMAX
	push	eax
	push	1
	push	eax
	push	1
	push	2
	push	OFFSET FLAT:_7S
	call	_make_info
	lea	eax, DWORD PTR _.T575_$[esp+120]
	add	esp, 24					; 00000018H
	push	eax
	push	8
	push	OFFSET FLAT:_7S
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _S, eax
	mov	eax, DWORD PTR _MMAX
	push	eax
	push	1
	push	eax
	push	1
	push	2
	push	OFFSET FLAT:_7H
	call	_make_info
	lea	eax, DWORD PTR _.T575_$[esp+120]
	add	esp, 24					; 00000018H
	push	eax
	push	8
	push	OFFSET FLAT:_7H
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _H, eax
	mov	eax, DWORD PTR _MMAX
	push	eax
	push	1
	push	1
	push	OFFSET FLAT:_7CO
	call	_make_info
	add	esp, 16					; 00000010H
	lea	eax, DWORD PTR _.T575_$[esp+96]
	push	eax
	push	8
	push	OFFSET FLAT:_7CO
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _CO, eax
	mov	eax, DWORD PTR _MMAX
	lea	eax, DWORD PTR [eax+eax*2]
	add	eax, eax
	push	eax
	push	1
	push	1
	push	OFFSET FLAT:_7GRAD
	call	_make_info
	lea	ecx, DWORD PTR _.T575_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7GRAD
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _MMAX
	mov	DWORD PTR _GRAD, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7T1
	call	_make_info
	lea	ecx, DWORD PTR _.T575_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7T1
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _MMAX
	mov	DWORD PTR _T1, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7T2
	call	_make_info
	lea	ecx, DWORD PTR _.T575_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7T2
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _MMAX
	mov	DWORD PTR _T2, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7T3
	call	_make_info
	lea	ecx, DWORD PTR _.T575_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7T3
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	ecx, DWORD PTR _MMAX
	mov	DWORD PTR _T3, eax
	npad	1
	push	ecx
	push	1
	push	1
	push	OFFSET FLAT:_7AA1
	call	_make_info
	lea	ecx, DWORD PTR _.T575_$[esp+112]
	add	esp, 16					; 00000010H
	push	ecx
	push	8
	push	OFFSET FLAT:_7AA1
	call	__Allocate_i
	mov	ecx, DWORD PTR _.T575_$[esp+108]
	add	esp, 12					; 0000000cH
	mov	edx, DWORD PTR _EVEC
	mov	DWORD PTR _AA1, eax
	push	edx
	mov	eax, DWORD PTR _B12
	push	eax
	mov	edx, DWORD PTR _RY2
	mov	eax, DWORD PTR _RX2
	mov	DWORD PTR _ERROR, ecx

; 918  : 	-			AA2(MMAX),RX2(MMAX),RY2(MMAX),RZ2(MMAX),
; 919  : 	-			B12(MMAX),EVEC(MMAX),S(MMAX,MMAX),H(MMAX,MMAX),
; 920  :      -			CO(MMAX),GRAD(MMAX*6),T1(MMAX),T2(MMAX),
; 921  :      -			T3(MMAX),STAT=ERROR)
; 922  : 
; 923  : 		CALL GUESS(METHODGEN,M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,

	mov	ecx, DWORD PTR _RZ2
	push	ecx
	push	edx
	push	eax
	mov	ecx, DWORD PTR _AA2
	push	ecx
	mov	edx, DWORD PTR _RZ1
	push	edx
	mov	eax, DWORD PTR _RY1
	push	eax
	mov	ecx, DWORD PTR _RX1
	push	ecx
	mov	edx, DWORD PTR _AA1
	push	edx
	push	OFFSET FLAT:_M
	push	OFFSET FLAT:_METHODGEN
	call	_GUESS@48

; 924  : 	-		EVEC)
; 925  :  		DO I=1,M

	mov	edx, DWORD PTR _M
	mov	ecx, DWORD PTR _T3
	mov	DWORD PTR _I$[esp+96], 1
	test	edx, edx
	jle	$L576
$L98336:

; 926  : 			T1(I)=AA1(I)	

	mov	edi, DWORD PTR _7AA1+8
	mov	eax, DWORD PTR _I$[esp+96]
	add	edi, eax
	mov	ebx, DWORD PTR _7T1+8
	mov	eax, DWORD PTR _AA1
	mov	ebp, DWORD PTR _I$[esp+96]
	add	ebx, ebp
	mov	esi, DWORD PTR _T1
	mov	ebp, DWORD PTR [eax+edi*8+4]
	dec	edx
	mov	eax, DWORD PTR [eax+edi*8]
	mov	DWORD PTR [esi+ebx*8+4], ebp
	mov	edi, DWORD PTR _T2
	mov	DWORD PTR [esi+ebx*8], eax

; 927  : 			T2(I)=AA2(I)	

	mov	ebp, DWORD PTR _7AA2+8
	mov	eax, DWORD PTR _I$[esp+96]
	add	ebp, eax
	mov	ebx, DWORD PTR _7T2+8
	mov	eax, DWORD PTR _AA2
	mov	esi, DWORD PTR _I$[esp+96]
	add	ebx, esi
	npad	1
	mov	esi, DWORD PTR [eax+ebp*8+4]
	mov	eax, DWORD PTR [eax+ebp*8]
	mov	DWORD PTR [edi+ebx*8+4], esi
	mov	DWORD PTR [edi+ebx*8], eax

; 928  : 			T3(I)=B12(I)	

	mov	edi, DWORD PTR _7B12+8
	mov	eax, DWORD PTR _I$[esp+96]
	mov	ebx, DWORD PTR _7T3+8
	add	edi, eax
	mov	ebp, DWORD PTR _I$[esp+96]
	mov	eax, DWORD PTR _B12
	add	ebx, ebp
	mov	esi, DWORD PTR [eax+edi*8+4]
	mov	eax, DWORD PTR [eax+edi*8]
	mov	DWORD PTR [ecx+ebx*8+4], esi
	mov	DWORD PTR [ecx+ebx*8], eax

; 924  : 	-		EVEC)
; 925  :  		DO I=1,M

	mov	ebp, DWORD PTR _I$[esp+96]
	inc	ebp
	test	edx, edx
	mov	DWORD PTR _I$[esp+96], ebp
	jne	$L98336
$L576:

; 929  : 		END DO
; 930  : C
; 931  : C
; 932  : 		IF (RELOADFLAG.EQ.1) THEN

	mov	eax, DWORD PTR _RELOADFLAG
	cmp	eax, 1
	jne	$L578

; 933  : 			WRITE(6,*)'HIT RELOAD!'

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	11					; 0000000bH
	push	OFFSET FLAT:_.SG580_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 934  : 			DO I=1,M-1

	mov	edx, DWORD PTR _M
	mov	ecx, DWORD PTR _T3
	mov	DWORD PTR _I$[esp+96], 1
	dec	edx
	test	edx, edx
	jle	$L581
$L98337:

; 935  : 				EVEC(I)=EVECSTORET(I)

	mov	ebp, DWORD PTR _7EVECSTORET+8
	mov	eax, DWORD PTR _I$[esp+96]
	add	ebp, eax
	mov	edi, DWORD PTR _7EVEC+8
	mov	eax, DWORD PTR _EVECSTORET
	mov	ebx, DWORD PTR _I$[esp+96]
	add	edi, ebx
	dec	edx
	mov	ebx, DWORD PTR _EVEC
	mov	esi, DWORD PTR [eax+ebp*8+4]
	mov	eax, DWORD PTR [eax+ebp*8]
	npad	1
	mov	DWORD PTR [ebx+edi*8+4], esi
	mov	esi, DWORD PTR _T1
	mov	DWORD PTR [ebx+edi*8], eax

; 936  : 				T1(I)=TSTORE1T(I)

	mov	eax, DWORD PTR _7TSTORE1T+8
	mov	ebx, DWORD PTR _I$[esp+96]
	mov	edi, DWORD PTR _TSTORE1T
	add	eax, ebx
	mov	ebp, DWORD PTR _I$[esp+96]
	mov	ebx, DWORD PTR _7T1+8
	add	ebx, ebp
	mov	ebp, DWORD PTR [edi+eax*8+4]
	mov	eax, DWORD PTR [edi+eax*8]
	mov	DWORD PTR [esi+ebx*8+4], ebp
	mov	DWORD PTR [esi+ebx*8], eax

; 937  : 				T2(I)=TSTORE2T(I)

	mov	eax, DWORD PTR _7TSTORE2T+8
	mov	ebx, DWORD PTR _I$[esp+96]
	mov	esi, DWORD PTR _TSTORE2T
	add	eax, ebx
	mov	edi, DWORD PTR _7T2+8
	add	edi, ebx
	mov	ebx, DWORD PTR _T2
	mov	ebp, DWORD PTR [esi+eax*8+4]
	mov	eax, DWORD PTR [esi+eax*8]
	mov	esi, DWORD PTR _TSTORE3T
	mov	DWORD PTR [ebx+edi*8+4], ebp
	mov	DWORD PTR [ebx+edi*8], eax

; 938  : 				T3(I)=TSTORE3T(I)

	mov	eax, DWORD PTR _7TSTORE3T+8
	mov	ebx, DWORD PTR _I$[esp+96]
	mov	ebp, DWORD PTR _I$[esp+96]
	add	eax, ebx
	mov	ebx, DWORD PTR _7T3+8
	add	ebx, ebp
	npad	1
	mov	edi, DWORD PTR [esi+eax*8+4]
	mov	eax, DWORD PTR [esi+eax*8]
	mov	esi, DWORD PTR _RZ1STORET
	mov	DWORD PTR [ecx+ebx*8+4], edi
	mov	DWORD PTR [ecx+ebx*8], eax

; 939  : 				RZ1(I)=RZ1STORET(I)

	mov	eax, DWORD PTR _I$[esp+96]
	mov	ebp, DWORD PTR _7RZ1STORET+8
	mov	edi, DWORD PTR _7RZ1+8
	add	ebp, eax
	add	edi, eax
	mov	ebx, DWORD PTR _RZ1
	npad	1
	mov	eax, DWORD PTR [esi+ebp*8+4]
	mov	ebp, DWORD PTR [esi+ebp*8]
	mov	DWORD PTR [ebx+edi*8+4], eax
	mov	DWORD PTR [ebx+edi*8], ebp

; 940  : 				RZ2(I)=RZ2STORET(I)

	mov	ebp, DWORD PTR _RZ2STORET
	mov	esi, DWORD PTR _7RZ2STORET+8
	mov	eax, DWORD PTR _I$[esp+96]
	add	esi, eax
	mov	edi, DWORD PTR _7RZ2+8
	add	edi, eax
	mov	ebx, DWORD PTR _RZ2
	mov	eax, DWORD PTR [ebp+esi*8+4]
	mov	ebp, DWORD PTR [ebp+esi*8]
	mov	DWORD PTR [ebx+edi*8+4], eax
	mov	DWORD PTR [ebx+edi*8], ebp

; 934  : 			DO I=1,M-1

	mov	eax, DWORD PTR _I$[esp+96]
	inc	eax
	test	edx, edx
	mov	DWORD PTR _I$[esp+96], eax
	jne	$L98337
$L581:

; 941  : 			END DO
; 942  : 			EVEC(M)=0.00D+00		

	mov	edx, DWORD PTR _7EVEC+8
	mov	eax, DWORD PTR _M
	add	edx, eax
	mov	eax, DWORD PTR _EVEC
	mov	DWORD PTR [eax+edx*8], 0
	mov	DWORD PTR [eax+edx*8+4], 0
$L578:

; 943  : 
; 944  : 		END IF
; 945  : 
; 946  : C
; 947  : 		CALL GRADOPT(M,T1,RX1,RY1,RZ1,T2,RX2,RY2,RZ2,T3,EVEC,GRAD)

	mov	eax, DWORD PTR _GRAD
	mov	edx, DWORD PTR _EVEC
	push	eax
	push	edx
	push	ecx
	mov	eax, DWORD PTR _RY2
	mov	ecx, DWORD PTR _RZ2
	mov	edx, DWORD PTR _RX2
	push	ecx
	push	eax
	push	edx
	mov	ecx, DWORD PTR _T2
	push	ecx
	mov	eax, DWORD PTR _RZ1
	push	eax
	mov	edx, DWORD PTR _RY1
	push	edx
	mov	ecx, DWORD PTR _RX1
	push	ecx
	mov	esi, DWORD PTR _T1
	push	esi
	push	OFFSET FLAT:_M
	call	_GRADOPT@48

; 948  : 		DO I=1,M

	mov	ecx, DWORD PTR _M
	mov	esi, DWORD PTR _T1
	mov	edi, DWORD PTR _T3
	mov	DWORD PTR -36+[esp+96], ecx
	mov	DWORD PTR _I$[esp+96], 1
	test	ecx, ecx
	jle	$L583
$L98338:
	mov	eax, DWORD PTR _7T1+8
	mov	ecx, DWORD PTR _I$[esp+96]
	add	eax, ecx
	mov	edx, DWORD PTR _I$[esp+96]
	mov	ecx, DWORD PTR _7T3+8
	add	ecx, edx
	lea	eax, DWORD PTR [esi+eax*8]

; 949  : 			AA1(I)=T1(I)**2 + T1(I)*T3(I)

	fld	QWORD PTR [edi+ecx*8]
	fadd	QWORD PTR [eax]
	mov	ecx, DWORD PTR _AA1
	fmul	QWORD PTR [eax]
	mov	eax, DWORD PTR _7AA1+8
	add	eax, edx
	fstp	QWORD PTR [ecx+eax*8]
	mov	eax, DWORD PTR _7T3+8
	mov	edx, DWORD PTR _I$[esp+96]
	add	eax, edx
	mov	ecx, DWORD PTR _7T1+8
	add	ecx, edx
	lea	eax, DWORD PTR [edi+eax*8]

; 950  : 			AA2(I)=T2(I)**2 + T3(I)**2 + T1(I)*T3(I)

	fld	QWORD PTR [esi+ecx*8]
	fadd	QWORD PTR [eax]
	mov	ecx, DWORD PTR _7T2+8
	add	ecx, edx
	fmul	QWORD PTR [eax]
	mov	eax, DWORD PTR _T2
	fld	QWORD PTR [eax+ecx*8]
	fld	QWORD PTR $T98441
	call	__FIfexp
	faddp	ST(1), ST(0)
	mov	ecx, DWORD PTR _7AA2+8
	mov	eax, DWORD PTR _I$[esp+96]
	add	ecx, eax
	mov	edx, DWORD PTR _AA2
	fstp	QWORD PTR [edx+ecx*8]

; 951  : 			B12(I)=-T1(I)*T3(I)

	mov	eax, DWORD PTR _7T3+8
	mov	ebx, DWORD PTR _I$[esp+96]
	add	eax, ebx
	mov	ecx, DWORD PTR _7B12+8
	add	ecx, ebx
	fld	QWORD PTR [edi+eax*8]
	mov	eax, DWORD PTR _7T1+8
	add	eax, ebx
	fmul	QWORD PTR [esi+eax*8]
	mov	eax, DWORD PTR _B12
	fchs
	fstp	QWORD PTR [eax+ecx*8]

; 952  : C
; 953  : 			EVECSTORE(I,K)=EVEC(I)

	mov	ebp, DWORD PTR _7EVEC+8
	mov	ebx, DWORD PTR _I$[esp+96]
	add	ebp, ebx
	mov	eax, DWORD PTR _EVEC
	mov	edx, DWORD PTR _7EVECSTORE+36
	mov	ecx, DWORD PTR _7EVECSTORE+8
	imul	edx, DWORD PTR _K$[esp+96]
	add	edx, ecx
	mov	ecx, DWORD PTR [eax+ebp*8+4]
	add	edx, ebx
	mov	eax, DWORD PTR [eax+ebp*8]
	mov	ebx, DWORD PTR _EVECSTORE
	mov	DWORD PTR [ebx+edx*8+4], ecx
	mov	DWORD PTR [ebx+edx*8], eax

; 954  : 			TSTORE1(I,K)=T1(I)

	mov	eax, DWORD PTR _I$[esp+96]
	mov	ecx, DWORD PTR _7T1+8
	mov	ebx, DWORD PTR _7TSTORE1+36
	add	ecx, eax
	mov	edx, DWORD PTR _TSTORE1
	imul	ebx, DWORD PTR _K$[esp+96]
	mov	eax, DWORD PTR _7TSTORE1+8
	mov	ebp, DWORD PTR [esi+ecx*8+4]
	add	ebx, eax
	mov	eax, DWORD PTR _I$[esp+96]
	add	ebx, eax
	mov	eax, DWORD PTR [esi+ecx*8]
	mov	DWORD PTR [edx+ebx*8+4], ebp
	mov	DWORD PTR [edx+ebx*8], eax

; 955  : 			TSTORE2(I,K)=T2(I)

	mov	eax, DWORD PTR _I$[esp+96]
	mov	ebp, DWORD PTR _7T2+8
	mov	ebx, DWORD PTR _7TSTORE2+36
	add	ebp, eax
	mov	ecx, DWORD PTR _7TSTORE2+8
	mov	eax, DWORD PTR _T2
	mov	edx, DWORD PTR _TSTORE2
	imul	ebx, DWORD PTR _K$[esp+96]
	add	ebx, ecx
	mov	ecx, DWORD PTR _I$[esp+96]
	add	ebx, ecx
	mov	ecx, DWORD PTR [eax+ebp*8+4]
	mov	eax, DWORD PTR [eax+ebp*8]
	mov	DWORD PTR [edx+ebx*8+4], ecx
	mov	DWORD PTR [edx+ebx*8], eax

; 956  : 			TSTORE3(I,K)=T3(I)

	mov	eax, DWORD PTR _I$[esp+96]
	mov	ecx, DWORD PTR _7T3+8
	mov	edx, DWORD PTR _7TSTORE3+36
	add	ecx, eax
	mov	ebx, DWORD PTR _TSTORE3
	imul	edx, DWORD PTR _K$[esp+96]
	add	edx, DWORD PTR _7TSTORE3+8
	mov	eax, DWORD PTR _I$[esp+96]
	add	edx, eax
	mov	ebp, DWORD PTR [edi+ecx*8+4]
	mov	eax, DWORD PTR [edi+ecx*8]
	mov	DWORD PTR [ebx+edx*8+4], ebp
	mov	DWORD PTR [ebx+edx*8], eax

; 957  : 			RZ1STORE(I,K)=RZ1(I)

	mov	eax, DWORD PTR _7RZ1+8
	mov	ecx, DWORD PTR _I$[esp+96]
	mov	ebx, DWORD PTR _7RZ1STORE+36
	add	eax, ecx
	mov	edx, DWORD PTR _7RZ1STORE+8
	imul	ebx, DWORD PTR _K$[esp+96]
	mov	ecx, DWORD PTR _RZ1
	add	ebx, edx
	add	ebx, DWORD PTR _I$[esp+96]
	mov	ebp, DWORD PTR [ecx+eax*8+4]
	mov	edx, DWORD PTR _RZ1STORE
	mov	eax, DWORD PTR [ecx+eax*8]
	mov	DWORD PTR [edx+ebx*8+4], ebp
	mov	DWORD PTR [edx+ebx*8], eax

; 958  : 			RZ2STORE(I,K)=RZ2(I)

	mov	eax, DWORD PTR _I$[esp+96]
	mov	ebp, DWORD PTR _7RZ2+8
	mov	ebx, DWORD PTR _K$[esp+96]
	add	ebp, eax
	mov	ecx, DWORD PTR _7RZ2STORE+8
	mov	eax, DWORD PTR _RZ2
	mov	edx, DWORD PTR _RZ2STORE
	imul	ebx, DWORD PTR _7RZ2STORE+36
	add	ebx, ecx
	mov	ecx, DWORD PTR _I$[esp+96]
	add	ebx, ecx
	mov	ecx, DWORD PTR [eax+ebp*8+4]
	mov	eax, DWORD PTR [eax+ebp*8]
	npad	1
	mov	DWORD PTR [edx+ebx*8+4], ecx
	mov	ecx, DWORD PTR -36+[esp+96]
	dec	ecx
	mov	DWORD PTR [edx+ebx*8], eax

; 948  : 		DO I=1,M

	mov	eax, DWORD PTR _I$[esp+96]
	mov	DWORD PTR -36+[esp+96], ecx
	inc	eax
	test	ecx, ecx
	mov	DWORD PTR _I$[esp+96], eax
	jne	$L98338
$L583:

; 959  : 		END DO
; 960  : 	
; 961  : 	CALL HAMILTON(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,EVEC,S,H,

	push	OFFSET FLAT:_DFLAG
	mov	eax, DWORD PTR _T2
	push	edi
	mov	ecx, DWORD PTR _GRAD
	push	eax
	mov	edx, DWORD PTR _CO
	push	esi
	mov	eax, DWORD PTR _H
	push	ecx
	push	OFFSET FLAT:_EMIN
	push	edx
	mov	ecx, DWORD PTR _S
	push	eax
	mov	edx, DWORD PTR _EVEC
	push	ecx
	mov	eax, DWORD PTR _B12
	push	edx
	mov	ecx, DWORD PTR _RZ2
	push	eax
	mov	edx, DWORD PTR _RY2
	push	ecx
	mov	eax, DWORD PTR _RX2
	push	edx
	mov	ecx, DWORD PTR _AA2
	push	eax
	mov	edx, DWORD PTR _RZ1
	push	ecx
	mov	eax, DWORD PTR _RY1
	push	edx
	mov	ecx, DWORD PTR _RX1
	push	eax
	mov	edx, DWORD PTR _AA1
	push	ecx
	push	edx
	push	OFFSET FLAT:_M
	call	_HAMILTON@80

; 962  : 	-				CO,EMIN,GRAD,T1,T2,T3,DFLAG)
; 963  : 		ENRGY(K)=EMIN

	mov	ecx, DWORD PTR _7ENRGY+8
	mov	eax, DWORD PTR _K$[esp+96]
	add	ecx, eax
	mov	edx, DWORD PTR _ENRGY
	mov	ebx, DWORD PTR _EMIN+4
	mov	ebp, DWORD PTR -24+[esp+96]
	dec	ebp
	mov	esi, DWORD PTR _K$[esp+96]
	inc	esi
	mov	DWORD PTR [edx+ecx*8+4], ebx
	mov	eax, DWORD PTR _EMIN
	mov	DWORD PTR -24+[esp+96], ebp
	mov	DWORD PTR [edx+ecx*8], eax
	mov	DWORD PTR _K$[esp+96], esi

; 910  : 	- BOX! ...ABORTING!'
; 911  : C
; 912  : 	DO K=1,LOOP

	test	ebp, ebp
	jne	$L98335
$L573:

; 964  : 	END DO			

	mov	edx, DWORD PTR _7ENRGY+8
	mov	eax, DWORD PTR _ENRGY
	mov	DWORD PTR _INDEX, 1
	mov	ecx, DWORD PTR [eax+edx*8+12]
	lea	edx, DWORD PTR [eax+edx*8]

; 965  : 	TENRGY=ENRGY(1)

	mov	esi, DWORD PTR _LOOP
	mov	DWORD PTR _TENRGY+4, ecx
	mov	eax, DWORD PTR [edx+8]
	dec	esi
	mov	ecx, 2
	mov	DWORD PTR _TENRGY, eax

; 966  : 	INDEX=1
; 967  :       DO K=2,LOOP

	test	esi, esi
	jle	SHORT $L587
	add	edx, 16					; 00000010H
$L98339:

; 968  : 		IF (TENRGY.GT.ENRGY(K)) THEN

	fld	QWORD PTR [edx]
	fcomp	QWORD PTR _TENRGY
	fwait
	fnstsw	ax
	test	ah, 1
	je	SHORT $L588

; 969  : 		INDEX=K

	mov	eax, DWORD PTR [edx+4]
	mov	ebx, DWORD PTR [edx]
	mov	DWORD PTR _INDEX, ecx
	mov	DWORD PTR _TENRGY+4, eax

; 970  : 		TENRGY=ENRGY(K)

	mov	DWORD PTR _TENRGY, ebx
$L588:

; 966  : 	INDEX=1
; 967  :       DO K=2,LOOP

	dec	esi
	add	edx, 8
	inc	ecx
	test	esi, esi
	jne	SHORT $L98339
$L587:

; 971  : 		END IF
; 972  : 	END DO
; 973  : 	WRITE(6,*)'LOWEST ENERGY = ',ENRGY(INDEX)

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	16					; 00000010H
	push	OFFSET FLAT:_.SG591_
	call	_lwrite_ch
	add	esp, 8
	mov	edx, DWORD PTR _7ENRGY+8
	add	edx, DWORD PTR _INDEX
	mov	eax, DWORD PTR _ENRGY
	mov	ecx, DWORD PTR [eax+edx*8+4]
	mov	edx, DWORD PTR [eax+edx*8]
	push	ecx
	push	edx
	call	_lwrite_d
	add	esp, 8
	call	_f90_endwrite

; 974  : 	RELOADFLAG=1

	mov	DWORD PTR _RELOADFLAG, 1

; 975  : 	OPEN (9, FILE = OPTIMIZEDWAVE, ACCESS = 'SEQUENTIAL',

	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	0
	push	10					; 0000000aH
	push	OFFSET FLAT:_.SG593_
	push	7
	push	OFFSET FLAT:_.SG592_
	push	66					; 00000042H
	push	OFFSET FLAT:_OPTIMIZEDWAVE
	push	9
	call	_f90_open
	add	esp, 124				; 0000007cH

; 976  :      -		STATUS = 'UNKNOWN')
; 977  : 	REWIND(9)

	push	0
	push	9
	call	_f90_rewind
	add	esp, 8

; 978  : 	WRITE(9,*)':BEGIN:'

	push	0
	push	0
	push	0
	push	0
	push	0
	push	9
	call	_f90_write
	add	esp, 24					; 00000018H
	push	7
	push	OFFSET FLAT:_.SG595_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 979  : 	DO I=1,M

	mov	DWORD PTR _I$[esp+96], 1
	mov	esi, DWORD PTR _M
	test	esi, esi
	jle	$L596
$L98340:

; 980  : 			EVECSTORET(I)=EVECSTORE(I,INDEX)

	mov	eax, DWORD PTR _INDEX
	mov	ecx, DWORD PTR _7EVECSTORE+8
	mov	edx, DWORD PTR _I$[esp+96]
	mov	edi, DWORD PTR _7EVECSTORET+8
	add	edi, edx
	mov	ebx, DWORD PTR _EVECSTORET
	imul	eax, DWORD PTR _7EVECSTORE+36
	add	eax, ecx
	push	0
	add	eax, edx
	mov	ecx, DWORD PTR _EVECSTORE
	push	0
	push	0
	mov	ebp, DWORD PTR [ecx+eax*8+4]
	push	0
	mov	eax, DWORD PTR [ecx+eax*8]
	push	0
	mov	ecx, DWORD PTR _INDEX
	push	9
	mov	DWORD PTR [ebx+edi*8+4], ebp
	mov	DWORD PTR [ebx+edi*8], eax

; 981  : 			TSTORE1T(I)=TSTORE1(I,INDEX)

	mov	eax, DWORD PTR _7TSTORE1+8
	mov	edx, DWORD PTR _I$[esp+120]
	mov	ebx, DWORD PTR _7TSTORE1T+8
	imul	ecx, DWORD PTR _7TSTORE1+36
	add	ecx, eax
	mov	eax, DWORD PTR _TSTORE1
	add	ecx, edx
	add	ebx, edx
	mov	edi, DWORD PTR _TSTORE1T
	mov	ebp, DWORD PTR [eax+ecx*8+4]
	mov	eax, DWORD PTR [eax+ecx*8]
	mov	DWORD PTR [edi+ebx*8+4], ebp
	mov	DWORD PTR [edi+ebx*8], eax

; 982  : 			TSTORE2T(I)=TSTORE2(I,INDEX)

	mov	eax, DWORD PTR _7TSTORE2+8
	mov	ecx, DWORD PTR _7TSTORE2+36
	mov	edx, DWORD PTR _I$[esp+120]
	mov	edi, DWORD PTR _7TSTORE2T+8
	mov	ebx, DWORD PTR _TSTORE2T
	imul	ecx, DWORD PTR _INDEX
	add	ecx, eax
	add	edi, edx
	add	ecx, edx
	mov	eax, DWORD PTR _TSTORE2
	mov	ebp, DWORD PTR [eax+ecx*8+4]
	mov	eax, DWORD PTR [eax+ecx*8]
	mov	DWORD PTR [ebx+edi*8+4], ebp
	mov	DWORD PTR [ebx+edi*8], eax

; 983  : 			TSTORE3T(I)=TSTORE3(I,INDEX)

	mov	eax, DWORD PTR _INDEX
	mov	ecx, DWORD PTR _7TSTORE3+8
	mov	edx, DWORD PTR _I$[esp+120]
	mov	edi, DWORD PTR _7TSTORE3T+8
	mov	ebx, DWORD PTR _TSTORE3T
	imul	eax, DWORD PTR _7TSTORE3+36
	add	eax, ecx
	add	edi, edx
	add	eax, edx
	mov	ecx, DWORD PTR _TSTORE3
	mov	ebp, DWORD PTR [ecx+eax*8+4]
	mov	eax, DWORD PTR [ecx+eax*8]
	mov	DWORD PTR [ebx+edi*8+4], ebp
	mov	DWORD PTR [ebx+edi*8], eax

; 984  : 			RZ1STORET(I)=RZ1STORE(I,INDEX)

	mov	eax, DWORD PTR _INDEX
	mov	ecx, DWORD PTR _7RZ1STORE+8
	mov	edx, DWORD PTR _I$[esp+120]
	mov	edi, DWORD PTR _7RZ1STORET+8
	mov	ebx, DWORD PTR _RZ1STORET
	imul	eax, DWORD PTR _7RZ1STORE+36
	add	eax, ecx
	add	edi, edx
	add	eax, edx
	mov	ecx, DWORD PTR _RZ1STORE
	mov	edx, DWORD PTR _INDEX
	npad	1
	mov	ebp, DWORD PTR [ecx+eax*8+4]
	mov	eax, DWORD PTR [ecx+eax*8]
	mov	ecx, DWORD PTR _RZ2STORE
	mov	DWORD PTR [ebx+edi*8+4], ebp
	mov	DWORD PTR [ebx+edi*8], eax

; 985  : 			RZ2STORET(I)=RZ2STORE(I,INDEX)

	mov	eax, DWORD PTR _7RZ2STORE+8
	mov	ebx, DWORD PTR _7RZ2STORET+8
	mov	edi, DWORD PTR _RZ2STORET
	imul	edx, DWORD PTR _7RZ2STORE+36
	add	edx, eax
	mov	eax, DWORD PTR _I$[esp+120]
	add	edx, eax
	add	ebx, eax
	mov	ebp, DWORD PTR [ecx+edx*8+4]
	mov	eax, DWORD PTR [ecx+edx*8]
	mov	DWORD PTR [edi+ebx*8+4], ebp
	mov	DWORD PTR [edi+ebx*8], eax

; 986  : 	WRITE(9,*) I,EVECSTORET(I),TSTORE1T(I),TSTORE2T(I),TSTORE3T(I),

	call	_f90_write
	mov	eax, DWORD PTR _I$[esp+120]
	add	esp, 24					; 00000018H
	push	eax
	call	_lwrite_i
	mov	eax, DWORD PTR _I$[esp+100]
	add	esp, 4
	mov	ecx, DWORD PTR _7EVECSTORET+8
	mov	edx, DWORD PTR _EVECSTORET
	add	ecx, eax
	mov	ebx, DWORD PTR [edx+ecx*8+4]
	mov	eax, DWORD PTR [edx+ecx*8]
	push	ebx
	dec	esi
	push	eax
	call	_lwrite_d
	mov	eax, DWORD PTR _I$[esp+104]
	add	esp, 8
	mov	ecx, DWORD PTR _7TSTORE1T+8
	mov	edx, DWORD PTR _TSTORE1T
	add	ecx, eax
	mov	ebx, DWORD PTR [edx+ecx*8+4]
	mov	eax, DWORD PTR [edx+ecx*8]
	push	ebx
	push	eax
	call	_lwrite_d
	mov	eax, DWORD PTR _I$[esp+104]
	add	esp, 8
	mov	ecx, DWORD PTR _7TSTORE2T+8
	mov	edx, DWORD PTR _TSTORE2T
	add	ecx, eax
	mov	ebx, DWORD PTR [edx+ecx*8+4]
	mov	eax, DWORD PTR [edx+ecx*8]
	push	ebx
	push	eax
	call	_lwrite_d
	mov	eax, DWORD PTR _I$[esp+104]
	add	esp, 8
	mov	ecx, DWORD PTR _7TSTORE3T+8
	mov	edx, DWORD PTR _TSTORE3T
	add	ecx, eax
	mov	ebx, DWORD PTR [edx+ecx*8+4]
	mov	eax, DWORD PTR [edx+ecx*8]
	push	ebx
	push	eax
	call	_lwrite_d
	mov	eax, DWORD PTR _I$[esp+104]
	add	esp, 8
	mov	edx, DWORD PTR _7RZ1STORET+8
	mov	ecx, DWORD PTR _RZ1STORET
	add	edx, eax
	mov	ebp, DWORD PTR [ecx+edx*8+4]
	mov	eax, DWORD PTR [ecx+edx*8]
	push	ebp
	push	eax
	call	_lwrite_d
	mov	ecx, DWORD PTR _I$[esp+104]
	add	esp, 8
	mov	eax, DWORD PTR _7RZ2STORET+8
	mov	edx, DWORD PTR _RZ2STORET
	add	eax, ecx
	mov	ebx, DWORD PTR [edx+eax*8+4]
	mov	eax, DWORD PTR [edx+eax*8]
	push	ebx
	push	eax
	call	_lwrite_d
	add	esp, 8
	call	_f90_endwrite

; 979  : 	DO I=1,M

	mov	eax, DWORD PTR _I$[esp+96]
	inc	eax
	test	esi, esi
	mov	DWORD PTR _I$[esp+96], eax
	jne	$L98340
$L596:

; 987  : 	-RZ1STORET(I),RZ2STORET(I)
; 988  : 	END DO
; 989  : 	WRITE(9,*)':END:'

	push	0
	push	0
	push	0
	push	0
	push	0
	push	9
	call	_f90_write
	add	esp, 24					; 00000018H
	push	5
	push	OFFSET FLAT:_.SG598_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 990  : 	WRITE(9,*)'ENERGY=',ENRGY(INDEX)

	push	0
	push	0
	push	0
	push	0
	push	0
	push	9
	call	_f90_write
	add	esp, 24					; 00000018H
	push	7
	push	OFFSET FLAT:_.SG599_
	call	_lwrite_ch
	add	esp, 8
	mov	eax, DWORD PTR _7ENRGY+8
	add	eax, DWORD PTR _INDEX
	mov	ecx, DWORD PTR _ENRGY
	mov	edx, DWORD PTR [ecx+eax*8+4]
	mov	eax, DWORD PTR [ecx+eax*8]
	push	edx
	push	eax
	call	_lwrite_d
	add	esp, 8
	call	_f90_endwrite

; 991  : 	CLOSE (9)

	push	0
	push	0
	push	0
	push	9
	call	_f90_close
	mov	ecx, DWORD PTR -8+[esp+112]
	mov	eax, DWORD PTR _J$[esp+112]
	add	esp, 16					; 00000010H
	dec	ecx

; 901  :      -			TSTORE2T(ENDBASIS),TSTORE3T(ENDBASIS),
; 902  :      -			RZ1STORET(ENDBASIS),RZ2STORET(ENDBASIS),
; 903  :      -			STAT=ERROR)
; 904  : 
; 905  : 		DO J=STARTBASIS,ENDBASIS

	inc	eax
	mov	esi, DWORD PTR _T1
	mov	edi, DWORD PTR _T3
	mov	DWORD PTR -8+[esp+96], ecx
	mov	DWORD PTR _J$[esp+96], eax
	test	ecx, ecx
	jne	$L98334
$L555:

; 992  : 		END DO
; 993  : 	END IF
; 994  : C
; 995  :       DO I=1,M

	mov	DWORD PTR _I$[esp+96], 1
	mov	ebx, DWORD PTR _M
	test	ebx, ebx
	jle	$L600
$L98341:
	mov	eax, DWORD PTR _7T1+8
	mov	ecx, DWORD PTR _I$[esp+96]
	add	eax, ecx
	mov	edx, DWORD PTR _I$[esp+96]
	mov	ecx, DWORD PTR _7T3+8
	dec	ebx
	lea	eax, DWORD PTR [esi+eax*8]
	add	ecx, edx

; 996  : 	    AA1(I)=T1(I)**2 + T1(I)*T3(I)

	fld	QWORD PTR [edi+ecx*8]
	fadd	QWORD PTR [eax]
	mov	ecx, DWORD PTR _AA1
	fmul	QWORD PTR [eax]
	mov	eax, DWORD PTR _7AA1+8
	add	eax, edx
	fstp	QWORD PTR [ecx+eax*8]
	mov	eax, DWORD PTR _7T3+8
	mov	edx, DWORD PTR _I$[esp+96]
	add	eax, edx
	mov	ecx, DWORD PTR _7T1+8
	add	ecx, edx
	lea	eax, DWORD PTR [edi+eax*8]

; 997  :           AA2(I)=T2(I)**2 + T3(I)**2 + T1(I)*T3(I)

	fld	QWORD PTR [esi+ecx*8]
	fadd	QWORD PTR [eax]
	mov	ecx, DWORD PTR _7T2+8
	add	ecx, edx
	fmul	QWORD PTR [eax]
	mov	eax, DWORD PTR _T2
	fld	QWORD PTR [eax+ecx*8]
	fld	QWORD PTR $T98441
	call	__FIfexp
	faddp	ST(1), ST(0)
	mov	ecx, DWORD PTR _7AA2+8
	mov	eax, DWORD PTR _I$[esp+96]
	add	ecx, eax
	mov	edx, DWORD PTR _AA2
	fstp	QWORD PTR [edx+ecx*8]

; 998  :           B12(I)=-T1(I)*T3(I)

	mov	eax, DWORD PTR _7T3+8
	mov	ebp, DWORD PTR _I$[esp+96]
	add	eax, ebp
	mov	ecx, DWORD PTR _7B12+8
	add	ecx, ebp
	fld	QWORD PTR [edi+eax*8]
	mov	eax, DWORD PTR _7T1+8
	add	eax, ebp
	fmul	QWORD PTR [esi+eax*8]
	mov	eax, DWORD PTR _B12
	fchs
	fstp	QWORD PTR [eax+ecx*8]

; 992  : 		END DO
; 993  : 	END IF
; 994  : C
; 995  :       DO I=1,M

	mov	ebp, DWORD PTR _I$[esp+96]
	inc	ebp
	test	ebx, ebx
	mov	DWORD PTR _I$[esp+96], ebp
	jne	$L98341
$L600:

; 999  : 	END DO
; 1000 : C
; 1001 : 	CALL BIOUT(' ')

	push	1
	push	OFFSET FLAT:_.SG602_
	call	_BIOUT@8

; 1002 : 	WRITE(6,*)' Square-Magnitude of the GRADIENT is',GG

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	36					; 00000024H
	push	OFFSET FLAT:_.SG603_
	call	_lwrite_ch
	add	esp, 8
	mov	eax, DWORD PTR _FINALGRAD+4
	mov	ecx, DWORD PTR _FINALGRAD
	npad	1
	push	eax
	push	ecx
	call	_lwrite_d
	add	esp, 8
	call	_f90_endwrite

; 1003 : 	WRITE(7,*)' Square-Magnitude of the GRADIENT is',GG

	push	0
	push	0
	push	0
	push	0
	push	0
	push	7
	call	_f90_write
	add	esp, 24					; 00000018H
	push	36					; 00000024H
	push	OFFSET FLAT:_.SG606_
	call	_lwrite_ch
	add	esp, 8
	mov	ecx, DWORD PTR _FINALGRAD+4
	mov	edx, DWORD PTR _FINALGRAD
	npad	1
	push	ecx
	push	edx
	mov	esi, 1
	call	_lwrite_d
	add	esp, 8
	call	_f90_endwrite

; 1004 : 	CALL BIOUT(' ')

	push	1
	push	OFFSET FLAT:_.SG607_
	call	_BIOUT@8

; 1005 : C
; 1006 : C
; 1007 : C
; 1008 :       J=1
; 1009 :       CALL BIOUT('     LinCoeff:.   AA1:        AA2:        B12:        

	push	69					; 00000045H
	push	OFFSET FLAT:_.SG608_
	call	_BIOUT@8

; 1010 : 	-RZ1:       RZ2:')
; 1011 :       CALL BIOUT('FNC: =================================================

	push	75					; 0000004bH
	push	OFFSET FLAT:_.SG609_
	call	_BIOUT@8

; 1012 :      -=====================')
; 1013 :  606  IF (J .LE. M) THEN

	cmp	DWORD PTR _M, esi
	jl	$L611
$L610:

; 1014 :           WRITE(6,1111) J,EVEC(J),AA1(J),AA2(J),B12(J),RZ1(J),RZ2(J)

	push	0
	push	0
	push	0
	push	-1
	push	OFFSET FLAT:_.SG613_
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	esi
	call	_write_i
	add	esp, 4
	mov	edx, DWORD PTR _7EVEC+8
	add	edx, esi
	mov	eax, DWORD PTR _EVEC
	mov	ecx, DWORD PTR [eax+edx*8+4]
	mov	edx, DWORD PTR [eax+edx*8]
	push	ecx
	push	edx
	call	_write_d
	add	esp, 8
	mov	eax, DWORD PTR _7AA1+8
	add	eax, esi
	mov	ecx, DWORD PTR _AA1
	mov	edx, DWORD PTR [ecx+eax*8+4]
	mov	eax, DWORD PTR [ecx+eax*8]
	push	edx
	push	eax
	call	_write_d
	add	esp, 8
	mov	eax, DWORD PTR _7AA2+8
	add	eax, esi
	mov	ecx, DWORD PTR _AA2
	mov	edx, DWORD PTR [ecx+eax*8+4]
	mov	eax, DWORD PTR [ecx+eax*8]
	push	edx
	push	eax
	call	_write_d
	add	esp, 8
	mov	ecx, DWORD PTR _7B12+8
	add	ecx, esi
	mov	eax, DWORD PTR _B12
	mov	edx, DWORD PTR [eax+ecx*8+4]
	mov	ecx, DWORD PTR [eax+ecx*8]
	push	edx
	push	ecx
	call	_write_d
	add	esp, 8
	mov	ecx, DWORD PTR _7RZ1+8
	add	ecx, esi
	mov	eax, DWORD PTR _RZ1
	mov	edx, DWORD PTR [eax+ecx*8+4]
	mov	ecx, DWORD PTR [eax+ecx*8]
	push	edx
	push	ecx
	call	_write_d
	add	esp, 8
	mov	ecx, DWORD PTR _7RZ2+8
	add	ecx, esi
	mov	eax, DWORD PTR _RZ2
	mov	edx, DWORD PTR [eax+ecx*8+4]
	mov	ecx, DWORD PTR [eax+ecx*8]
	push	edx
	push	ecx
	call	_write_d
	add	esp, 8
	call	_f90_endwrite

; 1015 :           WRITE(7,1111) J,EVEC(J),AA1(J),AA2(J),B12(J),RZ1(J),RZ2(J)

	push	0
	push	0
	push	0
	push	-1
	push	OFFSET FLAT:_.SG616_
	push	7
	call	_f90_write
	add	esp, 24					; 00000018H
	push	esi
	call	_write_i
	add	esp, 4
	mov	ecx, DWORD PTR _7EVEC+8
	add	ecx, esi
	mov	edx, DWORD PTR _EVEC
	mov	eax, DWORD PTR [edx+ecx*8+4]
	mov	ecx, DWORD PTR [edx+ecx*8]
	push	eax
	push	ecx
	call	_write_d
	add	esp, 8
	mov	ecx, DWORD PTR _7AA1+8
	add	ecx, esi
	mov	edx, DWORD PTR _AA1
	mov	eax, DWORD PTR [edx+ecx*8+4]
	mov	ecx, DWORD PTR [edx+ecx*8]
	push	eax
	push	ecx
	call	_write_d
	add	esp, 8
	mov	edx, DWORD PTR _7AA2+8
	add	edx, esi
	mov	ecx, DWORD PTR _AA2
	mov	eax, DWORD PTR [ecx+edx*8+4]
	push	eax
	mov	edx, DWORD PTR [ecx+edx*8]
	push	edx
	call	_write_d
	add	esp, 8
	mov	edx, DWORD PTR _7B12+8
	add	edx, esi
	mov	ecx, DWORD PTR _B12
	mov	eax, DWORD PTR [ecx+edx*8+4]
	mov	edx, DWORD PTR [ecx+edx*8]
	push	eax
	push	edx
	call	_write_d
	add	esp, 8
	mov	edx, DWORD PTR _7RZ1+8
	add	edx, esi
	mov	ecx, DWORD PTR _RZ1
	mov	eax, DWORD PTR [ecx+edx*8+4]
	mov	edx, DWORD PTR [ecx+edx*8]
	push	eax
	push	edx
	call	_write_d
	add	esp, 8
	mov	eax, DWORD PTR _7RZ2+8
	add	eax, esi
	mov	edx, DWORD PTR _RZ2
	inc	esi
	npad	1
	mov	ecx, DWORD PTR [edx+eax*8+4]
	mov	eax, DWORD PTR [edx+eax*8]
	push	ecx
	push	eax
	call	_write_d
	add	esp, 8
	call	_f90_endwrite

; 1016 :           J=J+1
; 1017 : 
; 1018 :       GOTO 606

	cmp	esi, DWORD PTR _M
	jle	$L610
$L611:

; 1019 :       END IF
; 1020 :       CALL BIOUT('------------------------------------------------------

	push	75					; 0000004bH
	push	OFFSET FLAT:_.SG617_
	call	_BIOUT@8

; 1021 :      1---------------------')
; 1022 : C
; 1023 : 	CALL BIOUT(' ')

	push	1
	push	OFFSET FLAT:_.SG618_
	call	_BIOUT@8

; 1024 : C
; 1025 :       DFLAG=0

	push	OFFSET FLAT:_DFLAG
	mov	eax, DWORD PTR _T3
	push	eax
	mov	ecx, DWORD PTR _T2
	push	ecx
	mov	edx, DWORD PTR _T1
	push	edx
	mov	eax, DWORD PTR _GRAD
	push	eax
	mov	ecx, DWORD PTR _CO
	push	OFFSET FLAT:_EMIN
	mov	edx, DWORD PTR _H
	push	ecx
	mov	eax, DWORD PTR _S
	push	edx
	mov	ecx, DWORD PTR _EVEC
	push	eax
	mov	edx, DWORD PTR _B12
	push	ecx
	mov	eax, DWORD PTR _RZ2
	push	edx
	mov	ecx, DWORD PTR _RY2
	push	eax
	mov	edx, DWORD PTR _RX2
	push	ecx
	mov	eax, DWORD PTR _AA2
	push	edx
	mov	ecx, DWORD PTR _RZ1
	push	eax
	mov	edx, DWORD PTR _RY1
	push	ecx
	mov	eax, DWORD PTR _RX1
	push	edx
	mov	ecx, DWORD PTR _AA1
	mov	DWORD PTR _DFLAG, 0

; 1026 : 	CALL HAMILTON(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,EVEC,S,H,

	push	eax
	push	ecx
	push	OFFSET FLAT:_M
	call	_HAMILTON@80

; 1027 : 	-				CO,EMIN,GRAD,T1,T2,T3,DFLAG)
; 1028 : C
; 1029 : 	EBEST=-1.1744758915D+00

	mov	DWORD PTR _EBEST, 998174039		; 3b7eed57H
	mov	DWORD PTR _EBEST+4, -1074607449		; bff2caa7H

; 1030 : C
; 1031 :       WRITE(6,*) '           OPTIMIZED ENERGY = ',EMIN

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	30					; 0000001eH
	push	OFFSET FLAT:_.SG620_
	call	_lwrite_ch
	add	esp, 8
	mov	ecx, DWORD PTR _EMIN+4
	mov	edx, DWORD PTR _EMIN
	npad	1
	push	ecx
	push	edx
	call	_lwrite_d
	add	esp, 8
	call	_f90_endwrite

; 1032 :       WRITE(7,*) '           OPTIMIZED ENERGY = ',EMIN

	push	0
	push	0
	push	0
	push	0
	push	0
	push	7
	call	_f90_write
	add	esp, 24					; 00000018H
	push	30					; 0000001eH
	push	OFFSET FLAT:_.SG621_
	call	_lwrite_ch
	add	esp, 8
	mov	ecx, DWORD PTR _EMIN+4
	mov	edx, DWORD PTR _EMIN
	npad	1
	push	ecx
	push	edx
	call	_lwrite_d
	add	esp, 8
	call	_f90_endwrite

; 1033 :       WRITE(6,*) '           BEST LIT. RESULT = ',EBEST

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	30					; 0000001eH
	push	OFFSET FLAT:_.SG622_
	call	_lwrite_ch
	add	esp, 8
	mov	ecx, DWORD PTR _EBEST+4
	mov	edx, DWORD PTR _EBEST
	npad	1
	push	ecx
	push	edx
	call	_lwrite_d
	add	esp, 8
	call	_f90_endwrite

; 1034 :       WRITE(7,*) '           BEST LIT. RESULT = ',EBEST

	push	0
	push	0
	push	0
	push	0
	push	0
	push	7
	call	_f90_write
	add	esp, 24					; 00000018H
	push	30					; 0000001eH
	push	OFFSET FLAT:_.SG623_
	call	_lwrite_ch
	add	esp, 8
	mov	ecx, DWORD PTR _EBEST+4
	mov	edx, DWORD PTR _EBEST
	npad	1
	push	ecx
	push	edx
	call	_lwrite_d
	add	esp, 8
	call	_f90_endwrite

; 1035 : 	WRITE(6,*) '           ------------------------------------------'

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	53					; 00000035H
	push	OFFSET FLAT:_.SG624_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 1036 : 	WRITE(7,*) '           ------------------------------------------'

	push	0
	push	0
	push	0
	push	0
	push	0
	push	7
	call	_f90_write
	add	esp, 24					; 00000018H
	push	53					; 00000035H
	push	OFFSET FLAT:_.SG625_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 1037 :       WRITE(6,*) '                 DIFFERENCE = ',EMIN-EBEST

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	30					; 0000001eH
	push	OFFSET FLAT:_.SG626_
	call	_lwrite_ch
	fld	QWORD PTR _EMIN
	fsub	QWORD PTR _EBEST
	add	esp, 8
	sub	esp, 8
	fstp	QWORD PTR [esp]
	call	_lwrite_d
	add	esp, 8
	call	_f90_endwrite

; 1038 :       WRITE(7,*) '                 DIFFERENCE = ',EMIN-EBEST

	push	0
	push	0
	push	0
	push	0
	push	0
	push	7
	call	_f90_write
	add	esp, 24					; 00000018H
	push	30					; 0000001eH
	push	OFFSET FLAT:_.SG627_
	call	_lwrite_ch
	fld	QWORD PTR _EMIN
	fsub	QWORD PTR _EBEST
	add	esp, 8
	sub	esp, 8
	fstp	QWORD PTR [esp]
	call	_lwrite_d
	add	esp, 8
	call	_f90_endwrite

; 1039 : C
; 1040 : C     -----------------------------------------------------
; 1041 : C     FORMAT(s):
; 1042 : C     -----------------------------------------------------
; 1043 : 1111  FORMAT(I3,6D12.4)
; 1044 : 1007  FORMAT(/) 
; 1045 : 9001  FORMAT ('     ',I2, '/', I2.2, '/', I4.4)  
; 1046 : 9002  FORMAT ('     ',I2, ':', I2.2, ':', I2.2, ':', I2.2, ' ', A, 'm')
; 1047 : 9003  FORMAT ('     ',I2, '/', I2.2, '/', I4.4)  
; 1048 : 9004  FORMAT ('     ',I2, ':', I2.2, ':', I2.2, ':', I2.2, ' ', A, 'm')
; 1049 : 9005  FORMAT (' Calculation of the FUNCTION requires ',I2, ':', I2.2, 
; 1050 : 	-':', I2.2, ' ', A, 'm')
; 1051 : 9006  FORMAT (' Calculation of the GRADIENT requires ',I2, ':', I2.2, 
; 1052 : 	-':', I2.2, ' ', A, 'm')
; 1053 : C     -----------------------------------------------------
; 1054 : C     TERMINATE PROGRAM EXECUTION!
; 1055 : C     -----------------------------------------------------
; 1056 :       CALL BIOUT(' ')

	push	1
	push	OFFSET FLAT:_.SG628_
	call	_BIOUT@8

; 1057 :       WRITE (7,*) ' PROGRAM FINISHES:  TIME STAMP...'

	push	0
	push	0
	push	0
	push	0
	push	0
	push	7
	call	_f90_write
	add	esp, 24					; 00000018H
	push	33					; 00000021H
	push	OFFSET FLAT:_.SG629_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 1058 :       WRITE (6,*) ' PROGRAM FINISHES:  TIME STAMP...'

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	33					; 00000021H
	push	OFFSET FLAT:_.SG630_
	call	_lwrite_ch
	add	esp, 8
	call	_f90_endwrite

; 1059 : C
; 1060 :       CALL GETDAT(tmpyear, tmpmonth, tmpday)

	push	OFFSET FLAT:_TMPDAY
	push	OFFSET FLAT:_TMPMONTH
	push	OFFSET FLAT:_TMPYEAR
	call	_GETDAT@12

; 1061 :       CALL GETTIM(tmphour, tmpminute, tmpsecond, tmphund)  

	push	OFFSET FLAT:_TMPHUND
	push	OFFSET FLAT:_TMPSECOND
	push	OFFSET FLAT:_TMPMINUTE
	push	OFFSET FLAT:_TMPHOUR
	call	_GETTIM@16

; 1062 : C
; 1063 :       WRITE (6, 9003) tmpmonth, tmpday, tmpyear

	push	0
	push	0
	push	0
	push	-1
	push	OFFSET FLAT:_.SG631_
	push	6
	call	_f90_write
	mov	cx, WORD PTR _TMPMONTH
	add	esp, 24					; 00000018H
	push	ecx
	call	_write_i2
	mov	cx, WORD PTR _TMPDAY
	add	esp, 4
	push	ecx
	call	_write_i2
	mov	cx, WORD PTR _TMPYEAR
	add	esp, 4
	push	ecx
	call	_write_i2
	add	esp, 4
	call	_f90_endwrite

; 1064 :       WRITE (6, 9004) tmphour, tmpminute, tmpsecond, tmphund

	push	0
	push	0
	push	0
	push	-1
	push	OFFSET FLAT:_.SG632_
	push	6
	call	_f90_write
	mov	cx, WORD PTR _TMPHOUR
	add	esp, 24					; 00000018H
	push	ecx
	call	_write_i2
	mov	cx, WORD PTR _TMPMINUTE
	add	esp, 4
	push	ecx
	call	_write_i2
	mov	cx, WORD PTR _TMPSECOND
	add	esp, 4
	push	ecx
	call	_write_i2
	mov	cx, WORD PTR _TMPHUND
	add	esp, 4
	push	ecx
	call	_write_i2
	add	esp, 4
	call	_f90_endwrite

; 1065 : C
; 1066 :       WRITE (7, 9003) tmpmonth, tmpday, tmpyear

	push	0
	push	0
	push	0
	push	-1
	push	OFFSET FLAT:_.SG633_
	push	7
	call	_f90_write
	mov	cx, WORD PTR _TMPMONTH
	add	esp, 24					; 00000018H
	push	ecx
	call	_write_i2
	mov	cx, WORD PTR _TMPDAY
	add	esp, 4
	push	ecx
	call	_write_i2
	mov	cx, WORD PTR _TMPYEAR
	add	esp, 4
	push	ecx
	call	_write_i2
	add	esp, 4
	call	_f90_endwrite

; 1067 :       WRITE (7, 9004) tmphour, tmpminute, tmpsecond, tmphund

	push	0
	push	0
	push	0
	push	-1
	push	OFFSET FLAT:_.SG634_
	push	7
	call	_f90_write
	mov	cx, WORD PTR _TMPHOUR
	add	esp, 24					; 00000018H
	push	ecx
	call	_write_i2
	mov	cx, WORD PTR _TMPMINUTE
	add	esp, 4
	push	ecx
	call	_write_i2
	mov	cx, WORD PTR _TMPSECOND
	add	esp, 4
	push	ecx
	call	_write_i2
	mov	cx, WORD PTR _TMPHUND
	add	esp, 4
	push	ecx
	call	_write_i2
	add	esp, 4
	call	_f90_endwrite

; 1068 : C
; 1069 : 	CLOSE (7)

	push	0
	push	0
	push	0
	push	7
	call	_f90_close
	add	esp, 16					; 00000010H

; 1070 : 	CLOSE (8)

	push	0
	push	0
	push	0
	push	8
	call	_f90_close
	add	esp, 16					; 00000010H

; 1071 : 	CLOSE (9)

	push	0
	push	0
	push	0
	push	9
	call	_f90_close
	add	esp, 16					; 00000010H
	xor	eax, eax

; 1072 :       END

	pop	ebp
	pop	edi
	pop	esi
	pop	ebx
	add	esp, 80					; 00000050H
	ret	0
	npad	4
_main	ENDP
_TEXT	ENDS
END
