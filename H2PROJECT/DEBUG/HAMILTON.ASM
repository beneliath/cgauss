	TITLE	C:\users\dgilmore\Graduate.School\Adamowicz.Research\H2ANL\HAMILTON.FOR
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
PUBLIC	_TNOVER
PUBLIC	_OVER
PUBLIC	_DOVN
PUBLIC	_DTK
PUBLIC	_DOV
PUBLIC	_DEN
PUBLIC	_DER
EXTRN	_ENERGY@24:NEAR
EXTRN	__FFsSTOP:NEAR
EXTRN	__FFljj:NEAR
EXTRN	_ubound_info:NEAR
EXTRN	_GEOM@24:NEAR
EXTRN	_lbound_info:NEAR
EXTRN	_make_info:NEAR
EXTRN	_f90_write:NEAR
EXTRN	_lwrite_ch:NEAR
EXTRN	_lwrite_d:NEAR
EXTRN	_f90_endwrite:NEAR
EXTRN	__Allocate_i:NEAR
EXTRN	_TK@84:NEAR
EXTRN	_OVERLP@88:NEAR
EXTRN	__Deallocate_i:NEAR
EXTRN	_ENUCL@100:NEAR
EXTRN	_lwrite_i:NEAR
EXTRN	_ELREP@88:NEAR
_DATA	SEGMENT
_TNOVER	DB	0cH DUP (?)
_OVER	DQ	01H DUP (?)
_DOVN	DB	028H DUP (?)
_DTK	DB	028H DUP (?)
_DOV	DB	028H DUP (?)
_DEN	DB	050H DUP (?)
_DER	DB	028H DUP (?)
	ORG $-260
_C:\users\dgilmore\Graduate.School\Adamowicz.Research\H2ANL\HAMILTON.FOR DB 
	DB	'C:\users\dgilmore\Graduate.School\Adamowicz.Research\H2ANL\HA'
	DB	'MILTON.FOR', 00H
_DATA	ENDS
PUBLIC	_HAMILTON@80
EXTRN	___FHlink@4:NEAR
EXTRN	___FHovfl@0:NEAR
EXTRN	___FHranger@0:NEAR
EXTRN	___FHunlink@0:NEAR
EXTRN	__fltused:NEAR
_BSS	SEGMENT
_DER1DA2 DQ	01H DUP (?)
_D2EN2DAZ DQ	01H DUP (?)
_D1T2DA1 DQ	01H DUP (?)
_D2EN2DBZ DQ	01H DUP (?)
_D1T2DA2 DQ	01H DUP (?)
_DER1DA5 DQ	01H DUP (?)
_EN11	DQ	01H DUP (?)
_D1T2DA5 DQ	01H DUP (?)
_GRADC	DQ	00H DUP (?)
_D2EN3DA1 DQ	01H DUP (?)
_EN2	DQ	01H DUP (?)
_D1T2DAZ DQ	01H DUP (?)
_DER1DAZ DQ	01H DUP (?)
_D2EN3DA2 DQ	01H DUP (?)
_D1EN2DA1 DQ	01H DUP (?)
_D1T2DBZ DQ	01H DUP (?)
_KIM	DD	01H DUP (?)
	ALIGN	8

_D2EN3DA5 DQ	01H DUP (?)
_D1EN2DA2 DQ	01H DUP (?)
_TEMP	DQ	00H DUP (?)
_SUMH	DQ	01H DUP (?)
_DER1DBZ DQ	01H DUP (?)
_D2EN3DAZ DQ	01H DUP (?)
_D1EN2DA5 DQ	01H DUP (?)
_SUMS	DQ	01H DUP (?)
_D2EN3DBZ DQ	01H DUP (?)
_D1EN2DAZ DQ	01H DUP (?)
_T6	DQ	01H DUP (?)
_ERROR	DD	01H DUP (?)
	ALIGN	8

_EN12	DQ	01H DUP (?)
_D1EN2DBZ DQ	01H DUP (?)
_D2T2DA1 DQ	01H DUP (?)
_D2EN4DA1 DQ	01H DUP (?)
_EN3	DQ	01H DUP (?)
_D2T2DA2 DQ	01H DUP (?)
_ER2	DQ	01H DUP (?)
_D2EN4DA2 DQ	01H DUP (?)
_D1EN3DA1 DQ	01H DUP (?)
_D2T2DA5 DQ	01H DUP (?)
_DER2DA1 DQ	01H DUP (?)
_D2EN4DA5 DQ	01H DUP (?)
_D1EN3DA2 DQ	01H DUP (?)
_D2T2DAZ DQ	01H DUP (?)
_DER2DA2 DQ	01H DUP (?)
_D2EN4DAZ DQ	01H DUP (?)
_D1EN3DA5 DQ	01H DUP (?)
_D2T2DBZ DQ	01H DUP (?)
_DER2DA5 DQ	01H DUP (?)
_D2EN4DBZ DQ	01H DUP (?)
_D1EN3DAZ DQ	01H DUP (?)
_T7	DQ	01H DUP (?)
_EX1	DQ	01H DUP (?)
_ENG	DQ	01H DUP (?)
_DER2DAZ DQ	01H DUP (?)
_EN13	DQ	01H DUP (?)
_D1EN3DBZ DQ	01H DUP (?)
_D3T2DA1 DQ	01H DUP (?)
_EY1	DQ	01H DUP (?)
_DER2DBZ DQ	01H DUP (?)
_D2EN5DA1 DQ	01H DUP (?)
_EN4	DQ	01H DUP (?)
_D3T2DA2 DQ	01H DUP (?)
_EZ1	DQ	01H DUP (?)
_ER3	DQ	01H DUP (?)
_D2EN5DA2 DQ	01H DUP (?)
_D1EN4DA1 DQ	01H DUP (?)
_D3T2DA5 DQ	01H DUP (?)
_EX2	DQ	01H DUP (?)
_DER3DA1 DQ	01H DUP (?)
_D2EN5DA5 DQ	01H DUP (?)
_D1EN4DA2 DQ	01H DUP (?)
_D3T2DAZ DQ	01H DUP (?)
_EY2	DQ	01H DUP (?)
_DER3DA2 DQ	01H DUP (?)
_D2EN5DAZ DQ	01H DUP (?)
_D1EN4DA5 DQ	01H DUP (?)
_D3T2DBZ DQ	01H DUP (?)
_EZ2	DQ	01H DUP (?)
_SCALE	DQ	01H DUP (?)
_DER3DA5 DQ	01H DUP (?)
_D2EN5DBZ DQ	01H DUP (?)
_D1EN4DAZ DQ	01H DUP (?)
_T8	DQ	01H DUP (?)
_DER3DAZ DQ	01H DUP (?)
_D1EN4DBZ DQ	01H DUP (?)
_D4T2DA1 DQ	01H DUP (?)
_MMAX	DD	01H DUP (?)
	ALIGN	8

_DER3DBZ DQ	01H DUP (?)
_EN5	DQ	01H DUP (?)
_D4T2DA2 DQ	01H DUP (?)
_ER4	DQ	01H DUP (?)
_EN14	DQ	01H DUP (?)
_D1EN5DA1 DQ	01H DUP (?)
_D4T2DA5 DQ	01H DUP (?)
_DER4DA1 DQ	01H DUP (?)
_D2EN6DA1 DQ	01H DUP (?)
_D4T2DAZ DQ	01H DUP (?)
_DHDA1	DQ	00H DUP (?)
_NN	DD	01H DUP (?)
	ALIGN	8

_DER4DA2 DQ	01H DUP (?)
_D2EN6DA2 DQ	01H DUP (?)
_D1EN5DA2 DQ	01H DUP (?)
_D4T2DBZ DQ	01H DUP (?)
_NB	DD	01H DUP (?)
	ALIGN	8

_DER4DA5 DQ	01H DUP (?)
_D2EN6DA5 DQ	01H DUP (?)
_DER4DAZ DQ	01H DUP (?)
_D2EN6DAZ DQ	01H DUP (?)
_D1EN5DA5 DQ	01H DUP (?)
_DER4DBZ DQ	01H DUP (?)
_D2EN6DBZ DQ	01H DUP (?)
_T	DQ	01H DUP (?)
_DHDA2	DQ	00H DUP (?)
_EN15	DQ	01H DUP (?)
_D1EN5DAZ DQ	01H DUP (?)
_D2EN7DA1 DQ	01H DUP (?)
_IC	DD	01H DUP (?)
	ALIGN	8

_EER1	DQ	01H DUP (?)
_D2EN7DA2 DQ	01H DUP (?)
_D1EN5DBZ DQ	01H DUP (?)
_T1	DQ	01H DUP (?)
_OV1	DQ	01H DUP (?)
_DHDA5	DQ	00H DUP (?)
_D2EN7DA5 DQ	01H DUP (?)
_D1T1DA1 DQ	01H DUP (?)
_KLM	DD	01H DUP (?)
	ALIGN	8

_D2EN7DAZ DQ	01H DUP (?)
_DOV1DA1 DQ	01H DUP (?)
_DEERDA1 DQ	01H DUP (?)
_D2EN7DBZ DQ	01H DUP (?)
_DHDAZ	DQ	00H DUP (?)
_DEERDA2 DQ	01H DUP (?)
_EN16	DQ	01H DUP (?)
_EN6	DQ	01H DUP (?)
_DEERDA5 DQ	01H DUP (?)
_D2EN8DA1 DQ	01H DUP (?)
_D1EN6DA1 DQ	01H DUP (?)
_D1T1DA2 DQ	01H DUP (?)
_DEERDAZ DQ	01H DUP (?)
_D2EN8DA2 DQ	01H DUP (?)
_D1EN6DA2 DQ	01H DUP (?)
_DOV1DA2 DQ	01H DUP (?)
_DHDBZ	DQ	00H DUP (?)
_DEERDBZ DQ	01H DUP (?)
_D2EN8DA5 DQ	01H DUP (?)
_D1EN6DA5 DQ	01H DUP (?)
_D1T1DA5 DQ	01H DUP (?)
_D2EN8DAZ DQ	01H DUP (?)
_D1EN6DAZ DQ	01H DUP (?)
_DTDA1	DQ	01H DUP (?)
_DOV1DA5 DQ	01H DUP (?)
_D2EN8DBZ DQ	01H DUP (?)
_D1EN6DBZ DQ	01H DUP (?)
_DTDA2	DQ	01H DUP (?)
_D1T1DAZ DQ	01H DUP (?)
_DSDA1	DQ	00H DUP (?)
_EN7	DQ	01H DUP (?)
_DTDA5	DQ	01H DUP (?)
_DOV1DAZ DQ	01H DUP (?)
_D1EN7DA1 DQ	01H DUP (?)
_DTDAZ	DQ	01H DUP (?)
_D1T1DBZ DQ	01H DUP (?)
_ENN1	DQ	01H DUP (?)
_D1EN7DA2 DQ	01H DUP (?)
_DTDBZ	DQ	01H DUP (?)
_DOV1DBZ DQ	01H DUP (?)
_DSDA2	DQ	00H DUP (?)
_DA1	DQ	01H DUP (?)
_D1EN7DA5 DQ	01H DUP (?)
_DA2	DQ	01H DUP (?)
_D1EN7DAZ DQ	01H DUP (?)
_DB	DQ	01H DUP (?)
_DENADA1 DQ	01H DUP (?)
_D1EN7DBZ DQ	01H DUP (?)
_T2	DQ	01H DUP (?)
_DSDA5	DQ	00H DUP (?)
_NX	DD	01H DUP (?)
	ALIGN	8

_DENADA2 DQ	01H DUP (?)
_EN8	DQ	01H DUP (?)
_D2T1DA1 DQ	01H DUP (?)
_OV2	DQ	01H DUP (?)
_DENADA5 DQ	01H DUP (?)
_D1EN8DA1 DQ	01H DUP (?)
_D2T1DA2 DQ	01H DUP (?)
_DOV2DA1 DQ	01H DUP (?)
_EN	DQ	01H DUP (?)
_DENADAZ DQ	01H DUP (?)
_D1EN8DA2 DQ	01H DUP (?)
_D2T1DA5 DQ	01H DUP (?)
_DOV2DA2 DQ	01H DUP (?)
_DSDAZ	DQ	00H DUP (?)
_DENADBZ DQ	01H DUP (?)
_D1EN8DA5 DQ	01H DUP (?)
_D2T1DAZ DQ	01H DUP (?)
_DOV2DA5 DQ	01H DUP (?)
_D1EN8DAZ DQ	01H DUP (?)
_D2T1DBZ DQ	01H DUP (?)
_DOV2DAZ DQ	01H DUP (?)
_D1EN8DBZ DQ	01H DUP (?)
_T3	DQ	01H DUP (?)
_DOV2DBZ DQ	01H DUP (?)
_DSDBZ	DQ	00H DUP (?)
_EN9	DQ	01H DUP (?)
_EN1	DQ	01H DUP (?)
_D3T1DA1 DQ	01H DUP (?)
_OV3	DQ	01H DUP (?)
_D2EN1DA1 DQ	01H DUP (?)
_D1EN1DA1 DQ	01H DUP (?)
_D3T1DA2 DQ	01H DUP (?)
_DOV3DA1 DQ	01H DUP (?)
_D2EN1DA2 DQ	01H DUP (?)
_D3T1DA5 DQ	01H DUP (?)
_DOV3DA2 DQ	01H DUP (?)
_GRADM	DQ	00H DUP (?)
_D2EN1DA5 DQ	01H DUP (?)
_D3T1DAZ DQ	01H DUP (?)
_DOV3DA5 DQ	01H DUP (?)
_D2EN1DAZ DQ	01H DUP (?)
_D3T1DBZ DQ	01H DUP (?)
_DOV3DAZ DQ	01H DUP (?)
_D2EN1DBZ DQ	01H DUP (?)
_D1EN1DA2 DQ	01H DUP (?)
_T4	DQ	01H DUP (?)
_DOV3DBZ DQ	01H DUP (?)
_DS	DQ	00H DUP (?)
_D4T1DA1 DQ	01H DUP (?)
_OV4	DQ	01H DUP (?)
_ER	DQ	01H DUP (?)
_ER1	DQ	01H DUP (?)
_D1EN1DA5 DQ	01H DUP (?)
_D4T1DA2 DQ	01H DUP (?)
_DOV4DA1 DQ	01H DUP (?)
_DER1DA1 DQ	01H DUP (?)
_EN10	DQ	01H DUP (?)
_D4T1DA5 DQ	01H DUP (?)
_DOV4DA2 DQ	01H DUP (?)
_DH	DQ	00H DUP (?)
_D2EN2DA1 DQ	01H DUP (?)
_D1EN1DAZ DQ	01H DUP (?)
_D4T1DAZ DQ	01H DUP (?)
_DOV4DA5 DQ	01H DUP (?)
_D2EN2DA2 DQ	01H DUP (?)
_D4T1DBZ DQ	01H DUP (?)
_DOV4DAZ DQ	01H DUP (?)
_D2EN2DA5 DQ	01H DUP (?)
_D1EN1DBZ DQ	01H DUP (?)
_T5	DQ	01H DUP (?)
_DOV4DBZ DQ	01H DUP (?)
_GRADA	DQ	00H DUP (?)
_BSS	ENDS
_DATA	SEGMENT
_7GRADC	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H
_7GRADA	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H
_7DH	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H
	ORG $+4
_7DS	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H
	ORG $+4
_7GRADM	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H
	ORG $+4
_7DSDBZ	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H
	ORG $+4
_7DSDAZ	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H
	ORG $+4
_7DSDA5	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H
	ORG $+4
_7DSDA2	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H
	ORG $+4
_7DSDA1	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H
	ORG $+4
_7DHDBZ	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H
	ORG $+4
_7DHDAZ	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H
	ORG $+4
_7DHDA5	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H
	ORG $+4
_7DHDA2	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H
	ORG $+4
_7DHDA1	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H
	ORG $+4
_7TEMP	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H, 00H
	DB	00H, 00H, 00H, 00H, 00H, 00H, 00H
_.SG139_ DB	'!NOT ENOUGH MEMORY SPACE AVAILABLE ON THIS BOX! ...ABORT'
	DB	'ING!'
_.SG596_ DB	'EN1='
_.SG633_ DB	'ER='
	ORG $+1
_.SG635_ DB	'EN='
	ORG $+1
$T98574	DB	'_HAMILTON@80', 00H
_DATA	ENDS
_CONST	SEGMENT
$T99039	DQ	04000000000000000r		; 2
$T99768	DQ	03fe6d6d5aa4ead2dr		; 0.713725
$T99889	DQ	00000000000000000r		; 0
$T99995	DQ	03ff0000000000000r		; 1
_CONST	ENDS
_TEXT	SEGMENT
; File C:\users\dgilmore\Graduate.School\Adamowicz.Research\H2ANL\HAMILTON.FOR
$T98496 = -1192
$T98432 = -680
$T98368 = -1072
_H$ = 56
$T98497 = -1200
$T98433 = -688
$T98369 = -1080
_7CO$ = -472
$T98498 = -1208
$T98434 = -696
$T98370 = -1088
_CO$ = 60
$T98499 = -1216
$T98435 = -704
$T98371 = -1096
_EMIN$ = 64
$T98500 = -1224
$T98436 = -712
$T98372 = -1104
$T98308 = -592
_7GRAD$ = -120
$T98501 = -1232
$T98437 = -720
$T98373 = -1112
$T98309 = -600
_GRAD$ = 68
$T98502 = -1240
$T98438 = -728
$T98374 = -1120
$T98310 = -608
_7TTT1$ = -504
$T98503 = -1248
$T98439 = -736
$T98375 = -1128
$T98311 = -616
_TTT1$ = 72
$T98504 = -1256
$T98440 = -744
$T98376 = -1136
$T98312 = -624
_J$ = -616
_7TTT2$ = -536
$T98505 = -1264
$T98441 = -752
$T98377 = -1144
$T98313 = -632
_TTT2$ = 76
$T98506 = -1272
$T98442 = -760
$T98378 = -1152
$T98314 = -640
_7TTT3$ = -568
$T98507 = -1280
$T98443 = -768
$T98379 = -1160
$T98315 = -648
_I$ = -616
_TTT3$ = 80
$T98508 = -1288
$T98444 = -776
$T98380 = -1168
$T98316 = -656
_DFLAG$ = 84
$T98573 = -584
$T98509 = -1296
$T98445 = -784
$T98381 = -1176
$T98317 = -664
$T98510 = -1304
$T98446 = -792
$T98382 = -1184
$T98318 = -672
$T98511 = -1312
$T98447 = -800
$T98383 = -1192
$T98319 = -680
$T98512 = -1320
$T98448 = -808
$T98384 = -1200
$T98320 = -688
$T98513 = -1328
$T98449 = -816
$T98385 = -1208
$T98321 = -696
$T98514 = -1336
$T98450 = -824
$T98386 = -1216
$T98322 = -704
$T98515 = -1344
$T98451 = -832
$T98387 = -1224
$T98323 = -712
$T98516 = -1352
$T98452 = -840
$T98388 = -1232
$T98324 = -720
$T98453 = -848
$T98389 = -1240
$T98325 = -728
_.T85_$ = -1368
$T98454 = -856
$T98390 = -1248
$T98326 = -736
$T98455 = -864
$T98391 = -1256
$T98327 = -744
$T98456 = -872
$T98392 = -1264
$T98328 = -752
$T98457 = -880
$T98393 = -1272
$T98329 = -760
$T98458 = -888
$T98394 = -1280
$T98330 = -768
$T98459 = -896
$T98395 = -1288
$T98331 = -776
$T98460 = -904
$T98396 = -1296
$T98332 = -784
$T98461 = -912
$T98397 = -1304
$T98333 = -792
$T98462 = -920
$T98398 = -1312
$T98334 = -800
$T98463 = -928
$T98399 = -1320
$T98335 = -808
$T98464 = -936
$T98400 = -1328
$T98336 = -816
$T98465 = -944
$T98401 = -1336
$T98337 = -824
$T98466 = -952
$T98402 = -1344
$T98338 = -832
$T98467 = -960
$T98403 = -1352
$T98339 = -840
$T98468 = -968
$T98340 = -848
$T98469 = -976
$T98341 = -856
$T98470 = -984
$T98342 = -864
$T98471 = -992
$T98343 = -872
$T98472 = -1000
$T98344 = -880
_M$ = 8
$T98473 = -1008
$T98345 = -888
_7AA1$ = -152
$T98474 = -1016
$T98346 = -896
_AA1$ = 12
$T98475 = -1024
$T98347 = -904
_7RX1$ = -184
$T98476 = -1032
$T98348 = -912
_RX1$ = 16
$T98477 = -1040
$T98349 = -920
_7RY1$ = -216
$T98478 = -1048
$T98350 = -928
_RY1$ = 20
$T98479 = -1056
$T98351 = -936
_7RZ1$ = -248
$T98480 = -1064
$T98352 = -944
_RZ1$ = 24
$T98481 = -1072
$T98353 = -952
_7AA2$ = -280
$T98482 = -1080
$T98354 = -960
_AA2$ = 28
$T98483 = -1088
$T98355 = -968
_7RX2$ = -312
$T98484 = -1096
$T98356 = -976
_RX2$ = 32
$T98485 = -1104
$T98421 = -592
$T98357 = -984
_7RY2$ = -344
$T98486 = -1112
$T98422 = -600
$T98358 = -992
_RY2$ = 36
$T98487 = -1120
$T98423 = -608
$T98359 = -1000
_7RZ2$ = -376
$T98488 = -1128
$T98424 = -616
$T98360 = -1008
_RZ2$ = 40
$T98489 = -1136
$T98425 = -624
$T98361 = -1016
_7B12$ = -408
$T98490 = -1144
$T98426 = -632
$T98362 = -1024
_B12$ = 44
$T98491 = -1152
$T98427 = -640
$T98363 = -1032
_7EVEC$ = -440
$T98492 = -1160
$T98428 = -648
$T98364 = -1040
_EVEC$ = 48
$T98493 = -1168
$T98429 = -656
$T98365 = -1048
_7S$ = -88
$T98494 = -1176
$T98430 = -664
$T98366 = -1056
_S$ = 52
$T98495 = -1184
$T98431 = -672
$T98367 = -1064
_7H$ = -44
_HAMILTON@80 PROC NEAR

; 3    :       SUBROUTINE HAMILTON(M,AA1,RX1,RY1,RZ1,AA2,RX2,RY2,RZ2,B12,

	sub	esp, 1376				; 00000560H
	lea	eax, DWORD PTR $T98573[esp+1376]
	push	ebx
	push	esi
	push	edi
	push	ebp
	mov	esi, 6
	push	eax
	call	___FHlink@4
	mov	DWORD PTR $T98573[esp+1396], OFFSET FLAT:_C:\users\dgilmore\Graduate.School\Adamowicz.Research\H2ANL\HAMILTON.FOR
	mov	DWORD PTR $T98573[esp+1400], OFFSET FLAT:$T98574
	mov	eax, DWORD PTR _M$[esp+1388]
	lea	ecx, DWORD PTR _7TTT3$[esp+1392]
	mov	DWORD PTR $T98573[esp+1404], esi

; 4    :      -                    EVEC,S,H,CO,EMIN,GRAD,TTT1,TTT2,TTT3,DFLAG)
; 5    : C     -----------------------------------------------------
; 6    :       IMPLICIT REAL*8 (A-H,O-Z)

	mov	edi, DWORD PTR [eax]
	push	edi
	push	1
	push	1
	push	ecx
	call	_make_info
	lea	eax, DWORD PTR _7TTT2$[esp+1408]
	add	esp, 16					; 00000010H
	push	edi
	push	1
	push	1
	push	eax
	call	_make_info
	lea	eax, DWORD PTR _7TTT1$[esp+1408]
	add	esp, 16					; 00000010H
	push	edi
	push	1
	push	1
	push	eax
	call	_make_info
	lea	eax, DWORD PTR _7CO$[esp+1408]
	add	esp, 16					; 00000010H
	push	edi
	push	1
	push	1
	push	eax
	call	_make_info
	lea	eax, DWORD PTR _7EVEC$[esp+1408]
	add	esp, 16					; 00000010H
	push	edi
	push	1
	push	1
	push	eax
	call	_make_info
	lea	eax, DWORD PTR _7B12$[esp+1408]
	add	esp, 16					; 00000010H
	push	edi
	push	1
	push	1
	push	eax
	call	_make_info
	lea	eax, DWORD PTR _7RZ2$[esp+1408]
	add	esp, 16					; 00000010H
	push	edi
	push	1
	push	1
	push	eax
	call	_make_info
	lea	eax, DWORD PTR _7RY2$[esp+1408]
	add	esp, 16					; 00000010H
	push	edi
	push	1
	push	1
	push	eax
	call	_make_info
	lea	eax, DWORD PTR _7RX2$[esp+1408]
	add	esp, 16					; 00000010H
	push	edi
	push	1
	push	1
	push	eax
	call	_make_info
	lea	eax, DWORD PTR _7AA2$[esp+1408]
	add	esp, 16					; 00000010H
	push	edi
	push	1
	push	1
	push	eax
	call	_make_info
	lea	eax, DWORD PTR _7RZ1$[esp+1408]
	add	esp, 16					; 00000010H
	push	edi
	push	1
	push	1
	push	eax
	call	_make_info
	lea	eax, DWORD PTR _7RY1$[esp+1408]
	add	esp, 16					; 00000010H
	push	edi
	push	1
	push	1
	push	eax
	call	_make_info
	add	esp, 16					; 00000010H
	push	edi
	push	1
	push	1
	lea	eax, DWORD PTR _7RX1$[esp+1404]
	push	eax
	call	_make_info
	lea	eax, DWORD PTR _7AA1$[esp+1408]
	add	esp, 16					; 00000010H
	push	edi
	push	1
	push	1
	push	eax
	call	_make_info
	mov	eax, DWORD PTR _M$[esp+1404]
	add	esp, 16					; 00000010H
	mov	edi, DWORD PTR [eax]
	imul	edi, esi
	jno	SHORT $L98576
	call	___FHovfl@0
$L98576:
	lea	eax, DWORD PTR _7GRAD$[esp+1392]
	push	edi
	push	1
	push	1
	push	eax
	call	_make_info
	mov	eax, DWORD PTR _M$[esp+1404]
	add	esp, 16					; 00000010H
	mov	esi, DWORD PTR [eax]
	push	esi
	mov	edi, esi
	push	1
	lea	eax, DWORD PTR _7H$[esp+1400]
	push	esi
	push	1
	push	2
	push	eax
	call	_make_info
	lea	eax, DWORD PTR _7S$[esp+1416]
	add	esp, 24					; 00000018H
	push	edi
	push	1
	push	esi
	push	1
	push	2
	push	eax
	mov	esi, 39					; 00000027H
	call	_make_info
	mov	DWORD PTR $T98573[esp+1428], 38		; 00000026H
	mov	ecx, DWORD PTR _M$[esp+1412]
	add	esp, 24					; 00000018H

; 7    :       INTEGER M,ERROR,DFLAG
; 8    :       REAL*8, DIMENSION(M) ::	AA1,RX1,RY1,RZ1,
; 9    : 	-						AA2,RX2,RY2,RZ2,
; 10   :      -						B12,EVEC,CO,TTT1,TTT2,TTT3
; 11   : 	REAL*8, DIMENSION(M*6) :: GRAD
; 12   : 	REAL*8, DIMENSION(M,M) :: S,H
; 13   : C
; 14   : 	REAL*8  TEMP[ALLOCATABLE](:)
; 15   : 
; 16   : 	REAL*8	DHDA1[ALLOCATABLE](:,:),DHDA2[ALLOCATABLE](:,:),
; 17   : 	-		DHDA5[ALLOCATABLE](:,:),DHDAZ[ALLOCATABLE](:,:),
; 18   :      -		DHDBZ[ALLOCATABLE](:,:)
; 19   : 	REAL*8	DSDA1[ALLOCATABLE](:,:),DSDA2[ALLOCATABLE](:,:),
; 20   : 	-		DSDA5[ALLOCATABLE](:,:),DSDAZ[ALLOCATABLE](:,:),
; 21   :      -		DSDBZ[ALLOCATABLE](:,:)
; 22   : 	REAL*8	GRADM[ALLOCATABLE](:,:),DS[ALLOCATABLE](:,:),
; 23   : 	-		DH[ALLOCATABLE](:,:),GRADA[ALLOCATABLE](:),
; 24   :      -		GRADC[ALLOCATABLE](:)
; 25   : C     -----------------------------------------------------
; 26   : 	COMMON /OVER/ L,K
; 27   : 	COMMON /TNOVER/ NXAC,XAC
; 28   : 	COMMON /DOV/ DOVDA1,DOVDA2,DOVDA5,DOVDAZ,DOVDBZ
; 29   : 	COMMON /DOVn/ DOVnDA1,DOVnDA2,DOVnDA5,DOVnDAZ,DOVnDBZ
; 30   : 	COMMON /DTK/ DTKDA1,DTKDA2,DTKDA5,DTKDAZ,DTKDBZ
; 31   : 	COMMON /DEN/ DENDA1,DENDA2,DENDA5,DENDAZ,DENDBZ,
; 32   : 	1			 DENDA1P,DENDA2P,DENDA5P,DENDAZP,DENDBZP
; 33   : 	COMMON /DER/ DERDA1,DERDA2,DERDA5,DERDAZ,DERDBZ
; 34   : C     -----------------------------------------------------
; 35   : 	PARAMETER(ISYMM=1)
; 36   : 	PARAMETER(IPERM=0)
; 37   : C     -----------------------------------------------------
; 38   :       MMAX=M

	mov	eax, DWORD PTR [ecx]
	mov	DWORD PTR $T98573[esp+1404], esi
	push	eax
	mov	DWORD PTR _MMAX, eax

; 39   :       ALLOCATE(	TEMP(MMAX),DHDA1(MMAX,MMAX),DHDA2(MMAX,MMAX),

	mov	DWORD PTR _.T85_$[esp+1396], 0
	push	1
	push	eax
	push	1
	push	2
	push	OFFSET FLAT:_7DHDA1
	call	_make_info
	add	esp, 24					; 00000018H
	mov	eax, DWORD PTR _DHDA1
	test	eax, eax
	je	SHORT $L88
	mov	DWORD PTR _.T85_$[esp+1392], 2
	jmp	SHORT $L86
$L88:
	mov	esi, 39					; 00000027H
	lea	eax, DWORD PTR _.T85_$[esp+1392]
	mov	DWORD PTR $T98573[esp+1404], esi
	push	eax
	push	8
	push	OFFSET FLAT:_7DHDA1
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _DHDA1, eax
$L86:
	mov	esi, 39					; 00000027H
	mov	eax, DWORD PTR _MMAX
	mov	DWORD PTR $T98573[esp+1404], esi
	push	eax
	push	1
	push	eax
	push	1
	push	2
	push	OFFSET FLAT:_7DHDA2
	call	_make_info
	add	esp, 24					; 00000018H
	mov	eax, DWORD PTR _DHDA2
	test	eax, eax
	je	SHORT $L92
	mov	DWORD PTR _.T85_$[esp+1392], 2
	jmp	SHORT $L90
$L92:
	mov	esi, 39					; 00000027H
	lea	eax, DWORD PTR _.T85_$[esp+1392]
	mov	DWORD PTR $T98573[esp+1404], esi
	push	eax
	push	8
	push	OFFSET FLAT:_7DHDA2
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _DHDA2, eax
$L90:
	mov	esi, 39					; 00000027H
	mov	eax, DWORD PTR _MMAX
	mov	DWORD PTR $T98573[esp+1404], esi
	push	eax
	push	1
	push	eax
	push	1
	push	2
	push	OFFSET FLAT:_7DHDA5
	call	_make_info
	add	esp, 24					; 00000018H
	mov	eax, DWORD PTR _DHDA5
	test	eax, eax
	je	SHORT $L95
	mov	DWORD PTR _.T85_$[esp+1392], 2
	jmp	SHORT $L93
$L95:
	mov	esi, 39					; 00000027H
	lea	eax, DWORD PTR _.T85_$[esp+1392]
	mov	DWORD PTR $T98573[esp+1404], esi
	push	eax
	push	8
	push	OFFSET FLAT:_7DHDA5
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _DHDA5, eax
$L93:
	mov	esi, 39					; 00000027H
	mov	eax, DWORD PTR _MMAX
	mov	DWORD PTR $T98573[esp+1404], esi
	push	eax
	push	1
	push	eax
	push	1
	push	2
	push	OFFSET FLAT:_7DHDAZ
	call	_make_info
	add	esp, 24					; 00000018H
	mov	eax, DWORD PTR _DHDAZ
	test	eax, eax
	je	SHORT $L98
	mov	DWORD PTR _.T85_$[esp+1392], 2
	jmp	SHORT $L96
$L98:
	mov	esi, 39					; 00000027H
	lea	eax, DWORD PTR _.T85_$[esp+1392]
	mov	DWORD PTR $T98573[esp+1404], esi
	push	eax
	push	8
	push	OFFSET FLAT:_7DHDAZ
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _DHDAZ, eax
$L96:
	mov	esi, 39					; 00000027H
	mov	eax, DWORD PTR _MMAX
	mov	DWORD PTR $T98573[esp+1404], esi
	push	eax
	push	1
	push	eax
	push	1
	push	2
	push	OFFSET FLAT:_7DHDBZ
	call	_make_info
	add	esp, 24					; 00000018H
	mov	eax, DWORD PTR _DHDBZ
	test	eax, eax
	je	SHORT $L101
	mov	DWORD PTR _.T85_$[esp+1392], 2
	jmp	SHORT $L99
$L101:
	mov	esi, 39					; 00000027H
	lea	eax, DWORD PTR _.T85_$[esp+1392]
	mov	DWORD PTR $T98573[esp+1404], esi
	push	eax
	push	8
	push	OFFSET FLAT:_7DHDBZ
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _DHDBZ, eax
$L99:
	mov	esi, 39					; 00000027H
	mov	eax, DWORD PTR _MMAX
	mov	DWORD PTR $T98573[esp+1404], esi
	push	eax
	push	1
	push	eax
	push	1
	push	2
	push	OFFSET FLAT:_7DSDA1
	call	_make_info
	add	esp, 24					; 00000018H
	mov	eax, DWORD PTR _DSDA1
	test	eax, eax
	je	SHORT $L104
	mov	DWORD PTR _.T85_$[esp+1392], 2
	jmp	SHORT $L102
$L104:
	mov	esi, 39					; 00000027H
	lea	eax, DWORD PTR _.T85_$[esp+1392]
	mov	DWORD PTR $T98573[esp+1404], esi
	push	eax
	push	8
	push	OFFSET FLAT:_7DSDA1
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _DSDA1, eax
$L102:
	mov	esi, 39					; 00000027H
	mov	eax, DWORD PTR _MMAX
	mov	DWORD PTR $T98573[esp+1404], esi
	push	eax
	push	1
	push	eax
	push	1
	push	2
	push	OFFSET FLAT:_7DSDA2
	call	_make_info
	add	esp, 24					; 00000018H
	mov	eax, DWORD PTR _DSDA2
	test	eax, eax
	je	SHORT $L107
	mov	DWORD PTR _.T85_$[esp+1392], 2
	jmp	SHORT $L105
$L107:
	mov	esi, 39					; 00000027H
	lea	eax, DWORD PTR _.T85_$[esp+1392]
	mov	DWORD PTR $T98573[esp+1404], esi
	push	eax
	push	8
	push	OFFSET FLAT:_7DSDA2
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _DSDA2, eax
$L105:
	mov	esi, 39					; 00000027H
	mov	eax, DWORD PTR _MMAX
	mov	DWORD PTR $T98573[esp+1404], esi
	push	eax
	push	1
	push	eax
	push	1
	push	2
	push	OFFSET FLAT:_7DSDA5
	call	_make_info
	add	esp, 24					; 00000018H
	mov	eax, DWORD PTR _DSDA5
	test	eax, eax
	je	SHORT $L110
	mov	DWORD PTR _.T85_$[esp+1392], 2
	jmp	SHORT $L108
$L110:
	mov	esi, 39					; 00000027H
	lea	eax, DWORD PTR _.T85_$[esp+1392]
	mov	DWORD PTR $T98573[esp+1404], esi
	push	eax
	push	8
	push	OFFSET FLAT:_7DSDA5
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _DSDA5, eax
$L108:
	mov	esi, 39					; 00000027H
	mov	eax, DWORD PTR _MMAX
	mov	DWORD PTR $T98573[esp+1404], esi
	push	eax
	push	1
	push	eax
	push	1
	push	2
	push	OFFSET FLAT:_7DSDAZ
	call	_make_info
	add	esp, 24					; 00000018H
	mov	eax, DWORD PTR _DSDAZ
	test	eax, eax
	je	SHORT $L113
	mov	DWORD PTR _.T85_$[esp+1392], 2
	jmp	SHORT $L111
$L113:
	mov	esi, 39					; 00000027H
	lea	eax, DWORD PTR _.T85_$[esp+1392]
	mov	DWORD PTR $T98573[esp+1404], esi
	push	eax
	push	8
	push	OFFSET FLAT:_7DSDAZ
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _DSDAZ, eax
$L111:
	mov	esi, 39					; 00000027H
	mov	eax, DWORD PTR _MMAX
	mov	DWORD PTR $T98573[esp+1404], esi
	push	eax
	push	1
	push	eax
	push	1
	push	2
	push	OFFSET FLAT:_7DSDBZ
	call	_make_info
	add	esp, 24					; 00000018H
	mov	eax, DWORD PTR _DSDBZ
	test	eax, eax
	je	SHORT $L116
	mov	DWORD PTR _.T85_$[esp+1392], 2
	jmp	SHORT $L114
$L116:
	mov	esi, 39					; 00000027H
	lea	eax, DWORD PTR _.T85_$[esp+1392]
	mov	DWORD PTR $T98573[esp+1404], esi
	push	eax
	push	8
	push	OFFSET FLAT:_7DSDBZ
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _DSDBZ, eax
$L114:
	mov	DWORD PTR $T98573[esp+1404], 39		; 00000027H
	mov	esi, DWORD PTR _MMAX
	mov	edi, esi
	imul	edi, 5
	jno	SHORT $L98578
	call	___FHovfl@0
$L98578:
	push	edi
	push	1
	push	esi
	push	1
	push	2
	push	OFFSET FLAT:_7GRADM
	call	_make_info
	add	esp, 24					; 00000018H
	mov	eax, DWORD PTR _GRADM
	test	eax, eax
	je	SHORT $L119
	mov	DWORD PTR _.T85_$[esp+1392], 2
	mov	eax, 39					; 00000027H
	mov	DWORD PTR $T98573[esp+1404], eax
	jmp	SHORT $L117
$L119:
	mov	esi, 39					; 00000027H
	lea	eax, DWORD PTR _.T85_$[esp+1392]
	mov	DWORD PTR $T98573[esp+1404], esi
	push	eax
	push	8
	push	OFFSET FLAT:_7GRADM
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _GRADM, eax
$L117:
	mov	DWORD PTR $T98573[esp+1404], 39		; 00000027H
	mov	esi, DWORD PTR _MMAX
	push	esi
	push	1
	imul	esi, 5
	jno	SHORT $L98579
	call	___FHovfl@0
$L98579:
	push	esi
	push	1
	push	2
	push	OFFSET FLAT:_7DS
	call	_make_info
	add	esp, 24					; 00000018H
	mov	eax, DWORD PTR _DS
	test	eax, eax
	je	SHORT $L122
	mov	DWORD PTR _.T85_$[esp+1392], 2
	mov	eax, 39					; 00000027H
	mov	DWORD PTR $T98573[esp+1404], eax
	jmp	SHORT $L120
$L122:
	mov	esi, 39					; 00000027H
	lea	eax, DWORD PTR _.T85_$[esp+1392]
	mov	DWORD PTR $T98573[esp+1404], esi
	push	eax
	push	8
	push	OFFSET FLAT:_7DS
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _DS, eax
$L120:
	mov	DWORD PTR $T98573[esp+1404], 39		; 00000027H
	mov	esi, DWORD PTR _MMAX
	push	esi
	push	1
	imul	esi, 5
	jno	SHORT $L98580
	call	___FHovfl@0
$L98580:
	push	esi
	push	1
	push	2
	push	OFFSET FLAT:_7DH
	call	_make_info
	add	esp, 24					; 00000018H
	mov	eax, DWORD PTR _DH
	test	eax, eax
	je	SHORT $L125
	mov	DWORD PTR _.T85_$[esp+1392], 2
	mov	eax, 39					; 00000027H
	mov	DWORD PTR $T98573[esp+1404], eax
	jmp	SHORT $L123
$L125:
	mov	esi, 39					; 00000027H
	lea	eax, DWORD PTR _.T85_$[esp+1392]
	mov	DWORD PTR $T98573[esp+1404], esi
	push	eax
	push	8
	push	OFFSET FLAT:_7DH
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _DH, eax
$L123:
	mov	DWORD PTR $T98573[esp+1404], 39		; 00000027H
	mov	esi, DWORD PTR _MMAX
	imul	esi, 5
	jno	SHORT $L98581
	call	___FHovfl@0
$L98581:
	push	esi
	push	1
	push	1
	push	OFFSET FLAT:_7GRADA
	call	_make_info
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _GRADA
	test	eax, eax
	je	SHORT $L128
	mov	DWORD PTR _.T85_$[esp+1392], 2
	mov	eax, 39					; 00000027H
	mov	DWORD PTR $T98573[esp+1404], eax
	jmp	SHORT $L126
$L128:
	mov	esi, 39					; 00000027H
	lea	eax, DWORD PTR _.T85_$[esp+1392]
	mov	DWORD PTR $T98573[esp+1404], esi
	push	eax
	push	8
	push	OFFSET FLAT:_7GRADA
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _GRADA, eax
$L126:
	mov	esi, 39					; 00000027H
	mov	eax, DWORD PTR _MMAX
	mov	DWORD PTR $T98573[esp+1404], esi
	push	eax
	push	1
	push	1
	push	OFFSET FLAT:_7GRADC
	call	_make_info
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _GRADC
	test	eax, eax
	je	SHORT $L131
	mov	DWORD PTR _.T85_$[esp+1392], 2
	jmp	SHORT $L129
$L131:
	mov	esi, 39					; 00000027H
	lea	eax, DWORD PTR _.T85_$[esp+1392]
	mov	DWORD PTR $T98573[esp+1404], esi
	push	eax
	push	8
	push	OFFSET FLAT:_7GRADC
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _GRADC, eax
$L129:
	mov	esi, 39					; 00000027H
	mov	eax, DWORD PTR _MMAX
	mov	DWORD PTR $T98573[esp+1404], esi
	push	eax
	push	1
	push	1
	push	OFFSET FLAT:_7TEMP
	call	_make_info
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _TEMP
	test	eax, eax
	je	SHORT $L134
	mov	DWORD PTR _.T85_$[esp+1392], 2
	jmp	SHORT $L132
$L134:
	mov	esi, 39					; 00000027H
	lea	eax, DWORD PTR _.T85_$[esp+1392]
	mov	DWORD PTR $T98573[esp+1404], esi
	push	eax
	push	8
	push	OFFSET FLAT:_7TEMP
	call	__Allocate_i
	add	esp, 12					; 0000000cH
	mov	DWORD PTR _TEMP, eax
$L132:
	mov	eax, DWORD PTR _.T85_$[esp+1392]
	mov	esi, 46					; 0000002eH
	mov	DWORD PTR _ERROR, eax
	test	eax, eax
	mov	DWORD PTR $T98573[esp+1404], 39		; 00000027H
	mov	DWORD PTR $T98573[esp+1404], esi

; 40   : 	-			DHDA5(MMAX,MMAX),DHDAZ(MMAX,MMAX),DHDBZ(MMAX,MMAX),
; 41   :      -			DSDA1(MMAX,MMAX),DSDA2(MMAX,MMAX),DSDA5(MMAX,MMAX),
; 42   :      -			DSDAZ(MMAX,MMAX),DSDBZ(MMAX,MMAX),GRADM(MMAX,MMAX*5),
; 43   :      -			DS(MMAX*5,MMAX),DH(MMAX*5,MMAX),GRADA(MMAX*5),
; 44   :      -			GRADC(MMAX),STAT=ERROR)
; 45   : C
; 46   : 	IF (ERROR .NE. 0) STOP '!NOT ENOUGH MEMORY SPACE AVAILABLE ON THIS

	je	SHORT $L136
	push	OFFSET FLAT:_.SG139_
	push	60					; 0000003cH
	call	__FFsSTOP
	add	esp, 8
$L136:
	mov	DWORD PTR $T98573[esp+1404], 51		; 00000033H

; 47   : 	- BOX! ...ABORTING!'
; 48   : C     -----------------------------------------------------
; 49   :       DATA ONE,TWO,FOUR,SIX,EIGHT,HALF/1.,2.,4.,6.,8.,.5/
; 50   : C     -----------------------------------------------------
; 51   :       CALL GEOM(EX1,EY1,EZ1,EX2,EY2,EZ2)

	push	OFFSET FLAT:_EZ2
	push	OFFSET FLAT:_EY2
	push	OFFSET FLAT:_EX2
	push	OFFSET FLAT:_EZ1
	push	OFFSET FLAT:_EY1
	push	OFFSET FLAT:_EX1
	call	_GEOM@24

; 52   : C     -----------------------------------------------------	  
; 53   : 	IF (DFLAG .EQ. 1) THEN

	mov	ecx, DWORD PTR _DFLAG$[esp+1388]
	mov	DWORD PTR $T98573[esp+1404], 53		; 00000035H
	cmp	DWORD PTR [ecx], 1
	jne	$L148
	mov	ecx, 54					; 00000036H
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	DWORD PTR _OVER, 1
	mov	edx, DWORD PTR [eax]
	mov	DWORD PTR $T98573[esp+1404], ecx

; 54   :       DO L=1,M

	mov	DWORD PTR -1364+[esp+1392], edx
	test	edx, edx
	jle	$L152
$L98306:
	mov	ecx, 55					; 00000037H
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	DWORD PTR _OVER+4, 1
	mov	eax, DWORD PTR [eax]
	mov	DWORD PTR $T98573[esp+1404], ecx

; 55   :           DO K=1,M

	mov	DWORD PTR -1372+[esp+1392], eax
	test	eax, eax
	jle	$L153
$L98307:
	mov	eax, DWORD PTR _M$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_OV1
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1412], 55		; 00000037H
	mov	DWORD PTR $T98573[esp+1412], 59		; 0000003bH
	mov	esi, DWORD PTR [eax]

; 56   : C     -----------------------------------------------------
; 57   : C     --- Calculate OVERLAP INTEGRAL ---
; 58   : C     -----------------------------------------------------
; 59   :              CALL OVERLP(M,DFLAG,RX1(L),RY1(L),RZ1(L),AA1(L),

	cmp	edi, esi
	jg	SHORT $L98582
	cmp	edi, 1
	jge	SHORT $L98583
$L98582:
	call	___FHranger@0
$L98583:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98585
	cmp	edi, 1
	jge	SHORT $L98586
$L98585:
	call	___FHranger@0
$L98586:
	mov	eax, DWORD PTR _B12$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98587
	cmp	edi, 1
	jge	SHORT $L98588
$L98587:
	call	___FHranger@0
$L98588:
	mov	eax, DWORD PTR _AA2$[esp+1404]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98589
	cmp	edi, 1
	jge	SHORT $L98590
$L98589:
	call	___FHranger@0
$L98590:
	mov	eax, DWORD PTR _RZ2$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98591
	cmp	edi, 1
	jge	SHORT $L98592
$L98591:
	call	___FHranger@0
$L98592:
	mov	eax, DWORD PTR _RY2$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98593
	cmp	edi, 1
	jge	SHORT $L98594
$L98593:
	call	___FHranger@0
$L98594:
	mov	eax, DWORD PTR _RX2$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98595
	cmp	edi, 1
	jge	SHORT $L98596
$L98595:
	call	___FHranger@0
$L98596:
	mov	eax, DWORD PTR _AA2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98597
	cmp	edi, 1
	jge	SHORT $L98598
$L98597:
	call	___FHranger@0
$L98598:
	mov	eax, DWORD PTR _RZ2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98599
	cmp	edi, 1
	jge	SHORT $L98600
$L98599:
	call	___FHranger@0
$L98600:
	mov	eax, DWORD PTR _RY2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98601
	cmp	edi, 1
	jge	SHORT $L98602
$L98601:
	call	___FHranger@0
$L98602:
	mov	eax, DWORD PTR _RX2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98603
	cmp	edi, 1
	jge	SHORT $L98604
$L98603:
	call	___FHranger@0
$L98604:
	mov	eax, DWORD PTR _AA1$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98605
	cmp	edi, 1
	jge	SHORT $L98606
$L98605:
	call	___FHranger@0
$L98606:
	mov	eax, DWORD PTR _RZ1$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98607
	cmp	edi, 1
	jge	SHORT $L98608
$L98607:
	call	___FHranger@0
$L98608:
	mov	eax, DWORD PTR _RY1$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98609
	cmp	edi, 1
	jge	SHORT $L98610
$L98609:
	call	___FHranger@0
$L98610:
	mov	eax, DWORD PTR _RX1$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98611
	cmp	edi, 1
	jge	SHORT $L98612
$L98611:
	call	___FHranger@0
$L98612:
	mov	eax, DWORD PTR _AA1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98613
	cmp	edi, 1
	jge	SHORT $L98614
$L98613:
	call	___FHranger@0
$L98614:
	mov	eax, DWORD PTR _RZ1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98615
	cmp	edi, 1
	jge	SHORT $L98616
$L98615:
	call	___FHranger@0
$L98616:
	mov	eax, DWORD PTR _RY1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98617
	cmp	edi, 1
	jge	SHORT $L98618
$L98617:
	call	___FHranger@0
$L98618:
	mov	eax, DWORD PTR _RX1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1464]
	sub	edi, 8
	push	edi
	push	eax
	mov	eax, DWORD PTR _M$[esp+1472]
	push	eax
	call	_OVERLP@88
	mov	eax, DWORD PTR _DOV+4
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_OV2
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1412], 65		; 00000041H

; 60   :      1            RX1(K),RY1(K),RZ1(K),AA1(K),                      
; 61   :      2            RX2(L),RY2(L),RZ2(L),AA2(L),
; 62   :      3            RX2(K),RY2(K),RZ2(K),AA2(K),
; 63   :      4            B12(L),B12(K),OV1,KLM)
; 64   : C
; 65   : 	DOV1DA1=DOVDA1

	mov	DWORD PTR $T98573[esp+1412], 66		; 00000042H
	mov	DWORD PTR $T98573[esp+1412], 67		; 00000043H
	mov	DWORD PTR $T98573[esp+1412], 68		; 00000044H
	mov	DWORD PTR $T98573[esp+1412], 69		; 00000045H
	mov	DWORD PTR $T98573[esp+1412], 72		; 00000048H
	mov	DWORD PTR _DOV1DA1+4, eax
	mov	eax, DWORD PTR _DOV
	mov	DWORD PTR _DOV1DA1, eax

; 66   : 	DOV1DA2=DOVDA2

	mov	eax, DWORD PTR _DOV+12
	mov	DWORD PTR _DOV1DA2+4, eax
	mov	eax, DWORD PTR _DOV+8
	mov	DWORD PTR _DOV1DA2, eax

; 67   : 	DOV1DA5=DOVDA5

	mov	eax, DWORD PTR _DOV+20
	mov	DWORD PTR _DOV1DA5+4, eax
	mov	eax, DWORD PTR _DOV+16
	mov	DWORD PTR _DOV1DA5, eax

; 68   : 	DOV1DAZ=DOVDAZ

	mov	eax, DWORD PTR _DOV+28
	mov	DWORD PTR _DOV1DAZ+4, eax
	mov	eax, DWORD PTR _DOV+24
	mov	DWORD PTR _DOV1DAZ, eax

; 69   : 	DOV1DBZ=DOVDBZ

	mov	eax, DWORD PTR _DOV+36
	mov	DWORD PTR _DOV1DBZ+4, eax
	mov	eax, DWORD PTR _DOV+32
	mov	DWORD PTR _DOV1DBZ, eax
	mov	eax, DWORD PTR _M$[esp+1396]
	mov	esi, DWORD PTR [eax]

; 70   : C
; 71   : 	IF(ISYMM) THEN
; 72   : 		CALL OVERLP(M,DFLAG,RX1(L),RY1(L),RZ1(L),AA1(L),

	cmp	edi, esi
	jg	SHORT $L98619
	cmp	edi, 1
	jge	SHORT $L98620
$L98619:
	call	___FHranger@0
$L98620:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98621
	cmp	edi, 1
	jge	SHORT $L98622
$L98621:
	call	___FHranger@0
$L98622:
	mov	eax, DWORD PTR _B12$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98623
	cmp	edi, 1
	jge	SHORT $L98624
$L98623:
	call	___FHranger@0
$L98624:
	mov	eax, DWORD PTR _AA1$[esp+1404]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98625
	cmp	edi, 1
	jge	SHORT $L98626
$L98625:
	call	___FHranger@0
$L98626:
	mov	eax, DWORD PTR _RZ1$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98627
	cmp	edi, 1
	jge	SHORT $L98628
$L98627:
	call	___FHranger@0
$L98628:
	mov	eax, DWORD PTR _RY1$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98629
	cmp	edi, 1
	jge	SHORT $L98630
$L98629:
	call	___FHranger@0
$L98630:
	mov	eax, DWORD PTR _RX1$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98631
	cmp	edi, 1
	jge	SHORT $L98632
$L98631:
	call	___FHranger@0
$L98632:
	mov	eax, DWORD PTR _AA2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98633
	cmp	edi, 1
	jge	SHORT $L98634
$L98633:
	call	___FHranger@0
$L98634:
	mov	eax, DWORD PTR _RZ2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98635
	cmp	edi, 1
	jge	SHORT $L98636
$L98635:
	call	___FHranger@0
$L98636:
	mov	eax, DWORD PTR _RY2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98637
	cmp	edi, 1
	jge	SHORT $L98638
$L98637:
	call	___FHranger@0
$L98638:
	mov	eax, DWORD PTR _RX2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98639
	cmp	edi, 1
	jge	SHORT $L98640
$L98639:
	call	___FHranger@0
$L98640:
	mov	eax, DWORD PTR _AA2$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98641
	cmp	edi, 1
	jge	SHORT $L98642
$L98641:
	call	___FHranger@0
$L98642:
	mov	eax, DWORD PTR _RZ2$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98643
	cmp	edi, 1
	jge	SHORT $L98644
$L98643:
	call	___FHranger@0
$L98644:
	mov	eax, DWORD PTR _RY2$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98645
	cmp	edi, 1
	jge	SHORT $L98646
$L98645:
	call	___FHranger@0
$L98646:
	mov	eax, DWORD PTR _RX2$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98647
	cmp	edi, 1
	jge	SHORT $L98648
$L98647:
	call	___FHranger@0
$L98648:
	mov	eax, DWORD PTR _AA1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98649
	cmp	edi, 1
	jge	SHORT $L98650
$L98649:
	call	___FHranger@0
$L98650:
	mov	eax, DWORD PTR _RZ1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98651
	cmp	edi, 1
	jge	SHORT $L98652
$L98651:
	call	___FHranger@0
$L98652:
	mov	eax, DWORD PTR _RY1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98653
	cmp	edi, 1
	jge	SHORT $L98654
$L98653:
	call	___FHranger@0
$L98654:
	mov	eax, DWORD PTR _RX1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1464]
	sub	edi, 8
	push	edi
	push	eax
	mov	eax, DWORD PTR _M$[esp+1472]
	push	eax
	call	_OVERLP@88
	mov	eax, DWORD PTR _DOV+4
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 78		; 0000004eH

; 73   :      1            RX2(K),RY2(K),RZ2(K),AA2(K),                      
; 74   :      2            RX2(L),RY2(L),RZ2(L),AA2(L),
; 75   :      3            RX1(K),RY1(K),RZ1(K),AA1(K),
; 76   :      4            B12(L),B12(K),OV2,KLM)
; 77   : C
; 78   : 	DOV2DA1=DOVDA1

	mov	DWORD PTR $T98573[esp+1404], 79		; 0000004fH
	mov	DWORD PTR $T98573[esp+1404], 80		; 00000050H
	mov	DWORD PTR $T98573[esp+1404], 81		; 00000051H
	mov	DWORD PTR $T98573[esp+1404], 82		; 00000052H
	mov	DWORD PTR $T98573[esp+1404], 86		; 00000056H
	mov	DWORD PTR _DOV2DA1+4, eax
	mov	eax, DWORD PTR _DOV
	mov	DWORD PTR _DOV2DA1, eax

; 79   : 	DOV2DA2=DOVDA2

	mov	eax, DWORD PTR _DOV+12
	mov	DWORD PTR _DOV2DA2+4, eax
	mov	eax, DWORD PTR _DOV+8
	mov	DWORD PTR _DOV2DA2, eax

; 80   : 	DOV2DA5=DOVDA5

	mov	eax, DWORD PTR _DOV+20
	mov	DWORD PTR _DOV2DA5+4, eax
	mov	eax, DWORD PTR _DOV+16
	mov	DWORD PTR _DOV2DA5, eax

; 81   : 	DOV2DAZ=DOVDAZ

	mov	eax, DWORD PTR _DOV+28
	mov	DWORD PTR _DOV2DAZ+4, eax
	mov	eax, DWORD PTR _DOV+24
	mov	DWORD PTR _DOV2DAZ, eax

; 82   : 	DOV2DBZ=DOVDBZ

	mov	eax, DWORD PTR _DOV+36
	mov	DWORD PTR _DOV2DBZ+4, eax
	mov	eax, DWORD PTR _DOV+32
	mov	DWORD PTR _DOV2DBZ, eax
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	esi, DWORD PTR [eax]

; 83   : C
; 84   : C     ...Nuclear Symmetrization Overlap Integral: P_[H1,H2]
; 85   : C
; 86   :               CALL OVERLP(M,DFLAG,RX1(L),RY1(L),RZ1(L),AA1(L),

	cmp	edi, esi
	jg	SHORT $L98655
	cmp	edi, 1
	jge	SHORT $L98656
$L98655:
	call	___FHranger@0
$L98656:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98308[esp+1392]
	cmp	edi, esi
	jg	SHORT $L98657
	cmp	edi, 1
	jge	SHORT $L98658
$L98657:
	call	___FHranger@0
$L98658:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98309[esp+1392]
	cmp	edi, esi
	jg	SHORT $L98659
	cmp	edi, 1
	jge	SHORT $L98660
$L98659:
	call	___FHranger@0
$L98660:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98310[esp+1392]
	cmp	edi, esi
	jg	SHORT $L98661
	cmp	edi, 1
	jge	SHORT $L98662
$L98661:
	call	___FHranger@0
$L98662:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98311[esp+1392]
	cmp	edi, esi
	jg	SHORT $L98663
	cmp	edi, 1
	jge	SHORT $L98664
$L98663:
	call	___FHranger@0
$L98664:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98312[esp+1392]
	cmp	edi, esi
	jg	SHORT $L98665
	cmp	edi, 1
	jge	SHORT $L98666
$L98665:
	call	___FHranger@0
$L98666:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_OV3
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98313[esp+1400]
	cmp	edi, esi
	jg	SHORT $L98667
	cmp	edi, 1
	jge	SHORT $L98668
$L98667:
	call	___FHranger@0
$L98668:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98669
	cmp	edi, 1
	jge	SHORT $L98670
$L98669:
	call	___FHranger@0
$L98670:
	mov	eax, DWORD PTR _B12$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98671
	cmp	edi, 1
	jge	SHORT $L98672
$L98671:
	call	___FHranger@0
$L98672:
	mov	eax, DWORD PTR _AA2$[esp+1404]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98308[esp+1408]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98309[esp+1416]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98310[esp+1420]
	push	eax
	jg	SHORT $L98673
	cmp	edi, 1
	jge	SHORT $L98674
$L98673:
	call	___FHranger@0
$L98674:
	mov	eax, DWORD PTR _AA2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98675
	cmp	edi, 1
	jge	SHORT $L98676
$L98675:
	call	___FHranger@0
$L98676:
	mov	eax, DWORD PTR _RZ2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98677
	cmp	edi, 1
	jge	SHORT $L98678
$L98677:
	call	___FHranger@0
$L98678:
	mov	eax, DWORD PTR _RY2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98679
	cmp	edi, 1
	jge	SHORT $L98680
$L98679:
	call	___FHranger@0
$L98680:
	mov	eax, DWORD PTR _RX2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98681
	cmp	edi, 1
	jge	SHORT $L98682
$L98681:
	call	___FHranger@0
$L98682:
	mov	eax, DWORD PTR _AA1$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98311[esp+1440]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98312[esp+1448]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98313[esp+1452]
	push	eax
	jg	SHORT $L98683
	cmp	edi, 1
	jge	SHORT $L98684
$L98683:
	call	___FHranger@0
$L98684:
	mov	eax, DWORD PTR _AA1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98685
	cmp	edi, 1
	jge	SHORT $L98686
$L98685:
	call	___FHranger@0
$L98686:
	mov	eax, DWORD PTR _RZ1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98687
	cmp	edi, 1
	jge	SHORT $L98688
$L98687:
	call	___FHranger@0
$L98688:
	mov	eax, DWORD PTR _RY1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98689
	cmp	edi, 1
	jge	SHORT $L98690
$L98689:
	call	___FHranger@0
$L98690:
	mov	eax, DWORD PTR _RX1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1464]
	sub	edi, 8
	push	edi
	push	eax
	mov	eax, DWORD PTR _M$[esp+1472]
	push	eax
	call	_OVERLP@88
	mov	eax, DWORD PTR _DOV+4
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 92		; 0000005cH

; 87   :      1            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
; 88   :      2            RX2(L),RY2(L),RZ2(L),AA2(L),
; 89   :      3            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
; 90   :      4            B12(L),B12(K),OV3,KLM)
; 91   : C
; 92   : 	DOV3DA1=DOVDA1

	mov	DWORD PTR $T98573[esp+1404], 93		; 0000005dH
	mov	DWORD PTR $T98573[esp+1404], 94		; 0000005eH
	mov	DWORD PTR $T98573[esp+1404], 95		; 0000005fH
	mov	DWORD PTR $T98573[esp+1404], 96		; 00000060H
	mov	DWORD PTR $T98573[esp+1404], 101	; 00000065H
	mov	DWORD PTR _DOV3DA1+4, eax
	mov	eax, DWORD PTR _DOV
	mov	DWORD PTR _DOV3DA1, eax

; 93   : 	DOV3DA2=DOVDA2

	mov	eax, DWORD PTR _DOV+12
	mov	DWORD PTR _DOV3DA2+4, eax
	mov	eax, DWORD PTR _DOV+8
	mov	DWORD PTR _DOV3DA2, eax

; 94   : 	DOV3DA5=DOVDA5

	mov	eax, DWORD PTR _DOV+20
	mov	DWORD PTR _DOV3DA5+4, eax
	mov	eax, DWORD PTR _DOV+16
	mov	DWORD PTR _DOV3DA5, eax

; 95   : 	DOV3DAZ=DOVDAZ

	mov	eax, DWORD PTR _DOV+28
	mov	DWORD PTR _DOV3DAZ+4, eax
	mov	eax, DWORD PTR _DOV+24
	mov	DWORD PTR _DOV3DAZ, eax

; 96   : 	DOV3DBZ=DOVDBZ

	mov	eax, DWORD PTR _DOV+36
	mov	DWORD PTR _DOV3DBZ+4, eax
	mov	eax, DWORD PTR _DOV+32
	mov	DWORD PTR _DOV3DBZ, eax
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	esi, DWORD PTR [eax]

; 97   : C
; 98   : C     ...Electro/Nuclear Symmetrization Overlap Integral:
; 99   : C         P_[elec.1,elec.2]*P_[H1,H2]
; 100  : C
; 101  :               CALL OVERLP(M,DFLAG,RX1(L),RY1(L),RZ1(L),AA1(L),

	cmp	edi, esi
	jg	SHORT $L98691
	cmp	edi, 1
	jge	SHORT $L98692
$L98691:
	call	___FHranger@0
$L98692:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98314[esp+1392]
	cmp	edi, esi
	jg	SHORT $L98693
	cmp	edi, 1
	jge	SHORT $L98694
$L98693:
	call	___FHranger@0
$L98694:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98315[esp+1392]
	cmp	edi, esi
	jg	SHORT $L98695
	cmp	edi, 1
	jge	SHORT $L98696
$L98695:
	call	___FHranger@0
$L98696:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98316[esp+1392]
	cmp	edi, esi
	jg	SHORT $L98697
	cmp	edi, 1
	jge	SHORT $L98698
$L98697:
	call	___FHranger@0
$L98698:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98317[esp+1392]
	cmp	edi, esi
	jg	SHORT $L98699
	cmp	edi, 1
	jge	SHORT $L98700
$L98699:
	call	___FHranger@0
$L98700:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98318[esp+1392]
	cmp	edi, esi
	jg	SHORT $L98701
	cmp	edi, 1
	jge	SHORT $L98702
$L98701:
	call	___FHranger@0
$L98702:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_OV4
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98319[esp+1400]
	cmp	edi, esi
	jg	SHORT $L98703
	cmp	edi, 1
	jge	SHORT $L98704
$L98703:
	call	___FHranger@0
$L98704:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98705
	cmp	edi, 1
	jge	SHORT $L98706
$L98705:
	call	___FHranger@0
$L98706:
	mov	eax, DWORD PTR _B12$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98707
	cmp	edi, 1
	jge	SHORT $L98708
$L98707:
	call	___FHranger@0
$L98708:
	mov	eax, DWORD PTR _AA1$[esp+1404]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98314[esp+1408]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98315[esp+1416]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98316[esp+1420]
	push	eax
	jg	SHORT $L98709
	cmp	edi, 1
	jge	SHORT $L98710
$L98709:
	call	___FHranger@0
$L98710:
	mov	eax, DWORD PTR _AA2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98711
	cmp	edi, 1
	jge	SHORT $L98712
$L98711:
	call	___FHranger@0
$L98712:
	mov	eax, DWORD PTR _RZ2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98713
	cmp	edi, 1
	jge	SHORT $L98714
$L98713:
	call	___FHranger@0
$L98714:
	mov	eax, DWORD PTR _RY2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98715
	cmp	edi, 1
	jge	SHORT $L98716
$L98715:
	call	___FHranger@0
$L98716:
	mov	eax, DWORD PTR _RX2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98717
	cmp	edi, 1
	jge	SHORT $L98718
$L98717:
	call	___FHranger@0
$L98718:
	mov	eax, DWORD PTR _AA2$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98317[esp+1440]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98318[esp+1448]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98319[esp+1452]
	push	eax
	jg	SHORT $L98719
	cmp	edi, 1
	jge	SHORT $L98720
$L98719:
	call	___FHranger@0
$L98720:
	mov	eax, DWORD PTR _AA1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98721
	cmp	edi, 1
	jge	SHORT $L98722
$L98721:
	call	___FHranger@0
$L98722:
	mov	eax, DWORD PTR _RZ1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98723
	cmp	edi, 1
	jge	SHORT $L98724
$L98723:
	call	___FHranger@0
$L98724:
	mov	eax, DWORD PTR _RY1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98725
	cmp	edi, 1
	jge	SHORT $L98726
$L98725:
	call	___FHranger@0
$L98726:
	mov	eax, DWORD PTR _RX1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1464]
	sub	edi, 8
	push	edi
	push	eax
	mov	eax, DWORD PTR _M$[esp+1472]
	push	eax
	call	_OVERLP@88
	mov	eax, DWORD PTR _DOV+4
	mov	edi, DWORD PTR _OVER
	mov	DWORD PTR $T98573[esp+1404], 107	; 0000006bH

; 102  :      1            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
; 103  :      2            RX2(L),RY2(L),RZ2(L),AA2(L),
; 104  :      3            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
; 105  :      4            B12(L),B12(K),OV4,KLM)
; 106  : C              
; 107  : 	DOV4DA1=DOVDA1

	mov	DWORD PTR $T98573[esp+1404], 108	; 0000006cH
	mov	DWORD PTR $T98573[esp+1404], 109	; 0000006dH
	mov	DWORD PTR $T98573[esp+1404], 110	; 0000006eH
	mov	DWORD PTR $T98573[esp+1404], 111	; 0000006fH
	mov	DWORD PTR $T98573[esp+1404], 115	; 00000073H
	mov	DWORD PTR _DOV4DA1+4, eax
	mov	eax, DWORD PTR _DOV
	mov	DWORD PTR _DOV4DA1, eax

; 108  : 	DOV4DA2=DOVDA2

	mov	eax, DWORD PTR _DOV+12
	mov	DWORD PTR _DOV4DA2+4, eax
	mov	eax, DWORD PTR _DOV+8
	mov	DWORD PTR _DOV4DA2, eax

; 109  : 	DOV4DA5=DOVDA5

	mov	eax, DWORD PTR _DOV+20
	mov	DWORD PTR _DOV4DA5+4, eax
	mov	eax, DWORD PTR _DOV+16
	mov	DWORD PTR _DOV4DA5, eax

; 110  : 	DOV4DAZ=DOVDAZ

	mov	eax, DWORD PTR _DOV+28
	mov	DWORD PTR _DOV4DAZ+4, eax
	mov	eax, DWORD PTR _DOV+24
	mov	DWORD PTR _DOV4DAZ, eax

; 111  : 	DOV4DBZ=DOVDBZ

	mov	eax, DWORD PTR _DOV+36
	mov	DWORD PTR _DOV4DBZ+4, eax
	mov	eax, DWORD PTR _DOV+32
	mov	DWORD PTR _DOV4DBZ, eax
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	esi, DWORD PTR [eax]

; 112  : 	END IF
; 113  : C     -----------------------------------------------------
; 114  : 	IF(ISYMM) THEN
; 115  :               S(L,K)=OV1+OV2+OV3+OV4

	cmp	edi, esi
	jg	SHORT $L98727
	cmp	edi, 1
	jge	SHORT $L98728
$L98727:
	call	___FHranger@0
$L98728:
	mov	ebx, DWORD PTR _OVER+4
	cmp	ebx, esi
	jg	SHORT $L98729
	cmp	ebx, 1
	jge	SHORT $L98730
$L98729:
	call	___FHranger@0
$L98730:
	fld	QWORD PTR _OV3
	fadd	QWORD PTR _OV1
	mov	ecx, DWORD PTR _7S$[esp+1428]
	mov	eax, DWORD PTR _S$[esp+1388]
	mov	esi, DWORD PTR _OVER
	push	1
	fadd	QWORD PTR _OV2
	imul	ecx, ebx
	fadd	QWORD PTR _OV4
	add	ecx, DWORD PTR _7S$[esp+1404]
	push	OFFSET FLAT:_7DSDA1
	mov	DWORD PTR $T98573[esp+1412], 133	; 00000085H
	add	ecx, edi
	fstp	QWORD PTR [eax+ecx*8]

; 116  : 		ELSE
; 117  :               S(L,K)=OV1
; 118  : 	END IF
; 119  : 	IF(IPERM.EQ.1) THEN
; 120  : 		S(L,K)=OV1
; 121  : 	END IF
; 122  : 	IF(IPERM.EQ.2) THEN
; 123  : 		S(L,K)=OV2
; 124  : 	END IF
; 125  : 	IF(IPERM.EQ.3) THEN
; 126  : 		S(L,K)=OV3
; 127  : 	END IF
; 128  : 	IF(IPERM.EQ.4) THEN
; 129  : 		S(L,K)=OV4
; 130  : 	END IF
; 131  : C
; 132  : 	IF(ISYMM) THEN
; 133  : 	DSDA1(L,K)=DOV1DA1+DOV2DA1+DOV3DA1+DOV4DA1

	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L98731
	push	1
	push	OFFSET FLAT:_7DSDA1
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L98732
$L98731:
	call	___FHranger@0
$L98732:
	mov	edi, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DSDA1
	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L98733
	push	2
	push	OFFSET FLAT:_7DSDA1
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L98734
$L98733:
	call	___FHranger@0
$L98734:
	fld	QWORD PTR _DOV3DA1
	fadd	QWORD PTR _DOV1DA1
	imul	edi, DWORD PTR _7DSDA1+36
	fadd	QWORD PTR _DOV2DA1
	add	edi, esi
	mov	eax, DWORD PTR _DSDA1
	add	edi, DWORD PTR _7DSDA1+8
	mov	esi, DWORD PTR _OVER
	fadd	QWORD PTR _DOV4DA1
	mov	DWORD PTR $T98573[esp+1404], 134	; 00000086H
	push	1
	push	OFFSET FLAT:_7DSDA2
	fstp	QWORD PTR [eax+edi*8]

; 134  : 	DSDA2(L,K)=DOV1DA2+DOV2DA2+DOV3DA2+DOV4DA2

	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L98735
	push	1
	push	OFFSET FLAT:_7DSDA2
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L98736
$L98735:
	call	___FHranger@0
$L98736:
	mov	edi, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DSDA2
	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L98737
	push	2
	push	OFFSET FLAT:_7DSDA2
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L98738
$L98737:
	call	___FHranger@0
$L98738:
	fld	QWORD PTR _DOV3DA2
	fadd	QWORD PTR _DOV1DA2
	imul	edi, DWORD PTR _7DSDA2+36
	fadd	QWORD PTR _DOV2DA2
	add	edi, esi
	mov	eax, DWORD PTR _DSDA2
	add	edi, DWORD PTR _7DSDA2+8
	mov	esi, DWORD PTR _OVER
	fadd	QWORD PTR _DOV4DA2
	mov	DWORD PTR $T98573[esp+1404], 135	; 00000087H
	push	1
	push	OFFSET FLAT:_7DSDA5
	fstp	QWORD PTR [eax+edi*8]

; 135  : 	DSDA5(L,K)=DOV1DA5+DOV2DA5+DOV3DA5+DOV4DA5

	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L98739
	push	1
	push	OFFSET FLAT:_7DSDA5
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L98740
$L98739:
	call	___FHranger@0
$L98740:
	mov	edi, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DSDA5
	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L98741
	push	2
	push	OFFSET FLAT:_7DSDA5
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L98742
$L98741:
	call	___FHranger@0
$L98742:
	fld	QWORD PTR _DOV1DA5
	fadd	QWORD PTR _DOV3DA5
	mov	ecx, DWORD PTR _7DSDA5+36
	mov	eax, DWORD PTR _DSDA5
	mov	DWORD PTR $T98573[esp+1404], 136	; 00000088H
	fadd	QWORD PTR _DOV2DA5
	imul	ecx, edi
	fadd	QWORD PTR _DOV4DA5
	add	ecx, DWORD PTR _7DSDA5+8
	push	1
	add	ecx, esi
	push	OFFSET FLAT:_7DSDAZ
	mov	esi, DWORD PTR _OVER
	fstp	QWORD PTR [eax+ecx*8]

; 136  : 	DSDAZ(L,K)=DOV1DAZ+DOV2DAZ+DOV3DAZ+DOV4DAZ

	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L98743
	push	1
	push	OFFSET FLAT:_7DSDAZ
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L98744
$L98743:
	call	___FHranger@0
$L98744:
	mov	edi, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DSDAZ
	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L98745
	push	2
	push	OFFSET FLAT:_7DSDAZ
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L98746
$L98745:
	call	___FHranger@0
$L98746:
	fld	QWORD PTR _DOV3DAZ
	fadd	QWORD PTR _DOV2DAZ
	mov	ecx, DWORD PTR _7DSDAZ+36
	mov	eax, DWORD PTR _DSDAZ
	mov	DWORD PTR $T98573[esp+1404], 137	; 00000089H
	fadd	QWORD PTR _DOV4DAZ
	imul	ecx, edi
	fadd	QWORD PTR _DOV1DAZ
	add	ecx, DWORD PTR _7DSDAZ+8
	push	1
	add	ecx, esi
	push	OFFSET FLAT:_7DSDBZ
	mov	esi, DWORD PTR _OVER
	fstp	QWORD PTR [eax+ecx*8]

; 137  : 	DSDBZ(L,K)=DOV1DBZ+DOV2DBZ+DOV3DBZ+DOV4DBZ

	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L98747
	push	1
	push	OFFSET FLAT:_7DSDBZ
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L98748
$L98747:
	call	___FHranger@0
$L98748:
	mov	edi, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DSDBZ
	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L98749
	push	2
	push	OFFSET FLAT:_7DSDBZ
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L98750
$L98749:
	call	___FHranger@0
$L98750:
	fld	QWORD PTR _DOV2DBZ
	fadd	QWORD PTR _DOV4DBZ
	mov	ecx, DWORD PTR _7DSDBZ+36
	mov	eax, DWORD PTR _DSDBZ
	mov	DWORD PTR $T98573[esp+1404], 176	; 000000b0H
	fadd	QWORD PTR _DOV1DBZ
	imul	ecx, edi
	fadd	QWORD PTR _DOV3DBZ
	add	ecx, DWORD PTR _7DSDBZ+8
	push	OFFSET FLAT:_T1
	add	ecx, esi
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1408], 177	; 000000b1H
	fstp	QWORD PTR [eax+ecx*8]

; 138  : 	ELSE
; 139  : 	DSDA1(L,K)=DOV1DA1
; 140  : 	DSDA2(L,K)=DOV1DA2
; 141  : 	DSDA5(L,K)=DOV1DA5
; 142  : 	DSDAZ(L,K)=DOV1DAZ
; 143  : 	DSDBZ(L,K)=DOV1DBZ
; 144  : 	END IF		
; 145  : 	IF(IPERM.EQ.1) THEN
; 146  : 	DSDA1(L,K)=DOV1DA1
; 147  : 	DSDA2(L,K)=DOV1DA2
; 148  : 	DSDA5(L,K)=DOV1DA5
; 149  : 	DSDAZ(L,K)=DOV1DAZ
; 150  : 	DSDBZ(L,K)=DOV1DBZ
; 151  : 	END IF
; 152  : 	IF(IPERM.EQ.2) THEN
; 153  : 	DSDA1(L,K)=DOV2DA1
; 154  : 	DSDA2(L,K)=DOV2DA2
; 155  : 	DSDA5(L,K)=DOV2DA5
; 156  : 	DSDAZ(L,K)=DOV2DAZ
; 157  : 	DSDBZ(L,K)=DOV2DBZ
; 158  : 	END IF
; 159  : 	IF(IPERM.EQ.3) THEN
; 160  : 	DSDA1(L,K)=DOV3DA1
; 161  : 	DSDA2(L,K)=DOV3DA2
; 162  : 	DSDA5(L,K)=DOV3DA5
; 163  : 	DSDAZ(L,K)=DOV3DAZ
; 164  : 	DSDBZ(L,K)=DOV3DBZ
; 165  : 	END IF
; 166  : 	IF(IPERM.EQ.4) THEN
; 167  : 	DSDA1(L,K)=DOV4DA1
; 168  : 	DSDA2(L,K)=DOV4DA2
; 169  : 	DSDA5(L,K)=DOV4DA5
; 170  : 	DSDAZ(L,K)=DOV4DAZ
; 171  : 	DSDBZ(L,K)=DOV4DBZ
; 172  : 	END IF
; 173  : C     -----------------------------------------------------
; 174  : C     --- Calculate KINETIC ENERGY OF ELECTRON 1 ---
; 175  : C     -----------------------------------------------------
; 176  : 	DOVnDA1=DOV1DA1

	mov	DWORD PTR $T98573[esp+1408], 178	; 000000b2H
	mov	DWORD PTR $T98573[esp+1408], 179	; 000000b3H
	mov	DWORD PTR $T98573[esp+1408], 180	; 000000b4H
	mov	DWORD PTR $T98573[esp+1408], 181	; 000000b5H
	mov	eax, DWORD PTR _DOV1DA1+4
	mov	DWORD PTR _DOVN+4, eax
	mov	eax, DWORD PTR _DOV1DA1
	mov	DWORD PTR _DOVN, eax

; 177  : 	DOVnDA2=DOV1DA2

	mov	eax, DWORD PTR _DOV1DA2+4
	mov	DWORD PTR _DOVN+12, eax
	mov	eax, DWORD PTR _DOV1DA2
	mov	DWORD PTR _DOVN+8, eax

; 178  : 	DOVnDA5=DOV1DA5

	mov	eax, DWORD PTR _DOV1DA5+4
	mov	DWORD PTR _DOVN+20, eax
	mov	eax, DWORD PTR _DOV1DA5
	mov	DWORD PTR _DOVN+16, eax

; 179  : 	DOVnDAZ=DOV1DAZ

	mov	eax, DWORD PTR _DOV1DAZ+4
	mov	DWORD PTR _DOVN+28, eax
	mov	eax, DWORD PTR _DOV1DAZ
	mov	DWORD PTR _DOVN+24, eax

; 180  : 	DOVnDBZ=DOV1DBZ

	mov	eax, DWORD PTR _DOV1DBZ+4
	mov	DWORD PTR _DOVN+36, eax
	mov	eax, DWORD PTR _DOV1DBZ
	mov	DWORD PTR _DOVN+32, eax
	mov	eax, DWORD PTR _M$[esp+1392]
	mov	esi, DWORD PTR [eax]

; 181  :               CALL TK(DFLAG,OV1,

	cmp	edi, esi
	jg	SHORT $L98751
	cmp	edi, 1
	jge	SHORT $L98752
$L98751:
	call	___FHranger@0
$L98752:
	mov	eax, DWORD PTR _B12$[esp+1392]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98753
	cmp	edi, 1
	jge	SHORT $L98754
$L98753:
	call	___FHranger@0
$L98754:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98755
	cmp	edi, 1
	jge	SHORT $L98756
$L98755:
	call	___FHranger@0
$L98756:
	mov	eax, DWORD PTR _AA2$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98757
	cmp	edi, 1
	jge	SHORT $L98758
$L98757:
	call	___FHranger@0
$L98758:
	mov	eax, DWORD PTR _RZ2$[esp+1404]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98759
	cmp	edi, 1
	jge	SHORT $L98760
$L98759:
	call	___FHranger@0
$L98760:
	mov	eax, DWORD PTR _RY2$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98761
	cmp	edi, 1
	jge	SHORT $L98762
$L98761:
	call	___FHranger@0
$L98762:
	mov	eax, DWORD PTR _RX2$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98763
	cmp	edi, 1
	jge	SHORT $L98764
$L98763:
	call	___FHranger@0
$L98764:
	mov	eax, DWORD PTR _AA2$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98765
	cmp	edi, 1
	jge	SHORT $L98766
$L98765:
	call	___FHranger@0
$L98766:
	mov	eax, DWORD PTR _RZ2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98767
	cmp	edi, 1
	jge	SHORT $L98768
$L98767:
	call	___FHranger@0
$L98768:
	mov	eax, DWORD PTR _RY2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98769
	cmp	edi, 1
	jge	SHORT $L98770
$L98769:
	call	___FHranger@0
$L98770:
	mov	eax, DWORD PTR _RX2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98771
	cmp	edi, 1
	jge	SHORT $L98772
$L98771:
	call	___FHranger@0
$L98772:
	mov	eax, DWORD PTR _AA1$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98773
	cmp	edi, 1
	jge	SHORT $L98774
$L98773:
	call	___FHranger@0
$L98774:
	mov	eax, DWORD PTR _RZ1$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98775
	cmp	edi, 1
	jge	SHORT $L98776
$L98775:
	call	___FHranger@0
$L98776:
	mov	eax, DWORD PTR _RY1$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98777
	cmp	edi, 1
	jge	SHORT $L98778
$L98777:
	call	___FHranger@0
$L98778:
	mov	eax, DWORD PTR _RX1$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98779
	cmp	edi, 1
	jge	SHORT $L98780
$L98779:
	call	___FHranger@0
$L98780:
	mov	eax, DWORD PTR _AA1$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98781
	cmp	edi, 1
	jge	SHORT $L98782
$L98781:
	call	___FHranger@0
$L98782:
	mov	eax, DWORD PTR _RZ1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98783
	cmp	edi, 1
	jge	SHORT $L98784
$L98783:
	call	___FHranger@0
$L98784:
	mov	eax, DWORD PTR _RY1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98785
	cmp	edi, 1
	jge	SHORT $L98786
$L98785:
	call	___FHranger@0
$L98786:
	mov	eax, DWORD PTR _RX1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1460]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV1
	push	eax
	call	_TK@84
	mov	eax, DWORD PTR _DTK+4
	push	OFFSET FLAT:_T2
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR _D1T1DA1+4, eax
	mov	DWORD PTR $T98573[esp+1408], 187	; 000000bbH

; 182  :      1            RX1(L),RY1(L),RZ1(L),AA1(L),
; 183  :      2            RX1(K),RY1(K),RZ1(K),AA1(K),                      
; 184  :      3            RX2(L),RY2(L),RZ2(L),AA2(L),
; 185  :      4            RX2(K),RY2(K),RZ2(K),AA2(K),
; 186  :      5            B12(L),B12(K),T1)  
; 187  : 	D1T1DA1=DTKDA1

	mov	DWORD PTR $T98573[esp+1408], 188	; 000000bcH
	mov	DWORD PTR $T98573[esp+1408], 189	; 000000bdH
	mov	DWORD PTR $T98573[esp+1408], 190	; 000000beH
	mov	DWORD PTR $T98573[esp+1408], 191	; 000000bfH
	mov	DWORD PTR $T98573[esp+1408], 195	; 000000c3H
	mov	DWORD PTR $T98573[esp+1408], 196	; 000000c4H
	mov	DWORD PTR $T98573[esp+1408], 197	; 000000c5H
	mov	DWORD PTR $T98573[esp+1408], 198	; 000000c6H
	mov	DWORD PTR $T98573[esp+1408], 199	; 000000c7H
	mov	DWORD PTR $T98573[esp+1408], 200	; 000000c8H
	mov	eax, DWORD PTR _DTK
	mov	DWORD PTR _D1T1DA1, eax

; 188  : 	D1T1DA2=DTKDA2

	mov	eax, DWORD PTR _DTK+12
	mov	DWORD PTR _D1T1DA2+4, eax
	mov	eax, DWORD PTR _DTK+8
	mov	DWORD PTR _D1T1DA2, eax

; 189  : 	D1T1DA5=DTKDA5

	mov	eax, DWORD PTR _DTK+20
	mov	DWORD PTR _D1T1DA5+4, eax
	mov	eax, DWORD PTR _DTK+16
	mov	DWORD PTR _D1T1DA5, eax

; 190  : 	D1T1DAZ=DTKDAZ

	mov	eax, DWORD PTR _DTK+28
	mov	DWORD PTR _D1T1DAZ+4, eax
	mov	eax, DWORD PTR _DTK+24
	mov	DWORD PTR _D1T1DAZ, eax

; 191  : 	D1T1DBZ=DTKDBZ

	mov	eax, DWORD PTR _DTK+36
	mov	DWORD PTR _D1T1DBZ+4, eax
	mov	eax, DWORD PTR _DTK+32
	mov	DWORD PTR _D1T1DBZ, eax

; 192  : C
; 193  : C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
; 194  : 	IF (ISYMM) THEN
; 195  : 	DOVnDA1=DOV2DA1

	mov	eax, DWORD PTR _DOV2DA1+4
	mov	DWORD PTR _DOVN+4, eax
	mov	eax, DWORD PTR _DOV2DA1
	mov	DWORD PTR _DOVN, eax

; 196  : 	DOVnDA2=DOV2DA2

	mov	eax, DWORD PTR _DOV2DA2+4
	mov	DWORD PTR _DOVN+12, eax
	mov	eax, DWORD PTR _DOV2DA2
	mov	DWORD PTR _DOVN+8, eax

; 197  : 	DOVnDA5=DOV2DA5

	mov	eax, DWORD PTR _DOV2DA5+4
	mov	DWORD PTR _DOVN+20, eax
	mov	eax, DWORD PTR _DOV2DA5
	mov	DWORD PTR _DOVN+16, eax

; 198  : 	DOVnDAZ=DOV2DAZ

	mov	eax, DWORD PTR _DOV2DAZ+4
	mov	DWORD PTR _DOVN+28, eax
	mov	eax, DWORD PTR _DOV2DAZ
	mov	DWORD PTR _DOVN+24, eax

; 199  : 	DOVnDBZ=DOV2DBZ

	mov	eax, DWORD PTR _DOV2DBZ+4
	mov	DWORD PTR _DOVN+36, eax
	mov	eax, DWORD PTR _DOV2DBZ
	mov	DWORD PTR _DOVN+32, eax
	mov	eax, DWORD PTR _M$[esp+1392]
	mov	esi, DWORD PTR [eax]

; 200  :               CALL TK(DFLAG,OV2,

	cmp	edi, esi
	jg	SHORT $L98787
	cmp	edi, 1
	jge	SHORT $L98788
$L98787:
	call	___FHranger@0
$L98788:
	mov	eax, DWORD PTR _B12$[esp+1392]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98789
	cmp	edi, 1
	jge	SHORT $L98790
$L98789:
	call	___FHranger@0
$L98790:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98791
	cmp	edi, 1
	jge	SHORT $L98792
$L98791:
	call	___FHranger@0
$L98792:
	mov	eax, DWORD PTR _AA1$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98793
	cmp	edi, 1
	jge	SHORT $L98794
$L98793:
	call	___FHranger@0
$L98794:
	mov	eax, DWORD PTR _RZ1$[esp+1404]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98795
	cmp	edi, 1
	jge	SHORT $L98796
$L98795:
	call	___FHranger@0
$L98796:
	mov	eax, DWORD PTR _RY1$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98797
	cmp	edi, 1
	jge	SHORT $L98798
$L98797:
	call	___FHranger@0
$L98798:
	mov	eax, DWORD PTR _RX1$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98799
	cmp	edi, 1
	jge	SHORT $L98800
$L98799:
	call	___FHranger@0
$L98800:
	mov	eax, DWORD PTR _AA2$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98801
	cmp	edi, 1
	jge	SHORT $L98802
$L98801:
	call	___FHranger@0
$L98802:
	mov	eax, DWORD PTR _RZ2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98803
	cmp	edi, 1
	jge	SHORT $L98804
$L98803:
	call	___FHranger@0
$L98804:
	mov	eax, DWORD PTR _RY2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98805
	cmp	edi, 1
	jge	SHORT $L98806
$L98805:
	call	___FHranger@0
$L98806:
	mov	eax, DWORD PTR _RX2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98807
	cmp	edi, 1
	jge	SHORT $L98808
$L98807:
	call	___FHranger@0
$L98808:
	mov	eax, DWORD PTR _AA2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98809
	cmp	edi, 1
	jge	SHORT $L98810
$L98809:
	call	___FHranger@0
$L98810:
	mov	eax, DWORD PTR _RZ2$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98811
	cmp	edi, 1
	jge	SHORT $L98812
$L98811:
	call	___FHranger@0
$L98812:
	mov	eax, DWORD PTR _RY2$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98813
	cmp	edi, 1
	jge	SHORT $L98814
$L98813:
	call	___FHranger@0
$L98814:
	mov	eax, DWORD PTR _RX2$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98815
	cmp	edi, 1
	jge	SHORT $L98816
$L98815:
	call	___FHranger@0
$L98816:
	mov	eax, DWORD PTR _AA1$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98817
	cmp	edi, 1
	jge	SHORT $L98818
$L98817:
	call	___FHranger@0
$L98818:
	mov	eax, DWORD PTR _RZ1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98819
	cmp	edi, 1
	jge	SHORT $L98820
$L98819:
	call	___FHranger@0
$L98820:
	mov	eax, DWORD PTR _RY1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98821
	cmp	edi, 1
	jge	SHORT $L98822
$L98821:
	call	___FHranger@0
$L98822:
	mov	eax, DWORD PTR _RX1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1460]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV2
	push	eax
	call	_TK@84
	mov	eax, DWORD PTR _DTK+4
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 206	; 000000ceH

; 201  :      1            RX1(L),RY1(L),RZ1(L),AA1(L),
; 202  :      2            RX2(K),RY2(K),RZ2(K),AA2(K),                      
; 203  :      3            RX2(L),RY2(L),RZ2(L),AA2(L),
; 204  :      4            RX1(K),RY1(K),RZ1(K),AA1(K),
; 205  :      5            B12(L),B12(K),T2)  
; 206  : 	D2T1DA1=DTKDA1

	mov	DWORD PTR $T98573[esp+1404], 207	; 000000cfH
	mov	DWORD PTR $T98573[esp+1404], 208	; 000000d0H
	mov	DWORD PTR $T98573[esp+1404], 209	; 000000d1H
	mov	DWORD PTR $T98573[esp+1404], 210	; 000000d2H
	mov	DWORD PTR $T98573[esp+1404], 212	; 000000d4H
	mov	DWORD PTR $T98573[esp+1404], 213	; 000000d5H
	mov	DWORD PTR $T98573[esp+1404], 214	; 000000d6H
	mov	DWORD PTR $T98573[esp+1404], 215	; 000000d7H
	mov	DWORD PTR $T98573[esp+1404], 216	; 000000d8H
	mov	DWORD PTR $T98573[esp+1404], 217	; 000000d9H
	mov	DWORD PTR _D2T1DA1+4, eax
	mov	eax, DWORD PTR _DTK
	mov	DWORD PTR _D2T1DA1, eax

; 207  : 	D2T1DA2=DTKDA2

	mov	eax, DWORD PTR _DTK+12
	mov	DWORD PTR _D2T1DA2+4, eax
	mov	eax, DWORD PTR _DTK+8
	mov	DWORD PTR _D2T1DA2, eax

; 208  : 	D2T1DA5=DTKDA5

	mov	eax, DWORD PTR _DTK+20
	mov	DWORD PTR _D2T1DA5+4, eax
	mov	eax, DWORD PTR _DTK+16
	mov	DWORD PTR _D2T1DA5, eax

; 209  : 	D2T1DAZ=DTKDAZ

	mov	eax, DWORD PTR _DTK+28
	mov	DWORD PTR _D2T1DAZ+4, eax
	mov	eax, DWORD PTR _DTK+24
	mov	DWORD PTR _D2T1DAZ, eax

; 210  : 	D2T1DBZ=DTKDBZ

	mov	eax, DWORD PTR _DTK+36
	mov	DWORD PTR _D2T1DBZ+4, eax
	mov	eax, DWORD PTR _DTK+32
	mov	DWORD PTR _D2T1DBZ, eax

; 211  : C     ...Nuclear Symmetrization Integral: P_[H1,H2]
; 212  : 	DOVnDA1=DOV3DA1

	mov	eax, DWORD PTR _DOV3DA1+4
	mov	DWORD PTR _DOVN+4, eax
	mov	eax, DWORD PTR _DOV3DA1
	mov	DWORD PTR _DOVN, eax

; 213  : 	DOVnDA2=DOV3DA2

	mov	eax, DWORD PTR _DOV3DA2+4
	mov	DWORD PTR _DOVN+12, eax
	mov	eax, DWORD PTR _DOV3DA2
	mov	DWORD PTR _DOVN+8, eax

; 214  : 	DOVnDA5=DOV3DA5

	mov	eax, DWORD PTR _DOV3DA5+4
	mov	DWORD PTR _DOVN+20, eax
	mov	eax, DWORD PTR _DOV3DA5
	mov	DWORD PTR _DOVN+16, eax

; 215  : 	DOVnDAZ=DOV3DAZ

	mov	eax, DWORD PTR _DOV3DAZ+4
	mov	DWORD PTR _DOVN+28, eax
	mov	eax, DWORD PTR _DOV3DAZ
	mov	DWORD PTR _DOVN+24, eax

; 216  : 	DOVnDBZ=DOV3DBZ

	mov	eax, DWORD PTR _DOV3DBZ+4
	mov	DWORD PTR _DOVN+36, eax
	mov	eax, DWORD PTR _DOV3DBZ
	mov	DWORD PTR _DOVN+32, eax
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	esi, DWORD PTR [eax]

; 217  :               CALL TK(DFLAG,OV3,

	cmp	edi, esi
	jg	SHORT $L98823
	cmp	edi, 1
	jge	SHORT $L98824
$L98823:
	call	___FHranger@0
$L98824:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98320[esp+1392]
	cmp	edi, esi
	jg	SHORT $L98825
	cmp	edi, 1
	jge	SHORT $L98826
$L98825:
	call	___FHranger@0
$L98826:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98321[esp+1392]
	cmp	edi, esi
	jg	SHORT $L98827
	cmp	edi, 1
	jge	SHORT $L98828
$L98827:
	call	___FHranger@0
$L98828:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98322[esp+1392]
	cmp	edi, esi
	jg	SHORT $L98829
	cmp	edi, 1
	jge	SHORT $L98830
$L98829:
	call	___FHranger@0
$L98830:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98323[esp+1392]
	cmp	edi, esi
	jg	SHORT $L98831
	cmp	edi, 1
	jge	SHORT $L98832
$L98831:
	call	___FHranger@0
$L98832:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98324[esp+1392]
	cmp	edi, esi
	jg	SHORT $L98833
	cmp	edi, 1
	jge	SHORT $L98834
$L98833:
	call	___FHranger@0
$L98834:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	push	OFFSET FLAT:_T3
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98325[esp+1396]
	cmp	edi, esi
	jg	SHORT $L98835
	cmp	edi, 1
	jge	SHORT $L98836
$L98835:
	call	___FHranger@0
$L98836:
	mov	eax, DWORD PTR _B12$[esp+1392]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98837
	cmp	edi, 1
	jge	SHORT $L98838
$L98837:
	call	___FHranger@0
$L98838:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98839
	cmp	edi, 1
	jge	SHORT $L98840
$L98839:
	call	___FHranger@0
$L98840:
	mov	eax, DWORD PTR _AA2$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98320[esp+1404]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98321[esp+1412]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98322[esp+1416]
	push	eax
	jg	SHORT $L98841
	cmp	edi, 1
	jge	SHORT $L98842
$L98841:
	call	___FHranger@0
$L98842:
	mov	eax, DWORD PTR _AA2$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98843
	cmp	edi, 1
	jge	SHORT $L98844
$L98843:
	call	___FHranger@0
$L98844:
	mov	eax, DWORD PTR _RZ2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98845
	cmp	edi, 1
	jge	SHORT $L98846
$L98845:
	call	___FHranger@0
$L98846:
	mov	eax, DWORD PTR _RY2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98847
	cmp	edi, 1
	jge	SHORT $L98848
$L98847:
	call	___FHranger@0
$L98848:
	mov	eax, DWORD PTR _RX2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98849
	cmp	edi, 1
	jge	SHORT $L98850
$L98849:
	call	___FHranger@0
$L98850:
	mov	eax, DWORD PTR _AA1$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98323[esp+1436]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98324[esp+1444]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98325[esp+1448]
	push	eax
	jg	SHORT $L98851
	cmp	edi, 1
	jge	SHORT $L98852
$L98851:
	call	___FHranger@0
$L98852:
	mov	eax, DWORD PTR _AA1$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98853
	cmp	edi, 1
	jge	SHORT $L98854
$L98853:
	call	___FHranger@0
$L98854:
	mov	eax, DWORD PTR _RZ1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98855
	cmp	edi, 1
	jge	SHORT $L98856
$L98855:
	call	___FHranger@0
$L98856:
	mov	eax, DWORD PTR _RY1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98857
	cmp	edi, 1
	jge	SHORT $L98858
$L98857:
	call	___FHranger@0
$L98858:
	mov	eax, DWORD PTR _RX1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1460]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV3
	push	eax
	call	_TK@84
	mov	eax, DWORD PTR _DTK+4
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 223	; 000000dfH

; 218  :      1            RX1(L),RY1(L),RZ1(L),AA1(L),
; 219  :      2            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
; 220  :      3            RX2(L),RY2(L),RZ2(L),AA2(L),
; 221  :      4            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
; 222  :      5            B12(L),B12(K),T3)
; 223  : 	D3T1DA1=DTKDA1

	mov	DWORD PTR $T98573[esp+1404], 224	; 000000e0H
	mov	DWORD PTR $T98573[esp+1404], 225	; 000000e1H
	mov	DWORD PTR $T98573[esp+1404], 226	; 000000e2H
	mov	DWORD PTR $T98573[esp+1404], 227	; 000000e3H
	mov	DWORD PTR $T98573[esp+1404], 230	; 000000e6H
	mov	DWORD PTR $T98573[esp+1404], 231	; 000000e7H
	mov	DWORD PTR $T98573[esp+1404], 232	; 000000e8H
	mov	DWORD PTR $T98573[esp+1404], 233	; 000000e9H
	mov	DWORD PTR $T98573[esp+1404], 234	; 000000eaH
	mov	DWORD PTR $T98573[esp+1404], 235	; 000000ebH
	mov	DWORD PTR _D3T1DA1+4, eax
	mov	eax, DWORD PTR _DTK
	mov	DWORD PTR _D3T1DA1, eax

; 224  : 	D3T1DA2=DTKDA2

	mov	eax, DWORD PTR _DTK+12
	mov	DWORD PTR _D3T1DA2+4, eax
	mov	eax, DWORD PTR _DTK+8
	mov	DWORD PTR _D3T1DA2, eax

; 225  : 	D3T1DA5=DTKDA5

	mov	eax, DWORD PTR _DTK+20
	mov	DWORD PTR _D3T1DA5+4, eax
	mov	eax, DWORD PTR _DTK+16
	mov	DWORD PTR _D3T1DA5, eax

; 226  : 	D3T1DAZ=DTKDAZ

	mov	eax, DWORD PTR _DTK+28
	mov	DWORD PTR _D3T1DAZ+4, eax
	mov	eax, DWORD PTR _DTK+24
	mov	DWORD PTR _D3T1DAZ, eax

; 227  : 	D3T1DBZ=DTKDBZ

	mov	eax, DWORD PTR _DTK+36
	mov	DWORD PTR _D3T1DBZ+4, eax
	mov	eax, DWORD PTR _DTK+32
	mov	DWORD PTR _D3T1DBZ, eax

; 228  : C     ...Electro/Nuclear Symmetrization Integral:
; 229  : C         P_[elec.1,elec.2]*P_[H1,H2]
; 230  : 	DOVnDA1=DOV4DA1

	mov	eax, DWORD PTR _DOV4DA1+4
	mov	DWORD PTR _DOVN+4, eax
	mov	eax, DWORD PTR _DOV4DA1
	mov	DWORD PTR _DOVN, eax

; 231  : 	DOVnDA2=DOV4DA2

	mov	eax, DWORD PTR _DOV4DA2+4
	mov	DWORD PTR _DOVN+12, eax
	mov	eax, DWORD PTR _DOV4DA2
	mov	DWORD PTR _DOVN+8, eax

; 232  : 	DOVnDA5=DOV4DA5

	mov	eax, DWORD PTR _DOV4DA5+4
	mov	DWORD PTR _DOVN+20, eax
	mov	eax, DWORD PTR _DOV4DA5
	mov	DWORD PTR _DOVN+16, eax

; 233  : 	DOVnDAZ=DOV4DAZ

	mov	eax, DWORD PTR _DOV4DAZ+4
	mov	DWORD PTR _DOVN+28, eax
	mov	eax, DWORD PTR _DOV4DAZ
	mov	DWORD PTR _DOVN+24, eax

; 234  : 	DOVnDBZ=DOV4DBZ

	mov	eax, DWORD PTR _DOV4DBZ+4
	mov	DWORD PTR _DOVN+36, eax
	mov	eax, DWORD PTR _DOV4DBZ
	mov	DWORD PTR _DOVN+32, eax
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	esi, DWORD PTR [eax]

; 235  :               CALL TK(DFLAG,OV4,

	cmp	edi, esi
	jg	SHORT $L98859
	cmp	edi, 1
	jge	SHORT $L98860
$L98859:
	call	___FHranger@0
$L98860:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98326[esp+1392]
	cmp	edi, esi
	jg	SHORT $L98861
	cmp	edi, 1
	jge	SHORT $L98862
$L98861:
	call	___FHranger@0
$L98862:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98327[esp+1392]
	cmp	edi, esi
	jg	SHORT $L98863
	cmp	edi, 1
	jge	SHORT $L98864
$L98863:
	call	___FHranger@0
$L98864:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98328[esp+1392]
	cmp	edi, esi
	jg	SHORT $L98865
	cmp	edi, 1
	jge	SHORT $L98866
$L98865:
	call	___FHranger@0
$L98866:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98329[esp+1392]
	cmp	edi, esi
	jg	SHORT $L98867
	cmp	edi, 1
	jge	SHORT $L98868
$L98867:
	call	___FHranger@0
$L98868:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98330[esp+1392]
	cmp	edi, esi
	jg	SHORT $L98869
	cmp	edi, 1
	jge	SHORT $L98870
$L98869:
	call	___FHranger@0
$L98870:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	push	OFFSET FLAT:_T4
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98331[esp+1396]
	cmp	edi, esi
	jg	SHORT $L98871
	cmp	edi, 1
	jge	SHORT $L98872
$L98871:
	call	___FHranger@0
$L98872:
	mov	eax, DWORD PTR _B12$[esp+1392]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98873
	cmp	edi, 1
	jge	SHORT $L98874
$L98873:
	call	___FHranger@0
$L98874:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98875
	cmp	edi, 1
	jge	SHORT $L98876
$L98875:
	call	___FHranger@0
$L98876:
	mov	eax, DWORD PTR _AA1$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98326[esp+1404]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98327[esp+1412]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98328[esp+1416]
	push	eax
	jg	SHORT $L98877
	cmp	edi, 1
	jge	SHORT $L98878
$L98877:
	call	___FHranger@0
$L98878:
	mov	eax, DWORD PTR _AA2$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98879
	cmp	edi, 1
	jge	SHORT $L98880
$L98879:
	call	___FHranger@0
$L98880:
	mov	eax, DWORD PTR _RZ2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98881
	cmp	edi, 1
	jge	SHORT $L98882
$L98881:
	call	___FHranger@0
$L98882:
	mov	eax, DWORD PTR _RY2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98883
	cmp	edi, 1
	jge	SHORT $L98884
$L98883:
	call	___FHranger@0
$L98884:
	mov	eax, DWORD PTR _RX2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98885
	cmp	edi, 1
	jge	SHORT $L98886
$L98885:
	call	___FHranger@0
$L98886:
	mov	eax, DWORD PTR _AA2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98329[esp+1436]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98330[esp+1444]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98331[esp+1448]
	push	eax
	jg	SHORT $L98887
	cmp	edi, 1
	jge	SHORT $L98888
$L98887:
	call	___FHranger@0
$L98888:
	mov	eax, DWORD PTR _AA1$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98889
	cmp	edi, 1
	jge	SHORT $L98890
$L98889:
	call	___FHranger@0
$L98890:
	mov	eax, DWORD PTR _RZ1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98891
	cmp	edi, 1
	jge	SHORT $L98892
$L98891:
	call	___FHranger@0
$L98892:
	mov	eax, DWORD PTR _RY1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98893
	cmp	edi, 1
	jge	SHORT $L98894
$L98893:
	call	___FHranger@0
$L98894:
	mov	eax, DWORD PTR _RX1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1460]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV4
	push	eax
	call	_TK@84
	mov	eax, DWORD PTR _DTK+4
	push	OFFSET FLAT:_T5
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR _D4T1DA1+4, eax
	mov	DWORD PTR $T98573[esp+1408], 241	; 000000f1H

; 236  :      1            RX1(L),RY1(L),RZ1(L),AA1(L),
; 237  :      2            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
; 238  :      3            RX2(L),RY2(L),RZ2(L),AA2(L),
; 239  :      4            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
; 240  :      5            B12(L),B12(K),T4)
; 241  : 	D4T1DA1=DTKDA1

	mov	DWORD PTR $T98573[esp+1408], 242	; 000000f2H
	mov	DWORD PTR $T98573[esp+1408], 243	; 000000f3H
	mov	DWORD PTR $T98573[esp+1408], 244	; 000000f4H
	mov	DWORD PTR $T98573[esp+1408], 245	; 000000f5H
	mov	DWORD PTR $T98573[esp+1408], 250	; 000000faH
	mov	DWORD PTR $T98573[esp+1408], 251	; 000000fbH
	mov	DWORD PTR $T98573[esp+1408], 252	; 000000fcH
	mov	DWORD PTR $T98573[esp+1408], 253	; 000000fdH
	mov	DWORD PTR $T98573[esp+1408], 254	; 000000feH
	mov	DWORD PTR $T98573[esp+1408], 255	; 000000ffH
	mov	eax, DWORD PTR _DTK
	mov	DWORD PTR _D4T1DA1, eax

; 242  : 	D4T1DA2=DTKDA2

	mov	eax, DWORD PTR _DTK+12
	mov	DWORD PTR _D4T1DA2+4, eax
	mov	eax, DWORD PTR _DTK+8
	mov	DWORD PTR _D4T1DA2, eax

; 243  : 	D4T1DA5=DTKDA5

	mov	eax, DWORD PTR _DTK+20
	mov	DWORD PTR _D4T1DA5+4, eax
	mov	eax, DWORD PTR _DTK+16
	mov	DWORD PTR _D4T1DA5, eax

; 244  : 	D4T1DAZ=DTKDAZ

	mov	eax, DWORD PTR _DTK+28
	mov	DWORD PTR _D4T1DAZ+4, eax
	mov	eax, DWORD PTR _DTK+24
	mov	DWORD PTR _D4T1DAZ, eax

; 245  : 	D4T1DBZ=DTKDBZ

	mov	eax, DWORD PTR _DTK+36
	mov	DWORD PTR _D4T1DBZ+4, eax
	mov	eax, DWORD PTR _DTK+32
	mov	DWORD PTR _D4T1DBZ, eax

; 246  : 	END IF
; 247  : C     -----------------------------------------------------
; 248  : C     --- Calculate KINETIC ENERGY OF ELECTRON 2 ---
; 249  : C     -----------------------------------------------------
; 250  : 	DOVnDA1=DOV1DA2

	mov	eax, DWORD PTR _DOV1DA2+4
	mov	DWORD PTR _DOVN+4, eax
	mov	eax, DWORD PTR _DOV1DA2
	mov	DWORD PTR _DOVN, eax

; 251  : 	DOVnDA2=DOV1DA1

	mov	eax, DWORD PTR _DOV1DA1+4
	mov	DWORD PTR _DOVN+12, eax
	mov	eax, DWORD PTR _DOV1DA1
	mov	DWORD PTR _DOVN+8, eax

; 252  : 	DOVnDA5=DOV1DA5

	mov	eax, DWORD PTR _DOV1DA5+4
	mov	DWORD PTR _DOVN+20, eax
	mov	eax, DWORD PTR _DOV1DA5
	mov	DWORD PTR _DOVN+16, eax

; 253  : 	DOVnDAZ=DOV1DBZ

	mov	eax, DWORD PTR _DOV1DBZ+4
	mov	DWORD PTR _DOVN+28, eax
	mov	eax, DWORD PTR _DOV1DBZ
	mov	DWORD PTR _DOVN+24, eax

; 254  : 	DOVnDBZ=DOV1DAZ

	mov	eax, DWORD PTR _DOV1DAZ+4
	mov	DWORD PTR _DOVN+36, eax
	mov	eax, DWORD PTR _DOV1DAZ
	mov	DWORD PTR _DOVN+32, eax
	mov	eax, DWORD PTR _M$[esp+1392]
	mov	esi, DWORD PTR [eax]

; 255  :              CALL TK(DFLAG,OV1,

	cmp	edi, esi
	jg	SHORT $L98895
	cmp	edi, 1
	jge	SHORT $L98896
$L98895:
	call	___FHranger@0
$L98896:
	mov	eax, DWORD PTR _B12$[esp+1392]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98897
	cmp	edi, 1
	jge	SHORT $L98898
$L98897:
	call	___FHranger@0
$L98898:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98899
	cmp	edi, 1
	jge	SHORT $L98900
$L98899:
	call	___FHranger@0
$L98900:
	mov	eax, DWORD PTR _AA1$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98901
	cmp	edi, 1
	jge	SHORT $L98902
$L98901:
	call	___FHranger@0
$L98902:
	mov	eax, DWORD PTR _RZ1$[esp+1404]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98903
	cmp	edi, 1
	jge	SHORT $L98904
$L98903:
	call	___FHranger@0
$L98904:
	mov	eax, DWORD PTR _RY1$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98905
	cmp	edi, 1
	jge	SHORT $L98906
$L98905:
	call	___FHranger@0
$L98906:
	mov	eax, DWORD PTR _RX1$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98907
	cmp	edi, 1
	jge	SHORT $L98908
$L98907:
	call	___FHranger@0
$L98908:
	mov	eax, DWORD PTR _AA1$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98909
	cmp	edi, 1
	jge	SHORT $L98910
$L98909:
	call	___FHranger@0
$L98910:
	mov	eax, DWORD PTR _RZ1$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98911
	cmp	edi, 1
	jge	SHORT $L98912
$L98911:
	call	___FHranger@0
$L98912:
	mov	eax, DWORD PTR _RY1$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98913
	cmp	edi, 1
	jge	SHORT $L98914
$L98913:
	call	___FHranger@0
$L98914:
	mov	eax, DWORD PTR _RX1$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98915
	cmp	edi, 1
	jge	SHORT $L98916
$L98915:
	call	___FHranger@0
$L98916:
	mov	eax, DWORD PTR _AA2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98917
	cmp	edi, 1
	jge	SHORT $L98918
$L98917:
	call	___FHranger@0
$L98918:
	mov	eax, DWORD PTR _RZ2$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98919
	cmp	edi, 1
	jge	SHORT $L98920
$L98919:
	call	___FHranger@0
$L98920:
	mov	eax, DWORD PTR _RY2$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98921
	cmp	edi, 1
	jge	SHORT $L98922
$L98921:
	call	___FHranger@0
$L98922:
	mov	eax, DWORD PTR _RX2$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98923
	cmp	edi, 1
	jge	SHORT $L98924
$L98923:
	call	___FHranger@0
$L98924:
	mov	eax, DWORD PTR _AA2$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98925
	cmp	edi, 1
	jge	SHORT $L98926
$L98925:
	call	___FHranger@0
$L98926:
	mov	eax, DWORD PTR _RZ2$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98927
	cmp	edi, 1
	jge	SHORT $L98928
$L98927:
	call	___FHranger@0
$L98928:
	mov	eax, DWORD PTR _RY2$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98929
	cmp	edi, 1
	jge	SHORT $L98930
$L98929:
	call	___FHranger@0
$L98930:
	mov	eax, DWORD PTR _RX2$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1460]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV1
	push	eax
	call	_TK@84
	mov	eax, DWORD PTR _DTK+12
	push	OFFSET FLAT:_T6
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR _D1T2DA1+4, eax
	mov	DWORD PTR $T98573[esp+1408], 261	; 00000105H

; 256  :      1            RX2(L),RY2(L),RZ2(L),AA2(L),
; 257  :      2            RX2(K),RY2(K),RZ2(K),AA2(K),                      
; 258  :      3            RX1(L),RY1(L),RZ1(L),AA1(L),
; 259  :      4            RX1(K),RY1(K),RZ1(K),AA1(K),
; 260  :      5            B12(L),B12(K),T5)  
; 261  : 	D1T2DA1=DTKDA2

	mov	DWORD PTR $T98573[esp+1408], 262	; 00000106H
	mov	DWORD PTR $T98573[esp+1408], 263	; 00000107H
	mov	DWORD PTR $T98573[esp+1408], 264	; 00000108H
	mov	DWORD PTR $T98573[esp+1408], 265	; 00000109H
	mov	DWORD PTR $T98573[esp+1408], 268	; 0000010cH
	mov	DWORD PTR $T98573[esp+1408], 269	; 0000010dH
	mov	DWORD PTR $T98573[esp+1408], 270	; 0000010eH
	mov	DWORD PTR $T98573[esp+1408], 271	; 0000010fH
	mov	DWORD PTR $T98573[esp+1408], 272	; 00000110H
	mov	DWORD PTR $T98573[esp+1408], 273	; 00000111H
	mov	eax, DWORD PTR _DTK+8
	mov	DWORD PTR _D1T2DA1, eax

; 262  : 	D1T2DA2=DTKDA1

	mov	eax, DWORD PTR _DTK+4
	mov	DWORD PTR _D1T2DA2+4, eax
	mov	eax, DWORD PTR _DTK
	mov	DWORD PTR _D1T2DA2, eax

; 263  : 	D1T2DA5=DTKDA5

	mov	eax, DWORD PTR _DTK+20
	mov	DWORD PTR _D1T2DA5+4, eax
	mov	eax, DWORD PTR _DTK+16
	mov	DWORD PTR _D1T2DA5, eax

; 264  : 	D1T2DAZ=DTKDBZ

	mov	eax, DWORD PTR _DTK+36
	mov	DWORD PTR _D1T2DAZ+4, eax
	mov	eax, DWORD PTR _DTK+32
	mov	DWORD PTR _D1T2DAZ, eax

; 265  : 	D1T2DBZ=DTKDAZ

	mov	eax, DWORD PTR _DTK+28
	mov	DWORD PTR _D1T2DBZ+4, eax
	mov	eax, DWORD PTR _DTK+24
	mov	DWORD PTR _D1T2DBZ, eax

; 266  : C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
; 267  : 	IF (ISYMM) THEN
; 268  : 	DOVnDA1=DOV2DA2

	mov	eax, DWORD PTR _DOV2DA2+4
	mov	DWORD PTR _DOVN+4, eax
	mov	eax, DWORD PTR _DOV2DA2
	mov	DWORD PTR _DOVN, eax

; 269  : 	DOVnDA2=DOV2DA1

	mov	eax, DWORD PTR _DOV2DA1+4
	mov	DWORD PTR _DOVN+12, eax
	mov	eax, DWORD PTR _DOV2DA1
	mov	DWORD PTR _DOVN+8, eax

; 270  : 	DOVnDA5=DOV2DA5

	mov	eax, DWORD PTR _DOV2DA5+4
	mov	DWORD PTR _DOVN+20, eax
	mov	eax, DWORD PTR _DOV2DA5
	mov	DWORD PTR _DOVN+16, eax

; 271  : 	DOVnDAZ=DOV2DBZ

	mov	eax, DWORD PTR _DOV2DBZ+4
	mov	DWORD PTR _DOVN+28, eax
	mov	eax, DWORD PTR _DOV2DBZ
	mov	DWORD PTR _DOVN+24, eax

; 272  : 	DOVnDBZ=DOV2DAZ

	mov	eax, DWORD PTR _DOV2DAZ+4
	mov	DWORD PTR _DOVN+36, eax
	mov	eax, DWORD PTR _DOV2DAZ
	mov	DWORD PTR _DOVN+32, eax
	mov	eax, DWORD PTR _M$[esp+1392]
	mov	esi, DWORD PTR [eax]

; 273  :               CALL TK(DFLAG,OV2,

	cmp	edi, esi
	jg	SHORT $L98931
	cmp	edi, 1
	jge	SHORT $L98932
$L98931:
	call	___FHranger@0
$L98932:
	mov	eax, DWORD PTR _B12$[esp+1392]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98933
	cmp	edi, 1
	jge	SHORT $L98934
$L98933:
	call	___FHranger@0
$L98934:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98935
	cmp	edi, 1
	jge	SHORT $L98936
$L98935:
	call	___FHranger@0
$L98936:
	mov	eax, DWORD PTR _AA2$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98937
	cmp	edi, 1
	jge	SHORT $L98938
$L98937:
	call	___FHranger@0
$L98938:
	mov	eax, DWORD PTR _RZ2$[esp+1404]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98939
	cmp	edi, 1
	jge	SHORT $L98940
$L98939:
	call	___FHranger@0
$L98940:
	mov	eax, DWORD PTR _RY2$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98941
	cmp	edi, 1
	jge	SHORT $L98942
$L98941:
	call	___FHranger@0
$L98942:
	mov	eax, DWORD PTR _RX2$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98943
	cmp	edi, 1
	jge	SHORT $L98944
$L98943:
	call	___FHranger@0
$L98944:
	mov	eax, DWORD PTR _AA1$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98945
	cmp	edi, 1
	jge	SHORT $L98946
$L98945:
	call	___FHranger@0
$L98946:
	mov	eax, DWORD PTR _RZ1$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98947
	cmp	edi, 1
	jge	SHORT $L98948
$L98947:
	call	___FHranger@0
$L98948:
	mov	eax, DWORD PTR _RY1$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98949
	cmp	edi, 1
	jge	SHORT $L98950
$L98949:
	call	___FHranger@0
$L98950:
	mov	eax, DWORD PTR _RX1$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98951
	cmp	edi, 1
	jge	SHORT $L98952
$L98951:
	call	___FHranger@0
$L98952:
	mov	eax, DWORD PTR _AA1$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98953
	cmp	edi, 1
	jge	SHORT $L98954
$L98953:
	call	___FHranger@0
$L98954:
	mov	eax, DWORD PTR _RZ1$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98955
	cmp	edi, 1
	jge	SHORT $L98956
$L98955:
	call	___FHranger@0
$L98956:
	mov	eax, DWORD PTR _RY1$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98957
	cmp	edi, 1
	jge	SHORT $L98958
$L98957:
	call	___FHranger@0
$L98958:
	mov	eax, DWORD PTR _RX1$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98959
	cmp	edi, 1
	jge	SHORT $L98960
$L98959:
	call	___FHranger@0
$L98960:
	mov	eax, DWORD PTR _AA2$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98961
	cmp	edi, 1
	jge	SHORT $L98962
$L98961:
	call	___FHranger@0
$L98962:
	mov	eax, DWORD PTR _RZ2$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98963
	cmp	edi, 1
	jge	SHORT $L98964
$L98963:
	call	___FHranger@0
$L98964:
	mov	eax, DWORD PTR _RY2$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98965
	cmp	edi, 1
	jge	SHORT $L98966
$L98965:
	call	___FHranger@0
$L98966:
	mov	eax, DWORD PTR _RX2$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1460]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV2
	push	eax
	call	_TK@84
	mov	eax, DWORD PTR _DTK+12
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 279	; 00000117H

; 274  :      1            RX2(L),RY2(L),RZ2(L),AA2(L),
; 275  :      2            RX1(K),RY1(K),RZ1(K),AA1(K),                      
; 276  :      3            RX1(L),RY1(L),RZ1(L),AA1(L),
; 277  :      4            RX2(K),RY2(K),RZ2(K),AA2(K),
; 278  :      5            B12(L),B12(K),T6) 
; 279  : 	D2T2DA1=DTKDA2

	mov	DWORD PTR $T98573[esp+1404], 280	; 00000118H
	mov	DWORD PTR $T98573[esp+1404], 281	; 00000119H
	mov	DWORD PTR $T98573[esp+1404], 282	; 0000011aH
	mov	DWORD PTR $T98573[esp+1404], 283	; 0000011bH
	mov	DWORD PTR $T98573[esp+1404], 285	; 0000011dH
	mov	DWORD PTR $T98573[esp+1404], 286	; 0000011eH
	mov	DWORD PTR $T98573[esp+1404], 287	; 0000011fH
	mov	DWORD PTR $T98573[esp+1404], 288	; 00000120H
	mov	DWORD PTR $T98573[esp+1404], 289	; 00000121H
	mov	DWORD PTR $T98573[esp+1404], 290	; 00000122H
	mov	DWORD PTR _D2T2DA1+4, eax
	mov	eax, DWORD PTR _DTK+8
	mov	DWORD PTR _D2T2DA1, eax

; 280  : 	D2T2DA2=DTKDA1

	mov	eax, DWORD PTR _DTK+4
	mov	DWORD PTR _D2T2DA2+4, eax
	mov	eax, DWORD PTR _DTK
	mov	DWORD PTR _D2T2DA2, eax

; 281  : 	D2T2DA5=DTKDA5

	mov	eax, DWORD PTR _DTK+20
	mov	DWORD PTR _D2T2DA5+4, eax
	mov	eax, DWORD PTR _DTK+16
	mov	DWORD PTR _D2T2DA5, eax

; 282  : 	D2T2DAZ=DTKDBZ

	mov	eax, DWORD PTR _DTK+36
	mov	DWORD PTR _D2T2DAZ+4, eax
	mov	eax, DWORD PTR _DTK+32
	mov	DWORD PTR _D2T2DAZ, eax

; 283  : 	D2T2DBZ=DTKDAZ

	mov	eax, DWORD PTR _DTK+28
	mov	DWORD PTR _D2T2DBZ+4, eax
	mov	eax, DWORD PTR _DTK+24
	mov	DWORD PTR _D2T2DBZ, eax

; 284  : C     ...Nuclear Symmetrization Integral: P_[H1,H2]
; 285  : 	DOVnDA1=DOV3DA2

	mov	eax, DWORD PTR _DOV3DA2+4
	mov	DWORD PTR _DOVN+4, eax
	mov	eax, DWORD PTR _DOV3DA2
	mov	DWORD PTR _DOVN, eax

; 286  : 	DOVnDA2=DOV3DA1

	mov	eax, DWORD PTR _DOV3DA1+4
	mov	DWORD PTR _DOVN+12, eax
	mov	eax, DWORD PTR _DOV3DA1
	mov	DWORD PTR _DOVN+8, eax

; 287  : 	DOVnDA5=DOV3DA5

	mov	eax, DWORD PTR _DOV3DA5+4
	mov	DWORD PTR _DOVN+20, eax
	mov	eax, DWORD PTR _DOV3DA5
	mov	DWORD PTR _DOVN+16, eax

; 288  : 	DOVnDAZ=DOV3DBZ

	mov	eax, DWORD PTR _DOV3DBZ+4
	mov	DWORD PTR _DOVN+28, eax
	mov	eax, DWORD PTR _DOV3DBZ
	mov	DWORD PTR _DOVN+24, eax

; 289  : 	DOVnDBZ=DOV3DAZ

	mov	eax, DWORD PTR _DOV3DAZ+4
	mov	DWORD PTR _DOVN+36, eax
	mov	eax, DWORD PTR _DOV3DAZ
	mov	DWORD PTR _DOVN+32, eax
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	esi, DWORD PTR [eax]

; 290  :               CALL TK(DFLAG,OV3,

	cmp	edi, esi
	jg	SHORT $L98967
	cmp	edi, 1
	jge	SHORT $L98968
$L98967:
	call	___FHranger@0
$L98968:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98332[esp+1392]
	cmp	edi, esi
	jg	SHORT $L98969
	cmp	edi, 1
	jge	SHORT $L98970
$L98969:
	call	___FHranger@0
$L98970:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98333[esp+1392]
	cmp	edi, esi
	jg	SHORT $L98971
	cmp	edi, 1
	jge	SHORT $L98972
$L98971:
	call	___FHranger@0
$L98972:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98334[esp+1392]
	cmp	edi, esi
	jg	SHORT $L98973
	cmp	edi, 1
	jge	SHORT $L98974
$L98973:
	call	___FHranger@0
$L98974:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98335[esp+1392]
	cmp	edi, esi
	jg	SHORT $L98975
	cmp	edi, 1
	jge	SHORT $L98976
$L98975:
	call	___FHranger@0
$L98976:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98336[esp+1392]
	cmp	edi, esi
	jg	SHORT $L98977
	cmp	edi, 1
	jge	SHORT $L98978
$L98977:
	call	___FHranger@0
$L98978:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	push	OFFSET FLAT:_T7
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98337[esp+1396]
	cmp	edi, esi
	jg	SHORT $L98979
	cmp	edi, 1
	jge	SHORT $L98980
$L98979:
	call	___FHranger@0
$L98980:
	mov	eax, DWORD PTR _B12$[esp+1392]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98981
	cmp	edi, 1
	jge	SHORT $L98982
$L98981:
	call	___FHranger@0
$L98982:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98983
	cmp	edi, 1
	jge	SHORT $L98984
$L98983:
	call	___FHranger@0
$L98984:
	mov	eax, DWORD PTR _AA1$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98332[esp+1404]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98333[esp+1412]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98334[esp+1416]
	push	eax
	jg	SHORT $L98985
	cmp	edi, 1
	jge	SHORT $L98986
$L98985:
	call	___FHranger@0
$L98986:
	mov	eax, DWORD PTR _AA1$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98987
	cmp	edi, 1
	jge	SHORT $L98988
$L98987:
	call	___FHranger@0
$L98988:
	mov	eax, DWORD PTR _RZ1$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98989
	cmp	edi, 1
	jge	SHORT $L98990
$L98989:
	call	___FHranger@0
$L98990:
	mov	eax, DWORD PTR _RY1$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98991
	cmp	edi, 1
	jge	SHORT $L98992
$L98991:
	call	___FHranger@0
$L98992:
	mov	eax, DWORD PTR _RX1$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L98993
	cmp	edi, 1
	jge	SHORT $L98994
$L98993:
	call	___FHranger@0
$L98994:
	mov	eax, DWORD PTR _AA2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98335[esp+1436]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98336[esp+1444]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98337[esp+1448]
	push	eax
	jg	SHORT $L98995
	cmp	edi, 1
	jge	SHORT $L98996
$L98995:
	call	___FHranger@0
$L98996:
	mov	eax, DWORD PTR _AA2$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98997
	cmp	edi, 1
	jge	SHORT $L98998
$L98997:
	call	___FHranger@0
$L98998:
	mov	eax, DWORD PTR _RZ2$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L98999
	cmp	edi, 1
	jge	SHORT $L99000
$L98999:
	call	___FHranger@0
$L99000:
	mov	eax, DWORD PTR _RY2$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99001
	cmp	edi, 1
	jge	SHORT $L99002
$L99001:
	call	___FHranger@0
$L99002:
	mov	eax, DWORD PTR _RX2$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1460]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV3
	push	eax
	call	_TK@84
	mov	eax, DWORD PTR _DTK+12
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 296	; 00000128H

; 291  :      1            RX2(L),RY2(L),RZ2(L),AA2(L),
; 292  :      2            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
; 293  :      3            RX1(L),RY1(L),RZ1(L),AA1(L),
; 294  :      4            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
; 295  :      5            B12(L),B12(K),T7)
; 296  : 	D3T2DA1=DTKDA2

	mov	DWORD PTR $T98573[esp+1404], 297	; 00000129H
	mov	DWORD PTR $T98573[esp+1404], 298	; 0000012aH
	mov	DWORD PTR $T98573[esp+1404], 299	; 0000012bH
	mov	DWORD PTR $T98573[esp+1404], 300	; 0000012cH
	mov	DWORD PTR $T98573[esp+1404], 303	; 0000012fH
	mov	DWORD PTR $T98573[esp+1404], 304	; 00000130H
	mov	DWORD PTR $T98573[esp+1404], 305	; 00000131H
	mov	DWORD PTR $T98573[esp+1404], 306	; 00000132H
	mov	DWORD PTR $T98573[esp+1404], 307	; 00000133H
	mov	DWORD PTR $T98573[esp+1404], 308	; 00000134H
	mov	DWORD PTR _D3T2DA1+4, eax
	mov	eax, DWORD PTR _DTK+8
	mov	DWORD PTR _D3T2DA1, eax

; 297  : 	D3T2DA2=DTKDA1

	mov	eax, DWORD PTR _DTK+4
	mov	DWORD PTR _D3T2DA2+4, eax
	mov	eax, DWORD PTR _DTK
	mov	DWORD PTR _D3T2DA2, eax

; 298  : 	D3T2DA5=DTKDA5

	mov	eax, DWORD PTR _DTK+20
	mov	DWORD PTR _D3T2DA5+4, eax
	mov	eax, DWORD PTR _DTK+16
	mov	DWORD PTR _D3T2DA5, eax

; 299  : 	D3T2DAZ=DTKDBZ

	mov	eax, DWORD PTR _DTK+36
	mov	DWORD PTR _D3T2DAZ+4, eax
	mov	eax, DWORD PTR _DTK+32
	mov	DWORD PTR _D3T2DAZ, eax

; 300  : 	D3T2DBZ=DTKDAZ

	mov	eax, DWORD PTR _DTK+28
	mov	DWORD PTR _D3T2DBZ+4, eax
	mov	eax, DWORD PTR _DTK+24
	mov	DWORD PTR _D3T2DBZ, eax

; 301  : C     ...Electro/Nuclear Symmetrization Integral:
; 302  : C         P_[elec.1,elec.2]*P_[H1,H2]
; 303  : 	DOVnDA1=DOV4DA2

	mov	eax, DWORD PTR _DOV4DA2+4
	mov	DWORD PTR _DOVN+4, eax
	mov	eax, DWORD PTR _DOV4DA2
	mov	DWORD PTR _DOVN, eax

; 304  : 	DOVnDA2=DOV4DA1

	mov	eax, DWORD PTR _DOV4DA1+4
	mov	DWORD PTR _DOVN+12, eax
	mov	eax, DWORD PTR _DOV4DA1
	mov	DWORD PTR _DOVN+8, eax

; 305  : 	DOVnDA5=DOV4DA5

	mov	eax, DWORD PTR _DOV4DA5+4
	mov	DWORD PTR _DOVN+20, eax
	mov	eax, DWORD PTR _DOV4DA5
	mov	DWORD PTR _DOVN+16, eax

; 306  : 	DOVnDAZ=DOV4DBZ

	mov	eax, DWORD PTR _DOV4DBZ+4
	mov	DWORD PTR _DOVN+28, eax
	mov	eax, DWORD PTR _DOV4DBZ
	mov	DWORD PTR _DOVN+24, eax

; 307  : 	DOVnDBZ=DOV4DAZ

	mov	eax, DWORD PTR _DOV4DAZ+4
	mov	DWORD PTR _DOVN+36, eax
	mov	eax, DWORD PTR _DOV4DAZ
	mov	DWORD PTR _DOVN+32, eax
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	esi, DWORD PTR [eax]

; 308  :               CALL TK(DFLAG,OV4,

	cmp	edi, esi
	jg	SHORT $L99003
	cmp	edi, 1
	jge	SHORT $L99004
$L99003:
	call	___FHranger@0
$L99004:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98338[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99005
	cmp	edi, 1
	jge	SHORT $L99006
$L99005:
	call	___FHranger@0
$L99006:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98339[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99007
	cmp	edi, 1
	jge	SHORT $L99008
$L99007:
	call	___FHranger@0
$L99008:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98340[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99009
	cmp	edi, 1
	jge	SHORT $L99010
$L99009:
	call	___FHranger@0
$L99010:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98341[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99011
	cmp	edi, 1
	jge	SHORT $L99012
$L99011:
	call	___FHranger@0
$L99012:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98342[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99013
	cmp	edi, 1
	jge	SHORT $L99014
$L99013:
	call	___FHranger@0
$L99014:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	push	OFFSET FLAT:_T8
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98343[esp+1396]
	cmp	edi, esi
	jg	SHORT $L99015
	cmp	edi, 1
	jge	SHORT $L99016
$L99015:
	call	___FHranger@0
$L99016:
	mov	eax, DWORD PTR _B12$[esp+1392]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99017
	cmp	edi, 1
	jge	SHORT $L99018
$L99017:
	call	___FHranger@0
$L99018:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99019
	cmp	edi, 1
	jge	SHORT $L99020
$L99019:
	call	___FHranger@0
$L99020:
	mov	eax, DWORD PTR _AA2$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98338[esp+1404]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98339[esp+1412]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98340[esp+1416]
	push	eax
	jg	SHORT $L99021
	cmp	edi, 1
	jge	SHORT $L99022
$L99021:
	call	___FHranger@0
$L99022:
	mov	eax, DWORD PTR _AA1$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99023
	cmp	edi, 1
	jge	SHORT $L99024
$L99023:
	call	___FHranger@0
$L99024:
	mov	eax, DWORD PTR _RZ1$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99025
	cmp	edi, 1
	jge	SHORT $L99026
$L99025:
	call	___FHranger@0
$L99026:
	mov	eax, DWORD PTR _RY1$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99027
	cmp	edi, 1
	jge	SHORT $L99028
$L99027:
	call	___FHranger@0
$L99028:
	mov	eax, DWORD PTR _RX1$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99029
	cmp	edi, 1
	jge	SHORT $L99030
$L99029:
	call	___FHranger@0
$L99030:
	mov	eax, DWORD PTR _AA1$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98341[esp+1436]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98342[esp+1444]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98343[esp+1448]
	push	eax
	jg	SHORT $L99031
	cmp	edi, 1
	jge	SHORT $L99032
$L99031:
	call	___FHranger@0
$L99032:
	mov	eax, DWORD PTR _AA2$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99033
	cmp	edi, 1
	jge	SHORT $L99034
$L99033:
	call	___FHranger@0
$L99034:
	mov	eax, DWORD PTR _RZ2$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99035
	cmp	edi, 1
	jge	SHORT $L99036
$L99035:
	call	___FHranger@0
$L99036:
	mov	eax, DWORD PTR _RY2$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99037
	cmp	edi, 1
	jge	SHORT $L99038
$L99037:
	call	___FHranger@0
$L99038:
	mov	eax, DWORD PTR _RX2$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1460]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV4
	push	eax
	call	_TK@84
	fld	QWORD PTR _T7
	fadd	QWORD PTR _T4
	fld	QWORD PTR _D3T2DA1
	fadd	QWORD PTR _D4T1DA1
	fxch	ST(1)
	fadd	QWORD PTR _T6
	fld	QWORD PTR _D4T1DA2
	fadd	QWORD PTR _D2T2DA2
	fxch	ST(1)
	fadd	QWORD PTR _T3
	fxch	ST(2)
	fadd	QWORD PTR _D2T2DA1
	fxch	ST(1)
	fadd	QWORD PTR _D3T1DA2
	fxch	ST(2)
	fadd	QWORD PTR _T1
	fxch	ST(1)
	fadd	QWORD PTR _D3T1DA1
	fxch	ST(2)
	fadd	QWORD PTR _D1T1DA2
	fxch	ST(1)
	fadd	QWORD PTR _T8
	fxch	ST(2)
	fadd	QWORD PTR _DTK+8
	fxch	ST(1)
	fadd	QWORD PTR _DTK
	fxch	ST(2)
	fadd	QWORD PTR _T5
	fxch	ST(1)
	fadd	QWORD PTR _D1T1DA1
	fxch	ST(2)
	fadd	QWORD PTR _D2T1DA2
	fxch	ST(1)
	fadd	QWORD PTR _T2
	fxch	ST(2)
	fadd	QWORD PTR _D2T1DA1
	fxch	ST(1)
	fadd	QWORD PTR _D1T2DA2
	fxch	ST(2)
	fdiv	QWORD PTR $T99039
	mov	eax, DWORD PTR _DTK+12
	mov	DWORD PTR $T98573[esp+1404], 314	; 0000013aH

; 309  :      1            RX2(L),RY2(L),RZ2(L),AA2(L),
; 310  :      2            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
; 311  :      3            RX1(L),RY1(L),RZ1(L),AA1(L),
; 312  :      4            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
; 313  :      5            B12(L),B12(K),T8)
; 314  : 	D4T2DA1=DTKDA2

	mov	DWORD PTR _D4T2DA1+4, eax
	mov	DWORD PTR $T98573[esp+1404], 315	; 0000013bH
	mov	eax, DWORD PTR _DTK+8
	mov	DWORD PTR $T98573[esp+1404], 316	; 0000013cH
	mov	DWORD PTR _D4T2DA1, eax

; 315  : 	D4T2DA2=DTKDA1

	mov	DWORD PTR $T98573[esp+1404], 317	; 0000013dH
	mov	eax, DWORD PTR _DTK+4
	mov	DWORD PTR $T98573[esp+1404], 318	; 0000013eH
	mov	DWORD PTR _D4T2DA2+4, eax
	mov	DWORD PTR $T98573[esp+1404], 322	; 00000142H
	mov	eax, DWORD PTR _DTK
	mov	DWORD PTR $T98573[esp+1404], 340	; 00000154H
	mov	DWORD PTR _D4T2DA2, eax

; 316  : 	D4T2DA5=DTKDA5

	mov	DWORD PTR $T98573[esp+1404], 342	; 00000156H
	mov	eax, DWORD PTR _DTK+20
	mov	DWORD PTR $T98573[esp+1404], 344	; 00000158H
	mov	DWORD PTR _D4T2DA5+4, eax
	mov	DWORD PTR $T98573[esp+1404], 346	; 0000015aH
	mov	eax, DWORD PTR _DTK+16
	mov	DWORD PTR $T98573[esp+1404], 348	; 0000015cH
	mov	DWORD PTR _D4T2DA5, eax

; 317  : 	D4T2DAZ=DTKDBZ

	mov	DWORD PTR $T98573[esp+1404], 388	; 00000184H
	mov	eax, DWORD PTR _DTK+36
	mov	DWORD PTR _D4T2DAZ+4, eax
	mov	eax, DWORD PTR _DTK+32
	mov	DWORD PTR _D4T2DAZ, eax

; 318  : 	D4T2DBZ=DTKDAZ

	mov	eax, DWORD PTR _DTK+28
	mov	DWORD PTR _D4T2DBZ+4, eax
	mov	eax, DWORD PTR _DTK+24
	mov	DWORD PTR _D4T2DBZ, eax

; 319  : 	END IF
; 320  : C
; 321  : 	IF (ISYMM) THEN
; 322  :             T=(T1+T2+T3+T4+T5+T6+T7+T8)/2.0D+00

	fxch	ST(1)
	fadd	QWORD PTR _D1T2DA1
	fxch	ST(2)
	fadd	QWORD PTR _D3T2DA2
	fxch	ST(1)
	fstp	QWORD PTR _T

; 323  : 	ELSE
; 324  : 		  T=(T1+T5)/2.0D+00
; 325  : 	END IF
; 326  : 	IF(IPERM.EQ.1) THEN
; 327  : 		T=(T1+T5)/2.0D+00
; 328  : 	END IF
; 329  : 	IF(IPERM.EQ.2) THEN
; 330  : 		T=(T2+T6)/2.0D+00
; 331  : 	END IF
; 332  : 	IF(IPERM.EQ.3) THEN
; 333  : 		T=(T3+T7)/2.0D+00
; 334  : 	END IF
; 335  : 	IF(IPERM.EQ.4) THEN
; 336  : 		T=(T4+T8)/2.0D+00
; 337  : 	END IF
; 338  : C
; 339  : 	IF (ISYMM) THEN
; 340  : 	DTDA1=(D1T1DA1+D2T1DA1+D3T1DA1+D4T1DA1+D1T2DA1+D2T2DA1+D3T2DA1+D4T

	fdiv	QWORD PTR $T99039
	fstp	QWORD PTR _DTDA2
	fdiv	QWORD PTR $T99039
	fstp	QWORD PTR _DTDA1

; 341  : 	-2DA1)/2.00D+00
; 342  : 	DTDA2=(D1T1DA2+D2T1DA2+D3T1DA2+D4T1DA2+D1T2DA2+D2T2DA2+D3T2DA2+D4T
; 343  : 	-2DA2)/2.00D+00
; 344  : 	DTDA5=(D1T1DA5+D2T1DA5+D3T1DA5+D4T1DA5+D1T2DA5+D2T2DA5+D3T2DA5+D4T

	fld	QWORD PTR _D4T1DA5
	fadd	QWORD PTR _D2T2DA5
	fld	QWORD PTR _D2T2DAZ
	fadd	QWORD PTR _D1T1DAZ
	fxch	ST(1)
	fadd	QWORD PTR _D1T1DA5
	fld	QWORD PTR _D1T1DBZ
	fadd	QWORD PTR _D2T2DBZ
	fxch	ST(1)
	fadd	QWORD PTR _D3T1DA5
	fxch	ST(2)
	fadd	QWORD PTR _D3T1DAZ
	fxch	ST(1)
	fadd	QWORD PTR _D3T1DBZ
	fxch	ST(2)
	fadd	QWORD PTR _DTK+16
	fxch	ST(1)
	fadd	QWORD PTR _DTK+32
	fxch	ST(2)
	fadd	QWORD PTR _DTK+24
	fxch	ST(1)
	fadd	QWORD PTR _D2T1DA5
	fxch	ST(2)
	fadd	QWORD PTR _D2T1DAZ
	fxch	ST(1)
	fadd	QWORD PTR _D2T1DBZ
	fxch	ST(2)
	fadd	QWORD PTR _D1T2DA5
	fxch	ST(1)
	fadd	QWORD PTR _D1T2DAZ
	fxch	ST(2)
	fadd	QWORD PTR _D1T2DBZ
	fxch	ST(1)
	fadd	QWORD PTR _D3T2DA5
	fxch	ST(2)
	fadd	QWORD PTR _D3T2DAZ
	fxch	ST(1)
	fadd	QWORD PTR _D3T2DBZ
	fxch	ST(2)
	fdiv	QWORD PTR $T99039
	fxch	ST(1)
	fadd	QWORD PTR _D4T1DAZ
	fxch	ST(2)
	fadd	QWORD PTR _D4T1DBZ
	fxch	ST(1)
	fstp	QWORD PTR _DTDA5

; 345  : 	-2DA5)/2.00D+00
; 346  : 	DTDAZ=(D1T1DAZ+D2T1DAZ+D3T1DAZ+D4T1DAZ+D1T2DAZ+D2T2DAZ+D3T2DAZ+D4T

	fdiv	QWORD PTR $T99039
	fstp	QWORD PTR _DTDBZ
	fdiv	QWORD PTR $T99039
	fstp	QWORD PTR _DTDAZ

; 347  : 	-2DAZ)/2.00D+00
; 348  : 	DTDBZ=(D1T1DBZ+D2T1DBZ+D3T1DBZ+D4T1DBZ+D1T2DBZ+D2T2DBZ+D3T2DBZ+D4T
; 349  : 	-2DBZ)/2.00D+00
; 350  : 	ELSE
; 351  : 	DTDA1=(D1T1DA1+D1T2DA1)/2.00D+00
; 352  : 	DTDA2=(D1T1DA2+D1T2DA2)/2.00D+00
; 353  : 	DTDA5=(D1T1DA5+D1T2DA5)/2.00D+00
; 354  : 	DTDAZ=(D1T1DAZ+D1T2DAZ)/2.00D+00
; 355  : 	DTDBZ=(D1T1DBZ+D1T2DBZ)/2.00D+00
; 356  : 	END IF
; 357  : 	IF(IPERM.EQ.1) THEN
; 358  : 	DTDA1=(D1T1DA1+D1T2DA1)/2.00D+00
; 359  : 	DTDA2=(D1T1DA2+D1T2DA2)/2.00D+00
; 360  : 	DTDA5=(D1T1DA5+D1T2DA5)/2.00D+00
; 361  : 	DTDAZ=(D1T1DAZ+D1T2DAZ)/2.00D+00
; 362  : 	DTDBZ=(D1T1DBZ+D1T2DBZ)/2.00D+00
; 363  : 	END IF
; 364  : 	IF(IPERM.EQ.2) THEN
; 365  : 	DTDA1=(D2T1DA1+D2T2DA1)/2.00D+00
; 366  : 	DTDA2=(D2T1DA2+D2T2DA2)/2.00D+00
; 367  : 	DTDA5=(D2T1DA5+D2T2DA5)/2.00D+00
; 368  : 	DTDAZ=(D2T1DAZ+D2T2DAZ)/2.00D+00
; 369  : 	DTDBZ=(D2T1DBZ+D2T2DBZ)/2.00D+00
; 370  : 	END IF
; 371  : 	IF(IPERM.EQ.3) THEN
; 372  : 	DTDA1=(D3T1DA1+D3T2DA1)/2.00D+00
; 373  : 	DTDA2=(D3T1DA2+D3T2DA2)/2.00D+00
; 374  : 	DTDA5=(D1T1DA5+D3T2DA5)/2.00D+00
; 375  : 	DTDAZ=(D3T1DAZ+D3T2DAZ)/2.00D+00
; 376  : 	DTDBZ=(D3T1DBZ+D3T2DBZ)/2.00D+00
; 377  : 	END IF
; 378  : 	IF(IPERM.EQ.4) THEN
; 379  : 	DTDA1=(D4T1DA1+D4T2DA1)/2.00D+00
; 380  : 	DTDA2=(D4T1DA2+D4T2DA2)/2.00D+00
; 381  : 	DTDA5=(D4T1DA5+D4T2DA5)/2.00D+00
; 382  : 	DTDAZ=(D4T1DAZ+D4T2DAZ)/2.00D+00
; 383  : 	DTDBZ=(D4T1DBZ+D4T2DBZ)/2.00D+00
; 384  : 	END IF
; 385  : C     -----------------------------------------------------
; 386  : C     --- Calculate Nuc./Elec. Attractn. FOR ELECTRON 1 ---
; 387  : C     -----------------------------------------------------
; 388  : 	DOVnDA1=DOV1DA1

	mov	eax, DWORD PTR _DOV1DA1+4
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN1
	mov	edi, DWORD PTR _OVER+4
	push	OFFSET FLAT:_EZ1
	mov	DWORD PTR _DOVN+4, eax
	mov	eax, DWORD PTR _DOV1DA1
	push	OFFSET FLAT:_EY1
	push	OFFSET FLAT:_EX1
	mov	DWORD PTR _DOVN, eax
	mov	DWORD PTR $T98573[esp+1424], 389	; 00000185H

; 389  : 	DOVnDA2=DOV1DA2

	mov	DWORD PTR $T98573[esp+1424], 390	; 00000186H
	mov	DWORD PTR $T98573[esp+1424], 391	; 00000187H
	mov	DWORD PTR $T98573[esp+1424], 392	; 00000188H
	mov	DWORD PTR $T98573[esp+1424], 393	; 00000189H
	mov	eax, DWORD PTR _DOV1DA2+4
	mov	DWORD PTR _DOVN+12, eax
	mov	eax, DWORD PTR _DOV1DA2
	mov	DWORD PTR _DOVN+8, eax

; 390  : 	DOVnDA5=DOV1DA5

	mov	eax, DWORD PTR _DOV1DA5+4
	mov	DWORD PTR _DOVN+20, eax
	mov	eax, DWORD PTR _DOV1DA5
	mov	DWORD PTR _DOVN+16, eax

; 391  : 	DOVnDAZ=DOV1DAZ

	mov	eax, DWORD PTR _DOV1DAZ+4
	mov	DWORD PTR _DOVN+28, eax
	mov	eax, DWORD PTR _DOV1DAZ
	mov	DWORD PTR _DOVN+24, eax

; 392  : 	DOVnDBZ=DOV1DBZ

	mov	eax, DWORD PTR _DOV1DBZ+4
	mov	DWORD PTR _DOVN+36, eax
	mov	eax, DWORD PTR _DOV1DBZ
	mov	DWORD PTR _DOVN+32, eax
	mov	eax, DWORD PTR _M$[esp+1408]
	mov	esi, DWORD PTR [eax]

; 393  :               CALL ENUCL(DFLAG,OV1,RX1(L),RY1(L),RZ1(L),AA1(L),

	cmp	edi, esi
	jg	SHORT $L99040
	cmp	edi, 1
	jge	SHORT $L99041
$L99040:
	call	___FHranger@0
$L99041:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99042
	cmp	edi, 1
	jge	SHORT $L99043
$L99042:
	call	___FHranger@0
$L99043:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99044
	cmp	edi, 1
	jge	SHORT $L99045
$L99044:
	call	___FHranger@0
$L99045:
	mov	eax, DWORD PTR _AA2$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99046
	cmp	edi, 1
	jge	SHORT $L99047
$L99046:
	call	___FHranger@0
$L99047:
	mov	eax, DWORD PTR _RZ2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99048
	cmp	edi, 1
	jge	SHORT $L99049
$L99048:
	call	___FHranger@0
$L99049:
	mov	eax, DWORD PTR _RY2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99050
	cmp	edi, 1
	jge	SHORT $L99051
$L99050:
	call	___FHranger@0
$L99051:
	mov	eax, DWORD PTR _RX2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99052
	cmp	edi, 1
	jge	SHORT $L99053
$L99052:
	call	___FHranger@0
$L99053:
	mov	eax, DWORD PTR _AA2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99054
	cmp	edi, 1
	jge	SHORT $L99055
$L99054:
	call	___FHranger@0
$L99055:
	mov	eax, DWORD PTR _RZ2$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99056
	cmp	edi, 1
	jge	SHORT $L99057
$L99056:
	call	___FHranger@0
$L99057:
	mov	eax, DWORD PTR _RY2$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99058
	cmp	edi, 1
	jge	SHORT $L99059
$L99058:
	call	___FHranger@0
$L99059:
	mov	eax, DWORD PTR _RX2$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99060
	cmp	edi, 1
	jge	SHORT $L99061
$L99060:
	call	___FHranger@0
$L99061:
	mov	eax, DWORD PTR _AA1$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99062
	cmp	edi, 1
	jge	SHORT $L99063
$L99062:
	call	___FHranger@0
$L99063:
	mov	eax, DWORD PTR _RZ1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99064
	cmp	edi, 1
	jge	SHORT $L99065
$L99064:
	call	___FHranger@0
$L99065:
	mov	eax, DWORD PTR _RY1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99066
	cmp	edi, 1
	jge	SHORT $L99067
$L99066:
	call	___FHranger@0
$L99067:
	mov	eax, DWORD PTR _RX1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99068
	cmp	edi, 1
	jge	SHORT $L99069
$L99068:
	call	___FHranger@0
$L99069:
	mov	eax, DWORD PTR _AA1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99070
	cmp	edi, 1
	jge	SHORT $L99071
$L99070:
	call	___FHranger@0
$L99071:
	mov	eax, DWORD PTR _RZ1$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99072
	cmp	edi, 1
	jge	SHORT $L99073
$L99072:
	call	___FHranger@0
$L99073:
	mov	eax, DWORD PTR _RY1$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99074
	cmp	edi, 1
	jge	SHORT $L99075
$L99074:
	call	___FHranger@0
$L99075:
	mov	eax, DWORD PTR _RX1$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV1
	push	eax
	call	_ENUCL@100
	mov	eax, DWORD PTR _DEN+4
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN2
	mov	edi, DWORD PTR _OVER+4
	push	OFFSET FLAT:_EZ1
	mov	DWORD PTR _D1EN1DA1+4, eax
	mov	eax, DWORD PTR _DEN
	push	OFFSET FLAT:_EY1
	push	OFFSET FLAT:_EX1
	mov	DWORD PTR _D1EN1DA1, eax
	mov	DWORD PTR $T98573[esp+1424], 399	; 0000018fH

; 394  :      1           RX1(K),RY1(K),RZ1(K),AA1(K),                      
; 395  :      2           RX2(L),RY2(L),RZ2(L),AA2(L),
; 396  :      3           RX2(K),RY2(K),RZ2(K),AA2(K),
; 397  :      4           B12(L),B12(K),
; 398  :      5           EX1,EY1,EZ1,EN1,KLM)
; 399  : 	D1EN1DA1=DENDA1

	mov	DWORD PTR $T98573[esp+1424], 400	; 00000190H

; 400  : 	D1EN1DA2=DENDA2

	mov	DWORD PTR $T98573[esp+1424], 401	; 00000191H
	mov	DWORD PTR $T98573[esp+1424], 402	; 00000192H
	mov	DWORD PTR $T98573[esp+1424], 403	; 00000193H
	mov	DWORD PTR $T98573[esp+1424], 408	; 00000198H
	mov	DWORD PTR $T98573[esp+1424], 409	; 00000199H
	mov	DWORD PTR $T98573[esp+1424], 410	; 0000019aH
	mov	DWORD PTR $T98573[esp+1424], 411	; 0000019bH
	mov	DWORD PTR $T98573[esp+1424], 412	; 0000019cH
	mov	DWORD PTR $T98573[esp+1424], 413	; 0000019dH
	mov	eax, DWORD PTR _DEN+12
	mov	DWORD PTR _D1EN1DA2+4, eax
	mov	eax, DWORD PTR _DEN+8
	mov	DWORD PTR _D1EN1DA2, eax

; 401  : 	D1EN1DA5=DENDA5

	mov	eax, DWORD PTR _DEN+20
	mov	DWORD PTR _D1EN1DA5+4, eax
	mov	eax, DWORD PTR _DEN+16
	mov	DWORD PTR _D1EN1DA5, eax

; 402  : 	D1EN1DAZ=DENDAZ

	mov	eax, DWORD PTR _DEN+28
	mov	DWORD PTR _D1EN1DAZ+4, eax
	mov	eax, DWORD PTR _DEN+24
	mov	DWORD PTR _D1EN1DAZ, eax

; 403  : 	D1EN1DBZ=DENDBZ

	mov	eax, DWORD PTR _DEN+36
	mov	DWORD PTR _D1EN1DBZ+4, eax
	mov	eax, DWORD PTR _DEN+32
	mov	DWORD PTR _D1EN1DBZ, eax

; 404  : C
; 405  : C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
; 406  : C
; 407  : 	IF (ISYMM) THEN
; 408  : 	DOVnDA1=DOV2DA1

	mov	eax, DWORD PTR _DOV2DA1+4
	mov	DWORD PTR _DOVN+4, eax
	mov	eax, DWORD PTR _DOV2DA1
	mov	DWORD PTR _DOVN, eax

; 409  : 	DOVnDA2=DOV2DA2

	mov	eax, DWORD PTR _DOV2DA2+4
	mov	DWORD PTR _DOVN+12, eax
	mov	eax, DWORD PTR _DOV2DA2
	mov	DWORD PTR _DOVN+8, eax

; 410  : 	DOVnDA5=DOV2DA5

	mov	eax, DWORD PTR _DOV2DA5+4
	mov	DWORD PTR _DOVN+20, eax
	mov	eax, DWORD PTR _DOV2DA5
	mov	DWORD PTR _DOVN+16, eax

; 411  : 	DOVnDAZ=DOV2DAZ

	mov	eax, DWORD PTR _DOV2DAZ+4
	mov	DWORD PTR _DOVN+28, eax
	mov	eax, DWORD PTR _DOV2DAZ
	mov	DWORD PTR _DOVN+24, eax

; 412  : 	DOVnDBZ=DOV2DBZ

	mov	eax, DWORD PTR _DOV2DBZ+4
	mov	DWORD PTR _DOVN+36, eax
	mov	eax, DWORD PTR _DOV2DBZ
	mov	DWORD PTR _DOVN+32, eax
	mov	eax, DWORD PTR _M$[esp+1408]
	mov	esi, DWORD PTR [eax]

; 413  :               CALL ENUCL(DFLAG,OV2,RX1(L),RY1(L),RZ1(L),AA1(L),

	cmp	edi, esi
	jg	SHORT $L99076
	cmp	edi, 1
	jge	SHORT $L99077
$L99076:
	call	___FHranger@0
$L99077:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99078
	cmp	edi, 1
	jge	SHORT $L99079
$L99078:
	call	___FHranger@0
$L99079:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99080
	cmp	edi, 1
	jge	SHORT $L99081
$L99080:
	call	___FHranger@0
$L99081:
	mov	eax, DWORD PTR _AA1$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99082
	cmp	edi, 1
	jge	SHORT $L99083
$L99082:
	call	___FHranger@0
$L99083:
	mov	eax, DWORD PTR _RZ1$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99084
	cmp	edi, 1
	jge	SHORT $L99085
$L99084:
	call	___FHranger@0
$L99085:
	mov	eax, DWORD PTR _RY1$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99086
	cmp	edi, 1
	jge	SHORT $L99087
$L99086:
	call	___FHranger@0
$L99087:
	mov	eax, DWORD PTR _RX1$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99088
	cmp	edi, 1
	jge	SHORT $L99089
$L99088:
	call	___FHranger@0
$L99089:
	mov	eax, DWORD PTR _AA2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99090
	cmp	edi, 1
	jge	SHORT $L99091
$L99090:
	call	___FHranger@0
$L99091:
	mov	eax, DWORD PTR _RZ2$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99092
	cmp	edi, 1
	jge	SHORT $L99093
$L99092:
	call	___FHranger@0
$L99093:
	mov	eax, DWORD PTR _RY2$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99094
	cmp	edi, 1
	jge	SHORT $L99095
$L99094:
	call	___FHranger@0
$L99095:
	mov	eax, DWORD PTR _RX2$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99096
	cmp	edi, 1
	jge	SHORT $L99097
$L99096:
	call	___FHranger@0
$L99097:
	mov	eax, DWORD PTR _AA2$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99098
	cmp	edi, 1
	jge	SHORT $L99099
$L99098:
	call	___FHranger@0
$L99099:
	mov	eax, DWORD PTR _RZ2$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99100
	cmp	edi, 1
	jge	SHORT $L99101
$L99100:
	call	___FHranger@0
$L99101:
	mov	eax, DWORD PTR _RY2$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99102
	cmp	edi, 1
	jge	SHORT $L99103
$L99102:
	call	___FHranger@0
$L99103:
	mov	eax, DWORD PTR _RX2$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99104
	cmp	edi, 1
	jge	SHORT $L99105
$L99104:
	call	___FHranger@0
$L99105:
	mov	eax, DWORD PTR _AA1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99106
	cmp	edi, 1
	jge	SHORT $L99107
$L99106:
	call	___FHranger@0
$L99107:
	mov	eax, DWORD PTR _RZ1$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99108
	cmp	edi, 1
	jge	SHORT $L99109
$L99108:
	call	___FHranger@0
$L99109:
	mov	eax, DWORD PTR _RY1$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99110
	cmp	edi, 1
	jge	SHORT $L99111
$L99110:
	call	___FHranger@0
$L99111:
	mov	eax, DWORD PTR _RX1$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV2
	push	eax
	call	_ENUCL@100
	mov	eax, DWORD PTR _DEN+4
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 419	; 000001a3H

; 414  :      1           RX2(K),RY2(K),RZ2(K),AA2(K),                      
; 415  :      2           RX2(L),RY2(L),RZ2(L),AA2(L),
; 416  :      3           RX1(K),RY1(K),RZ1(K),AA1(K),
; 417  :      4           B12(L),B12(K),
; 418  :      5           EX1,EY1,EZ1,EN2,KLM)                       
; 419  : 	D1EN2DA1=DENDA1

	mov	DWORD PTR $T98573[esp+1404], 420	; 000001a4H
	mov	DWORD PTR $T98573[esp+1404], 421	; 000001a5H
	mov	DWORD PTR $T98573[esp+1404], 422	; 000001a6H
	mov	DWORD PTR $T98573[esp+1404], 423	; 000001a7H
	mov	DWORD PTR $T98573[esp+1404], 427	; 000001abH
	mov	DWORD PTR $T98573[esp+1404], 428	; 000001acH
	mov	DWORD PTR $T98573[esp+1404], 429	; 000001adH
	mov	DWORD PTR $T98573[esp+1404], 430	; 000001aeH
	mov	DWORD PTR $T98573[esp+1404], 431	; 000001afH
	mov	DWORD PTR $T98573[esp+1404], 432	; 000001b0H
	mov	DWORD PTR _D1EN2DA1+4, eax
	mov	eax, DWORD PTR _DEN
	mov	DWORD PTR _D1EN2DA1, eax

; 420  : 	D1EN2DA2=DENDA2

	mov	eax, DWORD PTR _DEN+12
	mov	DWORD PTR _D1EN2DA2+4, eax
	mov	eax, DWORD PTR _DEN+8
	mov	DWORD PTR _D1EN2DA2, eax

; 421  : 	D1EN2DA5=DENDA5

	mov	eax, DWORD PTR _DEN+20
	mov	DWORD PTR _D1EN2DA5+4, eax
	mov	eax, DWORD PTR _DEN+16
	mov	DWORD PTR _D1EN2DA5, eax

; 422  : 	D1EN2DAZ=DENDAZ

	mov	eax, DWORD PTR _DEN+28
	mov	DWORD PTR _D1EN2DAZ+4, eax
	mov	eax, DWORD PTR _DEN+24
	mov	DWORD PTR _D1EN2DAZ, eax

; 423  : 	D1EN2DBZ=DENDBZ

	mov	eax, DWORD PTR _DEN+36
	mov	DWORD PTR _D1EN2DBZ+4, eax
	mov	eax, DWORD PTR _DEN+32
	mov	DWORD PTR _D1EN2DBZ, eax

; 424  : C
; 425  : C     ...Nuclear Symmetrization Integral: P_[H1,H2]
; 426  : C
; 427  : 	DOVnDA1=DOV3DA1

	mov	eax, DWORD PTR _DOV3DA1+4
	mov	DWORD PTR _DOVN+4, eax
	mov	eax, DWORD PTR _DOV3DA1
	mov	DWORD PTR _DOVN, eax

; 428  : 	DOVnDA2=DOV3DA2

	mov	eax, DWORD PTR _DOV3DA2+4
	mov	DWORD PTR _DOVN+12, eax
	mov	eax, DWORD PTR _DOV3DA2
	mov	DWORD PTR _DOVN+8, eax

; 429  : 	DOVnDA5=DOV3DA5

	mov	eax, DWORD PTR _DOV3DA5+4
	mov	DWORD PTR _DOVN+20, eax
	mov	eax, DWORD PTR _DOV3DA5
	mov	DWORD PTR _DOVN+16, eax

; 430  : 	DOVnDAZ=DOV3DAZ

	mov	eax, DWORD PTR _DOV3DAZ+4
	mov	DWORD PTR _DOVN+28, eax
	mov	eax, DWORD PTR _DOV3DAZ
	mov	DWORD PTR _DOVN+24, eax

; 431  : 	DOVnDBZ=DOV3DBZ

	mov	eax, DWORD PTR _DOV3DBZ+4
	mov	DWORD PTR _DOVN+36, eax
	mov	eax, DWORD PTR _DOV3DBZ
	mov	DWORD PTR _DOVN+32, eax
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	esi, DWORD PTR [eax]

; 432  :               CALL ENUCL(DFLAG,OV3,RX1(L),RY1(L),RZ1(L),AA1(L),

	cmp	edi, esi
	jg	SHORT $L99112
	cmp	edi, 1
	jge	SHORT $L99113
$L99112:
	call	___FHranger@0
$L99113:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98344[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99114
	cmp	edi, 1
	jge	SHORT $L99115
$L99114:
	call	___FHranger@0
$L99115:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98345[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99116
	cmp	edi, 1
	jge	SHORT $L99117
$L99116:
	call	___FHranger@0
$L99117:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98346[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99118
	cmp	edi, 1
	jge	SHORT $L99119
$L99118:
	call	___FHranger@0
$L99119:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98347[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99120
	cmp	edi, 1
	jge	SHORT $L99121
$L99120:
	call	___FHranger@0
$L99121:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98348[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99122
	cmp	edi, 1
	jge	SHORT $L99123
$L99122:
	call	___FHranger@0
$L99123:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN3
	push	OFFSET FLAT:_EZ1
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	push	OFFSET FLAT:_EY1
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98349[esp+1408]
	push	OFFSET FLAT:_EX1
	cmp	edi, esi
	jg	SHORT $L99124
	cmp	edi, 1
	jge	SHORT $L99125
$L99124:
	call	___FHranger@0
$L99125:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99126
	cmp	edi, 1
	jge	SHORT $L99127
$L99126:
	call	___FHranger@0
$L99127:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99128
	cmp	edi, 1
	jge	SHORT $L99129
$L99128:
	call	___FHranger@0
$L99129:
	mov	eax, DWORD PTR _AA2$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98344[esp+1420]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98345[esp+1428]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98346[esp+1432]
	push	eax
	jg	SHORT $L99130
	cmp	edi, 1
	jge	SHORT $L99131
$L99130:
	call	___FHranger@0
$L99131:
	mov	eax, DWORD PTR _AA2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99132
	cmp	edi, 1
	jge	SHORT $L99133
$L99132:
	call	___FHranger@0
$L99133:
	mov	eax, DWORD PTR _RZ2$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99134
	cmp	edi, 1
	jge	SHORT $L99135
$L99134:
	call	___FHranger@0
$L99135:
	mov	eax, DWORD PTR _RY2$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99136
	cmp	edi, 1
	jge	SHORT $L99137
$L99136:
	call	___FHranger@0
$L99137:
	mov	eax, DWORD PTR _RX2$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99138
	cmp	edi, 1
	jge	SHORT $L99139
$L99138:
	call	___FHranger@0
$L99139:
	mov	eax, DWORD PTR _AA1$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98347[esp+1452]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98348[esp+1460]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98349[esp+1464]
	push	eax
	jg	SHORT $L99140
	cmp	edi, 1
	jge	SHORT $L99141
$L99140:
	call	___FHranger@0
$L99141:
	mov	eax, DWORD PTR _AA1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99142
	cmp	edi, 1
	jge	SHORT $L99143
$L99142:
	call	___FHranger@0
$L99143:
	mov	eax, DWORD PTR _RZ1$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99144
	cmp	edi, 1
	jge	SHORT $L99145
$L99144:
	call	___FHranger@0
$L99145:
	mov	eax, DWORD PTR _RY1$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99146
	cmp	edi, 1
	jge	SHORT $L99147
$L99146:
	call	___FHranger@0
$L99147:
	mov	eax, DWORD PTR _RX1$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV3
	push	eax
	call	_ENUCL@100
	mov	eax, DWORD PTR _DEN+4
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 438	; 000001b6H

; 433  :      1            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
; 434  :      2            RX2(L),RY2(L),RZ2(L),AA2(L),
; 435  :      3            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
; 436  :      4            B12(L),B12(K),
; 437  :      5            EX1,EY1,EZ1,EN3,KLM)
; 438  : 	D1EN3DA1=DENDA1

	mov	DWORD PTR $T98573[esp+1404], 439	; 000001b7H
	mov	DWORD PTR $T98573[esp+1404], 440	; 000001b8H
	mov	DWORD PTR $T98573[esp+1404], 441	; 000001b9H
	mov	DWORD PTR $T98573[esp+1404], 442	; 000001baH
	mov	DWORD PTR $T98573[esp+1404], 447	; 000001bfH
	mov	DWORD PTR $T98573[esp+1404], 448	; 000001c0H
	mov	DWORD PTR $T98573[esp+1404], 449	; 000001c1H
	mov	DWORD PTR $T98573[esp+1404], 450	; 000001c2H
	mov	DWORD PTR $T98573[esp+1404], 451	; 000001c3H
	mov	DWORD PTR $T98573[esp+1404], 452	; 000001c4H
	mov	DWORD PTR _D1EN3DA1+4, eax
	mov	eax, DWORD PTR _DEN
	mov	DWORD PTR _D1EN3DA1, eax

; 439  : 	D1EN3DA2=DENDA2

	mov	eax, DWORD PTR _DEN+12
	mov	DWORD PTR _D1EN3DA2+4, eax
	mov	eax, DWORD PTR _DEN+8
	mov	DWORD PTR _D1EN3DA2, eax

; 440  : 	D1EN3DA5=DENDA5

	mov	eax, DWORD PTR _DEN+20
	mov	DWORD PTR _D1EN3DA5+4, eax
	mov	eax, DWORD PTR _DEN+16
	mov	DWORD PTR _D1EN3DA5, eax

; 441  : 	D1EN3DAZ=DENDAZ

	mov	eax, DWORD PTR _DEN+28
	mov	DWORD PTR _D1EN3DAZ+4, eax
	mov	eax, DWORD PTR _DEN+24
	mov	DWORD PTR _D1EN3DAZ, eax

; 442  : 	D1EN3DBZ=DENDBZ

	mov	eax, DWORD PTR _DEN+36
	mov	DWORD PTR _D1EN3DBZ+4, eax
	mov	eax, DWORD PTR _DEN+32
	mov	DWORD PTR _D1EN3DBZ, eax

; 443  : C
; 444  : C     ...Electro/Nuclear Symmetrization Integral:
; 445  : C         P_[elec.1,elec.2]*P_[H1,H2]
; 446  : C
; 447  : 	DOVnDA1=DOV4DA1

	mov	eax, DWORD PTR _DOV4DA1+4
	mov	DWORD PTR _DOVN+4, eax
	mov	eax, DWORD PTR _DOV4DA1
	mov	DWORD PTR _DOVN, eax

; 448  : 	DOVnDA2=DOV4DA2

	mov	eax, DWORD PTR _DOV4DA2+4
	mov	DWORD PTR _DOVN+12, eax
	mov	eax, DWORD PTR _DOV4DA2
	mov	DWORD PTR _DOVN+8, eax

; 449  : 	DOVnDA5=DOV4DA5

	mov	eax, DWORD PTR _DOV4DA5+4
	mov	DWORD PTR _DOVN+20, eax
	mov	eax, DWORD PTR _DOV4DA5
	mov	DWORD PTR _DOVN+16, eax

; 450  : 	DOVnDAZ=DOV4DAZ

	mov	eax, DWORD PTR _DOV4DAZ+4
	mov	DWORD PTR _DOVN+28, eax
	mov	eax, DWORD PTR _DOV4DAZ
	mov	DWORD PTR _DOVN+24, eax

; 451  : 	DOVnDBZ=DOV4DBZ

	mov	eax, DWORD PTR _DOV4DBZ+4
	mov	DWORD PTR _DOVN+36, eax
	mov	eax, DWORD PTR _DOV4DBZ
	mov	DWORD PTR _DOVN+32, eax
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	esi, DWORD PTR [eax]

; 452  :               CALL ENUCL(DFLAG,OV4,RX1(L),RY1(L),RZ1(L),AA1(L),

	cmp	edi, esi
	jg	SHORT $L99148
	cmp	edi, 1
	jge	SHORT $L99149
$L99148:
	call	___FHranger@0
$L99149:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98350[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99150
	cmp	edi, 1
	jge	SHORT $L99151
$L99150:
	call	___FHranger@0
$L99151:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98351[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99152
	cmp	edi, 1
	jge	SHORT $L99153
$L99152:
	call	___FHranger@0
$L99153:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98352[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99154
	cmp	edi, 1
	jge	SHORT $L99155
$L99154:
	call	___FHranger@0
$L99155:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98353[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99156
	cmp	edi, 1
	jge	SHORT $L99157
$L99156:
	call	___FHranger@0
$L99157:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98354[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99158
	cmp	edi, 1
	jge	SHORT $L99159
$L99158:
	call	___FHranger@0
$L99159:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN4
	push	OFFSET FLAT:_EZ1
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	push	OFFSET FLAT:_EY1
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98355[esp+1408]
	push	OFFSET FLAT:_EX1
	cmp	edi, esi
	jg	SHORT $L99160
	cmp	edi, 1
	jge	SHORT $L99161
$L99160:
	call	___FHranger@0
$L99161:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99162
	cmp	edi, 1
	jge	SHORT $L99163
$L99162:
	call	___FHranger@0
$L99163:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99164
	cmp	edi, 1
	jge	SHORT $L99165
$L99164:
	call	___FHranger@0
$L99165:
	mov	eax, DWORD PTR _AA1$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98350[esp+1420]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98351[esp+1428]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98352[esp+1432]
	push	eax
	jg	SHORT $L99166
	cmp	edi, 1
	jge	SHORT $L99167
$L99166:
	call	___FHranger@0
$L99167:
	mov	eax, DWORD PTR _AA2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99168
	cmp	edi, 1
	jge	SHORT $L99169
$L99168:
	call	___FHranger@0
$L99169:
	mov	eax, DWORD PTR _RZ2$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99170
	cmp	edi, 1
	jge	SHORT $L99171
$L99170:
	call	___FHranger@0
$L99171:
	mov	eax, DWORD PTR _RY2$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99172
	cmp	edi, 1
	jge	SHORT $L99173
$L99172:
	call	___FHranger@0
$L99173:
	mov	eax, DWORD PTR _RX2$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99174
	cmp	edi, 1
	jge	SHORT $L99175
$L99174:
	call	___FHranger@0
$L99175:
	mov	eax, DWORD PTR _AA2$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98353[esp+1452]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98354[esp+1460]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98355[esp+1464]
	push	eax
	jg	SHORT $L99176
	cmp	edi, 1
	jge	SHORT $L99177
$L99176:
	call	___FHranger@0
$L99177:
	mov	eax, DWORD PTR _AA1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99178
	cmp	edi, 1
	jge	SHORT $L99179
$L99178:
	call	___FHranger@0
$L99179:
	mov	eax, DWORD PTR _RZ1$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99180
	cmp	edi, 1
	jge	SHORT $L99181
$L99180:
	call	___FHranger@0
$L99181:
	mov	eax, DWORD PTR _RY1$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99182
	cmp	edi, 1
	jge	SHORT $L99183
$L99182:
	call	___FHranger@0
$L99183:
	mov	eax, DWORD PTR _RX1$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV4
	push	eax
	call	_ENUCL@100
	mov	eax, DWORD PTR _DEN+4
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN5
	mov	edi, DWORD PTR _OVER+4
	push	OFFSET FLAT:_EZ2
	mov	DWORD PTR _D1EN4DA1+4, eax
	mov	eax, DWORD PTR _DEN
	push	OFFSET FLAT:_EY2
	push	OFFSET FLAT:_EX2
	mov	DWORD PTR _D1EN4DA1, eax
	mov	DWORD PTR $T98573[esp+1424], 458	; 000001caH

; 453  :      1            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
; 454  :      2            RX2(L),RY2(L),RZ2(L),AA2(L),
; 455  :      3            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
; 456  :      4            B12(L),B12(K),
; 457  :      5            EX1,EY1,EZ1,EN4,KLM)
; 458  : 	D1EN4DA1=DENDA1

	mov	DWORD PTR $T98573[esp+1424], 459	; 000001cbH

; 459  : 	D1EN4DA2=DENDA2

	mov	DWORD PTR $T98573[esp+1424], 460	; 000001ccH
	mov	DWORD PTR $T98573[esp+1424], 461	; 000001cdH
	mov	DWORD PTR $T98573[esp+1424], 462	; 000001ceH
	mov	DWORD PTR $T98573[esp+1424], 466	; 000001d2H
	mov	DWORD PTR $T98573[esp+1424], 467	; 000001d3H
	mov	DWORD PTR $T98573[esp+1424], 468	; 000001d4H
	mov	DWORD PTR $T98573[esp+1424], 469	; 000001d5H
	mov	DWORD PTR $T98573[esp+1424], 470	; 000001d6H
	mov	DWORD PTR $T98573[esp+1424], 471	; 000001d7H
	mov	eax, DWORD PTR _DEN+12
	mov	DWORD PTR _D1EN4DA2+4, eax
	mov	eax, DWORD PTR _DEN+8
	mov	DWORD PTR _D1EN4DA2, eax

; 460  : 	D1EN4DA5=DENDA5

	mov	eax, DWORD PTR _DEN+20
	mov	DWORD PTR _D1EN4DA5+4, eax
	mov	eax, DWORD PTR _DEN+16
	mov	DWORD PTR _D1EN4DA5, eax

; 461  : 	D1EN4DAZ=DENDAZ

	mov	eax, DWORD PTR _DEN+28
	mov	DWORD PTR _D1EN4DAZ+4, eax
	mov	eax, DWORD PTR _DEN+24
	mov	DWORD PTR _D1EN4DAZ, eax

; 462  : 	D1EN4DBZ=DENDBZ

	mov	eax, DWORD PTR _DEN+36
	mov	DWORD PTR _D1EN4DBZ+4, eax
	mov	eax, DWORD PTR _DEN+32
	mov	DWORD PTR _D1EN4DBZ, eax

; 463  : 	END IF
; 464  : C  
; 465  : 
; 466  : 	DOVnDA1=DOV1DA1

	mov	eax, DWORD PTR _DOV1DA1+4
	mov	DWORD PTR _DOVN+4, eax
	mov	eax, DWORD PTR _DOV1DA1
	mov	DWORD PTR _DOVN, eax

; 467  : 	DOVnDA2=DOV1DA2

	mov	eax, DWORD PTR _DOV1DA2+4
	mov	DWORD PTR _DOVN+12, eax
	mov	eax, DWORD PTR _DOV1DA2
	mov	DWORD PTR _DOVN+8, eax

; 468  : 	DOVnDA5=DOV1DA5

	mov	eax, DWORD PTR _DOV1DA5+4
	mov	DWORD PTR _DOVN+20, eax
	mov	eax, DWORD PTR _DOV1DA5
	mov	DWORD PTR _DOVN+16, eax

; 469  : 	DOVnDAZ=DOV1DAZ

	mov	eax, DWORD PTR _DOV1DAZ+4
	mov	DWORD PTR _DOVN+28, eax
	mov	eax, DWORD PTR _DOV1DAZ
	mov	DWORD PTR _DOVN+24, eax

; 470  : 	DOVnDBZ=DOV1DBZ

	mov	eax, DWORD PTR _DOV1DBZ+4
	mov	DWORD PTR _DOVN+36, eax
	mov	eax, DWORD PTR _DOV1DBZ
	mov	DWORD PTR _DOVN+32, eax
	mov	eax, DWORD PTR _M$[esp+1408]
	mov	esi, DWORD PTR [eax]

; 471  :               CALL ENUCL(DFLAG,OV1,RX1(L),RY1(L),RZ1(L),AA1(L),

	cmp	edi, esi
	jg	SHORT $L99184
	cmp	edi, 1
	jge	SHORT $L99185
$L99184:
	call	___FHranger@0
$L99185:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99186
	cmp	edi, 1
	jge	SHORT $L99187
$L99186:
	call	___FHranger@0
$L99187:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99188
	cmp	edi, 1
	jge	SHORT $L99189
$L99188:
	call	___FHranger@0
$L99189:
	mov	eax, DWORD PTR _AA2$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99190
	cmp	edi, 1
	jge	SHORT $L99191
$L99190:
	call	___FHranger@0
$L99191:
	mov	eax, DWORD PTR _RZ2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99192
	cmp	edi, 1
	jge	SHORT $L99193
$L99192:
	call	___FHranger@0
$L99193:
	mov	eax, DWORD PTR _RY2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99194
	cmp	edi, 1
	jge	SHORT $L99195
$L99194:
	call	___FHranger@0
$L99195:
	mov	eax, DWORD PTR _RX2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99196
	cmp	edi, 1
	jge	SHORT $L99197
$L99196:
	call	___FHranger@0
$L99197:
	mov	eax, DWORD PTR _AA2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99198
	cmp	edi, 1
	jge	SHORT $L99199
$L99198:
	call	___FHranger@0
$L99199:
	mov	eax, DWORD PTR _RZ2$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99200
	cmp	edi, 1
	jge	SHORT $L99201
$L99200:
	call	___FHranger@0
$L99201:
	mov	eax, DWORD PTR _RY2$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99202
	cmp	edi, 1
	jge	SHORT $L99203
$L99202:
	call	___FHranger@0
$L99203:
	mov	eax, DWORD PTR _RX2$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99204
	cmp	edi, 1
	jge	SHORT $L99205
$L99204:
	call	___FHranger@0
$L99205:
	mov	eax, DWORD PTR _AA1$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99206
	cmp	edi, 1
	jge	SHORT $L99207
$L99206:
	call	___FHranger@0
$L99207:
	mov	eax, DWORD PTR _RZ1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99208
	cmp	edi, 1
	jge	SHORT $L99209
$L99208:
	call	___FHranger@0
$L99209:
	mov	eax, DWORD PTR _RY1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99210
	cmp	edi, 1
	jge	SHORT $L99211
$L99210:
	call	___FHranger@0
$L99211:
	mov	eax, DWORD PTR _RX1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99212
	cmp	edi, 1
	jge	SHORT $L99213
$L99212:
	call	___FHranger@0
$L99213:
	mov	eax, DWORD PTR _AA1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99214
	cmp	edi, 1
	jge	SHORT $L99215
$L99214:
	call	___FHranger@0
$L99215:
	mov	eax, DWORD PTR _RZ1$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99216
	cmp	edi, 1
	jge	SHORT $L99217
$L99216:
	call	___FHranger@0
$L99217:
	mov	eax, DWORD PTR _RY1$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99218
	cmp	edi, 1
	jge	SHORT $L99219
$L99218:
	call	___FHranger@0
$L99219:
	mov	eax, DWORD PTR _RX1$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV1
	push	eax
	call	_ENUCL@100
	mov	eax, DWORD PTR _DEN+44
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN6
	mov	edi, DWORD PTR _OVER+4
	push	OFFSET FLAT:_EZ2
	mov	DWORD PTR _D1EN5DA1+4, eax
	mov	eax, DWORD PTR _DEN+40
	push	OFFSET FLAT:_EY2
	push	OFFSET FLAT:_EX2
	mov	DWORD PTR _D1EN5DA1, eax
	mov	DWORD PTR $T98573[esp+1424], 477	; 000001ddH

; 472  :      1           RX1(K),RY1(K),RZ1(K),AA1(K),                      
; 473  :      2           RX2(L),RY2(L),RZ2(L),AA2(L),
; 474  :      3           RX2(K),RY2(K),RZ2(K),AA2(K),
; 475  :      4           B12(L),B12(K),
; 476  :      5           EX2,EY2,EZ2,EN5,KLM)
; 477  : 	D1EN5DA1=DENDA1P

	mov	DWORD PTR $T98573[esp+1424], 478	; 000001deH

; 478  : 	D1EN5DA2=DENDA2P

	mov	DWORD PTR $T98573[esp+1424], 479	; 000001dfH
	mov	DWORD PTR $T98573[esp+1424], 480	; 000001e0H
	mov	DWORD PTR $T98573[esp+1424], 481	; 000001e1H
	mov	DWORD PTR $T98573[esp+1424], 486	; 000001e6H
	mov	DWORD PTR $T98573[esp+1424], 487	; 000001e7H
	mov	DWORD PTR $T98573[esp+1424], 488	; 000001e8H
	mov	DWORD PTR $T98573[esp+1424], 489	; 000001e9H
	mov	DWORD PTR $T98573[esp+1424], 490	; 000001eaH
	mov	DWORD PTR $T98573[esp+1424], 491	; 000001ebH
	mov	eax, DWORD PTR _DEN+52
	mov	DWORD PTR _D1EN5DA2+4, eax
	mov	eax, DWORD PTR _DEN+48
	mov	DWORD PTR _D1EN5DA2, eax

; 479  : 	D1EN5DA5=DENDA5P

	mov	eax, DWORD PTR _DEN+60
	mov	DWORD PTR _D1EN5DA5+4, eax
	mov	eax, DWORD PTR _DEN+56
	mov	DWORD PTR _D1EN5DA5, eax

; 480  : 	D1EN5DAZ=DENDAZP

	mov	eax, DWORD PTR _DEN+68
	mov	DWORD PTR _D1EN5DAZ+4, eax
	mov	eax, DWORD PTR _DEN+64
	mov	DWORD PTR _D1EN5DAZ, eax

; 481  : 	D1EN5DBZ=DENDBZP

	mov	eax, DWORD PTR _DEN+76
	mov	DWORD PTR _D1EN5DBZ+4, eax
	mov	eax, DWORD PTR _DEN+72
	mov	DWORD PTR _D1EN5DBZ, eax

; 482  : C
; 483  : C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
; 484  : C
; 485  : 	IF (ISYMM) THEN
; 486  : 	DOVnDA1=DOV2DA1

	mov	eax, DWORD PTR _DOV2DA1+4
	mov	DWORD PTR _DOVN+4, eax
	mov	eax, DWORD PTR _DOV2DA1
	mov	DWORD PTR _DOVN, eax

; 487  : 	DOVnDA2=DOV2DA2

	mov	eax, DWORD PTR _DOV2DA2+4
	mov	DWORD PTR _DOVN+12, eax
	mov	eax, DWORD PTR _DOV2DA2
	mov	DWORD PTR _DOVN+8, eax

; 488  : 	DOVnDA5=DOV2DA5

	mov	eax, DWORD PTR _DOV2DA5+4
	mov	DWORD PTR _DOVN+20, eax
	mov	eax, DWORD PTR _DOV2DA5
	mov	DWORD PTR _DOVN+16, eax

; 489  : 	DOVnDAZ=DOV2DAZ

	mov	eax, DWORD PTR _DOV2DAZ+4
	mov	DWORD PTR _DOVN+28, eax
	mov	eax, DWORD PTR _DOV2DAZ
	mov	DWORD PTR _DOVN+24, eax

; 490  : 	DOVnDBZ=DOV2DBZ

	mov	eax, DWORD PTR _DOV2DBZ+4
	mov	DWORD PTR _DOVN+36, eax
	mov	eax, DWORD PTR _DOV2DBZ
	mov	DWORD PTR _DOVN+32, eax
	mov	eax, DWORD PTR _M$[esp+1408]
	mov	esi, DWORD PTR [eax]

; 491  :               CALL ENUCL(DFLAG,OV2,RX1(L),RY1(L),RZ1(L),AA1(L),

	cmp	edi, esi
	jg	SHORT $L99220
	cmp	edi, 1
	jge	SHORT $L99221
$L99220:
	call	___FHranger@0
$L99221:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99222
	cmp	edi, 1
	jge	SHORT $L99223
$L99222:
	call	___FHranger@0
$L99223:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99224
	cmp	edi, 1
	jge	SHORT $L99225
$L99224:
	call	___FHranger@0
$L99225:
	mov	eax, DWORD PTR _AA1$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99226
	cmp	edi, 1
	jge	SHORT $L99227
$L99226:
	call	___FHranger@0
$L99227:
	mov	eax, DWORD PTR _RZ1$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99228
	cmp	edi, 1
	jge	SHORT $L99229
$L99228:
	call	___FHranger@0
$L99229:
	mov	eax, DWORD PTR _RY1$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99230
	cmp	edi, 1
	jge	SHORT $L99231
$L99230:
	call	___FHranger@0
$L99231:
	mov	eax, DWORD PTR _RX1$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99232
	cmp	edi, 1
	jge	SHORT $L99233
$L99232:
	call	___FHranger@0
$L99233:
	mov	eax, DWORD PTR _AA2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99234
	cmp	edi, 1
	jge	SHORT $L99235
$L99234:
	call	___FHranger@0
$L99235:
	mov	eax, DWORD PTR _RZ2$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99236
	cmp	edi, 1
	jge	SHORT $L99237
$L99236:
	call	___FHranger@0
$L99237:
	mov	eax, DWORD PTR _RY2$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99238
	cmp	edi, 1
	jge	SHORT $L99239
$L99238:
	call	___FHranger@0
$L99239:
	mov	eax, DWORD PTR _RX2$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99240
	cmp	edi, 1
	jge	SHORT $L99241
$L99240:
	call	___FHranger@0
$L99241:
	mov	eax, DWORD PTR _AA2$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99242
	cmp	edi, 1
	jge	SHORT $L99243
$L99242:
	call	___FHranger@0
$L99243:
	mov	eax, DWORD PTR _RZ2$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99244
	cmp	edi, 1
	jge	SHORT $L99245
$L99244:
	call	___FHranger@0
$L99245:
	mov	eax, DWORD PTR _RY2$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99246
	cmp	edi, 1
	jge	SHORT $L99247
$L99246:
	call	___FHranger@0
$L99247:
	mov	eax, DWORD PTR _RX2$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99248
	cmp	edi, 1
	jge	SHORT $L99249
$L99248:
	call	___FHranger@0
$L99249:
	mov	eax, DWORD PTR _AA1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99250
	cmp	edi, 1
	jge	SHORT $L99251
$L99250:
	call	___FHranger@0
$L99251:
	mov	eax, DWORD PTR _RZ1$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99252
	cmp	edi, 1
	jge	SHORT $L99253
$L99252:
	call	___FHranger@0
$L99253:
	mov	eax, DWORD PTR _RY1$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99254
	cmp	edi, 1
	jge	SHORT $L99255
$L99254:
	call	___FHranger@0
$L99255:
	mov	eax, DWORD PTR _RX1$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV2
	push	eax
	call	_ENUCL@100
	mov	eax, DWORD PTR _DEN+44
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 497	; 000001f1H

; 492  :      1           RX2(K),RY2(K),RZ2(K),AA2(K),                      
; 493  :      2           RX2(L),RY2(L),RZ2(L),AA2(L),
; 494  :      3           RX1(K),RY1(K),RZ1(K),AA1(K),
; 495  :      4           B12(L),B12(K),
; 496  :      5           EX2,EY2,EZ2,EN6,KLM)                       
; 497  : 	D1EN6DA1=DENDA1P

	mov	DWORD PTR $T98573[esp+1404], 498	; 000001f2H
	mov	DWORD PTR $T98573[esp+1404], 499	; 000001f3H
	mov	DWORD PTR $T98573[esp+1404], 500	; 000001f4H
	mov	DWORD PTR $T98573[esp+1404], 501	; 000001f5H
	mov	DWORD PTR $T98573[esp+1404], 505	; 000001f9H
	mov	DWORD PTR $T98573[esp+1404], 506	; 000001faH
	mov	DWORD PTR $T98573[esp+1404], 507	; 000001fbH
	mov	DWORD PTR $T98573[esp+1404], 508	; 000001fcH
	mov	DWORD PTR $T98573[esp+1404], 509	; 000001fdH
	mov	DWORD PTR $T98573[esp+1404], 510	; 000001feH
	mov	DWORD PTR _D1EN6DA1+4, eax
	mov	eax, DWORD PTR _DEN+40
	mov	DWORD PTR _D1EN6DA1, eax

; 498  : 	D1EN6DA2=DENDA2P

	mov	eax, DWORD PTR _DEN+52
	mov	DWORD PTR _D1EN6DA2+4, eax
	mov	eax, DWORD PTR _DEN+48
	mov	DWORD PTR _D1EN6DA2, eax

; 499  : 	D1EN6DA5=DENDA5P

	mov	eax, DWORD PTR _DEN+60
	mov	DWORD PTR _D1EN6DA5+4, eax
	mov	eax, DWORD PTR _DEN+56
	mov	DWORD PTR _D1EN6DA5, eax

; 500  : 	D1EN6DAZ=DENDAZP

	mov	eax, DWORD PTR _DEN+68
	mov	DWORD PTR _D1EN6DAZ+4, eax
	mov	eax, DWORD PTR _DEN+64
	mov	DWORD PTR _D1EN6DAZ, eax

; 501  : 	D1EN6DBZ=DENDBZP

	mov	eax, DWORD PTR _DEN+76
	mov	DWORD PTR _D1EN6DBZ+4, eax
	mov	eax, DWORD PTR _DEN+72
	mov	DWORD PTR _D1EN6DBZ, eax

; 502  : C
; 503  : C     ...Nuclear Symmetrization Integral: P_[H1,H2]
; 504  : C
; 505  : 	DOVnDA1=DOV3DA1

	mov	eax, DWORD PTR _DOV3DA1+4
	mov	DWORD PTR _DOVN+4, eax
	mov	eax, DWORD PTR _DOV3DA1
	mov	DWORD PTR _DOVN, eax

; 506  : 	DOVnDA2=DOV3DA2

	mov	eax, DWORD PTR _DOV3DA2+4
	mov	DWORD PTR _DOVN+12, eax
	mov	eax, DWORD PTR _DOV3DA2
	mov	DWORD PTR _DOVN+8, eax

; 507  : 	DOVnDA5=DOV3DA5

	mov	eax, DWORD PTR _DOV3DA5+4
	mov	DWORD PTR _DOVN+20, eax
	mov	eax, DWORD PTR _DOV3DA5
	mov	DWORD PTR _DOVN+16, eax

; 508  : 	DOVnDAZ=DOV3DAZ

	mov	eax, DWORD PTR _DOV3DAZ+4
	mov	DWORD PTR _DOVN+28, eax
	mov	eax, DWORD PTR _DOV3DAZ
	mov	DWORD PTR _DOVN+24, eax

; 509  : 	DOVnDBZ=DOV3DBZ

	mov	eax, DWORD PTR _DOV3DBZ+4
	mov	DWORD PTR _DOVN+36, eax
	mov	eax, DWORD PTR _DOV3DBZ
	mov	DWORD PTR _DOVN+32, eax
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	esi, DWORD PTR [eax]

; 510  :               CALL ENUCL(DFLAG,OV3,RX1(L),RY1(L),RZ1(L),AA1(L),

	cmp	edi, esi
	jg	SHORT $L99256
	cmp	edi, 1
	jge	SHORT $L99257
$L99256:
	call	___FHranger@0
$L99257:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98356[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99258
	cmp	edi, 1
	jge	SHORT $L99259
$L99258:
	call	___FHranger@0
$L99259:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98357[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99260
	cmp	edi, 1
	jge	SHORT $L99261
$L99260:
	call	___FHranger@0
$L99261:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98358[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99262
	cmp	edi, 1
	jge	SHORT $L99263
$L99262:
	call	___FHranger@0
$L99263:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98359[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99264
	cmp	edi, 1
	jge	SHORT $L99265
$L99264:
	call	___FHranger@0
$L99265:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98360[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99266
	cmp	edi, 1
	jge	SHORT $L99267
$L99266:
	call	___FHranger@0
$L99267:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN7
	push	OFFSET FLAT:_EZ2
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	push	OFFSET FLAT:_EY2
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98361[esp+1408]
	push	OFFSET FLAT:_EX2
	cmp	edi, esi
	jg	SHORT $L99268
	cmp	edi, 1
	jge	SHORT $L99269
$L99268:
	call	___FHranger@0
$L99269:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99270
	cmp	edi, 1
	jge	SHORT $L99271
$L99270:
	call	___FHranger@0
$L99271:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99272
	cmp	edi, 1
	jge	SHORT $L99273
$L99272:
	call	___FHranger@0
$L99273:
	mov	eax, DWORD PTR _AA2$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98356[esp+1420]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98357[esp+1428]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98358[esp+1432]
	push	eax
	jg	SHORT $L99274
	cmp	edi, 1
	jge	SHORT $L99275
$L99274:
	call	___FHranger@0
$L99275:
	mov	eax, DWORD PTR _AA2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99276
	cmp	edi, 1
	jge	SHORT $L99277
$L99276:
	call	___FHranger@0
$L99277:
	mov	eax, DWORD PTR _RZ2$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99278
	cmp	edi, 1
	jge	SHORT $L99279
$L99278:
	call	___FHranger@0
$L99279:
	mov	eax, DWORD PTR _RY2$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99280
	cmp	edi, 1
	jge	SHORT $L99281
$L99280:
	call	___FHranger@0
$L99281:
	mov	eax, DWORD PTR _RX2$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99282
	cmp	edi, 1
	jge	SHORT $L99283
$L99282:
	call	___FHranger@0
$L99283:
	mov	eax, DWORD PTR _AA1$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98359[esp+1452]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98360[esp+1460]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98361[esp+1464]
	push	eax
	jg	SHORT $L99284
	cmp	edi, 1
	jge	SHORT $L99285
$L99284:
	call	___FHranger@0
$L99285:
	mov	eax, DWORD PTR _AA1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99286
	cmp	edi, 1
	jge	SHORT $L99287
$L99286:
	call	___FHranger@0
$L99287:
	mov	eax, DWORD PTR _RZ1$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99288
	cmp	edi, 1
	jge	SHORT $L99289
$L99288:
	call	___FHranger@0
$L99289:
	mov	eax, DWORD PTR _RY1$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99290
	cmp	edi, 1
	jge	SHORT $L99291
$L99290:
	call	___FHranger@0
$L99291:
	mov	eax, DWORD PTR _RX1$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV3
	push	eax
	call	_ENUCL@100
	mov	eax, DWORD PTR _DEN+44
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 516	; 00000204H

; 511  :      1            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
; 512  :      2            RX2(L),RY2(L),RZ2(L),AA2(L),
; 513  :      3            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
; 514  :      4            B12(L),B12(K),
; 515  :      5            EX2,EY2,EZ2,EN7,KLM)
; 516  : 	D1EN7DA1=DENDA1P

	mov	DWORD PTR $T98573[esp+1404], 517	; 00000205H
	mov	DWORD PTR $T98573[esp+1404], 518	; 00000206H
	mov	DWORD PTR $T98573[esp+1404], 519	; 00000207H
	mov	DWORD PTR $T98573[esp+1404], 520	; 00000208H
	mov	DWORD PTR $T98573[esp+1404], 525	; 0000020dH
	mov	DWORD PTR $T98573[esp+1404], 526	; 0000020eH
	mov	DWORD PTR $T98573[esp+1404], 527	; 0000020fH
	mov	DWORD PTR $T98573[esp+1404], 528	; 00000210H
	mov	DWORD PTR $T98573[esp+1404], 529	; 00000211H
	mov	DWORD PTR $T98573[esp+1404], 530	; 00000212H
	mov	DWORD PTR _D1EN7DA1+4, eax
	mov	eax, DWORD PTR _DEN+40
	mov	DWORD PTR _D1EN7DA1, eax

; 517  : 	D1EN7DA2=DENDA2P

	mov	eax, DWORD PTR _DEN+52
	mov	DWORD PTR _D1EN7DA2+4, eax
	mov	eax, DWORD PTR _DEN+48
	mov	DWORD PTR _D1EN7DA2, eax

; 518  : 	D1EN7DA5=DENDA5P

	mov	eax, DWORD PTR _DEN+60
	mov	DWORD PTR _D1EN7DA5+4, eax
	mov	eax, DWORD PTR _DEN+56
	mov	DWORD PTR _D1EN7DA5, eax

; 519  : 	D1EN7DAZ=DENDAZP

	mov	eax, DWORD PTR _DEN+68
	mov	DWORD PTR _D1EN7DAZ+4, eax
	mov	eax, DWORD PTR _DEN+64
	mov	DWORD PTR _D1EN7DAZ, eax

; 520  : 	D1EN7DBZ=DENDBZP

	mov	eax, DWORD PTR _DEN+76
	mov	DWORD PTR _D1EN7DBZ+4, eax
	mov	eax, DWORD PTR _DEN+72
	mov	DWORD PTR _D1EN7DBZ, eax

; 521  : C
; 522  : C     ...Electro/Nuclear Symmetrization Integral:
; 523  : C         P_[elec.1,elec.2]*P_[H1,H2]
; 524  : C
; 525  : 	DOVnDA1=DOV4DA1

	mov	eax, DWORD PTR _DOV4DA1+4
	mov	DWORD PTR _DOVN+4, eax
	mov	eax, DWORD PTR _DOV4DA1
	mov	DWORD PTR _DOVN, eax

; 526  : 	DOVnDA2=DOV4DA2

	mov	eax, DWORD PTR _DOV4DA2+4
	mov	DWORD PTR _DOVN+12, eax
	mov	eax, DWORD PTR _DOV4DA2
	mov	DWORD PTR _DOVN+8, eax

; 527  : 	DOVnDA5=DOV4DA5

	mov	eax, DWORD PTR _DOV4DA5+4
	mov	DWORD PTR _DOVN+20, eax
	mov	eax, DWORD PTR _DOV4DA5
	mov	DWORD PTR _DOVN+16, eax

; 528  : 	DOVnDAZ=DOV4DAZ

	mov	eax, DWORD PTR _DOV4DAZ+4
	mov	DWORD PTR _DOVN+28, eax
	mov	eax, DWORD PTR _DOV4DAZ
	mov	DWORD PTR _DOVN+24, eax

; 529  : 	DOVnDBZ=DOV4DBZ

	mov	eax, DWORD PTR _DOV4DBZ+4
	mov	DWORD PTR _DOVN+36, eax
	mov	eax, DWORD PTR _DOV4DBZ
	mov	DWORD PTR _DOVN+32, eax
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	esi, DWORD PTR [eax]

; 530  :               CALL ENUCL(DFLAG,OV4,RX1(L),RY1(L),RZ1(L),AA1(L),

	cmp	edi, esi
	jg	SHORT $L99292
	cmp	edi, 1
	jge	SHORT $L99293
$L99292:
	call	___FHranger@0
$L99293:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98362[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99294
	cmp	edi, 1
	jge	SHORT $L99295
$L99294:
	call	___FHranger@0
$L99295:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98363[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99296
	cmp	edi, 1
	jge	SHORT $L99297
$L99296:
	call	___FHranger@0
$L99297:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98364[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99298
	cmp	edi, 1
	jge	SHORT $L99299
$L99298:
	call	___FHranger@0
$L99299:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98365[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99300
	cmp	edi, 1
	jge	SHORT $L99301
$L99300:
	call	___FHranger@0
$L99301:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98366[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99302
	cmp	edi, 1
	jge	SHORT $L99303
$L99302:
	call	___FHranger@0
$L99303:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN8
	push	OFFSET FLAT:_EZ2
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	push	OFFSET FLAT:_EY2
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98367[esp+1408]
	push	OFFSET FLAT:_EX2
	cmp	edi, esi
	jg	SHORT $L99304
	cmp	edi, 1
	jge	SHORT $L99305
$L99304:
	call	___FHranger@0
$L99305:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99306
	cmp	edi, 1
	jge	SHORT $L99307
$L99306:
	call	___FHranger@0
$L99307:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99308
	cmp	edi, 1
	jge	SHORT $L99309
$L99308:
	call	___FHranger@0
$L99309:
	mov	eax, DWORD PTR _AA1$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98362[esp+1420]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98363[esp+1428]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98364[esp+1432]
	push	eax
	jg	SHORT $L99310
	cmp	edi, 1
	jge	SHORT $L99311
$L99310:
	call	___FHranger@0
$L99311:
	mov	eax, DWORD PTR _AA2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99312
	cmp	edi, 1
	jge	SHORT $L99313
$L99312:
	call	___FHranger@0
$L99313:
	mov	eax, DWORD PTR _RZ2$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99314
	cmp	edi, 1
	jge	SHORT $L99315
$L99314:
	call	___FHranger@0
$L99315:
	mov	eax, DWORD PTR _RY2$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99316
	cmp	edi, 1
	jge	SHORT $L99317
$L99316:
	call	___FHranger@0
$L99317:
	mov	eax, DWORD PTR _RX2$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99318
	cmp	edi, 1
	jge	SHORT $L99319
$L99318:
	call	___FHranger@0
$L99319:
	mov	eax, DWORD PTR _AA2$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98365[esp+1452]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98366[esp+1460]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98367[esp+1464]
	push	eax
	jg	SHORT $L99320
	cmp	edi, 1
	jge	SHORT $L99321
$L99320:
	call	___FHranger@0
$L99321:
	mov	eax, DWORD PTR _AA1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99322
	cmp	edi, 1
	jge	SHORT $L99323
$L99322:
	call	___FHranger@0
$L99323:
	mov	eax, DWORD PTR _RZ1$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99324
	cmp	edi, 1
	jge	SHORT $L99325
$L99324:
	call	___FHranger@0
$L99325:
	mov	eax, DWORD PTR _RY1$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99326
	cmp	edi, 1
	jge	SHORT $L99327
$L99326:
	call	___FHranger@0
$L99327:
	mov	eax, DWORD PTR _RX1$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV4
	push	eax
	call	_ENUCL@100
	mov	eax, DWORD PTR _DEN+44
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN9
	mov	edi, DWORD PTR _OVER+4
	push	OFFSET FLAT:_EZ1
	mov	DWORD PTR _D1EN8DA1+4, eax
	mov	eax, DWORD PTR _DEN+40
	push	OFFSET FLAT:_EY1
	push	OFFSET FLAT:_EX1
	mov	DWORD PTR _D1EN8DA1, eax
	mov	DWORD PTR $T98573[esp+1424], 536	; 00000218H

; 531  :      1            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
; 532  :      2            RX2(L),RY2(L),RZ2(L),AA2(L),
; 533  :      3            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
; 534  :      4            B12(L),B12(K),
; 535  :      5            EX2,EY2,EZ2,EN8,KLM)
; 536  : 	D1EN8DA1=DENDA1P

	mov	DWORD PTR $T98573[esp+1424], 537	; 00000219H

; 537  : 	D1EN8DA2=DENDA2P

	mov	DWORD PTR $T98573[esp+1424], 538	; 0000021aH
	mov	DWORD PTR $T98573[esp+1424], 539	; 0000021bH
	mov	DWORD PTR $T98573[esp+1424], 540	; 0000021cH
	mov	DWORD PTR $T98573[esp+1424], 546	; 00000222H
	mov	DWORD PTR $T98573[esp+1424], 547	; 00000223H
	mov	DWORD PTR $T98573[esp+1424], 548	; 00000224H
	mov	DWORD PTR $T98573[esp+1424], 549	; 00000225H
	mov	DWORD PTR $T98573[esp+1424], 550	; 00000226H
	mov	DWORD PTR $T98573[esp+1424], 551	; 00000227H
	mov	eax, DWORD PTR _DEN+52
	mov	DWORD PTR _D1EN8DA2+4, eax
	mov	eax, DWORD PTR _DEN+48
	mov	DWORD PTR _D1EN8DA2, eax

; 538  : 	D1EN8DA5=DENDA5P

	mov	eax, DWORD PTR _DEN+60
	mov	DWORD PTR _D1EN8DA5+4, eax
	mov	eax, DWORD PTR _DEN+56
	mov	DWORD PTR _D1EN8DA5, eax

; 539  : 	D1EN8DAZ=DENDAZP

	mov	eax, DWORD PTR _DEN+68
	mov	DWORD PTR _D1EN8DAZ+4, eax
	mov	eax, DWORD PTR _DEN+64
	mov	DWORD PTR _D1EN8DAZ, eax

; 540  : 	D1EN8DBZ=DENDBZP

	mov	eax, DWORD PTR _DEN+76
	mov	DWORD PTR _D1EN8DBZ+4, eax
	mov	eax, DWORD PTR _DEN+72
	mov	DWORD PTR _D1EN8DBZ, eax

; 541  : 	END IF
; 542  : C 
; 543  : C     -----------------------------------------------------
; 544  : C     --- Calculate Nuc./Elec. Attractn. FOR ELECTRON 2 ---
; 545  : C     -----------------------------------------------------
; 546  : 	DOVnDA1=DOV1DA2

	mov	eax, DWORD PTR _DOV1DA2+4
	mov	DWORD PTR _DOVN+4, eax
	mov	eax, DWORD PTR _DOV1DA2
	mov	DWORD PTR _DOVN, eax

; 547  : 	DOVnDA2=DOV1DA1

	mov	eax, DWORD PTR _DOV1DA1+4
	mov	DWORD PTR _DOVN+12, eax
	mov	eax, DWORD PTR _DOV1DA1
	mov	DWORD PTR _DOVN+8, eax

; 548  : 	DOVnDA5=DOV1DA5

	mov	eax, DWORD PTR _DOV1DA5+4
	mov	DWORD PTR _DOVN+20, eax
	mov	eax, DWORD PTR _DOV1DA5
	mov	DWORD PTR _DOVN+16, eax

; 549  : 	DOVnDAZ=DOV1DBZ

	mov	eax, DWORD PTR _DOV1DBZ+4
	mov	DWORD PTR _DOVN+28, eax
	mov	eax, DWORD PTR _DOV1DBZ
	mov	DWORD PTR _DOVN+24, eax

; 550  : 	DOVnDBZ=DOV1DAZ

	mov	eax, DWORD PTR _DOV1DAZ+4
	mov	DWORD PTR _DOVN+36, eax
	mov	eax, DWORD PTR _DOV1DAZ
	mov	DWORD PTR _DOVN+32, eax
	mov	eax, DWORD PTR _M$[esp+1408]
	mov	esi, DWORD PTR [eax]

; 551  :               CALL ENUCL(DFLAG,OV1,RX2(L),RY2(L),RZ2(L),AA2(L),

	cmp	edi, esi
	jg	SHORT $L99328
	cmp	edi, 1
	jge	SHORT $L99329
$L99328:
	call	___FHranger@0
$L99329:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99330
	cmp	edi, 1
	jge	SHORT $L99331
$L99330:
	call	___FHranger@0
$L99331:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99332
	cmp	edi, 1
	jge	SHORT $L99333
$L99332:
	call	___FHranger@0
$L99333:
	mov	eax, DWORD PTR _AA1$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99334
	cmp	edi, 1
	jge	SHORT $L99335
$L99334:
	call	___FHranger@0
$L99335:
	mov	eax, DWORD PTR _RZ1$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99336
	cmp	edi, 1
	jge	SHORT $L99337
$L99336:
	call	___FHranger@0
$L99337:
	mov	eax, DWORD PTR _RY1$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99338
	cmp	edi, 1
	jge	SHORT $L99339
$L99338:
	call	___FHranger@0
$L99339:
	mov	eax, DWORD PTR _RX1$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99340
	cmp	edi, 1
	jge	SHORT $L99341
$L99340:
	call	___FHranger@0
$L99341:
	mov	eax, DWORD PTR _AA1$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99342
	cmp	edi, 1
	jge	SHORT $L99343
$L99342:
	call	___FHranger@0
$L99343:
	mov	eax, DWORD PTR _RZ1$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99344
	cmp	edi, 1
	jge	SHORT $L99345
$L99344:
	call	___FHranger@0
$L99345:
	mov	eax, DWORD PTR _RY1$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99346
	cmp	edi, 1
	jge	SHORT $L99347
$L99346:
	call	___FHranger@0
$L99347:
	mov	eax, DWORD PTR _RX1$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99348
	cmp	edi, 1
	jge	SHORT $L99349
$L99348:
	call	___FHranger@0
$L99349:
	mov	eax, DWORD PTR _AA2$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99350
	cmp	edi, 1
	jge	SHORT $L99351
$L99350:
	call	___FHranger@0
$L99351:
	mov	eax, DWORD PTR _RZ2$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99352
	cmp	edi, 1
	jge	SHORT $L99353
$L99352:
	call	___FHranger@0
$L99353:
	mov	eax, DWORD PTR _RY2$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99354
	cmp	edi, 1
	jge	SHORT $L99355
$L99354:
	call	___FHranger@0
$L99355:
	mov	eax, DWORD PTR _RX2$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99356
	cmp	edi, 1
	jge	SHORT $L99357
$L99356:
	call	___FHranger@0
$L99357:
	mov	eax, DWORD PTR _AA2$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99358
	cmp	edi, 1
	jge	SHORT $L99359
$L99358:
	call	___FHranger@0
$L99359:
	mov	eax, DWORD PTR _RZ2$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99360
	cmp	edi, 1
	jge	SHORT $L99361
$L99360:
	call	___FHranger@0
$L99361:
	mov	eax, DWORD PTR _RY2$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99362
	cmp	edi, 1
	jge	SHORT $L99363
$L99362:
	call	___FHranger@0
$L99363:
	mov	eax, DWORD PTR _RX2$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV1
	push	eax
	call	_ENUCL@100
	mov	eax, DWORD PTR _DEN+12
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN10
	mov	edi, DWORD PTR _OVER+4
	push	OFFSET FLAT:_EZ1
	mov	DWORD PTR _D2EN1DA1+4, eax
	mov	eax, DWORD PTR _DEN+8
	push	OFFSET FLAT:_EY1
	push	OFFSET FLAT:_EX1
	mov	DWORD PTR _D2EN1DA1, eax
	mov	DWORD PTR $T98573[esp+1424], 557	; 0000022dH

; 552  :      1           RX2(K),RY2(K),RZ2(K),AA2(K),                      
; 553  :      2           RX1(L),RY1(L),RZ1(L),AA1(L),
; 554  :      3           RX1(K),RY1(K),RZ1(K),AA1(K),
; 555  :      4           B12(L),B12(K),
; 556  :      5           EX1,EY1,EZ1,EN9,KLM)                       
; 557  : 	D2EN1DA1=DENDA2

	mov	DWORD PTR $T98573[esp+1424], 558	; 0000022eH

; 558  : 	D2EN1DA2=DENDA1

	mov	DWORD PTR $T98573[esp+1424], 559	; 0000022fH
	mov	DWORD PTR $T98573[esp+1424], 560	; 00000230H
	mov	DWORD PTR $T98573[esp+1424], 561	; 00000231H
	mov	DWORD PTR $T98573[esp+1424], 566	; 00000236H
	mov	DWORD PTR $T98573[esp+1424], 567	; 00000237H
	mov	DWORD PTR $T98573[esp+1424], 568	; 00000238H
	mov	DWORD PTR $T98573[esp+1424], 569	; 00000239H
	mov	DWORD PTR $T98573[esp+1424], 570	; 0000023aH
	mov	DWORD PTR $T98573[esp+1424], 571	; 0000023bH
	mov	eax, DWORD PTR _DEN+4
	mov	DWORD PTR _D2EN1DA2+4, eax
	mov	eax, DWORD PTR _DEN
	mov	DWORD PTR _D2EN1DA2, eax

; 559  : 	D2EN1DA5=DENDA5

	mov	eax, DWORD PTR _DEN+20
	mov	DWORD PTR _D2EN1DA5+4, eax
	mov	eax, DWORD PTR _DEN+16
	mov	DWORD PTR _D2EN1DA5, eax

; 560  : 	D2EN1DAZ=DENDBZ

	mov	eax, DWORD PTR _DEN+36
	mov	DWORD PTR _D2EN1DAZ+4, eax
	mov	eax, DWORD PTR _DEN+32
	mov	DWORD PTR _D2EN1DAZ, eax

; 561  : 	D2EN1DBZ=DENDAZ

	mov	eax, DWORD PTR _DEN+28
	mov	DWORD PTR _D2EN1DBZ+4, eax
	mov	eax, DWORD PTR _DEN+24
	mov	DWORD PTR _D2EN1DBZ, eax

; 562  : C
; 563  : C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
; 564  : C
; 565  : 	IF (ISYMM) THEN
; 566  : 	DOVnDA1=DOV2DA2

	mov	eax, DWORD PTR _DOV2DA2+4
	mov	DWORD PTR _DOVN+4, eax
	mov	eax, DWORD PTR _DOV2DA2
	mov	DWORD PTR _DOVN, eax

; 567  : 	DOVnDA2=DOV2DA1

	mov	eax, DWORD PTR _DOV2DA1+4
	mov	DWORD PTR _DOVN+12, eax
	mov	eax, DWORD PTR _DOV2DA1
	mov	DWORD PTR _DOVN+8, eax

; 568  : 	DOVnDA5=DOV2DA5

	mov	eax, DWORD PTR _DOV2DA5+4
	mov	DWORD PTR _DOVN+20, eax
	mov	eax, DWORD PTR _DOV2DA5
	mov	DWORD PTR _DOVN+16, eax

; 569  : 	DOVnDAZ=DOV2DBZ

	mov	eax, DWORD PTR _DOV2DBZ+4
	mov	DWORD PTR _DOVN+28, eax
	mov	eax, DWORD PTR _DOV2DBZ
	mov	DWORD PTR _DOVN+24, eax

; 570  : 	DOVnDBZ=DOV2DAZ

	mov	eax, DWORD PTR _DOV2DAZ+4
	mov	DWORD PTR _DOVN+36, eax
	mov	eax, DWORD PTR _DOV2DAZ
	mov	DWORD PTR _DOVN+32, eax
	mov	eax, DWORD PTR _M$[esp+1408]
	mov	esi, DWORD PTR [eax]

; 571  :               CALL ENUCL(DFLAG,OV2,RX2(L),RY2(L),RZ2(L),AA2(L),

	cmp	edi, esi
	jg	SHORT $L99364
	cmp	edi, 1
	jge	SHORT $L99365
$L99364:
	call	___FHranger@0
$L99365:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99366
	cmp	edi, 1
	jge	SHORT $L99367
$L99366:
	call	___FHranger@0
$L99367:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99368
	cmp	edi, 1
	jge	SHORT $L99369
$L99368:
	call	___FHranger@0
$L99369:
	mov	eax, DWORD PTR _AA2$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99370
	cmp	edi, 1
	jge	SHORT $L99371
$L99370:
	call	___FHranger@0
$L99371:
	mov	eax, DWORD PTR _RZ2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99372
	cmp	edi, 1
	jge	SHORT $L99373
$L99372:
	call	___FHranger@0
$L99373:
	mov	eax, DWORD PTR _RY2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99374
	cmp	edi, 1
	jge	SHORT $L99375
$L99374:
	call	___FHranger@0
$L99375:
	mov	eax, DWORD PTR _RX2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99376
	cmp	edi, 1
	jge	SHORT $L99377
$L99376:
	call	___FHranger@0
$L99377:
	mov	eax, DWORD PTR _AA1$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99378
	cmp	edi, 1
	jge	SHORT $L99379
$L99378:
	call	___FHranger@0
$L99379:
	mov	eax, DWORD PTR _RZ1$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99380
	cmp	edi, 1
	jge	SHORT $L99381
$L99380:
	call	___FHranger@0
$L99381:
	mov	eax, DWORD PTR _RY1$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99382
	cmp	edi, 1
	jge	SHORT $L99383
$L99382:
	call	___FHranger@0
$L99383:
	mov	eax, DWORD PTR _RX1$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99384
	cmp	edi, 1
	jge	SHORT $L99385
$L99384:
	call	___FHranger@0
$L99385:
	mov	eax, DWORD PTR _AA1$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99386
	cmp	edi, 1
	jge	SHORT $L99387
$L99386:
	call	___FHranger@0
$L99387:
	mov	eax, DWORD PTR _RZ1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99388
	cmp	edi, 1
	jge	SHORT $L99389
$L99388:
	call	___FHranger@0
$L99389:
	mov	eax, DWORD PTR _RY1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99390
	cmp	edi, 1
	jge	SHORT $L99391
$L99390:
	call	___FHranger@0
$L99391:
	mov	eax, DWORD PTR _RX1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99392
	cmp	edi, 1
	jge	SHORT $L99393
$L99392:
	call	___FHranger@0
$L99393:
	mov	eax, DWORD PTR _AA2$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99394
	cmp	edi, 1
	jge	SHORT $L99395
$L99394:
	call	___FHranger@0
$L99395:
	mov	eax, DWORD PTR _RZ2$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99396
	cmp	edi, 1
	jge	SHORT $L99397
$L99396:
	call	___FHranger@0
$L99397:
	mov	eax, DWORD PTR _RY2$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99398
	cmp	edi, 1
	jge	SHORT $L99399
$L99398:
	call	___FHranger@0
$L99399:
	mov	eax, DWORD PTR _RX2$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV2
	push	eax
	call	_ENUCL@100
	mov	eax, DWORD PTR _DEN+12
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 577	; 00000241H

; 572  :      1           RX1(K),RY1(K),RZ1(K),AA1(K),                      
; 573  :      2           RX1(L),RY1(L),RZ1(L),AA1(L),
; 574  :      3           RX2(K),RY2(K),RZ2(K),AA2(K),
; 575  :      4           B12(L),B12(K),
; 576  :      5           EX1,EY1,EZ1,EN10,KLM)                       
; 577  : 	D2EN2DA1=DENDA2

	mov	DWORD PTR $T98573[esp+1404], 578	; 00000242H
	mov	DWORD PTR $T98573[esp+1404], 579	; 00000243H
	mov	DWORD PTR $T98573[esp+1404], 580	; 00000244H
	mov	DWORD PTR $T98573[esp+1404], 581	; 00000245H
	mov	DWORD PTR $T98573[esp+1404], 585	; 00000249H
	mov	DWORD PTR $T98573[esp+1404], 586	; 0000024aH
	mov	DWORD PTR $T98573[esp+1404], 587	; 0000024bH
	mov	DWORD PTR $T98573[esp+1404], 588	; 0000024cH
	mov	DWORD PTR $T98573[esp+1404], 589	; 0000024dH
	mov	DWORD PTR $T98573[esp+1404], 590	; 0000024eH
	mov	DWORD PTR _D2EN2DA1+4, eax
	mov	eax, DWORD PTR _DEN+8
	mov	DWORD PTR _D2EN2DA1, eax

; 578  : 	D2EN2DA2=DENDA1

	mov	eax, DWORD PTR _DEN+4
	mov	DWORD PTR _D2EN2DA2+4, eax
	mov	eax, DWORD PTR _DEN
	mov	DWORD PTR _D2EN2DA2, eax

; 579  : 	D2EN2DA5=DENDA5

	mov	eax, DWORD PTR _DEN+20
	mov	DWORD PTR _D2EN2DA5+4, eax
	mov	eax, DWORD PTR _DEN+16
	mov	DWORD PTR _D2EN2DA5, eax

; 580  : 	D2EN2DAZ=DENDBZ

	mov	eax, DWORD PTR _DEN+36
	mov	DWORD PTR _D2EN2DAZ+4, eax
	mov	eax, DWORD PTR _DEN+32
	mov	DWORD PTR _D2EN2DAZ, eax

; 581  : 	D2EN2DBZ=DENDAZ

	mov	eax, DWORD PTR _DEN+28
	mov	DWORD PTR _D2EN2DBZ+4, eax
	mov	eax, DWORD PTR _DEN+24
	mov	DWORD PTR _D2EN2DBZ, eax

; 582  : C
; 583  : C     ...Nuclear Symmetrization Integral: P_[H1,H2]
; 584  : C
; 585  : 	DOVnDA1=DOV3DA2

	mov	eax, DWORD PTR _DOV3DA2+4
	mov	DWORD PTR _DOVN+4, eax
	mov	eax, DWORD PTR _DOV3DA2
	mov	DWORD PTR _DOVN, eax

; 586  : 	DOVnDA2=DOV3DA1

	mov	eax, DWORD PTR _DOV3DA1+4
	mov	DWORD PTR _DOVN+12, eax
	mov	eax, DWORD PTR _DOV3DA1
	mov	DWORD PTR _DOVN+8, eax

; 587  : 	DOVnDA5=DOV3DA5

	mov	eax, DWORD PTR _DOV3DA5+4
	mov	DWORD PTR _DOVN+20, eax
	mov	eax, DWORD PTR _DOV3DA5
	mov	DWORD PTR _DOVN+16, eax

; 588  : 	DOVnDAZ=DOV3DBZ

	mov	eax, DWORD PTR _DOV3DBZ+4
	mov	DWORD PTR _DOVN+28, eax
	mov	eax, DWORD PTR _DOV3DBZ
	mov	DWORD PTR _DOVN+24, eax

; 589  : 	DOVnDBZ=DOV3DAZ

	mov	eax, DWORD PTR _DOV3DAZ+4
	mov	DWORD PTR _DOVN+36, eax
	mov	eax, DWORD PTR _DOV3DAZ
	mov	DWORD PTR _DOVN+32, eax
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	esi, DWORD PTR [eax]

; 590  :               CALL ENUCL(DFLAG,OV3,RX2(L),RY2(L),RZ2(L),AA2(L),

	cmp	edi, esi
	jg	SHORT $L99400
	cmp	edi, 1
	jge	SHORT $L99401
$L99400:
	call	___FHranger@0
$L99401:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98368[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99402
	cmp	edi, 1
	jge	SHORT $L99403
$L99402:
	call	___FHranger@0
$L99403:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98369[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99404
	cmp	edi, 1
	jge	SHORT $L99405
$L99404:
	call	___FHranger@0
$L99405:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98370[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99406
	cmp	edi, 1
	jge	SHORT $L99407
$L99406:
	call	___FHranger@0
$L99407:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98371[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99408
	cmp	edi, 1
	jge	SHORT $L99409
$L99408:
	call	___FHranger@0
$L99409:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98372[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99410
	cmp	edi, 1
	jge	SHORT $L99411
$L99410:
	call	___FHranger@0
$L99411:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN11
	push	OFFSET FLAT:_EZ1
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	push	OFFSET FLAT:_EY1
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98373[esp+1408]
	push	OFFSET FLAT:_EX1
	cmp	edi, esi
	jg	SHORT $L99412
	cmp	edi, 1
	jge	SHORT $L99413
$L99412:
	call	___FHranger@0
$L99413:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99414
	cmp	edi, 1
	jge	SHORT $L99415
$L99414:
	call	___FHranger@0
$L99415:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99416
	cmp	edi, 1
	jge	SHORT $L99417
$L99416:
	call	___FHranger@0
$L99417:
	mov	eax, DWORD PTR _AA1$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98368[esp+1420]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98369[esp+1428]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98370[esp+1432]
	push	eax
	jg	SHORT $L99418
	cmp	edi, 1
	jge	SHORT $L99419
$L99418:
	call	___FHranger@0
$L99419:
	mov	eax, DWORD PTR _AA1$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99420
	cmp	edi, 1
	jge	SHORT $L99421
$L99420:
	call	___FHranger@0
$L99421:
	mov	eax, DWORD PTR _RZ1$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99422
	cmp	edi, 1
	jge	SHORT $L99423
$L99422:
	call	___FHranger@0
$L99423:
	mov	eax, DWORD PTR _RY1$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99424
	cmp	edi, 1
	jge	SHORT $L99425
$L99424:
	call	___FHranger@0
$L99425:
	mov	eax, DWORD PTR _RX1$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99426
	cmp	edi, 1
	jge	SHORT $L99427
$L99426:
	call	___FHranger@0
$L99427:
	mov	eax, DWORD PTR _AA2$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98371[esp+1452]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98372[esp+1460]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98373[esp+1464]
	push	eax
	jg	SHORT $L99428
	cmp	edi, 1
	jge	SHORT $L99429
$L99428:
	call	___FHranger@0
$L99429:
	mov	eax, DWORD PTR _AA2$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99430
	cmp	edi, 1
	jge	SHORT $L99431
$L99430:
	call	___FHranger@0
$L99431:
	mov	eax, DWORD PTR _RZ2$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99432
	cmp	edi, 1
	jge	SHORT $L99433
$L99432:
	call	___FHranger@0
$L99433:
	mov	eax, DWORD PTR _RY2$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99434
	cmp	edi, 1
	jge	SHORT $L99435
$L99434:
	call	___FHranger@0
$L99435:
	mov	eax, DWORD PTR _RX2$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV3
	push	eax
	call	_ENUCL@100
	mov	eax, DWORD PTR _DEN+12
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 596	; 00000254H

; 591  :      1            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
; 592  :      2            RX1(L),RY1(L),RZ1(L),AA1(L),
; 593  :      3            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
; 594  :      4            B12(L),B12(K),
; 595  :      5            EX1,EY1,EZ1,EN11,KLM)
; 596  : 	D2EN3DA1=DENDA2

	mov	DWORD PTR $T98573[esp+1404], 597	; 00000255H
	mov	DWORD PTR $T98573[esp+1404], 598	; 00000256H
	mov	DWORD PTR $T98573[esp+1404], 599	; 00000257H
	mov	DWORD PTR $T98573[esp+1404], 600	; 00000258H
	mov	DWORD PTR $T98573[esp+1404], 605	; 0000025dH
	mov	DWORD PTR $T98573[esp+1404], 606	; 0000025eH
	mov	DWORD PTR $T98573[esp+1404], 607	; 0000025fH
	mov	DWORD PTR $T98573[esp+1404], 608	; 00000260H
	mov	DWORD PTR $T98573[esp+1404], 609	; 00000261H
	mov	DWORD PTR $T98573[esp+1404], 610	; 00000262H
	mov	DWORD PTR _D2EN3DA1+4, eax
	mov	eax, DWORD PTR _DEN+8
	mov	DWORD PTR _D2EN3DA1, eax

; 597  : 	D2EN3DA2=DENDA1

	mov	eax, DWORD PTR _DEN+4
	mov	DWORD PTR _D2EN3DA2+4, eax
	mov	eax, DWORD PTR _DEN
	mov	DWORD PTR _D2EN3DA2, eax

; 598  : 	D2EN3DA5=DENDA5

	mov	eax, DWORD PTR _DEN+20
	mov	DWORD PTR _D2EN3DA5+4, eax
	mov	eax, DWORD PTR _DEN+16
	mov	DWORD PTR _D2EN3DA5, eax

; 599  : 	D2EN3DAZ=DENDBZ

	mov	eax, DWORD PTR _DEN+36
	mov	DWORD PTR _D2EN3DAZ+4, eax
	mov	eax, DWORD PTR _DEN+32
	mov	DWORD PTR _D2EN3DAZ, eax

; 600  : 	D2EN3DBZ=DENDAZ

	mov	eax, DWORD PTR _DEN+28
	mov	DWORD PTR _D2EN3DBZ+4, eax
	mov	eax, DWORD PTR _DEN+24
	mov	DWORD PTR _D2EN3DBZ, eax

; 601  : C
; 602  : C     ...Electro/Nuclear Symmetrization Integral:
; 603  : C         P_[elec.1,elec.2]*P_[H1,H2]
; 604  : C
; 605  : 	DOVnDA1=DOV4DA2

	mov	eax, DWORD PTR _DOV4DA2+4
	mov	DWORD PTR _DOVN+4, eax
	mov	eax, DWORD PTR _DOV4DA2
	mov	DWORD PTR _DOVN, eax

; 606  : 	DOVnDA2=DOV4DA1

	mov	eax, DWORD PTR _DOV4DA1+4
	mov	DWORD PTR _DOVN+12, eax
	mov	eax, DWORD PTR _DOV4DA1
	mov	DWORD PTR _DOVN+8, eax

; 607  : 	DOVnDA5=DOV4DA5

	mov	eax, DWORD PTR _DOV4DA5+4
	mov	DWORD PTR _DOVN+20, eax
	mov	eax, DWORD PTR _DOV4DA5
	mov	DWORD PTR _DOVN+16, eax

; 608  : 	DOVnDAZ=DOV4DBZ

	mov	eax, DWORD PTR _DOV4DBZ+4
	mov	DWORD PTR _DOVN+28, eax
	mov	eax, DWORD PTR _DOV4DBZ
	mov	DWORD PTR _DOVN+24, eax

; 609  : 	DOVnDBZ=DOV4DAZ

	mov	eax, DWORD PTR _DOV4DAZ+4
	mov	DWORD PTR _DOVN+36, eax
	mov	eax, DWORD PTR _DOV4DAZ
	mov	DWORD PTR _DOVN+32, eax
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	esi, DWORD PTR [eax]

; 610  :               CALL ENUCL(DFLAG,OV4,RX2(L),RY2(L),RZ2(L),AA2(L),

	cmp	edi, esi
	jg	SHORT $L99436
	cmp	edi, 1
	jge	SHORT $L99437
$L99436:
	call	___FHranger@0
$L99437:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98374[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99438
	cmp	edi, 1
	jge	SHORT $L99439
$L99438:
	call	___FHranger@0
$L99439:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98375[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99440
	cmp	edi, 1
	jge	SHORT $L99441
$L99440:
	call	___FHranger@0
$L99441:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98376[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99442
	cmp	edi, 1
	jge	SHORT $L99443
$L99442:
	call	___FHranger@0
$L99443:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98377[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99444
	cmp	edi, 1
	jge	SHORT $L99445
$L99444:
	call	___FHranger@0
$L99445:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98378[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99446
	cmp	edi, 1
	jge	SHORT $L99447
$L99446:
	call	___FHranger@0
$L99447:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN12
	push	OFFSET FLAT:_EZ1
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	push	OFFSET FLAT:_EY1
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98379[esp+1408]
	push	OFFSET FLAT:_EX1
	cmp	edi, esi
	jg	SHORT $L99448
	cmp	edi, 1
	jge	SHORT $L99449
$L99448:
	call	___FHranger@0
$L99449:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99450
	cmp	edi, 1
	jge	SHORT $L99451
$L99450:
	call	___FHranger@0
$L99451:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99452
	cmp	edi, 1
	jge	SHORT $L99453
$L99452:
	call	___FHranger@0
$L99453:
	mov	eax, DWORD PTR _AA2$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98374[esp+1420]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98375[esp+1428]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98376[esp+1432]
	push	eax
	jg	SHORT $L99454
	cmp	edi, 1
	jge	SHORT $L99455
$L99454:
	call	___FHranger@0
$L99455:
	mov	eax, DWORD PTR _AA1$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99456
	cmp	edi, 1
	jge	SHORT $L99457
$L99456:
	call	___FHranger@0
$L99457:
	mov	eax, DWORD PTR _RZ1$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99458
	cmp	edi, 1
	jge	SHORT $L99459
$L99458:
	call	___FHranger@0
$L99459:
	mov	eax, DWORD PTR _RY1$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99460
	cmp	edi, 1
	jge	SHORT $L99461
$L99460:
	call	___FHranger@0
$L99461:
	mov	eax, DWORD PTR _RX1$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99462
	cmp	edi, 1
	jge	SHORT $L99463
$L99462:
	call	___FHranger@0
$L99463:
	mov	eax, DWORD PTR _AA1$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98377[esp+1452]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98378[esp+1460]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98379[esp+1464]
	push	eax
	jg	SHORT $L99464
	cmp	edi, 1
	jge	SHORT $L99465
$L99464:
	call	___FHranger@0
$L99465:
	mov	eax, DWORD PTR _AA2$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99466
	cmp	edi, 1
	jge	SHORT $L99467
$L99466:
	call	___FHranger@0
$L99467:
	mov	eax, DWORD PTR _RZ2$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99468
	cmp	edi, 1
	jge	SHORT $L99469
$L99468:
	call	___FHranger@0
$L99469:
	mov	eax, DWORD PTR _RY2$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99470
	cmp	edi, 1
	jge	SHORT $L99471
$L99470:
	call	___FHranger@0
$L99471:
	mov	eax, DWORD PTR _RX2$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV4
	push	eax
	call	_ENUCL@100
	mov	eax, DWORD PTR _DEN+12
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN13
	mov	edi, DWORD PTR _OVER+4
	push	OFFSET FLAT:_EZ2
	mov	DWORD PTR _D2EN4DA1+4, eax
	mov	eax, DWORD PTR _DEN+8
	push	OFFSET FLAT:_EY2
	push	OFFSET FLAT:_EX2
	mov	DWORD PTR _D2EN4DA1, eax
	mov	DWORD PTR $T98573[esp+1424], 616	; 00000268H

; 611  :      1            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
; 612  :      2            RX1(L),RY1(L),RZ1(L),AA1(L),
; 613  :      3            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
; 614  :      4            B12(L),B12(K),
; 615  :      5            EX1,EY1,EZ1,EN12,KLM)
; 616  : 	D2EN4DA1=DENDA2

	mov	DWORD PTR $T98573[esp+1424], 617	; 00000269H

; 617  : 	D2EN4DA2=DENDA1

	mov	DWORD PTR $T98573[esp+1424], 618	; 0000026aH
	mov	DWORD PTR $T98573[esp+1424], 619	; 0000026bH
	mov	DWORD PTR $T98573[esp+1424], 620	; 0000026cH
	mov	DWORD PTR $T98573[esp+1424], 625	; 00000271H
	mov	DWORD PTR $T98573[esp+1424], 626	; 00000272H
	mov	DWORD PTR $T98573[esp+1424], 627	; 00000273H
	mov	DWORD PTR $T98573[esp+1424], 628	; 00000274H
	mov	DWORD PTR $T98573[esp+1424], 629	; 00000275H
	mov	DWORD PTR $T98573[esp+1424], 630	; 00000276H
	mov	eax, DWORD PTR _DEN+4
	mov	DWORD PTR _D2EN4DA2+4, eax
	mov	eax, DWORD PTR _DEN
	mov	DWORD PTR _D2EN4DA2, eax

; 618  : 	D2EN4DA5=DENDA5

	mov	eax, DWORD PTR _DEN+20
	mov	DWORD PTR _D2EN4DA5+4, eax
	mov	eax, DWORD PTR _DEN+16
	mov	DWORD PTR _D2EN4DA5, eax

; 619  : 	D2EN4DAZ=DENDBZ

	mov	eax, DWORD PTR _DEN+36
	mov	DWORD PTR _D2EN4DAZ+4, eax
	mov	eax, DWORD PTR _DEN+32
	mov	DWORD PTR _D2EN4DAZ, eax

; 620  : 	D2EN4DBZ=DENDAZ

	mov	eax, DWORD PTR _DEN+28
	mov	DWORD PTR _D2EN4DBZ+4, eax
	mov	eax, DWORD PTR _DEN+24
	mov	DWORD PTR _D2EN4DBZ, eax

; 621  : 	END IF
; 622  : C 
; 623  : 
; 624  :  
; 625  : 	DOVnDA1=DOV1DA2

	mov	eax, DWORD PTR _DOV1DA2+4
	mov	DWORD PTR _DOVN+4, eax
	mov	eax, DWORD PTR _DOV1DA2
	mov	DWORD PTR _DOVN, eax

; 626  : 	DOVnDA2=DOV1DA1

	mov	eax, DWORD PTR _DOV1DA1+4
	mov	DWORD PTR _DOVN+12, eax
	mov	eax, DWORD PTR _DOV1DA1
	mov	DWORD PTR _DOVN+8, eax

; 627  : 	DOVnDA5=DOV1DA5

	mov	eax, DWORD PTR _DOV1DA5+4
	mov	DWORD PTR _DOVN+20, eax
	mov	eax, DWORD PTR _DOV1DA5
	mov	DWORD PTR _DOVN+16, eax

; 628  : 	DOVnDAZ=DOV1DBZ

	mov	eax, DWORD PTR _DOV1DBZ+4
	mov	DWORD PTR _DOVN+28, eax
	mov	eax, DWORD PTR _DOV1DBZ
	mov	DWORD PTR _DOVN+24, eax

; 629  : 	DOVnDBZ=DOV1DAZ

	mov	eax, DWORD PTR _DOV1DAZ+4
	mov	DWORD PTR _DOVN+36, eax
	mov	eax, DWORD PTR _DOV1DAZ
	mov	DWORD PTR _DOVN+32, eax
	mov	eax, DWORD PTR _M$[esp+1408]
	mov	esi, DWORD PTR [eax]

; 630  :               CALL ENUCL(DFLAG,OV1,RX2(L),RY2(L),RZ2(L),AA2(L),

	cmp	edi, esi
	jg	SHORT $L99472
	cmp	edi, 1
	jge	SHORT $L99473
$L99472:
	call	___FHranger@0
$L99473:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99474
	cmp	edi, 1
	jge	SHORT $L99475
$L99474:
	call	___FHranger@0
$L99475:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99476
	cmp	edi, 1
	jge	SHORT $L99477
$L99476:
	call	___FHranger@0
$L99477:
	mov	eax, DWORD PTR _AA1$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99478
	cmp	edi, 1
	jge	SHORT $L99479
$L99478:
	call	___FHranger@0
$L99479:
	mov	eax, DWORD PTR _RZ1$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99480
	cmp	edi, 1
	jge	SHORT $L99481
$L99480:
	call	___FHranger@0
$L99481:
	mov	eax, DWORD PTR _RY1$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99482
	cmp	edi, 1
	jge	SHORT $L99483
$L99482:
	call	___FHranger@0
$L99483:
	mov	eax, DWORD PTR _RX1$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99484
	cmp	edi, 1
	jge	SHORT $L99485
$L99484:
	call	___FHranger@0
$L99485:
	mov	eax, DWORD PTR _AA1$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99486
	cmp	edi, 1
	jge	SHORT $L99487
$L99486:
	call	___FHranger@0
$L99487:
	mov	eax, DWORD PTR _RZ1$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99488
	cmp	edi, 1
	jge	SHORT $L99489
$L99488:
	call	___FHranger@0
$L99489:
	mov	eax, DWORD PTR _RY1$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99490
	cmp	edi, 1
	jge	SHORT $L99491
$L99490:
	call	___FHranger@0
$L99491:
	mov	eax, DWORD PTR _RX1$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99492
	cmp	edi, 1
	jge	SHORT $L99493
$L99492:
	call	___FHranger@0
$L99493:
	mov	eax, DWORD PTR _AA2$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99494
	cmp	edi, 1
	jge	SHORT $L99495
$L99494:
	call	___FHranger@0
$L99495:
	mov	eax, DWORD PTR _RZ2$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99496
	cmp	edi, 1
	jge	SHORT $L99497
$L99496:
	call	___FHranger@0
$L99497:
	mov	eax, DWORD PTR _RY2$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99498
	cmp	edi, 1
	jge	SHORT $L99499
$L99498:
	call	___FHranger@0
$L99499:
	mov	eax, DWORD PTR _RX2$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99500
	cmp	edi, 1
	jge	SHORT $L99501
$L99500:
	call	___FHranger@0
$L99501:
	mov	eax, DWORD PTR _AA2$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99502
	cmp	edi, 1
	jge	SHORT $L99503
$L99502:
	call	___FHranger@0
$L99503:
	mov	eax, DWORD PTR _RZ2$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99504
	cmp	edi, 1
	jge	SHORT $L99505
$L99504:
	call	___FHranger@0
$L99505:
	mov	eax, DWORD PTR _RY2$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99506
	cmp	edi, 1
	jge	SHORT $L99507
$L99506:
	call	___FHranger@0
$L99507:
	mov	eax, DWORD PTR _RX2$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV1
	push	eax
	call	_ENUCL@100
	mov	eax, DWORD PTR _DEN+52
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN14
	mov	edi, DWORD PTR _OVER+4
	push	OFFSET FLAT:_EZ2
	mov	DWORD PTR _D2EN5DA1+4, eax
	mov	eax, DWORD PTR _DEN+48
	push	OFFSET FLAT:_EY2
	push	OFFSET FLAT:_EX2
	mov	DWORD PTR _D2EN5DA1, eax
	mov	DWORD PTR $T98573[esp+1424], 636	; 0000027cH

; 631  :      1           RX2(K),RY2(K),RZ2(K),AA2(K),                      
; 632  :      2           RX1(L),RY1(L),RZ1(L),AA1(L),
; 633  :      3           RX1(K),RY1(K),RZ1(K),AA1(K),
; 634  :      4           B12(L),B12(K),
; 635  :      5           EX2,EY2,EZ2,EN13,KLM)                       
; 636  : 	D2EN5DA1=DENDA2P

	mov	DWORD PTR $T98573[esp+1424], 637	; 0000027dH

; 637  : 	D2EN5DA2=DENDA1P

	mov	DWORD PTR $T98573[esp+1424], 638	; 0000027eH
	mov	DWORD PTR $T98573[esp+1424], 639	; 0000027fH
	mov	DWORD PTR $T98573[esp+1424], 640	; 00000280H
	mov	DWORD PTR $T98573[esp+1424], 645	; 00000285H
	mov	DWORD PTR $T98573[esp+1424], 646	; 00000286H
	mov	DWORD PTR $T98573[esp+1424], 647	; 00000287H
	mov	DWORD PTR $T98573[esp+1424], 648	; 00000288H
	mov	DWORD PTR $T98573[esp+1424], 649	; 00000289H
	mov	DWORD PTR $T98573[esp+1424], 650	; 0000028aH
	mov	eax, DWORD PTR _DEN+44
	mov	DWORD PTR _D2EN5DA2+4, eax
	mov	eax, DWORD PTR _DEN+40
	mov	DWORD PTR _D2EN5DA2, eax

; 638  : 	D2EN5DA5=DENDA5P

	mov	eax, DWORD PTR _DEN+60
	mov	DWORD PTR _D2EN5DA5+4, eax
	mov	eax, DWORD PTR _DEN+56
	mov	DWORD PTR _D2EN5DA5, eax

; 639  : 	D2EN5DAZ=DENDBZP

	mov	eax, DWORD PTR _DEN+76
	mov	DWORD PTR _D2EN5DAZ+4, eax
	mov	eax, DWORD PTR _DEN+72
	mov	DWORD PTR _D2EN5DAZ, eax

; 640  : 	D2EN5DBZ=DENDAZP

	mov	eax, DWORD PTR _DEN+68
	mov	DWORD PTR _D2EN5DBZ+4, eax
	mov	eax, DWORD PTR _DEN+64
	mov	DWORD PTR _D2EN5DBZ, eax

; 641  : C
; 642  : C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
; 643  : C
; 644  : 	IF (ISYMM) THEN
; 645  : 	DOVnDA1=DOV2DA2

	mov	eax, DWORD PTR _DOV2DA2+4
	mov	DWORD PTR _DOVN+4, eax
	mov	eax, DWORD PTR _DOV2DA2
	mov	DWORD PTR _DOVN, eax

; 646  : 	DOVnDA2=DOV2DA1

	mov	eax, DWORD PTR _DOV2DA1+4
	mov	DWORD PTR _DOVN+12, eax
	mov	eax, DWORD PTR _DOV2DA1
	mov	DWORD PTR _DOVN+8, eax

; 647  : 	DOVnDA5=DOV2DA5

	mov	eax, DWORD PTR _DOV2DA5+4
	mov	DWORD PTR _DOVN+20, eax
	mov	eax, DWORD PTR _DOV2DA5
	mov	DWORD PTR _DOVN+16, eax

; 648  : 	DOVnDAZ=DOV2DBZ

	mov	eax, DWORD PTR _DOV2DBZ+4
	mov	DWORD PTR _DOVN+28, eax
	mov	eax, DWORD PTR _DOV2DBZ
	mov	DWORD PTR _DOVN+24, eax

; 649  : 	DOVnDBZ=DOV2DAZ

	mov	eax, DWORD PTR _DOV2DAZ+4
	mov	DWORD PTR _DOVN+36, eax
	mov	eax, DWORD PTR _DOV2DAZ
	mov	DWORD PTR _DOVN+32, eax
	mov	eax, DWORD PTR _M$[esp+1408]
	mov	esi, DWORD PTR [eax]

; 650  :               CALL ENUCL(DFLAG,OV2,RX2(L),RY2(L),RZ2(L),AA2(L),

	cmp	edi, esi
	jg	SHORT $L99508
	cmp	edi, 1
	jge	SHORT $L99509
$L99508:
	call	___FHranger@0
$L99509:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99510
	cmp	edi, 1
	jge	SHORT $L99511
$L99510:
	call	___FHranger@0
$L99511:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99512
	cmp	edi, 1
	jge	SHORT $L99513
$L99512:
	call	___FHranger@0
$L99513:
	mov	eax, DWORD PTR _AA2$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99514
	cmp	edi, 1
	jge	SHORT $L99515
$L99514:
	call	___FHranger@0
$L99515:
	mov	eax, DWORD PTR _RZ2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99516
	cmp	edi, 1
	jge	SHORT $L99517
$L99516:
	call	___FHranger@0
$L99517:
	mov	eax, DWORD PTR _RY2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99518
	cmp	edi, 1
	jge	SHORT $L99519
$L99518:
	call	___FHranger@0
$L99519:
	mov	eax, DWORD PTR _RX2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99520
	cmp	edi, 1
	jge	SHORT $L99521
$L99520:
	call	___FHranger@0
$L99521:
	mov	eax, DWORD PTR _AA1$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99522
	cmp	edi, 1
	jge	SHORT $L99523
$L99522:
	call	___FHranger@0
$L99523:
	mov	eax, DWORD PTR _RZ1$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99524
	cmp	edi, 1
	jge	SHORT $L99525
$L99524:
	call	___FHranger@0
$L99525:
	mov	eax, DWORD PTR _RY1$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99526
	cmp	edi, 1
	jge	SHORT $L99527
$L99526:
	call	___FHranger@0
$L99527:
	mov	eax, DWORD PTR _RX1$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99528
	cmp	edi, 1
	jge	SHORT $L99529
$L99528:
	call	___FHranger@0
$L99529:
	mov	eax, DWORD PTR _AA1$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99530
	cmp	edi, 1
	jge	SHORT $L99531
$L99530:
	call	___FHranger@0
$L99531:
	mov	eax, DWORD PTR _RZ1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99532
	cmp	edi, 1
	jge	SHORT $L99533
$L99532:
	call	___FHranger@0
$L99533:
	mov	eax, DWORD PTR _RY1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99534
	cmp	edi, 1
	jge	SHORT $L99535
$L99534:
	call	___FHranger@0
$L99535:
	mov	eax, DWORD PTR _RX1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99536
	cmp	edi, 1
	jge	SHORT $L99537
$L99536:
	call	___FHranger@0
$L99537:
	mov	eax, DWORD PTR _AA2$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99538
	cmp	edi, 1
	jge	SHORT $L99539
$L99538:
	call	___FHranger@0
$L99539:
	mov	eax, DWORD PTR _RZ2$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99540
	cmp	edi, 1
	jge	SHORT $L99541
$L99540:
	call	___FHranger@0
$L99541:
	mov	eax, DWORD PTR _RY2$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99542
	cmp	edi, 1
	jge	SHORT $L99543
$L99542:
	call	___FHranger@0
$L99543:
	mov	eax, DWORD PTR _RX2$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV2
	push	eax
	call	_ENUCL@100
	mov	eax, DWORD PTR _DEN+52
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 656	; 00000290H

; 651  :      1           RX1(K),RY1(K),RZ1(K),AA1(K),                      
; 652  :      2           RX1(L),RY1(L),RZ1(L),AA1(L),
; 653  :      3           RX2(K),RY2(K),RZ2(K),AA2(K),
; 654  :      4           B12(L),B12(K),
; 655  :      5           EX2,EY2,EZ2,EN14,KLM)  
; 656  : 	D2EN6DA1=DENDA2P

	mov	DWORD PTR $T98573[esp+1404], 657	; 00000291H
	mov	DWORD PTR $T98573[esp+1404], 658	; 00000292H
	mov	DWORD PTR $T98573[esp+1404], 659	; 00000293H
	mov	DWORD PTR $T98573[esp+1404], 660	; 00000294H
	mov	DWORD PTR $T98573[esp+1404], 664	; 00000298H
	mov	DWORD PTR $T98573[esp+1404], 665	; 00000299H
	mov	DWORD PTR $T98573[esp+1404], 666	; 0000029aH
	mov	DWORD PTR $T98573[esp+1404], 667	; 0000029bH
	mov	DWORD PTR $T98573[esp+1404], 668	; 0000029cH
	mov	DWORD PTR $T98573[esp+1404], 669	; 0000029dH
	mov	DWORD PTR _D2EN6DA1+4, eax
	mov	eax, DWORD PTR _DEN+48
	mov	DWORD PTR _D2EN6DA1, eax

; 657  : 	D2EN6DA2=DENDA1P

	mov	eax, DWORD PTR _DEN+44
	mov	DWORD PTR _D2EN6DA2+4, eax
	mov	eax, DWORD PTR _DEN+40
	mov	DWORD PTR _D2EN6DA2, eax

; 658  : 	D2EN6DA5=DENDA5P

	mov	eax, DWORD PTR _DEN+60
	mov	DWORD PTR _D2EN6DA5+4, eax
	mov	eax, DWORD PTR _DEN+56
	mov	DWORD PTR _D2EN6DA5, eax

; 659  : 	D2EN6DAZ=DENDBZP

	mov	eax, DWORD PTR _DEN+76
	mov	DWORD PTR _D2EN6DAZ+4, eax
	mov	eax, DWORD PTR _DEN+72
	mov	DWORD PTR _D2EN6DAZ, eax

; 660  : 	D2EN6DBZ=DENDAZP

	mov	eax, DWORD PTR _DEN+68
	mov	DWORD PTR _D2EN6DBZ+4, eax
	mov	eax, DWORD PTR _DEN+64
	mov	DWORD PTR _D2EN6DBZ, eax

; 661  : C
; 662  : C     ...Nuclear Symmetrization Integral: P_[H1,H2]
; 663  : C
; 664  : 	DOVnDA1=DOV3DA2

	mov	eax, DWORD PTR _DOV3DA2+4
	mov	DWORD PTR _DOVN+4, eax
	mov	eax, DWORD PTR _DOV3DA2
	mov	DWORD PTR _DOVN, eax

; 665  : 	DOVnDA2=DOV3DA1

	mov	eax, DWORD PTR _DOV3DA1+4
	mov	DWORD PTR _DOVN+12, eax
	mov	eax, DWORD PTR _DOV3DA1
	mov	DWORD PTR _DOVN+8, eax

; 666  : 	DOVnDA5=DOV3DA5

	mov	eax, DWORD PTR _DOV3DA5+4
	mov	DWORD PTR _DOVN+20, eax
	mov	eax, DWORD PTR _DOV3DA5
	mov	DWORD PTR _DOVN+16, eax

; 667  : 	DOVnDAZ=DOV3DBZ

	mov	eax, DWORD PTR _DOV3DBZ+4
	mov	DWORD PTR _DOVN+28, eax
	mov	eax, DWORD PTR _DOV3DBZ
	mov	DWORD PTR _DOVN+24, eax

; 668  : 	DOVnDBZ=DOV3DAZ

	mov	eax, DWORD PTR _DOV3DAZ+4
	mov	DWORD PTR _DOVN+36, eax
	mov	eax, DWORD PTR _DOV3DAZ
	mov	DWORD PTR _DOVN+32, eax
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	esi, DWORD PTR [eax]

; 669  :               CALL ENUCL(DFLAG,OV3,RX2(L),RY2(L),RZ2(L),AA2(L),

	cmp	edi, esi
	jg	SHORT $L99544
	cmp	edi, 1
	jge	SHORT $L99545
$L99544:
	call	___FHranger@0
$L99545:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98380[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99546
	cmp	edi, 1
	jge	SHORT $L99547
$L99546:
	call	___FHranger@0
$L99547:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98381[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99548
	cmp	edi, 1
	jge	SHORT $L99549
$L99548:
	call	___FHranger@0
$L99549:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98382[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99550
	cmp	edi, 1
	jge	SHORT $L99551
$L99550:
	call	___FHranger@0
$L99551:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98383[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99552
	cmp	edi, 1
	jge	SHORT $L99553
$L99552:
	call	___FHranger@0
$L99553:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98384[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99554
	cmp	edi, 1
	jge	SHORT $L99555
$L99554:
	call	___FHranger@0
$L99555:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN15
	push	OFFSET FLAT:_EZ2
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	push	OFFSET FLAT:_EY2
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98385[esp+1408]
	push	OFFSET FLAT:_EX2
	cmp	edi, esi
	jg	SHORT $L99556
	cmp	edi, 1
	jge	SHORT $L99557
$L99556:
	call	___FHranger@0
$L99557:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99558
	cmp	edi, 1
	jge	SHORT $L99559
$L99558:
	call	___FHranger@0
$L99559:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99560
	cmp	edi, 1
	jge	SHORT $L99561
$L99560:
	call	___FHranger@0
$L99561:
	mov	eax, DWORD PTR _AA1$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98380[esp+1420]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98381[esp+1428]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98382[esp+1432]
	push	eax
	jg	SHORT $L99562
	cmp	edi, 1
	jge	SHORT $L99563
$L99562:
	call	___FHranger@0
$L99563:
	mov	eax, DWORD PTR _AA1$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99564
	cmp	edi, 1
	jge	SHORT $L99565
$L99564:
	call	___FHranger@0
$L99565:
	mov	eax, DWORD PTR _RZ1$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99566
	cmp	edi, 1
	jge	SHORT $L99567
$L99566:
	call	___FHranger@0
$L99567:
	mov	eax, DWORD PTR _RY1$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99568
	cmp	edi, 1
	jge	SHORT $L99569
$L99568:
	call	___FHranger@0
$L99569:
	mov	eax, DWORD PTR _RX1$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99570
	cmp	edi, 1
	jge	SHORT $L99571
$L99570:
	call	___FHranger@0
$L99571:
	mov	eax, DWORD PTR _AA2$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98383[esp+1452]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98384[esp+1460]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98385[esp+1464]
	push	eax
	jg	SHORT $L99572
	cmp	edi, 1
	jge	SHORT $L99573
$L99572:
	call	___FHranger@0
$L99573:
	mov	eax, DWORD PTR _AA2$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99574
	cmp	edi, 1
	jge	SHORT $L99575
$L99574:
	call	___FHranger@0
$L99575:
	mov	eax, DWORD PTR _RZ2$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99576
	cmp	edi, 1
	jge	SHORT $L99577
$L99576:
	call	___FHranger@0
$L99577:
	mov	eax, DWORD PTR _RY2$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99578
	cmp	edi, 1
	jge	SHORT $L99579
$L99578:
	call	___FHranger@0
$L99579:
	mov	eax, DWORD PTR _RX2$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV3
	push	eax
	call	_ENUCL@100
	mov	eax, DWORD PTR _DEN+52
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 675	; 000002a3H

; 670  :      1            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
; 671  :      2            RX1(L),RY1(L),RZ1(L),AA1(L),
; 672  :      3            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
; 673  :      4            B12(L),B12(K),
; 674  :      5            EX2,EY2,EZ2,EN15,KLM)
; 675  : 	D2EN7DA1=DENDA2P

	mov	DWORD PTR $T98573[esp+1404], 676	; 000002a4H
	mov	DWORD PTR $T98573[esp+1404], 677	; 000002a5H
	mov	DWORD PTR $T98573[esp+1404], 678	; 000002a6H
	mov	DWORD PTR $T98573[esp+1404], 679	; 000002a7H
	mov	DWORD PTR $T98573[esp+1404], 684	; 000002acH
	mov	DWORD PTR $T98573[esp+1404], 685	; 000002adH
	mov	DWORD PTR $T98573[esp+1404], 686	; 000002aeH
	mov	DWORD PTR $T98573[esp+1404], 687	; 000002afH
	mov	DWORD PTR $T98573[esp+1404], 688	; 000002b0H
	mov	DWORD PTR $T98573[esp+1404], 689	; 000002b1H
	mov	DWORD PTR _D2EN7DA1+4, eax
	mov	eax, DWORD PTR _DEN+48
	mov	DWORD PTR _D2EN7DA1, eax

; 676  : 	D2EN7DA2=DENDA1P

	mov	eax, DWORD PTR _DEN+44
	mov	DWORD PTR _D2EN7DA2+4, eax
	mov	eax, DWORD PTR _DEN+40
	mov	DWORD PTR _D2EN7DA2, eax

; 677  : 	D2EN7DA5=DENDA5P

	mov	eax, DWORD PTR _DEN+60
	mov	DWORD PTR _D2EN7DA5+4, eax
	mov	eax, DWORD PTR _DEN+56
	mov	DWORD PTR _D2EN7DA5, eax

; 678  : 	D2EN7DAZ=DENDBZP

	mov	eax, DWORD PTR _DEN+76
	mov	DWORD PTR _D2EN7DAZ+4, eax
	mov	eax, DWORD PTR _DEN+72
	mov	DWORD PTR _D2EN7DAZ, eax

; 679  : 	D2EN7DBZ=DENDAZP

	mov	eax, DWORD PTR _DEN+68
	mov	DWORD PTR _D2EN7DBZ+4, eax
	mov	eax, DWORD PTR _DEN+64
	mov	DWORD PTR _D2EN7DBZ, eax

; 680  : C
; 681  : C     ...Electro/Nuclear Symmetrization Integral:
; 682  : C         P_[elec.1,elec.2]*P_[H1,H2]
; 683  : C
; 684  : 	DOVnDA1=DOV4DA2

	mov	eax, DWORD PTR _DOV4DA2+4
	mov	DWORD PTR _DOVN+4, eax
	mov	eax, DWORD PTR _DOV4DA2
	mov	DWORD PTR _DOVN, eax

; 685  : 	DOVnDA2=DOV4DA1

	mov	eax, DWORD PTR _DOV4DA1+4
	mov	DWORD PTR _DOVN+12, eax
	mov	eax, DWORD PTR _DOV4DA1
	mov	DWORD PTR _DOVN+8, eax

; 686  : 	DOVnDA5=DOV4DA5

	mov	eax, DWORD PTR _DOV4DA5+4
	mov	DWORD PTR _DOVN+20, eax
	mov	eax, DWORD PTR _DOV4DA5
	mov	DWORD PTR _DOVN+16, eax

; 687  : 	DOVnDAZ=DOV4DBZ

	mov	eax, DWORD PTR _DOV4DBZ+4
	mov	DWORD PTR _DOVN+28, eax
	mov	eax, DWORD PTR _DOV4DBZ
	mov	DWORD PTR _DOVN+24, eax

; 688  : 	DOVnDBZ=DOV4DAZ

	mov	eax, DWORD PTR _DOV4DAZ+4
	mov	DWORD PTR _DOVN+36, eax
	mov	eax, DWORD PTR _DOV4DAZ
	mov	DWORD PTR _DOVN+32, eax
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	esi, DWORD PTR [eax]

; 689  :               CALL ENUCL(DFLAG,OV4,RX2(L),RY2(L),RZ2(L),AA2(L),

	cmp	edi, esi
	jg	SHORT $L99580
	cmp	edi, 1
	jge	SHORT $L99581
$L99580:
	call	___FHranger@0
$L99581:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98386[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99582
	cmp	edi, 1
	jge	SHORT $L99583
$L99582:
	call	___FHranger@0
$L99583:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98387[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99584
	cmp	edi, 1
	jge	SHORT $L99585
$L99584:
	call	___FHranger@0
$L99585:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98388[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99586
	cmp	edi, 1
	jge	SHORT $L99587
$L99586:
	call	___FHranger@0
$L99587:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98389[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99588
	cmp	edi, 1
	jge	SHORT $L99589
$L99588:
	call	___FHranger@0
$L99589:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98390[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99590
	cmp	edi, 1
	jge	SHORT $L99591
$L99590:
	call	___FHranger@0
$L99591:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN16
	push	OFFSET FLAT:_EZ2
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	push	OFFSET FLAT:_EY2
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98391[esp+1408]
	push	OFFSET FLAT:_EX2
	cmp	edi, esi
	jg	SHORT $L99592
	cmp	edi, 1
	jge	SHORT $L99593
$L99592:
	call	___FHranger@0
$L99593:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99594
	cmp	edi, 1
	jge	SHORT $L99595
$L99594:
	call	___FHranger@0
$L99595:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99596
	cmp	edi, 1
	jge	SHORT $L99597
$L99596:
	call	___FHranger@0
$L99597:
	mov	eax, DWORD PTR _AA2$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98386[esp+1420]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98387[esp+1428]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98388[esp+1432]
	push	eax
	jg	SHORT $L99598
	cmp	edi, 1
	jge	SHORT $L99599
$L99598:
	call	___FHranger@0
$L99599:
	mov	eax, DWORD PTR _AA1$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99600
	cmp	edi, 1
	jge	SHORT $L99601
$L99600:
	call	___FHranger@0
$L99601:
	mov	eax, DWORD PTR _RZ1$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99602
	cmp	edi, 1
	jge	SHORT $L99603
$L99602:
	call	___FHranger@0
$L99603:
	mov	eax, DWORD PTR _RY1$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99604
	cmp	edi, 1
	jge	SHORT $L99605
$L99604:
	call	___FHranger@0
$L99605:
	mov	eax, DWORD PTR _RX1$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99606
	cmp	edi, 1
	jge	SHORT $L99607
$L99606:
	call	___FHranger@0
$L99607:
	mov	eax, DWORD PTR _AA1$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98389[esp+1452]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98390[esp+1460]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98391[esp+1464]
	push	eax
	jg	SHORT $L99608
	cmp	edi, 1
	jge	SHORT $L99609
$L99608:
	call	___FHranger@0
$L99609:
	mov	eax, DWORD PTR _AA2$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99610
	cmp	edi, 1
	jge	SHORT $L99611
$L99610:
	call	___FHranger@0
$L99611:
	mov	eax, DWORD PTR _RZ2$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99612
	cmp	edi, 1
	jge	SHORT $L99613
$L99612:
	call	___FHranger@0
$L99613:
	mov	eax, DWORD PTR _RY2$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99614
	cmp	edi, 1
	jge	SHORT $L99615
$L99614:
	call	___FHranger@0
$L99615:
	mov	eax, DWORD PTR _RX2$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV4
	push	eax
	call	_ENUCL@100
	fld	QWORD PTR _EN9
	fadd	QWORD PTR _EN14
	fld	QWORD PTR _D2EN1DA1
	fadd	QWORD PTR _D2EN6DA1
	fxch	ST(1)
	fadd	QWORD PTR _EN6
	fld	QWORD PTR _D1EN6DA2
	fadd	QWORD PTR _D1EN4DA2
	fxch	ST(1)
	fadd	QWORD PTR _EN4
	fxch	ST(2)
	fadd	QWORD PTR _D1EN6DA1
	fxch	ST(1)
	fadd	QWORD PTR _D2EN3DA2
	fxch	ST(2)
	fadd	QWORD PTR _EN11
	fxch	ST(1)
	fadd	QWORD PTR _D1EN4DA1
	fxch	ST(2)
	fadd	QWORD PTR _D1EN8DA2
	fxch	ST(1)
	fadd	QWORD PTR _EN8
	fxch	ST(2)
	fadd	QWORD PTR _D2EN3DA1
	fxch	ST(1)
	fadd	QWORD PTR _D2EN5DA2
	fxch	ST(2)
	fadd	QWORD PTR _EN16
	fxch	ST(1)
	fadd	QWORD PTR _D1EN8DA1
	fxch	ST(2)
	fadd	QWORD PTR _D1EN3DA2
	fxch	ST(1)
	fadd	QWORD PTR _EN13
	fxch	ST(2)
	fadd	QWORD PTR _D2EN5DA1
	fxch	ST(1)
	fadd	QWORD PTR _D2EN2DA2
	fxch	ST(2)
	fadd	QWORD PTR _EN3
	fxch	ST(1)
	fadd	QWORD PTR _D1EN3DA1
	fxch	ST(2)
	fadd	QWORD PTR _D1EN1DA2
	fxch	ST(1)
	fadd	QWORD PTR _EN10
	fxch	ST(2)
	fadd	QWORD PTR _D2EN2DA1
	fxch	ST(1)
	fadd	QWORD PTR _DEN+40
	fxch	ST(2)
	fadd	QWORD PTR _EN15
	fxch	ST(1)
	fadd	QWORD PTR _DEN+48
	fxch	ST(2)
	fadd	QWORD PTR _D2EN7DA2
	fxch	ST(1)
	fadd	QWORD PTR _EN7
	fxch	ST(2)
	fadd	QWORD PTR _D2EN7DA1
	fxch	ST(1)
	fadd	QWORD PTR _D1EN7DA2
	fxch	ST(2)
	fadd	QWORD PTR _EN1
	fxch	ST(1)
	fadd	QWORD PTR _D1EN7DA1
	fxch	ST(2)
	fadd	QWORD PTR _D1EN5DA2
	fxch	ST(1)
	fadd	QWORD PTR _EN5
	fxch	ST(2)
	fadd	QWORD PTR _D1EN1DA1
	fxch	ST(1)
	fadd	QWORD PTR _D2EN4DA2
	fxch	ST(2)
	fadd	QWORD PTR _EN12
	fxch	ST(1)
	fadd	QWORD PTR _D1EN5DA1
	fxch	ST(2)
	fadd	QWORD PTR _D1EN2DA2
	fxch	ST(1)
	fadd	QWORD PTR _EN2
	fxch	ST(2)
	fadd	QWORD PTR _D2EN4DA1
	fxch	ST(1)
	fadd	QWORD PTR _D2EN1DA2
	fld	QWORD PTR _D1EN4DA5
	fadd	QWORD PTR _D2EN3DA5
	fxch	ST(2)
	fadd	QWORD PTR _D1EN2DA1
	fxch	ST(1)
	fadd	QWORD PTR _D2EN6DA2
	fxch	ST(2)
	fadd	QWORD PTR _D1EN8DA5
	fxch	ST(3)
	fstp	QWORD PTR _ENN1
	fstp	QWORD PTR _DENADA1
	fstp	QWORD PTR _DENADA2
	fadd	QWORD PTR _D2EN5DA5
	mov	DWORD PTR $T98573[esp+1404], 695	; 000002b7H

; 690  :      1            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
; 691  :      2            RX1(L),RY1(L),RZ1(L),AA1(L),
; 692  :      3            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
; 693  :      4            B12(L),B12(K),
; 694  :      5            EX2,EY2,EZ2,EN16,KLM)
; 695  : 	D2EN8DA1=DENDA2P

	mov	DWORD PTR $T98573[esp+1404], 696	; 000002b8H
	fadd	QWORD PTR _D1EN1DA5
	mov	DWORD PTR $T98573[esp+1404], 697	; 000002b9H
	mov	DWORD PTR $T98573[esp+1404], 698	; 000002baH
	fadd	QWORD PTR _D1EN3DA5
	mov	DWORD PTR $T98573[esp+1404], 699	; 000002bbH
	mov	DWORD PTR $T98573[esp+1404], 704	; 000002c0H
	fadd	QWORD PTR _D2EN2DA5
	mov	DWORD PTR $T98573[esp+1404], 711	; 000002c7H
	mov	DWORD PTR $T98573[esp+1404], 714	; 000002caH
	fadd	QWORD PTR _D2EN7DA5
	mov	DWORD PTR $T98573[esp+1404], 717	; 000002cdH
	mov	eax, DWORD PTR _DEN+52
	fadd	QWORD PTR _DEN+56
	mov	DWORD PTR _D2EN8DA1+4, eax
	mov	eax, DWORD PTR _DEN+48
	mov	DWORD PTR _D2EN8DA1, eax

; 696  : 	D2EN8DA2=DENDA1P

	mov	eax, DWORD PTR _DEN+44
	fadd	QWORD PTR _D1EN5DA5
	mov	DWORD PTR _D2EN8DA2+4, eax
	mov	eax, DWORD PTR _DEN+40
	mov	DWORD PTR _D2EN8DA2, eax

; 697  : 	D2EN8DA5=DENDA5P

	mov	eax, DWORD PTR _DEN+60
	fadd	QWORD PTR _D1EN7DA5
	mov	DWORD PTR _D2EN8DA5+4, eax
	mov	eax, DWORD PTR _DEN+56
	mov	DWORD PTR _D2EN8DA5, eax

; 698  : 	D2EN8DAZ=DENDBZP

	mov	eax, DWORD PTR _DEN+76
	fadd	QWORD PTR _D2EN4DA5
	mov	DWORD PTR _D2EN8DAZ+4, eax
	mov	eax, DWORD PTR _DEN+72
	mov	DWORD PTR _D2EN8DAZ, eax

; 699  : 	D2EN8DBZ=DENDAZP

	mov	eax, DWORD PTR _DEN+68
	fadd	QWORD PTR _D1EN2DA5
	mov	DWORD PTR _D2EN8DBZ+4, eax
	mov	eax, DWORD PTR _DEN+64
	mov	DWORD PTR _D2EN8DBZ, eax

; 700  : 	END IF
; 701  : C 
; 702  : C
; 703  : 	IF (ISYMM) THEN
; 704  :               ENN1=EN1+EN2+EN3+EN4+EN5+EN6+EN7+EN8
; 705  :      1            +EN9+EN10+EN11+EN12+EN13+EN14+EN15+EN16
; 706  : 	ELSE
; 707  :               ENN1=EN1+EN5+EN9+EN13
; 708  : 	END IF
; 709  : C
; 710  : 	IF (ISYMM) THEN
; 711  : 	DENADA1=D1EN1DA1+D1EN2DA1+D1EN3DA1+D1EN4DA1+D1EN5DA1+D1EN6DA1+D1EN
; 712  : 	-7DA1+D1EN8DA1+D2EN1DA1+D2EN2DA1+D2EN3DA1+D2EN4DA1+D2EN5DA1+D2EN6DA
; 713  :      -1+D2EN7DA1+D2EN8DA1
; 714  : 	DENADA2=D1EN1DA2+D1EN2DA2+D1EN3DA2+D1EN4DA2+D1EN5DA2+D1EN6DA2+D1EN
; 715  : 	-7DA2+D1EN8DA2+D2EN1DA2+D2EN2DA2+D2EN3DA2+D2EN4DA2+D2EN5DA2+D2EN6DA
; 716  :      -2+D2EN7DA2+D2EN8DA2
; 717  : 	DENADA5=D1EN1DA5+D1EN2DA5+D1EN3DA5+D1EN4DA5+D1EN5DA5+D1EN6DA5+D1EN

	fadd	QWORD PTR _D2EN1DA5
	fadd	QWORD PTR _D2EN6DA5
	fadd	QWORD PTR _D1EN6DA5
	fstp	QWORD PTR _DENADA5
	mov	DWORD PTR $T98573[esp+1404], 720	; 000002d0H

; 718  : 	-7DA5+D1EN8DA5+D2EN1DA5+D2EN2DA5+D2EN3DA5+D2EN4DA5+D2EN5DA5+D2EN6DA
; 719  :      -5+D2EN7DA5+D2EN8DA5
; 720  : 	DENADAZ=D1EN1DAZ+D1EN2DAZ+D1EN3DAZ+D1EN4DAZ+D1EN5DAZ+D1EN6DAZ+D1EN

	fld	QWORD PTR _D1EN4DAZ
	fadd	QWORD PTR _D2EN3DAZ
	fld	QWORD PTR _D1EN1DBZ
	fadd	QWORD PTR _D1EN8DBZ
	fxch	ST(1)
	fadd	QWORD PTR _D1EN8DAZ
	fxch	ST(1)
	fadd	QWORD PTR _D1EN3DBZ
	fxch	ST(1)
	fadd	QWORD PTR _D1EN1DAZ
	fxch	ST(1)
	fadd	QWORD PTR _D2EN2DBZ
	fxch	ST(1)
	fadd	QWORD PTR _D1EN3DAZ
	fxch	ST(1)
	fadd	QWORD PTR _D2EN5DBZ
	fxch	ST(1)
	fadd	QWORD PTR _D2EN2DAZ
	fxch	ST(1)
	fadd	QWORD PTR _D1EN5DBZ
	fxch	ST(1)
	fadd	QWORD PTR _D2EN5DAZ
	fxch	ST(1)
	fadd	QWORD PTR _D2EN7DBZ
	fxch	ST(1)
	fadd	QWORD PTR _D1EN5DAZ
	fxch	ST(1)
	fadd	QWORD PTR _D1EN7DBZ
	fxch	ST(1)
	fadd	QWORD PTR _D2EN7DAZ
	fxch	ST(1)
	fadd	QWORD PTR _DEN+64
	fxch	ST(1)
	fadd	QWORD PTR _D1EN7DAZ
	fxch	ST(1)
	fadd	QWORD PTR _D2EN4DBZ
	fxch	ST(1)
	fadd	QWORD PTR _DEN+72
	fxch	ST(1)
	fadd	QWORD PTR _D1EN2DBZ
	fxch	ST(1)
	fadd	QWORD PTR _D2EN4DAZ
	fxch	ST(1)
	fadd	QWORD PTR _D2EN1DBZ
	fxch	ST(1)
	fadd	QWORD PTR _D1EN2DAZ
	fxch	ST(1)
	fadd	QWORD PTR _D2EN6DBZ
	fxch	ST(1)
	fadd	QWORD PTR _D2EN1DAZ
	fxch	ST(1)
	fadd	QWORD PTR _D1EN6DBZ
	fxch	ST(1)
	fadd	QWORD PTR _D2EN6DAZ
	fxch	ST(1)
	fadd	QWORD PTR _D1EN4DBZ
	fxch	ST(1)
	fadd	QWORD PTR _D1EN6DAZ
	fxch	ST(1)
	fadd	QWORD PTR _D2EN3DBZ
	mov	eax, DWORD PTR _DOV1DA1+4
	push	OFFSET FLAT:_KLM
	fxch	ST(1)
	fstp	QWORD PTR _DENADAZ
	fstp	QWORD PTR _DENADBZ
	push	OFFSET FLAT:_ER1
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1412], 723	; 000002d3H

; 721  : 	-7DAZ+D1EN8DAZ+D2EN1DAZ+D2EN2DAZ+D2EN3DAZ+D2EN4DAZ+D2EN5DAZ+D2EN6DA
; 722  :      -Z+D2EN7DAZ+D2EN8DAZ
; 723  : 	DENADBZ=D1EN1DBZ+D1EN2DBZ+D1EN3DBZ+D1EN4DBZ+D1EN5DBZ+D1EN6DBZ+D1EN

	mov	DWORD PTR $T98573[esp+1412], 764	; 000002fcH

; 724  : 	-7DBZ+D1EN8DBZ+D2EN1DBZ+D2EN2DBZ+D2EN3DBZ+D2EN4DBZ+D2EN5DBZ+D2EN6DB
; 725  :      -Z+D2EN7DBZ+D2EN8DBZ
; 726  : 	ELSE
; 727  : 	DENADA1=D1EN1DA1+D1EN5DA1+D2EN1DA1+D2EN5DA1
; 728  : 	DENADA2=D1EN1DA2+D1EN5DA2+D2EN1DA2+D2EN5DA2
; 729  : 	DENADA5=D1EN1DA5+D1EN5DA5+D2EN1DA5+D2EN5DA5
; 730  : 	DENADAZ=D1EN1DAZ+D1EN5DAZ+D2EN1DAZ+D2EN5DAZ
; 731  : 	DENADBZ=D1EN1DBZ+D1EN5DBZ+D2EN1DBZ+D2EN5DBZ
; 732  : 	END IF
; 733  : 	IF(IPERM.EQ.1) THEN
; 734  : 	DENADA1=D1EN1DA1+D1EN5DA1+D2EN1DA1+D2EN5DA1
; 735  : 	DENADA2=D1EN1DA2+D1EN5DA2+D2EN1DA2+D2EN5DA2
; 736  : 	DENADA5=D1EN1DA5+D1EN5DA5+D2EN1DA5+D2EN5DA5
; 737  : 	DENADAZ=D1EN1DAZ+D1EN5DAZ+D2EN1DAZ+D2EN5DAZ
; 738  : 	DENADBZ=D1EN1DBZ+D1EN5DBZ+D2EN1DBZ+D2EN5DBZ
; 739  : 	END IF
; 740  : 	IF(IPERM.EQ.2) THEN
; 741  : 	DENADA1=D1EN2DA1+D1EN6DA1+D2EN2DA1+D2EN6DA1
; 742  : 	DENADA2=D1EN2DA2+D1EN6DA2+D2EN2DA2+D2EN6DA2
; 743  : 	DENADA5=D1EN2DA5+D1EN6DA5+D2EN2DA5+D2EN6DA5
; 744  : 	DENADAZ=D1EN2DAZ+D1EN6DAZ+D2EN2DAZ+D2EN6DAZ
; 745  : 	DENADBZ=D1EN2DBZ+D1EN6DBZ+D2EN2DBZ+D2EN6DBZ
; 746  : 	END IF
; 747  : 	IF(IPERM.EQ.3) THEN
; 748  : 	DENADA1=D1EN3DA1+D1EN7DA1+D2EN3DA1+D2EN7DA1
; 749  : 	DENADA2=D1EN3DA2+D1EN7DA2+D2EN3DA2+D2EN7DA2
; 750  : 	DENADA5=D1EN3DA5+D1EN7DA5+D2EN3DA5+D2EN7DA5
; 751  : 	DENADAZ=D1EN3DAZ+D1EN7DAZ+D2EN3DAZ+D2EN7DAZ
; 752  : 	DENADBZ=D1EN3DBZ+D1EN7DBZ+D2EN3DBZ+D2EN7DBZ
; 753  : 	END IF
; 754  : 	IF(IPERM.EQ.4) THEN
; 755  : 	DENADA1=D1EN4DA1+D1EN8DA1+D2EN4DA1+D2EN8DA1
; 756  : 	DENADA2=D1EN4DA2+D1EN8DA2+D2EN4DA2+D2EN8DA2
; 757  : 	DENADA5=D1EN4DA5+D1EN8DA5+D2EN4DA5+D2EN8DA5
; 758  : 	DENADAZ=D1EN4DAZ+D1EN8DAZ+D2EN4DAZ+D2EN8DAZ
; 759  : 	DENADBZ=D1EN4DBZ+D1EN8DBZ+D2EN4DBZ+D2EN8DBZ
; 760  : 	END IF
; 761  : C     -----------------------------------------------------
; 762  : C     --- Calculate the Elec./Elec. Repulsion ---
; 763  : C     -----------------------------------------------------
; 764  : 	DOVnDA1=DOV1DA1

	mov	DWORD PTR $T98573[esp+1412], 765	; 000002fdH
	mov	DWORD PTR $T98573[esp+1412], 766	; 000002feH
	mov	DWORD PTR $T98573[esp+1412], 767	; 000002ffH
	mov	DWORD PTR $T98573[esp+1412], 768	; 00000300H
	mov	DWORD PTR $T98573[esp+1412], 769	; 00000301H
	mov	DWORD PTR _DOVN+4, eax
	mov	eax, DWORD PTR _DOV1DA1
	mov	DWORD PTR _DOVN, eax

; 765  : 	DOVnDA2=DOV1DA2

	mov	eax, DWORD PTR _DOV1DA2+4
	mov	DWORD PTR _DOVN+12, eax
	mov	eax, DWORD PTR _DOV1DA2
	mov	DWORD PTR _DOVN+8, eax

; 766  : 	DOVnDA5=DOV1DA5

	mov	eax, DWORD PTR _DOV1DA5+4
	mov	DWORD PTR _DOVN+20, eax
	mov	eax, DWORD PTR _DOV1DA5
	mov	DWORD PTR _DOVN+16, eax

; 767  : 	DOVnDAZ=DOV1DAZ

	mov	eax, DWORD PTR _DOV1DAZ+4
	mov	DWORD PTR _DOVN+28, eax
	mov	eax, DWORD PTR _DOV1DAZ
	mov	DWORD PTR _DOVN+24, eax

; 768  : 	DOVnDBZ=DOV1DBZ

	mov	eax, DWORD PTR _DOV1DBZ+4
	mov	DWORD PTR _DOVN+36, eax
	mov	eax, DWORD PTR _DOV1DBZ
	mov	DWORD PTR _DOVN+32, eax
	mov	eax, DWORD PTR _M$[esp+1396]
	mov	esi, DWORD PTR [eax]

; 769  :               CALL ELREP(DFLAG,OV1,RX1(L),RY1(L),RZ1(L),AA1(L),

	cmp	edi, esi
	jg	SHORT $L99616
	cmp	edi, 1
	jge	SHORT $L99617
$L99616:
	call	___FHranger@0
$L99617:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99618
	cmp	edi, 1
	jge	SHORT $L99619
$L99618:
	call	___FHranger@0
$L99619:
	mov	eax, DWORD PTR _B12$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99620
	cmp	edi, 1
	jge	SHORT $L99621
$L99620:
	call	___FHranger@0
$L99621:
	mov	eax, DWORD PTR _AA2$[esp+1404]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99622
	cmp	edi, 1
	jge	SHORT $L99623
$L99622:
	call	___FHranger@0
$L99623:
	mov	eax, DWORD PTR _RZ2$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99624
	cmp	edi, 1
	jge	SHORT $L99625
$L99624:
	call	___FHranger@0
$L99625:
	mov	eax, DWORD PTR _RY2$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99626
	cmp	edi, 1
	jge	SHORT $L99627
$L99626:
	call	___FHranger@0
$L99627:
	mov	eax, DWORD PTR _RX2$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99628
	cmp	edi, 1
	jge	SHORT $L99629
$L99628:
	call	___FHranger@0
$L99629:
	mov	eax, DWORD PTR _AA2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99630
	cmp	edi, 1
	jge	SHORT $L99631
$L99630:
	call	___FHranger@0
$L99631:
	mov	eax, DWORD PTR _RZ2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99632
	cmp	edi, 1
	jge	SHORT $L99633
$L99632:
	call	___FHranger@0
$L99633:
	mov	eax, DWORD PTR _RY2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99634
	cmp	edi, 1
	jge	SHORT $L99635
$L99634:
	call	___FHranger@0
$L99635:
	mov	eax, DWORD PTR _RX2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99636
	cmp	edi, 1
	jge	SHORT $L99637
$L99636:
	call	___FHranger@0
$L99637:
	mov	eax, DWORD PTR _AA1$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99638
	cmp	edi, 1
	jge	SHORT $L99639
$L99638:
	call	___FHranger@0
$L99639:
	mov	eax, DWORD PTR _RZ1$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99640
	cmp	edi, 1
	jge	SHORT $L99641
$L99640:
	call	___FHranger@0
$L99641:
	mov	eax, DWORD PTR _RY1$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99642
	cmp	edi, 1
	jge	SHORT $L99643
$L99642:
	call	___FHranger@0
$L99643:
	mov	eax, DWORD PTR _RX1$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99644
	cmp	edi, 1
	jge	SHORT $L99645
$L99644:
	call	___FHranger@0
$L99645:
	mov	eax, DWORD PTR _AA1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99646
	cmp	edi, 1
	jge	SHORT $L99647
$L99646:
	call	___FHranger@0
$L99647:
	mov	eax, DWORD PTR _RZ1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99648
	cmp	edi, 1
	jge	SHORT $L99649
$L99648:
	call	___FHranger@0
$L99649:
	mov	eax, DWORD PTR _RY1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99650
	cmp	edi, 1
	jge	SHORT $L99651
$L99650:
	call	___FHranger@0
$L99651:
	mov	eax, DWORD PTR _RX1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1464]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV1
	push	eax
	call	_ELREP@88
	mov	eax, DWORD PTR _DER+4
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_ER2
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1412], 774	; 00000306H

; 770  :      1           RX1(K),RY1(K),RZ1(K),AA1(K),                      
; 771  :      2           RX2(L),RY2(L),RZ2(L),AA2(L),
; 772  :      3           RX2(K),RY2(K),RZ2(K),AA2(K),
; 773  :      4           B12(L),B12(K),ER1,KLM)  
; 774  : 	DER1DA1=DERDA1

	mov	DWORD PTR $T98573[esp+1412], 775	; 00000307H
	mov	DWORD PTR $T98573[esp+1412], 776	; 00000308H
	mov	DWORD PTR $T98573[esp+1412], 777	; 00000309H
	mov	DWORD PTR $T98573[esp+1412], 778	; 0000030aH
	mov	DWORD PTR $T98573[esp+1412], 783	; 0000030fH
	mov	DWORD PTR $T98573[esp+1412], 784	; 00000310H
	mov	DWORD PTR $T98573[esp+1412], 785	; 00000311H
	mov	DWORD PTR $T98573[esp+1412], 786	; 00000312H
	mov	DWORD PTR $T98573[esp+1412], 787	; 00000313H
	mov	DWORD PTR $T98573[esp+1412], 788	; 00000314H
	mov	DWORD PTR _DER1DA1+4, eax
	mov	eax, DWORD PTR _DER
	mov	DWORD PTR _DER1DA1, eax

; 775  : 	DER1DA2=DERDA2

	mov	eax, DWORD PTR _DER+12
	mov	DWORD PTR _DER1DA2+4, eax
	mov	eax, DWORD PTR _DER+8
	mov	DWORD PTR _DER1DA2, eax

; 776  : 	DER1DA5=DERDA5

	mov	eax, DWORD PTR _DER+20
	mov	DWORD PTR _DER1DA5+4, eax
	mov	eax, DWORD PTR _DER+16
	mov	DWORD PTR _DER1DA5, eax

; 777  : 	DER1DAZ=DERDAZ

	mov	eax, DWORD PTR _DER+28
	mov	DWORD PTR _DER1DAZ+4, eax
	mov	eax, DWORD PTR _DER+24
	mov	DWORD PTR _DER1DAZ, eax

; 778  : 	DER1DBZ=DERDBZ

	mov	eax, DWORD PTR _DER+36
	mov	DWORD PTR _DER1DBZ+4, eax
	mov	eax, DWORD PTR _DER+32
	mov	DWORD PTR _DER1DBZ, eax

; 779  : C
; 780  : C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
; 781  : C
; 782  : 	IF (ISYMM) THEN
; 783  : 	DOVnDA1=DOV2DA1

	mov	eax, DWORD PTR _DOV2DA1+4
	mov	DWORD PTR _DOVN+4, eax
	mov	eax, DWORD PTR _DOV2DA1
	mov	DWORD PTR _DOVN, eax

; 784  : 	DOVnDA2=DOV2DA2

	mov	eax, DWORD PTR _DOV2DA2+4
	mov	DWORD PTR _DOVN+12, eax
	mov	eax, DWORD PTR _DOV2DA2
	mov	DWORD PTR _DOVN+8, eax

; 785  : 	DOVnDA5=DOV2DA5

	mov	eax, DWORD PTR _DOV2DA5+4
	mov	DWORD PTR _DOVN+20, eax
	mov	eax, DWORD PTR _DOV2DA5
	mov	DWORD PTR _DOVN+16, eax

; 786  : 	DOVnDAZ=DOV2DAZ

	mov	eax, DWORD PTR _DOV2DAZ+4
	mov	DWORD PTR _DOVN+28, eax
	mov	eax, DWORD PTR _DOV2DAZ
	mov	DWORD PTR _DOVN+24, eax

; 787  : 	DOVnDBZ=DOV2DBZ

	mov	eax, DWORD PTR _DOV2DBZ+4
	mov	DWORD PTR _DOVN+36, eax
	mov	eax, DWORD PTR _DOV2DBZ
	mov	DWORD PTR _DOVN+32, eax
	mov	eax, DWORD PTR _M$[esp+1396]
	mov	esi, DWORD PTR [eax]

; 788  :               CALL ELREP(DFLAG,OV2,RX1(L),RY1(L),RZ1(L),AA1(L),

	cmp	edi, esi
	jg	SHORT $L99652
	cmp	edi, 1
	jge	SHORT $L99653
$L99652:
	call	___FHranger@0
$L99653:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99654
	cmp	edi, 1
	jge	SHORT $L99655
$L99654:
	call	___FHranger@0
$L99655:
	mov	eax, DWORD PTR _B12$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99656
	cmp	edi, 1
	jge	SHORT $L99657
$L99656:
	call	___FHranger@0
$L99657:
	mov	eax, DWORD PTR _AA1$[esp+1404]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99658
	cmp	edi, 1
	jge	SHORT $L99659
$L99658:
	call	___FHranger@0
$L99659:
	mov	eax, DWORD PTR _RZ1$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99660
	cmp	edi, 1
	jge	SHORT $L99661
$L99660:
	call	___FHranger@0
$L99661:
	mov	eax, DWORD PTR _RY1$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99662
	cmp	edi, 1
	jge	SHORT $L99663
$L99662:
	call	___FHranger@0
$L99663:
	mov	eax, DWORD PTR _RX1$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99664
	cmp	edi, 1
	jge	SHORT $L99665
$L99664:
	call	___FHranger@0
$L99665:
	mov	eax, DWORD PTR _AA2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99666
	cmp	edi, 1
	jge	SHORT $L99667
$L99666:
	call	___FHranger@0
$L99667:
	mov	eax, DWORD PTR _RZ2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99668
	cmp	edi, 1
	jge	SHORT $L99669
$L99668:
	call	___FHranger@0
$L99669:
	mov	eax, DWORD PTR _RY2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99670
	cmp	edi, 1
	jge	SHORT $L99671
$L99670:
	call	___FHranger@0
$L99671:
	mov	eax, DWORD PTR _RX2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99672
	cmp	edi, 1
	jge	SHORT $L99673
$L99672:
	call	___FHranger@0
$L99673:
	mov	eax, DWORD PTR _AA2$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99674
	cmp	edi, 1
	jge	SHORT $L99675
$L99674:
	call	___FHranger@0
$L99675:
	mov	eax, DWORD PTR _RZ2$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99676
	cmp	edi, 1
	jge	SHORT $L99677
$L99676:
	call	___FHranger@0
$L99677:
	mov	eax, DWORD PTR _RY2$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99678
	cmp	edi, 1
	jge	SHORT $L99679
$L99678:
	call	___FHranger@0
$L99679:
	mov	eax, DWORD PTR _RX2$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99680
	cmp	edi, 1
	jge	SHORT $L99681
$L99680:
	call	___FHranger@0
$L99681:
	mov	eax, DWORD PTR _AA1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99682
	cmp	edi, 1
	jge	SHORT $L99683
$L99682:
	call	___FHranger@0
$L99683:
	mov	eax, DWORD PTR _RZ1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99684
	cmp	edi, 1
	jge	SHORT $L99685
$L99684:
	call	___FHranger@0
$L99685:
	mov	eax, DWORD PTR _RY1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99686
	cmp	edi, 1
	jge	SHORT $L99687
$L99686:
	call	___FHranger@0
$L99687:
	mov	eax, DWORD PTR _RX1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1464]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV2
	push	eax
	call	_ELREP@88
	mov	eax, DWORD PTR _DER+4
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 793	; 00000319H

; 789  :      1           RX2(K),RY2(K),RZ2(K),AA2(K),                      
; 790  :      2           RX2(L),RY2(L),RZ2(L),AA2(L),
; 791  :      3           RX1(K),RY1(K),RZ1(K),AA1(K),
; 792  :      4           B12(L),B12(K),ER2,KLM)  
; 793  : 	DER2DA1=DERDA1

	mov	DWORD PTR $T98573[esp+1404], 794	; 0000031aH
	mov	DWORD PTR $T98573[esp+1404], 795	; 0000031bH
	mov	DWORD PTR $T98573[esp+1404], 796	; 0000031cH
	mov	DWORD PTR $T98573[esp+1404], 797	; 0000031dH
	mov	DWORD PTR $T98573[esp+1404], 801	; 00000321H
	mov	DWORD PTR $T98573[esp+1404], 802	; 00000322H
	mov	DWORD PTR $T98573[esp+1404], 803	; 00000323H
	mov	DWORD PTR $T98573[esp+1404], 804	; 00000324H
	mov	DWORD PTR $T98573[esp+1404], 805	; 00000325H
	mov	DWORD PTR $T98573[esp+1404], 806	; 00000326H
	mov	DWORD PTR _DER2DA1+4, eax
	mov	eax, DWORD PTR _DER
	mov	DWORD PTR _DER2DA1, eax

; 794  : 	DER2DA2=DERDA2

	mov	eax, DWORD PTR _DER+12
	mov	DWORD PTR _DER2DA2+4, eax
	mov	eax, DWORD PTR _DER+8
	mov	DWORD PTR _DER2DA2, eax

; 795  : 	DER2DA5=DERDA5

	mov	eax, DWORD PTR _DER+20
	mov	DWORD PTR _DER2DA5+4, eax
	mov	eax, DWORD PTR _DER+16
	mov	DWORD PTR _DER2DA5, eax

; 796  : 	DER2DAZ=DERDAZ

	mov	eax, DWORD PTR _DER+28
	mov	DWORD PTR _DER2DAZ+4, eax
	mov	eax, DWORD PTR _DER+24
	mov	DWORD PTR _DER2DAZ, eax

; 797  : 	DER2DBZ=DERDBZ

	mov	eax, DWORD PTR _DER+36
	mov	DWORD PTR _DER2DBZ+4, eax
	mov	eax, DWORD PTR _DER+32
	mov	DWORD PTR _DER2DBZ, eax

; 798  : C
; 799  : C     ...Nuclear Symmetrization Overlap Integral: P_[H1,H2]
; 800  : C
; 801  : 	DOVnDA1=DOV3DA1

	mov	eax, DWORD PTR _DOV3DA1+4
	mov	DWORD PTR _DOVN+4, eax
	mov	eax, DWORD PTR _DOV3DA1
	mov	DWORD PTR _DOVN, eax

; 802  : 	DOVnDA2=DOV3DA2

	mov	eax, DWORD PTR _DOV3DA2+4
	mov	DWORD PTR _DOVN+12, eax
	mov	eax, DWORD PTR _DOV3DA2
	mov	DWORD PTR _DOVN+8, eax

; 803  : 	DOVnDA5=DOV3DA5

	mov	eax, DWORD PTR _DOV3DA5+4
	mov	DWORD PTR _DOVN+20, eax
	mov	eax, DWORD PTR _DOV3DA5
	mov	DWORD PTR _DOVN+16, eax

; 804  : 	DOVnDAZ=DOV3DAZ

	mov	eax, DWORD PTR _DOV3DAZ+4
	mov	DWORD PTR _DOVN+28, eax
	mov	eax, DWORD PTR _DOV3DAZ
	mov	DWORD PTR _DOVN+24, eax

; 805  : 	DOVnDBZ=DOV3DBZ

	mov	eax, DWORD PTR _DOV3DBZ+4
	mov	DWORD PTR _DOVN+36, eax
	mov	eax, DWORD PTR _DOV3DBZ
	mov	DWORD PTR _DOVN+32, eax
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	esi, DWORD PTR [eax]

; 806  :               CALL ELREP(DFLAG,OV3,RX1(L),RY1(L),RZ1(L),AA1(L),

	cmp	edi, esi
	jg	SHORT $L99688
	cmp	edi, 1
	jge	SHORT $L99689
$L99688:
	call	___FHranger@0
$L99689:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98392[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99690
	cmp	edi, 1
	jge	SHORT $L99691
$L99690:
	call	___FHranger@0
$L99691:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98393[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99692
	cmp	edi, 1
	jge	SHORT $L99693
$L99692:
	call	___FHranger@0
$L99693:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98394[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99694
	cmp	edi, 1
	jge	SHORT $L99695
$L99694:
	call	___FHranger@0
$L99695:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98395[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99696
	cmp	edi, 1
	jge	SHORT $L99697
$L99696:
	call	___FHranger@0
$L99697:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98396[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99698
	cmp	edi, 1
	jge	SHORT $L99699
$L99698:
	call	___FHranger@0
$L99699:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_ER3
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98397[esp+1400]
	cmp	edi, esi
	jg	SHORT $L99700
	cmp	edi, 1
	jge	SHORT $L99701
$L99700:
	call	___FHranger@0
$L99701:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99702
	cmp	edi, 1
	jge	SHORT $L99703
$L99702:
	call	___FHranger@0
$L99703:
	mov	eax, DWORD PTR _B12$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99704
	cmp	edi, 1
	jge	SHORT $L99705
$L99704:
	call	___FHranger@0
$L99705:
	mov	eax, DWORD PTR _AA2$[esp+1404]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98392[esp+1408]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98393[esp+1416]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98394[esp+1420]
	push	eax
	jg	SHORT $L99706
	cmp	edi, 1
	jge	SHORT $L99707
$L99706:
	call	___FHranger@0
$L99707:
	mov	eax, DWORD PTR _AA2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99708
	cmp	edi, 1
	jge	SHORT $L99709
$L99708:
	call	___FHranger@0
$L99709:
	mov	eax, DWORD PTR _RZ2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99710
	cmp	edi, 1
	jge	SHORT $L99711
$L99710:
	call	___FHranger@0
$L99711:
	mov	eax, DWORD PTR _RY2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99712
	cmp	edi, 1
	jge	SHORT $L99713
$L99712:
	call	___FHranger@0
$L99713:
	mov	eax, DWORD PTR _RX2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99714
	cmp	edi, 1
	jge	SHORT $L99715
$L99714:
	call	___FHranger@0
$L99715:
	mov	eax, DWORD PTR _AA1$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98395[esp+1440]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98396[esp+1448]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98397[esp+1452]
	push	eax
	jg	SHORT $L99716
	cmp	edi, 1
	jge	SHORT $L99717
$L99716:
	call	___FHranger@0
$L99717:
	mov	eax, DWORD PTR _AA1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99718
	cmp	edi, 1
	jge	SHORT $L99719
$L99718:
	call	___FHranger@0
$L99719:
	mov	eax, DWORD PTR _RZ1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99720
	cmp	edi, 1
	jge	SHORT $L99721
$L99720:
	call	___FHranger@0
$L99721:
	mov	eax, DWORD PTR _RY1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99722
	cmp	edi, 1
	jge	SHORT $L99723
$L99722:
	call	___FHranger@0
$L99723:
	mov	eax, DWORD PTR _RX1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1464]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV3
	push	eax
	call	_ELREP@88
	mov	eax, DWORD PTR _DER+4
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 811	; 0000032bH

; 807  :      1            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
; 808  :      2            RX2(L),RY2(L),RZ2(L),AA2(L),
; 809  :      3            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
; 810  :      4            B12(L),B12(K),ER3,KLM)
; 811  : 	DER3DA1=DERDA1

	mov	DWORD PTR $T98573[esp+1404], 812	; 0000032cH
	mov	DWORD PTR $T98573[esp+1404], 813	; 0000032dH
	mov	DWORD PTR $T98573[esp+1404], 814	; 0000032eH
	mov	DWORD PTR $T98573[esp+1404], 815	; 0000032fH
	mov	DWORD PTR $T98573[esp+1404], 820	; 00000334H
	mov	DWORD PTR $T98573[esp+1404], 821	; 00000335H
	mov	DWORD PTR $T98573[esp+1404], 822	; 00000336H
	mov	DWORD PTR $T98573[esp+1404], 823	; 00000337H
	mov	DWORD PTR $T98573[esp+1404], 824	; 00000338H
	mov	DWORD PTR $T98573[esp+1404], 825	; 00000339H
	mov	DWORD PTR _DER3DA1+4, eax
	mov	eax, DWORD PTR _DER
	mov	DWORD PTR _DER3DA1, eax

; 812  : 	DER3DA2=DERDA2

	mov	eax, DWORD PTR _DER+12
	mov	DWORD PTR _DER3DA2+4, eax
	mov	eax, DWORD PTR _DER+8
	mov	DWORD PTR _DER3DA2, eax

; 813  : 	DER3DA5=DERDA5

	mov	eax, DWORD PTR _DER+20
	mov	DWORD PTR _DER3DA5+4, eax
	mov	eax, DWORD PTR _DER+16
	mov	DWORD PTR _DER3DA5, eax

; 814  : 	DER3DAZ=DERDAZ

	mov	eax, DWORD PTR _DER+28
	mov	DWORD PTR _DER3DAZ+4, eax
	mov	eax, DWORD PTR _DER+24
	mov	DWORD PTR _DER3DAZ, eax

; 815  : 	DER3DBZ=DERDBZ

	mov	eax, DWORD PTR _DER+36
	mov	DWORD PTR _DER3DBZ+4, eax
	mov	eax, DWORD PTR _DER+32
	mov	DWORD PTR _DER3DBZ, eax

; 816  : C
; 817  : C     ...Electro/Nuclear Symmetrization Overlap Integral:
; 818  : C         P_[elec.1,elec.2]*P_[H1,H2]
; 819  : C
; 820  : 	DOVnDA1=DOV4DA1

	mov	eax, DWORD PTR _DOV4DA1+4
	mov	DWORD PTR _DOVN+4, eax
	mov	eax, DWORD PTR _DOV4DA1
	mov	DWORD PTR _DOVN, eax

; 821  : 	DOVnDA2=DOV4DA2

	mov	eax, DWORD PTR _DOV4DA2+4
	mov	DWORD PTR _DOVN+12, eax
	mov	eax, DWORD PTR _DOV4DA2
	mov	DWORD PTR _DOVN+8, eax

; 822  : 	DOVnDA5=DOV4DA5

	mov	eax, DWORD PTR _DOV4DA5+4
	mov	DWORD PTR _DOVN+20, eax
	mov	eax, DWORD PTR _DOV4DA5
	mov	DWORD PTR _DOVN+16, eax

; 823  : 	DOVnDAZ=DOV4DAZ

	mov	eax, DWORD PTR _DOV4DAZ+4
	mov	DWORD PTR _DOVN+28, eax
	mov	eax, DWORD PTR _DOV4DAZ
	mov	DWORD PTR _DOVN+24, eax

; 824  : 	DOVnDBZ=DOV4DBZ

	mov	eax, DWORD PTR _DOV4DBZ+4
	mov	DWORD PTR _DOVN+36, eax
	mov	eax, DWORD PTR _DOV4DBZ
	mov	DWORD PTR _DOVN+32, eax
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	esi, DWORD PTR [eax]

; 825  :               CALL ELREP(DFLAG,OV4,RX1(L),RY1(L),RZ1(L),AA1(L),

	cmp	edi, esi
	jg	SHORT $L99724
	cmp	edi, 1
	jge	SHORT $L99725
$L99724:
	call	___FHranger@0
$L99725:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98398[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99726
	cmp	edi, 1
	jge	SHORT $L99727
$L99726:
	call	___FHranger@0
$L99727:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98399[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99728
	cmp	edi, 1
	jge	SHORT $L99729
$L99728:
	call	___FHranger@0
$L99729:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98400[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99730
	cmp	edi, 1
	jge	SHORT $L99731
$L99730:
	call	___FHranger@0
$L99731:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98401[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99732
	cmp	edi, 1
	jge	SHORT $L99733
$L99732:
	call	___FHranger@0
$L99733:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98402[esp+1392]
	cmp	edi, esi
	jg	SHORT $L99734
	cmp	edi, 1
	jge	SHORT $L99735
$L99734:
	call	___FHranger@0
$L99735:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_ER4
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98403[esp+1400]
	cmp	edi, esi
	jg	SHORT $L99736
	cmp	edi, 1
	jge	SHORT $L99737
$L99736:
	call	___FHranger@0
$L99737:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99738
	cmp	edi, 1
	jge	SHORT $L99739
$L99738:
	call	___FHranger@0
$L99739:
	mov	eax, DWORD PTR _B12$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99740
	cmp	edi, 1
	jge	SHORT $L99741
$L99740:
	call	___FHranger@0
$L99741:
	mov	eax, DWORD PTR _AA1$[esp+1404]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98398[esp+1408]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98399[esp+1416]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98400[esp+1420]
	push	eax
	jg	SHORT $L99742
	cmp	edi, 1
	jge	SHORT $L99743
$L99742:
	call	___FHranger@0
$L99743:
	mov	eax, DWORD PTR _AA2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99744
	cmp	edi, 1
	jge	SHORT $L99745
$L99744:
	call	___FHranger@0
$L99745:
	mov	eax, DWORD PTR _RZ2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99746
	cmp	edi, 1
	jge	SHORT $L99747
$L99746:
	call	___FHranger@0
$L99747:
	mov	eax, DWORD PTR _RY2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99748
	cmp	edi, 1
	jge	SHORT $L99749
$L99748:
	call	___FHranger@0
$L99749:
	mov	eax, DWORD PTR _RX2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L99750
	cmp	edi, 1
	jge	SHORT $L99751
$L99750:
	call	___FHranger@0
$L99751:
	mov	eax, DWORD PTR _AA2$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98401[esp+1440]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98402[esp+1448]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98403[esp+1452]
	push	eax
	jg	SHORT $L99752
	cmp	edi, 1
	jge	SHORT $L99753
$L99752:
	call	___FHranger@0
$L99753:
	mov	eax, DWORD PTR _AA1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99754
	cmp	edi, 1
	jge	SHORT $L99755
$L99754:
	call	___FHranger@0
$L99755:
	mov	eax, DWORD PTR _RZ1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99756
	cmp	edi, 1
	jge	SHORT $L99757
$L99756:
	call	___FHranger@0
$L99757:
	mov	eax, DWORD PTR _RY1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L99758
	cmp	edi, 1
	jge	SHORT $L99759
$L99758:
	call	___FHranger@0
$L99759:
	mov	eax, DWORD PTR _RX1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1464]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV4
	push	eax
	call	_ELREP@88
	fld	QWORD PTR _ER2
	fadd	QWORD PTR _ER4
	fld	QWORD PTR _DER3DA1
	fadd	QWORD PTR _DER
	fxch	ST(1)
	fadd	QWORD PTR _ER3
	fld	QWORD PTR _DER1DA2
	fadd	QWORD PTR _DER3DA2
	fxch	ST(1)
	fadd	QWORD PTR _ER1
	fxch	ST(2)
	fadd	QWORD PTR _DER1DA1
	fxch	ST(1)
	fadd	QWORD PTR _DER+8
	fld	QWORD PTR _DER1DA5
	fadd	QWORD PTR _DER3DA5
	fxch	ST(2)
	fadd	QWORD PTR _DER2DA1
	fxch	ST(1)
	fadd	QWORD PTR _DER2DA2
	fxch	ST(2)
	fadd	QWORD PTR _DER+16
	fxch	ST(3)
	fstp	QWORD PTR -1360+[esp+1392]
	fstp	QWORD PTR _DEERDA1
	fstp	QWORD PTR _DEERDA2
	fadd	QWORD PTR _DER2DA5
	fld	QWORD PTR _DER1DAZ
	fadd	QWORD PTR _DER3DAZ
	fld	QWORD PTR _DER1DBZ
	fadd	QWORD PTR _DER3DBZ
	fxch	ST(1)
	fadd	QWORD PTR _DER+24
	fxch	ST(2)
	fstp	QWORD PTR _DEERDA5
	fadd	QWORD PTR _DER+32
	fxch	ST(1)
	fadd	QWORD PTR _DER2DAZ
	fld	QWORD PTR _EN9
	fadd	QWORD PTR _EN14
	fxch	ST(2)
	fadd	QWORD PTR _DER2DBZ
	fxch	ST(1)
	fstp	QWORD PTR _DEERDAZ
	fstp	QWORD PTR _DEERDBZ
	fadd	QWORD PTR _EN6
	mov	DWORD PTR $T98573[esp+1404], 830	; 0000033eH

; 826  :      1            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
; 827  :      2            RX2(L),RY2(L),RZ2(L),AA2(L),
; 828  :      3            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
; 829  :      4            B12(L),B12(K),ER4,KLM)
; 830  : 	DER4DA1=DERDA1

	mov	DWORD PTR $T98573[esp+1404], 831	; 0000033fH
	fadd	QWORD PTR _EN4
	mov	DWORD PTR $T98573[esp+1404], 832	; 00000340H
	mov	DWORD PTR $T98573[esp+1404], 833	; 00000341H
	fadd	QWORD PTR _EN11
	mov	DWORD PTR $T98573[esp+1404], 834	; 00000342H
	mov	DWORD PTR $T98573[esp+1404], 839	; 00000347H
	fadd	QWORD PTR _EN8
	mov	DWORD PTR $T98573[esp+1404], 845	; 0000034dH
	mov	DWORD PTR $T98573[esp+1404], 846	; 0000034eH
	fadd	QWORD PTR _EN16
	mov	DWORD PTR $T98573[esp+1404], 847	; 0000034fH
	mov	DWORD PTR $T98573[esp+1404], 848	; 00000350H
	fadd	QWORD PTR _EN13
	mov	DWORD PTR $T98573[esp+1404], 849	; 00000351H
	mov	DWORD PTR $T98573[esp+1404], 887	; 00000377H
	fadd	QWORD PTR _EN3
	mov	DWORD PTR $T98573[esp+1404], 907	; 0000038bH
	mov	DWORD PTR $T98573[esp+1404], 925	; 0000039dH
	fadd	QWORD PTR _EN10
	mov	eax, DWORD PTR _DER+4
	mov	DWORD PTR _DER4DA1+4, eax
	mov	eax, DWORD PTR _DER
	fadd	QWORD PTR _EN15
	mov	DWORD PTR _DER4DA1, eax

; 831  : 	DER4DA2=DERDA2

	mov	eax, DWORD PTR _DER+12
	mov	DWORD PTR _DER4DA2+4, eax
	mov	eax, DWORD PTR _DER+8
	fadd	QWORD PTR _EN7
	mov	DWORD PTR _DER4DA2, eax

; 832  : 	DER4DA5=DERDA5

	mov	eax, DWORD PTR _DER+20
	mov	DWORD PTR _DER4DA5+4, eax
	mov	eax, DWORD PTR _DER+16
	fadd	QWORD PTR _EN1
	mov	DWORD PTR _DER4DA5, eax

; 833  : 	DER4DAZ=DERDAZ

	mov	eax, DWORD PTR _DER+28
	mov	DWORD PTR _DER4DAZ+4, eax
	mov	eax, DWORD PTR _DER+24
	fadd	QWORD PTR _EN5
	mov	DWORD PTR _DER4DAZ, eax

; 834  : 	DER4DBZ=DERDBZ

	mov	eax, DWORD PTR _DER+36
	mov	DWORD PTR _DER4DBZ+4, eax
	mov	eax, DWORD PTR _DER+32
	fadd	QWORD PTR _EN12
	mov	DWORD PTR _DER4DBZ, eax

; 835  : 	END IF
; 836  : 
; 837  : C              
; 838  : 	IF (ISYMM) THEN
; 839  : 		EER1=ER1+ER2+ER3+ER4

	mov	eax, DWORD PTR -1360+[esp+1396]
	mov	DWORD PTR _EER1+4, eax
	mov	eax, DWORD PTR -1360+[esp+1392]
	fadd	QWORD PTR _EN2
	mov	DWORD PTR _EER1, eax

; 840  : 	ELSE
; 841  : 		EER1=ER1
; 842  : 	END IF
; 843  : C
; 844  : 	IF (ISYMM) THEN
; 845  : 	DEERDA1=DER1DA1+DER2DA1+DER3DA1+DER4DA1
; 846  : 	DEERDA2=DER1DA2+DER2DA2+DER3DA2+DER4DA2
; 847  : 	DEERDA5=DER1DA5+DER2DA5+DER3DA5+DER4DA5
; 848  : 	DEERDAZ=DER1DAZ+DER2DAZ+DER3DAZ+DER4DAZ
; 849  : 	DEERDBZ=DER1DBZ+DER2DBZ+DER3DBZ+DER4DBZ
; 850  : 	ELSE
; 851  : 	DEERDA1=DER1DA1
; 852  : 	DEERDA2=DER1DA2
; 853  : 	DEERDA5=DER1DA5
; 854  : 	DEERDAZ=DER1DAZ
; 855  : 	DEERDBZ=DER1DBZ
; 856  : 	END IF
; 857  : 	IF(IPERM.EQ.1) THEN
; 858  : 	DEERDA1=DER1DA1
; 859  : 	DEERDA2=DER1DA2
; 860  : 	DEERDA5=DER1DA5
; 861  : 	DEERDAZ=DER1DAZ
; 862  : 	DEERDBZ=DER1DBZ
; 863  : 	END IF
; 864  : 	IF(IPERM.EQ.2) THEN
; 865  : 	DEERDA1=DER2DA1
; 866  : 	DEERDA2=DER2DA2
; 867  : 	DEERDA5=DER2DA5
; 868  : 	DEERDAZ=DER2DAZ
; 869  : 	DEERDBZ=DER2DBZ
; 870  : 	END IF
; 871  : 	IF(IPERM.EQ.3) THEN
; 872  : 	DEERDA1=DER3DA1
; 873  : 	DEERDA2=DER3DA2
; 874  : 	DEERDA5=DER3DA5
; 875  : 	DEERDAZ=DER3DAZ
; 876  : 	DEERDBZ=DER3DBZ
; 877  : 	END IF
; 878  : 	IF(IPERM.EQ.4) THEN
; 879  : 	DEERDA1=DER4DA1
; 880  : 	DEERDA2=DER4DA2
; 881  : 	DEERDA5=DER4DA5
; 882  : 	DEERDAZ=DER4DAZ
; 883  : 	DEERDBZ=DER4DBZ
; 884  : 	END IF
; 885  : C
; 886  : 	IF (ISYMM) THEN
; 887  :               EN=EN1+EN2+EN3+EN4+EN5+EN6+EN7+EN8

	mov	eax, DWORD PTR -1360+[esp+1396]
	mov	DWORD PTR _ER+4, eax
	mov	eax, DWORD PTR -1360+[esp+1392]
	mov	DWORD PTR _ER, eax
	mov	eax, DWORD PTR _M$[esp+1388]
	fstp	QWORD PTR _EN

; 888  :      1            +EN9+EN10+EN11+EN12+EN13+EN14+EN15+EN16
; 889  : 	ELSE
; 890  :               EN=EN1+EN5+EN9+EN13
; 891  : 	END IF
; 892  : C
; 893  : 	IF(IPERM.EQ.1) THEN
; 894  :               EN=EN1+EN5+EN9+EN13
; 895  : 	END IF
; 896  : 	IF(IPERM.EQ.2) THEN
; 897  :               EN=EN2+EN6+EN10+EN14
; 898  : 	END IF
; 899  : 	IF(IPERM.EQ.3) THEN
; 900  :               EN=EN3+EN7+EN11+EN15
; 901  : 	END IF
; 902  : 	IF(IPERM.EQ.4) THEN
; 903  :               EN=EN4+EN8+EN12+EN16
; 904  : 	END IF
; 905  : C
; 906  : 	IF (ISYMM) THEN
; 907  :               ER=ER1+ER2+ER3+ER4

	mov	eax, DWORD PTR [eax]
	mov	DWORD PTR -1376+[esp+1392], eax

; 908  : 	ELSE
; 909  : 			ER=ER1
; 910  : 	END IF
; 911  : C
; 912  : 	IF(IPERM.EQ.1) THEN
; 913  : 			ER=ER1
; 914  : 	END IF
; 915  : 	IF(IPERM.EQ.2) THEN
; 916  : 			ER=ER2
; 917  : 	END IF
; 918  : 	IF(IPERM.EQ.3) THEN
; 919  : 			ER=ER3
; 920  : 	END IF
; 921  : 	IF(IPERM.EQ.4) THEN
; 922  : 			ER=ER4
; 923  : 	END IF
; 924  : C
; 925  :               H(L,K)=ER+T-EN+(1.0D+00/1.4011D+00)*S(L,K)

	mov	edi, DWORD PTR _OVER
	cmp	edi, DWORD PTR -1376+[esp+1392]
	jg	SHORT $L99760
	cmp	edi, 1
	jge	SHORT $L99761
$L99760:
	call	___FHranger@0
$L99761:
	mov	esi, DWORD PTR _OVER+4
	cmp	esi, DWORD PTR -1376+[esp+1392]
	jg	SHORT $L99762
	cmp	esi, 1
	jge	SHORT $L99763
$L99762:
	call	___FHranger@0
$L99763:
	mov	ebx, DWORD PTR _OVER+4
	cmp	ebx, DWORD PTR -1376+[esp+1392]
	jg	SHORT $L99764
	cmp	ebx, 1
	jge	SHORT $L99765
$L99764:
	call	___FHranger@0
$L99765:
	imul	ebx, DWORD PTR _7S$[esp+1428]
	mov	ebp, DWORD PTR _OVER
	cmp	ebp, DWORD PTR -1376+[esp+1392]
	jg	SHORT $L99766
	cmp	ebp, 1
	jge	SHORT $L99767
$L99766:
	call	___FHranger@0
$L99767:
	add	ebx, ebp
	mov	eax, DWORD PTR _S$[esp+1388]
	add	ebx, DWORD PTR _7S$[esp+1400]
	mov	ecx, DWORD PTR _7H$[esp+1428]
	imul	ecx, esi
	fld	QWORD PTR [eax+ebx*8]
	fmul	QWORD PTR $T99768
	add	ecx, DWORD PTR _7H$[esp+1400]
	mov	eax, DWORD PTR _H$[esp+1388]
	fsub	QWORD PTR _EN
	add	ecx, edi
	mov	esi, DWORD PTR _OVER
	mov	DWORD PTR $T98573[esp+1404], 927	; 0000039fH
	fadd	QWORD PTR _T
	push	1
	push	OFFSET FLAT:_7DHDA1
	fadd	QWORD PTR -1360+[esp+1400]
	fstp	QWORD PTR [eax+ecx*8]

; 926  : C
; 927  : 	DHDA1(L,K)=DEERDA1+DTDA1-DENADA1+(1.0D+00/1.4011D+00)*DSDA1(L,K)

	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L99769
	push	1
	push	OFFSET FLAT:_7DHDA1
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L99770
$L99769:
	call	___FHranger@0
$L99770:
	mov	edi, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DHDA1
	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L99771
	push	2
	push	OFFSET FLAT:_7DHDA1
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L99772
$L99771:
	call	___FHranger@0
$L99772:
	mov	ebx, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DSDA1
	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L99773
	push	2
	push	OFFSET FLAT:_7DSDA1
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L99774
$L99773:
	call	___FHranger@0
$L99774:
	mov	ebp, DWORD PTR _OVER
	push	1
	imul	ebx, DWORD PTR _7DSDA1+36
	push	OFFSET FLAT:_7DSDA1
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L99775
	push	1
	push	OFFSET FLAT:_7DSDA1
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L99776
$L99775:
	call	___FHranger@0
$L99776:
	add	ebx, ebp
	mov	eax, DWORD PTR _DSDA1
	add	ebx, DWORD PTR _7DSDA1+8
	mov	ecx, DWORD PTR _DHDA1
	mov	DWORD PTR $T98573[esp+1404], 928	; 000003a0H
	fld	QWORD PTR [eax+ebx*8]
	fmul	QWORD PTR $T99768
	mov	eax, DWORD PTR _7DHDA1+36
	push	1
	push	OFFSET FLAT:_7DHDA2
	fsub	QWORD PTR _DENADA1
	imul	eax, edi
	fadd	QWORD PTR _DEERDA1
	add	eax, DWORD PTR _7DHDA1+8
	fadd	QWORD PTR _DTDA1
	add	eax, esi
	mov	esi, DWORD PTR _OVER
	fstp	QWORD PTR [ecx+eax*8]

; 928  : 	DHDA2(L,K)=DEERDA2+DTDA2-DENADA2+(1.0D+00/1.4011D+00)*DSDA2(L,K)

	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L99777
	push	1
	push	OFFSET FLAT:_7DHDA2
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L99778
$L99777:
	call	___FHranger@0
$L99778:
	mov	edi, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DHDA2
	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L99779
	push	2
	push	OFFSET FLAT:_7DHDA2
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L99780
$L99779:
	call	___FHranger@0
$L99780:
	mov	ebx, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DSDA2
	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L99781
	push	2
	push	OFFSET FLAT:_7DSDA2
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L99782
$L99781:
	call	___FHranger@0
$L99782:
	mov	ebp, DWORD PTR _OVER
	push	1
	imul	ebx, DWORD PTR _7DSDA2+36
	push	OFFSET FLAT:_7DSDA2
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L99783
	push	1
	push	OFFSET FLAT:_7DSDA2
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L99784
$L99783:
	call	___FHranger@0
$L99784:
	add	ebx, ebp
	mov	eax, DWORD PTR _DSDA2
	add	ebx, DWORD PTR _7DSDA2+8
	mov	ecx, DWORD PTR _7DHDA2+36
	imul	ecx, edi
	fld	QWORD PTR [eax+ebx*8]
	fmul	QWORD PTR $T99768
	add	ecx, DWORD PTR _7DHDA2+8
	mov	eax, DWORD PTR _DHDA2
	fsub	QWORD PTR _DENADA2
	add	ecx, esi
	mov	edi, DWORD PTR _OVER
	mov	DWORD PTR $T98573[esp+1404], 929	; 000003a1H
	fadd	QWORD PTR _DTDA2
	push	1
	push	OFFSET FLAT:_7DHDA5
	fadd	QWORD PTR _DEERDA2
	fstp	QWORD PTR [eax+ecx*8]

; 929  : 	DHDA5(L,K)=DEERDA5+DTDA5-DENADA5+(1.0D+00/1.4011D+00)*DSDA5(L,K)

	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L99785
	push	1
	push	OFFSET FLAT:_7DHDA5
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L99786
$L99785:
	call	___FHranger@0
$L99786:
	mov	esi, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DHDA5
	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L99787
	push	2
	push	OFFSET FLAT:_7DHDA5
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L99788
$L99787:
	call	___FHranger@0
$L99788:
	mov	ebx, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DSDA5
	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L99789
	push	2
	push	OFFSET FLAT:_7DSDA5
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L99790
$L99789:
	call	___FHranger@0
$L99790:
	mov	ebp, DWORD PTR _OVER
	push	1
	imul	ebx, DWORD PTR _7DSDA5+36
	push	OFFSET FLAT:_7DSDA5
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L99791
	push	1
	push	OFFSET FLAT:_7DSDA5
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L99792
$L99791:
	call	___FHranger@0
$L99792:
	add	ebx, ebp
	mov	eax, DWORD PTR _DSDA5
	add	ebx, DWORD PTR _7DSDA5+8
	mov	ecx, DWORD PTR _7DHDA5+36
	imul	ecx, esi
	fld	QWORD PTR [eax+ebx*8]
	fmul	QWORD PTR $T99768
	add	ecx, DWORD PTR _7DHDA5+8
	mov	eax, DWORD PTR _DHDA5
	fsub	QWORD PTR _DENADA5
	add	ecx, edi
	mov	esi, DWORD PTR _OVER
	mov	DWORD PTR $T98573[esp+1404], 930	; 000003a2H
	fadd	QWORD PTR _DTDA5
	push	1
	push	OFFSET FLAT:_7DHDAZ
	fadd	QWORD PTR _DEERDA5
	fstp	QWORD PTR [eax+ecx*8]

; 930  : 	DHDAZ(L,K)=DEERDAZ+DTDAZ-DENADAZ+(1.0D+00/1.4011D+00)*DSDAZ(L,K)

	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L99793
	push	1
	push	OFFSET FLAT:_7DHDAZ
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L99794
$L99793:
	call	___FHranger@0
$L99794:
	mov	edi, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DHDAZ
	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L99795
	push	2
	push	OFFSET FLAT:_7DHDAZ
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L99796
$L99795:
	call	___FHranger@0
$L99796:
	mov	ebx, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DSDAZ
	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L99797
	push	2
	push	OFFSET FLAT:_7DSDAZ
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L99798
$L99797:
	call	___FHranger@0
$L99798:
	mov	ebp, DWORD PTR _OVER
	push	1
	imul	ebx, DWORD PTR _7DSDAZ+36
	push	OFFSET FLAT:_7DSDAZ
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L99799
	push	1
	push	OFFSET FLAT:_7DSDAZ
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L99800
$L99799:
	call	___FHranger@0
$L99800:
	add	ebx, ebp
	mov	eax, DWORD PTR _DSDAZ
	add	ebx, DWORD PTR _7DSDAZ+8
	push	1
	imul	edi, DWORD PTR _7DHDAZ+36
	fld	QWORD PTR [eax+ebx*8]
	fmul	QWORD PTR $T99768
	add	edi, esi
	mov	eax, DWORD PTR _DHDAZ
	add	edi, DWORD PTR _7DHDAZ+8
	mov	esi, DWORD PTR _OVER
	fsub	QWORD PTR _DENADAZ
	mov	DWORD PTR $T98573[esp+1408], 931	; 000003a3H
	push	OFFSET FLAT:_7DHDBZ
	fadd	QWORD PTR _DTDAZ
	fadd	QWORD PTR _DEERDAZ
	fstp	QWORD PTR [eax+edi*8]

; 931  : 	DHDBZ(L,K)=DEERDBZ+DTDBZ-DENADBZ+(1.0D+00/1.4011D+00)*DSDBZ(L,K)

	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L99801
	push	1
	push	OFFSET FLAT:_7DHDBZ
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L99802
$L99801:
	call	___FHranger@0
$L99802:
	mov	edi, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DHDBZ
	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L99803
	push	2
	push	OFFSET FLAT:_7DHDBZ
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L99804
$L99803:
	call	___FHranger@0
$L99804:
	mov	ebx, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DSDBZ
	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L99805
	push	2
	push	OFFSET FLAT:_7DSDBZ
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L99806
$L99805:
	call	___FHranger@0
$L99806:
	mov	ebp, DWORD PTR _OVER
	push	1
	imul	ebx, DWORD PTR _7DSDBZ+36
	push	OFFSET FLAT:_7DSDBZ
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L99807
	push	1
	push	OFFSET FLAT:_7DSDBZ
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L99808
$L99807:
	call	___FHranger@0
$L99808:
	add	ebx, ebp
	mov	eax, DWORD PTR _DSDBZ
	mov	DWORD PTR $T98573[esp+1404], 933	; 000003a5H
	imul	edi, DWORD PTR _7DHDBZ+36
	add	ebx, DWORD PTR _7DSDBZ+8
	add	edi, esi
	add	edi, DWORD PTR _7DHDBZ+8
	fld	QWORD PTR [eax+ebx*8]
	fmul	QWORD PTR $T99768
	mov	eax, DWORD PTR _DHDBZ
	fsub	QWORD PTR _DENADBZ
	fadd	QWORD PTR _DTDBZ
	fadd	QWORD PTR _DEERDBZ
	fstp	QWORD PTR [eax+edi*8]

; 932  : C
; 933  : 	IF (L.EQ.K) THEN

	mov	eax, DWORD PTR _OVER
	cmp	DWORD PTR _OVER+4, eax
	jne	$L156
	mov	esi, DWORD PTR _OVER
	push	1
	mov	DWORD PTR $T98573[esp+1408], 934	; 000003a6H

; 934  : 		DHDA1(L,K)=2.00D+00*DHDA1(L,K)

	push	OFFSET FLAT:_7DHDA1
	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L99809
	push	1
	push	OFFSET FLAT:_7DHDA1
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L99810
$L99809:
	call	___FHranger@0
$L99810:
	mov	edi, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DHDA1
	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L99811
	push	2
	push	OFFSET FLAT:_7DHDA1
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L99812
$L99811:
	call	___FHranger@0
$L99812:
	mov	ebx, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DHDA1
	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L99813
	push	2
	push	OFFSET FLAT:_7DHDA1
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L99814
$L99813:
	call	___FHranger@0
$L99814:
	mov	ebp, DWORD PTR _OVER
	push	1
	imul	ebx, DWORD PTR _7DHDA1+36
	push	OFFSET FLAT:_7DHDA1
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L99815
	push	1
	push	OFFSET FLAT:_7DHDA1
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L99816
$L99815:
	call	___FHranger@0
$L99816:
	add	ebx, ebp
	mov	eax, DWORD PTR _7DHDA1+8
	add	ebx, eax
	mov	ecx, DWORD PTR _DHDA1
	mov	edx, DWORD PTR _DHDA1
	push	1
	fld	QWORD PTR [ecx+ebx*8]
	fmul	QWORD PTR $T99039
	mov	ecx, DWORD PTR _7DHDA1+36
	push	OFFSET FLAT:_7DHDA2
	mov	DWORD PTR $T98573[esp+1412], 935	; 000003a7H
	imul	ecx, edi
	add	ecx, eax
	mov	edi, DWORD PTR _OVER
	add	ecx, esi
	fstp	QWORD PTR [edx+ecx*8]

; 935  : 		DHDA2(L,K)=2.00D+00*DHDA2(L,K)

	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L99817
	push	1
	push	OFFSET FLAT:_7DHDA2
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L99818
$L99817:
	call	___FHranger@0
$L99818:
	mov	esi, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DHDA2
	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L99819
	push	2
	push	OFFSET FLAT:_7DHDA2
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L99820
$L99819:
	call	___FHranger@0
$L99820:
	mov	ebx, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DHDA2
	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L99821
	push	2
	push	OFFSET FLAT:_7DHDA2
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L99822
$L99821:
	call	___FHranger@0
$L99822:
	mov	ebp, DWORD PTR _OVER
	push	1
	imul	ebx, DWORD PTR _7DHDA2+36
	push	OFFSET FLAT:_7DHDA2
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L99823
	push	1
	push	OFFSET FLAT:_7DHDA2
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L99824
$L99823:
	call	___FHranger@0
$L99824:
	add	ebx, ebp
	mov	eax, DWORD PTR _7DHDA2+8
	add	ebx, eax
	mov	ecx, DWORD PTR _DHDA2
	mov	edx, DWORD PTR _DHDA2
	push	1
	fld	QWORD PTR [ecx+ebx*8]
	fmul	QWORD PTR $T99039
	mov	ecx, DWORD PTR _7DHDA2+36
	push	OFFSET FLAT:_7DHDA5
	mov	DWORD PTR $T98573[esp+1412], 936	; 000003a8H
	imul	ecx, esi
	add	ecx, eax
	add	ecx, edi
	mov	edi, DWORD PTR _OVER
	fstp	QWORD PTR [edx+ecx*8]

; 936  : 		DHDA5(L,K)=2.00D+00*DHDA5(L,K)

	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L99825
	push	1
	push	OFFSET FLAT:_7DHDA5
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L99826
$L99825:
	call	___FHranger@0
$L99826:
	mov	esi, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DHDA5
	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L99827
	push	2
	push	OFFSET FLAT:_7DHDA5
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L99828
$L99827:
	call	___FHranger@0
$L99828:
	mov	ebx, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DHDA5
	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L99829
	push	2
	push	OFFSET FLAT:_7DHDA5
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L99830
$L99829:
	call	___FHranger@0
$L99830:
	mov	ebp, DWORD PTR _OVER
	push	1
	imul	ebx, DWORD PTR _7DHDA5+36
	push	OFFSET FLAT:_7DHDA5
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L99831
	push	1
	push	OFFSET FLAT:_7DHDA5
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L99832
$L99831:
	call	___FHranger@0
$L99832:
	add	ebx, ebp
	mov	eax, DWORD PTR _7DHDA5+8
	add	ebx, eax
	mov	ecx, DWORD PTR _DHDA5
	mov	eax, DWORD PTR _7DHDA5+36
	mov	edx, DWORD PTR _7DHDA5+8
	fld	QWORD PTR [ecx+ebx*8]
	fmul	QWORD PTR $T99039
	imul	eax, esi
	add	eax, edx
	mov	esi, DWORD PTR _OVER
	add	eax, edi
	push	1
	mov	DWORD PTR $T98573[esp+1408], 937	; 000003a9H
	fstp	QWORD PTR [ecx+eax*8]

; 937  : 		DHDAZ(L,K)=2.00D+00*DHDAZ(L,K)

	push	OFFSET FLAT:_7DHDAZ
	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L99833
	push	1
	push	OFFSET FLAT:_7DHDAZ
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L99834
$L99833:
	call	___FHranger@0
$L99834:
	mov	edi, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DHDAZ
	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L99835
	push	2
	push	OFFSET FLAT:_7DHDAZ
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L99836
$L99835:
	call	___FHranger@0
$L99836:
	mov	ebx, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DHDAZ
	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L99837
	push	2
	push	OFFSET FLAT:_7DHDAZ
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L99838
$L99837:
	call	___FHranger@0
$L99838:
	mov	ebp, DWORD PTR _OVER
	push	1
	imul	ebx, DWORD PTR _7DHDAZ+36
	push	OFFSET FLAT:_7DHDAZ
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L99839
	push	1
	push	OFFSET FLAT:_7DHDAZ
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L99840
$L99839:
	call	___FHranger@0
$L99840:
	add	ebx, ebp
	mov	eax, DWORD PTR _7DHDAZ+8
	add	ebx, eax
	mov	ecx, DWORD PTR _DHDAZ
	mov	DWORD PTR $T98573[esp+1404], 938	; 000003aaH
	fld	QWORD PTR [ecx+ebx*8]
	fmul	QWORD PTR $T99039
	imul	edi, DWORD PTR _7DHDAZ+36
	add	edi, esi
	push	1
	add	edi, eax
	mov	esi, DWORD PTR _OVER
	push	OFFSET FLAT:_7DHDBZ
	fstp	QWORD PTR [ecx+edi*8]

; 938  : 		DHDBZ(L,K)=2.00D+00*DHDBZ(L,K)

	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L99841
	push	1
	push	OFFSET FLAT:_7DHDBZ
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L99842
$L99841:
	call	___FHranger@0
$L99842:
	mov	edi, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DHDBZ
	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L99843
	push	2
	push	OFFSET FLAT:_7DHDBZ
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L99844
$L99843:
	call	___FHranger@0
$L99844:
	mov	ebx, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DHDBZ
	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L99845
	push	2
	push	OFFSET FLAT:_7DHDBZ
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L99846
$L99845:
	call	___FHranger@0
$L99846:
	mov	ebp, DWORD PTR _OVER
	push	1
	imul	ebx, DWORD PTR _7DHDBZ+36
	push	OFFSET FLAT:_7DHDBZ
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L99847
	push	1
	push	OFFSET FLAT:_7DHDBZ
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L99848
$L99847:
	call	___FHranger@0
$L99848:
	add	ebx, ebp
	mov	eax, DWORD PTR _7DHDBZ+8
	add	ebx, eax
	mov	ecx, DWORD PTR _DHDBZ
	mov	DWORD PTR $T98573[esp+1404], 940	; 000003acH
	fld	QWORD PTR [ecx+ebx*8]
	fmul	QWORD PTR $T99039
	imul	edi, DWORD PTR _7DHDBZ+36
	add	edi, esi
	push	1
	add	edi, eax
	push	OFFSET FLAT:_7DSDA1
	fstp	QWORD PTR [ecx+edi*8]

; 939  : C
; 940  : 		DSDA1(L,K)=2.00D+00*DSDA1(L,K)

	mov	edi, DWORD PTR _OVER
	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L99849
	push	1
	push	OFFSET FLAT:_7DSDA1
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L99850
$L99849:
	call	___FHranger@0
$L99850:
	mov	esi, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DSDA1
	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L99851
	push	2
	push	OFFSET FLAT:_7DSDA1
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L99852
$L99851:
	call	___FHranger@0
$L99852:
	mov	ebx, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DSDA1
	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L99853
	push	2
	push	OFFSET FLAT:_7DSDA1
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L99854
$L99853:
	call	___FHranger@0
$L99854:
	mov	ebp, DWORD PTR _OVER
	push	1
	imul	ebx, DWORD PTR _7DSDA1+36
	push	OFFSET FLAT:_7DSDA1
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L99855
	push	1
	push	OFFSET FLAT:_7DSDA1
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L99856
$L99855:
	call	___FHranger@0
$L99856:
	add	ebx, ebp
	mov	eax, DWORD PTR _7DSDA1+8
	add	ebx, eax
	mov	ecx, DWORD PTR _DSDA1
	mov	DWORD PTR $T98573[esp+1404], 941	; 000003adH
	fld	QWORD PTR [ecx+ebx*8]
	fmul	QWORD PTR $T99039
	imul	esi, DWORD PTR _7DSDA1+36
	add	esi, edi
	push	1
	add	esi, eax
	push	OFFSET FLAT:_7DSDA2
	fstp	QWORD PTR [ecx+esi*8]

; 941  : 		DSDA2(L,K)=2.00D+00*DSDA2(L,K)

	mov	esi, DWORD PTR _OVER
	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L99857
	push	1
	push	OFFSET FLAT:_7DSDA2
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L99858
$L99857:
	call	___FHranger@0
$L99858:
	mov	edi, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DSDA2
	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L99859
	push	2
	push	OFFSET FLAT:_7DSDA2
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L99860
$L99859:
	call	___FHranger@0
$L99860:
	mov	ebx, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DSDA2
	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L99861
	push	2
	push	OFFSET FLAT:_7DSDA2
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L99862
$L99861:
	call	___FHranger@0
$L99862:
	mov	ebp, DWORD PTR _OVER
	push	1
	imul	ebx, DWORD PTR _7DSDA2+36
	push	OFFSET FLAT:_7DSDA2
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L99863
	push	1
	push	OFFSET FLAT:_7DSDA2
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L99864
$L99863:
	call	___FHranger@0
$L99864:
	add	ebx, ebp
	mov	eax, DWORD PTR _7DSDA2+8
	add	ebx, eax
	mov	ecx, DWORD PTR _DSDA2
	mov	DWORD PTR $T98573[esp+1404], 942	; 000003aeH
	fld	QWORD PTR [ecx+ebx*8]
	fmul	QWORD PTR $T99039
	imul	edi, DWORD PTR _7DSDA2+36
	add	edi, esi
	push	1
	add	edi, eax
	push	OFFSET FLAT:_7DSDA5
	fstp	QWORD PTR [ecx+edi*8]

; 942  : 		DSDA5(L,K)=2.00D+00*DSDA5(L,K)

	mov	edi, DWORD PTR _OVER
	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L99865
	push	1
	push	OFFSET FLAT:_7DSDA5
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L99866
$L99865:
	call	___FHranger@0
$L99866:
	mov	esi, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DSDA5
	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L99867
	push	2
	push	OFFSET FLAT:_7DSDA5
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L99868
$L99867:
	call	___FHranger@0
$L99868:
	mov	ebx, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DSDA5
	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L99869
	push	2
	push	OFFSET FLAT:_7DSDA5
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L99870
$L99869:
	call	___FHranger@0
$L99870:
	mov	ebp, DWORD PTR _OVER
	push	1
	imul	ebx, DWORD PTR _7DSDA5+36
	push	OFFSET FLAT:_7DSDA5
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L99871
	push	1
	push	OFFSET FLAT:_7DSDA5
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L99872
$L99871:
	call	___FHranger@0
$L99872:
	add	ebx, ebp
	mov	eax, DWORD PTR _7DSDA5+8
	add	ebx, eax
	mov	ecx, DWORD PTR _DSDA5
	mov	edx, DWORD PTR _DSDA5
	push	1
	fld	QWORD PTR [ecx+ebx*8]
	fmul	QWORD PTR $T99039
	mov	ecx, DWORD PTR _7DSDA5+36
	push	OFFSET FLAT:_7DSDAZ
	mov	DWORD PTR $T98573[esp+1412], 943	; 000003afH
	imul	ecx, esi
	add	ecx, eax
	add	ecx, edi
	mov	edi, DWORD PTR _OVER
	fstp	QWORD PTR [edx+ecx*8]

; 943  : 		DSDAZ(L,K)=2.00D+00*DSDAZ(L,K)

	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L99873
	push	1
	push	OFFSET FLAT:_7DSDAZ
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L99874
$L99873:
	call	___FHranger@0
$L99874:
	mov	esi, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DSDAZ
	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L99875
	push	2
	push	OFFSET FLAT:_7DSDAZ
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L99876
$L99875:
	call	___FHranger@0
$L99876:
	mov	ebx, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DSDAZ
	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L99877
	push	2
	push	OFFSET FLAT:_7DSDAZ
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L99878
$L99877:
	call	___FHranger@0
$L99878:
	mov	ebp, DWORD PTR _OVER
	push	1
	imul	ebx, DWORD PTR _7DSDAZ+36
	push	OFFSET FLAT:_7DSDAZ
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L99879
	push	1
	push	OFFSET FLAT:_7DSDAZ
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L99880
$L99879:
	call	___FHranger@0
$L99880:
	add	ebx, ebp
	mov	eax, DWORD PTR _7DSDAZ+8
	add	ebx, eax
	mov	ecx, DWORD PTR _DSDAZ
	mov	edx, DWORD PTR _DSDAZ
	push	1
	fld	QWORD PTR [ecx+ebx*8]
	fmul	QWORD PTR $T99039
	mov	ecx, DWORD PTR _7DSDAZ+36
	push	OFFSET FLAT:_7DSDBZ
	mov	DWORD PTR $T98573[esp+1412], 944	; 000003b0H
	imul	ecx, esi
	add	ecx, eax
	add	ecx, edi
	mov	edi, DWORD PTR _OVER
	fstp	QWORD PTR [edx+ecx*8]

; 944  : 		DSDBZ(L,K)=2.00D+00*DSDBZ(L,K)

	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L99881
	push	1
	push	OFFSET FLAT:_7DSDBZ
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L99882
$L99881:
	call	___FHranger@0
$L99882:
	mov	esi, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DSDBZ
	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L99883
	push	2
	push	OFFSET FLAT:_7DSDBZ
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L99884
$L99883:
	call	___FHranger@0
$L99884:
	mov	ebx, DWORD PTR _OVER+4
	push	2
	push	OFFSET FLAT:_7DSDBZ
	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L99885
	push	2
	push	OFFSET FLAT:_7DSDBZ
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L99886
$L99885:
	call	___FHranger@0
$L99886:
	mov	ebp, DWORD PTR _OVER
	push	1
	imul	ebx, DWORD PTR _7DSDBZ+36
	push	OFFSET FLAT:_7DSDBZ
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L99887
	push	1
	push	OFFSET FLAT:_7DSDBZ
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L99888
$L99887:
	call	___FHranger@0
$L99888:
	add	ebx, ebp
	mov	eax, DWORD PTR _7DSDBZ+8
	add	ebx, eax
	mov	ecx, DWORD PTR _DSDBZ
	mov	DWORD PTR $T98573[esp+1404], 55		; 00000037H
	fld	QWORD PTR [ecx+ebx*8]
	fmul	QWORD PTR $T99039
	mov	ecx, DWORD PTR _7DSDBZ+36
	mov	edx, DWORD PTR _DSDBZ
	imul	ecx, esi
	add	ecx, eax
	add	ecx, edi
	fstp	QWORD PTR [edx+ecx*8]
$L156:

; 55   :           DO K=1,M

	mov	eax, DWORD PTR -1372+[esp+1392]
	mov	ecx, DWORD PTR _OVER+4
	dec	eax
	inc	ecx
	mov	DWORD PTR -1372+[esp+1392], eax
	mov	DWORD PTR _OVER+4, ecx
	mov	DWORD PTR $T98573[esp+1404], 55		; 00000037H
	test	eax, eax
	jne	$L98307
	mov	DWORD PTR $T98573[esp+1404], 54		; 00000036H
$L153:

; 54   :       DO L=1,M

	mov	eax, DWORD PTR -1364+[esp+1392]
	mov	ecx, DWORD PTR _OVER
	dec	eax
	inc	ecx
	mov	DWORD PTR -1364+[esp+1392], eax
	mov	DWORD PTR _OVER, ecx
	mov	DWORD PTR $T98573[esp+1404], 54		; 00000036H
	test	eax, eax
	jne	$L98306
$L152:

; 945  : 	END IF
; 946  : C     ----------------------------------------------------- 
; 947  :           END DO 
; 948  :       END DO

	mov	DWORD PTR $T98573[esp+1404], 955	; 000003bbH

; 949  : C     ----------------------------------------------------- 
; 950  : C      CALL TAB(S,M,M,MMAX,MMAX)                  
; 951  : C      CALL TAB(H,M,M,MMAX,MMAX)
; 952  : C      CALL ENERGY(M,H,S,EMIN,CO,KIM)
; 953  : C		STOP
; 954  : C     -----------------------------------------------------      
; 955  :       SUMH=0.00D+00

	xor	ecx, ecx
	mov	edx, DWORD PTR _M$[esp+1388]
	mov	eax, 957				; 000003bdH
	mov	DWORD PTR _SUMH, ecx
	mov	ebp, DWORD PTR [edx]
	mov	DWORD PTR _SUMH+4, ecx
	mov	edi, 1
	mov	DWORD PTR _SUMS, ecx
	mov	DWORD PTR _SUMS+4, ecx
	mov	DWORD PTR -592+[esp+1392], ebp
	mov	DWORD PTR $T98573[esp+1404], 956	; 000003bcH

; 956  : 	SUMS=0.00D+00

	mov	DWORD PTR $T98573[esp+1404], eax
	cmp	ebp, ecx

; 957  : 	DO J=1,M

	jle	$L521
$L98404:
	mov	ecx, 958				; 000003beH
	mov	esi, 1
	mov	DWORD PTR $T98573[esp+1404], ecx
	mov	DWORD PTR -600+[esp+1392], ebp

; 958  : 		DO I=1,M

	test	ebp, ebp
	jle	$L522
$L98405:
	cmp	edi, ebp
	mov	ebx, edi
	mov	DWORD PTR $T98573[esp+1404], 959	; 000003bfH

; 959  : 			SUMH=SUMH+EVEC(I)*H(I,J)*EVEC(J)

	jg	SHORT $L99890
	cmp	ebx, 1
	jge	SHORT $L99891
$L99890:
	call	___FHranger@0
$L99891:
	imul	ebx, DWORD PTR _7H$[esp+1428]
	cmp	esi, ebp
	jg	SHORT $L99892
	cmp	esi, 1
	jge	SHORT $L99893
$L99892:
	call	___FHranger@0
$L99893:
	add	ebx, esi
	mov	eax, DWORD PTR _H$[esp+1388]
	add	ebx, DWORD PTR _7H$[esp+1400]
	cmp	edi, ebp
	fld	QWORD PTR [eax+ebx*8]
	jg	SHORT $L99894
	cmp	edi, 1
	jge	SHORT $L99895
$L99894:
	call	___FHranger@0
$L99895:
	mov	eax, DWORD PTR _EVEC$[esp+1388]
	cmp	esi, ebp
	fmul	QWORD PTR [eax+edi*8-8]
	jg	SHORT $L99896
	cmp	esi, 1
	jge	SHORT $L99897
$L99896:
	call	___FHranger@0
$L99897:
	cmp	edi, ebp
	mov	eax, DWORD PTR _EVEC$[esp+1388]
	mov	DWORD PTR $T98573[esp+1404], 960	; 000003c0H
	fmul	QWORD PTR [eax+esi*8-8]
	fadd	QWORD PTR _SUMH
	fstp	QWORD PTR _SUMH

; 960  : 			SUMS=SUMS+EVEC(I)*S(I,J)*EVEC(J)

	mov	ebx, edi
	jg	SHORT $L99898
	cmp	ebx, 1
	jge	SHORT $L99899
$L99898:
	call	___FHranger@0
$L99899:
	imul	ebx, DWORD PTR _7S$[esp+1428]
	cmp	esi, ebp
	jg	SHORT $L99900
	cmp	esi, 1
	jge	SHORT $L99901
$L99900:
	call	___FHranger@0
$L99901:
	add	ebx, esi
	mov	eax, DWORD PTR _S$[esp+1388]
	add	ebx, DWORD PTR _7S$[esp+1400]
	cmp	edi, ebp
	fld	QWORD PTR [eax+ebx*8]
	jg	SHORT $L99902
	cmp	edi, 1
	jge	SHORT $L99903
$L99902:
	call	___FHranger@0
$L99903:
	mov	eax, DWORD PTR _EVEC$[esp+1388]
	cmp	esi, ebp
	fmul	QWORD PTR [eax+edi*8-8]
	jg	SHORT $L99904
	cmp	esi, 1
	jge	SHORT $L99905
$L99904:
	call	___FHranger@0
$L99905:
	mov	eax, DWORD PTR _EVEC$[esp+1388]
	inc	esi
	mov	DWORD PTR $T98573[esp+1404], 958	; 000003beH
	fmul	QWORD PTR [eax+esi*8-16]
	fadd	QWORD PTR _SUMS
	fstp	QWORD PTR _SUMS

; 958  : 		DO I=1,M

	dec	DWORD PTR -600+[esp+1392]
	cmp	DWORD PTR -600+[esp+1392], 0
	jne	$L98405
	mov	DWORD PTR $T98573[esp+1404], 957	; 000003bdH
$L522:

; 957  : 	DO J=1,M

	mov	eax, DWORD PTR -592+[esp+1392]
	inc	edi
	mov	DWORD PTR $T98573[esp+1404], 957	; 000003bdH
	dec	eax
	mov	DWORD PTR -592+[esp+1392], eax
	test	eax, eax
	jne	$L98404
$L521:

; 961  : 		END DO
; 962  : 	END DO

	fld	QWORD PTR _SUMH
	fdiv	QWORD PTR _SUMS
	mov	ecx, DWORD PTR _EMIN$[esp+1388]
	mov	DWORD PTR $T98573[esp+1404], 964	; 000003c4H

; 963  : C
; 964  : 	ENG=SUMH/SUMS

	mov	DWORD PTR $T98573[esp+1404], 965	; 000003c5H
	fstp	QWORD PTR _ENG

; 965  : 	EMIN=ENG

	mov	edx, DWORD PTR _ENG+4
	mov	eax, DWORD PTR _ENG
	mov	DWORD PTR [ecx+4], edx
	mov	DWORD PTR -592+[esp+1392], ebp
	mov	DWORD PTR [ecx], eax
	test	ebp, ebp
	mov	DWORD PTR _I$[esp+1392], 1
	mov	eax, 967				; 000003c7H
	mov	DWORD PTR $T98573[esp+1404], eax

; 966  : C
; 967  : 	DO I=1,M

	jle	$L527
$L98406:
	mov	esi, 1
	mov	ecx, 969				; 000003c9H
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	DWORD PTR _OVER+4, esi
	mov	DWORD PTR $T98573[esp+1404], 968	; 000003c8H

; 968  : 		K=1

	mov	eax, DWORD PTR [eax]
	mov	DWORD PTR $T98573[esp+1404], ecx

; 969  : 		DO J=1,M

	mov	DWORD PTR -600+[esp+1392], eax
	test	eax, eax
	jle	$L528
$L98407:
	mov	ebp, DWORD PTR _OVER+4
	push	1
	mov	DWORD PTR $T98573[esp+1408], 970	; 000003caH

; 970  : 			DH(K,I)=DHDA1(J,I)

	push	OFFSET FLAT:_7DH
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L99906
	push	1
	push	OFFSET FLAT:_7DH
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L99907
$L99906:
	call	___FHranger@0
$L99907:
	mov	edi, DWORD PTR _I$[esp+1392]
	push	2
	push	OFFSET FLAT:_7DH
	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L99908
	push	2
	push	OFFSET FLAT:_7DH
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L99909
$L99908:
	call	___FHranger@0
$L99909:
	mov	ebx, DWORD PTR _I$[esp+1392]
	push	2
	push	OFFSET FLAT:_7DHDA1
	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L99910
	push	2
	push	OFFSET FLAT:_7DHDA1
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L99911
$L99910:
	call	___FHranger@0
$L99911:
	imul	ebx, DWORD PTR _7DHDA1+36
	push	1
	push	OFFSET FLAT:_7DHDA1
	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L99912
	push	1
	push	OFFSET FLAT:_7DHDA1
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L99913
$L99912:
	call	___FHranger@0
$L99913:
	add	ebx, esi
	mov	edx, DWORD PTR _DHDA1
	add	ebx, DWORD PTR _7DHDA1+8
	mov	ecx, DWORD PTR _7DH+36
	imul	ecx, edi
	add	ecx, DWORD PTR _7DH+8
	mov	eax, DWORD PTR [edx+ebx*8+4]
	add	ecx, ebp
	mov	edi, DWORD PTR _DH
	mov	DWORD PTR $T98573[esp+1404], 971	; 000003cbH
	mov	ebp, DWORD PTR _OVER+4
	mov	DWORD PTR [edi+ecx*8+4], eax
	mov	eax, DWORD PTR [edx+ebx*8]
	inc	ebp
	mov	DWORD PTR [edi+ecx*8], eax

; 971  : 			DH(K+1,I)=DHDA2(J,I)

	jno	SHORT $L99914
	call	___FHovfl@0
$L99914:
	push	1
	push	OFFSET FLAT:_7DH
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L99915
	push	1
	push	OFFSET FLAT:_7DH
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L99916
$L99915:
	call	___FHranger@0
$L99916:
	mov	edi, DWORD PTR _I$[esp+1392]
	push	2
	push	OFFSET FLAT:_7DH
	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L99917
	push	2
	push	OFFSET FLAT:_7DH
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L99918
$L99917:
	call	___FHranger@0
$L99918:
	mov	ebx, DWORD PTR _I$[esp+1392]
	push	2
	push	OFFSET FLAT:_7DHDA2
	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L99919
	push	2
	push	OFFSET FLAT:_7DHDA2
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L99920
$L99919:
	call	___FHranger@0
$L99920:
	imul	ebx, DWORD PTR _7DHDA2+36
	push	1
	push	OFFSET FLAT:_7DHDA2
	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L99921
	push	1
	push	OFFSET FLAT:_7DHDA2
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L99922
$L99921:
	call	___FHranger@0
$L99922:
	add	ebx, esi
	mov	edx, DWORD PTR _DHDA2
	add	ebx, DWORD PTR _7DHDA2+8
	mov	ecx, DWORD PTR _7DH+36
	imul	ecx, edi
	add	ecx, DWORD PTR _7DH+8
	mov	eax, DWORD PTR [edx+ebx*8+4]
	add	ecx, ebp
	mov	edi, DWORD PTR _DH
	mov	DWORD PTR $T98573[esp+1404], 972	; 000003ccH
	mov	ebp, DWORD PTR _OVER+4
	mov	DWORD PTR [edi+ecx*8+4], eax
	mov	eax, DWORD PTR [edx+ebx*8]
	add	ebp, 2
	mov	DWORD PTR [edi+ecx*8], eax

; 972  : 			DH(K+2,I)=DHDA5(J,I)

	jno	SHORT $L99923
	call	___FHovfl@0
$L99923:
	push	1
	push	OFFSET FLAT:_7DH
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L99924
	push	1
	push	OFFSET FLAT:_7DH
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L99925
$L99924:
	call	___FHranger@0
$L99925:
	mov	edi, DWORD PTR _I$[esp+1392]
	push	2
	push	OFFSET FLAT:_7DH
	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L99926
	push	2
	push	OFFSET FLAT:_7DH
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L99927
$L99926:
	call	___FHranger@0
$L99927:
	mov	ebx, DWORD PTR _I$[esp+1392]
	push	2
	push	OFFSET FLAT:_7DHDA5
	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L99928
	push	2
	push	OFFSET FLAT:_7DHDA5
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L99929
$L99928:
	call	___FHranger@0
$L99929:
	imul	ebx, DWORD PTR _7DHDA5+36
	push	1
	push	OFFSET FLAT:_7DHDA5
	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L99930
	push	1
	push	OFFSET FLAT:_7DHDA5
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L99931
$L99930:
	call	___FHranger@0
$L99931:
	add	ebx, esi
	mov	edx, DWORD PTR _DHDA5
	add	ebx, DWORD PTR _7DHDA5+8
	mov	ecx, DWORD PTR _7DH+36
	imul	ecx, edi
	add	ecx, DWORD PTR _7DH+8
	mov	eax, DWORD PTR [edx+ebx*8+4]
	add	ecx, ebp
	mov	edi, DWORD PTR _DH
	mov	DWORD PTR $T98573[esp+1404], 973	; 000003cdH
	mov	ebp, DWORD PTR _OVER+4
	mov	DWORD PTR [edi+ecx*8+4], eax
	mov	eax, DWORD PTR [edx+ebx*8]
	add	ebp, 3
	mov	DWORD PTR [edi+ecx*8], eax

; 973  : 			DH(K+3,I)=DHDAZ(J,I)

	jno	SHORT $L99932
	call	___FHovfl@0
$L99932:
	push	1
	push	OFFSET FLAT:_7DH
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L99933
	push	1
	push	OFFSET FLAT:_7DH
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L99934
$L99933:
	call	___FHranger@0
$L99934:
	mov	edi, DWORD PTR _I$[esp+1392]
	push	2
	push	OFFSET FLAT:_7DH
	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L99935
	push	2
	push	OFFSET FLAT:_7DH
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L99936
$L99935:
	call	___FHranger@0
$L99936:
	mov	ebx, DWORD PTR _I$[esp+1392]
	push	2
	push	OFFSET FLAT:_7DHDAZ
	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L99937
	push	2
	push	OFFSET FLAT:_7DHDAZ
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L99938
$L99937:
	call	___FHranger@0
$L99938:
	imul	ebx, DWORD PTR _7DHDAZ+36
	push	1
	push	OFFSET FLAT:_7DHDAZ
	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L99939
	push	1
	push	OFFSET FLAT:_7DHDAZ
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L99940
$L99939:
	call	___FHranger@0
$L99940:
	add	ebx, esi
	mov	edx, DWORD PTR _DHDAZ
	add	ebx, DWORD PTR _7DHDAZ+8
	mov	ecx, DWORD PTR _7DH+36
	imul	ecx, edi
	add	ecx, DWORD PTR _7DH+8
	mov	eax, DWORD PTR [edx+ebx*8+4]
	add	ecx, ebp
	mov	edi, DWORD PTR _DH
	mov	DWORD PTR $T98573[esp+1404], 974	; 000003ceH
	mov	ebp, DWORD PTR _OVER+4
	mov	DWORD PTR [edi+ecx*8+4], eax
	mov	eax, DWORD PTR [edx+ebx*8]
	add	ebp, 4
	mov	DWORD PTR [edi+ecx*8], eax

; 974  : 			DH(K+4,I)=DHDBZ(J,I)

	jno	SHORT $L99941
	call	___FHovfl@0
$L99941:
	push	1
	push	OFFSET FLAT:_7DH
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L99942
	push	1
	push	OFFSET FLAT:_7DH
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L99943
$L99942:
	call	___FHranger@0
$L99943:
	mov	edi, DWORD PTR _I$[esp+1392]
	push	2
	push	OFFSET FLAT:_7DH
	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L99944
	push	2
	push	OFFSET FLAT:_7DH
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L99945
$L99944:
	call	___FHranger@0
$L99945:
	mov	ebx, DWORD PTR _I$[esp+1392]
	push	2
	push	OFFSET FLAT:_7DHDBZ
	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L99946
	push	2
	push	OFFSET FLAT:_7DHDBZ
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L99947
$L99946:
	call	___FHranger@0
$L99947:
	imul	ebx, DWORD PTR _7DHDBZ+36
	push	1
	push	OFFSET FLAT:_7DHDBZ
	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L99948
	push	1
	push	OFFSET FLAT:_7DHDBZ
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L99949
$L99948:
	call	___FHranger@0
$L99949:
	add	ebx, esi
	mov	edx, DWORD PTR _DHDBZ
	add	ebx, DWORD PTR _7DHDBZ+8
	mov	ecx, DWORD PTR _7DH+36
	imul	ecx, edi
	add	ecx, DWORD PTR _7DH+8
	mov	eax, DWORD PTR [edx+ebx*8+4]
	add	ecx, ebp
	mov	edi, DWORD PTR _DH
	mov	ebp, DWORD PTR _OVER+4
	push	1
	push	OFFSET FLAT:_7DS
	mov	DWORD PTR [edi+ecx*8+4], eax
	mov	DWORD PTR $T98573[esp+1412], 976	; 000003d0H
	mov	eax, DWORD PTR [edx+ebx*8]
	mov	DWORD PTR [edi+ecx*8], eax

; 975  : C
; 976  : 			DS(K,I)=DSDA1(J,I)

	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L99950
	push	1
	push	OFFSET FLAT:_7DS
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L99951
$L99950:
	call	___FHranger@0
$L99951:
	mov	edi, DWORD PTR _I$[esp+1392]
	push	2
	push	OFFSET FLAT:_7DS
	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L99952
	push	2
	push	OFFSET FLAT:_7DS
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L99953
$L99952:
	call	___FHranger@0
$L99953:
	mov	ebx, DWORD PTR _I$[esp+1392]
	push	2
	push	OFFSET FLAT:_7DSDA1
	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L99954
	push	2
	push	OFFSET FLAT:_7DSDA1
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L99955
$L99954:
	call	___FHranger@0
$L99955:
	imul	ebx, DWORD PTR _7DSDA1+36
	push	1
	push	OFFSET FLAT:_7DSDA1
	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L99956
	push	1
	push	OFFSET FLAT:_7DSDA1
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L99957
$L99956:
	call	___FHranger@0
$L99957:
	add	ebx, esi
	mov	ecx, DWORD PTR _7DS+36
	add	ebx, DWORD PTR _7DSDA1+8
	mov	eax, DWORD PTR _DSDA1
	imul	ecx, edi
	add	ecx, DWORD PTR _7DS+8
	mov	eax, DWORD PTR [eax+ebx*8+4]
	add	ecx, ebp
	mov	edx, DWORD PTR _DS
	mov	ebp, DWORD PTR _OVER+4
	mov	DWORD PTR [edx+ecx*8+4], eax
	mov	eax, DWORD PTR _DSDA1
	inc	ebp
	mov	DWORD PTR $T98573[esp+1404], 977	; 000003d1H
	mov	eax, DWORD PTR [eax+ebx*8]
	mov	DWORD PTR [edx+ecx*8], eax

; 977  : 			DS(K+1,I)=DSDA2(J,I)

	jno	SHORT $L99958
	call	___FHovfl@0
$L99958:
	push	1
	push	OFFSET FLAT:_7DS
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L99959
	push	1
	push	OFFSET FLAT:_7DS
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L99960
$L99959:
	call	___FHranger@0
$L99960:
	mov	edi, DWORD PTR _I$[esp+1392]
	push	2
	push	OFFSET FLAT:_7DS
	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L99961
	push	2
	push	OFFSET FLAT:_7DS
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L99962
$L99961:
	call	___FHranger@0
$L99962:
	mov	ebx, DWORD PTR _I$[esp+1392]
	push	2
	push	OFFSET FLAT:_7DSDA2
	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L99963
	push	2
	push	OFFSET FLAT:_7DSDA2
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L99964
$L99963:
	call	___FHranger@0
$L99964:
	imul	ebx, DWORD PTR _7DSDA2+36
	push	1
	push	OFFSET FLAT:_7DSDA2
	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L99965
	push	1
	push	OFFSET FLAT:_7DSDA2
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L99966
$L99965:
	call	___FHranger@0
$L99966:
	add	ebx, esi
	mov	edx, DWORD PTR _DSDA2
	add	ebx, DWORD PTR _7DSDA2+8
	mov	ecx, DWORD PTR _7DS+36
	imul	ecx, edi
	add	ecx, DWORD PTR _7DS+8
	mov	eax, DWORD PTR [edx+ebx*8+4]
	add	ecx, ebp
	mov	edi, DWORD PTR _DS
	mov	DWORD PTR $T98573[esp+1404], 978	; 000003d2H
	mov	ebp, DWORD PTR _OVER+4
	mov	DWORD PTR [edi+ecx*8+4], eax
	mov	eax, DWORD PTR [edx+ebx*8]
	add	ebp, 2
	mov	DWORD PTR [edi+ecx*8], eax

; 978  : 			DS(K+2,I)=DSDA5(J,I)

	jno	SHORT $L99967
	call	___FHovfl@0
$L99967:
	push	1
	push	OFFSET FLAT:_7DS
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L99968
	push	1
	push	OFFSET FLAT:_7DS
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L99969
$L99968:
	call	___FHranger@0
$L99969:
	mov	edi, DWORD PTR _I$[esp+1392]
	push	2
	push	OFFSET FLAT:_7DS
	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L99970
	push	2
	push	OFFSET FLAT:_7DS
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L99971
$L99970:
	call	___FHranger@0
$L99971:
	mov	ebx, DWORD PTR _I$[esp+1392]
	push	2
	push	OFFSET FLAT:_7DSDA5
	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L99972
	push	2
	push	OFFSET FLAT:_7DSDA5
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L99973
$L99972:
	call	___FHranger@0
$L99973:
	imul	ebx, DWORD PTR _7DSDA5+36
	push	1
	push	OFFSET FLAT:_7DSDA5
	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L99974
	push	1
	push	OFFSET FLAT:_7DSDA5
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L99975
$L99974:
	call	___FHranger@0
$L99975:
	add	ebx, esi
	mov	edx, DWORD PTR _DSDA5
	add	ebx, DWORD PTR _7DSDA5+8
	mov	ecx, DWORD PTR _7DS+36
	imul	ecx, edi
	add	ecx, DWORD PTR _7DS+8
	mov	eax, DWORD PTR [edx+ebx*8+4]
	add	ecx, ebp
	mov	edi, DWORD PTR _DS
	mov	DWORD PTR $T98573[esp+1404], 979	; 000003d3H
	mov	ebp, DWORD PTR _OVER+4
	mov	DWORD PTR [edi+ecx*8+4], eax
	mov	eax, DWORD PTR [edx+ebx*8]
	add	ebp, 3
	mov	DWORD PTR [edi+ecx*8], eax

; 979  : 			DS(K+3,I)=DSDAZ(J,I)

	jno	SHORT $L99976
	call	___FHovfl@0
$L99976:
	push	1
	push	OFFSET FLAT:_7DS
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L99977
	push	1
	push	OFFSET FLAT:_7DS
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L99978
$L99977:
	call	___FHranger@0
$L99978:
	mov	edi, DWORD PTR _I$[esp+1392]
	push	2
	push	OFFSET FLAT:_7DS
	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L99979
	push	2
	push	OFFSET FLAT:_7DS
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L99980
$L99979:
	call	___FHranger@0
$L99980:
	mov	ebx, DWORD PTR _I$[esp+1392]
	push	2
	push	OFFSET FLAT:_7DSDAZ
	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L99981
	push	2
	push	OFFSET FLAT:_7DSDAZ
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L99982
$L99981:
	call	___FHranger@0
$L99982:
	imul	ebx, DWORD PTR _7DSDAZ+36
	push	1
	push	OFFSET FLAT:_7DSDAZ
	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L99983
	push	1
	push	OFFSET FLAT:_7DSDAZ
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L99984
$L99983:
	call	___FHranger@0
$L99984:
	add	ebx, esi
	mov	ecx, DWORD PTR _7DS+36
	add	ebx, DWORD PTR _7DSDAZ+8
	mov	eax, DWORD PTR _DSDAZ
	imul	ecx, edi
	add	ecx, DWORD PTR _7DS+8
	mov	eax, DWORD PTR [eax+ebx*8+4]
	add	ecx, ebp
	mov	edx, DWORD PTR _DS
	mov	ebp, DWORD PTR _OVER+4
	mov	DWORD PTR [edx+ecx*8+4], eax
	mov	eax, DWORD PTR _DSDAZ
	add	ebp, 4
	mov	DWORD PTR $T98573[esp+1404], 980	; 000003d4H
	mov	eax, DWORD PTR [eax+ebx*8]
	mov	DWORD PTR [edx+ecx*8], eax

; 980  : 			DS(K+4,I)=DSDBZ(J,I)

	jno	SHORT $L99985
	call	___FHovfl@0
$L99985:
	push	1
	push	OFFSET FLAT:_7DS
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L99986
	push	1
	push	OFFSET FLAT:_7DS
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L99987
$L99986:
	call	___FHranger@0
$L99987:
	mov	edi, DWORD PTR _I$[esp+1392]
	push	2
	push	OFFSET FLAT:_7DS
	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L99988
	push	2
	push	OFFSET FLAT:_7DS
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L99989
$L99988:
	call	___FHranger@0
$L99989:
	mov	ebx, DWORD PTR _I$[esp+1392]
	push	2
	push	OFFSET FLAT:_7DSDBZ
	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L99990
	push	2
	push	OFFSET FLAT:_7DSDBZ
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L99991
$L99990:
	call	___FHranger@0
$L99991:
	imul	ebx, DWORD PTR _7DSDBZ+36
	push	1
	push	OFFSET FLAT:_7DSDBZ
	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L99992
	push	1
	push	OFFSET FLAT:_7DSDBZ
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L99993
$L99992:
	call	___FHranger@0
$L99993:
	add	ebx, esi
	mov	ecx, DWORD PTR _7DS+36
	add	ebx, DWORD PTR _7DSDBZ+8
	mov	eax, DWORD PTR _DSDBZ
	imul	ecx, edi
	add	ecx, DWORD PTR _7DS+8
	mov	eax, DWORD PTR [eax+ebx*8+4]
	add	ecx, ebp
	mov	edx, DWORD PTR _DS
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR [edx+ecx*8+4], eax
	mov	eax, DWORD PTR _DSDBZ
	add	edi, 5
	mov	DWORD PTR $T98573[esp+1404], 982	; 000003d6H
	mov	eax, DWORD PTR [eax+ebx*8]
	mov	DWORD PTR [edx+ecx*8], eax

; 981  : C
; 982  : 			K=K+5

	jno	SHORT $L99994
	call	___FHovfl@0
$L99994:
	inc	esi
	mov	DWORD PTR _OVER+4, edi
	mov	DWORD PTR $T98573[esp+1404], 969	; 000003c9H

; 969  : 		DO J=1,M

	dec	DWORD PTR -600+[esp+1392]
	cmp	DWORD PTR -600+[esp+1392], 0
	jne	$L98407
	mov	DWORD PTR $T98573[esp+1404], 967	; 000003c7H
$L528:

; 966  : C
; 967  : 	DO I=1,M

	mov	eax, DWORD PTR -592+[esp+1392]
	mov	ecx, DWORD PTR _I$[esp+1392]
	dec	eax
	inc	ecx
	mov	DWORD PTR -592+[esp+1392], eax
	mov	DWORD PTR _I$[esp+1392], ecx
	mov	DWORD PTR $T98573[esp+1404], 967	; 000003c7H
	test	eax, eax
	jne	$L98406
	mov	DWORD PTR $T98573[esp+1404], 984	; 000003d8H
$L527:

; 983  : 		END DO
; 984  : 	END DO

	fld	QWORD PTR $T99995
	fdiv	QWORD PTR _SUMS
	mov	edx, DWORD PTR _M$[esp+1388]
	mov	DWORD PTR $T98573[esp+1404], 986	; 000003daH

; 985  : C
; 986  : 	SCALE=1.00D+00/SUMS

	mov	DWORD PTR _I$[esp+1392], 1
	mov	eax, DWORD PTR [edx]
	mov	ecx, 988				; 000003dcH
	test	eax, eax
	mov	DWORD PTR $T98573[esp+1404], ecx
	mov	DWORD PTR -600+[esp+1392], eax
	fstp	QWORD PTR _SCALE

; 987  : C
; 988  : 	DO I=1,M

	jle	$L532
$L98408:

; 989  : 		DO J=1,M*5

	mov	eax, DWORD PTR _M$[esp+1388]
	mov	DWORD PTR $T98573[esp+1404], 989	; 000003ddH
	mov	edi, DWORD PTR [eax]
	imul	edi, 5
	jno	SHORT $L99996
	call	___FHovfl@0
$L99996:
	mov	esi, 1
	mov	DWORD PTR -608+[esp+1392], edi
	test	edi, edi
	jle	$L533
	mov	DWORD PTR $T98573[esp+1404], 989	; 000003ddH
$L98409:
	mov	ebx, DWORD PTR _I$[esp+1392]
	push	1
	mov	DWORD PTR $T98573[esp+1408], 990	; 000003deH

; 990  : 			GRADM(I,J)=DH(J,I)-ENG*DS(J,I)

	push	OFFSET FLAT:_7GRADM
	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L99997
	push	1
	push	OFFSET FLAT:_7GRADM
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L99998
$L99997:
	call	___FHranger@0
$L99998:
	mov	edi, esi
	push	2
	push	OFFSET FLAT:_7GRADM
	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L99999
	push	2
	push	OFFSET FLAT:_7GRADM
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L100000
$L99999:
	call	___FHranger@0
$L100000:
	mov	ebp, DWORD PTR _I$[esp+1392]
	push	2
	push	OFFSET FLAT:_7DS
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L100001
	push	2
	push	OFFSET FLAT:_7DS
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L100002
$L100001:
	call	___FHranger@0
$L100002:
	imul	ebp, DWORD PTR _7DS+36
	push	1
	push	OFFSET FLAT:_7DS
	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L100003
	push	1
	push	OFFSET FLAT:_7DS
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L100004
$L100003:
	call	___FHranger@0
$L100004:
	add	ebp, esi
	mov	eax, DWORD PTR _DS
	add	ebp, DWORD PTR _7DS+8
	push	2
	push	OFFSET FLAT:_7DH
	fld	QWORD PTR [eax+ebp*8]
	fmul	QWORD PTR _ENG
	mov	ebp, DWORD PTR _I$[esp+1400]
	fstp	QWORD PTR -592+[esp+1400]
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L100005
	push	2
	push	OFFSET FLAT:_7DH
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L100006
$L100005:
	call	___FHranger@0
$L100006:
	imul	ebp, DWORD PTR _7DH+36
	push	1
	push	OFFSET FLAT:_7DH
	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L100007
	push	1
	push	OFFSET FLAT:_7DH
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L100008
$L100007:
	call	___FHranger@0
$L100008:
	add	ebp, esi
	mov	eax, DWORD PTR _DH
	add	ebp, DWORD PTR _7DH+8
	mov	ecx, DWORD PTR _7GRADM+36
	imul	ecx, edi
	fld	QWORD PTR [eax+ebp*8]
	fsub	QWORD PTR -592+[esp+1392]
	add	ecx, DWORD PTR _7GRADM+8
	mov	eax, DWORD PTR _GRADM
	add	ecx, ebx
	push	1
	mov	ebx, DWORD PTR _I$[esp+1396]
	push	OFFSET FLAT:_7GRADM
	fstp	QWORD PTR [eax+ecx*8]
	mov	DWORD PTR $T98573[esp+1412], 991	; 000003dfH

; 991  : 			GRADM(I,J)=SCALE*GRADM(I,J)

	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L100009
	push	1
	push	OFFSET FLAT:_7GRADM
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L100010
$L100009:
	call	___FHranger@0
$L100010:
	mov	edi, esi
	push	2
	push	OFFSET FLAT:_7GRADM
	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L100011
	push	2
	push	OFFSET FLAT:_7GRADM
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L100012
$L100011:
	call	___FHranger@0
$L100012:
	mov	ebp, esi
	push	2
	push	OFFSET FLAT:_7GRADM
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L100013
	push	2
	push	OFFSET FLAT:_7GRADM
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L100014
$L100013:
	call	___FHranger@0
$L100014:
	mov	eax, DWORD PTR _I$[esp+1392]
	push	1
	mov	DWORD PTR -592+[esp+1396], eax
	push	OFFSET FLAT:_7GRADM
	imul	ebp, DWORD PTR _7GRADM+36
	call	_ubound_info
	add	esp, 8
	cmp	DWORD PTR -592+[esp+1392], eax
	jg	SHORT $L100015
	push	1
	push	OFFSET FLAT:_7GRADM
	call	_lbound_info
	add	esp, 8
	cmp	DWORD PTR -592+[esp+1392], eax
	jge	SHORT $L100016
$L100015:
	call	___FHranger@0
$L100016:
	add	ebp, DWORD PTR -592+[esp+1392]
	mov	ecx, DWORD PTR _GRADM
	add	ebp, DWORD PTR _7GRADM+8
	mov	eax, DWORD PTR _7GRADM+36
	imul	eax, edi
	fld	QWORD PTR [ecx+ebp*8]
	fmul	QWORD PTR _SCALE
	add	eax, DWORD PTR _7GRADM+8
	inc	esi
	add	eax, ebx
	dec	DWORD PTR -608+[esp+1392]
	mov	DWORD PTR $T98573[esp+1404], 989	; 000003ddH
	fstp	QWORD PTR [ecx+eax*8]

; 989  : 		DO J=1,M*5

	cmp	DWORD PTR -608+[esp+1392], 0
	jne	$L98409
	mov	DWORD PTR $T98573[esp+1404], 988	; 000003dcH
$L533:

; 987  : C
; 988  : 	DO I=1,M

	mov	eax, DWORD PTR -600+[esp+1392]
	mov	ecx, DWORD PTR _I$[esp+1392]
	dec	eax
	inc	ecx
	mov	DWORD PTR -600+[esp+1392], eax
	mov	DWORD PTR _I$[esp+1392], ecx
	mov	DWORD PTR $T98573[esp+1404], 988	; 000003dcH
	test	eax, eax
	jne	$L98408
	mov	DWORD PTR $T98573[esp+1404], 993	; 000003e1H
$L532:

; 992  : 		END DO
; 993  : 	END DO

	mov	DWORD PTR $T98573[esp+1404], 995	; 000003e3H

; 994  : C
; 995  : 	NN=5

	mov	DWORD PTR _NN, 5

; 996  : 	NB=M

	mov	ecx, DWORD PTR _M$[esp+1388]
	mov	edi, 1
	mov	DWORD PTR $T98573[esp+1404], 996	; 000003e4H
	mov	edx, DWORD PTR [ecx]
	mov	ecx, 997				; 000003e5H
	mov	DWORD PTR _NB, edx
	mov	DWORD PTR $T98573[esp+1404], ecx

; 997  : 	DO J=1,NB

	mov	DWORD PTR -600+[esp+1392], edx
	test	edx, edx
	jle	$L538
$L98410:
	mov	eax, 998				; 000003e6H
	mov	ecx, DWORD PTR _NN
	mov	DWORD PTR $T98573[esp+1404], eax
	mov	DWORD PTR -608+[esp+1392], ecx

; 998  : 		DO K=1,NN

	mov	DWORD PTR _OVER+4, 1
	test	ecx, ecx
	jle	$L539
$L98411:
	mov	DWORD PTR $T98573[esp+1404], 999	; 000003e7H

; 999  : 			IC=(J-1)*NN+K

	mov	esi, edi
	dec	esi
	jno	SHORT $L100017
	call	___FHovfl@0
$L100017:
	imul	esi, DWORD PTR _NN
	jno	SHORT $L100018
	call	___FHovfl@0
$L100018:
	add	esi, DWORD PTR _OVER+4
	jno	SHORT $L100019
	call	___FHovfl@0
$L100019:
	push	1
	mov	DWORD PTR _IC, esi
	mov	DWORD PTR $T98573[esp+1408], 1000	; 000003e8H

; 1000 : 			GRADA(IC)=0.00D+00

	push	OFFSET FLAT:_7GRADA
	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L100020
	push	1
	push	OFFSET FLAT:_7GRADA
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L100021
$L100020:
	call	___FHranger@0
$L100021:
	mov	ecx, DWORD PTR _7GRADA+8
	mov	eax, DWORD PTR _GRADA
	add	ecx, esi
	mov	esi, 1
	mov	DWORD PTR [eax+ecx*8], 0
	mov	DWORD PTR [eax+ecx*8+4], 0
	mov	ecx, 1001				; 000003e9H
	mov	eax, DWORD PTR _NB
	mov	DWORD PTR $T98573[esp+1404], ecx
	mov	DWORD PTR -616+[esp+1392], eax

; 1001 : 			DO I=1,NB

	test	eax, eax
	jle	$L541
$L98412:
	mov	DWORD PTR $T98573[esp+1404], 1002	; 000003eaH

; 1002 : 				IF(I.NE.J) THEN

	cmp	esi, edi
	je	$L546
	mov	ebx, DWORD PTR _IC
	push	1
	mov	DWORD PTR $T98573[esp+1408], 1003	; 000003ebH

; 1003 : 					GRADA(IC)=GRADA(IC)+2.00D+00*EVEC(I)*EVEC(J)*

	push	OFFSET FLAT:_7GRADA
	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L100022
	push	1
	push	OFFSET FLAT:_7GRADA
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L100023
$L100022:
	call	___FHranger@0
$L100023:
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	ebp, DWORD PTR _IC
	push	2
	push	OFFSET FLAT:_7GRADM
	mov	ecx, DWORD PTR [eax]
	mov	DWORD PTR -592+[esp+1400], ecx
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L100024
	push	2
	push	OFFSET FLAT:_7GRADM
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L100025
$L100024:
	call	___FHranger@0
$L100025:
	imul	ebp, DWORD PTR _7GRADM+36
	push	1
	push	OFFSET FLAT:_7GRADM
	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L100026
	push	1
	push	OFFSET FLAT:_7GRADM
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L100027
$L100026:
	call	___FHranger@0
$L100027:
	add	ebp, esi
	mov	eax, DWORD PTR _7GRADM+8
	add	ebp, eax
	mov	ecx, DWORD PTR _GRADM
	mov	edx, DWORD PTR -592+[esp+1392]
	fld	QWORD PTR [ecx+ebp*8]
	cmp	edi, edx
	jg	SHORT $L100028
	cmp	edi, 1
	jge	SHORT $L100029
$L100028:
	call	___FHranger@0
$L100029:
	mov	eax, DWORD PTR _EVEC$[esp+1388]
	cmp	esi, DWORD PTR -592+[esp+1392]
	fmul	QWORD PTR [eax+edi*8-8]
	jg	SHORT $L100030
	cmp	esi, 1
	jge	SHORT $L100031
$L100030:
	call	___FHranger@0
$L100031:
	mov	eax, DWORD PTR _EVEC$[esp+1388]
	mov	ebp, DWORD PTR _IC
	push	1
	push	OFFSET FLAT:_7GRADA
	fmul	QWORD PTR [eax+esi*8-8]
	fmul	QWORD PTR $T99039
	fstp	QWORD PTR -592+[esp+1400]
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L100032
	push	1
	push	OFFSET FLAT:_7GRADA
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L100033
$L100032:
	call	___FHranger@0
$L100033:
	mov	eax, DWORD PTR _7GRADA+8
	mov	ecx, DWORD PTR _GRADA
	add	ebp, eax
	add	eax, ebx
	fld	QWORD PTR [ecx+ebp*8]
	fadd	QWORD PTR -592+[esp+1392]
	fstp	QWORD PTR [ecx+eax*8]
	mov	eax, 1005				; 000003edH
	mov	DWORD PTR $T98573[esp+1404], eax

; 1004 : 	-						GRADM(I,IC)
; 1005 : 				ELSE

	jmp	$L544
$L546:
	mov	ebx, DWORD PTR _IC
	push	1
	mov	DWORD PTR $T98573[esp+1408], 1006	; 000003eeH

; 1006 : 					GRADA(IC)=GRADA(IC)+EVEC(I)*EVEC(J)*GRADM(I,IC)

	push	OFFSET FLAT:_7GRADA
	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L100034
	push	1
	push	OFFSET FLAT:_7GRADA
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L100035
$L100034:
	call	___FHranger@0
$L100035:
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	ebp, DWORD PTR _IC
	push	2
	push	OFFSET FLAT:_7GRADM
	mov	ecx, DWORD PTR [eax]
	mov	DWORD PTR -592+[esp+1400], ecx
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L100036
	push	2
	push	OFFSET FLAT:_7GRADM
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L100037
$L100036:
	call	___FHranger@0
$L100037:
	imul	ebp, DWORD PTR _7GRADM+36
	push	1
	push	OFFSET FLAT:_7GRADM
	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L100038
	push	1
	push	OFFSET FLAT:_7GRADM
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L100039
$L100038:
	call	___FHranger@0
$L100039:
	add	ebp, esi
	mov	eax, DWORD PTR _7GRADM+8
	add	ebp, eax
	mov	ecx, DWORD PTR _GRADM
	mov	edx, DWORD PTR -592+[esp+1392]
	fld	QWORD PTR [ecx+ebp*8]
	cmp	edi, edx
	jg	SHORT $L100040
	cmp	edi, 1
	jge	SHORT $L100041
$L100040:
	call	___FHranger@0
$L100041:
	mov	eax, DWORD PTR _EVEC$[esp+1388]
	cmp	esi, DWORD PTR -592+[esp+1392]
	fmul	QWORD PTR [eax+edi*8-8]
	jg	SHORT $L100042
	cmp	esi, 1
	jge	SHORT $L100043
$L100042:
	call	___FHranger@0
$L100043:
	mov	eax, DWORD PTR _EVEC$[esp+1388]
	mov	ebp, DWORD PTR _IC
	push	1
	push	OFFSET FLAT:_7GRADA
	fmul	QWORD PTR [eax+esi*8-8]
	fstp	QWORD PTR -592+[esp+1400]
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L100044
	push	1
	push	OFFSET FLAT:_7GRADA
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L100045
$L100044:
	call	___FHranger@0
$L100045:
	mov	eax, DWORD PTR _7GRADA+8
	mov	ecx, DWORD PTR _GRADA
	mov	DWORD PTR $T98573[esp+1404], 1001	; 000003e9H
	add	ebp, eax
	add	eax, ebx
	fld	QWORD PTR [ecx+ebp*8]
	fadd	QWORD PTR -592+[esp+1392]
	fstp	QWORD PTR [ecx+eax*8]
$L544:

; 1001 : 			DO I=1,NB

	mov	eax, DWORD PTR -616+[esp+1392]
	inc	esi
	mov	DWORD PTR $T98573[esp+1404], 1001	; 000003e9H
	dec	eax
	mov	DWORD PTR -616+[esp+1392], eax
	test	eax, eax
	jne	$L98412
	mov	DWORD PTR $T98573[esp+1404], 998	; 000003e6H
$L541:

; 998  : 		DO K=1,NN

	mov	eax, DWORD PTR -608+[esp+1392]
	mov	ecx, DWORD PTR _OVER+4
	dec	eax
	inc	ecx
	mov	DWORD PTR -608+[esp+1392], eax
	mov	DWORD PTR _OVER+4, ecx
	mov	DWORD PTR $T98573[esp+1404], 998	; 000003e6H
	test	eax, eax
	jne	$L98411
	mov	DWORD PTR $T98573[esp+1404], 997	; 000003e5H
$L539:

; 997  : 	DO J=1,NB

	mov	eax, DWORD PTR -600+[esp+1392]
	inc	edi
	mov	DWORD PTR $T98573[esp+1404], 997	; 000003e5H
	dec	eax
	mov	DWORD PTR -600+[esp+1392], eax
	test	eax, eax
	jne	$L98410
$L538:

; 1007 : 				END IF
; 1008 : 			END DO
; 1009 : 		END DO
; 1010 : 	END DO

	mov	eax, 1012				; 000003f4H
	mov	edi, DWORD PTR _NB
	mov	esi, 1
	mov	DWORD PTR $T98573[esp+1404], eax

; 1011 : C
; 1012 :       DO I=1,NB

	test	edi, edi
	jle	SHORT $L547
$L98413:
	mov	DWORD PTR $T98573[esp+1404], 1013	; 000003f5H

; 1013 : 		GRADC(I) = 0.00D+00

	push	1
	push	OFFSET FLAT:_7GRADC
	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L100046
	push	1
	push	OFFSET FLAT:_7GRADC
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L100047
$L100046:
	call	___FHranger@0
$L100047:
	mov	ecx, DWORD PTR _7GRADC+8
	mov	eax, DWORD PTR _GRADC
	add	ecx, esi
	dec	edi
	mov	DWORD PTR $T98573[esp+1404], 1012	; 000003f4H
	inc	esi
	mov	DWORD PTR [eax+ecx*8], 0
	mov	DWORD PTR [eax+ecx*8+4], 0

; 1011 : C
; 1012 :       DO I=1,NB

	test	edi, edi
	jne	SHORT $L98413
	mov	DWORD PTR $T98573[esp+1404], 1014	; 000003f6H
$L547:

; 1014 : 	END DO

	mov	ecx, 1016				; 000003f8H
	mov	eax, DWORD PTR _NB
	mov	DWORD PTR $T98573[esp+1404], ecx
	mov	DWORD PTR -600+[esp+1392], eax

; 1015 : 	
; 1016 :       DO J=1,NB

	mov	DWORD PTR _J$[esp+1392], 1
	test	eax, eax
	jle	$L549
$L98414:
	mov	ecx, 1017				; 000003f9H
	mov	eax, DWORD PTR _NB
	mov	esi, 1
	mov	DWORD PTR $T98573[esp+1404], ecx

; 1017 : 		DO I=1,NB

	mov	DWORD PTR -608+[esp+1392], eax
	test	eax, eax
	jle	$L550
$L98415:
	mov	ebp, esi
	push	1
	mov	DWORD PTR $T98573[esp+1408], 1018	; 000003faH

; 1018 : 		  GRADC(I)=GRADC(I)+2.00D+00*SCALE*( H(I,J)-ENG*S(I,J) )*

	push	OFFSET FLAT:_7GRADC
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L100048
	push	1
	push	OFFSET FLAT:_7GRADC
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L100049
$L100048:
	call	___FHranger@0
$L100049:
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	ebx, DWORD PTR _J$[esp+1392]
	mov	edi, DWORD PTR [eax]
	cmp	ebx, edi
	jg	SHORT $L100050
	cmp	ebx, 1
	jge	SHORT $L100051
$L100050:
	call	___FHranger@0
$L100051:
	imul	ebx, DWORD PTR _7S$[esp+1428]
	cmp	esi, edi
	jg	SHORT $L100052
	cmp	esi, 1
	jge	SHORT $L100053
$L100052:
	call	___FHranger@0
$L100053:
	add	ebx, esi
	mov	eax, DWORD PTR _S$[esp+1388]
	add	ebx, DWORD PTR _7S$[esp+1400]
	fld	QWORD PTR [eax+ebx*8]
	fmul	QWORD PTR _ENG
	mov	ebx, DWORD PTR _J$[esp+1392]
	cmp	ebx, edi
	jg	SHORT $L100054
	cmp	ebx, 1
	jge	SHORT $L100055
$L100054:
	call	___FHranger@0
$L100055:
	imul	ebx, DWORD PTR _7H$[esp+1428]
	cmp	esi, edi
	jg	SHORT $L100056
	cmp	esi, 1
	jge	SHORT $L100057
$L100056:
	call	___FHranger@0
$L100057:
	add	ebx, esi
	mov	eax, DWORD PTR _H$[esp+1388]
	add	ebx, DWORD PTR _7H$[esp+1400]
	fsubr	QWORD PTR [eax+ebx*8]
	mov	ebx, DWORD PTR _J$[esp+1392]
	cmp	ebx, edi
	jg	SHORT $L100058
	cmp	ebx, 1
	jge	SHORT $L100059
$L100058:
	call	___FHranger@0
$L100059:
	mov	eax, DWORD PTR _EVEC$[esp+1388]
	mov	edi, esi
	push	1
	push	OFFSET FLAT:_7GRADC
	fmul	QWORD PTR [eax+ebx*8-8]
	fmul	QWORD PTR _SCALE
	fmul	QWORD PTR $T99039
	fstp	QWORD PTR -592+[esp+1400]
	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L100060
	push	1
	push	OFFSET FLAT:_7GRADC
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L100061
$L100060:
	call	___FHranger@0
$L100061:
	add	edi, DWORD PTR _7GRADC+8
	mov	eax, DWORD PTR _GRADC
	mov	ecx, DWORD PTR _7GRADC+8
	inc	esi
	fld	QWORD PTR [eax+edi*8]
	fadd	QWORD PTR -592+[esp+1392]
	add	ecx, ebp
	dec	DWORD PTR -608+[esp+1392]
	mov	DWORD PTR $T98573[esp+1404], 1017	; 000003f9H
	cmp	DWORD PTR -608+[esp+1392], 0
	fstp	QWORD PTR [eax+ecx*8]

; 1017 : 		DO I=1,NB

	jne	$L98415
	mov	DWORD PTR $T98573[esp+1404], 1016	; 000003f8H
$L550:

; 1015 : 	
; 1016 :       DO J=1,NB

	mov	eax, DWORD PTR -600+[esp+1392]
	mov	ecx, DWORD PTR _J$[esp+1392]
	dec	eax
	inc	ecx
	mov	DWORD PTR -600+[esp+1392], eax
	mov	DWORD PTR _J$[esp+1392], ecx
	mov	DWORD PTR $T98573[esp+1404], 1016	; 000003f8H
	test	eax, eax
	jne	$L98414
	mov	DWORD PTR $T98573[esp+1404], 1021	; 000003fdH
$L549:

; 1019 : 	-				EVEC(J)
; 1020 : 		END DO
; 1021 : 	END DO

	mov	DWORD PTR $T98573[esp+1404], 1023	; 000003ffH

; 1022 : 
; 1023 : 	K=1

	mov	DWORD PTR _OVER+4, 1
	mov	DWORD PTR $T98573[esp+1404], 1024	; 00000400H

; 1024 : 	DO I=1,NB*5,5

	mov	esi, DWORD PTR _NB
	imul	esi, 5
	jno	SHORT $L100062
	call	___FHovfl@0
$L100062:
	lea	eax, DWORD PTR [esi+4]
	mov	ecx, 5
	cdq
	idiv	ecx
	mov	edi, 1
	mov	DWORD PTR -592+[esp+1392], eax
	test	eax, eax
	jle	$L553
	mov	DWORD PTR $T98573[esp+1404], 1024	; 00000400H
$L98416:
	mov	esi, edi
	push	1
	mov	DWORD PTR $T98573[esp+1408], 1025	; 00000401H

; 1025 : 		DA1=GRADA(I)

	push	OFFSET FLAT:_7GRADA
	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L100063
	push	1
	push	OFFSET FLAT:_7GRADA
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L100064
$L100063:
	call	___FHranger@0
$L100064:
	add	esi, DWORD PTR _7GRADA+8
	mov	eax, DWORD PTR _GRADA
	mov	eax, DWORD PTR [eax+esi*8+4]
	mov	ebx, edi
	inc	ebx
	mov	DWORD PTR _DA1+4, eax
	mov	eax, DWORD PTR _GRADA
	mov	DWORD PTR $T98573[esp+1404], 1026	; 00000402H
	mov	eax, DWORD PTR [eax+esi*8]
	mov	esi, DWORD PTR _GRADA
	mov	DWORD PTR _DA1, eax

; 1026 : 		DA2=GRADA(I+1)

	jno	SHORT $L100065
	call	___FHovfl@0
$L100065:
	push	1
	push	OFFSET FLAT:_7GRADA
	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L100066
	push	1
	push	OFFSET FLAT:_7GRADA
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L100067
$L100066:
	call	___FHranger@0
$L100067:
	add	ebx, DWORD PTR _7GRADA+8
	mov	DWORD PTR $T98573[esp+1404], 1027	; 00000403H
	mov	eax, DWORD PTR [esi+ebx*8+4]
	mov	DWORD PTR _DA2+4, eax
	mov	eax, DWORD PTR [esi+ebx*8]
	mov	ebx, edi
	mov	esi, DWORD PTR _GRADA
	add	ebx, 2
	mov	DWORD PTR _DA2, eax

; 1027 : 		DB=GRADA(I+2)

	jno	SHORT $L100068
	call	___FHovfl@0
$L100068:
	push	1
	push	OFFSET FLAT:_7GRADA
	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L100069
	push	1
	push	OFFSET FLAT:_7GRADA
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L100070
$L100069:
	call	___FHranger@0
$L100070:
	add	ebx, DWORD PTR _7GRADA+8
	push	1
	mov	DWORD PTR $T98573[esp+1408], 1028	; 00000404H
	mov	eax, DWORD PTR [esi+ebx*8+4]
	push	OFFSET FLAT:_7GRADA
	mov	DWORD PTR _DB+4, eax
	mov	eax, DWORD PTR [esi+ebx*8]
	mov	DWORD PTR _DB, eax

; 1028 : 	GRADA(I)=DA1*(2.00D+00*TTT1(K)+TTT3(K))+DA2*TTT3(K)-DB*TTT3(K)

	call	_ubound_info
	add	esp, 8
	cmp	edi, eax
	jg	SHORT $L100071
	push	1
	push	OFFSET FLAT:_7GRADA
	call	_lbound_info
	add	esp, 8
	cmp	edi, eax
	jge	SHORT $L100072
$L100071:
	call	___FHranger@0
$L100072:
	mov	ebp, edi
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	ebx, DWORD PTR _OVER+4
	mov	esi, DWORD PTR [eax]
	cmp	ebx, esi
	jg	SHORT $L100073
	cmp	ebx, 1
	jge	SHORT $L100074
$L100073:
	call	___FHranger@0
$L100074:
	mov	eax, DWORD PTR _TTT1$[esp+1388]
	fld	QWORD PTR [eax+ebx*8-8]
	fmul	QWORD PTR $T99039
	mov	ebx, DWORD PTR _OVER+4
	cmp	ebx, esi
	jg	SHORT $L100075
	cmp	ebx, 1
	jge	SHORT $L100076
$L100075:
	call	___FHranger@0
$L100076:
	mov	eax, DWORD PTR _TTT3$[esp+1388]
	fadd	QWORD PTR [eax+ebx*8-8]
	fmul	QWORD PTR _DA1
	mov	ebx, DWORD PTR _OVER+4
	cmp	ebx, esi
	jg	SHORT $L100077
	cmp	ebx, 1
	jge	SHORT $L100078
$L100077:
	call	___FHranger@0
$L100078:
	mov	eax, DWORD PTR _TTT3$[esp+1388]
	fld	QWORD PTR [eax+ebx*8-8]
	fmul	QWORD PTR _DB
	mov	ebx, DWORD PTR _OVER+4
	cmp	ebx, esi
	fsubp	ST(1), ST(0)
	jg	SHORT $L100079
	cmp	ebx, 1
	jge	SHORT $L100080
$L100079:
	call	___FHranger@0
$L100080:
	mov	eax, DWORD PTR _TTT3$[esp+1388]
	mov	ecx, DWORD PTR _7GRADA+8
	add	ecx, ebp
	mov	esi, edi
	fld	QWORD PTR [eax+ebx*8-8]
	fmul	QWORD PTR _DA2
	faddp	ST(1), ST(0)
	mov	eax, DWORD PTR _GRADA
	inc	esi
	mov	DWORD PTR $T98573[esp+1404], 1029	; 00000405H
	fstp	QWORD PTR [eax+ecx*8]

; 1029 : 	GRADA(I+1)=2.00D+00*TTT2(K)*DA2

	jno	SHORT $L100081
	call	___FHovfl@0
$L100081:
	push	1
	push	OFFSET FLAT:_7GRADA
	call	_ubound_info
	add	esp, 8
	cmp	esi, eax
	jg	SHORT $L100082
	push	1
	push	OFFSET FLAT:_7GRADA
	call	_lbound_info
	add	esp, 8
	cmp	esi, eax
	jge	SHORT $L100083
$L100082:
	call	___FHranger@0
$L100083:
	mov	ebx, DWORD PTR _OVER+4
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	ebp, DWORD PTR _GRADA
	npad	1
	cmp	ebx, DWORD PTR [eax]
	jg	SHORT $L100084
	cmp	ebx, 1
	jge	SHORT $L100085
$L100084:
	call	___FHranger@0
$L100085:
	mov	DWORD PTR $T98573[esp+1404], 1030	; 00000406H
	mov	eax, DWORD PTR _TTT2$[esp+1388]
	fld	QWORD PTR [eax+ebx*8-8]
	fmul	QWORD PTR _DA2
	mov	eax, DWORD PTR _7GRADA+8
	mov	ebx, edi
	add	eax, esi
	add	ebx, 2
	fmul	QWORD PTR $T99039
	fstp	QWORD PTR [ebp+eax*8]

; 1030 : 	GRADA(I+2)=DA1*TTT1(K)+(2.00D+00*TTT3(K)+TTT1(K))*DA2-DB*TTT1(K)

	jno	SHORT $L100086
	call	___FHovfl@0
$L100086:
	push	1
	push	OFFSET FLAT:_7GRADA
	call	_ubound_info
	add	esp, 8
	cmp	ebx, eax
	jg	SHORT $L100087
	push	1
	push	OFFSET FLAT:_7GRADA
	call	_lbound_info
	add	esp, 8
	cmp	ebx, eax
	jge	SHORT $L100088
$L100087:
	call	___FHranger@0
$L100088:
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	ebp, DWORD PTR _OVER+4
	mov	esi, DWORD PTR [eax]
	cmp	ebp, esi
	jg	SHORT $L100089
	cmp	ebp, 1
	jge	SHORT $L100090
$L100089:
	call	___FHranger@0
$L100090:
	mov	eax, DWORD PTR _TTT3$[esp+1388]
	fld	QWORD PTR [eax+ebp*8-8]
	fmul	QWORD PTR $T99039
	mov	ebp, DWORD PTR _OVER+4
	cmp	ebp, esi
	jg	SHORT $L100091
	cmp	ebp, 1
	jge	SHORT $L100092
$L100091:
	call	___FHranger@0
$L100092:
	mov	eax, DWORD PTR _TTT1$[esp+1388]
	fadd	QWORD PTR [eax+ebp*8-8]
	fmul	QWORD PTR _DA2
	mov	ebp, DWORD PTR _OVER+4
	cmp	ebp, esi
	jg	SHORT $L100093
	cmp	ebp, 1
	jge	SHORT $L100094
$L100093:
	call	___FHranger@0
$L100094:
	mov	eax, DWORD PTR _TTT1$[esp+1388]
	fld	QWORD PTR [eax+ebp*8-8]
	fmul	QWORD PTR _DB
	mov	ebp, DWORD PTR _OVER+4
	cmp	ebp, esi
	fsubp	ST(1), ST(0)
	jg	SHORT $L100095
	cmp	ebp, 1
	jge	SHORT $L100096
$L100095:
	call	___FHranger@0
$L100096:
	mov	eax, DWORD PTR _TTT1$[esp+1388]
	mov	ecx, DWORD PTR _7GRADA+8
	add	ecx, ebx
	mov	esi, DWORD PTR _OVER+4
	fld	QWORD PTR [eax+ebp*8-8]
	fmul	QWORD PTR _DA1
	faddp	ST(1), ST(0)
	mov	eax, DWORD PTR _GRADA
	inc	esi
	mov	DWORD PTR $T98573[esp+1404], 1031	; 00000407H
	fstp	QWORD PTR [eax+ecx*8]

; 1031 : 	K=K+1

	jno	SHORT $L100097
	call	___FHovfl@0
$L100097:
	add	edi, 5
	mov	DWORD PTR _OVER+4, esi
	mov	DWORD PTR $T98573[esp+1404], 1024	; 00000400H

; 1024 : 	DO I=1,NB*5,5

	dec	DWORD PTR -592+[esp+1392]
	cmp	DWORD PTR -592+[esp+1392], 0
	jne	$L98416
	mov	DWORD PTR $T98573[esp+1404], 1032	; 00000408H
$L553:

; 1032 :       END DO
; 1033 : C
; 1034 :       NX=M*5

	mov	ecx, DWORD PTR _M$[esp+1388]
	mov	DWORD PTR $T98573[esp+1404], 1034	; 0000040aH
	mov	esi, DWORD PTR [ecx]
	imul	esi, 5
	jno	SHORT $L100098
	call	___FHovfl@0
$L100098:
	mov	eax, 1035				; 0000040bH
	mov	edi, esi
	mov	DWORD PTR _NX, esi
	mov	DWORD PTR $T98573[esp+1404], eax

; 1035 : 	DO I=1,NX

	mov	DWORD PTR _I$[esp+1392], 1
	test	esi, esi
	mov	esi, DWORD PTR _GRAD$[esp+1388]
	jle	$L98570
$L98417:
	mov	ebx, DWORD PTR _I$[esp+1392]
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	DWORD PTR $T98573[esp+1404], 1036	; 0000040cH

; 1036 : 		GRAD(I) = GRADA(I)

	mov	ebp, DWORD PTR [eax]
	imul	ebp, 6
	jno	SHORT $L100099
	call	___FHovfl@0
$L100099:
	cmp	ebx, ebp
	jg	SHORT $L100100
	cmp	ebx, 1
	jge	SHORT $L100101
$L100100:
	call	___FHranger@0
$L100101:
	mov	ebp, DWORD PTR _I$[esp+1392]
	push	1
	push	OFFSET FLAT:_7GRADA
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L100102
	push	1
	push	OFFSET FLAT:_7GRADA
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L100103
$L100102:
	call	___FHranger@0
$L100103:
	add	ebp, DWORD PTR _7GRADA+8
	mov	eax, DWORD PTR _GRADA
	dec	edi
	inc	DWORD PTR _I$[esp+1392]
	mov	eax, DWORD PTR [eax+ebp*8+4]
	test	edi, edi
	mov	DWORD PTR [esi+ebx*8-4], eax
	mov	eax, DWORD PTR _GRADA
	mov	DWORD PTR $T98573[esp+1404], 1035	; 0000040bH
	mov	eax, DWORD PTR [eax+ebp*8]
	mov	DWORD PTR [esi+ebx*8-8], eax

; 1035 : 	DO I=1,NX

	jne	$L98417
$L98570:

; 1037 : 	END DO

	mov	eax, 1039				; 0000040fH
	mov	edi, DWORD PTR _NB
	mov	DWORD PTR $T98573[esp+1404], eax
	test	edi, edi

; 1038 : 	
; 1039 :       DO I=1,NB

	mov	DWORD PTR _I$[esp+1392], 1
	jle	$L561
$L98418:
	mov	DWORD PTR $T98573[esp+1404], 1040	; 00000410H

; 1040 : 		GRAD(I+NX) = GRADC(I)

	mov	ebx, DWORD PTR _I$[esp+1392]
	add	ebx, DWORD PTR _NX
	jno	SHORT $L100104
	call	___FHovfl@0
$L100104:
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	ebp, DWORD PTR [eax]
	imul	ebp, 6
	jno	SHORT $L100105
	call	___FHovfl@0
$L100105:
	cmp	ebx, ebp
	jg	SHORT $L100106
	cmp	ebx, 1
	jge	SHORT $L100107
$L100106:
	call	___FHranger@0
$L100107:
	mov	ebp, DWORD PTR _I$[esp+1392]
	push	1
	push	OFFSET FLAT:_7GRADC
	call	_ubound_info
	add	esp, 8
	cmp	ebp, eax
	jg	SHORT $L100108
	push	1
	push	OFFSET FLAT:_7GRADC
	call	_lbound_info
	add	esp, 8
	cmp	ebp, eax
	jge	SHORT $L100109
$L100108:
	call	___FHranger@0
$L100109:
	add	ebp, DWORD PTR _7GRADC+8
	mov	ecx, DWORD PTR _GRADC
	dec	edi
	inc	DWORD PTR _I$[esp+1392]
	mov	eax, DWORD PTR [ecx+ebp*8+4]
	test	edi, edi
	mov	DWORD PTR $T98573[esp+1404], 1039	; 0000040fH
	mov	DWORD PTR [esi+ebx*8-4], eax
	mov	eax, DWORD PTR [ecx+ebp*8]
	mov	DWORD PTR [esi+ebx*8-8], eax

; 1038 : 	
; 1039 :       DO I=1,NB

	jne	$L98418
	mov	DWORD PTR $T98573[esp+1404], 1041	; 00000411H
$L561:

; 1041 : 	END DO

	xor	esi, esi
	mov	ecx, DWORD PTR _DHDA1
	mov	DWORD PTR $T98573[esp+1404], 1043	; 00000413H

; 1042 : C     -----------------------------------------------------      
; 1043 :       DEALLOCATE(	TEMP,DHDA1,DHDA2,DHDA5,DHDAZ,DHDBZ,DSDA1,DSDA2,

	push	esi
	push	8
	push	OFFSET FLAT:_7DHDA1
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _DHDA2
	mov	DWORD PTR _DHDA1, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7DHDA2
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _DHDA5
	mov	DWORD PTR _DHDA2, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7DHDA5
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _DHDAZ
	mov	DWORD PTR _DHDA5, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7DHDAZ
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _DHDBZ
	mov	DWORD PTR _DHDAZ, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7DHDBZ
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _DSDA1
	mov	DWORD PTR _DHDBZ, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7DSDA1
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _DSDA2
	mov	DWORD PTR _DSDA1, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7DSDA2
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _DSDA5
	mov	DWORD PTR _DSDA2, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7DSDA5
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _DSDAZ
	mov	DWORD PTR _DSDA5, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7DSDAZ
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _DSDBZ
	mov	DWORD PTR _DSDAZ, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7DSDBZ
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _GRADM
	mov	DWORD PTR _DSDBZ, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7GRADM
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _DS
	mov	DWORD PTR _GRADM, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7DS
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	DWORD PTR _DS, esi
	push	esi
	push	8
	push	OFFSET FLAT:_7DH
	mov	eax, DWORD PTR _DH
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _GRADA
	mov	DWORD PTR _DH, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7GRADA
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _GRADC
	mov	DWORD PTR _GRADA, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7GRADC
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _TEMP
	mov	DWORD PTR _GRADC, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7TEMP
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, 1048				; 00000418H
	mov	DWORD PTR _TEMP, esi
	npad	1
	mov	DWORD PTR $T98573[esp+1404], eax

; 1044 : 	-			DSDA5,DSDAZ,DSDBZ,GRADM,DS,DH,GRADA,GRADC)
; 1045 : C
; 1046 : C
; 1047 : C
; 1048 : 	ELSE

	jmp	$L652
$L148:
	mov	ecx, 1052				; 0000041cH
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	DWORD PTR _OVER, 1
	mov	edx, DWORD PTR [eax]
	mov	DWORD PTR $T98573[esp+1404], ecx

; 1049 : C
; 1050 : C
; 1051 : C
; 1052 :       DO L=1,M

	mov	DWORD PTR -1364+[esp+1392], edx
	test	edx, edx
	jle	$L564
$L98419:
	mov	eax, 1053				; 0000041dH
	mov	ecx, DWORD PTR _OVER
	mov	DWORD PTR $T98573[esp+1404], eax
	mov	DWORD PTR -1372+[esp+1392], ecx

; 1053 :           DO K=1,L

	mov	DWORD PTR _OVER+4, 1
	test	ecx, ecx
	jle	$L565
$L98420:
	mov	eax, DWORD PTR _M$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_OV1
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1412], 1053	; 0000041dH
	mov	DWORD PTR $T98573[esp+1412], 1057	; 00000421H
	mov	esi, DWORD PTR [eax]

; 1054 : C     -----------------------------------------------------
; 1055 : C     --- Calculate OVERLAP INTEGRAL ---
; 1056 : C     -----------------------------------------------------
; 1057 :              CALL OVERLP(M,DFLAG,RX1(L),RY1(L),RZ1(L),AA1(L),

	cmp	edi, esi
	jg	SHORT $L100110
	cmp	edi, 1
	jge	SHORT $L100111
$L100110:
	call	___FHranger@0
$L100111:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100112
	cmp	edi, 1
	jge	SHORT $L100113
$L100112:
	call	___FHranger@0
$L100113:
	mov	eax, DWORD PTR _B12$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100114
	cmp	edi, 1
	jge	SHORT $L100115
$L100114:
	call	___FHranger@0
$L100115:
	mov	eax, DWORD PTR _AA2$[esp+1404]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100116
	cmp	edi, 1
	jge	SHORT $L100117
$L100116:
	call	___FHranger@0
$L100117:
	mov	eax, DWORD PTR _RZ2$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100118
	cmp	edi, 1
	jge	SHORT $L100119
$L100118:
	call	___FHranger@0
$L100119:
	mov	eax, DWORD PTR _RY2$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100120
	cmp	edi, 1
	jge	SHORT $L100121
$L100120:
	call	___FHranger@0
$L100121:
	mov	eax, DWORD PTR _RX2$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100122
	cmp	edi, 1
	jge	SHORT $L100123
$L100122:
	call	___FHranger@0
$L100123:
	mov	eax, DWORD PTR _AA2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100124
	cmp	edi, 1
	jge	SHORT $L100125
$L100124:
	call	___FHranger@0
$L100125:
	mov	eax, DWORD PTR _RZ2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100126
	cmp	edi, 1
	jge	SHORT $L100127
$L100126:
	call	___FHranger@0
$L100127:
	mov	eax, DWORD PTR _RY2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100128
	cmp	edi, 1
	jge	SHORT $L100129
$L100128:
	call	___FHranger@0
$L100129:
	mov	eax, DWORD PTR _RX2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100130
	cmp	edi, 1
	jge	SHORT $L100131
$L100130:
	call	___FHranger@0
$L100131:
	mov	eax, DWORD PTR _AA1$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100132
	cmp	edi, 1
	jge	SHORT $L100133
$L100132:
	call	___FHranger@0
$L100133:
	mov	eax, DWORD PTR _RZ1$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100134
	cmp	edi, 1
	jge	SHORT $L100135
$L100134:
	call	___FHranger@0
$L100135:
	mov	eax, DWORD PTR _RY1$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100136
	cmp	edi, 1
	jge	SHORT $L100137
$L100136:
	call	___FHranger@0
$L100137:
	mov	eax, DWORD PTR _RX1$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100138
	cmp	edi, 1
	jge	SHORT $L100139
$L100138:
	call	___FHranger@0
$L100139:
	mov	eax, DWORD PTR _AA1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100140
	cmp	edi, 1
	jge	SHORT $L100141
$L100140:
	call	___FHranger@0
$L100141:
	mov	eax, DWORD PTR _RZ1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100142
	cmp	edi, 1
	jge	SHORT $L100143
$L100142:
	call	___FHranger@0
$L100143:
	mov	eax, DWORD PTR _RY1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100144
	cmp	edi, 1
	jge	SHORT $L100145
$L100144:
	call	___FHranger@0
$L100145:
	mov	eax, DWORD PTR _RX1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1464]
	sub	edi, 8
	push	edi
	push	eax
	mov	eax, DWORD PTR _M$[esp+1472]
	push	eax
	call	_OVERLP@88
	mov	eax, DWORD PTR _M$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_OV2
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1412], 1065	; 00000429H
	mov	esi, DWORD PTR [eax]

; 1058 :      1            RX1(K),RY1(K),RZ1(K),AA1(K),                      
; 1059 :      2            RX2(L),RY2(L),RZ2(L),AA2(L),
; 1060 :      3            RX2(K),RY2(K),RZ2(K),AA2(K),
; 1061 :      4            B12(L),B12(K),OV1,KLM)
; 1062 : C
; 1063 : C
; 1064 : 	IF(ISYMM) THEN
; 1065 : 		CALL OVERLP(M,DFLAG,RX1(L),RY1(L),RZ1(L),AA1(L),

	cmp	edi, esi
	jg	SHORT $L100146
	cmp	edi, 1
	jge	SHORT $L100147
$L100146:
	call	___FHranger@0
$L100147:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100148
	cmp	edi, 1
	jge	SHORT $L100149
$L100148:
	call	___FHranger@0
$L100149:
	mov	eax, DWORD PTR _B12$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100150
	cmp	edi, 1
	jge	SHORT $L100151
$L100150:
	call	___FHranger@0
$L100151:
	mov	eax, DWORD PTR _AA1$[esp+1404]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100152
	cmp	edi, 1
	jge	SHORT $L100153
$L100152:
	call	___FHranger@0
$L100153:
	mov	eax, DWORD PTR _RZ1$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100154
	cmp	edi, 1
	jge	SHORT $L100155
$L100154:
	call	___FHranger@0
$L100155:
	mov	eax, DWORD PTR _RY1$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100156
	cmp	edi, 1
	jge	SHORT $L100157
$L100156:
	call	___FHranger@0
$L100157:
	mov	eax, DWORD PTR _RX1$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100158
	cmp	edi, 1
	jge	SHORT $L100159
$L100158:
	call	___FHranger@0
$L100159:
	mov	eax, DWORD PTR _AA2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100160
	cmp	edi, 1
	jge	SHORT $L100161
$L100160:
	call	___FHranger@0
$L100161:
	mov	eax, DWORD PTR _RZ2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100162
	cmp	edi, 1
	jge	SHORT $L100163
$L100162:
	call	___FHranger@0
$L100163:
	mov	eax, DWORD PTR _RY2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100164
	cmp	edi, 1
	jge	SHORT $L100165
$L100164:
	call	___FHranger@0
$L100165:
	mov	eax, DWORD PTR _RX2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100166
	cmp	edi, 1
	jge	SHORT $L100167
$L100166:
	call	___FHranger@0
$L100167:
	mov	eax, DWORD PTR _AA2$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100168
	cmp	edi, 1
	jge	SHORT $L100169
$L100168:
	call	___FHranger@0
$L100169:
	mov	eax, DWORD PTR _RZ2$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100170
	cmp	edi, 1
	jge	SHORT $L100171
$L100170:
	call	___FHranger@0
$L100171:
	mov	eax, DWORD PTR _RY2$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100172
	cmp	edi, 1
	jge	SHORT $L100173
$L100172:
	call	___FHranger@0
$L100173:
	mov	eax, DWORD PTR _RX2$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100174
	cmp	edi, 1
	jge	SHORT $L100175
$L100174:
	call	___FHranger@0
$L100175:
	mov	eax, DWORD PTR _AA1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100176
	cmp	edi, 1
	jge	SHORT $L100177
$L100176:
	call	___FHranger@0
$L100177:
	mov	eax, DWORD PTR _RZ1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100178
	cmp	edi, 1
	jge	SHORT $L100179
$L100178:
	call	___FHranger@0
$L100179:
	mov	eax, DWORD PTR _RY1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100180
	cmp	edi, 1
	jge	SHORT $L100181
$L100180:
	call	___FHranger@0
$L100181:
	mov	eax, DWORD PTR _RX1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1464]
	sub	edi, 8
	push	edi
	push	eax
	mov	eax, DWORD PTR _M$[esp+1472]
	push	eax
	call	_OVERLP@88
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 1074	; 00000432H
	mov	esi, DWORD PTR [eax]

; 1066 :      1            RX2(K),RY2(K),RZ2(K),AA2(K),                      
; 1067 :      2            RX2(L),RY2(L),RZ2(L),AA2(L),
; 1068 :      3            RX1(K),RY1(K),RZ1(K),AA1(K),
; 1069 :      4            B12(L),B12(K),OV2,KLM)
; 1070 : C
; 1071 : C
; 1072 : C     ...Nuclear Symmetrization Overlap Integral: P_[H1,H2]
; 1073 : C
; 1074 :               CALL OVERLP(M,DFLAG,RX1(L),RY1(L),RZ1(L),AA1(L),

	cmp	edi, esi
	jg	SHORT $L100182
	cmp	edi, 1
	jge	SHORT $L100183
$L100182:
	call	___FHranger@0
$L100183:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98421[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100184
	cmp	edi, 1
	jge	SHORT $L100185
$L100184:
	call	___FHranger@0
$L100185:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98422[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100186
	cmp	edi, 1
	jge	SHORT $L100187
$L100186:
	call	___FHranger@0
$L100187:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98423[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100188
	cmp	edi, 1
	jge	SHORT $L100189
$L100188:
	call	___FHranger@0
$L100189:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98424[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100190
	cmp	edi, 1
	jge	SHORT $L100191
$L100190:
	call	___FHranger@0
$L100191:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98425[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100192
	cmp	edi, 1
	jge	SHORT $L100193
$L100192:
	call	___FHranger@0
$L100193:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_OV3
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98426[esp+1400]
	cmp	edi, esi
	jg	SHORT $L100194
	cmp	edi, 1
	jge	SHORT $L100195
$L100194:
	call	___FHranger@0
$L100195:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100196
	cmp	edi, 1
	jge	SHORT $L100197
$L100196:
	call	___FHranger@0
$L100197:
	mov	eax, DWORD PTR _B12$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100198
	cmp	edi, 1
	jge	SHORT $L100199
$L100198:
	call	___FHranger@0
$L100199:
	mov	eax, DWORD PTR _AA2$[esp+1404]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98421[esp+1408]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98422[esp+1416]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98423[esp+1420]
	push	eax
	jg	SHORT $L100200
	cmp	edi, 1
	jge	SHORT $L100201
$L100200:
	call	___FHranger@0
$L100201:
	mov	eax, DWORD PTR _AA2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100202
	cmp	edi, 1
	jge	SHORT $L100203
$L100202:
	call	___FHranger@0
$L100203:
	mov	eax, DWORD PTR _RZ2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100204
	cmp	edi, 1
	jge	SHORT $L100205
$L100204:
	call	___FHranger@0
$L100205:
	mov	eax, DWORD PTR _RY2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100206
	cmp	edi, 1
	jge	SHORT $L100207
$L100206:
	call	___FHranger@0
$L100207:
	mov	eax, DWORD PTR _RX2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100208
	cmp	edi, 1
	jge	SHORT $L100209
$L100208:
	call	___FHranger@0
$L100209:
	mov	eax, DWORD PTR _AA1$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98424[esp+1440]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98425[esp+1448]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98426[esp+1452]
	push	eax
	jg	SHORT $L100210
	cmp	edi, 1
	jge	SHORT $L100211
$L100210:
	call	___FHranger@0
$L100211:
	mov	eax, DWORD PTR _AA1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100212
	cmp	edi, 1
	jge	SHORT $L100213
$L100212:
	call	___FHranger@0
$L100213:
	mov	eax, DWORD PTR _RZ1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100214
	cmp	edi, 1
	jge	SHORT $L100215
$L100214:
	call	___FHranger@0
$L100215:
	mov	eax, DWORD PTR _RY1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100216
	cmp	edi, 1
	jge	SHORT $L100217
$L100216:
	call	___FHranger@0
$L100217:
	mov	eax, DWORD PTR _RX1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1464]
	sub	edi, 8
	push	edi
	push	eax
	mov	eax, DWORD PTR _M$[esp+1472]
	push	eax
	call	_OVERLP@88
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 1084	; 0000043cH
	mov	esi, DWORD PTR [eax]

; 1075 :      1            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
; 1076 :      2            RX2(L),RY2(L),RZ2(L),AA2(L),
; 1077 :      3            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
; 1078 :      4            B12(L),B12(K),OV3,KLM)
; 1079 : C
; 1080 : C
; 1081 : C     ...Electro/Nuclear Symmetrization Overlap Integral:
; 1082 : C         P_[elec.1,elec.2]*P_[H1,H2]
; 1083 : C
; 1084 :               CALL OVERLP(M,DFLAG,RX1(L),RY1(L),RZ1(L),AA1(L),

	cmp	edi, esi
	jg	SHORT $L100218
	cmp	edi, 1
	jge	SHORT $L100219
$L100218:
	call	___FHranger@0
$L100219:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98427[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100220
	cmp	edi, 1
	jge	SHORT $L100221
$L100220:
	call	___FHranger@0
$L100221:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98428[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100222
	cmp	edi, 1
	jge	SHORT $L100223
$L100222:
	call	___FHranger@0
$L100223:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98429[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100224
	cmp	edi, 1
	jge	SHORT $L100225
$L100224:
	call	___FHranger@0
$L100225:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98430[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100226
	cmp	edi, 1
	jge	SHORT $L100227
$L100226:
	call	___FHranger@0
$L100227:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98431[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100228
	cmp	edi, 1
	jge	SHORT $L100229
$L100228:
	call	___FHranger@0
$L100229:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_OV4
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98432[esp+1400]
	cmp	edi, esi
	jg	SHORT $L100230
	cmp	edi, 1
	jge	SHORT $L100231
$L100230:
	call	___FHranger@0
$L100231:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100232
	cmp	edi, 1
	jge	SHORT $L100233
$L100232:
	call	___FHranger@0
$L100233:
	mov	eax, DWORD PTR _B12$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100234
	cmp	edi, 1
	jge	SHORT $L100235
$L100234:
	call	___FHranger@0
$L100235:
	mov	eax, DWORD PTR _AA1$[esp+1404]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98427[esp+1408]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98428[esp+1416]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98429[esp+1420]
	push	eax
	jg	SHORT $L100236
	cmp	edi, 1
	jge	SHORT $L100237
$L100236:
	call	___FHranger@0
$L100237:
	mov	eax, DWORD PTR _AA2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100238
	cmp	edi, 1
	jge	SHORT $L100239
$L100238:
	call	___FHranger@0
$L100239:
	mov	eax, DWORD PTR _RZ2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100240
	cmp	edi, 1
	jge	SHORT $L100241
$L100240:
	call	___FHranger@0
$L100241:
	mov	eax, DWORD PTR _RY2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100242
	cmp	edi, 1
	jge	SHORT $L100243
$L100242:
	call	___FHranger@0
$L100243:
	mov	eax, DWORD PTR _RX2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100244
	cmp	edi, 1
	jge	SHORT $L100245
$L100244:
	call	___FHranger@0
$L100245:
	mov	eax, DWORD PTR _AA2$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98430[esp+1440]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98431[esp+1448]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98432[esp+1452]
	push	eax
	jg	SHORT $L100246
	cmp	edi, 1
	jge	SHORT $L100247
$L100246:
	call	___FHranger@0
$L100247:
	mov	eax, DWORD PTR _AA1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100248
	cmp	edi, 1
	jge	SHORT $L100249
$L100248:
	call	___FHranger@0
$L100249:
	mov	eax, DWORD PTR _RZ1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100250
	cmp	edi, 1
	jge	SHORT $L100251
$L100250:
	call	___FHranger@0
$L100251:
	mov	eax, DWORD PTR _RY1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100252
	cmp	edi, 1
	jge	SHORT $L100253
$L100252:
	call	___FHranger@0
$L100253:
	mov	eax, DWORD PTR _RX1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1464]
	sub	edi, 8
	push	edi
	push	eax
	mov	eax, DWORD PTR _M$[esp+1472]
	push	eax
	call	_OVERLP@88
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	edi, DWORD PTR _OVER
	mov	DWORD PTR $T98573[esp+1404], 1093	; 00000445H
	mov	esi, DWORD PTR [eax]

; 1085 :      1            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
; 1086 :      2            RX2(L),RY2(L),RZ2(L),AA2(L),
; 1087 :      3            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
; 1088 :      4            B12(L),B12(K),OV4,KLM)
; 1089 : C              
; 1090 : 	END IF
; 1091 : C     -----------------------------------------------------
; 1092 : 	IF(ISYMM) THEN
; 1093 :               S(L,K)=OV1+OV2+OV3+OV4

	cmp	edi, esi
	jg	SHORT $L100254
	cmp	edi, 1
	jge	SHORT $L100255
$L100254:
	call	___FHranger@0
$L100255:
	mov	ebx, DWORD PTR _OVER+4
	cmp	ebx, esi
	jg	SHORT $L100256
	cmp	ebx, 1
	jge	SHORT $L100257
$L100256:
	call	___FHranger@0
$L100257:
	fld	QWORD PTR _OV3
	fadd	QWORD PTR _OV1
	mov	ecx, DWORD PTR _7S$[esp+1428]
	mov	eax, DWORD PTR _S$[esp+1388]
	mov	DWORD PTR $T98573[esp+1404], 1116	; 0000045cH
	fadd	QWORD PTR _OV2
	imul	ecx, ebx
	fadd	QWORD PTR _OV4
	add	ecx, DWORD PTR _7S$[esp+1400]
	push	OFFSET FLAT:_T1
	add	ecx, edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	fstp	QWORD PTR [eax+ecx*8]

; 1094 : 		ELSE
; 1095 :               S(L,K)=OV1
; 1096 : 	END IF
; 1097 : 	IF(IPERM.EQ.1) THEN
; 1098 : 		S(L,K)=OV1
; 1099 : 	END IF
; 1100 : 	IF(IPERM.EQ.2) THEN
; 1101 : 		S(L,K)=OV2
; 1102 : 	END IF
; 1103 : 	IF(IPERM.EQ.3) THEN
; 1104 : 		S(L,K)=OV3
; 1105 : 	END IF
; 1106 : 	IF(IPERM.EQ.4) THEN
; 1107 : 		S(L,K)=OV4
; 1108 : 	END IF
; 1109 : C
; 1110 : C	IF (K.EQ.1) THEN
; 1111 : C		WRITE(6,*)'S(',L,',1)=',S(L,K)
; 1112 : C	END IF
; 1113 : C     -----------------------------------------------------
; 1114 : C     --- Calculate KINETIC ENERGY OF ELECTRON 1 ---
; 1115 : C     -----------------------------------------------------
; 1116 :               CALL TK(DFLAG,OV1,

	jg	SHORT $L100258
	cmp	edi, 1
	jge	SHORT $L100259
$L100258:
	call	___FHranger@0
$L100259:
	mov	eax, DWORD PTR _B12$[esp+1392]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100260
	cmp	edi, 1
	jge	SHORT $L100261
$L100260:
	call	___FHranger@0
$L100261:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100262
	cmp	edi, 1
	jge	SHORT $L100263
$L100262:
	call	___FHranger@0
$L100263:
	mov	eax, DWORD PTR _AA2$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100264
	cmp	edi, 1
	jge	SHORT $L100265
$L100264:
	call	___FHranger@0
$L100265:
	mov	eax, DWORD PTR _RZ2$[esp+1404]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100266
	cmp	edi, 1
	jge	SHORT $L100267
$L100266:
	call	___FHranger@0
$L100267:
	mov	eax, DWORD PTR _RY2$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100268
	cmp	edi, 1
	jge	SHORT $L100269
$L100268:
	call	___FHranger@0
$L100269:
	mov	eax, DWORD PTR _RX2$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100270
	cmp	edi, 1
	jge	SHORT $L100271
$L100270:
	call	___FHranger@0
$L100271:
	mov	eax, DWORD PTR _AA2$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100272
	cmp	edi, 1
	jge	SHORT $L100273
$L100272:
	call	___FHranger@0
$L100273:
	mov	eax, DWORD PTR _RZ2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100274
	cmp	edi, 1
	jge	SHORT $L100275
$L100274:
	call	___FHranger@0
$L100275:
	mov	eax, DWORD PTR _RY2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100276
	cmp	edi, 1
	jge	SHORT $L100277
$L100276:
	call	___FHranger@0
$L100277:
	mov	eax, DWORD PTR _RX2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100278
	cmp	edi, 1
	jge	SHORT $L100279
$L100278:
	call	___FHranger@0
$L100279:
	mov	eax, DWORD PTR _AA1$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100280
	cmp	edi, 1
	jge	SHORT $L100281
$L100280:
	call	___FHranger@0
$L100281:
	mov	eax, DWORD PTR _RZ1$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100282
	cmp	edi, 1
	jge	SHORT $L100283
$L100282:
	call	___FHranger@0
$L100283:
	mov	eax, DWORD PTR _RY1$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100284
	cmp	edi, 1
	jge	SHORT $L100285
$L100284:
	call	___FHranger@0
$L100285:
	mov	eax, DWORD PTR _RX1$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100286
	cmp	edi, 1
	jge	SHORT $L100287
$L100286:
	call	___FHranger@0
$L100287:
	mov	eax, DWORD PTR _AA1$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100288
	cmp	edi, 1
	jge	SHORT $L100289
$L100288:
	call	___FHranger@0
$L100289:
	mov	eax, DWORD PTR _RZ1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100290
	cmp	edi, 1
	jge	SHORT $L100291
$L100290:
	call	___FHranger@0
$L100291:
	mov	eax, DWORD PTR _RY1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100292
	cmp	edi, 1
	jge	SHORT $L100293
$L100292:
	call	___FHranger@0
$L100293:
	mov	eax, DWORD PTR _RX1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1460]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV1
	push	eax
	call	_TK@84
	mov	eax, DWORD PTR _M$[esp+1388]
	push	OFFSET FLAT:_T2
	mov	DWORD PTR $T98573[esp+1408], 1125	; 00000465H
	mov	esi, DWORD PTR [eax]
	mov	edi, DWORD PTR _OVER+4

; 1117 :      1            RX1(L),RY1(L),RZ1(L),AA1(L),
; 1118 :      2            RX1(K),RY1(K),RZ1(K),AA1(K),                      
; 1119 :      3            RX2(L),RY2(L),RZ2(L),AA2(L),
; 1120 :      4            RX2(K),RY2(K),RZ2(K),AA2(K),
; 1121 :      5            B12(L),B12(K),T1)  
; 1122 : C
; 1123 : C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
; 1124 : 	IF (ISYMM) THEN
; 1125 :               CALL TK(DFLAG,OV2,

	cmp	edi, esi
	jg	SHORT $L100294
	cmp	edi, 1
	jge	SHORT $L100295
$L100294:
	call	___FHranger@0
$L100295:
	mov	eax, DWORD PTR _B12$[esp+1392]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100296
	cmp	edi, 1
	jge	SHORT $L100297
$L100296:
	call	___FHranger@0
$L100297:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100298
	cmp	edi, 1
	jge	SHORT $L100299
$L100298:
	call	___FHranger@0
$L100299:
	mov	eax, DWORD PTR _AA1$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100300
	cmp	edi, 1
	jge	SHORT $L100301
$L100300:
	call	___FHranger@0
$L100301:
	mov	eax, DWORD PTR _RZ1$[esp+1404]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100302
	cmp	edi, 1
	jge	SHORT $L100303
$L100302:
	call	___FHranger@0
$L100303:
	mov	eax, DWORD PTR _RY1$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100304
	cmp	edi, 1
	jge	SHORT $L100305
$L100304:
	call	___FHranger@0
$L100305:
	mov	eax, DWORD PTR _RX1$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100306
	cmp	edi, 1
	jge	SHORT $L100307
$L100306:
	call	___FHranger@0
$L100307:
	mov	eax, DWORD PTR _AA2$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100308
	cmp	edi, 1
	jge	SHORT $L100309
$L100308:
	call	___FHranger@0
$L100309:
	mov	eax, DWORD PTR _RZ2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100310
	cmp	edi, 1
	jge	SHORT $L100311
$L100310:
	call	___FHranger@0
$L100311:
	mov	eax, DWORD PTR _RY2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100312
	cmp	edi, 1
	jge	SHORT $L100313
$L100312:
	call	___FHranger@0
$L100313:
	mov	eax, DWORD PTR _RX2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100314
	cmp	edi, 1
	jge	SHORT $L100315
$L100314:
	call	___FHranger@0
$L100315:
	mov	eax, DWORD PTR _AA2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100316
	cmp	edi, 1
	jge	SHORT $L100317
$L100316:
	call	___FHranger@0
$L100317:
	mov	eax, DWORD PTR _RZ2$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100318
	cmp	edi, 1
	jge	SHORT $L100319
$L100318:
	call	___FHranger@0
$L100319:
	mov	eax, DWORD PTR _RY2$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100320
	cmp	edi, 1
	jge	SHORT $L100321
$L100320:
	call	___FHranger@0
$L100321:
	mov	eax, DWORD PTR _RX2$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100322
	cmp	edi, 1
	jge	SHORT $L100323
$L100322:
	call	___FHranger@0
$L100323:
	mov	eax, DWORD PTR _AA1$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100324
	cmp	edi, 1
	jge	SHORT $L100325
$L100324:
	call	___FHranger@0
$L100325:
	mov	eax, DWORD PTR _RZ1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100326
	cmp	edi, 1
	jge	SHORT $L100327
$L100326:
	call	___FHranger@0
$L100327:
	mov	eax, DWORD PTR _RY1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100328
	cmp	edi, 1
	jge	SHORT $L100329
$L100328:
	call	___FHranger@0
$L100329:
	mov	eax, DWORD PTR _RX1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1460]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV2
	push	eax
	call	_TK@84
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 1132	; 0000046cH
	mov	esi, DWORD PTR [eax]

; 1126 :      1            RX1(L),RY1(L),RZ1(L),AA1(L),
; 1127 :      2            RX2(K),RY2(K),RZ2(K),AA2(K),                      
; 1128 :      3            RX2(L),RY2(L),RZ2(L),AA2(L),
; 1129 :      4            RX1(K),RY1(K),RZ1(K),AA1(K),
; 1130 :      5            B12(L),B12(K),T2)  
; 1131 : C     ...Nuclear Symmetrization Integral: P_[H1,H2]
; 1132 :               CALL TK(DFLAG,OV3,

	cmp	edi, esi
	jg	SHORT $L100330
	cmp	edi, 1
	jge	SHORT $L100331
$L100330:
	call	___FHranger@0
$L100331:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98433[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100332
	cmp	edi, 1
	jge	SHORT $L100333
$L100332:
	call	___FHranger@0
$L100333:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98434[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100334
	cmp	edi, 1
	jge	SHORT $L100335
$L100334:
	call	___FHranger@0
$L100335:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98435[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100336
	cmp	edi, 1
	jge	SHORT $L100337
$L100336:
	call	___FHranger@0
$L100337:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98436[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100338
	cmp	edi, 1
	jge	SHORT $L100339
$L100338:
	call	___FHranger@0
$L100339:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98437[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100340
	cmp	edi, 1
	jge	SHORT $L100341
$L100340:
	call	___FHranger@0
$L100341:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	push	OFFSET FLAT:_T3
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98438[esp+1396]
	cmp	edi, esi
	jg	SHORT $L100342
	cmp	edi, 1
	jge	SHORT $L100343
$L100342:
	call	___FHranger@0
$L100343:
	mov	eax, DWORD PTR _B12$[esp+1392]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100344
	cmp	edi, 1
	jge	SHORT $L100345
$L100344:
	call	___FHranger@0
$L100345:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100346
	cmp	edi, 1
	jge	SHORT $L100347
$L100346:
	call	___FHranger@0
$L100347:
	mov	eax, DWORD PTR _AA2$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98433[esp+1404]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98434[esp+1412]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98435[esp+1416]
	push	eax
	jg	SHORT $L100348
	cmp	edi, 1
	jge	SHORT $L100349
$L100348:
	call	___FHranger@0
$L100349:
	mov	eax, DWORD PTR _AA2$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100350
	cmp	edi, 1
	jge	SHORT $L100351
$L100350:
	call	___FHranger@0
$L100351:
	mov	eax, DWORD PTR _RZ2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100352
	cmp	edi, 1
	jge	SHORT $L100353
$L100352:
	call	___FHranger@0
$L100353:
	mov	eax, DWORD PTR _RY2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100354
	cmp	edi, 1
	jge	SHORT $L100355
$L100354:
	call	___FHranger@0
$L100355:
	mov	eax, DWORD PTR _RX2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100356
	cmp	edi, 1
	jge	SHORT $L100357
$L100356:
	call	___FHranger@0
$L100357:
	mov	eax, DWORD PTR _AA1$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98436[esp+1436]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98437[esp+1444]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98438[esp+1448]
	push	eax
	jg	SHORT $L100358
	cmp	edi, 1
	jge	SHORT $L100359
$L100358:
	call	___FHranger@0
$L100359:
	mov	eax, DWORD PTR _AA1$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100360
	cmp	edi, 1
	jge	SHORT $L100361
$L100360:
	call	___FHranger@0
$L100361:
	mov	eax, DWORD PTR _RZ1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100362
	cmp	edi, 1
	jge	SHORT $L100363
$L100362:
	call	___FHranger@0
$L100363:
	mov	eax, DWORD PTR _RY1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100364
	cmp	edi, 1
	jge	SHORT $L100365
$L100364:
	call	___FHranger@0
$L100365:
	mov	eax, DWORD PTR _RX1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1460]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV3
	push	eax
	call	_TK@84
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 1140	; 00000474H
	mov	esi, DWORD PTR [eax]

; 1133 :      1            RX1(L),RY1(L),RZ1(L),AA1(L),
; 1134 :      2            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
; 1135 :      3            RX2(L),RY2(L),RZ2(L),AA2(L),
; 1136 :      4            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
; 1137 :      5            B12(L),B12(K),T3)
; 1138 : C     ...Electro/Nuclear Symmetrization Integral:
; 1139 : C         P_[elec.1,elec.2]*P_[H1,H2]
; 1140 :               CALL TK(DFLAG,OV4,

	cmp	edi, esi
	jg	SHORT $L100366
	cmp	edi, 1
	jge	SHORT $L100367
$L100366:
	call	___FHranger@0
$L100367:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98439[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100368
	cmp	edi, 1
	jge	SHORT $L100369
$L100368:
	call	___FHranger@0
$L100369:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98440[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100370
	cmp	edi, 1
	jge	SHORT $L100371
$L100370:
	call	___FHranger@0
$L100371:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98441[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100372
	cmp	edi, 1
	jge	SHORT $L100373
$L100372:
	call	___FHranger@0
$L100373:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98442[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100374
	cmp	edi, 1
	jge	SHORT $L100375
$L100374:
	call	___FHranger@0
$L100375:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98443[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100376
	cmp	edi, 1
	jge	SHORT $L100377
$L100376:
	call	___FHranger@0
$L100377:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	push	OFFSET FLAT:_T4
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98444[esp+1396]
	cmp	edi, esi
	jg	SHORT $L100378
	cmp	edi, 1
	jge	SHORT $L100379
$L100378:
	call	___FHranger@0
$L100379:
	mov	eax, DWORD PTR _B12$[esp+1392]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100380
	cmp	edi, 1
	jge	SHORT $L100381
$L100380:
	call	___FHranger@0
$L100381:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100382
	cmp	edi, 1
	jge	SHORT $L100383
$L100382:
	call	___FHranger@0
$L100383:
	mov	eax, DWORD PTR _AA1$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98439[esp+1404]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98440[esp+1412]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98441[esp+1416]
	push	eax
	jg	SHORT $L100384
	cmp	edi, 1
	jge	SHORT $L100385
$L100384:
	call	___FHranger@0
$L100385:
	mov	eax, DWORD PTR _AA2$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100386
	cmp	edi, 1
	jge	SHORT $L100387
$L100386:
	call	___FHranger@0
$L100387:
	mov	eax, DWORD PTR _RZ2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100388
	cmp	edi, 1
	jge	SHORT $L100389
$L100388:
	call	___FHranger@0
$L100389:
	mov	eax, DWORD PTR _RY2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100390
	cmp	edi, 1
	jge	SHORT $L100391
$L100390:
	call	___FHranger@0
$L100391:
	mov	eax, DWORD PTR _RX2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100392
	cmp	edi, 1
	jge	SHORT $L100393
$L100392:
	call	___FHranger@0
$L100393:
	mov	eax, DWORD PTR _AA2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98442[esp+1436]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98443[esp+1444]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98444[esp+1448]
	push	eax
	jg	SHORT $L100394
	cmp	edi, 1
	jge	SHORT $L100395
$L100394:
	call	___FHranger@0
$L100395:
	mov	eax, DWORD PTR _AA1$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100396
	cmp	edi, 1
	jge	SHORT $L100397
$L100396:
	call	___FHranger@0
$L100397:
	mov	eax, DWORD PTR _RZ1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100398
	cmp	edi, 1
	jge	SHORT $L100399
$L100398:
	call	___FHranger@0
$L100399:
	mov	eax, DWORD PTR _RY1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100400
	cmp	edi, 1
	jge	SHORT $L100401
$L100400:
	call	___FHranger@0
$L100401:
	mov	eax, DWORD PTR _RX1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1460]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV4
	push	eax
	call	_TK@84
	mov	eax, DWORD PTR _M$[esp+1388]
	push	OFFSET FLAT:_T5
	mov	DWORD PTR $T98573[esp+1408], 1150	; 0000047eH
	mov	esi, DWORD PTR [eax]
	mov	edi, DWORD PTR _OVER+4

; 1141 :      1            RX1(L),RY1(L),RZ1(L),AA1(L),
; 1142 :      2            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
; 1143 :      3            RX2(L),RY2(L),RZ2(L),AA2(L),
; 1144 :      4            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
; 1145 :      5            B12(L),B12(K),T4)
; 1146 : 	END IF
; 1147 : C     -----------------------------------------------------
; 1148 : C     --- Calculate KINETIC ENERGY OF ELECTRON 2 ---
; 1149 : C     -----------------------------------------------------
; 1150 :              CALL TK(DFLAG,OV1,

	cmp	edi, esi
	jg	SHORT $L100402
	cmp	edi, 1
	jge	SHORT $L100403
$L100402:
	call	___FHranger@0
$L100403:
	mov	eax, DWORD PTR _B12$[esp+1392]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100404
	cmp	edi, 1
	jge	SHORT $L100405
$L100404:
	call	___FHranger@0
$L100405:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100406
	cmp	edi, 1
	jge	SHORT $L100407
$L100406:
	call	___FHranger@0
$L100407:
	mov	eax, DWORD PTR _AA1$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100408
	cmp	edi, 1
	jge	SHORT $L100409
$L100408:
	call	___FHranger@0
$L100409:
	mov	eax, DWORD PTR _RZ1$[esp+1404]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100410
	cmp	edi, 1
	jge	SHORT $L100411
$L100410:
	call	___FHranger@0
$L100411:
	mov	eax, DWORD PTR _RY1$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100412
	cmp	edi, 1
	jge	SHORT $L100413
$L100412:
	call	___FHranger@0
$L100413:
	mov	eax, DWORD PTR _RX1$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100414
	cmp	edi, 1
	jge	SHORT $L100415
$L100414:
	call	___FHranger@0
$L100415:
	mov	eax, DWORD PTR _AA1$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100416
	cmp	edi, 1
	jge	SHORT $L100417
$L100416:
	call	___FHranger@0
$L100417:
	mov	eax, DWORD PTR _RZ1$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100418
	cmp	edi, 1
	jge	SHORT $L100419
$L100418:
	call	___FHranger@0
$L100419:
	mov	eax, DWORD PTR _RY1$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100420
	cmp	edi, 1
	jge	SHORT $L100421
$L100420:
	call	___FHranger@0
$L100421:
	mov	eax, DWORD PTR _RX1$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100422
	cmp	edi, 1
	jge	SHORT $L100423
$L100422:
	call	___FHranger@0
$L100423:
	mov	eax, DWORD PTR _AA2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100424
	cmp	edi, 1
	jge	SHORT $L100425
$L100424:
	call	___FHranger@0
$L100425:
	mov	eax, DWORD PTR _RZ2$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100426
	cmp	edi, 1
	jge	SHORT $L100427
$L100426:
	call	___FHranger@0
$L100427:
	mov	eax, DWORD PTR _RY2$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100428
	cmp	edi, 1
	jge	SHORT $L100429
$L100428:
	call	___FHranger@0
$L100429:
	mov	eax, DWORD PTR _RX2$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100430
	cmp	edi, 1
	jge	SHORT $L100431
$L100430:
	call	___FHranger@0
$L100431:
	mov	eax, DWORD PTR _AA2$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100432
	cmp	edi, 1
	jge	SHORT $L100433
$L100432:
	call	___FHranger@0
$L100433:
	mov	eax, DWORD PTR _RZ2$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100434
	cmp	edi, 1
	jge	SHORT $L100435
$L100434:
	call	___FHranger@0
$L100435:
	mov	eax, DWORD PTR _RY2$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100436
	cmp	edi, 1
	jge	SHORT $L100437
$L100436:
	call	___FHranger@0
$L100437:
	mov	eax, DWORD PTR _RX2$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1460]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV1
	push	eax
	call	_TK@84
	mov	eax, DWORD PTR _M$[esp+1388]
	push	OFFSET FLAT:_T6
	mov	DWORD PTR $T98573[esp+1408], 1158	; 00000486H
	mov	esi, DWORD PTR [eax]
	mov	edi, DWORD PTR _OVER+4

; 1151 :      1            RX2(L),RY2(L),RZ2(L),AA2(L),
; 1152 :      2            RX2(K),RY2(K),RZ2(K),AA2(K),                      
; 1153 :      3            RX1(L),RY1(L),RZ1(L),AA1(L),
; 1154 :      4            RX1(K),RY1(K),RZ1(K),AA1(K),
; 1155 :      5            B12(L),B12(K),T5)  
; 1156 : C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
; 1157 : 	IF (ISYMM) THEN
; 1158 :               CALL TK(DFLAG,OV2,

	cmp	edi, esi
	jg	SHORT $L100438
	cmp	edi, 1
	jge	SHORT $L100439
$L100438:
	call	___FHranger@0
$L100439:
	mov	eax, DWORD PTR _B12$[esp+1392]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100440
	cmp	edi, 1
	jge	SHORT $L100441
$L100440:
	call	___FHranger@0
$L100441:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100442
	cmp	edi, 1
	jge	SHORT $L100443
$L100442:
	call	___FHranger@0
$L100443:
	mov	eax, DWORD PTR _AA2$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100444
	cmp	edi, 1
	jge	SHORT $L100445
$L100444:
	call	___FHranger@0
$L100445:
	mov	eax, DWORD PTR _RZ2$[esp+1404]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100446
	cmp	edi, 1
	jge	SHORT $L100447
$L100446:
	call	___FHranger@0
$L100447:
	mov	eax, DWORD PTR _RY2$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100448
	cmp	edi, 1
	jge	SHORT $L100449
$L100448:
	call	___FHranger@0
$L100449:
	mov	eax, DWORD PTR _RX2$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100450
	cmp	edi, 1
	jge	SHORT $L100451
$L100450:
	call	___FHranger@0
$L100451:
	mov	eax, DWORD PTR _AA1$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100452
	cmp	edi, 1
	jge	SHORT $L100453
$L100452:
	call	___FHranger@0
$L100453:
	mov	eax, DWORD PTR _RZ1$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100454
	cmp	edi, 1
	jge	SHORT $L100455
$L100454:
	call	___FHranger@0
$L100455:
	mov	eax, DWORD PTR _RY1$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100456
	cmp	edi, 1
	jge	SHORT $L100457
$L100456:
	call	___FHranger@0
$L100457:
	mov	eax, DWORD PTR _RX1$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100458
	cmp	edi, 1
	jge	SHORT $L100459
$L100458:
	call	___FHranger@0
$L100459:
	mov	eax, DWORD PTR _AA1$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100460
	cmp	edi, 1
	jge	SHORT $L100461
$L100460:
	call	___FHranger@0
$L100461:
	mov	eax, DWORD PTR _RZ1$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100462
	cmp	edi, 1
	jge	SHORT $L100463
$L100462:
	call	___FHranger@0
$L100463:
	mov	eax, DWORD PTR _RY1$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100464
	cmp	edi, 1
	jge	SHORT $L100465
$L100464:
	call	___FHranger@0
$L100465:
	mov	eax, DWORD PTR _RX1$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100466
	cmp	edi, 1
	jge	SHORT $L100467
$L100466:
	call	___FHranger@0
$L100467:
	mov	eax, DWORD PTR _AA2$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100468
	cmp	edi, 1
	jge	SHORT $L100469
$L100468:
	call	___FHranger@0
$L100469:
	mov	eax, DWORD PTR _RZ2$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100470
	cmp	edi, 1
	jge	SHORT $L100471
$L100470:
	call	___FHranger@0
$L100471:
	mov	eax, DWORD PTR _RY2$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100472
	cmp	edi, 1
	jge	SHORT $L100473
$L100472:
	call	___FHranger@0
$L100473:
	mov	eax, DWORD PTR _RX2$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1460]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV2
	push	eax
	call	_TK@84
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 1165	; 0000048dH
	mov	esi, DWORD PTR [eax]

; 1159 :      1            RX2(L),RY2(L),RZ2(L),AA2(L),
; 1160 :      2            RX1(K),RY1(K),RZ1(K),AA1(K),                      
; 1161 :      3            RX1(L),RY1(L),RZ1(L),AA1(L),
; 1162 :      4            RX2(K),RY2(K),RZ2(K),AA2(K),
; 1163 :      5            B12(L),B12(K),T6) 
; 1164 : C     ...Nuclear Symmetrization Integral: P_[H1,H2]
; 1165 :               CALL TK(DFLAG,OV3,

	cmp	edi, esi
	jg	SHORT $L100474
	cmp	edi, 1
	jge	SHORT $L100475
$L100474:
	call	___FHranger@0
$L100475:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98445[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100476
	cmp	edi, 1
	jge	SHORT $L100477
$L100476:
	call	___FHranger@0
$L100477:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98446[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100478
	cmp	edi, 1
	jge	SHORT $L100479
$L100478:
	call	___FHranger@0
$L100479:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98447[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100480
	cmp	edi, 1
	jge	SHORT $L100481
$L100480:
	call	___FHranger@0
$L100481:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98448[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100482
	cmp	edi, 1
	jge	SHORT $L100483
$L100482:
	call	___FHranger@0
$L100483:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98449[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100484
	cmp	edi, 1
	jge	SHORT $L100485
$L100484:
	call	___FHranger@0
$L100485:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	push	OFFSET FLAT:_T7
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98450[esp+1396]
	cmp	edi, esi
	jg	SHORT $L100486
	cmp	edi, 1
	jge	SHORT $L100487
$L100486:
	call	___FHranger@0
$L100487:
	mov	eax, DWORD PTR _B12$[esp+1392]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100488
	cmp	edi, 1
	jge	SHORT $L100489
$L100488:
	call	___FHranger@0
$L100489:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100490
	cmp	edi, 1
	jge	SHORT $L100491
$L100490:
	call	___FHranger@0
$L100491:
	mov	eax, DWORD PTR _AA1$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98445[esp+1404]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98446[esp+1412]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98447[esp+1416]
	push	eax
	jg	SHORT $L100492
	cmp	edi, 1
	jge	SHORT $L100493
$L100492:
	call	___FHranger@0
$L100493:
	mov	eax, DWORD PTR _AA1$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100494
	cmp	edi, 1
	jge	SHORT $L100495
$L100494:
	call	___FHranger@0
$L100495:
	mov	eax, DWORD PTR _RZ1$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100496
	cmp	edi, 1
	jge	SHORT $L100497
$L100496:
	call	___FHranger@0
$L100497:
	mov	eax, DWORD PTR _RY1$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100498
	cmp	edi, 1
	jge	SHORT $L100499
$L100498:
	call	___FHranger@0
$L100499:
	mov	eax, DWORD PTR _RX1$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100500
	cmp	edi, 1
	jge	SHORT $L100501
$L100500:
	call	___FHranger@0
$L100501:
	mov	eax, DWORD PTR _AA2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98448[esp+1436]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98449[esp+1444]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98450[esp+1448]
	push	eax
	jg	SHORT $L100502
	cmp	edi, 1
	jge	SHORT $L100503
$L100502:
	call	___FHranger@0
$L100503:
	mov	eax, DWORD PTR _AA2$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100504
	cmp	edi, 1
	jge	SHORT $L100505
$L100504:
	call	___FHranger@0
$L100505:
	mov	eax, DWORD PTR _RZ2$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100506
	cmp	edi, 1
	jge	SHORT $L100507
$L100506:
	call	___FHranger@0
$L100507:
	mov	eax, DWORD PTR _RY2$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100508
	cmp	edi, 1
	jge	SHORT $L100509
$L100508:
	call	___FHranger@0
$L100509:
	mov	eax, DWORD PTR _RX2$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1460]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV3
	push	eax
	call	_TK@84
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 1173	; 00000495H
	mov	esi, DWORD PTR [eax]

; 1166 :      1            RX2(L),RY2(L),RZ2(L),AA2(L),
; 1167 :      2            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
; 1168 :      3            RX1(L),RY1(L),RZ1(L),AA1(L),
; 1169 :      4            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
; 1170 :      5            B12(L),B12(K),T7)
; 1171 : C     ...Electro/Nuclear Symmetrization Integral:
; 1172 : C         P_[elec.1,elec.2]*P_[H1,H2]
; 1173 :               CALL TK(DFLAG,OV4,

	cmp	edi, esi
	jg	SHORT $L100510
	cmp	edi, 1
	jge	SHORT $L100511
$L100510:
	call	___FHranger@0
$L100511:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98451[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100512
	cmp	edi, 1
	jge	SHORT $L100513
$L100512:
	call	___FHranger@0
$L100513:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98452[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100514
	cmp	edi, 1
	jge	SHORT $L100515
$L100514:
	call	___FHranger@0
$L100515:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98453[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100516
	cmp	edi, 1
	jge	SHORT $L100517
$L100516:
	call	___FHranger@0
$L100517:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98454[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100518
	cmp	edi, 1
	jge	SHORT $L100519
$L100518:
	call	___FHranger@0
$L100519:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98455[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100520
	cmp	edi, 1
	jge	SHORT $L100521
$L100520:
	call	___FHranger@0
$L100521:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	push	OFFSET FLAT:_T8
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98456[esp+1396]
	cmp	edi, esi
	jg	SHORT $L100522
	cmp	edi, 1
	jge	SHORT $L100523
$L100522:
	call	___FHranger@0
$L100523:
	mov	eax, DWORD PTR _B12$[esp+1392]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100524
	cmp	edi, 1
	jge	SHORT $L100525
$L100524:
	call	___FHranger@0
$L100525:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100526
	cmp	edi, 1
	jge	SHORT $L100527
$L100526:
	call	___FHranger@0
$L100527:
	mov	eax, DWORD PTR _AA2$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98451[esp+1404]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98452[esp+1412]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98453[esp+1416]
	push	eax
	jg	SHORT $L100528
	cmp	edi, 1
	jge	SHORT $L100529
$L100528:
	call	___FHranger@0
$L100529:
	mov	eax, DWORD PTR _AA1$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100530
	cmp	edi, 1
	jge	SHORT $L100531
$L100530:
	call	___FHranger@0
$L100531:
	mov	eax, DWORD PTR _RZ1$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100532
	cmp	edi, 1
	jge	SHORT $L100533
$L100532:
	call	___FHranger@0
$L100533:
	mov	eax, DWORD PTR _RY1$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100534
	cmp	edi, 1
	jge	SHORT $L100535
$L100534:
	call	___FHranger@0
$L100535:
	mov	eax, DWORD PTR _RX1$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100536
	cmp	edi, 1
	jge	SHORT $L100537
$L100536:
	call	___FHranger@0
$L100537:
	mov	eax, DWORD PTR _AA1$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98454[esp+1436]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98455[esp+1444]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98456[esp+1448]
	push	eax
	jg	SHORT $L100538
	cmp	edi, 1
	jge	SHORT $L100539
$L100538:
	call	___FHranger@0
$L100539:
	mov	eax, DWORD PTR _AA2$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100540
	cmp	edi, 1
	jge	SHORT $L100541
$L100540:
	call	___FHranger@0
$L100541:
	mov	eax, DWORD PTR _RZ2$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100542
	cmp	edi, 1
	jge	SHORT $L100543
$L100542:
	call	___FHranger@0
$L100543:
	mov	eax, DWORD PTR _RY2$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100544
	cmp	edi, 1
	jge	SHORT $L100545
$L100544:
	call	___FHranger@0
$L100545:
	mov	eax, DWORD PTR _RX2$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1460]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV4
	push	eax
	call	_TK@84
	fld	QWORD PTR _T7
	fadd	QWORD PTR _T4
	mov	eax, DWORD PTR _M$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN1
	mov	edi, DWORD PTR _OVER+4
	fadd	QWORD PTR _T6
	mov	esi, DWORD PTR [eax]
	push	OFFSET FLAT:_EZ1
	push	OFFSET FLAT:_EY1
	cmp	edi, esi
	fadd	QWORD PTR _T3
	mov	DWORD PTR $T98573[esp+1420], 1182	; 0000049eH

; 1174 :      1            RX2(L),RY2(L),RZ2(L),AA2(L),
; 1175 :      2            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
; 1176 :      3            RX1(L),RY1(L),RZ1(L),AA1(L),
; 1177 :      4            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
; 1178 :      5            B12(L),B12(K),T8)
; 1179 : 	END IF
; 1180 : C
; 1181 : 	IF (ISYMM) THEN
; 1182 :             T=(T1+T2+T3+T4+T5+T6+T7+T8)/2.0D+00

	mov	DWORD PTR $T98573[esp+1420], 1203	; 000004b3H
	fadd	QWORD PTR _T1
	push	OFFSET FLAT:_EX1
	fadd	QWORD PTR _T8
	fadd	QWORD PTR _T5
	fadd	QWORD PTR _T2
	fdiv	QWORD PTR $T99039
	fstp	QWORD PTR _T

; 1183 : 	ELSE
; 1184 : 		  T=(T1+T5)/2.0D+00
; 1185 : 	END IF
; 1186 : 	IF(IPERM.EQ.1) THEN
; 1187 : 		T=(T1+T5)/2.0D+00
; 1188 : 	END IF
; 1189 : 	IF(IPERM.EQ.2) THEN
; 1190 : 		T=(T2+T6)/2.0D+00
; 1191 : 	END IF
; 1192 : 	IF(IPERM.EQ.3) THEN
; 1193 : 		T=(T3+T7)/2.0D+00
; 1194 : 	END IF
; 1195 : 	IF(IPERM.EQ.4) THEN
; 1196 : 		T=(T4+T8)/2.0D+00
; 1197 : 	END IF
; 1198 : 
; 1199 : C	write(6,*)L,K,T
; 1200 : C     -----------------------------------------------------
; 1201 : C     --- Calculate Nuc./Elec. Attractn. FOR ELECTRON 1 ---
; 1202 : C     -----------------------------------------------------
; 1203 :               CALL ENUCL(DFLAG,OV1,RX1(L),RY1(L),RZ1(L),AA1(L),

	jg	SHORT $L100546
	cmp	edi, 1
	jge	SHORT $L100547
$L100546:
	call	___FHranger@0
$L100547:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100548
	cmp	edi, 1
	jge	SHORT $L100549
$L100548:
	call	___FHranger@0
$L100549:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100550
	cmp	edi, 1
	jge	SHORT $L100551
$L100550:
	call	___FHranger@0
$L100551:
	mov	eax, DWORD PTR _AA2$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100552
	cmp	edi, 1
	jge	SHORT $L100553
$L100552:
	call	___FHranger@0
$L100553:
	mov	eax, DWORD PTR _RZ2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100554
	cmp	edi, 1
	jge	SHORT $L100555
$L100554:
	call	___FHranger@0
$L100555:
	mov	eax, DWORD PTR _RY2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100556
	cmp	edi, 1
	jge	SHORT $L100557
$L100556:
	call	___FHranger@0
$L100557:
	mov	eax, DWORD PTR _RX2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100558
	cmp	edi, 1
	jge	SHORT $L100559
$L100558:
	call	___FHranger@0
$L100559:
	mov	eax, DWORD PTR _AA2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100560
	cmp	edi, 1
	jge	SHORT $L100561
$L100560:
	call	___FHranger@0
$L100561:
	mov	eax, DWORD PTR _RZ2$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100562
	cmp	edi, 1
	jge	SHORT $L100563
$L100562:
	call	___FHranger@0
$L100563:
	mov	eax, DWORD PTR _RY2$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100564
	cmp	edi, 1
	jge	SHORT $L100565
$L100564:
	call	___FHranger@0
$L100565:
	mov	eax, DWORD PTR _RX2$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100566
	cmp	edi, 1
	jge	SHORT $L100567
$L100566:
	call	___FHranger@0
$L100567:
	mov	eax, DWORD PTR _AA1$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100568
	cmp	edi, 1
	jge	SHORT $L100569
$L100568:
	call	___FHranger@0
$L100569:
	mov	eax, DWORD PTR _RZ1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100570
	cmp	edi, 1
	jge	SHORT $L100571
$L100570:
	call	___FHranger@0
$L100571:
	mov	eax, DWORD PTR _RY1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100572
	cmp	edi, 1
	jge	SHORT $L100573
$L100572:
	call	___FHranger@0
$L100573:
	mov	eax, DWORD PTR _RX1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100574
	cmp	edi, 1
	jge	SHORT $L100575
$L100574:
	call	___FHranger@0
$L100575:
	mov	eax, DWORD PTR _AA1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100576
	cmp	edi, 1
	jge	SHORT $L100577
$L100576:
	call	___FHranger@0
$L100577:
	mov	eax, DWORD PTR _RZ1$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100578
	cmp	edi, 1
	jge	SHORT $L100579
$L100578:
	call	___FHranger@0
$L100579:
	mov	eax, DWORD PTR _RY1$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100580
	cmp	edi, 1
	jge	SHORT $L100581
$L100580:
	call	___FHranger@0
$L100581:
	mov	eax, DWORD PTR _RX1$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV1
	push	eax
	call	_ENUCL@100
	mov	DWORD PTR $T98573[esp+1404], 1211	; 000004bbH

; 1204 :      1           RX1(K),RY1(K),RZ1(K),AA1(K),                      
; 1205 :      2           RX2(L),RY2(L),RZ2(L),AA2(L),
; 1206 :      3           RX2(K),RY2(K),RZ2(K),AA2(K),
; 1207 :      4           B12(L),B12(K),
; 1208 :      5           EX1,EY1,EZ1,EN1,KLM)
; 1209 : C
; 1210 : 
; 1211 : 	WRITE(6,*)'EN1=',EN1

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	4
	push	OFFSET FLAT:_.SG596_
	call	_lwrite_ch
	add	esp, 8
	mov	eax, DWORD PTR _EN1+4
	push	eax
	mov	eax, DWORD PTR _EN1
	push	eax
	call	_lwrite_d
	add	esp, 8
	call	_f90_endwrite
	mov	eax, DWORD PTR _M$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN2
	mov	edi, DWORD PTR _OVER+4
	mov	esi, DWORD PTR [eax]
	push	OFFSET FLAT:_EZ1
	push	OFFSET FLAT:_EY1
	cmp	edi, esi
	mov	DWORD PTR $T98573[esp+1420], 1216	; 000004c0H

; 1212 : 
; 1213 : C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
; 1214 : C
; 1215 : 	IF (ISYMM) THEN
; 1216 :               CALL ENUCL(DFLAG,OV2,RX1(L),RY1(L),RZ1(L),AA1(L),

	push	OFFSET FLAT:_EX1
	jg	SHORT $L100582
	cmp	edi, 1
	jge	SHORT $L100583
$L100582:
	call	___FHranger@0
$L100583:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100584
	cmp	edi, 1
	jge	SHORT $L100585
$L100584:
	call	___FHranger@0
$L100585:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100586
	cmp	edi, 1
	jge	SHORT $L100587
$L100586:
	call	___FHranger@0
$L100587:
	mov	eax, DWORD PTR _AA1$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100588
	cmp	edi, 1
	jge	SHORT $L100589
$L100588:
	call	___FHranger@0
$L100589:
	mov	eax, DWORD PTR _RZ1$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100590
	cmp	edi, 1
	jge	SHORT $L100591
$L100590:
	call	___FHranger@0
$L100591:
	mov	eax, DWORD PTR _RY1$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100592
	cmp	edi, 1
	jge	SHORT $L100593
$L100592:
	call	___FHranger@0
$L100593:
	mov	eax, DWORD PTR _RX1$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100594
	cmp	edi, 1
	jge	SHORT $L100595
$L100594:
	call	___FHranger@0
$L100595:
	mov	eax, DWORD PTR _AA2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100596
	cmp	edi, 1
	jge	SHORT $L100597
$L100596:
	call	___FHranger@0
$L100597:
	mov	eax, DWORD PTR _RZ2$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100598
	cmp	edi, 1
	jge	SHORT $L100599
$L100598:
	call	___FHranger@0
$L100599:
	mov	eax, DWORD PTR _RY2$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100600
	cmp	edi, 1
	jge	SHORT $L100601
$L100600:
	call	___FHranger@0
$L100601:
	mov	eax, DWORD PTR _RX2$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100602
	cmp	edi, 1
	jge	SHORT $L100603
$L100602:
	call	___FHranger@0
$L100603:
	mov	eax, DWORD PTR _AA2$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100604
	cmp	edi, 1
	jge	SHORT $L100605
$L100604:
	call	___FHranger@0
$L100605:
	mov	eax, DWORD PTR _RZ2$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100606
	cmp	edi, 1
	jge	SHORT $L100607
$L100606:
	call	___FHranger@0
$L100607:
	mov	eax, DWORD PTR _RY2$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100608
	cmp	edi, 1
	jge	SHORT $L100609
$L100608:
	call	___FHranger@0
$L100609:
	mov	eax, DWORD PTR _RX2$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100610
	cmp	edi, 1
	jge	SHORT $L100611
$L100610:
	call	___FHranger@0
$L100611:
	mov	eax, DWORD PTR _AA1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100612
	cmp	edi, 1
	jge	SHORT $L100613
$L100612:
	call	___FHranger@0
$L100613:
	mov	eax, DWORD PTR _RZ1$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100614
	cmp	edi, 1
	jge	SHORT $L100615
$L100614:
	call	___FHranger@0
$L100615:
	mov	eax, DWORD PTR _RY1$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100616
	cmp	edi, 1
	jge	SHORT $L100617
$L100616:
	call	___FHranger@0
$L100617:
	mov	eax, DWORD PTR _RX1$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV2
	push	eax
	call	_ENUCL@100
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 1225	; 000004c9H
	mov	esi, DWORD PTR [eax]

; 1217 :      1           RX2(K),RY2(K),RZ2(K),AA2(K),                      
; 1218 :      2           RX2(L),RY2(L),RZ2(L),AA2(L),
; 1219 :      3           RX1(K),RY1(K),RZ1(K),AA1(K),
; 1220 :      4           B12(L),B12(K),
; 1221 :      5           EX1,EY1,EZ1,EN2,KLM)                       
; 1222 : C
; 1223 : C     ...Nuclear Symmetrization Integral: P_[H1,H2]
; 1224 : C
; 1225 :               CALL ENUCL(DFLAG,OV3,RX1(L),RY1(L),RZ1(L),AA1(L),

	cmp	edi, esi
	jg	SHORT $L100618
	cmp	edi, 1
	jge	SHORT $L100619
$L100618:
	call	___FHranger@0
$L100619:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98457[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100620
	cmp	edi, 1
	jge	SHORT $L100621
$L100620:
	call	___FHranger@0
$L100621:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98458[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100622
	cmp	edi, 1
	jge	SHORT $L100623
$L100622:
	call	___FHranger@0
$L100623:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98459[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100624
	cmp	edi, 1
	jge	SHORT $L100625
$L100624:
	call	___FHranger@0
$L100625:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98460[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100626
	cmp	edi, 1
	jge	SHORT $L100627
$L100626:
	call	___FHranger@0
$L100627:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98461[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100628
	cmp	edi, 1
	jge	SHORT $L100629
$L100628:
	call	___FHranger@0
$L100629:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN3
	push	OFFSET FLAT:_EZ1
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	push	OFFSET FLAT:_EY1
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98462[esp+1408]
	push	OFFSET FLAT:_EX1
	cmp	edi, esi
	jg	SHORT $L100630
	cmp	edi, 1
	jge	SHORT $L100631
$L100630:
	call	___FHranger@0
$L100631:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100632
	cmp	edi, 1
	jge	SHORT $L100633
$L100632:
	call	___FHranger@0
$L100633:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100634
	cmp	edi, 1
	jge	SHORT $L100635
$L100634:
	call	___FHranger@0
$L100635:
	mov	eax, DWORD PTR _AA2$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98457[esp+1420]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98458[esp+1428]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98459[esp+1432]
	push	eax
	jg	SHORT $L100636
	cmp	edi, 1
	jge	SHORT $L100637
$L100636:
	call	___FHranger@0
$L100637:
	mov	eax, DWORD PTR _AA2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100638
	cmp	edi, 1
	jge	SHORT $L100639
$L100638:
	call	___FHranger@0
$L100639:
	mov	eax, DWORD PTR _RZ2$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100640
	cmp	edi, 1
	jge	SHORT $L100641
$L100640:
	call	___FHranger@0
$L100641:
	mov	eax, DWORD PTR _RY2$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100642
	cmp	edi, 1
	jge	SHORT $L100643
$L100642:
	call	___FHranger@0
$L100643:
	mov	eax, DWORD PTR _RX2$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100644
	cmp	edi, 1
	jge	SHORT $L100645
$L100644:
	call	___FHranger@0
$L100645:
	mov	eax, DWORD PTR _AA1$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98460[esp+1452]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98461[esp+1460]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98462[esp+1464]
	push	eax
	jg	SHORT $L100646
	cmp	edi, 1
	jge	SHORT $L100647
$L100646:
	call	___FHranger@0
$L100647:
	mov	eax, DWORD PTR _AA1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100648
	cmp	edi, 1
	jge	SHORT $L100649
$L100648:
	call	___FHranger@0
$L100649:
	mov	eax, DWORD PTR _RZ1$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100650
	cmp	edi, 1
	jge	SHORT $L100651
$L100650:
	call	___FHranger@0
$L100651:
	mov	eax, DWORD PTR _RY1$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100652
	cmp	edi, 1
	jge	SHORT $L100653
$L100652:
	call	___FHranger@0
$L100653:
	mov	eax, DWORD PTR _RX1$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV3
	push	eax
	call	_ENUCL@100
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 1235	; 000004d3H
	mov	esi, DWORD PTR [eax]

; 1226 :      1            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
; 1227 :      2            RX2(L),RY2(L),RZ2(L),AA2(L),
; 1228 :      3            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
; 1229 :      4            B12(L),B12(K),
; 1230 :      5            EX1,EY1,EZ1,EN3,KLM)
; 1231 : C
; 1232 : C     ...Electro/Nuclear Symmetrization Integral:
; 1233 : C         P_[elec.1,elec.2]*P_[H1,H2]
; 1234 : C
; 1235 :               CALL ENUCL(DFLAG,OV4,RX1(L),RY1(L),RZ1(L),AA1(L),

	cmp	edi, esi
	jg	SHORT $L100654
	cmp	edi, 1
	jge	SHORT $L100655
$L100654:
	call	___FHranger@0
$L100655:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98463[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100656
	cmp	edi, 1
	jge	SHORT $L100657
$L100656:
	call	___FHranger@0
$L100657:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98464[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100658
	cmp	edi, 1
	jge	SHORT $L100659
$L100658:
	call	___FHranger@0
$L100659:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98465[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100660
	cmp	edi, 1
	jge	SHORT $L100661
$L100660:
	call	___FHranger@0
$L100661:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98466[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100662
	cmp	edi, 1
	jge	SHORT $L100663
$L100662:
	call	___FHranger@0
$L100663:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98467[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100664
	cmp	edi, 1
	jge	SHORT $L100665
$L100664:
	call	___FHranger@0
$L100665:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN4
	push	OFFSET FLAT:_EZ1
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	push	OFFSET FLAT:_EY1
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98468[esp+1408]
	push	OFFSET FLAT:_EX1
	cmp	edi, esi
	jg	SHORT $L100666
	cmp	edi, 1
	jge	SHORT $L100667
$L100666:
	call	___FHranger@0
$L100667:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100668
	cmp	edi, 1
	jge	SHORT $L100669
$L100668:
	call	___FHranger@0
$L100669:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100670
	cmp	edi, 1
	jge	SHORT $L100671
$L100670:
	call	___FHranger@0
$L100671:
	mov	eax, DWORD PTR _AA1$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98463[esp+1420]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98464[esp+1428]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98465[esp+1432]
	push	eax
	jg	SHORT $L100672
	cmp	edi, 1
	jge	SHORT $L100673
$L100672:
	call	___FHranger@0
$L100673:
	mov	eax, DWORD PTR _AA2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100674
	cmp	edi, 1
	jge	SHORT $L100675
$L100674:
	call	___FHranger@0
$L100675:
	mov	eax, DWORD PTR _RZ2$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100676
	cmp	edi, 1
	jge	SHORT $L100677
$L100676:
	call	___FHranger@0
$L100677:
	mov	eax, DWORD PTR _RY2$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100678
	cmp	edi, 1
	jge	SHORT $L100679
$L100678:
	call	___FHranger@0
$L100679:
	mov	eax, DWORD PTR _RX2$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100680
	cmp	edi, 1
	jge	SHORT $L100681
$L100680:
	call	___FHranger@0
$L100681:
	mov	eax, DWORD PTR _AA2$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98466[esp+1452]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98467[esp+1460]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98468[esp+1464]
	push	eax
	jg	SHORT $L100682
	cmp	edi, 1
	jge	SHORT $L100683
$L100682:
	call	___FHranger@0
$L100683:
	mov	eax, DWORD PTR _AA1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100684
	cmp	edi, 1
	jge	SHORT $L100685
$L100684:
	call	___FHranger@0
$L100685:
	mov	eax, DWORD PTR _RZ1$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100686
	cmp	edi, 1
	jge	SHORT $L100687
$L100686:
	call	___FHranger@0
$L100687:
	mov	eax, DWORD PTR _RY1$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100688
	cmp	edi, 1
	jge	SHORT $L100689
$L100688:
	call	___FHranger@0
$L100689:
	mov	eax, DWORD PTR _RX1$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV4
	push	eax
	call	_ENUCL@100
	mov	eax, DWORD PTR _M$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN5
	mov	edi, DWORD PTR _OVER+4
	mov	esi, DWORD PTR [eax]
	push	OFFSET FLAT:_EZ2
	push	OFFSET FLAT:_EY2
	cmp	edi, esi
	mov	DWORD PTR $T98573[esp+1420], 1243	; 000004dbH

; 1236 :      1            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
; 1237 :      2            RX2(L),RY2(L),RZ2(L),AA2(L),
; 1238 :      3            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
; 1239 :      4            B12(L),B12(K),
; 1240 :      5            EX1,EY1,EZ1,EN4,KLM)
; 1241 : 	END IF
; 1242 : C  
; 1243 :               CALL ENUCL(DFLAG,OV1,RX1(L),RY1(L),RZ1(L),AA1(L),

	push	OFFSET FLAT:_EX2
	jg	SHORT $L100690
	cmp	edi, 1
	jge	SHORT $L100691
$L100690:
	call	___FHranger@0
$L100691:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100692
	cmp	edi, 1
	jge	SHORT $L100693
$L100692:
	call	___FHranger@0
$L100693:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100694
	cmp	edi, 1
	jge	SHORT $L100695
$L100694:
	call	___FHranger@0
$L100695:
	mov	eax, DWORD PTR _AA2$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100696
	cmp	edi, 1
	jge	SHORT $L100697
$L100696:
	call	___FHranger@0
$L100697:
	mov	eax, DWORD PTR _RZ2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100698
	cmp	edi, 1
	jge	SHORT $L100699
$L100698:
	call	___FHranger@0
$L100699:
	mov	eax, DWORD PTR _RY2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100700
	cmp	edi, 1
	jge	SHORT $L100701
$L100700:
	call	___FHranger@0
$L100701:
	mov	eax, DWORD PTR _RX2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100702
	cmp	edi, 1
	jge	SHORT $L100703
$L100702:
	call	___FHranger@0
$L100703:
	mov	eax, DWORD PTR _AA2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100704
	cmp	edi, 1
	jge	SHORT $L100705
$L100704:
	call	___FHranger@0
$L100705:
	mov	eax, DWORD PTR _RZ2$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100706
	cmp	edi, 1
	jge	SHORT $L100707
$L100706:
	call	___FHranger@0
$L100707:
	mov	eax, DWORD PTR _RY2$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100708
	cmp	edi, 1
	jge	SHORT $L100709
$L100708:
	call	___FHranger@0
$L100709:
	mov	eax, DWORD PTR _RX2$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100710
	cmp	edi, 1
	jge	SHORT $L100711
$L100710:
	call	___FHranger@0
$L100711:
	mov	eax, DWORD PTR _AA1$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100712
	cmp	edi, 1
	jge	SHORT $L100713
$L100712:
	call	___FHranger@0
$L100713:
	mov	eax, DWORD PTR _RZ1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100714
	cmp	edi, 1
	jge	SHORT $L100715
$L100714:
	call	___FHranger@0
$L100715:
	mov	eax, DWORD PTR _RY1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100716
	cmp	edi, 1
	jge	SHORT $L100717
$L100716:
	call	___FHranger@0
$L100717:
	mov	eax, DWORD PTR _RX1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100718
	cmp	edi, 1
	jge	SHORT $L100719
$L100718:
	call	___FHranger@0
$L100719:
	mov	eax, DWORD PTR _AA1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100720
	cmp	edi, 1
	jge	SHORT $L100721
$L100720:
	call	___FHranger@0
$L100721:
	mov	eax, DWORD PTR _RZ1$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100722
	cmp	edi, 1
	jge	SHORT $L100723
$L100722:
	call	___FHranger@0
$L100723:
	mov	eax, DWORD PTR _RY1$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100724
	cmp	edi, 1
	jge	SHORT $L100725
$L100724:
	call	___FHranger@0
$L100725:
	mov	eax, DWORD PTR _RX1$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV1
	push	eax
	call	_ENUCL@100
	mov	eax, DWORD PTR _M$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN6
	mov	edi, DWORD PTR _OVER+4
	mov	esi, DWORD PTR [eax]
	push	OFFSET FLAT:_EZ2
	push	OFFSET FLAT:_EY2
	cmp	edi, esi
	mov	DWORD PTR $T98573[esp+1420], 1253	; 000004e5H

; 1244 :      1           RX1(K),RY1(K),RZ1(K),AA1(K),                      
; 1245 :      2           RX2(L),RY2(L),RZ2(L),AA2(L),
; 1246 :      3           RX2(K),RY2(K),RZ2(K),AA2(K),
; 1247 :      4           B12(L),B12(K),
; 1248 :      5           EX2,EY2,EZ2,EN5,KLM)
; 1249 : C
; 1250 : C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
; 1251 : C
; 1252 : 	IF (ISYMM) THEN
; 1253 :               CALL ENUCL(DFLAG,OV2,RX1(L),RY1(L),RZ1(L),AA1(L),

	push	OFFSET FLAT:_EX2
	jg	SHORT $L100726
	cmp	edi, 1
	jge	SHORT $L100727
$L100726:
	call	___FHranger@0
$L100727:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100728
	cmp	edi, 1
	jge	SHORT $L100729
$L100728:
	call	___FHranger@0
$L100729:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100730
	cmp	edi, 1
	jge	SHORT $L100731
$L100730:
	call	___FHranger@0
$L100731:
	mov	eax, DWORD PTR _AA1$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100732
	cmp	edi, 1
	jge	SHORT $L100733
$L100732:
	call	___FHranger@0
$L100733:
	mov	eax, DWORD PTR _RZ1$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100734
	cmp	edi, 1
	jge	SHORT $L100735
$L100734:
	call	___FHranger@0
$L100735:
	mov	eax, DWORD PTR _RY1$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100736
	cmp	edi, 1
	jge	SHORT $L100737
$L100736:
	call	___FHranger@0
$L100737:
	mov	eax, DWORD PTR _RX1$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100738
	cmp	edi, 1
	jge	SHORT $L100739
$L100738:
	call	___FHranger@0
$L100739:
	mov	eax, DWORD PTR _AA2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100740
	cmp	edi, 1
	jge	SHORT $L100741
$L100740:
	call	___FHranger@0
$L100741:
	mov	eax, DWORD PTR _RZ2$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100742
	cmp	edi, 1
	jge	SHORT $L100743
$L100742:
	call	___FHranger@0
$L100743:
	mov	eax, DWORD PTR _RY2$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100744
	cmp	edi, 1
	jge	SHORT $L100745
$L100744:
	call	___FHranger@0
$L100745:
	mov	eax, DWORD PTR _RX2$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100746
	cmp	edi, 1
	jge	SHORT $L100747
$L100746:
	call	___FHranger@0
$L100747:
	mov	eax, DWORD PTR _AA2$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100748
	cmp	edi, 1
	jge	SHORT $L100749
$L100748:
	call	___FHranger@0
$L100749:
	mov	eax, DWORD PTR _RZ2$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100750
	cmp	edi, 1
	jge	SHORT $L100751
$L100750:
	call	___FHranger@0
$L100751:
	mov	eax, DWORD PTR _RY2$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100752
	cmp	edi, 1
	jge	SHORT $L100753
$L100752:
	call	___FHranger@0
$L100753:
	mov	eax, DWORD PTR _RX2$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100754
	cmp	edi, 1
	jge	SHORT $L100755
$L100754:
	call	___FHranger@0
$L100755:
	mov	eax, DWORD PTR _AA1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100756
	cmp	edi, 1
	jge	SHORT $L100757
$L100756:
	call	___FHranger@0
$L100757:
	mov	eax, DWORD PTR _RZ1$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100758
	cmp	edi, 1
	jge	SHORT $L100759
$L100758:
	call	___FHranger@0
$L100759:
	mov	eax, DWORD PTR _RY1$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100760
	cmp	edi, 1
	jge	SHORT $L100761
$L100760:
	call	___FHranger@0
$L100761:
	mov	eax, DWORD PTR _RX1$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV2
	push	eax
	call	_ENUCL@100
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 1262	; 000004eeH
	mov	esi, DWORD PTR [eax]

; 1254 :      1           RX2(K),RY2(K),RZ2(K),AA2(K),                      
; 1255 :      2           RX2(L),RY2(L),RZ2(L),AA2(L),
; 1256 :      3           RX1(K),RY1(K),RZ1(K),AA1(K),
; 1257 :      4           B12(L),B12(K),
; 1258 :      5           EX2,EY2,EZ2,EN6,KLM)                       
; 1259 : C
; 1260 : C     ...Nuclear Symmetrization Integral: P_[H1,H2]
; 1261 : C
; 1262 :               CALL ENUCL(DFLAG,OV3,RX1(L),RY1(L),RZ1(L),AA1(L),

	cmp	edi, esi
	jg	SHORT $L100762
	cmp	edi, 1
	jge	SHORT $L100763
$L100762:
	call	___FHranger@0
$L100763:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98469[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100764
	cmp	edi, 1
	jge	SHORT $L100765
$L100764:
	call	___FHranger@0
$L100765:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98470[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100766
	cmp	edi, 1
	jge	SHORT $L100767
$L100766:
	call	___FHranger@0
$L100767:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98471[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100768
	cmp	edi, 1
	jge	SHORT $L100769
$L100768:
	call	___FHranger@0
$L100769:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98472[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100770
	cmp	edi, 1
	jge	SHORT $L100771
$L100770:
	call	___FHranger@0
$L100771:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98473[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100772
	cmp	edi, 1
	jge	SHORT $L100773
$L100772:
	call	___FHranger@0
$L100773:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN7
	push	OFFSET FLAT:_EZ2
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	push	OFFSET FLAT:_EY2
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98474[esp+1408]
	push	OFFSET FLAT:_EX2
	cmp	edi, esi
	jg	SHORT $L100774
	cmp	edi, 1
	jge	SHORT $L100775
$L100774:
	call	___FHranger@0
$L100775:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100776
	cmp	edi, 1
	jge	SHORT $L100777
$L100776:
	call	___FHranger@0
$L100777:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100778
	cmp	edi, 1
	jge	SHORT $L100779
$L100778:
	call	___FHranger@0
$L100779:
	mov	eax, DWORD PTR _AA2$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98469[esp+1420]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98470[esp+1428]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98471[esp+1432]
	push	eax
	jg	SHORT $L100780
	cmp	edi, 1
	jge	SHORT $L100781
$L100780:
	call	___FHranger@0
$L100781:
	mov	eax, DWORD PTR _AA2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100782
	cmp	edi, 1
	jge	SHORT $L100783
$L100782:
	call	___FHranger@0
$L100783:
	mov	eax, DWORD PTR _RZ2$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100784
	cmp	edi, 1
	jge	SHORT $L100785
$L100784:
	call	___FHranger@0
$L100785:
	mov	eax, DWORD PTR _RY2$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100786
	cmp	edi, 1
	jge	SHORT $L100787
$L100786:
	call	___FHranger@0
$L100787:
	mov	eax, DWORD PTR _RX2$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100788
	cmp	edi, 1
	jge	SHORT $L100789
$L100788:
	call	___FHranger@0
$L100789:
	mov	eax, DWORD PTR _AA1$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98472[esp+1452]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98473[esp+1460]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98474[esp+1464]
	push	eax
	jg	SHORT $L100790
	cmp	edi, 1
	jge	SHORT $L100791
$L100790:
	call	___FHranger@0
$L100791:
	mov	eax, DWORD PTR _AA1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100792
	cmp	edi, 1
	jge	SHORT $L100793
$L100792:
	call	___FHranger@0
$L100793:
	mov	eax, DWORD PTR _RZ1$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100794
	cmp	edi, 1
	jge	SHORT $L100795
$L100794:
	call	___FHranger@0
$L100795:
	mov	eax, DWORD PTR _RY1$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100796
	cmp	edi, 1
	jge	SHORT $L100797
$L100796:
	call	___FHranger@0
$L100797:
	mov	eax, DWORD PTR _RX1$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV3
	push	eax
	call	_ENUCL@100
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 1272	; 000004f8H
	mov	esi, DWORD PTR [eax]

; 1263 :      1            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
; 1264 :      2            RX2(L),RY2(L),RZ2(L),AA2(L),
; 1265 :      3            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
; 1266 :      4            B12(L),B12(K),
; 1267 :      5            EX2,EY2,EZ2,EN7,KLM)
; 1268 : C
; 1269 : C     ...Electro/Nuclear Symmetrization Integral:
; 1270 : C         P_[elec.1,elec.2]*P_[H1,H2]
; 1271 : C
; 1272 :               CALL ENUCL(DFLAG,OV4,RX1(L),RY1(L),RZ1(L),AA1(L),

	cmp	edi, esi
	jg	SHORT $L100798
	cmp	edi, 1
	jge	SHORT $L100799
$L100798:
	call	___FHranger@0
$L100799:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98475[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100800
	cmp	edi, 1
	jge	SHORT $L100801
$L100800:
	call	___FHranger@0
$L100801:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98476[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100802
	cmp	edi, 1
	jge	SHORT $L100803
$L100802:
	call	___FHranger@0
$L100803:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98477[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100804
	cmp	edi, 1
	jge	SHORT $L100805
$L100804:
	call	___FHranger@0
$L100805:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98478[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100806
	cmp	edi, 1
	jge	SHORT $L100807
$L100806:
	call	___FHranger@0
$L100807:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98479[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100808
	cmp	edi, 1
	jge	SHORT $L100809
$L100808:
	call	___FHranger@0
$L100809:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN8
	push	OFFSET FLAT:_EZ2
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	push	OFFSET FLAT:_EY2
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98480[esp+1408]
	push	OFFSET FLAT:_EX2
	cmp	edi, esi
	jg	SHORT $L100810
	cmp	edi, 1
	jge	SHORT $L100811
$L100810:
	call	___FHranger@0
$L100811:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100812
	cmp	edi, 1
	jge	SHORT $L100813
$L100812:
	call	___FHranger@0
$L100813:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100814
	cmp	edi, 1
	jge	SHORT $L100815
$L100814:
	call	___FHranger@0
$L100815:
	mov	eax, DWORD PTR _AA1$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98475[esp+1420]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98476[esp+1428]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98477[esp+1432]
	push	eax
	jg	SHORT $L100816
	cmp	edi, 1
	jge	SHORT $L100817
$L100816:
	call	___FHranger@0
$L100817:
	mov	eax, DWORD PTR _AA2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100818
	cmp	edi, 1
	jge	SHORT $L100819
$L100818:
	call	___FHranger@0
$L100819:
	mov	eax, DWORD PTR _RZ2$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100820
	cmp	edi, 1
	jge	SHORT $L100821
$L100820:
	call	___FHranger@0
$L100821:
	mov	eax, DWORD PTR _RY2$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100822
	cmp	edi, 1
	jge	SHORT $L100823
$L100822:
	call	___FHranger@0
$L100823:
	mov	eax, DWORD PTR _RX2$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100824
	cmp	edi, 1
	jge	SHORT $L100825
$L100824:
	call	___FHranger@0
$L100825:
	mov	eax, DWORD PTR _AA2$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98478[esp+1452]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98479[esp+1460]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98480[esp+1464]
	push	eax
	jg	SHORT $L100826
	cmp	edi, 1
	jge	SHORT $L100827
$L100826:
	call	___FHranger@0
$L100827:
	mov	eax, DWORD PTR _AA1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100828
	cmp	edi, 1
	jge	SHORT $L100829
$L100828:
	call	___FHranger@0
$L100829:
	mov	eax, DWORD PTR _RZ1$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100830
	cmp	edi, 1
	jge	SHORT $L100831
$L100830:
	call	___FHranger@0
$L100831:
	mov	eax, DWORD PTR _RY1$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100832
	cmp	edi, 1
	jge	SHORT $L100833
$L100832:
	call	___FHranger@0
$L100833:
	mov	eax, DWORD PTR _RX1$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV4
	push	eax
	call	_ENUCL@100
	mov	eax, DWORD PTR _M$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN9
	mov	edi, DWORD PTR _OVER+4
	mov	esi, DWORD PTR [eax]
	push	OFFSET FLAT:_EZ1
	push	OFFSET FLAT:_EY1
	cmp	edi, esi
	mov	DWORD PTR $T98573[esp+1420], 1283	; 00000503H

; 1273 :      1            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
; 1274 :      2            RX2(L),RY2(L),RZ2(L),AA2(L),
; 1275 :      3            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
; 1276 :      4            B12(L),B12(K),
; 1277 :      5            EX2,EY2,EZ2,EN8,KLM)
; 1278 : 	END IF
; 1279 : C 
; 1280 : C     -----------------------------------------------------
; 1281 : C     --- Calculate Nuc./Elec. Attractn. FOR ELECTRON 2 ---
; 1282 : C     -----------------------------------------------------
; 1283 :               CALL ENUCL(DFLAG,OV1,RX2(L),RY2(L),RZ2(L),AA2(L),

	push	OFFSET FLAT:_EX1
	jg	SHORT $L100834
	cmp	edi, 1
	jge	SHORT $L100835
$L100834:
	call	___FHranger@0
$L100835:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100836
	cmp	edi, 1
	jge	SHORT $L100837
$L100836:
	call	___FHranger@0
$L100837:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100838
	cmp	edi, 1
	jge	SHORT $L100839
$L100838:
	call	___FHranger@0
$L100839:
	mov	eax, DWORD PTR _AA1$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100840
	cmp	edi, 1
	jge	SHORT $L100841
$L100840:
	call	___FHranger@0
$L100841:
	mov	eax, DWORD PTR _RZ1$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100842
	cmp	edi, 1
	jge	SHORT $L100843
$L100842:
	call	___FHranger@0
$L100843:
	mov	eax, DWORD PTR _RY1$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100844
	cmp	edi, 1
	jge	SHORT $L100845
$L100844:
	call	___FHranger@0
$L100845:
	mov	eax, DWORD PTR _RX1$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100846
	cmp	edi, 1
	jge	SHORT $L100847
$L100846:
	call	___FHranger@0
$L100847:
	mov	eax, DWORD PTR _AA1$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100848
	cmp	edi, 1
	jge	SHORT $L100849
$L100848:
	call	___FHranger@0
$L100849:
	mov	eax, DWORD PTR _RZ1$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100850
	cmp	edi, 1
	jge	SHORT $L100851
$L100850:
	call	___FHranger@0
$L100851:
	mov	eax, DWORD PTR _RY1$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100852
	cmp	edi, 1
	jge	SHORT $L100853
$L100852:
	call	___FHranger@0
$L100853:
	mov	eax, DWORD PTR _RX1$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100854
	cmp	edi, 1
	jge	SHORT $L100855
$L100854:
	call	___FHranger@0
$L100855:
	mov	eax, DWORD PTR _AA2$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100856
	cmp	edi, 1
	jge	SHORT $L100857
$L100856:
	call	___FHranger@0
$L100857:
	mov	eax, DWORD PTR _RZ2$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100858
	cmp	edi, 1
	jge	SHORT $L100859
$L100858:
	call	___FHranger@0
$L100859:
	mov	eax, DWORD PTR _RY2$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100860
	cmp	edi, 1
	jge	SHORT $L100861
$L100860:
	call	___FHranger@0
$L100861:
	mov	eax, DWORD PTR _RX2$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100862
	cmp	edi, 1
	jge	SHORT $L100863
$L100862:
	call	___FHranger@0
$L100863:
	mov	eax, DWORD PTR _AA2$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100864
	cmp	edi, 1
	jge	SHORT $L100865
$L100864:
	call	___FHranger@0
$L100865:
	mov	eax, DWORD PTR _RZ2$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100866
	cmp	edi, 1
	jge	SHORT $L100867
$L100866:
	call	___FHranger@0
$L100867:
	mov	eax, DWORD PTR _RY2$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100868
	cmp	edi, 1
	jge	SHORT $L100869
$L100868:
	call	___FHranger@0
$L100869:
	mov	eax, DWORD PTR _RX2$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV1
	push	eax
	call	_ENUCL@100
	mov	eax, DWORD PTR _M$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN10
	mov	edi, DWORD PTR _OVER+4
	mov	esi, DWORD PTR [eax]
	push	OFFSET FLAT:_EZ1
	push	OFFSET FLAT:_EY1
	cmp	edi, esi
	mov	DWORD PTR $T98573[esp+1420], 1293	; 0000050dH

; 1284 :      1           RX2(K),RY2(K),RZ2(K),AA2(K),                      
; 1285 :      2           RX1(L),RY1(L),RZ1(L),AA1(L),
; 1286 :      3           RX1(K),RY1(K),RZ1(K),AA1(K),
; 1287 :      4           B12(L),B12(K),
; 1288 :      5           EX1,EY1,EZ1,EN9,KLM)                       
; 1289 : C
; 1290 : C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
; 1291 : C
; 1292 : 	IF (ISYMM) THEN
; 1293 :               CALL ENUCL(DFLAG,OV2,RX2(L),RY2(L),RZ2(L),AA2(L),

	push	OFFSET FLAT:_EX1
	jg	SHORT $L100870
	cmp	edi, 1
	jge	SHORT $L100871
$L100870:
	call	___FHranger@0
$L100871:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100872
	cmp	edi, 1
	jge	SHORT $L100873
$L100872:
	call	___FHranger@0
$L100873:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100874
	cmp	edi, 1
	jge	SHORT $L100875
$L100874:
	call	___FHranger@0
$L100875:
	mov	eax, DWORD PTR _AA2$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100876
	cmp	edi, 1
	jge	SHORT $L100877
$L100876:
	call	___FHranger@0
$L100877:
	mov	eax, DWORD PTR _RZ2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100878
	cmp	edi, 1
	jge	SHORT $L100879
$L100878:
	call	___FHranger@0
$L100879:
	mov	eax, DWORD PTR _RY2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100880
	cmp	edi, 1
	jge	SHORT $L100881
$L100880:
	call	___FHranger@0
$L100881:
	mov	eax, DWORD PTR _RX2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100882
	cmp	edi, 1
	jge	SHORT $L100883
$L100882:
	call	___FHranger@0
$L100883:
	mov	eax, DWORD PTR _AA1$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100884
	cmp	edi, 1
	jge	SHORT $L100885
$L100884:
	call	___FHranger@0
$L100885:
	mov	eax, DWORD PTR _RZ1$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100886
	cmp	edi, 1
	jge	SHORT $L100887
$L100886:
	call	___FHranger@0
$L100887:
	mov	eax, DWORD PTR _RY1$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100888
	cmp	edi, 1
	jge	SHORT $L100889
$L100888:
	call	___FHranger@0
$L100889:
	mov	eax, DWORD PTR _RX1$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100890
	cmp	edi, 1
	jge	SHORT $L100891
$L100890:
	call	___FHranger@0
$L100891:
	mov	eax, DWORD PTR _AA1$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100892
	cmp	edi, 1
	jge	SHORT $L100893
$L100892:
	call	___FHranger@0
$L100893:
	mov	eax, DWORD PTR _RZ1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100894
	cmp	edi, 1
	jge	SHORT $L100895
$L100894:
	call	___FHranger@0
$L100895:
	mov	eax, DWORD PTR _RY1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100896
	cmp	edi, 1
	jge	SHORT $L100897
$L100896:
	call	___FHranger@0
$L100897:
	mov	eax, DWORD PTR _RX1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100898
	cmp	edi, 1
	jge	SHORT $L100899
$L100898:
	call	___FHranger@0
$L100899:
	mov	eax, DWORD PTR _AA2$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100900
	cmp	edi, 1
	jge	SHORT $L100901
$L100900:
	call	___FHranger@0
$L100901:
	mov	eax, DWORD PTR _RZ2$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100902
	cmp	edi, 1
	jge	SHORT $L100903
$L100902:
	call	___FHranger@0
$L100903:
	mov	eax, DWORD PTR _RY2$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100904
	cmp	edi, 1
	jge	SHORT $L100905
$L100904:
	call	___FHranger@0
$L100905:
	mov	eax, DWORD PTR _RX2$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV2
	push	eax
	call	_ENUCL@100
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 1302	; 00000516H
	mov	esi, DWORD PTR [eax]

; 1294 :      1           RX1(K),RY1(K),RZ1(K),AA1(K),                      
; 1295 :      2           RX1(L),RY1(L),RZ1(L),AA1(L),
; 1296 :      3           RX2(K),RY2(K),RZ2(K),AA2(K),
; 1297 :      4           B12(L),B12(K),
; 1298 :      5           EX1,EY1,EZ1,EN10,KLM)                       
; 1299 : C
; 1300 : C     ...Nuclear Symmetrization Integral: P_[H1,H2]
; 1301 : C
; 1302 :               CALL ENUCL(DFLAG,OV3,RX2(L),RY2(L),RZ2(L),AA2(L),

	cmp	edi, esi
	jg	SHORT $L100906
	cmp	edi, 1
	jge	SHORT $L100907
$L100906:
	call	___FHranger@0
$L100907:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98481[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100908
	cmp	edi, 1
	jge	SHORT $L100909
$L100908:
	call	___FHranger@0
$L100909:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98482[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100910
	cmp	edi, 1
	jge	SHORT $L100911
$L100910:
	call	___FHranger@0
$L100911:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98483[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100912
	cmp	edi, 1
	jge	SHORT $L100913
$L100912:
	call	___FHranger@0
$L100913:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98484[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100914
	cmp	edi, 1
	jge	SHORT $L100915
$L100914:
	call	___FHranger@0
$L100915:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98485[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100916
	cmp	edi, 1
	jge	SHORT $L100917
$L100916:
	call	___FHranger@0
$L100917:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN11
	push	OFFSET FLAT:_EZ1
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	push	OFFSET FLAT:_EY1
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98486[esp+1408]
	push	OFFSET FLAT:_EX1
	cmp	edi, esi
	jg	SHORT $L100918
	cmp	edi, 1
	jge	SHORT $L100919
$L100918:
	call	___FHranger@0
$L100919:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100920
	cmp	edi, 1
	jge	SHORT $L100921
$L100920:
	call	___FHranger@0
$L100921:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100922
	cmp	edi, 1
	jge	SHORT $L100923
$L100922:
	call	___FHranger@0
$L100923:
	mov	eax, DWORD PTR _AA1$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98481[esp+1420]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98482[esp+1428]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98483[esp+1432]
	push	eax
	jg	SHORT $L100924
	cmp	edi, 1
	jge	SHORT $L100925
$L100924:
	call	___FHranger@0
$L100925:
	mov	eax, DWORD PTR _AA1$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100926
	cmp	edi, 1
	jge	SHORT $L100927
$L100926:
	call	___FHranger@0
$L100927:
	mov	eax, DWORD PTR _RZ1$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100928
	cmp	edi, 1
	jge	SHORT $L100929
$L100928:
	call	___FHranger@0
$L100929:
	mov	eax, DWORD PTR _RY1$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100930
	cmp	edi, 1
	jge	SHORT $L100931
$L100930:
	call	___FHranger@0
$L100931:
	mov	eax, DWORD PTR _RX1$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100932
	cmp	edi, 1
	jge	SHORT $L100933
$L100932:
	call	___FHranger@0
$L100933:
	mov	eax, DWORD PTR _AA2$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98484[esp+1452]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98485[esp+1460]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98486[esp+1464]
	push	eax
	jg	SHORT $L100934
	cmp	edi, 1
	jge	SHORT $L100935
$L100934:
	call	___FHranger@0
$L100935:
	mov	eax, DWORD PTR _AA2$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100936
	cmp	edi, 1
	jge	SHORT $L100937
$L100936:
	call	___FHranger@0
$L100937:
	mov	eax, DWORD PTR _RZ2$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100938
	cmp	edi, 1
	jge	SHORT $L100939
$L100938:
	call	___FHranger@0
$L100939:
	mov	eax, DWORD PTR _RY2$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100940
	cmp	edi, 1
	jge	SHORT $L100941
$L100940:
	call	___FHranger@0
$L100941:
	mov	eax, DWORD PTR _RX2$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV3
	push	eax
	call	_ENUCL@100
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 1312	; 00000520H
	mov	esi, DWORD PTR [eax]

; 1303 :      1            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
; 1304 :      2            RX1(L),RY1(L),RZ1(L),AA1(L),
; 1305 :      3            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
; 1306 :      4            B12(L),B12(K),
; 1307 :      5            EX1,EY1,EZ1,EN11,KLM)
; 1308 : C
; 1309 : C     ...Electro/Nuclear Symmetrization Integral:
; 1310 : C         P_[elec.1,elec.2]*P_[H1,H2]
; 1311 : C
; 1312 :               CALL ENUCL(DFLAG,OV4,RX2(L),RY2(L),RZ2(L),AA2(L),

	cmp	edi, esi
	jg	SHORT $L100942
	cmp	edi, 1
	jge	SHORT $L100943
$L100942:
	call	___FHranger@0
$L100943:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98487[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100944
	cmp	edi, 1
	jge	SHORT $L100945
$L100944:
	call	___FHranger@0
$L100945:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98488[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100946
	cmp	edi, 1
	jge	SHORT $L100947
$L100946:
	call	___FHranger@0
$L100947:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98489[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100948
	cmp	edi, 1
	jge	SHORT $L100949
$L100948:
	call	___FHranger@0
$L100949:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98490[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100950
	cmp	edi, 1
	jge	SHORT $L100951
$L100950:
	call	___FHranger@0
$L100951:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98491[esp+1392]
	cmp	edi, esi
	jg	SHORT $L100952
	cmp	edi, 1
	jge	SHORT $L100953
$L100952:
	call	___FHranger@0
$L100953:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN12
	push	OFFSET FLAT:_EZ1
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	push	OFFSET FLAT:_EY1
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98492[esp+1408]
	push	OFFSET FLAT:_EX1
	cmp	edi, esi
	jg	SHORT $L100954
	cmp	edi, 1
	jge	SHORT $L100955
$L100954:
	call	___FHranger@0
$L100955:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100956
	cmp	edi, 1
	jge	SHORT $L100957
$L100956:
	call	___FHranger@0
$L100957:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100958
	cmp	edi, 1
	jge	SHORT $L100959
$L100958:
	call	___FHranger@0
$L100959:
	mov	eax, DWORD PTR _AA2$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98487[esp+1420]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98488[esp+1428]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98489[esp+1432]
	push	eax
	jg	SHORT $L100960
	cmp	edi, 1
	jge	SHORT $L100961
$L100960:
	call	___FHranger@0
$L100961:
	mov	eax, DWORD PTR _AA1$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100962
	cmp	edi, 1
	jge	SHORT $L100963
$L100962:
	call	___FHranger@0
$L100963:
	mov	eax, DWORD PTR _RZ1$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100964
	cmp	edi, 1
	jge	SHORT $L100965
$L100964:
	call	___FHranger@0
$L100965:
	mov	eax, DWORD PTR _RY1$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100966
	cmp	edi, 1
	jge	SHORT $L100967
$L100966:
	call	___FHranger@0
$L100967:
	mov	eax, DWORD PTR _RX1$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100968
	cmp	edi, 1
	jge	SHORT $L100969
$L100968:
	call	___FHranger@0
$L100969:
	mov	eax, DWORD PTR _AA1$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98490[esp+1452]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98491[esp+1460]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98492[esp+1464]
	push	eax
	jg	SHORT $L100970
	cmp	edi, 1
	jge	SHORT $L100971
$L100970:
	call	___FHranger@0
$L100971:
	mov	eax, DWORD PTR _AA2$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100972
	cmp	edi, 1
	jge	SHORT $L100973
$L100972:
	call	___FHranger@0
$L100973:
	mov	eax, DWORD PTR _RZ2$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100974
	cmp	edi, 1
	jge	SHORT $L100975
$L100974:
	call	___FHranger@0
$L100975:
	mov	eax, DWORD PTR _RY2$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100976
	cmp	edi, 1
	jge	SHORT $L100977
$L100976:
	call	___FHranger@0
$L100977:
	mov	eax, DWORD PTR _RX2$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV4
	push	eax
	call	_ENUCL@100
	mov	eax, DWORD PTR _M$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN13
	mov	edi, DWORD PTR _OVER+4
	mov	esi, DWORD PTR [eax]
	push	OFFSET FLAT:_EZ2
	push	OFFSET FLAT:_EY2
	cmp	edi, esi
	mov	DWORD PTR $T98573[esp+1420], 1320	; 00000528H

; 1313 :      1            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
; 1314 :      2            RX1(L),RY1(L),RZ1(L),AA1(L),
; 1315 :      3            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
; 1316 :      4            B12(L),B12(K),
; 1317 :      5            EX1,EY1,EZ1,EN12,KLM)
; 1318 : 	END IF
; 1319 : C 
; 1320 :               CALL ENUCL(DFLAG,OV1,RX2(L),RY2(L),RZ2(L),AA2(L),

	push	OFFSET FLAT:_EX2
	jg	SHORT $L100978
	cmp	edi, 1
	jge	SHORT $L100979
$L100978:
	call	___FHranger@0
$L100979:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100980
	cmp	edi, 1
	jge	SHORT $L100981
$L100980:
	call	___FHranger@0
$L100981:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100982
	cmp	edi, 1
	jge	SHORT $L100983
$L100982:
	call	___FHranger@0
$L100983:
	mov	eax, DWORD PTR _AA1$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100984
	cmp	edi, 1
	jge	SHORT $L100985
$L100984:
	call	___FHranger@0
$L100985:
	mov	eax, DWORD PTR _RZ1$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100986
	cmp	edi, 1
	jge	SHORT $L100987
$L100986:
	call	___FHranger@0
$L100987:
	mov	eax, DWORD PTR _RY1$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100988
	cmp	edi, 1
	jge	SHORT $L100989
$L100988:
	call	___FHranger@0
$L100989:
	mov	eax, DWORD PTR _RX1$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100990
	cmp	edi, 1
	jge	SHORT $L100991
$L100990:
	call	___FHranger@0
$L100991:
	mov	eax, DWORD PTR _AA1$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100992
	cmp	edi, 1
	jge	SHORT $L100993
$L100992:
	call	___FHranger@0
$L100993:
	mov	eax, DWORD PTR _RZ1$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100994
	cmp	edi, 1
	jge	SHORT $L100995
$L100994:
	call	___FHranger@0
$L100995:
	mov	eax, DWORD PTR _RY1$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L100996
	cmp	edi, 1
	jge	SHORT $L100997
$L100996:
	call	___FHranger@0
$L100997:
	mov	eax, DWORD PTR _RX1$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L100998
	cmp	edi, 1
	jge	SHORT $L100999
$L100998:
	call	___FHranger@0
$L100999:
	mov	eax, DWORD PTR _AA2$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101000
	cmp	edi, 1
	jge	SHORT $L101001
$L101000:
	call	___FHranger@0
$L101001:
	mov	eax, DWORD PTR _RZ2$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101002
	cmp	edi, 1
	jge	SHORT $L101003
$L101002:
	call	___FHranger@0
$L101003:
	mov	eax, DWORD PTR _RY2$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101004
	cmp	edi, 1
	jge	SHORT $L101005
$L101004:
	call	___FHranger@0
$L101005:
	mov	eax, DWORD PTR _RX2$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101006
	cmp	edi, 1
	jge	SHORT $L101007
$L101006:
	call	___FHranger@0
$L101007:
	mov	eax, DWORD PTR _AA2$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101008
	cmp	edi, 1
	jge	SHORT $L101009
$L101008:
	call	___FHranger@0
$L101009:
	mov	eax, DWORD PTR _RZ2$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101010
	cmp	edi, 1
	jge	SHORT $L101011
$L101010:
	call	___FHranger@0
$L101011:
	mov	eax, DWORD PTR _RY2$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101012
	cmp	edi, 1
	jge	SHORT $L101013
$L101012:
	call	___FHranger@0
$L101013:
	mov	eax, DWORD PTR _RX2$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV1
	push	eax
	call	_ENUCL@100
	mov	eax, DWORD PTR _M$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN14
	mov	edi, DWORD PTR _OVER+4
	mov	esi, DWORD PTR [eax]
	push	OFFSET FLAT:_EZ2
	push	OFFSET FLAT:_EY2
	cmp	edi, esi
	mov	DWORD PTR $T98573[esp+1420], 1330	; 00000532H

; 1321 :      1           RX2(K),RY2(K),RZ2(K),AA2(K),                      
; 1322 :      2           RX1(L),RY1(L),RZ1(L),AA1(L),
; 1323 :      3           RX1(K),RY1(K),RZ1(K),AA1(K),
; 1324 :      4           B12(L),B12(K),
; 1325 :      5           EX2,EY2,EZ2,EN13,KLM)                       
; 1326 : C
; 1327 : C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
; 1328 : C
; 1329 : 	IF (ISYMM) THEN
; 1330 :               CALL ENUCL(DFLAG,OV2,RX2(L),RY2(L),RZ2(L),AA2(L),

	push	OFFSET FLAT:_EX2
	jg	SHORT $L101014
	cmp	edi, 1
	jge	SHORT $L101015
$L101014:
	call	___FHranger@0
$L101015:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101016
	cmp	edi, 1
	jge	SHORT $L101017
$L101016:
	call	___FHranger@0
$L101017:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101018
	cmp	edi, 1
	jge	SHORT $L101019
$L101018:
	call	___FHranger@0
$L101019:
	mov	eax, DWORD PTR _AA2$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101020
	cmp	edi, 1
	jge	SHORT $L101021
$L101020:
	call	___FHranger@0
$L101021:
	mov	eax, DWORD PTR _RZ2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101022
	cmp	edi, 1
	jge	SHORT $L101023
$L101022:
	call	___FHranger@0
$L101023:
	mov	eax, DWORD PTR _RY2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101024
	cmp	edi, 1
	jge	SHORT $L101025
$L101024:
	call	___FHranger@0
$L101025:
	mov	eax, DWORD PTR _RX2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101026
	cmp	edi, 1
	jge	SHORT $L101027
$L101026:
	call	___FHranger@0
$L101027:
	mov	eax, DWORD PTR _AA1$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101028
	cmp	edi, 1
	jge	SHORT $L101029
$L101028:
	call	___FHranger@0
$L101029:
	mov	eax, DWORD PTR _RZ1$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101030
	cmp	edi, 1
	jge	SHORT $L101031
$L101030:
	call	___FHranger@0
$L101031:
	mov	eax, DWORD PTR _RY1$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101032
	cmp	edi, 1
	jge	SHORT $L101033
$L101032:
	call	___FHranger@0
$L101033:
	mov	eax, DWORD PTR _RX1$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101034
	cmp	edi, 1
	jge	SHORT $L101035
$L101034:
	call	___FHranger@0
$L101035:
	mov	eax, DWORD PTR _AA1$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101036
	cmp	edi, 1
	jge	SHORT $L101037
$L101036:
	call	___FHranger@0
$L101037:
	mov	eax, DWORD PTR _RZ1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101038
	cmp	edi, 1
	jge	SHORT $L101039
$L101038:
	call	___FHranger@0
$L101039:
	mov	eax, DWORD PTR _RY1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101040
	cmp	edi, 1
	jge	SHORT $L101041
$L101040:
	call	___FHranger@0
$L101041:
	mov	eax, DWORD PTR _RX1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101042
	cmp	edi, 1
	jge	SHORT $L101043
$L101042:
	call	___FHranger@0
$L101043:
	mov	eax, DWORD PTR _AA2$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101044
	cmp	edi, 1
	jge	SHORT $L101045
$L101044:
	call	___FHranger@0
$L101045:
	mov	eax, DWORD PTR _RZ2$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101046
	cmp	edi, 1
	jge	SHORT $L101047
$L101046:
	call	___FHranger@0
$L101047:
	mov	eax, DWORD PTR _RY2$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101048
	cmp	edi, 1
	jge	SHORT $L101049
$L101048:
	call	___FHranger@0
$L101049:
	mov	eax, DWORD PTR _RX2$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV2
	push	eax
	call	_ENUCL@100
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 1339	; 0000053bH
	mov	esi, DWORD PTR [eax]

; 1331 :      1           RX1(K),RY1(K),RZ1(K),AA1(K),                      
; 1332 :      2           RX1(L),RY1(L),RZ1(L),AA1(L),
; 1333 :      3           RX2(K),RY2(K),RZ2(K),AA2(K),
; 1334 :      4           B12(L),B12(K),
; 1335 :      5           EX2,EY2,EZ2,EN14,KLM)  
; 1336 : C
; 1337 : C     ...Nuclear Symmetrization Integral: P_[H1,H2]
; 1338 : C
; 1339 :               CALL ENUCL(DFLAG,OV3,RX2(L),RY2(L),RZ2(L),AA2(L),

	cmp	edi, esi
	jg	SHORT $L101050
	cmp	edi, 1
	jge	SHORT $L101051
$L101050:
	call	___FHranger@0
$L101051:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98493[esp+1392]
	cmp	edi, esi
	jg	SHORT $L101052
	cmp	edi, 1
	jge	SHORT $L101053
$L101052:
	call	___FHranger@0
$L101053:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98494[esp+1392]
	cmp	edi, esi
	jg	SHORT $L101054
	cmp	edi, 1
	jge	SHORT $L101055
$L101054:
	call	___FHranger@0
$L101055:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98495[esp+1392]
	cmp	edi, esi
	jg	SHORT $L101056
	cmp	edi, 1
	jge	SHORT $L101057
$L101056:
	call	___FHranger@0
$L101057:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98496[esp+1392]
	cmp	edi, esi
	jg	SHORT $L101058
	cmp	edi, 1
	jge	SHORT $L101059
$L101058:
	call	___FHranger@0
$L101059:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98497[esp+1392]
	cmp	edi, esi
	jg	SHORT $L101060
	cmp	edi, 1
	jge	SHORT $L101061
$L101060:
	call	___FHranger@0
$L101061:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN15
	push	OFFSET FLAT:_EZ2
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	push	OFFSET FLAT:_EY2
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98498[esp+1408]
	push	OFFSET FLAT:_EX2
	cmp	edi, esi
	jg	SHORT $L101062
	cmp	edi, 1
	jge	SHORT $L101063
$L101062:
	call	___FHranger@0
$L101063:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101064
	cmp	edi, 1
	jge	SHORT $L101065
$L101064:
	call	___FHranger@0
$L101065:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101066
	cmp	edi, 1
	jge	SHORT $L101067
$L101066:
	call	___FHranger@0
$L101067:
	mov	eax, DWORD PTR _AA1$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98493[esp+1420]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98494[esp+1428]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98495[esp+1432]
	push	eax
	jg	SHORT $L101068
	cmp	edi, 1
	jge	SHORT $L101069
$L101068:
	call	___FHranger@0
$L101069:
	mov	eax, DWORD PTR _AA1$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101070
	cmp	edi, 1
	jge	SHORT $L101071
$L101070:
	call	___FHranger@0
$L101071:
	mov	eax, DWORD PTR _RZ1$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101072
	cmp	edi, 1
	jge	SHORT $L101073
$L101072:
	call	___FHranger@0
$L101073:
	mov	eax, DWORD PTR _RY1$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101074
	cmp	edi, 1
	jge	SHORT $L101075
$L101074:
	call	___FHranger@0
$L101075:
	mov	eax, DWORD PTR _RX1$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101076
	cmp	edi, 1
	jge	SHORT $L101077
$L101076:
	call	___FHranger@0
$L101077:
	mov	eax, DWORD PTR _AA2$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98496[esp+1452]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98497[esp+1460]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98498[esp+1464]
	push	eax
	jg	SHORT $L101078
	cmp	edi, 1
	jge	SHORT $L101079
$L101078:
	call	___FHranger@0
$L101079:
	mov	eax, DWORD PTR _AA2$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101080
	cmp	edi, 1
	jge	SHORT $L101081
$L101080:
	call	___FHranger@0
$L101081:
	mov	eax, DWORD PTR _RZ2$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101082
	cmp	edi, 1
	jge	SHORT $L101083
$L101082:
	call	___FHranger@0
$L101083:
	mov	eax, DWORD PTR _RY2$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101084
	cmp	edi, 1
	jge	SHORT $L101085
$L101084:
	call	___FHranger@0
$L101085:
	mov	eax, DWORD PTR _RX2$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV3
	push	eax
	call	_ENUCL@100
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 1349	; 00000545H
	mov	esi, DWORD PTR [eax]

; 1340 :      1            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
; 1341 :      2            RX1(L),RY1(L),RZ1(L),AA1(L),
; 1342 :      3            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
; 1343 :      4            B12(L),B12(K),
; 1344 :      5            EX2,EY2,EZ2,EN15,KLM)
; 1345 : C
; 1346 : C     ...Electro/Nuclear Symmetrization Integral:
; 1347 : C         P_[elec.1,elec.2]*P_[H1,H2]
; 1348 : C
; 1349 :               CALL ENUCL(DFLAG,OV4,RX2(L),RY2(L),RZ2(L),AA2(L),

	cmp	edi, esi
	jg	SHORT $L101086
	cmp	edi, 1
	jge	SHORT $L101087
$L101086:
	call	___FHranger@0
$L101087:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98499[esp+1392]
	cmp	edi, esi
	jg	SHORT $L101088
	cmp	edi, 1
	jge	SHORT $L101089
$L101088:
	call	___FHranger@0
$L101089:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98500[esp+1392]
	cmp	edi, esi
	jg	SHORT $L101090
	cmp	edi, 1
	jge	SHORT $L101091
$L101090:
	call	___FHranger@0
$L101091:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98501[esp+1392]
	cmp	edi, esi
	jg	SHORT $L101092
	cmp	edi, 1
	jge	SHORT $L101093
$L101092:
	call	___FHranger@0
$L101093:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98502[esp+1392]
	cmp	edi, esi
	jg	SHORT $L101094
	cmp	edi, 1
	jge	SHORT $L101095
$L101094:
	call	___FHranger@0
$L101095:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98503[esp+1392]
	cmp	edi, esi
	jg	SHORT $L101096
	cmp	edi, 1
	jge	SHORT $L101097
$L101096:
	call	___FHranger@0
$L101097:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_EN16
	push	OFFSET FLAT:_EZ2
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	push	OFFSET FLAT:_EY2
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98504[esp+1408]
	push	OFFSET FLAT:_EX2
	cmp	edi, esi
	jg	SHORT $L101098
	cmp	edi, 1
	jge	SHORT $L101099
$L101098:
	call	___FHranger@0
$L101099:
	mov	eax, DWORD PTR _B12$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101100
	cmp	edi, 1
	jge	SHORT $L101101
$L101100:
	call	___FHranger@0
$L101101:
	mov	eax, DWORD PTR _B12$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101102
	cmp	edi, 1
	jge	SHORT $L101103
$L101102:
	call	___FHranger@0
$L101103:
	mov	eax, DWORD PTR _AA2$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98499[esp+1420]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98500[esp+1428]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98501[esp+1432]
	push	eax
	jg	SHORT $L101104
	cmp	edi, 1
	jge	SHORT $L101105
$L101104:
	call	___FHranger@0
$L101105:
	mov	eax, DWORD PTR _AA1$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101106
	cmp	edi, 1
	jge	SHORT $L101107
$L101106:
	call	___FHranger@0
$L101107:
	mov	eax, DWORD PTR _RZ1$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101108
	cmp	edi, 1
	jge	SHORT $L101109
$L101108:
	call	___FHranger@0
$L101109:
	mov	eax, DWORD PTR _RY1$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101110
	cmp	edi, 1
	jge	SHORT $L101111
$L101110:
	call	___FHranger@0
$L101111:
	mov	eax, DWORD PTR _RX1$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101112
	cmp	edi, 1
	jge	SHORT $L101113
$L101112:
	call	___FHranger@0
$L101113:
	mov	eax, DWORD PTR _AA1$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98502[esp+1452]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98503[esp+1460]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98504[esp+1464]
	push	eax
	jg	SHORT $L101114
	cmp	edi, 1
	jge	SHORT $L101115
$L101114:
	call	___FHranger@0
$L101115:
	mov	eax, DWORD PTR _AA2$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101116
	cmp	edi, 1
	jge	SHORT $L101117
$L101116:
	call	___FHranger@0
$L101117:
	mov	eax, DWORD PTR _RZ2$[esp+1468]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101118
	cmp	edi, 1
	jge	SHORT $L101119
$L101118:
	call	___FHranger@0
$L101119:
	mov	eax, DWORD PTR _RY2$[esp+1472]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101120
	cmp	edi, 1
	jge	SHORT $L101121
$L101120:
	call	___FHranger@0
$L101121:
	mov	eax, DWORD PTR _RX2$[esp+1476]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1476]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV4
	push	eax
	call	_ENUCL@100
	fld	QWORD PTR _EN9
	fadd	QWORD PTR _EN14
	mov	eax, DWORD PTR _M$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_ER1
	mov	edi, DWORD PTR _OVER+4
	fadd	QWORD PTR _EN6
	mov	DWORD PTR $T98573[esp+1412], 1359	; 0000054fH

; 1350 :      1            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
; 1351 :      2            RX1(L),RY1(L),RZ1(L),AA1(L),
; 1352 :      3            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
; 1353 :      4            B12(L),B12(K),
; 1354 :      5            EX2,EY2,EZ2,EN16,KLM)
; 1355 : 	END IF
; 1356 : C 
; 1357 : C
; 1358 : 	IF (ISYMM) THEN
; 1359 :               ENN1=EN1+EN2+EN3+EN4+EN5+EN6+EN7+EN8

	mov	DWORD PTR $T98573[esp+1412], 1367	; 00000557H
	fadd	QWORD PTR _EN4
	mov	esi, DWORD PTR [eax]
	cmp	edi, esi
	fadd	QWORD PTR _EN11
	fadd	QWORD PTR _EN8
	fadd	QWORD PTR _EN16
	fadd	QWORD PTR _EN13
	fadd	QWORD PTR _EN3
	fadd	QWORD PTR _EN10
	fadd	QWORD PTR _EN15
	fadd	QWORD PTR _EN7
	fadd	QWORD PTR _EN1
	fadd	QWORD PTR _EN5
	fadd	QWORD PTR _EN12
	fadd	QWORD PTR _EN2
	fstp	QWORD PTR _ENN1

; 1360 :      1            +EN9+EN10+EN11+EN12+EN13+EN14+EN15+EN16
; 1361 : 	ELSE
; 1362 :               ENN1=EN1+EN5+EN9+EN13
; 1363 : 	END IF
; 1364 : C     -----------------------------------------------------
; 1365 : C     --- Calculate the Elec./Elec. Repulsion ---
; 1366 : C     -----------------------------------------------------
; 1367 :               CALL ELREP(DFLAG,OV1,RX1(L),RY1(L),RZ1(L),AA1(L),

	jg	SHORT $L101122
	cmp	edi, 1
	jge	SHORT $L101123
$L101122:
	call	___FHranger@0
$L101123:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101124
	cmp	edi, 1
	jge	SHORT $L101125
$L101124:
	call	___FHranger@0
$L101125:
	mov	eax, DWORD PTR _B12$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101126
	cmp	edi, 1
	jge	SHORT $L101127
$L101126:
	call	___FHranger@0
$L101127:
	mov	eax, DWORD PTR _AA2$[esp+1404]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101128
	cmp	edi, 1
	jge	SHORT $L101129
$L101128:
	call	___FHranger@0
$L101129:
	mov	eax, DWORD PTR _RZ2$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101130
	cmp	edi, 1
	jge	SHORT $L101131
$L101130:
	call	___FHranger@0
$L101131:
	mov	eax, DWORD PTR _RY2$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101132
	cmp	edi, 1
	jge	SHORT $L101133
$L101132:
	call	___FHranger@0
$L101133:
	mov	eax, DWORD PTR _RX2$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101134
	cmp	edi, 1
	jge	SHORT $L101135
$L101134:
	call	___FHranger@0
$L101135:
	mov	eax, DWORD PTR _AA2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101136
	cmp	edi, 1
	jge	SHORT $L101137
$L101136:
	call	___FHranger@0
$L101137:
	mov	eax, DWORD PTR _RZ2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101138
	cmp	edi, 1
	jge	SHORT $L101139
$L101138:
	call	___FHranger@0
$L101139:
	mov	eax, DWORD PTR _RY2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101140
	cmp	edi, 1
	jge	SHORT $L101141
$L101140:
	call	___FHranger@0
$L101141:
	mov	eax, DWORD PTR _RX2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101142
	cmp	edi, 1
	jge	SHORT $L101143
$L101142:
	call	___FHranger@0
$L101143:
	mov	eax, DWORD PTR _AA1$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101144
	cmp	edi, 1
	jge	SHORT $L101145
$L101144:
	call	___FHranger@0
$L101145:
	mov	eax, DWORD PTR _RZ1$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101146
	cmp	edi, 1
	jge	SHORT $L101147
$L101146:
	call	___FHranger@0
$L101147:
	mov	eax, DWORD PTR _RY1$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101148
	cmp	edi, 1
	jge	SHORT $L101149
$L101148:
	call	___FHranger@0
$L101149:
	mov	eax, DWORD PTR _RX1$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101150
	cmp	edi, 1
	jge	SHORT $L101151
$L101150:
	call	___FHranger@0
$L101151:
	mov	eax, DWORD PTR _AA1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101152
	cmp	edi, 1
	jge	SHORT $L101153
$L101152:
	call	___FHranger@0
$L101153:
	mov	eax, DWORD PTR _RZ1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101154
	cmp	edi, 1
	jge	SHORT $L101155
$L101154:
	call	___FHranger@0
$L101155:
	mov	eax, DWORD PTR _RY1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101156
	cmp	edi, 1
	jge	SHORT $L101157
$L101156:
	call	___FHranger@0
$L101157:
	mov	eax, DWORD PTR _RX1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1464]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV1
	push	eax
	call	_ELREP@88
	mov	eax, DWORD PTR _M$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_ER2
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1412], 1376	; 00000560H
	mov	esi, DWORD PTR [eax]

; 1368 :      1           RX1(K),RY1(K),RZ1(K),AA1(K),                      
; 1369 :      2           RX2(L),RY2(L),RZ2(L),AA2(L),
; 1370 :      3           RX2(K),RY2(K),RZ2(K),AA2(K),
; 1371 :      4           B12(L),B12(K),ER1,KLM)  
; 1372 : C
; 1373 : C     ...Electron Symmetrization Integral: P_[elec.1,elec.2]
; 1374 : C
; 1375 : 	IF (ISYMM) THEN
; 1376 :               CALL ELREP(DFLAG,OV2,RX1(L),RY1(L),RZ1(L),AA1(L),

	cmp	edi, esi
	jg	SHORT $L101158
	cmp	edi, 1
	jge	SHORT $L101159
$L101158:
	call	___FHranger@0
$L101159:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101160
	cmp	edi, 1
	jge	SHORT $L101161
$L101160:
	call	___FHranger@0
$L101161:
	mov	eax, DWORD PTR _B12$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101162
	cmp	edi, 1
	jge	SHORT $L101163
$L101162:
	call	___FHranger@0
$L101163:
	mov	eax, DWORD PTR _AA1$[esp+1404]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101164
	cmp	edi, 1
	jge	SHORT $L101165
$L101164:
	call	___FHranger@0
$L101165:
	mov	eax, DWORD PTR _RZ1$[esp+1408]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101166
	cmp	edi, 1
	jge	SHORT $L101167
$L101166:
	call	___FHranger@0
$L101167:
	mov	eax, DWORD PTR _RY1$[esp+1412]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101168
	cmp	edi, 1
	jge	SHORT $L101169
$L101168:
	call	___FHranger@0
$L101169:
	mov	eax, DWORD PTR _RX1$[esp+1416]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101170
	cmp	edi, 1
	jge	SHORT $L101171
$L101170:
	call	___FHranger@0
$L101171:
	mov	eax, DWORD PTR _AA2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101172
	cmp	edi, 1
	jge	SHORT $L101173
$L101172:
	call	___FHranger@0
$L101173:
	mov	eax, DWORD PTR _RZ2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101174
	cmp	edi, 1
	jge	SHORT $L101175
$L101174:
	call	___FHranger@0
$L101175:
	mov	eax, DWORD PTR _RY2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101176
	cmp	edi, 1
	jge	SHORT $L101177
$L101176:
	call	___FHranger@0
$L101177:
	mov	eax, DWORD PTR _RX2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101178
	cmp	edi, 1
	jge	SHORT $L101179
$L101178:
	call	___FHranger@0
$L101179:
	mov	eax, DWORD PTR _AA2$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101180
	cmp	edi, 1
	jge	SHORT $L101181
$L101180:
	call	___FHranger@0
$L101181:
	mov	eax, DWORD PTR _RZ2$[esp+1440]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101182
	cmp	edi, 1
	jge	SHORT $L101183
$L101182:
	call	___FHranger@0
$L101183:
	mov	eax, DWORD PTR _RY2$[esp+1444]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101184
	cmp	edi, 1
	jge	SHORT $L101185
$L101184:
	call	___FHranger@0
$L101185:
	mov	eax, DWORD PTR _RX2$[esp+1448]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101186
	cmp	edi, 1
	jge	SHORT $L101187
$L101186:
	call	___FHranger@0
$L101187:
	mov	eax, DWORD PTR _AA1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101188
	cmp	edi, 1
	jge	SHORT $L101189
$L101188:
	call	___FHranger@0
$L101189:
	mov	eax, DWORD PTR _RZ1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101190
	cmp	edi, 1
	jge	SHORT $L101191
$L101190:
	call	___FHranger@0
$L101191:
	mov	eax, DWORD PTR _RY1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101192
	cmp	edi, 1
	jge	SHORT $L101193
$L101192:
	call	___FHranger@0
$L101193:
	mov	eax, DWORD PTR _RX1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1464]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV2
	push	eax
	call	_ELREP@88
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 1384	; 00000568H
	mov	esi, DWORD PTR [eax]

; 1377 :      1           RX2(K),RY2(K),RZ2(K),AA2(K),                      
; 1378 :      2           RX2(L),RY2(L),RZ2(L),AA2(L),
; 1379 :      3           RX1(K),RY1(K),RZ1(K),AA1(K),
; 1380 :      4           B12(L),B12(K),ER2,KLM)  
; 1381 : C
; 1382 : C     ...Nuclear Symmetrization Overlap Integral: P_[H1,H2]
; 1383 : C
; 1384 :               CALL ELREP(DFLAG,OV3,RX1(L),RY1(L),RZ1(L),AA1(L),

	cmp	edi, esi
	jg	SHORT $L101194
	cmp	edi, 1
	jge	SHORT $L101195
$L101194:
	call	___FHranger@0
$L101195:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98505[esp+1392]
	cmp	edi, esi
	jg	SHORT $L101196
	cmp	edi, 1
	jge	SHORT $L101197
$L101196:
	call	___FHranger@0
$L101197:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98506[esp+1392]
	cmp	edi, esi
	jg	SHORT $L101198
	cmp	edi, 1
	jge	SHORT $L101199
$L101198:
	call	___FHranger@0
$L101199:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98507[esp+1392]
	cmp	edi, esi
	jg	SHORT $L101200
	cmp	edi, 1
	jge	SHORT $L101201
$L101200:
	call	___FHranger@0
$L101201:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98508[esp+1392]
	cmp	edi, esi
	jg	SHORT $L101202
	cmp	edi, 1
	jge	SHORT $L101203
$L101202:
	call	___FHranger@0
$L101203:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98509[esp+1392]
	cmp	edi, esi
	jg	SHORT $L101204
	cmp	edi, 1
	jge	SHORT $L101205
$L101204:
	call	___FHranger@0
$L101205:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_ER3
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98510[esp+1400]
	cmp	edi, esi
	jg	SHORT $L101206
	cmp	edi, 1
	jge	SHORT $L101207
$L101206:
	call	___FHranger@0
$L101207:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101208
	cmp	edi, 1
	jge	SHORT $L101209
$L101208:
	call	___FHranger@0
$L101209:
	mov	eax, DWORD PTR _B12$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101210
	cmp	edi, 1
	jge	SHORT $L101211
$L101210:
	call	___FHranger@0
$L101211:
	mov	eax, DWORD PTR _AA2$[esp+1404]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98505[esp+1408]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98506[esp+1416]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98507[esp+1420]
	push	eax
	jg	SHORT $L101212
	cmp	edi, 1
	jge	SHORT $L101213
$L101212:
	call	___FHranger@0
$L101213:
	mov	eax, DWORD PTR _AA2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101214
	cmp	edi, 1
	jge	SHORT $L101215
$L101214:
	call	___FHranger@0
$L101215:
	mov	eax, DWORD PTR _RZ2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101216
	cmp	edi, 1
	jge	SHORT $L101217
$L101216:
	call	___FHranger@0
$L101217:
	mov	eax, DWORD PTR _RY2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101218
	cmp	edi, 1
	jge	SHORT $L101219
$L101218:
	call	___FHranger@0
$L101219:
	mov	eax, DWORD PTR _RX2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101220
	cmp	edi, 1
	jge	SHORT $L101221
$L101220:
	call	___FHranger@0
$L101221:
	mov	eax, DWORD PTR _AA1$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98508[esp+1440]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98509[esp+1448]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98510[esp+1452]
	push	eax
	jg	SHORT $L101222
	cmp	edi, 1
	jge	SHORT $L101223
$L101222:
	call	___FHranger@0
$L101223:
	mov	eax, DWORD PTR _AA1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101224
	cmp	edi, 1
	jge	SHORT $L101225
$L101224:
	call	___FHranger@0
$L101225:
	mov	eax, DWORD PTR _RZ1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101226
	cmp	edi, 1
	jge	SHORT $L101227
$L101226:
	call	___FHranger@0
$L101227:
	mov	eax, DWORD PTR _RY1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101228
	cmp	edi, 1
	jge	SHORT $L101229
$L101228:
	call	___FHranger@0
$L101229:
	mov	eax, DWORD PTR _RX1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1464]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV3
	push	eax
	call	_ELREP@88
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	edi, DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 1393	; 00000571H
	mov	esi, DWORD PTR [eax]

; 1385 :      1            -RX1(K),-RY1(K),-RZ1(K),AA1(K),                      
; 1386 :      2            RX2(L),RY2(L),RZ2(L),AA2(L),
; 1387 :      3            -RX2(K),-RY2(K),-RZ2(K),AA2(K),
; 1388 :      4            B12(L),B12(K),ER3,KLM)
; 1389 : C
; 1390 : C     ...Electro/Nuclear Symmetrization Overlap Integral:
; 1391 : C         P_[elec.1,elec.2]*P_[H1,H2]
; 1392 : C
; 1393 :               CALL ELREP(DFLAG,OV4,RX1(L),RY1(L),RZ1(L),AA1(L),

	cmp	edi, esi
	jg	SHORT $L101230
	cmp	edi, 1
	jge	SHORT $L101231
$L101230:
	call	___FHranger@0
$L101231:
	mov	eax, DWORD PTR _RZ1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98511[esp+1392]
	cmp	edi, esi
	jg	SHORT $L101232
	cmp	edi, 1
	jge	SHORT $L101233
$L101232:
	call	___FHranger@0
$L101233:
	mov	eax, DWORD PTR _RY1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98512[esp+1392]
	cmp	edi, esi
	jg	SHORT $L101234
	cmp	edi, 1
	jge	SHORT $L101235
$L101234:
	call	___FHranger@0
$L101235:
	mov	eax, DWORD PTR _RX1$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98513[esp+1392]
	cmp	edi, esi
	jg	SHORT $L101236
	cmp	edi, 1
	jge	SHORT $L101237
$L101236:
	call	___FHranger@0
$L101237:
	mov	eax, DWORD PTR _RZ2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98514[esp+1392]
	cmp	edi, esi
	jg	SHORT $L101238
	cmp	edi, 1
	jge	SHORT $L101239
$L101238:
	call	___FHranger@0
$L101239:
	mov	eax, DWORD PTR _RY2$[esp+1388]
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98515[esp+1392]
	cmp	edi, esi
	jg	SHORT $L101240
	cmp	edi, 1
	jge	SHORT $L101241
$L101240:
	call	___FHranger@0
$L101241:
	mov	eax, DWORD PTR _RX2$[esp+1388]
	push	OFFSET FLAT:_KLM
	push	OFFSET FLAT:_ER4
	fld	QWORD PTR [eax+edi*8-8]
	fchs
	mov	edi, DWORD PTR _OVER+4
	fstp	QWORD PTR $T98516[esp+1400]
	cmp	edi, esi
	jg	SHORT $L101242
	cmp	edi, 1
	jge	SHORT $L101243
$L101242:
	call	___FHranger@0
$L101243:
	mov	eax, DWORD PTR _B12$[esp+1396]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101244
	cmp	edi, 1
	jge	SHORT $L101245
$L101244:
	call	___FHranger@0
$L101245:
	mov	eax, DWORD PTR _B12$[esp+1400]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101246
	cmp	edi, 1
	jge	SHORT $L101247
$L101246:
	call	___FHranger@0
$L101247:
	mov	eax, DWORD PTR _AA1$[esp+1404]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98511[esp+1408]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98512[esp+1416]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98513[esp+1420]
	push	eax
	jg	SHORT $L101248
	cmp	edi, 1
	jge	SHORT $L101249
$L101248:
	call	___FHranger@0
$L101249:
	mov	eax, DWORD PTR _AA2$[esp+1420]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101250
	cmp	edi, 1
	jge	SHORT $L101251
$L101250:
	call	___FHranger@0
$L101251:
	mov	eax, DWORD PTR _RZ2$[esp+1424]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101252
	cmp	edi, 1
	jge	SHORT $L101253
$L101252:
	call	___FHranger@0
$L101253:
	mov	eax, DWORD PTR _RY2$[esp+1428]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101254
	cmp	edi, 1
	jge	SHORT $L101255
$L101254:
	call	___FHranger@0
$L101255:
	mov	eax, DWORD PTR _RX2$[esp+1432]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER+4
	cmp	edi, esi
	jg	SHORT $L101256
	cmp	edi, 1
	jge	SHORT $L101257
$L101256:
	call	___FHranger@0
$L101257:
	mov	eax, DWORD PTR _AA2$[esp+1436]
	lea	edi, DWORD PTR [eax+edi*8]
	lea	eax, DWORD PTR $T98514[esp+1440]
	sub	edi, 8
	push	edi
	push	eax
	lea	eax, DWORD PTR $T98515[esp+1448]
	mov	edi, DWORD PTR _OVER
	push	eax
	cmp	edi, esi
	lea	eax, DWORD PTR $T98516[esp+1452]
	push	eax
	jg	SHORT $L101258
	cmp	edi, 1
	jge	SHORT $L101259
$L101258:
	call	___FHranger@0
$L101259:
	mov	eax, DWORD PTR _AA1$[esp+1452]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101260
	cmp	edi, 1
	jge	SHORT $L101261
$L101260:
	call	___FHranger@0
$L101261:
	mov	eax, DWORD PTR _RZ1$[esp+1456]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101262
	cmp	edi, 1
	jge	SHORT $L101263
$L101262:
	call	___FHranger@0
$L101263:
	mov	eax, DWORD PTR _RY1$[esp+1460]
	lea	edi, DWORD PTR [eax+edi*8]
	sub	edi, 8
	push	edi
	mov	edi, DWORD PTR _OVER
	cmp	edi, esi
	jg	SHORT $L101264
	cmp	edi, 1
	jge	SHORT $L101265
$L101264:
	call	___FHranger@0
$L101265:
	mov	eax, DWORD PTR _RX1$[esp+1464]
	lea	edi, DWORD PTR [eax+edi*8]
	mov	eax, DWORD PTR _DFLAG$[esp+1464]
	sub	edi, 8
	push	edi
	push	OFFSET FLAT:_OV4
	push	eax
	call	_ELREP@88
	fld	QWORD PTR _ER2
	fadd	QWORD PTR _ER4
	fld	QWORD PTR _EN9
	fadd	QWORD PTR _EN14
	fxch	ST(1)
	fadd	QWORD PTR _ER3
	fxch	ST(1)
	fadd	QWORD PTR _EN6
	fxch	ST(1)
	fadd	QWORD PTR _ER1
	fxch	ST(1)
	fadd	QWORD PTR _EN4
	mov	DWORD PTR $T98573[esp+1404], 1402	; 0000057aH
	mov	DWORD PTR $T98573[esp+1404], 1409	; 00000581H
	fadd	QWORD PTR _EN11
	fxch	ST(1)
	fstp	QWORD PTR -1360+[esp+1392]

; 1394 :      1            -RX2(K),-RY2(K),-RZ2(K),AA2(K),                      
; 1395 :      2            RX2(L),RY2(L),RZ2(L),AA2(L),
; 1396 :      3            -RX1(K),-RY1(K),-RZ1(K),AA1(K),
; 1397 :      4            B12(L),B12(K),ER4,KLM)
; 1398 : 	END IF
; 1399 : 
; 1400 : C              
; 1401 : 	IF (ISYMM) THEN
; 1402 : 		EER1=ER1+ER2+ER3+ER4

	fadd	QWORD PTR _EN8
	mov	eax, DWORD PTR -1360+[esp+1396]
	mov	edi, DWORD PTR _OVER
	mov	DWORD PTR $T98573[esp+1404], 1429	; 00000595H
	fadd	QWORD PTR _EN16
	mov	DWORD PTR $T98573[esp+1404], 1447	; 000005a7H
	mov	DWORD PTR _EER1+4, eax
	mov	eax, DWORD PTR -1360+[esp+1392]
	fadd	QWORD PTR _EN13
	mov	DWORD PTR _EER1, eax

; 1403 : 	ELSE
; 1404 : 		EER1=ER1
; 1405 : 	END IF
; 1406 : C
; 1407 : C
; 1408 : 	IF (ISYMM) THEN
; 1409 :               EN=EN1+EN2+EN3+EN4+EN5+EN6+EN7+EN8

	mov	eax, DWORD PTR -1360+[esp+1396]
	mov	DWORD PTR _ER+4, eax
	mov	eax, DWORD PTR -1360+[esp+1392]
	fadd	QWORD PTR _EN3
	mov	DWORD PTR _ER, eax
	mov	eax, DWORD PTR _M$[esp+1388]
	fadd	QWORD PTR _EN10
	mov	eax, DWORD PTR [eax]
	mov	DWORD PTR -1376+[esp+1392], eax
	fadd	QWORD PTR _EN15
	cmp	edi, DWORD PTR -1376+[esp+1392]
	fadd	QWORD PTR _EN7
	fadd	QWORD PTR _EN1
	fadd	QWORD PTR _EN5
	fadd	QWORD PTR _EN12
	fadd	QWORD PTR _EN2
	fstp	QWORD PTR _EN

; 1410 :      1            +EN9+EN10+EN11+EN12+EN13+EN14+EN15+EN16
; 1411 : 	ELSE
; 1412 :               EN=EN1+EN5+EN9+EN13
; 1413 : 	END IF
; 1414 : C
; 1415 : 	IF(IPERM.EQ.1) THEN
; 1416 :               EN=EN1+EN5+EN9+EN13
; 1417 : 	END IF
; 1418 : 	IF(IPERM.EQ.2) THEN
; 1419 :               EN=EN2+EN6+EN10+EN14
; 1420 : 	END IF
; 1421 : 	IF(IPERM.EQ.3) THEN
; 1422 :               EN=EN3+EN7+EN11+EN15
; 1423 : 	END IF
; 1424 : 	IF(IPERM.EQ.4) THEN
; 1425 :               EN=EN4+EN8+EN12+EN16
; 1426 : 	END IF
; 1427 : C
; 1428 : 	IF (ISYMM) THEN
; 1429 :               ER=ER1+ER2+ER3+ER4
; 1430 : 	ELSE
; 1431 : 			ER=ER1
; 1432 : 	END IF
; 1433 : C
; 1434 : 	IF(IPERM.EQ.1) THEN
; 1435 : 			ER=ER1
; 1436 : 	END IF
; 1437 : 	IF(IPERM.EQ.2) THEN
; 1438 : 			ER=ER2
; 1439 : 	END IF
; 1440 : 	IF(IPERM.EQ.3) THEN
; 1441 : 			ER=ER3
; 1442 : 	END IF
; 1443 : 	IF(IPERM.EQ.4) THEN
; 1444 : 			ER=ER4
; 1445 : 	END IF
; 1446 : C
; 1447 :               H(L,K)=ER+T-EN+(1.0D+00/1.4011D+00)*S(L,K)

	jg	SHORT $L101266
	cmp	edi, 1
	jge	SHORT $L101267
$L101266:
	call	___FHranger@0
$L101267:
	mov	esi, DWORD PTR _OVER+4
	cmp	esi, DWORD PTR -1376+[esp+1392]
	jg	SHORT $L101268
	cmp	esi, 1
	jge	SHORT $L101269
$L101268:
	call	___FHranger@0
$L101269:
	mov	ebx, DWORD PTR _OVER+4
	cmp	ebx, DWORD PTR -1376+[esp+1392]
	jg	SHORT $L101270
	cmp	ebx, 1
	jge	SHORT $L101271
$L101270:
	call	___FHranger@0
$L101271:
	imul	ebx, DWORD PTR _7S$[esp+1428]
	mov	ebp, DWORD PTR _OVER
	cmp	ebp, DWORD PTR -1376+[esp+1392]
	jg	SHORT $L101272
	cmp	ebp, 1
	jge	SHORT $L101273
$L101272:
	call	___FHranger@0
$L101273:
	add	ebx, ebp
	mov	eax, DWORD PTR _S$[esp+1388]
	add	ebx, DWORD PTR _7S$[esp+1400]
	mov	ecx, DWORD PTR _7H$[esp+1428]
	imul	ecx, esi
	fld	QWORD PTR [eax+ebx*8]
	fmul	QWORD PTR $T99768
	add	ecx, DWORD PTR _7H$[esp+1400]
	mov	eax, DWORD PTR _H$[esp+1388]
	fsub	QWORD PTR _EN
	add	ecx, edi
	push	0
	mov	DWORD PTR $T98573[esp+1408], 1449	; 000005a9H
	fadd	QWORD PTR _T
	push	0
	push	0
	push	0
	push	0
	fadd	QWORD PTR -1360+[esp+1412]
	push	6
	fstp	QWORD PTR [eax+ecx*8]

; 1448 : C
; 1449 : 	WRITE(6,*)'ER=',ER,L,K

	call	_f90_write
	add	esp, 24					; 00000018H
	push	3
	push	OFFSET FLAT:_.SG633_
	call	_lwrite_ch
	add	esp, 8
	mov	eax, DWORD PTR _ER+4
	push	eax
	mov	eax, DWORD PTR _ER
	push	eax
	call	_lwrite_d
	add	esp, 8
	mov	eax, DWORD PTR _OVER
	push	eax
	call	_lwrite_i
	add	esp, 4
	mov	eax, DWORD PTR _OVER+4
	push	eax
	call	_lwrite_i
	add	esp, 4
	call	_f90_endwrite
	mov	DWORD PTR $T98573[esp+1404], 1450	; 000005aaH

; 1450 : 	WRITE(6,*)'EN=',EN,L,K

	push	0
	push	0
	push	0
	push	0
	push	0
	push	6
	call	_f90_write
	add	esp, 24					; 00000018H
	push	3
	push	OFFSET FLAT:_.SG635_
	call	_lwrite_ch
	add	esp, 8
	mov	eax, DWORD PTR _EN+4
	push	eax
	mov	eax, DWORD PTR _EN
	push	eax
	call	_lwrite_d
	add	esp, 8
	mov	eax, DWORD PTR _OVER
	push	eax
	call	_lwrite_i
	add	esp, 4
	mov	eax, DWORD PTR _OVER+4
	push	eax
	call	_lwrite_i
	add	esp, 4
	call	_f90_endwrite
	dec	DWORD PTR -1372+[esp+1392]
	inc	DWORD PTR _OVER+4
	mov	DWORD PTR $T98573[esp+1404], 1053	; 0000041dH

; 1053 :           DO K=1,L

	cmp	DWORD PTR -1372+[esp+1392], 0
	jne	$L98420
	mov	DWORD PTR $T98573[esp+1404], 1052	; 0000041cH
$L565:

; 1049 : C
; 1050 : C
; 1051 : C
; 1052 :       DO L=1,M

	mov	eax, DWORD PTR -1364+[esp+1392]
	mov	ecx, DWORD PTR _OVER
	dec	eax
	inc	ecx
	mov	DWORD PTR -1364+[esp+1392], eax
	mov	DWORD PTR _OVER, ecx
	mov	DWORD PTR $T98573[esp+1404], 1052	; 0000041cH
	test	eax, eax
	jne	$L98419
$L564:

; 1451 : C     ----------------------------------------------------- 
; 1452 :           END DO 
; 1453 :       END DO

	mov	DWORD PTR $T98573[esp+1404], 1453	; 000005adH
	mov	eax, DWORD PTR _M$[esp+1388]
	mov	DWORD PTR _OVER, 1
	mov	ecx, DWORD PTR [eax]
	mov	DWORD PTR -608+[esp+1392], ecx
	mov	ecx, 1454				; 000005aeH
	mov	eax, DWORD PTR -608+[esp+1392]
	mov	DWORD PTR $T98573[esp+1404], ecx

; 1454 : 	DO L=1,M

	mov	DWORD PTR -592+[esp+1392], eax
	test	eax, eax
	jle	$L636
$L98517:
	mov	ecx, 1455				; 000005afH
	mov	eax, DWORD PTR _OVER
	mov	DWORD PTR $T98573[esp+1404], ecx
	mov	DWORD PTR -600+[esp+1392], eax

; 1455 : 		DO K=1,L

	mov	DWORD PTR _OVER+4, 1
	test	eax, eax
	jle	$L637
$L98518:
	mov	DWORD PTR $T98573[esp+1404], 1456	; 000005b0H

; 1456 : 			IF (K .NE. L) THEN

	mov	eax, DWORD PTR _OVER
	cmp	DWORD PTR _OVER+4, eax
	je	$L639
	mov	ebx, DWORD PTR _OVER+4
	mov	ecx, DWORD PTR -608+[esp+1392]
	mov	DWORD PTR $T98573[esp+1404], 1457	; 000005b1H

; 1457 : 				S(K,L)=S(L,K)

	cmp	ebx, ecx
	jg	SHORT $L101274
	cmp	ebx, 1
	jge	SHORT $L101275
$L101274:
	call	___FHranger@0
$L101275:
	mov	esi, DWORD PTR _OVER
	mov	eax, DWORD PTR -608+[esp+1392]
	cmp	esi, eax
	jg	SHORT $L101276
	cmp	esi, 1
	jge	SHORT $L101277
$L101276:
	call	___FHranger@0
$L101277:
	mov	edi, DWORD PTR _OVER+4
	mov	eax, DWORD PTR -608+[esp+1392]
	cmp	edi, eax
	jg	SHORT $L101278
	cmp	edi, 1
	jge	SHORT $L101279
$L101278:
	call	___FHranger@0
$L101279:
	mov	ebp, DWORD PTR _OVER
	mov	eax, DWORD PTR -608+[esp+1392]
	imul	edi, DWORD PTR _7S$[esp+1428]
	cmp	ebp, eax
	jg	SHORT $L101280
	cmp	ebp, 1
	jge	SHORT $L101281
$L101280:
	call	___FHranger@0
$L101281:
	add	edi, ebp
	mov	eax, DWORD PTR _7S$[esp+1400]
	add	edi, eax
	mov	edx, DWORD PTR _7S$[esp+1428]
	imul	edx, esi
	add	edx, eax
	mov	eax, DWORD PTR -608+[esp+1392]
	add	edx, ebx
	mov	ebx, DWORD PTR _S$[esp+1388]
	mov	DWORD PTR $T98573[esp+1404], 1458	; 000005b2H
	mov	ecx, DWORD PTR [ebx+edi*8+4]
	mov	DWORD PTR [ebx+edx*8+4], ecx
	mov	ecx, DWORD PTR [ebx+edi*8]
	mov	DWORD PTR [ebx+edx*8], ecx

; 1458 : 				H(K,L)=H(L,K)

	mov	ebx, DWORD PTR _OVER+4
	cmp	ebx, eax
	jg	SHORT $L101282
	cmp	ebx, 1
	jge	SHORT $L101283
$L101282:
	call	___FHranger@0
$L101283:
	mov	esi, DWORD PTR _OVER
	mov	eax, DWORD PTR -608+[esp+1392]
	cmp	esi, eax
	jg	SHORT $L101284
	cmp	esi, 1
	jge	SHORT $L101285
$L101284:
	call	___FHranger@0
$L101285:
	mov	edi, DWORD PTR _OVER+4
	mov	eax, DWORD PTR -608+[esp+1392]
	cmp	edi, eax
	jg	SHORT $L101286
	cmp	edi, 1
	jge	SHORT $L101287
$L101286:
	call	___FHranger@0
$L101287:
	mov	ebp, DWORD PTR _OVER
	mov	eax, DWORD PTR -608+[esp+1392]
	imul	edi, DWORD PTR _7H$[esp+1428]
	cmp	ebp, eax
	jg	SHORT $L101288
	cmp	ebp, 1
	jge	SHORT $L101289
$L101288:
	call	___FHranger@0
$L101289:
	add	edi, ebp
	mov	eax, DWORD PTR _7H$[esp+1400]
	add	edi, eax
	mov	edx, DWORD PTR _7H$[esp+1428]
	imul	edx, esi
	add	edx, eax
	add	edx, ebx
	mov	ebx, DWORD PTR _H$[esp+1388]
	mov	DWORD PTR $T98573[esp+1404], 1455	; 000005afH
	mov	ecx, DWORD PTR [ebx+edi*8+4]
	mov	DWORD PTR [ebx+edx*8+4], ecx
	mov	ecx, DWORD PTR [ebx+edi*8]
	mov	DWORD PTR [ebx+edx*8], ecx
$L639:

; 1455 : 		DO K=1,L

	mov	eax, DWORD PTR -600+[esp+1392]
	mov	ecx, DWORD PTR _OVER+4
	dec	eax
	inc	ecx
	mov	DWORD PTR -600+[esp+1392], eax
	mov	DWORD PTR _OVER+4, ecx
	mov	DWORD PTR $T98573[esp+1404], 1455	; 000005afH
	test	eax, eax
	jne	$L98518
	mov	DWORD PTR $T98573[esp+1404], 1454	; 000005aeH
$L637:

; 1454 : 	DO L=1,M

	mov	eax, DWORD PTR -592+[esp+1392]
	mov	ecx, DWORD PTR _OVER
	dec	eax
	inc	ecx
	mov	DWORD PTR -592+[esp+1392], eax
	mov	DWORD PTR _OVER, ecx
	mov	DWORD PTR $T98573[esp+1404], 1454	; 000005aeH
	test	eax, eax
	jne	$L98517
	mov	DWORD PTR $T98573[esp+1404], 1461	; 000005b5H
$L636:

; 1459 : 			END IF
; 1460 : 		END DO
; 1461 : 	END DO

	mov	ecx, 1465				; 000005b9H
	mov	eax, DWORD PTR _DFLAG$[esp+1388]
	mov	DWORD PTR $T98573[esp+1404], ecx

; 1462 : C     ----------------------------------------------------- 
; 1463 : C      CALL TAB(S,M,M,MMAX,MMAX)                  
; 1464 : C      CALL TAB(H,M,M,MMAX,MMAX)
; 1465 : 	IF (DFLAG .NE. 3) CALL ENERGY(M,H,S,EMIN,CO,KIM)

	mov	edx, DWORD PTR [eax]
	cmp	edx, 3
	je	SHORT $L642
	mov	eax, DWORD PTR _CO$[esp+1388]
	push	OFFSET FLAT:_KIM
	mov	ecx, DWORD PTR _EMIN$[esp+1392]
	push	eax
	mov	edx, DWORD PTR _S$[esp+1396]
	push	ecx
	mov	eax, DWORD PTR _H$[esp+1400]
	push	edx
	mov	ecx, DWORD PTR _M$[esp+1404]
	push	eax
	push	ecx
	call	_ENERGY@24
$L642:

; 1466 : C		STOP
; 1467 : C     -----------------------------------------------------      
; 1468 : 	IF (DFLAG .EQ. 3) THEN

	mov	ecx, DWORD PTR _DFLAG$[esp+1388]
	mov	DWORD PTR $T98573[esp+1404], 1468	; 000005bcH
	cmp	DWORD PTR [ecx], 3
	jne	$L646
	mov	DWORD PTR $T98573[esp+1404], 1469	; 000005bdH

; 1469 :       SUMH=0.00D+00

	xor	ecx, ecx
	mov	edx, DWORD PTR _M$[esp+1388]
	mov	eax, 1471				; 000005bfH
	mov	DWORD PTR _SUMH, ecx
	mov	edi, DWORD PTR [edx]
	mov	DWORD PTR _SUMH+4, ecx
	mov	esi, 1
	mov	DWORD PTR _SUMS, ecx
	mov	DWORD PTR _SUMS+4, ecx
	mov	DWORD PTR -592+[esp+1392], edi
	mov	DWORD PTR $T98573[esp+1404], 1470	; 000005beH

; 1470 : 	SUMS=0.00D+00

	mov	DWORD PTR $T98573[esp+1404], eax
	cmp	edi, ecx

; 1471 : 	DO J=1,M

	jle	$L648
$L98519:
	mov	ecx, 1472				; 000005c0H
	mov	ebx, 1
	mov	DWORD PTR $T98573[esp+1404], ecx
	mov	DWORD PTR -600+[esp+1392], edi

; 1472 : 		DO I=1,M

	test	edi, edi
	jle	$L649
$L98520:
	cmp	esi, edi
	mov	ebp, esi
	mov	DWORD PTR $T98573[esp+1404], 1473	; 000005c1H

; 1473 : 			SUMH=SUMH+EVEC(I)*H(I,J)*EVEC(J)

	jg	SHORT $L101290
	cmp	ebp, 1
	jge	SHORT $L101291
$L101290:
	call	___FHranger@0
$L101291:
	imul	ebp, DWORD PTR _7H$[esp+1428]
	cmp	ebx, edi
	jg	SHORT $L101292
	cmp	ebx, 1
	jge	SHORT $L101293
$L101292:
	call	___FHranger@0
$L101293:
	add	ebp, ebx
	mov	eax, DWORD PTR _H$[esp+1388]
	add	ebp, DWORD PTR _7H$[esp+1400]
	cmp	esi, edi
	fld	QWORD PTR [eax+ebp*8]
	jg	SHORT $L101294
	cmp	esi, 1
	jge	SHORT $L101295
$L101294:
	call	___FHranger@0
$L101295:
	mov	eax, DWORD PTR _EVEC$[esp+1388]
	cmp	ebx, edi
	fmul	QWORD PTR [eax+esi*8-8]
	jg	SHORT $L101296
	cmp	ebx, 1
	jge	SHORT $L101297
$L101296:
	call	___FHranger@0
$L101297:
	cmp	esi, edi
	mov	eax, DWORD PTR _EVEC$[esp+1388]
	mov	DWORD PTR $T98573[esp+1404], 1474	; 000005c2H
	fmul	QWORD PTR [eax+ebx*8-8]
	fadd	QWORD PTR _SUMH
	fstp	QWORD PTR _SUMH

; 1474 : 			SUMS=SUMS+EVEC(I)*S(I,J)*EVEC(J)

	mov	ebp, esi
	jg	SHORT $L101298
	cmp	ebp, 1
	jge	SHORT $L101299
$L101298:
	call	___FHranger@0
$L101299:
	imul	ebp, DWORD PTR _7S$[esp+1428]
	cmp	ebx, edi
	jg	SHORT $L101300
	cmp	ebx, 1
	jge	SHORT $L101301
$L101300:
	call	___FHranger@0
$L101301:
	add	ebp, ebx
	mov	eax, DWORD PTR _S$[esp+1388]
	add	ebp, DWORD PTR _7S$[esp+1400]
	cmp	esi, edi
	fld	QWORD PTR [eax+ebp*8]
	jg	SHORT $L101302
	cmp	esi, 1
	jge	SHORT $L101303
$L101302:
	call	___FHranger@0
$L101303:
	mov	eax, DWORD PTR _EVEC$[esp+1388]
	cmp	ebx, edi
	fmul	QWORD PTR [eax+esi*8-8]
	jg	SHORT $L101304
	cmp	ebx, 1
	jge	SHORT $L101305
$L101304:
	call	___FHranger@0
$L101305:
	mov	eax, DWORD PTR _EVEC$[esp+1388]
	inc	ebx
	mov	DWORD PTR $T98573[esp+1404], 1472	; 000005c0H
	fmul	QWORD PTR [eax+ebx*8-16]
	fadd	QWORD PTR _SUMS
	fstp	QWORD PTR _SUMS

; 1472 : 		DO I=1,M

	dec	DWORD PTR -600+[esp+1392]
	cmp	DWORD PTR -600+[esp+1392], 0
	jne	$L98520
	mov	DWORD PTR $T98573[esp+1404], 1471	; 000005bfH
$L649:

; 1471 : 	DO J=1,M

	mov	eax, DWORD PTR -592+[esp+1392]
	inc	esi
	mov	DWORD PTR $T98573[esp+1404], 1471	; 000005bfH
	dec	eax
	mov	DWORD PTR -592+[esp+1392], eax
	test	eax, eax
	jne	$L98519
	mov	DWORD PTR $T98573[esp+1404], 1476	; 000005c4H
$L648:

; 1475 : 		END DO
; 1476 : 	END DO

	fld	QWORD PTR _SUMH
	fdiv	QWORD PTR _SUMS
	mov	ecx, DWORD PTR _EMIN$[esp+1388]
	mov	DWORD PTR $T98573[esp+1404], 1478	; 000005c6H

; 1477 : C
; 1478 : 	ENG=SUMH/SUMS

	mov	DWORD PTR $T98573[esp+1404], 1479	; 000005c7H
	mov	DWORD PTR $T98573[esp+1404], 1480	; 000005c8H
	fstp	QWORD PTR _ENG

; 1479 : 	EMIN=ENG

	mov	edx, DWORD PTR _ENG+4
	mov	eax, DWORD PTR _ENG
	mov	DWORD PTR [ecx+4], edx
	mov	DWORD PTR [ecx], eax
$L646:

; 1480 : 	END IF

	xor	esi, esi
	mov	ecx, DWORD PTR _DHDA1
	mov	DWORD PTR $T98573[esp+1404], 1483	; 000005cbH

; 1481 : C
; 1482 : C     -----------------------------------------------------      
; 1483 :       DEALLOCATE(	TEMP,DHDA1,DHDA2,DHDA5,DHDAZ,DHDBZ,DSDA1,DSDA2,

	push	esi
	push	8
	push	OFFSET FLAT:_7DHDA1
	push	ecx
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _DHDA2
	mov	DWORD PTR _DHDA1, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7DHDA2
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _DHDA5
	mov	DWORD PTR _DHDA2, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7DHDA5
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _DHDAZ
	mov	DWORD PTR _DHDA5, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7DHDAZ
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _DHDBZ
	mov	DWORD PTR _DHDAZ, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7DHDBZ
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _DSDA1
	mov	DWORD PTR _DHDBZ, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7DSDA1
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _DSDA2
	mov	DWORD PTR _DSDA1, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7DSDA2
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _DSDA5
	mov	DWORD PTR _DSDA2, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7DSDA5
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _DSDAZ
	mov	DWORD PTR _DSDA5, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7DSDAZ
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _DSDBZ
	mov	DWORD PTR _DSDAZ, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7DSDBZ
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _GRADM
	mov	DWORD PTR _DSDBZ, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7GRADM
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _DS
	mov	DWORD PTR _GRADM, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7DS
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	DWORD PTR _DS, esi
	push	esi
	push	8
	push	OFFSET FLAT:_7DH
	mov	eax, DWORD PTR _DH
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _GRADA
	mov	DWORD PTR _DH, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7GRADA
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _GRADC
	mov	DWORD PTR _GRADA, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7GRADC
	push	eax
	call	__Deallocate_i
	add	esp, 16					; 00000010H
	mov	eax, DWORD PTR _TEMP
	mov	DWORD PTR _GRADC, esi
	npad	1
	push	esi
	push	8
	push	OFFSET FLAT:_7TEMP
	push	eax
	call	__Deallocate_i
	mov	DWORD PTR $T98573[esp+1420], 1488	; 000005d0H
	add	esp, 16					; 00000010H
	mov	DWORD PTR _TEMP, esi
$L652:

; 1484 : 	-			DSDA5,DSDAZ,DSDBZ,GRADM,DS,DH,GRADA,GRADC)
; 1485 : 	END IF
; 1486 : 
; 1487 :       RETURN
; 1488 :       END

	mov	DWORD PTR $T98573[esp+1404], 1488	; 000005d0H
	call	___FHunlink@0
	pop	ebp
	pop	edi
	pop	esi
	pop	ebx
	add	esp, 1376				; 00000560H
	ret	80					; 00000050H
	npad	7
	npad	7
_HAMILTON@80 ENDP
_TEXT	ENDS
END
