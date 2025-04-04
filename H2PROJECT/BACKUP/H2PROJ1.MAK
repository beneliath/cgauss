# Microsoft Developer Studio Generated NMAKE File, Format Version 4.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

!IF "$(CFG)" == ""
CFG=H2PROJ - Win32 Debug
!MESSAGE No configuration specified.  Defaulting to H2PROJ - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "H2PROJ - Win32 Release" && "$(CFG)" != "H2PROJ - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE on this makefile
!MESSAGE by defining the macro CFG on the command line.  For example:
!MESSAGE 
!MESSAGE NMAKE /f "H2PROJ1.MAK" CFG="H2PROJ - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "H2PROJ - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "H2PROJ - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 
################################################################################
# Begin Project
# PROP Target_Last_Scanned "H2PROJ - Win32 Debug"
F90=fl32.exe
RSC=rc.exe

!IF  "$(CFG)" == "H2PROJ - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
OUTDIR=.\Release
INTDIR=.\Release

ALL : "$(OUTDIR)\H2PROJ1.exe" "$(OUTDIR)\H2PROJ1.bsc"

CLEAN : 
	-@erase ".\Release\H2PROJ1.bsc"
	-@erase ".\Release\trudge.sbr"
	-@erase ".\Release\OVERLP.SBR"
	-@erase ".\Release\ELREP.SBR"
	-@erase ".\Release\FU.SBR"
	-@erase ".\Release\DLINMIN.sbr"
	-@erase ".\Release\KFACT.SBR"
	-@erase ".\Release\eign.sbr"
	-@erase ".\Release\VEC.SBR"
	-@erase ".\Release\DF1DIM.sbr"
	-@erase ".\Release\GEOM.SBR"
	-@erase ".\Release\GUESS.SBR"
	-@erase ".\Release\YACOBI.SBR"
	-@erase ".\Release\SMMULT.SBR"
	-@erase ".\Release\H2.SBR"
	-@erase ".\Release\F0.SBR"
	-@erase ".\Release\GRNUE.SBR"
	-@erase ".\Release\BLANKM.SBR"
	-@erase ".\Release\DBRENT.sbr"
	-@erase ".\Release\search.sbr"
	-@erase ".\Release\TK.SBR"
	-@erase ".\Release\TKJOHN.SBR"
	-@erase ".\Release\VECH.SBR"
	-@erase ".\Release\HAMILTON.SBR"
	-@erase ".\Release\TAB.SBR"
	-@erase ".\Release\ENERGY.SBR"
	-@erase ".\Release\MAIN.SBR"
	-@erase ".\Release\FUNC.SBR"
	-@erase ".\Release\GRELE.SBR"
	-@erase ".\Release\BIOUT.SBR"
	-@erase ".\Release\F.SBR"
	-@erase ".\Release\GRAD.SBR"
	-@erase ".\Release\BLANKDIR.SBR"
	-@erase ".\Release\ENUCL.SBR"
	-@erase ".\Release\DEVEC.SBR"
	-@erase ".\Release\DK.SBR"
	-@erase ".\Release\DMNBRAK.SBR"
	-@erase ".\Release\GRKIN.SBR"
	-@erase ".\Release\MATRIXM.SBR"
	-@erase ".\Release\GEMPRT.SBR"
	-@erase ".\Release\DPOWELL.sbr"
	-@erase ".\Release\MATRIXCOPY.SBR"
	-@erase ".\Release\GROVER.SBR"
	-@erase ".\Release\TRANSPOSE.SBR"
	-@erase ".\Release\H2PROJ1.exe"
	-@erase ".\Release\DEVEC.OBJ"
	-@erase ".\Release\DK.OBJ"
	-@erase ".\Release\DMNBRAK.OBJ"
	-@erase ".\Release\GRKIN.OBJ"
	-@erase ".\Release\MATRIXM.OBJ"
	-@erase ".\Release\GEMPRT.OBJ"
	-@erase ".\Release\DPOWELL.obj"
	-@erase ".\Release\MATRIXCOPY.OBJ"
	-@erase ".\Release\GROVER.OBJ"
	-@erase ".\Release\TRANSPOSE.OBJ"
	-@erase ".\Release\trudge.obj"
	-@erase ".\Release\OVERLP.OBJ"
	-@erase ".\Release\ELREP.OBJ"
	-@erase ".\Release\FU.OBJ"
	-@erase ".\Release\DLINMIN.obj"
	-@erase ".\Release\KFACT.OBJ"
	-@erase ".\Release\eign.obj"
	-@erase ".\Release\VEC.OBJ"
	-@erase ".\Release\DF1DIM.obj"
	-@erase ".\Release\GEOM.OBJ"
	-@erase ".\Release\GUESS.OBJ"
	-@erase ".\Release\YACOBI.OBJ"
	-@erase ".\Release\SMMULT.OBJ"
	-@erase ".\Release\H2.OBJ"
	-@erase ".\Release\F0.OBJ"
	-@erase ".\Release\GRNUE.OBJ"
	-@erase ".\Release\BLANKM.OBJ"
	-@erase ".\Release\DBRENT.obj"
	-@erase ".\Release\search.obj"
	-@erase ".\Release\TK.OBJ"
	-@erase ".\Release\TKJOHN.OBJ"
	-@erase ".\Release\VECH.OBJ"
	-@erase ".\Release\HAMILTON.OBJ"
	-@erase ".\Release\TAB.OBJ"
	-@erase ".\Release\ENERGY.OBJ"
	-@erase ".\Release\MAIN.OBJ"
	-@erase ".\Release\FUNC.OBJ"
	-@erase ".\Release\GRELE.OBJ"
	-@erase ".\Release\BIOUT.OBJ"
	-@erase ".\Release\F.OBJ"
	-@erase ".\Release\GRAD.OBJ"
	-@erase ".\Release\BLANKDIR.OBJ"
	-@erase ".\Release\ENUCL.OBJ"
	-@erase ".\Release\H2PROJ1.map"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE F90 /Ox /I "Release/" /c /nologo
# ADD F90 /Ox /Op /G5 /4Yb /W0 /FR /Zi /I "Release/" /c /nologo /MT
F90_PROJ=/Ox /Op /G5 /4Yb /W0 /FR"Release/" /Zi /I "Release/" /c /nologo /MT\
 /Fo"Release/" /Fd"Release/H2PROJ1.pdb" 
F90_OBJS=.\Release/
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/H2PROJ1.bsc" 
BSC32_SBRS= \
	"$(INTDIR)/trudge.sbr" \
	"$(INTDIR)/OVERLP.SBR" \
	"$(INTDIR)/ELREP.SBR" \
	"$(INTDIR)/FU.SBR" \
	"$(INTDIR)/DLINMIN.sbr" \
	"$(INTDIR)/KFACT.SBR" \
	"$(INTDIR)/eign.sbr" \
	"$(INTDIR)/VEC.SBR" \
	"$(INTDIR)/DF1DIM.sbr" \
	"$(INTDIR)/GEOM.SBR" \
	"$(INTDIR)/GUESS.SBR" \
	"$(INTDIR)/YACOBI.SBR" \
	"$(INTDIR)/SMMULT.SBR" \
	"$(INTDIR)/H2.SBR" \
	"$(INTDIR)/F0.SBR" \
	"$(INTDIR)/GRNUE.SBR" \
	"$(INTDIR)/BLANKM.SBR" \
	"$(INTDIR)/DBRENT.sbr" \
	"$(INTDIR)/search.sbr" \
	"$(INTDIR)/TK.SBR" \
	"$(INTDIR)/TKJOHN.SBR" \
	"$(INTDIR)/VECH.SBR" \
	"$(INTDIR)/HAMILTON.SBR" \
	"$(INTDIR)/TAB.SBR" \
	"$(INTDIR)/ENERGY.SBR" \
	"$(INTDIR)/MAIN.SBR" \
	"$(INTDIR)/FUNC.SBR" \
	"$(INTDIR)/GRELE.SBR" \
	"$(INTDIR)/BIOUT.SBR" \
	"$(INTDIR)/F.SBR" \
	"$(INTDIR)/GRAD.SBR" \
	"$(INTDIR)/BLANKDIR.SBR" \
	"$(INTDIR)/ENUCL.SBR" \
	"$(INTDIR)/DEVEC.SBR" \
	"$(INTDIR)/DK.SBR" \
	"$(INTDIR)/DMNBRAK.SBR" \
	"$(INTDIR)/GRKIN.SBR" \
	"$(INTDIR)/MATRIXM.SBR" \
	"$(INTDIR)/GEMPRT.SBR" \
	"$(INTDIR)/DPOWELL.sbr" \
	"$(INTDIR)/MATRIXCOPY.SBR" \
	"$(INTDIR)/GROVER.SBR" \
	"$(INTDIR)/TRANSPOSE.SBR"

"$(OUTDIR)\H2PROJ1.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /profile /map /debug /debugtype:both /machine:I386
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /profile\
 /map:"$(INTDIR)/H2PROJ1.map" /debug /debugtype:both /machine:I386\
 /out:"$(OUTDIR)/H2PROJ1.exe" 
LINK32_OBJS= \
	"$(INTDIR)/DEVEC.OBJ" \
	"$(INTDIR)/DK.OBJ" \
	"$(INTDIR)/DMNBRAK.OBJ" \
	"$(INTDIR)/GRKIN.OBJ" \
	"$(INTDIR)/MATRIXM.OBJ" \
	"$(INTDIR)/GEMPRT.OBJ" \
	"$(INTDIR)/DPOWELL.obj" \
	"$(INTDIR)/MATRIXCOPY.OBJ" \
	"$(INTDIR)/GROVER.OBJ" \
	"$(INTDIR)/TRANSPOSE.OBJ" \
	"$(INTDIR)/trudge.obj" \
	"$(INTDIR)/OVERLP.OBJ" \
	"$(INTDIR)/ELREP.OBJ" \
	"$(INTDIR)/FU.OBJ" \
	"$(INTDIR)/DLINMIN.obj" \
	"$(INTDIR)/KFACT.OBJ" \
	"$(INTDIR)/eign.obj" \
	"$(INTDIR)/VEC.OBJ" \
	"$(INTDIR)/DF1DIM.obj" \
	"$(INTDIR)/GEOM.OBJ" \
	"$(INTDIR)/GUESS.OBJ" \
	"$(INTDIR)/YACOBI.OBJ" \
	"$(INTDIR)/SMMULT.OBJ" \
	"$(INTDIR)/H2.OBJ" \
	"$(INTDIR)/F0.OBJ" \
	"$(INTDIR)/GRNUE.OBJ" \
	"$(INTDIR)/BLANKM.OBJ" \
	"$(INTDIR)/DBRENT.obj" \
	"$(INTDIR)/search.obj" \
	"$(INTDIR)/TK.OBJ" \
	"$(INTDIR)/TKJOHN.OBJ" \
	"$(INTDIR)/VECH.OBJ" \
	"$(INTDIR)/HAMILTON.OBJ" \
	"$(INTDIR)/TAB.OBJ" \
	"$(INTDIR)/ENERGY.OBJ" \
	"$(INTDIR)/MAIN.OBJ" \
	"$(INTDIR)/FUNC.OBJ" \
	"$(INTDIR)/GRELE.OBJ" \
	"$(INTDIR)/BIOUT.OBJ" \
	"$(INTDIR)/F.OBJ" \
	"$(INTDIR)/GRAD.OBJ" \
	"$(INTDIR)/BLANKDIR.OBJ" \
	"$(INTDIR)/ENUCL.OBJ"

"$(OUTDIR)\H2PROJ1.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "H2PROJ - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
OUTDIR=.\Debug
INTDIR=.\Debug

ALL : "$(OUTDIR)\H2PROJ1.exe" "$(OUTDIR)\H2PROJ1.bsc"

CLEAN : 
	-@erase ".\Debug\H2PROJ1.bsc"
	-@erase ".\Debug\YACOBI.SBR"
	-@erase ".\Debug\GRNUE.SBR"
	-@erase ".\Debug\BLANKM.SBR"
	-@erase ".\Debug\TRANSPOSE.SBR"
	-@erase ".\Debug\DBRENT.sbr"
	-@erase ".\Debug\search.sbr"
	-@erase ".\Debug\ENUCL.SBR"
	-@erase ".\Debug\DK.SBR"
	-@erase ".\Debug\TKJOHN.SBR"
	-@erase ".\Debug\GRAD.SBR"
	-@erase ".\Debug\TK.SBR"
	-@erase ".\Debug\GEMPRT.SBR"
	-@erase ".\Debug\MATRIXM.SBR"
	-@erase ".\Debug\DPOWELL.sbr"
	-@erase ".\Debug\F.SBR"
	-@erase ".\Debug\FU.SBR"
	-@erase ".\Debug\OVERLP.SBR"
	-@erase ".\Debug\BIOUT.SBR"
	-@erase ".\Debug\DLINMIN.sbr"
	-@erase ".\Debug\HAMILTON.SBR"
	-@erase ".\Debug\eign.sbr"
	-@erase ".\Debug\H2.SBR"
	-@erase ".\Debug\F0.SBR"
	-@erase ".\Debug\GEOM.SBR"
	-@erase ".\Debug\VEC.SBR"
	-@erase ".\Debug\SMMULT.SBR"
	-@erase ".\Debug\DEVEC.SBR"
	-@erase ".\Debug\MATRIXCOPY.SBR"
	-@erase ".\Debug\GRKIN.SBR"
	-@erase ".\Debug\GROVER.SBR"
	-@erase ".\Debug\DMNBRAK.SBR"
	-@erase ".\Debug\BLANKDIR.SBR"
	-@erase ".\Debug\VECH.SBR"
	-@erase ".\Debug\ENERGY.SBR"
	-@erase ".\Debug\ELREP.SBR"
	-@erase ".\Debug\trudge.sbr"
	-@erase ".\Debug\KFACT.SBR"
	-@erase ".\Debug\GRELE.SBR"
	-@erase ".\Debug\TAB.SBR"
	-@erase ".\Debug\MAIN.SBR"
	-@erase ".\Debug\FUNC.SBR"
	-@erase ".\Debug\DF1DIM.sbr"
	-@erase ".\Debug\GUESS.SBR"
	-@erase ".\Debug\H2PROJ1.exe"
	-@erase ".\Debug\DMNBRAK.OBJ"
	-@erase ".\Debug\BLANKDIR.OBJ"
	-@erase ".\Debug\VECH.OBJ"
	-@erase ".\Debug\ENERGY.OBJ"
	-@erase ".\Debug\ELREP.OBJ"
	-@erase ".\Debug\trudge.obj"
	-@erase ".\Debug\KFACT.OBJ"
	-@erase ".\Debug\GRELE.OBJ"
	-@erase ".\Debug\TAB.OBJ"
	-@erase ".\Debug\MAIN.OBJ"
	-@erase ".\Debug\FUNC.OBJ"
	-@erase ".\Debug\DF1DIM.obj"
	-@erase ".\Debug\GUESS.OBJ"
	-@erase ".\Debug\YACOBI.OBJ"
	-@erase ".\Debug\GRNUE.OBJ"
	-@erase ".\Debug\BLANKM.OBJ"
	-@erase ".\Debug\TRANSPOSE.OBJ"
	-@erase ".\Debug\DBRENT.obj"
	-@erase ".\Debug\search.obj"
	-@erase ".\Debug\ENUCL.OBJ"
	-@erase ".\Debug\DK.OBJ"
	-@erase ".\Debug\TKJOHN.OBJ"
	-@erase ".\Debug\GRAD.OBJ"
	-@erase ".\Debug\TK.OBJ"
	-@erase ".\Debug\GEMPRT.OBJ"
	-@erase ".\Debug\MATRIXM.OBJ"
	-@erase ".\Debug\DPOWELL.obj"
	-@erase ".\Debug\F.OBJ"
	-@erase ".\Debug\FU.OBJ"
	-@erase ".\Debug\OVERLP.OBJ"
	-@erase ".\Debug\BIOUT.OBJ"
	-@erase ".\Debug\DLINMIN.obj"
	-@erase ".\Debug\HAMILTON.OBJ"
	-@erase ".\Debug\eign.obj"
	-@erase ".\Debug\H2.OBJ"
	-@erase ".\Debug\F0.OBJ"
	-@erase ".\Debug\GEOM.OBJ"
	-@erase ".\Debug\VEC.OBJ"
	-@erase ".\Debug\SMMULT.OBJ"
	-@erase ".\Debug\DEVEC.OBJ"
	-@erase ".\Debug\MATRIXCOPY.OBJ"
	-@erase ".\Debug\GRKIN.OBJ"
	-@erase ".\Debug\GROVER.OBJ"
	-@erase ".\Debug\H2PROJ1.map"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE F90 /Zi /I "Debug/" /c /nologo
# ADD F90 /Ox /Op /G5 /4Yb /W0 /FR /Zi /I "Debug/" /c /nologo /MT
F90_PROJ=/Ox /Op /G5 /4Yb /W0 /FR"Debug/" /Zi /I "Debug/" /c /nologo /MT\
 /Fo"Debug/" /Fd"Debug/H2PROJ1.pdb" 
F90_OBJS=.\Debug/
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/H2PROJ1.bsc" 
BSC32_SBRS= \
	"$(INTDIR)/YACOBI.SBR" \
	"$(INTDIR)/GRNUE.SBR" \
	"$(INTDIR)/BLANKM.SBR" \
	"$(INTDIR)/TRANSPOSE.SBR" \
	"$(INTDIR)/DBRENT.sbr" \
	"$(INTDIR)/search.sbr" \
	"$(INTDIR)/ENUCL.SBR" \
	"$(INTDIR)/DK.SBR" \
	"$(INTDIR)/TKJOHN.SBR" \
	"$(INTDIR)/GRAD.SBR" \
	"$(INTDIR)/TK.SBR" \
	"$(INTDIR)/GEMPRT.SBR" \
	"$(INTDIR)/MATRIXM.SBR" \
	"$(INTDIR)/DPOWELL.sbr" \
	"$(INTDIR)/F.SBR" \
	"$(INTDIR)/FU.SBR" \
	"$(INTDIR)/OVERLP.SBR" \
	"$(INTDIR)/BIOUT.SBR" \
	"$(INTDIR)/DLINMIN.sbr" \
	"$(INTDIR)/HAMILTON.SBR" \
	"$(INTDIR)/eign.sbr" \
	"$(INTDIR)/H2.SBR" \
	"$(INTDIR)/F0.SBR" \
	"$(INTDIR)/GEOM.SBR" \
	"$(INTDIR)/VEC.SBR" \
	"$(INTDIR)/SMMULT.SBR" \
	"$(INTDIR)/DEVEC.SBR" \
	"$(INTDIR)/MATRIXCOPY.SBR" \
	"$(INTDIR)/GRKIN.SBR" \
	"$(INTDIR)/GROVER.SBR" \
	"$(INTDIR)/DMNBRAK.SBR" \
	"$(INTDIR)/BLANKDIR.SBR" \
	"$(INTDIR)/VECH.SBR" \
	"$(INTDIR)/ENERGY.SBR" \
	"$(INTDIR)/ELREP.SBR" \
	"$(INTDIR)/trudge.sbr" \
	"$(INTDIR)/KFACT.SBR" \
	"$(INTDIR)/GRELE.SBR" \
	"$(INTDIR)/TAB.SBR" \
	"$(INTDIR)/MAIN.SBR" \
	"$(INTDIR)/FUNC.SBR" \
	"$(INTDIR)/DF1DIM.sbr" \
	"$(INTDIR)/GUESS.SBR"

"$(OUTDIR)\H2PROJ1.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /profile /map /debug /debugtype:both /machine:I386
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /profile\
 /map:"$(INTDIR)/H2PROJ1.map" /debug /debugtype:both /machine:I386\
 /out:"$(OUTDIR)/H2PROJ1.exe" 
LINK32_OBJS= \
	"$(INTDIR)/DMNBRAK.OBJ" \
	"$(INTDIR)/BLANKDIR.OBJ" \
	"$(INTDIR)/VECH.OBJ" \
	"$(INTDIR)/ENERGY.OBJ" \
	"$(INTDIR)/ELREP.OBJ" \
	"$(INTDIR)/trudge.obj" \
	"$(INTDIR)/KFACT.OBJ" \
	"$(INTDIR)/GRELE.OBJ" \
	"$(INTDIR)/TAB.OBJ" \
	"$(INTDIR)/MAIN.OBJ" \
	"$(INTDIR)/FUNC.OBJ" \
	"$(INTDIR)/DF1DIM.obj" \
	"$(INTDIR)/GUESS.OBJ" \
	"$(INTDIR)/YACOBI.OBJ" \
	"$(INTDIR)/GRNUE.OBJ" \
	"$(INTDIR)/BLANKM.OBJ" \
	"$(INTDIR)/TRANSPOSE.OBJ" \
	"$(INTDIR)/DBRENT.obj" \
	"$(INTDIR)/search.obj" \
	"$(INTDIR)/ENUCL.OBJ" \
	"$(INTDIR)/DK.OBJ" \
	"$(INTDIR)/TKJOHN.OBJ" \
	"$(INTDIR)/GRAD.OBJ" \
	"$(INTDIR)/TK.OBJ" \
	"$(INTDIR)/GEMPRT.OBJ" \
	"$(INTDIR)/MATRIXM.OBJ" \
	"$(INTDIR)/DPOWELL.obj" \
	"$(INTDIR)/F.OBJ" \
	"$(INTDIR)/FU.OBJ" \
	"$(INTDIR)/OVERLP.OBJ" \
	"$(INTDIR)/BIOUT.OBJ" \
	"$(INTDIR)/DLINMIN.obj" \
	"$(INTDIR)/HAMILTON.OBJ" \
	"$(INTDIR)/eign.obj" \
	"$(INTDIR)/H2.OBJ" \
	"$(INTDIR)/F0.OBJ" \
	"$(INTDIR)/GEOM.OBJ" \
	"$(INTDIR)/VEC.OBJ" \
	"$(INTDIR)/SMMULT.OBJ" \
	"$(INTDIR)/DEVEC.OBJ" \
	"$(INTDIR)/MATRIXCOPY.OBJ" \
	"$(INTDIR)/GRKIN.OBJ" \
	"$(INTDIR)/GROVER.OBJ"

"$(OUTDIR)\H2PROJ1.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ENDIF 

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

################################################################################
# Begin Target

# Name "H2PROJ - Win32 Release"
# Name "H2PROJ - Win32 Debug"

!IF  "$(CFG)" == "H2PROJ - Win32 Release"

!ELSEIF  "$(CFG)" == "H2PROJ - Win32 Debug"

!ENDIF 

################################################################################
# Begin Source File

SOURCE=.\BIOUT.FOR

"$(INTDIR)\BIOUT.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\BIOUT.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\BLANKDIR.FOR

"$(INTDIR)\BLANKDIR.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\BLANKDIR.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\BLANKM.FOR

"$(INTDIR)\BLANKM.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\BLANKM.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\DEVEC.FOR

"$(INTDIR)\DEVEC.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\DEVEC.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\DK.FOR

"$(INTDIR)\DK.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\DK.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\eign.for

"$(INTDIR)\eign.obj" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\eign.sbr" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\ELREP.FOR

"$(INTDIR)\ELREP.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\ELREP.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\ENERGY.FOR

"$(INTDIR)\ENERGY.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\ENERGY.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\ENUCL.FOR

"$(INTDIR)\ENUCL.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\ENUCL.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\F.FOR

"$(INTDIR)\F.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\F.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\F0.FOR

"$(INTDIR)\F0.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\F0.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\FUNC.FOR

"$(INTDIR)\FUNC.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\FUNC.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\GEMPRT.FOR

"$(INTDIR)\GEMPRT.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\GEMPRT.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\GEOM.FOR

"$(INTDIR)\GEOM.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\GEOM.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\GRAD.FOR

"$(INTDIR)\GRAD.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\GRAD.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\GRELE.FOR

"$(INTDIR)\GRELE.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\GRELE.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\GRKIN.FOR

"$(INTDIR)\GRKIN.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\GRKIN.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\GRNUE.FOR

"$(INTDIR)\GRNUE.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\GRNUE.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\GROVER.FOR

"$(INTDIR)\GROVER.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\GROVER.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\GUESS.FOR

"$(INTDIR)\GUESS.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\GUESS.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\H2.FOR

"$(INTDIR)\H2.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\H2.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\HAMILTON.FOR

"$(INTDIR)\HAMILTON.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\HAMILTON.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\KFACT.FOR

"$(INTDIR)\KFACT.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\KFACT.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\MAIN.FOR

"$(INTDIR)\MAIN.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\MAIN.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\MATRIXCOPY.FOR

"$(INTDIR)\MATRIXCOPY.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\MATRIXCOPY.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\MATRIXM.FOR

"$(INTDIR)\MATRIXM.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\MATRIXM.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\OVERLP.FOR

"$(INTDIR)\OVERLP.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\OVERLP.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\search.for

"$(INTDIR)\search.obj" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\search.sbr" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\SMMULT.FOR

"$(INTDIR)\SMMULT.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\SMMULT.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\TAB.FOR

"$(INTDIR)\TAB.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\TAB.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\TK.FOR

"$(INTDIR)\TK.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\TK.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\TKJOHN.FOR

"$(INTDIR)\TKJOHN.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\TKJOHN.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\TRANSPOSE.FOR

"$(INTDIR)\TRANSPOSE.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\TRANSPOSE.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\trudge.for

"$(INTDIR)\trudge.obj" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\trudge.sbr" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\VEC.FOR

"$(INTDIR)\VEC.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\VEC.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\VECH.FOR

"$(INTDIR)\VECH.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\VECH.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\YACOBI.FOR

"$(INTDIR)\YACOBI.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\YACOBI.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=D:\APPS\WIN\FPSNT\RECIPES\double\DBRENT.for

BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\DBRENT.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\DBRENT.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=D:\APPS\WIN\FPSNT\RECIPES\double\DF1DIM.for

BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\DF1DIM.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\DF1DIM.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=D:\APPS\WIN\FPSNT\RECIPES\double\DLINMIN.for

BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\DLINMIN.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\DLINMIN.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=D:\APPS\WIN\FPSNT\RECIPES\double\DPOWELL.for

BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\DPOWELL.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\DPOWELL.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=D:\APPS\WIN\FPSNT\RECIPES\double\DMNBRAK.FOR

BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\DMNBRAK.OBJ" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\DMNBRAK.SBR" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

# End Source File
################################################################################
# Begin Source File

SOURCE=.\FU.FOR

"$(INTDIR)\FU.OBJ" : $(SOURCE) "$(INTDIR)"

"$(INTDIR)\FU.SBR" : $(SOURCE) "$(INTDIR)"


# End Source File
# End Target
# End Project
################################################################################
