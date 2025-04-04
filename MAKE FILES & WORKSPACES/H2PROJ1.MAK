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
	-@erase ".\Release\search.sbr"
	-@erase ".\Release\GEOM.SBR"
	-@erase ".\Release\TKJOHN.SBR"
	-@erase ".\Release\DMNBRAK.SBR"
	-@erase ".\Release\GRKIN.SBR"
	-@erase ".\Release\MATRIXM.SBR"
	-@erase ".\Release\TRANSPOSE.SBR"
	-@erase ".\Release\GEMPRT.SBR"
	-@erase ".\Release\DK.SBR"
	-@erase ".\Release\BLANKDIR.SBR"
	-@erase ".\Release\KFACT.SBR"
	-@erase ".\Release\DLINMIN.sbr"
	-@erase ".\Release\VECH.SBR"
	-@erase ".\Release\OVERLP.SBR"
	-@erase ".\Release\GUESS.SBR"
	-@erase ".\Release\FUNC.SBR"
	-@erase ".\Release\DF1DIM.sbr"
	-@erase ".\Release\eign.sbr"
	-@erase ".\Release\GRNUE.SBR"
	-@erase ".\Release\FU.SBR"
	-@erase ".\Release\SMMULT.SBR"
	-@erase ".\Release\DBRENT.sbr"
	-@erase ".\Release\DPOWELL.sbr"
	-@erase ".\Release\H2.SBR"
	-@erase ".\Release\F0.SBR"
	-@erase ".\Release\GRAD.SBR"
	-@erase ".\Release\TAB.SBR"
	-@erase ".\Release\GROVER.SBR"
	-@erase ".\Release\ELREP.SBR"
	-@erase ".\Release\TK.SBR"
	-@erase ".\Release\ENERGY.SBR"
	-@erase ".\Release\GRELE.SBR"
	-@erase ".\Release\VEC.SBR"
	-@erase ".\Release\trudge.sbr"
	-@erase ".\Release\BIOUT.SBR"
	-@erase ".\Release\MAIN.SBR"
	-@erase ".\Release\F.SBR"
	-@erase ".\Release\HAMILTON.SBR"
	-@erase ".\Release\YACOBI.SBR"
	-@erase ".\Release\BLANKM.SBR"
	-@erase ".\Release\MATRIXCOPY.SBR"
	-@erase ".\Release\ENUCL.SBR"
	-@erase ".\Release\DEVEC.SBR"
	-@erase ".\Release\H2PROJ1.exe"
	-@erase ".\Release\BIOUT.OBJ"
	-@erase ".\Release\MAIN.OBJ"
	-@erase ".\Release\F.OBJ"
	-@erase ".\Release\HAMILTON.OBJ"
	-@erase ".\Release\YACOBI.OBJ"
	-@erase ".\Release\BLANKM.OBJ"
	-@erase ".\Release\MATRIXCOPY.OBJ"
	-@erase ".\Release\ENUCL.OBJ"
	-@erase ".\Release\DEVEC.OBJ"
	-@erase ".\Release\search.obj"
	-@erase ".\Release\GEOM.OBJ"
	-@erase ".\Release\TKJOHN.OBJ"
	-@erase ".\Release\DMNBRAK.OBJ"
	-@erase ".\Release\GRKIN.OBJ"
	-@erase ".\Release\MATRIXM.OBJ"
	-@erase ".\Release\TRANSPOSE.OBJ"
	-@erase ".\Release\GEMPRT.OBJ"
	-@erase ".\Release\DK.OBJ"
	-@erase ".\Release\BLANKDIR.OBJ"
	-@erase ".\Release\KFACT.OBJ"
	-@erase ".\Release\DLINMIN.obj"
	-@erase ".\Release\VECH.OBJ"
	-@erase ".\Release\OVERLP.OBJ"
	-@erase ".\Release\GUESS.OBJ"
	-@erase ".\Release\FUNC.OBJ"
	-@erase ".\Release\DF1DIM.obj"
	-@erase ".\Release\eign.obj"
	-@erase ".\Release\GRNUE.OBJ"
	-@erase ".\Release\FU.OBJ"
	-@erase ".\Release\SMMULT.OBJ"
	-@erase ".\Release\DBRENT.obj"
	-@erase ".\Release\DPOWELL.obj"
	-@erase ".\Release\H2.OBJ"
	-@erase ".\Release\F0.OBJ"
	-@erase ".\Release\GRAD.OBJ"
	-@erase ".\Release\TAB.OBJ"
	-@erase ".\Release\GROVER.OBJ"
	-@erase ".\Release\ELREP.OBJ"
	-@erase ".\Release\TK.OBJ"
	-@erase ".\Release\ENERGY.OBJ"
	-@erase ".\Release\GRELE.OBJ"
	-@erase ".\Release\VEC.OBJ"
	-@erase ".\Release\trudge.obj"
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
	".\Release\search.sbr" \
	".\Release\GEOM.SBR" \
	".\Release\TKJOHN.SBR" \
	".\Release\DMNBRAK.SBR" \
	".\Release\GRKIN.SBR" \
	".\Release\MATRIXM.SBR" \
	".\Release\TRANSPOSE.SBR" \
	".\Release\GEMPRT.SBR" \
	".\Release\DK.SBR" \
	".\Release\BLANKDIR.SBR" \
	".\Release\KFACT.SBR" \
	".\Release\DLINMIN.sbr" \
	".\Release\VECH.SBR" \
	".\Release\OVERLP.SBR" \
	".\Release\GUESS.SBR" \
	".\Release\FUNC.SBR" \
	".\Release\DF1DIM.sbr" \
	".\Release\eign.sbr" \
	".\Release\GRNUE.SBR" \
	".\Release\FU.SBR" \
	".\Release\SMMULT.SBR" \
	".\Release\DBRENT.sbr" \
	".\Release\DPOWELL.sbr" \
	".\Release\H2.SBR" \
	".\Release\F0.SBR" \
	".\Release\GRAD.SBR" \
	".\Release\TAB.SBR" \
	".\Release\GROVER.SBR" \
	".\Release\ELREP.SBR" \
	".\Release\TK.SBR" \
	".\Release\ENERGY.SBR" \
	".\Release\GRELE.SBR" \
	".\Release\VEC.SBR" \
	".\Release\trudge.sbr" \
	".\Release\BIOUT.SBR" \
	".\Release\MAIN.SBR" \
	".\Release\F.SBR" \
	".\Release\HAMILTON.SBR" \
	".\Release\YACOBI.SBR" \
	".\Release\BLANKM.SBR" \
	".\Release\MATRIXCOPY.SBR" \
	".\Release\ENUCL.SBR" \
	".\Release\DEVEC.SBR"

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
	".\Release\BIOUT.OBJ" \
	".\Release\MAIN.OBJ" \
	".\Release\F.OBJ" \
	".\Release\HAMILTON.OBJ" \
	".\Release\YACOBI.OBJ" \
	".\Release\BLANKM.OBJ" \
	".\Release\MATRIXCOPY.OBJ" \
	".\Release\ENUCL.OBJ" \
	".\Release\DEVEC.OBJ" \
	".\Release\search.obj" \
	".\Release\GEOM.OBJ" \
	".\Release\TKJOHN.OBJ" \
	".\Release\DMNBRAK.OBJ" \
	".\Release\GRKIN.OBJ" \
	".\Release\MATRIXM.OBJ" \
	".\Release\TRANSPOSE.OBJ" \
	".\Release\GEMPRT.OBJ" \
	".\Release\DK.OBJ" \
	".\Release\BLANKDIR.OBJ" \
	".\Release\KFACT.OBJ" \
	".\Release\DLINMIN.obj" \
	".\Release\VECH.OBJ" \
	".\Release\OVERLP.OBJ" \
	".\Release\GUESS.OBJ" \
	".\Release\FUNC.OBJ" \
	".\Release\DF1DIM.obj" \
	".\Release\eign.obj" \
	".\Release\GRNUE.OBJ" \
	".\Release\FU.OBJ" \
	".\Release\SMMULT.OBJ" \
	".\Release\DBRENT.obj" \
	".\Release\DPOWELL.obj" \
	".\Release\H2.OBJ" \
	".\Release\F0.OBJ" \
	".\Release\GRAD.OBJ" \
	".\Release\TAB.OBJ" \
	".\Release\GROVER.OBJ" \
	".\Release\ELREP.OBJ" \
	".\Release\TK.OBJ" \
	".\Release\ENERGY.OBJ" \
	".\Release\GRELE.OBJ" \
	".\Release\VEC.OBJ" \
	".\Release\trudge.obj"

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
	-@erase ".\Debug\search.sbr"
	-@erase ".\Debug\ELREP.SBR"
	-@erase ".\Debug\TAB.SBR"
	-@erase ".\Debug\TKJOHN.SBR"
	-@erase ".\Debug\KFACT.SBR"
	-@erase ".\Debug\GUESS.SBR"
	-@erase ".\Debug\VEC.SBR"
	-@erase ".\Debug\HAMILTON.SBR"
	-@erase ".\Debug\DK.SBR"
	-@erase ".\Debug\VECH.SBR"
	-@erase ".\Debug\GRNUE.SBR"
	-@erase ".\Debug\TRANSPOSE.SBR"
	-@erase ".\Debug\MAIN.SBR"
	-@erase ".\Debug\FUNC.SBR"
	-@erase ".\Debug\ENUCL.SBR"
	-@erase ".\Debug\DMNBRAK.SBR"
	-@erase ".\Debug\FU.SBR"
	-@erase ".\Debug\MATRIXM.SBR"
	-@erase ".\Debug\SMMULT.SBR"
	-@erase ".\Debug\DPOWELL.sbr"
	-@erase ".\Debug\GEOM.SBR"
	-@erase ".\Debug\BLANKDIR.SBR"
	-@erase ".\Debug\GRELE.SBR"
	-@erase ".\Debug\GRAD.SBR"
	-@erase ".\Debug\GEMPRT.SBR"
	-@erase ".\Debug\H2.SBR"
	-@erase ".\Debug\BIOUT.SBR"
	-@erase ".\Debug\F0.SBR"
	-@erase ".\Debug\GROVER.SBR"
	-@erase ".\Debug\DLINMIN.sbr"
	-@erase ".\Debug\MATRIXCOPY.SBR"
	-@erase ".\Debug\TK.SBR"
	-@erase ".\Debug\ENERGY.SBR"
	-@erase ".\Debug\trudge.sbr"
	-@erase ".\Debug\OVERLP.SBR"
	-@erase ".\Debug\DEVEC.SBR"
	-@erase ".\Debug\DF1DIM.sbr"
	-@erase ".\Debug\F.SBR"
	-@erase ".\Debug\GRKIN.SBR"
	-@erase ".\Debug\eign.sbr"
	-@erase ".\Debug\YACOBI.SBR"
	-@erase ".\Debug\BLANKM.SBR"
	-@erase ".\Debug\DBRENT.sbr"
	-@erase ".\Debug\H2PROJ1.exe"
	-@erase ".\Debug\DF1DIM.obj"
	-@erase ".\Debug\F.OBJ"
	-@erase ".\Debug\GRKIN.OBJ"
	-@erase ".\Debug\eign.obj"
	-@erase ".\Debug\YACOBI.OBJ"
	-@erase ".\Debug\BLANKM.OBJ"
	-@erase ".\Debug\DBRENT.obj"
	-@erase ".\Debug\search.obj"
	-@erase ".\Debug\ELREP.OBJ"
	-@erase ".\Debug\TAB.OBJ"
	-@erase ".\Debug\TKJOHN.OBJ"
	-@erase ".\Debug\KFACT.OBJ"
	-@erase ".\Debug\GUESS.OBJ"
	-@erase ".\Debug\VEC.OBJ"
	-@erase ".\Debug\HAMILTON.OBJ"
	-@erase ".\Debug\DK.OBJ"
	-@erase ".\Debug\VECH.OBJ"
	-@erase ".\Debug\GRNUE.OBJ"
	-@erase ".\Debug\TRANSPOSE.OBJ"
	-@erase ".\Debug\MAIN.OBJ"
	-@erase ".\Debug\FUNC.OBJ"
	-@erase ".\Debug\ENUCL.OBJ"
	-@erase ".\Debug\DMNBRAK.OBJ"
	-@erase ".\Debug\FU.OBJ"
	-@erase ".\Debug\MATRIXM.OBJ"
	-@erase ".\Debug\SMMULT.OBJ"
	-@erase ".\Debug\DPOWELL.obj"
	-@erase ".\Debug\GEOM.OBJ"
	-@erase ".\Debug\BLANKDIR.OBJ"
	-@erase ".\Debug\GRELE.OBJ"
	-@erase ".\Debug\GRAD.OBJ"
	-@erase ".\Debug\GEMPRT.OBJ"
	-@erase ".\Debug\H2.OBJ"
	-@erase ".\Debug\BIOUT.OBJ"
	-@erase ".\Debug\F0.OBJ"
	-@erase ".\Debug\GROVER.OBJ"
	-@erase ".\Debug\DLINMIN.obj"
	-@erase ".\Debug\MATRIXCOPY.OBJ"
	-@erase ".\Debug\TK.OBJ"
	-@erase ".\Debug\ENERGY.OBJ"
	-@erase ".\Debug\trudge.obj"
	-@erase ".\Debug\OVERLP.OBJ"
	-@erase ".\Debug\DEVEC.OBJ"
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
	".\Debug\search.sbr" \
	".\Debug\ELREP.SBR" \
	".\Debug\TAB.SBR" \
	".\Debug\TKJOHN.SBR" \
	".\Debug\KFACT.SBR" \
	".\Debug\GUESS.SBR" \
	".\Debug\VEC.SBR" \
	".\Debug\HAMILTON.SBR" \
	".\Debug\DK.SBR" \
	".\Debug\VECH.SBR" \
	".\Debug\GRNUE.SBR" \
	".\Debug\TRANSPOSE.SBR" \
	".\Debug\MAIN.SBR" \
	".\Debug\FUNC.SBR" \
	".\Debug\ENUCL.SBR" \
	".\Debug\DMNBRAK.SBR" \
	".\Debug\FU.SBR" \
	".\Debug\MATRIXM.SBR" \
	".\Debug\SMMULT.SBR" \
	".\Debug\DPOWELL.sbr" \
	".\Debug\GEOM.SBR" \
	".\Debug\BLANKDIR.SBR" \
	".\Debug\GRELE.SBR" \
	".\Debug\GRAD.SBR" \
	".\Debug\GEMPRT.SBR" \
	".\Debug\H2.SBR" \
	".\Debug\BIOUT.SBR" \
	".\Debug\F0.SBR" \
	".\Debug\GROVER.SBR" \
	".\Debug\DLINMIN.sbr" \
	".\Debug\MATRIXCOPY.SBR" \
	".\Debug\TK.SBR" \
	".\Debug\ENERGY.SBR" \
	".\Debug\trudge.sbr" \
	".\Debug\OVERLP.SBR" \
	".\Debug\DEVEC.SBR" \
	".\Debug\DF1DIM.sbr" \
	".\Debug\F.SBR" \
	".\Debug\GRKIN.SBR" \
	".\Debug\eign.sbr" \
	".\Debug\YACOBI.SBR" \
	".\Debug\BLANKM.SBR" \
	".\Debug\DBRENT.sbr"

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
	".\Debug\DF1DIM.obj" \
	".\Debug\F.OBJ" \
	".\Debug\GRKIN.OBJ" \
	".\Debug\eign.obj" \
	".\Debug\YACOBI.OBJ" \
	".\Debug\BLANKM.OBJ" \
	".\Debug\DBRENT.obj" \
	".\Debug\search.obj" \
	".\Debug\ELREP.OBJ" \
	".\Debug\TAB.OBJ" \
	".\Debug\TKJOHN.OBJ" \
	".\Debug\KFACT.OBJ" \
	".\Debug\GUESS.OBJ" \
	".\Debug\VEC.OBJ" \
	".\Debug\HAMILTON.OBJ" \
	".\Debug\DK.OBJ" \
	".\Debug\VECH.OBJ" \
	".\Debug\GRNUE.OBJ" \
	".\Debug\TRANSPOSE.OBJ" \
	".\Debug\MAIN.OBJ" \
	".\Debug\FUNC.OBJ" \
	".\Debug\ENUCL.OBJ" \
	".\Debug\DMNBRAK.OBJ" \
	".\Debug\FU.OBJ" \
	".\Debug\MATRIXM.OBJ" \
	".\Debug\SMMULT.OBJ" \
	".\Debug\DPOWELL.obj" \
	".\Debug\GEOM.OBJ" \
	".\Debug\BLANKDIR.OBJ" \
	".\Debug\GRELE.OBJ" \
	".\Debug\GRAD.OBJ" \
	".\Debug\GEMPRT.OBJ" \
	".\Debug\H2.OBJ" \
	".\Debug\BIOUT.OBJ" \
	".\Debug\F0.OBJ" \
	".\Debug\GROVER.OBJ" \
	".\Debug\DLINMIN.obj" \
	".\Debug\MATRIXCOPY.OBJ" \
	".\Debug\TK.OBJ" \
	".\Debug\ENERGY.OBJ" \
	".\Debug\trudge.obj" \
	".\Debug\OVERLP.OBJ" \
	".\Debug\DEVEC.OBJ"

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
