# Microsoft Developer Studio Generated NMAKE File, Format Version 40001
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

!IF "$(CFG)" == ""
CFG=cgauss - Win32 Debug
!MESSAGE No configuration specified.  Defaulting to cgauss - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "cgauss - Win32 Release" && "$(CFG)" != "cgauss - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE on this makefile
!MESSAGE by defining the macro CFG on the command line.  For example:
!MESSAGE 
!MESSAGE NMAKE /f "cgauss.mak" CFG="cgauss - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "cgauss - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "cgauss - Win32 Debug" (based on "Win32 (x86) Console Application")
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
# PROP Target_Last_Scanned "cgauss - Win32 Debug"
RSC=rc.exe
CPP=cl.exe
F90=fl32.exe

!IF  "$(CFG)" == "cgauss - Win32 Release"

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

ALL : "$(OUTDIR)\cgauss.exe" "$(OUTDIR)\random_numerical_generation.mod"

CLEAN : 
	-@erase ".\Release\random_numerical_generation.mod"
	-@erase ".\Release\cgauss.exe"
	-@erase ".\Release\shared strings.obj"
	-@erase ".\Release\gradient optimization.obj"
	-@erase ".\Release\scalars_to_be_shared.mod"
	-@erase ".\Release\vectors_to_be_shared.mod"
	-@erase ".\Release\matrices_to_be_shared.mod"
	-@erase ".\Release\default program settings.obj"
	-@erase ".\Release\strings_to_be_shared.mod"
	-@erase ".\Release\program_settings_to_be_shared.mod"
	-@erase ".\Release\jacobi matrix diagonalization.obj"
	-@erase ".\Release\matrix print (TAB).obj"
	-@erase ".\Release\guess random starting data.obj"
	-@erase ".\Release\overlap element.obj"
	-@erase ".\Release\electro-nuclear attraction element.obj"
	-@erase ".\Release\shared settings.obj"
	-@erase ".\Release\load restart point.obj"
	-@erase ".\Release\get current date & time.obj"
	-@erase ".\Release\calculate hamiltonian.obj"
	-@erase ".\Release\return value of the function.obj"
	-@erase ".\Release\lower triangle duplication to matrix.obj"
	-@erase ".\Release\cholesky decondition.obj"
	-@erase ".\Release\print wavefunction stats.obj"
	-@erase ".\Release\matrix blank.obj"
	-@erase ".\Release\error function approximation.obj"
	-@erase ".\Release\load default testing parameters.obj"
	-@erase ".\Release\construction of S & H matrices.obj"
	-@erase ".\Release\solve eigenvalue problem.obj"
	-@erase ".\Release\print preliminary matter.obj"
	-@erase ".\Release\screen & general file text output generation.obj"
	-@erase ".\Release\cholesky condition.obj"
	-@erase ".\Release\prep default testing parameters.obj"
	-@erase ".\Release\return value of the gradient.obj"
	-@erase ".\Release\piped input processor.obj"
	-@erase ".\Release\shared matrices.obj"
	-@erase ".\Release\calculate rayleigh quotient.obj"
	-@erase ".\Release\electron repulsion element.obj"
	-@erase ".\Release\define geometry.obj"
	-@erase ".\Release\calculate kinetic energy.obj"
	-@erase ".\Release\shared vectors.obj"
	-@erase ".\Release\solve 2 by 2.obj"
	-@erase ".\Release\delectro-nuclear attraction element.obj"
	-@erase ".\Release\energy diagnosis.obj"
	-@erase ".\Release\construction of gradient.obj"
	-@erase ".\Release\calculate hamiltonian gradient.obj"
	-@erase ".\Release\delectron repulsion element.obj"
	-@erase ".\Release\doverlap element.obj"
	-@erase ".\Release\program core.obj"
	-@erase ".\Release\dkinetic energy element.obj"
	-@erase ".\Release\kinetic energy element.obj"
	-@erase ".\Release\shared scalars.obj"
	-@erase ".\Release\shared random number generation.obj"
	-@erase ".\Release\simple set printout.obj"
	-@erase ".\Release\construct vector-half.obj"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE F90 /Ox /I "Release/" /c /nologo
# ADD F90 /Ox /Op /G5 /4R8 /I "Release/" /c /nologo /Ob2 /MT /4Mf
# SUBTRACT F90 /FR
F90_PROJ=/Ox /Op /G5 /4R8 /I "Release/" /c /nologo /Ob2 /MT /Fo"Release/" /4Mf 
F90_OBJS=.\Release/
CPP_OBJS=.\Release/
CPP_SBRS=
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/cgauss.bsc" 
BSC32_SBRS=
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /pdb:none /machine:I386
# SUBTRACT LINK32 /profile /map /debug
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /pdb:none /machine:I386\
 /out:"$(OUTDIR)/cgauss.exe" 
LINK32_OBJS= \
	"$(INTDIR)/shared strings.obj" \
	"$(INTDIR)/gradient optimization.obj" \
	"$(INTDIR)/default program settings.obj" \
	"$(INTDIR)/jacobi matrix diagonalization.obj" \
	"$(INTDIR)/matrix print (TAB).obj" \
	"$(INTDIR)/guess random starting data.obj" \
	"$(INTDIR)/overlap element.obj" \
	"$(INTDIR)/electro-nuclear attraction element.obj" \
	"$(INTDIR)/shared settings.obj" \
	"$(INTDIR)/load restart point.obj" \
	"$(INTDIR)/get current date & time.obj" \
	"$(INTDIR)/calculate hamiltonian.obj" \
	"$(INTDIR)/return value of the function.obj" \
	"$(INTDIR)/lower triangle duplication to matrix.obj" \
	"$(INTDIR)/cholesky decondition.obj" \
	"$(INTDIR)/print wavefunction stats.obj" \
	"$(INTDIR)/matrix blank.obj" \
	"$(INTDIR)/error function approximation.obj" \
	"$(INTDIR)/load default testing parameters.obj" \
	"$(INTDIR)/construction of S & H matrices.obj" \
	"$(INTDIR)/solve eigenvalue problem.obj" \
	"$(INTDIR)/print preliminary matter.obj" \
	"$(INTDIR)/screen & general file text output generation.obj" \
	"$(INTDIR)/cholesky condition.obj" \
	"$(INTDIR)/prep default testing parameters.obj" \
	"$(INTDIR)/return value of the gradient.obj" \
	"$(INTDIR)/piped input processor.obj" \
	"$(INTDIR)/shared matrices.obj" \
	"$(INTDIR)/calculate rayleigh quotient.obj" \
	"$(INTDIR)/electron repulsion element.obj" \
	"$(INTDIR)/define geometry.obj" \
	"$(INTDIR)/calculate kinetic energy.obj" \
	"$(INTDIR)/shared vectors.obj" \
	"$(INTDIR)/solve 2 by 2.obj" \
	"$(INTDIR)/delectro-nuclear attraction element.obj" \
	"$(INTDIR)/energy diagnosis.obj" \
	"$(INTDIR)/construction of gradient.obj" \
	"$(INTDIR)/calculate hamiltonian gradient.obj" \
	"$(INTDIR)/delectron repulsion element.obj" \
	"$(INTDIR)/doverlap element.obj" \
	"$(INTDIR)/program core.obj" \
	"$(INTDIR)/dkinetic energy element.obj" \
	"$(INTDIR)/kinetic energy element.obj" \
	"$(INTDIR)/shared scalars.obj" \
	"$(INTDIR)/shared random number generation.obj" \
	"$(INTDIR)/simple set printout.obj" \
	"$(INTDIR)/construct vector-half.obj" \
	"..\..\..\..\..\..\MSDEV\LIB\MATHD.LIB"

"$(OUTDIR)\cgauss.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

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

ALL : "$(OUTDIR)\cgauss.exe" "$(OUTDIR)\random_numerical_generation.mod"\
 "$(OUTDIR)\cgauss.bsc"

CLEAN : 
	-@erase ".\Debug\random_numerical_generation.mod"
	-@erase ".\Debug\cgauss.bsc"
	-@erase ".\Debug\define geometry.sbr"
	-@erase ".\Debug\shared random number generation.sbr"
	-@erase ".\Debug\prep default testing parameters.sbr"
	-@erase ".\Debug\scalars_to_be_shared.mod"
	-@erase ".\Debug\load restart point.sbr"
	-@erase ".\Debug\vectors_to_be_shared.mod"
	-@erase ".\Debug\matrices_to_be_shared.mod"
	-@erase ".\Debug\strings_to_be_shared.mod"
	-@erase ".\Debug\lower triangle duplication to matrix.sbr"
	-@erase ".\Debug\construct vector-half.sbr"
	-@erase ".\Debug\default program settings.sbr"
	-@erase ".\Debug\program_settings_to_be_shared.mod"
	-@erase ".\Debug\matrix blank.sbr"
	-@erase ".\Debug\gradient optimization.sbr"
	-@erase ".\Debug\construction of gradient.sbr"
	-@erase ".\Debug\jacobi matrix diagonalization.sbr"
	-@erase ".\Debug\cholesky condition.sbr"
	-@erase ".\Debug\guess random starting data.sbr"
	-@erase ".\Debug\calculate rayleigh quotient.sbr"
	-@erase ".\Debug\get current date & time.sbr"
	-@erase ".\Debug\construction of S & H matrices.sbr"
	-@erase ".\Debug\return value of the gradient.sbr"
	-@erase ".\Debug\simple set printout.sbr"
	-@erase ".\Debug\shared matrices.sbr"
	-@erase ".\Debug\electro-nuclear attraction element.sbr"
	-@erase ".\Debug\shared vectors.sbr"
	-@erase ".\Debug\print wavefunction stats.sbr"
	-@erase ".\Debug\calculate hamiltonian.sbr"
	-@erase ".\Debug\delectron repulsion element.sbr"
	-@erase ".\Debug\load default testing parameters.sbr"
	-@erase ".\Debug\solve 2 by 2.sbr"
	-@erase ".\Debug\solve eigenvalue problem.sbr"
	-@erase ".\Debug\print preliminary matter.sbr"
	-@erase ".\Debug\energy diagnosis.sbr"
	-@erase ".\Debug\shared settings.sbr"
	-@erase ".\Debug\screen & general file text output generation.sbr"
	-@erase ".\Debug\calculate hamiltonian gradient.sbr"
	-@erase ".\Debug\program core.sbr"
	-@erase ".\Debug\doverlap element.sbr"
	-@erase ".\Debug\shared scalars.sbr"
	-@erase ".\Debug\piped input processor.sbr"
	-@erase ".\Debug\return value of the function.sbr"
	-@erase ".\Debug\cholesky decondition.sbr"
	-@erase ".\Debug\calculate kinetic energy.sbr"
	-@erase ".\Debug\error function approximation.sbr"
	-@erase ".\Debug\electron repulsion element.sbr"
	-@erase ".\Debug\shared strings.sbr"
	-@erase ".\Debug\dkinetic energy element.sbr"
	-@erase ".\Debug\kinetic energy element.sbr"
	-@erase ".\Debug\delectro-nuclear attraction element.sbr"
	-@erase ".\Debug\matrix print (TAB).sbr"
	-@erase ".\Debug\overlap element.sbr"
	-@erase ".\Debug\cgauss.exe"
	-@erase ".\Debug\error function approximation.obj"
	-@erase ".\Debug\electron repulsion element.obj"
	-@erase ".\Debug\shared strings.obj"
	-@erase ".\Debug\dkinetic energy element.obj"
	-@erase ".\Debug\kinetic energy element.obj"
	-@erase ".\Debug\delectro-nuclear attraction element.obj"
	-@erase ".\Debug\matrix print (TAB).obj"
	-@erase ".\Debug\overlap element.obj"
	-@erase ".\Debug\define geometry.obj"
	-@erase ".\Debug\shared random number generation.obj"
	-@erase ".\Debug\prep default testing parameters.obj"
	-@erase ".\Debug\load restart point.obj"
	-@erase ".\Debug\lower triangle duplication to matrix.obj"
	-@erase ".\Debug\construct vector-half.obj"
	-@erase ".\Debug\default program settings.obj"
	-@erase ".\Debug\matrix blank.obj"
	-@erase ".\Debug\gradient optimization.obj"
	-@erase ".\Debug\construction of gradient.obj"
	-@erase ".\Debug\jacobi matrix diagonalization.obj"
	-@erase ".\Debug\cholesky condition.obj"
	-@erase ".\Debug\guess random starting data.obj"
	-@erase ".\Debug\calculate rayleigh quotient.obj"
	-@erase ".\Debug\get current date & time.obj"
	-@erase ".\Debug\construction of S & H matrices.obj"
	-@erase ".\Debug\return value of the gradient.obj"
	-@erase ".\Debug\simple set printout.obj"
	-@erase ".\Debug\shared matrices.obj"
	-@erase ".\Debug\electro-nuclear attraction element.obj"
	-@erase ".\Debug\shared vectors.obj"
	-@erase ".\Debug\print wavefunction stats.obj"
	-@erase ".\Debug\calculate hamiltonian.obj"
	-@erase ".\Debug\delectron repulsion element.obj"
	-@erase ".\Debug\load default testing parameters.obj"
	-@erase ".\Debug\solve 2 by 2.obj"
	-@erase ".\Debug\solve eigenvalue problem.obj"
	-@erase ".\Debug\print preliminary matter.obj"
	-@erase ".\Debug\energy diagnosis.obj"
	-@erase ".\Debug\shared settings.obj"
	-@erase ".\Debug\screen & general file text output generation.obj"
	-@erase ".\Debug\calculate hamiltonian gradient.obj"
	-@erase ".\Debug\program core.obj"
	-@erase ".\Debug\doverlap element.obj"
	-@erase ".\Debug\shared scalars.obj"
	-@erase ".\Debug\piped input processor.obj"
	-@erase ".\Debug\return value of the function.obj"
	-@erase ".\Debug\cholesky decondition.obj"
	-@erase ".\Debug\calculate kinetic energy.obj"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE F90 /Zi /I "Debug/" /c /nologo
# ADD F90 /Ox /Op /G5 /4R8 /4Yb /FR /Zi /I "Debug/" /c /nologo /Ob2 /4Mf
F90_PROJ=/Ox /Op /G5 /4R8 /4Yb /FR"Debug/" /Zi /I "Debug/" /c /nologo /Ob2\
 /Fo"Debug/" /Fd"Debug/cgauss.pdb" /4Mf 
F90_OBJS=.\Debug/
CPP_OBJS=.\Debug/
CPP_SBRS=
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/cgauss.bsc" 
BSC32_SBRS= \
	"$(INTDIR)/define geometry.sbr" \
	"$(INTDIR)/shared random number generation.sbr" \
	"$(INTDIR)/prep default testing parameters.sbr" \
	"$(INTDIR)/load restart point.sbr" \
	"$(INTDIR)/lower triangle duplication to matrix.sbr" \
	"$(INTDIR)/construct vector-half.sbr" \
	"$(INTDIR)/default program settings.sbr" \
	"$(INTDIR)/matrix blank.sbr" \
	"$(INTDIR)/gradient optimization.sbr" \
	"$(INTDIR)/construction of gradient.sbr" \
	"$(INTDIR)/jacobi matrix diagonalization.sbr" \
	"$(INTDIR)/cholesky condition.sbr" \
	"$(INTDIR)/guess random starting data.sbr" \
	"$(INTDIR)/calculate rayleigh quotient.sbr" \
	"$(INTDIR)/get current date & time.sbr" \
	"$(INTDIR)/construction of S & H matrices.sbr" \
	"$(INTDIR)/return value of the gradient.sbr" \
	"$(INTDIR)/simple set printout.sbr" \
	"$(INTDIR)/shared matrices.sbr" \
	"$(INTDIR)/electro-nuclear attraction element.sbr" \
	"$(INTDIR)/shared vectors.sbr" \
	"$(INTDIR)/print wavefunction stats.sbr" \
	"$(INTDIR)/calculate hamiltonian.sbr" \
	"$(INTDIR)/delectron repulsion element.sbr" \
	"$(INTDIR)/load default testing parameters.sbr" \
	"$(INTDIR)/solve 2 by 2.sbr" \
	"$(INTDIR)/solve eigenvalue problem.sbr" \
	"$(INTDIR)/print preliminary matter.sbr" \
	"$(INTDIR)/energy diagnosis.sbr" \
	"$(INTDIR)/shared settings.sbr" \
	"$(INTDIR)/screen & general file text output generation.sbr" \
	"$(INTDIR)/calculate hamiltonian gradient.sbr" \
	"$(INTDIR)/program core.sbr" \
	"$(INTDIR)/doverlap element.sbr" \
	"$(INTDIR)/shared scalars.sbr" \
	"$(INTDIR)/piped input processor.sbr" \
	"$(INTDIR)/return value of the function.sbr" \
	"$(INTDIR)/cholesky decondition.sbr" \
	"$(INTDIR)/calculate kinetic energy.sbr" \
	"$(INTDIR)/error function approximation.sbr" \
	"$(INTDIR)/electron repulsion element.sbr" \
	"$(INTDIR)/shared strings.sbr" \
	"$(INTDIR)/dkinetic energy element.sbr" \
	"$(INTDIR)/kinetic energy element.sbr" \
	"$(INTDIR)/delectro-nuclear attraction element.sbr" \
	"$(INTDIR)/matrix print (TAB).sbr" \
	"$(INTDIR)/overlap element.sbr"

"$(OUTDIR)\cgauss.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /pdb:none /debug /machine:I386
# SUBTRACT LINK32 /profile /map
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /pdb:none /debug\
 /machine:I386 /out:"$(OUTDIR)/cgauss.exe" 
LINK32_OBJS= \
	"$(INTDIR)/error function approximation.obj" \
	"$(INTDIR)/electron repulsion element.obj" \
	"$(INTDIR)/shared strings.obj" \
	"$(INTDIR)/dkinetic energy element.obj" \
	"$(INTDIR)/kinetic energy element.obj" \
	"$(INTDIR)/delectro-nuclear attraction element.obj" \
	"$(INTDIR)/matrix print (TAB).obj" \
	"$(INTDIR)/overlap element.obj" \
	"$(INTDIR)/define geometry.obj" \
	"$(INTDIR)/shared random number generation.obj" \
	"$(INTDIR)/prep default testing parameters.obj" \
	"$(INTDIR)/load restart point.obj" \
	"$(INTDIR)/lower triangle duplication to matrix.obj" \
	"$(INTDIR)/construct vector-half.obj" \
	"$(INTDIR)/default program settings.obj" \
	"$(INTDIR)/matrix blank.obj" \
	"$(INTDIR)/gradient optimization.obj" \
	"$(INTDIR)/construction of gradient.obj" \
	"$(INTDIR)/jacobi matrix diagonalization.obj" \
	"$(INTDIR)/cholesky condition.obj" \
	"$(INTDIR)/guess random starting data.obj" \
	"$(INTDIR)/calculate rayleigh quotient.obj" \
	"$(INTDIR)/get current date & time.obj" \
	"$(INTDIR)/construction of S & H matrices.obj" \
	"$(INTDIR)/return value of the gradient.obj" \
	"$(INTDIR)/simple set printout.obj" \
	"$(INTDIR)/shared matrices.obj" \
	"$(INTDIR)/electro-nuclear attraction element.obj" \
	"$(INTDIR)/shared vectors.obj" \
	"$(INTDIR)/print wavefunction stats.obj" \
	"$(INTDIR)/calculate hamiltonian.obj" \
	"$(INTDIR)/delectron repulsion element.obj" \
	"$(INTDIR)/load default testing parameters.obj" \
	"$(INTDIR)/solve 2 by 2.obj" \
	"$(INTDIR)/solve eigenvalue problem.obj" \
	"$(INTDIR)/print preliminary matter.obj" \
	"$(INTDIR)/energy diagnosis.obj" \
	"$(INTDIR)/shared settings.obj" \
	"$(INTDIR)/screen & general file text output generation.obj" \
	"$(INTDIR)/calculate hamiltonian gradient.obj" \
	"$(INTDIR)/program core.obj" \
	"$(INTDIR)/doverlap element.obj" \
	"$(INTDIR)/shared scalars.obj" \
	"$(INTDIR)/piped input processor.obj" \
	"$(INTDIR)/return value of the function.obj" \
	"$(INTDIR)/cholesky decondition.obj" \
	"$(INTDIR)/calculate kinetic energy.obj" \
	"..\..\..\..\..\..\MSDEV\LIB\MATHD.LIB"

"$(OUTDIR)\cgauss.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ENDIF 

CPP_PROJ=/ML /Fo"$(INTDIR)/" 

.c{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.cpp{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.cxx{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.c{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

.cpp{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

.cxx{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

################################################################################
# Begin Target

# Name "cgauss - Win32 Release"
# Name "cgauss - Win32 Debug"

!IF  "$(CFG)" == "cgauss - Win32 Release"

!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

!ENDIF 

################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\shared matrices.f90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

F90_MODOUT=\
	"matrices_to_be_shared"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\shared matrices.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\matrices_to_be_shared.mod" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

F90_MODOUT=\
	"matrices_to_be_shared"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\shared matrices.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\shared matrices.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\matrices_to_be_shared.mod" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\shared vectors.f90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

F90_MODOUT=\
	"vectors_to_be_shared"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\shared vectors.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\vectors_to_be_shared.mod" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

F90_MODOUT=\
	"vectors_to_be_shared"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\shared vectors.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\shared vectors.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\vectors_to_be_shared.mod" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\program core.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"


"$(INTDIR)\program core.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\program core.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\program core.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\shared scalars.f90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

F90_MODOUT=\
	"scalars_to_be_shared"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\shared scalars.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\scalars_to_be_shared.mod" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

F90_MODOUT=\
	"scalars_to_be_shared"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\shared scalars.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\shared scalars.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\scalars_to_be_shared.mod" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\shared settings.f90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

F90_MODOUT=\
	"program_settings_to_be_shared"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\shared settings.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\program_settings_to_be_shared.mod" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

F90_MODOUT=\
	"program_settings_to_be_shared"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\shared settings.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\shared settings.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\program_settings_to_be_shared.mod" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\piped input processor.f90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

NODEP_F90_PIPED=\
	".\Release\program_settings_to_be_shared.mod"\
	

"$(INTDIR)\piped input processor.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\program_settings_to_be_shared.mod"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

NODEP_F90_PIPED=\
	".\Debug\program_settings_to_be_shared.mod"\
	

BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\piped input processor.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\program_settings_to_be_shared.mod"
   $(BuildCmds)

"$(INTDIR)\piped input processor.sbr" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\program_settings_to_be_shared.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\shared strings.f90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

F90_MODOUT=\
	"strings_to_be_shared"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\shared strings.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\strings_to_be_shared.mod" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

F90_MODOUT=\
	"strings_to_be_shared"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\shared strings.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\shared strings.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\strings_to_be_shared.mod" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\guess random starting data.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"


"$(INTDIR)\guess random starting data.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\guess random starting data.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\guess random starting data.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\screen & general file text output generation.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"


"$(INTDIR)\screen & general file text output generation.obj" : $(SOURCE)\
 "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\screen & general file text output generation.obj" : $(SOURCE)\
 "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\screen & general file text output generation.sbr" : $(SOURCE)\
 "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\construction of S & H matrices.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

NODEP_F90_CONST=\
	".\Release\scalars_to_be_shared.mod"\
	

"$(INTDIR)\construction of S & H matrices.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

DEP_F90_CONST=\
	".\Debug\scalars_to_be_shared.mod"\
	

BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\construction of S & H matrices.obj" : $(SOURCE) $(DEP_F90_CONST)\
 "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"
   $(BuildCmds)

"$(INTDIR)\construction of S & H matrices.sbr" : $(SOURCE) $(DEP_F90_CONST)\
 "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\overlap element.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

NODEP_F90_OVERL=\
	".\Release\scalars_to_be_shared.mod"\
	".\Release\vectors_to_be_shared.mod"\
	".\Release\matrices_to_be_shared.mod"\
	".\Release\program_settings_to_be_shared.mod"\
	

"$(INTDIR)\overlap element.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod" "$(INTDIR)\vectors_to_be_shared.mod"\
 "$(INTDIR)\matrices_to_be_shared.mod"\
 "$(INTDIR)\program_settings_to_be_shared.mod"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

DEP_F90_OVERL=\
	".\Debug\scalars_to_be_shared.mod"\
	".\Debug\vectors_to_be_shared.mod"\
	".\Debug\matrices_to_be_shared.mod"\
	".\Debug\program_settings_to_be_shared.mod"\
	

BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\overlap element.obj" : $(SOURCE) $(DEP_F90_OVERL) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod" "$(INTDIR)\vectors_to_be_shared.mod"\
 "$(INTDIR)\matrices_to_be_shared.mod"\
 "$(INTDIR)\program_settings_to_be_shared.mod"
   $(BuildCmds)

"$(INTDIR)\overlap element.sbr" : $(SOURCE) $(DEP_F90_OVERL) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod" "$(INTDIR)\vectors_to_be_shared.mod"\
 "$(INTDIR)\matrices_to_be_shared.mod"\
 "$(INTDIR)\program_settings_to_be_shared.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\matrix print (TAB).F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"


"$(INTDIR)\matrix print (TAB).obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\matrix print (TAB).obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\matrix print (TAB).sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\calculate kinetic energy.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"


"$(INTDIR)\calculate kinetic energy.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\calculate kinetic energy.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\calculate kinetic energy.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\kinetic energy element.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"


"$(INTDIR)\kinetic energy element.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\kinetic energy element.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\kinetic energy element.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\define geometry.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"


"$(INTDIR)\define geometry.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\define geometry.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\define geometry.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\electro-nuclear attraction element.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"


"$(INTDIR)\electro-nuclear attraction element.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\electro-nuclear attraction element.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\electro-nuclear attraction element.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\error function approximation.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"


"$(INTDIR)\error function approximation.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\error function approximation.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\error function approximation.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\electron repulsion element.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

NODEP_F90_ELECT=\
	".\Release\scalars_to_be_shared.mod"\
	

"$(INTDIR)\electron repulsion element.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

DEP_F90_ELECT=\
	".\Debug\scalars_to_be_shared.mod"\
	

BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\electron repulsion element.obj" : $(SOURCE) $(DEP_F90_ELECT)\
 "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"
   $(BuildCmds)

"$(INTDIR)\electron repulsion element.sbr" : $(SOURCE) $(DEP_F90_ELECT)\
 "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\solve eigenvalue problem.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"


"$(INTDIR)\solve eigenvalue problem.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\solve eigenvalue problem.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\solve eigenvalue problem.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\construct vector-half.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"


"$(INTDIR)\construct vector-half.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\construct vector-half.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\construct vector-half.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\jacobi matrix diagonalization.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"


"$(INTDIR)\jacobi matrix diagonalization.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\jacobi matrix diagonalization.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\jacobi matrix diagonalization.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\lower triangle duplication to matrix.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"


"$(INTDIR)\lower triangle duplication to matrix.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\lower triangle duplication to matrix.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\lower triangle duplication to matrix.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\matrix blank.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"


"$(INTDIR)\matrix blank.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\matrix blank.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\matrix blank.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\print preliminary matter.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

NODEP_F90_PRINT=\
	".\Release\scalars_to_be_shared.mod"\
	".\Release\vectors_to_be_shared.mod"\
	".\Release\matrices_to_be_shared.mod"\
	".\Release\strings_to_be_shared.mod"\
	".\Release\program_settings_to_be_shared.mod"\
	

"$(INTDIR)\print preliminary matter.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod" "$(INTDIR)\vectors_to_be_shared.mod"\
 "$(INTDIR)\matrices_to_be_shared.mod" "$(INTDIR)\strings_to_be_shared.mod"\
 "$(INTDIR)\program_settings_to_be_shared.mod"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

DEP_F90_PRINT=\
	".\Debug\scalars_to_be_shared.mod"\
	".\Debug\vectors_to_be_shared.mod"\
	".\Debug\matrices_to_be_shared.mod"\
	".\Debug\strings_to_be_shared.mod"\
	".\Debug\program_settings_to_be_shared.mod"\
	

BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\print preliminary matter.obj" : $(SOURCE) $(DEP_F90_PRINT)\
 "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"\
 "$(INTDIR)\vectors_to_be_shared.mod" "$(INTDIR)\matrices_to_be_shared.mod"\
 "$(INTDIR)\strings_to_be_shared.mod"\
 "$(INTDIR)\program_settings_to_be_shared.mod"
   $(BuildCmds)

"$(INTDIR)\print preliminary matter.sbr" : $(SOURCE) $(DEP_F90_PRINT)\
 "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"\
 "$(INTDIR)\vectors_to_be_shared.mod" "$(INTDIR)\matrices_to_be_shared.mod"\
 "$(INTDIR)\strings_to_be_shared.mod"\
 "$(INTDIR)\program_settings_to_be_shared.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\load default testing parameters.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

NODEP_F90_LOAD_=\
	".\Release\vectors_to_be_shared.mod"\
	

"$(INTDIR)\load default testing parameters.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\vectors_to_be_shared.mod"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

DEP_F90_LOAD_=\
	".\Debug\vectors_to_be_shared.mod"\
	

BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\load default testing parameters.obj" : $(SOURCE) $(DEP_F90_LOAD_)\
 "$(INTDIR)" "$(INTDIR)\vectors_to_be_shared.mod"
   $(BuildCmds)

"$(INTDIR)\load default testing parameters.sbr" : $(SOURCE) $(DEP_F90_LOAD_)\
 "$(INTDIR)" "$(INTDIR)\vectors_to_be_shared.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\prep default testing parameters.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

NODEP_F90_PREP_=\
	".\Release\scalars_to_be_shared.mod"\
	

"$(INTDIR)\prep default testing parameters.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

DEP_F90_PREP_=\
	".\Debug\scalars_to_be_shared.mod"\
	

BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\prep default testing parameters.obj" : $(SOURCE) $(DEP_F90_PREP_)\
 "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"
   $(BuildCmds)

"$(INTDIR)\prep default testing parameters.sbr" : $(SOURCE) $(DEP_F90_PREP_)\
 "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\load restart point.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

NODEP_F90_LOAD_R=\
	".\Release\scalars_to_be_shared.mod"\
	".\Release\vectors_to_be_shared.mod"\
	".\Release\matrices_to_be_shared.mod"\
	".\Release\strings_to_be_shared.mod"\
	

"$(INTDIR)\load restart point.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod" "$(INTDIR)\vectors_to_be_shared.mod"\
 "$(INTDIR)\matrices_to_be_shared.mod" "$(INTDIR)\strings_to_be_shared.mod"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

DEP_F90_LOAD_R=\
	".\Debug\scalars_to_be_shared.mod"\
	".\Debug\vectors_to_be_shared.mod"\
	".\Debug\matrices_to_be_shared.mod"\
	".\Debug\strings_to_be_shared.mod"\
	

BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\load restart point.obj" : $(SOURCE) $(DEP_F90_LOAD_R) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod" "$(INTDIR)\vectors_to_be_shared.mod"\
 "$(INTDIR)\matrices_to_be_shared.mod" "$(INTDIR)\strings_to_be_shared.mod"
   $(BuildCmds)

"$(INTDIR)\load restart point.sbr" : $(SOURCE) $(DEP_F90_LOAD_R) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod" "$(INTDIR)\vectors_to_be_shared.mod"\
 "$(INTDIR)\matrices_to_be_shared.mod" "$(INTDIR)\strings_to_be_shared.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\cholesky condition.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

NODEP_F90_CHOLE=\
	".\Release\scalars_to_be_shared.mod"\
	".\Release\vectors_to_be_shared.mod"\
	

"$(INTDIR)\cholesky condition.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod" "$(INTDIR)\vectors_to_be_shared.mod"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

DEP_F90_CHOLE=\
	".\Debug\scalars_to_be_shared.mod"\
	".\Debug\vectors_to_be_shared.mod"\
	

BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\cholesky condition.obj" : $(SOURCE) $(DEP_F90_CHOLE) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod" "$(INTDIR)\vectors_to_be_shared.mod"
   $(BuildCmds)

"$(INTDIR)\cholesky condition.sbr" : $(SOURCE) $(DEP_F90_CHOLE) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod" "$(INTDIR)\vectors_to_be_shared.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\print wavefunction stats.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

NODEP_F90_PRINT_=\
	".\Release\scalars_to_be_shared.mod"\
	".\Release\vectors_to_be_shared.mod"\
	".\Release\matrices_to_be_shared.mod"\
	

"$(INTDIR)\print wavefunction stats.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod" "$(INTDIR)\vectors_to_be_shared.mod"\
 "$(INTDIR)\matrices_to_be_shared.mod"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

DEP_F90_PRINT_=\
	".\Debug\scalars_to_be_shared.mod"\
	".\Debug\vectors_to_be_shared.mod"\
	".\Debug\matrices_to_be_shared.mod"\
	

BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\print wavefunction stats.obj" : $(SOURCE) $(DEP_F90_PRINT_)\
 "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"\
 "$(INTDIR)\vectors_to_be_shared.mod" "$(INTDIR)\matrices_to_be_shared.mod"
   $(BuildCmds)

"$(INTDIR)\print wavefunction stats.sbr" : $(SOURCE) $(DEP_F90_PRINT_)\
 "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"\
 "$(INTDIR)\vectors_to_be_shared.mod" "$(INTDIR)\matrices_to_be_shared.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\gradient optimization.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

NODEP_F90_GRADI=\
	".\Release\scalars_to_be_shared.mod"\
	".\Release\vectors_to_be_shared.mod"\
	".\Release\matrices_to_be_shared.mod"\
	

"$(INTDIR)\gradient optimization.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod" "$(INTDIR)\vectors_to_be_shared.mod"\
 "$(INTDIR)\matrices_to_be_shared.mod"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

DEP_F90_GRADI=\
	".\Debug\scalars_to_be_shared.mod"\
	".\Debug\vectors_to_be_shared.mod"\
	".\Debug\matrices_to_be_shared.mod"\
	

BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\gradient optimization.obj" : $(SOURCE) $(DEP_F90_GRADI) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod" "$(INTDIR)\vectors_to_be_shared.mod"\
 "$(INTDIR)\matrices_to_be_shared.mod"
   $(BuildCmds)

"$(INTDIR)\gradient optimization.sbr" : $(SOURCE) $(DEP_F90_GRADI) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod" "$(INTDIR)\vectors_to_be_shared.mod"\
 "$(INTDIR)\matrices_to_be_shared.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\return value of the function.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

NODEP_F90_RETUR=\
	".\Release\scalars_to_be_shared.mod"\
	".\Release\vectors_to_be_shared.mod"\
	".\Release\matrices_to_be_shared.mod"\
	

"$(INTDIR)\return value of the function.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod" "$(INTDIR)\vectors_to_be_shared.mod"\
 "$(INTDIR)\matrices_to_be_shared.mod"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

DEP_F90_RETUR=\
	".\Debug\scalars_to_be_shared.mod"\
	".\Debug\vectors_to_be_shared.mod"\
	".\Debug\matrices_to_be_shared.mod"\
	

BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\return value of the function.obj" : $(SOURCE) $(DEP_F90_RETUR)\
 "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"\
 "$(INTDIR)\vectors_to_be_shared.mod" "$(INTDIR)\matrices_to_be_shared.mod"
   $(BuildCmds)

"$(INTDIR)\return value of the function.sbr" : $(SOURCE) $(DEP_F90_RETUR)\
 "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"\
 "$(INTDIR)\vectors_to_be_shared.mod" "$(INTDIR)\matrices_to_be_shared.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\shared random number generation.f90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

F90_MODOUT=\
	"random_numerical_generation"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\shared random number generation.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\random_numerical_generation.mod" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

F90_MODOUT=\
	"random_numerical_generation"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\shared random number generation.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\shared random number generation.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\random_numerical_generation.mod" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\get current date & time.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"


"$(INTDIR)\get current date & time.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\get current date & time.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\get current date & time.sbr" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\doverlap element.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

NODEP_F90_DOVER=\
	".\Release\scalars_to_be_shared.mod"\
	

"$(INTDIR)\doverlap element.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

DEP_F90_DOVER=\
	".\Debug\scalars_to_be_shared.mod"\
	

BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\doverlap element.obj" : $(SOURCE) $(DEP_F90_DOVER) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod"
   $(BuildCmds)

"$(INTDIR)\doverlap element.sbr" : $(SOURCE) $(DEP_F90_DOVER) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\construction of gradient.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

NODEP_F90_CONSTR=\
	".\Release\scalars_to_be_shared.mod"\
	".\Release\matrices_to_be_shared.mod"\
	

"$(INTDIR)\construction of gradient.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod" "$(INTDIR)\matrices_to_be_shared.mod"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

DEP_F90_CONSTR=\
	".\Debug\scalars_to_be_shared.mod"\
	".\Debug\matrices_to_be_shared.mod"\
	

BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\construction of gradient.obj" : $(SOURCE) $(DEP_F90_CONSTR)\
 "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"\
 "$(INTDIR)\matrices_to_be_shared.mod"
   $(BuildCmds)

"$(INTDIR)\construction of gradient.sbr" : $(SOURCE) $(DEP_F90_CONSTR)\
 "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"\
 "$(INTDIR)\matrices_to_be_shared.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\calculate rayleigh quotient.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

NODEP_F90_CALCU=\
	".\Release\scalars_to_be_shared.mod"\
	".\Release\vectors_to_be_shared.mod"\
	".\Release\matrices_to_be_shared.mod"\
	

"$(INTDIR)\calculate rayleigh quotient.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod" "$(INTDIR)\vectors_to_be_shared.mod"\
 "$(INTDIR)\matrices_to_be_shared.mod"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

DEP_F90_CALCU=\
	".\Debug\scalars_to_be_shared.mod"\
	".\Debug\vectors_to_be_shared.mod"\
	".\Debug\matrices_to_be_shared.mod"\
	

BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\calculate rayleigh quotient.obj" : $(SOURCE) $(DEP_F90_CALCU)\
 "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"\
 "$(INTDIR)\vectors_to_be_shared.mod" "$(INTDIR)\matrices_to_be_shared.mod"
   $(BuildCmds)

"$(INTDIR)\calculate rayleigh quotient.sbr" : $(SOURCE) $(DEP_F90_CALCU)\
 "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"\
 "$(INTDIR)\vectors_to_be_shared.mod" "$(INTDIR)\matrices_to_be_shared.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\calculate hamiltonian gradient.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

NODEP_F90_CALCUL=\
	".\Release\scalars_to_be_shared.mod"\
	".\Release\vectors_to_be_shared.mod"\
	".\Release\matrices_to_be_shared.mod"\
	

"$(INTDIR)\calculate hamiltonian gradient.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod" "$(INTDIR)\vectors_to_be_shared.mod"\
 "$(INTDIR)\matrices_to_be_shared.mod"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

DEP_F90_CALCUL=\
	".\Debug\scalars_to_be_shared.mod"\
	".\Debug\vectors_to_be_shared.mod"\
	".\Debug\matrices_to_be_shared.mod"\
	

BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\calculate hamiltonian gradient.obj" : $(SOURCE) $(DEP_F90_CALCUL)\
 "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"\
 "$(INTDIR)\vectors_to_be_shared.mod" "$(INTDIR)\matrices_to_be_shared.mod"
   $(BuildCmds)

"$(INTDIR)\calculate hamiltonian gradient.sbr" : $(SOURCE) $(DEP_F90_CALCUL)\
 "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"\
 "$(INTDIR)\vectors_to_be_shared.mod" "$(INTDIR)\matrices_to_be_shared.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\dkinetic energy element.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

NODEP_F90_DKINE=\
	".\Release\scalars_to_be_shared.mod"\
	

"$(INTDIR)\dkinetic energy element.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

DEP_F90_DKINE=\
	".\Debug\scalars_to_be_shared.mod"\
	

BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\dkinetic energy element.obj" : $(SOURCE) $(DEP_F90_DKINE)\
 "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"
   $(BuildCmds)

"$(INTDIR)\dkinetic energy element.sbr" : $(SOURCE) $(DEP_F90_DKINE)\
 "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\delectron repulsion element.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

NODEP_F90_DELEC=\
	".\Release\scalars_to_be_shared.mod"\
	

"$(INTDIR)\delectron repulsion element.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

DEP_F90_DELEC=\
	".\Debug\scalars_to_be_shared.mod"\
	

BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\delectron repulsion element.obj" : $(SOURCE) $(DEP_F90_DELEC)\
 "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"
   $(BuildCmds)

"$(INTDIR)\delectron repulsion element.sbr" : $(SOURCE) $(DEP_F90_DELEC)\
 "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\delectro-nuclear attraction element.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

NODEP_F90_DELECT=\
	".\Release\scalars_to_be_shared.mod"\
	

"$(INTDIR)\delectro-nuclear attraction element.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

DEP_F90_DELECT=\
	".\Debug\scalars_to_be_shared.mod"\
	

BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\delectro-nuclear attraction element.obj" : $(SOURCE)\
 $(DEP_F90_DELECT) "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"
   $(BuildCmds)

"$(INTDIR)\delectro-nuclear attraction element.sbr" : $(SOURCE)\
 $(DEP_F90_DELECT) "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\return value of the gradient.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

NODEP_F90_RETURN=\
	".\Release\scalars_to_be_shared.mod"\
	".\Release\vectors_to_be_shared.mod"\
	".\Release\matrices_to_be_shared.mod"\
	".\Release\strings_to_be_shared.mod"\
	

"$(INTDIR)\return value of the gradient.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod" "$(INTDIR)\vectors_to_be_shared.mod"\
 "$(INTDIR)\matrices_to_be_shared.mod" "$(INTDIR)\strings_to_be_shared.mod"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

DEP_F90_RETURN=\
	".\Debug\scalars_to_be_shared.mod"\
	".\Debug\vectors_to_be_shared.mod"\
	".\Debug\matrices_to_be_shared.mod"\
	".\Debug\strings_to_be_shared.mod"\
	

BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\return value of the gradient.obj" : $(SOURCE) $(DEP_F90_RETURN)\
 "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"\
 "$(INTDIR)\vectors_to_be_shared.mod" "$(INTDIR)\matrices_to_be_shared.mod"\
 "$(INTDIR)\strings_to_be_shared.mod"
   $(BuildCmds)

"$(INTDIR)\return value of the gradient.sbr" : $(SOURCE) $(DEP_F90_RETURN)\
 "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"\
 "$(INTDIR)\vectors_to_be_shared.mod" "$(INTDIR)\matrices_to_be_shared.mod"\
 "$(INTDIR)\strings_to_be_shared.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\energy diagnosis.f90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

DEP_F90_ENERG=\
	".\Release\scalars_to_be_shared.mod"\
	".\Release\vectors_to_be_shared.mod"\
	".\Release\matrices_to_be_shared.mod"\
	

"$(INTDIR)\energy diagnosis.obj" : $(SOURCE) $(DEP_F90_ENERG) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod" "$(INTDIR)\vectors_to_be_shared.mod"\
 "$(INTDIR)\matrices_to_be_shared.mod"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

DEP_F90_ENERG=\
	".\Debug\scalars_to_be_shared.mod"\
	".\Debug\vectors_to_be_shared.mod"\
	".\Debug\matrices_to_be_shared.mod"\
	

BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\energy diagnosis.obj" : $(SOURCE) $(DEP_F90_ENERG) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod" "$(INTDIR)\vectors_to_be_shared.mod"\
 "$(INTDIR)\matrices_to_be_shared.mod"
   $(BuildCmds)

"$(INTDIR)\energy diagnosis.sbr" : $(SOURCE) $(DEP_F90_ENERG) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod" "$(INTDIR)\vectors_to_be_shared.mod"\
 "$(INTDIR)\matrices_to_be_shared.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\solve 2 by 2.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

DEP_F90_SOLVE=\
	".\Release\scalars_to_be_shared.mod"\
	".\Release\vectors_to_be_shared.mod"\
	

"$(INTDIR)\solve 2 by 2.obj" : $(SOURCE) $(DEP_F90_SOLVE) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod" "$(INTDIR)\vectors_to_be_shared.mod"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

DEP_F90_SOLVE=\
	".\Debug\scalars_to_be_shared.mod"\
	".\Debug\vectors_to_be_shared.mod"\
	

BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\solve 2 by 2.obj" : $(SOURCE) $(DEP_F90_SOLVE) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod" "$(INTDIR)\vectors_to_be_shared.mod"
   $(BuildCmds)

"$(INTDIR)\solve 2 by 2.sbr" : $(SOURCE) $(DEP_F90_SOLVE) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod" "$(INTDIR)\vectors_to_be_shared.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\MSDEV\LIB\MATHD.LIB

!IF  "$(CFG)" == "cgauss - Win32 Release"

!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\simple set printout.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

DEP_F90_SIMPL=\
	".\Release\scalars_to_be_shared.mod"\
	".\Release\vectors_to_be_shared.mod"\
	

"$(INTDIR)\simple set printout.obj" : $(SOURCE) $(DEP_F90_SIMPL) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod" "$(INTDIR)\vectors_to_be_shared.mod"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

DEP_F90_SIMPL=\
	".\Debug\scalars_to_be_shared.mod"\
	".\Debug\vectors_to_be_shared.mod"\
	

BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\simple set printout.obj" : $(SOURCE) $(DEP_F90_SIMPL) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod" "$(INTDIR)\vectors_to_be_shared.mod"
   $(BuildCmds)

"$(INTDIR)\simple set printout.sbr" : $(SOURCE) $(DEP_F90_SIMPL) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod" "$(INTDIR)\vectors_to_be_shared.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\cholesky decondition.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

DEP_F90_CHOLES=\
	".\Release\scalars_to_be_shared.mod"\
	".\Release\vectors_to_be_shared.mod"\
	

"$(INTDIR)\cholesky decondition.obj" : $(SOURCE) $(DEP_F90_CHOLES) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod" "$(INTDIR)\vectors_to_be_shared.mod"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

DEP_F90_CHOLES=\
	".\Debug\scalars_to_be_shared.mod"\
	".\Debug\vectors_to_be_shared.mod"\
	

BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\cholesky decondition.obj" : $(SOURCE) $(DEP_F90_CHOLES) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod" "$(INTDIR)\vectors_to_be_shared.mod"
   $(BuildCmds)

"$(INTDIR)\cholesky decondition.sbr" : $(SOURCE) $(DEP_F90_CHOLES) "$(INTDIR)"\
 "$(INTDIR)\scalars_to_be_shared.mod" "$(INTDIR)\vectors_to_be_shared.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\calculate hamiltonian.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

DEP_F90_CALCULA=\
	".\Release\scalars_to_be_shared.mod"\
	".\Release\vectors_to_be_shared.mod"\
	".\Release\matrices_to_be_shared.mod"\
	

"$(INTDIR)\calculate hamiltonian.obj" : $(SOURCE) $(DEP_F90_CALCULA)\
 "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"\
 "$(INTDIR)\vectors_to_be_shared.mod" "$(INTDIR)\matrices_to_be_shared.mod"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

DEP_F90_CALCULA=\
	".\Debug\scalars_to_be_shared.mod"\
	".\Debug\vectors_to_be_shared.mod"\
	".\Debug\matrices_to_be_shared.mod"\
	

BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\calculate hamiltonian.obj" : $(SOURCE) $(DEP_F90_CALCULA)\
 "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"\
 "$(INTDIR)\vectors_to_be_shared.mod" "$(INTDIR)\matrices_to_be_shared.mod"
   $(BuildCmds)

"$(INTDIR)\calculate hamiltonian.sbr" : $(SOURCE) $(DEP_F90_CALCULA)\
 "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"\
 "$(INTDIR)\vectors_to_be_shared.mod" "$(INTDIR)\matrices_to_be_shared.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\
"\users\dgilmore\Graduate.School\Adamowicz.Research\cgauss\default program settings.F90"

!IF  "$(CFG)" == "cgauss - Win32 Release"

DEP_F90_DEFAU=\
	".\Release\scalars_to_be_shared.mod"\
	".\Release\strings_to_be_shared.mod"\
	".\Release\program_settings_to_be_shared.mod"\
	

"$(INTDIR)\default program settings.obj" : $(SOURCE) $(DEP_F90_DEFAU)\
 "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"\
 "$(INTDIR)\strings_to_be_shared.mod"\
 "$(INTDIR)\program_settings_to_be_shared.mod"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "cgauss - Win32 Debug"

DEP_F90_DEFAU=\
	".\Debug\scalars_to_be_shared.mod"\
	".\Debug\strings_to_be_shared.mod"\
	".\Debug\program_settings_to_be_shared.mod"\
	

BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\default program settings.obj" : $(SOURCE) $(DEP_F90_DEFAU)\
 "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"\
 "$(INTDIR)\strings_to_be_shared.mod"\
 "$(INTDIR)\program_settings_to_be_shared.mod"
   $(BuildCmds)

"$(INTDIR)\default program settings.sbr" : $(SOURCE) $(DEP_F90_DEFAU)\
 "$(INTDIR)" "$(INTDIR)\scalars_to_be_shared.mod"\
 "$(INTDIR)\strings_to_be_shared.mod"\
 "$(INTDIR)\program_settings_to_be_shared.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
# End Target
# End Project
################################################################################
