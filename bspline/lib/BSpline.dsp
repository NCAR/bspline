# Microsoft Developer Studio Project File - Name="BSpline" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Dynamic-Link Library" 0x0102

CFG=BSpline - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "BSpline.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "BSpline.mak" CFG="BSpline - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "BSpline - Win32 Release" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE "BSpline - Win32 Debug" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
MTL=midl.exe
RSC=rc.exe

!IF  "$(CFG)" == "BSpline - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "BSpline_"
# PROP BASE Intermediate_Dir "BSpline_"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "c:/Binaries/BSpline/Release/Lib"
# PROP Intermediate_Dir "c:/Binaries/BSpline/Release/Obj"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MT /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /YX /FD /c
# ADD CPP /nologo /MD /W3 /GX /I "BSpline" /I "../BSpline" /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /YX /FD /c
# ADD BASE MTL /nologo /D "NDEBUG" /mktyplib203 /o "NUL" /win32
# ADD MTL /nologo /D "NDEBUG" /mktyplib203 /o "NUL" /win32
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /dll /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /dll /machine:I386
# Begin Special Build Tool
TargetPath=\Binaries\BSpline\Release\Lib\BSpline.dll
SOURCE="$(InputPath)"
PostBuild_Cmds=copy $(TargetPath)  c:\winnt\system32
# End Special Build Tool

!ELSEIF  "$(CFG)" == "BSpline - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "BSpline0"
# PROP BASE Intermediate_Dir "BSpline0"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "c:/Binaries/BSpline/Debug/Lib"
# PROP Intermediate_Dir "c:/Binaries/BSpline/Debug/Obj"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MTd /W3 /Gm /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /FD /c
# ADD CPP /MDd /W3 /Gm /GX /ZI /Od /I "BSpline" /I "../BSpline" /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /FD /c
# SUBTRACT CPP /nologo
# ADD BASE MTL /nologo /D "_DEBUG" /mktyplib203 /o "NUL" /win32
# ADD MTL /D "_DEBUG" /mktyplib203 /o "NUL" /win32
# SUBTRACT MTL /nologo
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /dll /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /subsystem:windows /dll /debug /machine:I386 /pdbtype:sept
# SUBTRACT LINK32 /nologo /verbose
# Begin Special Build Tool
TargetPath=\Binaries\BSpline\Debug\Lib\BSpline.dll
SOURCE="$(InputPath)"
PostBuild_Desc=Copy dll to system folder
PostBuild_Cmds=copy $(TargetPath)  c:\winnt\system32
# End Special Build Tool

!ENDIF 

# Begin Target

# Name "BSpline - Win32 Release"
# Name "BSpline - Win32 Debug"
# Begin Source File

SOURCE=.\BandedMatrix.h
# End Source File
# Begin Source File

SOURCE=.\BSpline.cpp
# End Source File
# Begin Source File

SOURCE=..\BSpline\BSpline.h
# End Source File
# Begin Source File

SOURCE=.\BSplineD.cpp
# End Source File
# Begin Source File

SOURCE=.\BSplineF.cpp
# End Source File
# End Target
# End Project
