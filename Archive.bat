@ECHO OFF

@REM================================================================================================================================
:SetVars

@SET ARCHPATH=Archive
@SET PROGNAME=TurbSim
@SET ARCHNAME=TurbSim_v%1

::=======================================================================================================
IF "%COMPUTERNAME%"=="APLATT-21846S" GOTO APLATT-21846S
IF "%COMPUTERNAME%"=="BJONKMAN-23080S" GOTO BJONKMAN-23080S

:APLATT-21846S
@SET WINZIP="C:\Program Files (x86)\WinZip\WZZip"
@SET WINZIPSE="C:\Program Files (x86)\WinZip Self-Extractor\wzipse32.exe"
@SET SEVENZIP="C:\Program Files\7-Zip\7z.exe"
GOTO CheckSyntax

:BJONKMAN-23080S
@SET WINZIP="C:\Program Files (x86)\WinZip\WZZip"
@SET WINZIPSE="C:\Program Files (x86)\WinZip Self-Extractor\WZIPSE22\wzipse32.exe"
@SET SEVENZIP="C:\Program Files\7-Zip\7z.exe"
GOTO CheckSyntax

::=======================================================================================================


:CheckSyntax
IF NOT "%1"==""  GOTO DeleteOld

@ECHO 
@ECHO  The syntax for creating an archive is "Archive <version>"
@ECHO.
@ECHO  Example:  "archive 2.01.00"

@GOTO Done

:DeleteOld
@IF EXIST ARCHTMP.zip DEL ARCHTMP.zip
@IF EXIST ARCHTMP.exe DEL ARCHTMP.exe
@IF EXIST ARCHTMP.tar DEL ARCHTMP.tar
@IF EXIST ARCHTMP.tar.gz DEL ARCHTMP.tar.gz

rem CALL CopyFilesForRelease.bat

:DoIt
@ECHO.
@ECHO --------------------------------------------------------------------------------------
@ECHO Archiving %PROGNAME% for general distribution on Windows.
@ECHO --------------------------------------------------------------------------------------
@ECHO.

@%WINZIP% -a -o -P ARCHTMP @ArcFiles.txt @ArcWin.txt
@%WINZIPSE% ARCHTMP.zip -d. -y -win32 -le -overwrite -st"Unzipping %PROGNAME%" -m Disclaimer.txt

@COPY ARCHTMP.exe %ARCHPATH%\%ARCHNAME%.exe
@DEL ARCHTMP.zip, ARCHTMP.exe

:: @ECHO.
:: @ECHO -------------------------------------------------------------------------------------
:: @ECHO Archiving %PROGNAME% for maintenance and internal use (including certification tests).
:: @ECHO --------------------------------------------------------------------------------------
:: @ECHO.
::
:: @%WINZIP% -a -o -P ARCHTMP @ArcFiles.txt @ArcWin.txt @ArcMaint.txt
:: @%WINZIPSE% ARCHTMP.zip -d. -y -win32 -le -overwrite -st"Unzipping %PROGNAME%" -m Disclaimer.txt
::
:: @COPY ARCHTMP.exe %ARCHPATH%\%ARCHNAME%_all.exe
:: @DEL ARCHTMP.zip, ARCHTMP.exe


@ECHO --------------------------------------------------------------------------------------
@ECHO Archiving %PROGNAME% for general distribution (tar.gz).
@ECHO --------------------------------------------------------------------------------------
@ECHO.
@rem first create a tar file, then compress it (gzip allows only one file)
@%SEVENZIP% a -ttar ARCHTMP @ArcFiles.txt
@%SEVENZIP% a -tgzip ARCHTMP.tar.gz ARCHTMP.tar
@COPY ARCHTMP.tar.gz %ARCHPATH%\%ARCHNAME%.tar.gz
@DEL ARCHTMP.tar, ARCHTMP.tar.gz



@REM================================================================================================================================
:Done
@REM================================================================================================================================
@SET ARCHNAME=
@SET ARCHPATH=
@SET PROGNAME=
@SET WINZIP=
@SET WINZIPSE=
@SET SEVENZIP=
