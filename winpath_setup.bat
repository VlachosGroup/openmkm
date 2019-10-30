@echo OFF
echo Setting up environment for OpenMKM
set src=..\Data
for /f "tokens=*" %%a in ("%~dp0.\%src%") do set "src=%%~fa"
( ENDLOCAL & REM RETURN VALUES
    set cdata=%src%
)

%comspec% /k "set CANTERA_DATA=%cdata%;&& cd /d %HOMEPATH% && set PATH=%PATH%;%~dp0"

