@echo off
setlocal
if defined VCVARS64 (
    call "%VCVARS64%"
    if errorlevel 1 exit /b 1
)
where cl.exe >nul 2>nul
if errorlevel 1 (
    echo MSVC not found. Use an x64 Native Tools Command Prompt or set VCVARS64.
    exit /b 1
)
if /I not "%VSCMD_ARG_TGT_ARCH%"=="x64" (
    echo Use the x64 MSVC environment.
    exit /b 1
)
pushd "%~dp0"
if errorlevel 1 exit /b 1
if not exist build mkdir build
if not exist build (
    popd
    exit /b 1
)
cl /nologo /LD /O2 /fp:strict /EHsc /W4 /utf-8 /Fo:build\ /Fe:build\recovered_primitives.dll src\recovered_primitives.cpp src\recovered_motor.cpp src\recovered_initialization.cpp src\recovered_grid.cpp src\recovered_commands.cpp src\recovered_parameters.cpp src\recovered_controller.cpp /link /IMPLIB:build\recovered_primitives.lib
set "BUILD_RESULT=%errorlevel%"
popd
exit /b %BUILD_RESULT%
