@echo off
setlocal
rem ============================================================================================
rem  One-click build: configure + compile (single-config generator, Ninja).
rem
rem  Prerequisites (either one):
rem    1) Run this script from "x64 Native Tools Command Prompt for VS 2022" -- cl.exe already in PATH;
rem    2) or set VCVARS to the path of vcvars64.bat first, e.g.:
rem         set VCVARS=C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat
rem
rem  The SGGK SDK is located automatically by cmake/FindSGGK.cmake
rem  (-DSGGK_SDK_DIR= / env SGGK_ROOT / ../SGGK / third_party/SGGK / ./SGGK; see README).
rem  Arguments after this script name are forwarded to cmake as-is, e.g.:
rem         build.bat -DSGGK_SDK_DIR=D:\somewhere\SGGK
rem
rem  NOTE: keep this file ASCII-only with CRLF line endings -- cmd.exe parses batch files using
rem        the system ANSI code page, so UTF-8 text or LF-only endings break command parsing.
rem ============================================================================================

set "ROOT=%~dp0"
if "%ROOT:~-1%"=="\" set "ROOT=%ROOT:~0,-1%"

where cl.exe >nul 2>nul
if errorlevel 1 (
    if defined VCVARS (
        echo [build] cl.exe not in PATH; calling "%VCVARS%"
        call "%VCVARS%"
        if errorlevel 1 (
            echo [build] failed to initialize MSVC environment via VCVARS
            exit /b 1
        )
    ) else (
        echo [build] cl.exe not found and VCVARS is not set.
        echo [build] Run this script from "x64 Native Tools Command Prompt for VS 2022",
        echo [build] or set VCVARS to the path of vcvars64.bat and retry.
        exit /b 1
    )
)

cmake -S "%ROOT%" -B "%ROOT%\build" -G Ninja -DCMAKE_BUILD_TYPE=Release ^
      -DCMAKE_CXX_FLAGS="/utf-8 /EHsc /MD" %*
if errorlevel 1 (
    echo [build] cmake configure failed
    exit /b 1
)

cmake --build "%ROOT%\build"
if errorlevel 1 (
    echo [build] cmake build failed
    exit /b 1
)

exit /b 0
