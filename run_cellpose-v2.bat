@echo off
setlocal enabledelayedexpansion

::———————— Settings ————————
:: (1) Where your raw images live:
set INPUT_DIR=C:\Users\user\Desktop\input_cellpose

:: (2) Where masks + overlays will get written:
set MASK_DIR=%INPUT_DIR%\masks

:: (3) Path to your Python .py quantification script:
set PYTHON_EXE=python
set QUANT_SCRIPT=C:\Users\user\Desktop\scripts_and_notes_cellpose\cellpose_quantify.py

:: (4) Path to your log file (captures everything):
set LOG_FILE=%INPUT_DIR%\cellpose_run_log.txt

:: Make sure the log file exists (or gets created):
if not exist "%LOG_FILE%" (
    rem Create an empty file and then echo a header
    > "%LOG_FILE%" echo ===== Cellpose run log started: %DATE% %TIME% =====
)

::———————— Interactive Prompts ————————

:: (A) Which channel to segment on?
echo.
echo Which channel do you want to use for segmentation?
echo    [0] Grayscale RGB
echo    [1] Red
echo    [2] Green
echo    [3] Blue
set /p channel="Enter channel number (0-3): "

:: (B) Which Cellpose model?
echo.
echo Which Cellpose model do you want to use?
echo    [1] cyto
echo    [2] nuclei
set /p model_choice="Enter model number (1 or 2): "
if "%model_choice%"=="1" (
    set MODEL=cyto
) else (
    set MODEL=nuclei
)

:: (C) Auto vs. manual diameter
echo.
echo Do you want Cellpose to auto-estimate diameter, or specify it yourself?
echo    [A] Auto (default auto-diameter)
echo    [M] Manual (you will be asked to type in a pixel value)
:diam_choice
set /p diam_choice="Enter A or M: "
if /I "%diam_choice%"=="A" (
    set DIAMETER=0
) else (
    if /I "%diam_choice%"=="M" (
        set /p DIAMETER="Enter expected diameter (in pixels): "
    ) else (
        echo Invalid choice—please type A or M.
        goto diam_choice
    )
)

::———————— Make sure output folder exists ————————
if not exist "%MASK_DIR%" (
    mkdir "%MASK_DIR%"
    echo Created directory "%MASK_DIR%" >> "%LOG_FILE%" 2>&1
)

::———————— Run Cellpose with redirection ————————
echo. >> "%LOG_FILE%" 2>&1
echo Running Cellpose on %INPUT_DIR% with: >> "%LOG_FILE%" 2>&1
echo    channel = %channel%           >> "%LOG_FILE%" 2>&1
echo    model   = %MODEL%             >> "%LOG_FILE%" 2>&1
echo    diameter= %DIAMETER%          >> "%LOG_FILE%" 2>&1
echo --------------------------------------- >> "%LOG_FILE%" 2>&1

cellpose --dir "%INPUT_DIR%" --pretrained_model %MODEL% --chan %channel% --diameter %DIAMETER% --savedir "%MASK_DIR%" --save_png --verbose >> "%LOG_FILE%" 2>&1
if errorlevel 1 (
    echo [ERROR] Cellpose returned a non-zero exit code. Check the log. >> "%LOG_FILE%" 2>&1
    goto end
)

::———————— Run Python quantification (append to same log) ————————
echo.                                   >> "%LOG_FILE%" 2>&1
echo Running quantification script...    >> "%LOG_FILE%" 2>&1
echo --------------------------------------- >> "%LOG_FILE%" 2>&1

:: FIX: Change the second argument to %INPUT_DIR%
%PYTHON_EXE% "%QUANT_SCRIPT%" "%INPUT_DIR%" "%INPUT_DIR%" >> "%LOG_FILE%" 2>&1
if errorlevel 1 (
    echo [ERROR] Python quantification returned a non-zero exit code. Check the log. >> "%LOG_FILE%" 2>&1
    goto end
)

::———————— Done —————————————
:end
echo.                                   >> "%LOG_FILE%" 2>&1
echo ===== Completed at %DATE% %TIME% ===== >> "%LOG_FILE%" 2>&1

endlocal