@echo off
setlocal enabledelayedexpansion

:: ─────────────────── Settings ───────────────────
:: (1) Where your raw images live:
set "INPUT_DIR=C:\Users\user\Desktop\input_cellpose_v5"

:: (2) Where output folder and masks will get written:
set "OUTPUT_DIR=%INPUT_DIR%\output"
set "MASK_DIR=%OUTPUT_DIR%\masks"
set "SARCASM_DIR=%OUTPUT_DIR%\masks_for_sarcasm"

:: (3) Path to your Python .py quantification script:
set "PYTHON_EXE=python"
set "QUANT_SCRIPT=C:\Users\user\Desktop\scripts_and_notes_cellpose\cellpose_quantify_v5.py"

:: (4) Path to your log file (captures everything):
set "LOG_FILE=%OUTPUT_DIR%\cellpose_run_log.txt"

:: (5) Path to save previous settings (in same folder as Python script):
for %%F in ("%QUANT_SCRIPT%") do set "SCRIPT_DIR=%%~dpF"
set "SETTINGS_FILE=%SCRIPT_DIR%cellpose_previous_settings.txt"

:: ─────────────────── Make sure the output directories exist ───────────────────
for %%D in ("%OUTPUT_DIR%" "%MASK_DIR%" "%SARCASM_DIR%") do (
    if not exist "%%~D" (
        mkdir "%%~D"
    )
)

:: ─────────────────── Make sure the log file exists ───────────────────
if not exist "%LOG_FILE%" (
    > "%LOG_FILE%" echo ===== Cellpose run log started: %DATE% %TIME% =====
)

:: ─────────────────── Interactive Prompts ───────────────────
:get_channel
echo.
echo Which channel do you want to use for segmentation?
echo    [0] Grayscale RGB
echo    [1] Red
echo    [2] Green
echo    [3] Blue
echo    [4] Re-run with previous settings
set /p channel="Enter choice (0-4): "

:: Re-use previous settings?  ───────────────────
if "%channel%"=="4" (
    if exist "%SETTINGS_FILE%" (
        echo.
        echo Loading previous settings...
        for /f "usebackq tokens=1,2 delims==" %%a in ("%SETTINGS_FILE%") do (
            if "%%a"=="channel"    set channel=%%b
            if "%%a"=="model"      set MODEL=%%b
            if "%%a"=="diameter"   set DIAMETER=%%b
            if "%%a"=="diam_mode"  set DIAM_MODE=%%b
            if "%%a"=="lwr_filter" set LWR_FILTER=%%b
            if "%%a"=="lwr_cutoff" set LWR_CUTOFF=%%b
        )
        goto show_loaded
    ) else (
        echo.
        echo No previous settings found.  Please enter settings manually.
        goto get_channel
    )
)

:: Channel entered manually ───────────────────
if "%channel%"=="0" (set CH_DESC=Grayscale RGB)
if "%channel%"=="1" (set CH_DESC=Red)
if "%channel%"=="2" (set CH_DESC=Green)
if "%channel%"=="3" (set CH_DESC=Blue)

:: ─────────────────── Choose model ───────────────────
:choose_model
echo.
echo Which Cellpose model?
echo    [1] cyto
echo    [2] nuclei
set /p model_choice="Enter 1 or 2: "
if "%model_choice%"=="1" (set MODEL=cyto) else if "%model_choice%"=="2" (set MODEL=nuclei) else (
    echo Invalid choice.
    goto choose_model
)

:: ─────────────────── Diameter ───────────────────
:diam_choice
echo.
echo Diameter options:
echo    [1] Auto-estimate
echo    [2] Manual
set /p diam_choice="Enter 1 or 2: "
if "%diam_choice%"=="1" (
    set DIAM_MODE=auto
    set DIAMETER=0
) else if "%diam_choice%"=="2" (
    set DIAM_MODE=manual
    set /p DIAMETER="Enter expected diameter (pixels): "
) else (
    echo Invalid choice.
    goto diam_choice
)

:: ─────────────────── Length/Width filter ───────────────────
:lwr_choice
echo.
echo Apply length/width ratio filter?
echo    [1] Yes
echo    [2] No
set /p lwr_choice="Enter 1 or 2: "
if "%lwr_choice%"=="1" (
    set LWR_FILTER=1
    set /p LWR_CUTOFF="Enter cutoff (e.g., 3 for adult CMs): "
    if "%LWR_CUTOFF%"=="" set LWR_CUTOFF=3
) else if "%lwr_choice%"=="2" (
    set LWR_FILTER=0
    set LWR_CUTOFF=0
) else (
    echo Invalid choice.
    goto lwr_choice
)

:show_loaded
echo.
echo --------------------- Settings to be used ---------------------------
echo Channel       : %channel%  (%CH_DESC%)
echo Model         : %MODEL%
echo Diameter mode : %DIAM_MODE%
echo Diameter      : %DIAMETER%
echo L/W filter    : %LWR_FILTER%  (cutoff=%LWR_CUTOFF%)
echo ---------------------------------------------------------------------

:: ─────────────────── Save current settings for next run ───────────────────
> "%SETTINGS_FILE%" (
    echo channel=%channel%
    echo model=%MODEL%
    echo diameter=%DIAMETER%
    echo diam_mode=%DIAM_MODE%
    echo lwr_filter=%LWR_FILTER%
    echo lwr_cutoff=%LWR_CUTOFF%
)

:: ─────────────────── Run Cellpose ───────────────────
echo.
echo Running Cellpose segmentation...
echo --------------------------------------- >> "%LOG_FILE%" 2>&1
echo Starting Cellpose with channel=%channel% model=%MODEL% >> "%LOG_FILE%" 2>&1

cellpose ^
    --dir "%INPUT_DIR%" ^
    --pretrained_model %MODEL% ^
    --chan %channel% ^
    --diameter %DIAMETER% ^
    --savedir "%MASK_DIR%" ^
    --save_png ^
    --verbose >> "%LOG_FILE%" 2>&1

if errorlevel 1 (
    echo [ERROR] Cellpose exited with a non-zero code. Check the log.
    goto wrap_up
)

:: If Cellpose put *_seg.npy files back in INPUT_DIR, move them:
move "%INPUT_DIR%\*_seg.npy" "%MASK_DIR%\" >nul 2>&1

:: ─────────────────── Run Python post-processing ───────────────────
echo.
echo Running Python quantification...
echo --------------------------------------- >> "%LOG_FILE%" 2>&1

"%PYTHON_EXE%" "%QUANT_SCRIPT%" ^
    "%INPUT_DIR%" "%OUTPUT_DIR%" ^
    %LWR_FILTER% %LWR_CUTOFF% %channel% %MODEL% %DIAM_MODE% %DIAMETER% ^
    >> "%LOG_FILE%" 2>&1

if errorlevel 1 (
    echo [ERROR] Python script returned a non-zero exit code. Check the log.
)

:wrap_up
echo. >> "%LOG_FILE%" 2>&1
echo ===== Completed at %DATE% %TIME% ===== >> "%LOG_FILE%" 2>&1
echo Done.
endlocal
pause
