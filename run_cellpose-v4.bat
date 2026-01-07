@echo off
setlocal enabledelayedexpansion

:: ————————————— Settings —————————————
:: (1) Where your raw images live:
set INPUT_DIR=C:\Users\user\Desktop\input_cellpose_v4

:: (2) Where output folder and masks will get written:
set OUTPUT_DIR=%INPUT_DIR%\output
set MASK_DIR=%OUTPUT_DIR%\masks

:: (3) Path to your Python .py quantification script:
set PYTHON_EXE=python
set QUANT_SCRIPT=C:\Users\user\Desktop\scripts_and_notes_cellpose\cellpose_quantify_v4.py

:: (4) Location of your custom Cellpose models (for menu building):
set MODEL_DIR=C:\Users\user\.cellpose\models

:: (5) Path to your log file (captures everything):
set LOG_FILE=%OUTPUT_DIR%\cellpose_run_log.txt

:: (6) Path to save previous settings (in same folder as Python script):
for %%F in ("%QUANT_SCRIPT%") do set SCRIPT_DIR=%%~dpF
set SETTINGS_FILE=%SCRIPT_DIR%cellpose_previous_settings.txt

:: Make sure the output and masks directories exist:
if not exist "%OUTPUT_DIR%" (
    mkdir "%OUTPUT_DIR%"
    echo Created directory "%OUTPUT_DIR%" >> "%LOG_FILE%" 2>&1
)
if not exist "%MASK_DIR%" (
    mkdir "%MASK_DIR%"
    echo Created directory "%MASK_DIR%" >> "%LOG_FILE%" 2>&1
)

:: Make sure the log file exists (or gets created):
if not exist "%LOG_FILE%" (
    > "%LOG_FILE%" echo ===== Cellpose run log started: %DATE% %TIME% =====
)

:: Build menu of models (defaults + any files found in MODEL_DIR)
set "MODEL_COUNT=0"
call :AddModel cyto
call :AddModel nuclei

if exist "%MODEL_DIR%" (
    if exist "%MODEL_DIR%\*.npy" (
        for %%F in ("%MODEL_DIR%\*.npy") do call :AddModel "%%~fF" "%%~nF"
    )
    if exist "%MODEL_DIR%\*.pt" (
        for %%F in ("%MODEL_DIR%\*.pt") do call :AddModel "%%~fF" "%%~nF"
    )
    if exist "%MODEL_DIR%\*.pth" (
        for %%F in ("%MODEL_DIR%\*.pth") do call :AddModel "%%~fF" "%%~nF"
    )
) else (
    echo [WARN] Model directory "%MODEL_DIR%" not found. Only default models will be listed.
)

if "!MODEL_COUNT!"=="0" (
    echo [WARN] No models detected. Adding cyto as fallback.
    call :AddModel cyto
)

:: ————————————— Interactive Prompts —————————————

:: (A) Which channel to segment on? (with option to reuse previous settings)
echo.
echo Which channel do you want to use for segmentation?
echo    [0] Grayscale RGB
echo    [1] Red
echo    [2] Green
echo    [3] Blue
echo     or
echo    [4] Re-run analysis with all the same settings used in the previous run
set /p channel="Enter a channel number (0-3), or enter (4) to re-run analysis with the same settings: "

:: Check if user wants to reuse previous settings
if "%channel%"=="4" (
    if exist "%SETTINGS_FILE%" (
        echo.
        echo Loading previous settings...
        echo ===============================================
        
        :: Read previous settings from file
        for /f "usebackq tokens=1,2 delims==" %%a in ("%SETTINGS_FILE%") do (
            if "%%a"=="channel" set channel=%%b
            if "%%a"=="model" set MODEL=%%b
            if "%%a"=="diameter" set DIAMETER=%%b
            if "%%a"=="diam_mode" set DIAM_MODE=%%b
            if "%%a"=="lwr_filter" set LWR_FILTER=%%b
            if "%%a"=="lwr_cutoff" set LWR_CUTOFF=%%b
        )
        
        :: Display loaded settings
        echo Previous settings loaded:
        echo    Channel: %channel%
        echo    Model: %MODEL%
        echo    Diameter mode: %DIAM_MODE%
        echo    Diameter: %DIAMETER%
        echo    L/W filter: %LWR_FILTER%
        echo    L/W cutoff: %LWR_CUTOFF%
        echo ===============================================
        echo.
        
        :: Log the reuse of previous settings
        echo. >> "%LOG_FILE%" 2>&1
        echo =============================================== >> "%LOG_FILE%" 2>&1
        echo REUSING PREVIOUS SETTINGS: >> "%LOG_FILE%" 2>&1
        echo    Channel: %channel% >> "%LOG_FILE%" 2>&1
        echo    Model: %MODEL% >> "%LOG_FILE%" 2>&1
        echo    Diameter mode: %DIAM_MODE% >> "%LOG_FILE%" 2>&1
        echo    Diameter: %DIAMETER% >> "%LOG_FILE%" 2>&1
        echo    L/W filter: %LWR_FILTER% >> "%LOG_FILE%" 2>&1
        echo    L/W cutoff: %LWR_CUTOFF% >> "%LOG_FILE%" 2>&1
        echo =============================================== >> "%LOG_FILE%" 2>&1
        
        goto run_cellpose
    ) else (
        echo.
        echo No previous settings found. Please enter settings manually.
        echo.
        goto get_channel
    )
) else (
    goto continue_prompts
)

:get_channel
echo Which channel do you want to use for segmentation?
echo    [0] Grayscale RGB
echo    [1] Red
echo    [2] Green
echo    [3] Blue
set /p channel="Enter channel number (0-3): "

:continue_prompts
:: (B) Which Cellpose model?
echo.
echo Available Cellpose models (including any custom models in "%MODEL_DIR%"):
for /L %%I in (1,1,!MODEL_COUNT!) do (
    echo    [%%I] !MODEL_LABEL_%%I!
)
:choose_model
set /p model_choice="Enter a model number (1-!MODEL_COUNT!): "
if not defined model_choice (
    echo Please enter a number between 1 and !MODEL_COUNT!.
    goto choose_model
)
for /f "delims=0123456789" %%a in ("%model_choice%") do (
    echo Please enter digits only.
    goto choose_model
)
if %model_choice% lss 1 (
    echo Invalid choice—please type a number between 1 and !MODEL_COUNT!.
    goto choose_model
)
if %model_choice% gtr !MODEL_COUNT! (
    echo Invalid choice—please type a number between 1 and !MODEL_COUNT!.
    goto choose_model
)
set MODEL=!MODEL_VALUE_%model_choice%!
set MODEL_LABEL=!MODEL_LABEL_%model_choice%!
echo Selected model: !MODEL_LABEL!

:: (C) Auto vs. manual diameter  (changed to 1/2)
:diam_choice
echo.
echo Do you want Cellpose to auto-estimate diameter, or specify it yourself?
echo    [1] Auto (default auto-diameter)
echo    [2] Manual (you will be asked to type in a pixel value)
set /p diam_choice="Enter 1 or 2: "
if "%diam_choice%"=="1" (
    set DIAMETER=0
    set DIAM_MODE=auto
) else if "%diam_choice%"=="2" (
    set DIAM_MODE=manual
    set /p DIAMETER="Enter expected diameter (in pixels): "
) else (
    echo Invalid choice—please type 1 or 2.
    goto diam_choice
)

:: (D) Length/Width ratio filter
:lwr_choice
echo.
echo Should data generation be limited to objects that exceed a particular length to width ratio?
echo    [1] Yes
echo    [2] No
set /p lwr_choice="Enter 1 or 2: "
if "%lwr_choice%"=="1" (
    set LWR_FILTER=1
    set /p LWR_CUTOFF="Only objects with length/width exceeding x will be included in the data output. Enter a value for x (3 is typical for adult cardiomyocytes; press Enter for 3): "
    if "%LWR_CUTOFF%"=="" set LWR_CUTOFF=3
) else if "%lwr_choice%"=="2" (
    set LWR_FILTER=0
    set LWR_CUTOFF=0
) else (
    echo Invalid choice—please type 1 or 2.
    goto lwr_choice
)

:: ————————————— Save Current Settings —————————————
echo Saving current settings...
> "%SETTINGS_FILE%" echo channel=%channel%
>> "%SETTINGS_FILE%" echo model=%MODEL%
>> "%SETTINGS_FILE%" echo diameter=%DIAMETER%
>> "%SETTINGS_FILE%" echo diam_mode=%DIAM_MODE%
>> "%SETTINGS_FILE%" echo lwr_filter=%LWR_FILTER%
>> "%SETTINGS_FILE%" echo lwr_cutoff=%LWR_CUTOFF%

:run_cellpose
:: ————————————— Run Cellpose with redirection —————————————
echo Running segmentation...
echo. >> "%LOG_FILE%" 2>&1
echo Running Cellpose on %INPUT_DIR% with: >> "%LOG_FILE%" 2>&1
echo    channel = %channel%           >> "%LOG_FILE%" 2>&1
echo    model   = %MODEL%             >> "%LOG_FILE%" 2>&1
echo    diam_mode= %DIAM_MODE%        >> "%LOG_FILE%" 2>&1
echo    diameter= %DIAMETER%          >> "%LOG_FILE%" 2>&1
echo    L/W filter? %LWR_FILTER% (cutoff=%LWR_CUTOFF%) >> "%LOG_FILE%" 2>&1
echo --------------------------------------- >> "%LOG_FILE%" 2>&1

cellpose --dir "%INPUT_DIR%" --pretrained_model "%MODEL%" --chan %channel% --diameter %DIAMETER% --savedir "%MASK_DIR%" --save_png --verbose >> "%LOG_FILE%" 2>&1
if errorlevel 1 (
    echo [ERROR] Cellpose returned a non-zero exit code. Check the log. >> "%LOG_FILE%" 2>&1
    goto end
)

move "%INPUT_DIR%\*_seg.npy" "%MASK_DIR%\"

:: ————————————— Run Python quantification (append to same log) —————————————
echo Running quantification script...
echo.                                   >> "%LOG_FILE%" 2>&1
echo Running quantification script...    >> "%LOG_FILE%" 2>&1
echo --------------------------------------- >> "%LOG_FILE%" 2>&1

:: Pass run parameters so they can be logged in the CSV.
:: Args: image_dir output_dir LWR_FILTER LWR_CUTOFF channel MODEL DIAM_MODE DIAMETER
%PYTHON_EXE% "%QUANT_SCRIPT%" "%INPUT_DIR%" "%OUTPUT_DIR%" %LWR_FILTER% %LWR_CUTOFF% %channel% %MODEL% %DIAM_MODE% %DIAMETER% >> "%LOG_FILE%" 2>&1
if errorlevel 1 (
    echo [ERROR] Python quantification returned a non-zero exit code. Check the log. >> "%LOG_FILE%" 2>&1
    goto end
)

:: ————————————— Done —————————————
:end
echo.                                   >> "%LOG_FILE%" 2>&1
echo ===== Completed at %DATE% %TIME% ===== >> "%LOG_FILE%" 2>&1
echo Finished
endlocal
pause

goto :EOF

:AddModel
set "NEW_MODEL_VALUE=%~1"
set "NEW_MODEL_LABEL=%~2"
if "!NEW_MODEL_VALUE!"=="" goto :EOF
if "!NEW_MODEL_LABEL!"=="" set "NEW_MODEL_LABEL=%~1"

for /L %%I in (1,1,!MODEL_COUNT!) do (
    if /I "!MODEL_VALUE_%%I!"=="!NEW_MODEL_VALUE!" goto :EOF
)
set /a MODEL_COUNT+=1
set "MODEL_VALUE_!MODEL_COUNT!=!NEW_MODEL_VALUE!"
set "MODEL_LABEL_!MODEL_COUNT!=!NEW_MODEL_LABEL!"
goto :EOF
