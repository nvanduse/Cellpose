@echo off
setlocal EnableExtensions EnableDelayedExpansion

:: ================================================
:: Cellpose Runner (Cells + Nuclei) -> Cytoplasm metrics
:: - Segments CELLS (model=cyto*/nuclei selectable)
:: - Segments NUCLEI (usually model=nuclei)
:: - Subtracts nuclei from cells; measures intensities on cytoplasm only
:: ================================================

:: --------- USER SETTINGS ---------
set "INPUT_DIR=C:\Users\user\Desktop\input_cellpose_Nppa"            
set "OUTPUT_DIR=%INPUT_DIR%\output"

:: Python exe and script location
set "PYTHON_EXE=python"
set "QUANT_SCRIPT=C:\Users\user\Desktop\scripts_and_notes_cellpose\cellpose_quantify_Nppa_v1.py"

:: Defaults (you can change at prompts)
set "CELL_MODEL=cyto3"
set "NUC_MODEL=nuclei"

echo.
echo ==== CELL SEGMENTATION PARAMETERS ====
set /p CELL_CH="Cell channel (R/G/B or 0/1/2; blank if grayscale) [G]: "
if "%CELL_CH%"=="" set "CELL_CH=G"
set /p CELL_DIAM="Cell diameter in pixels or 'auto' [auto]: "
if "%CELL_DIAM%"=="" set "CELL_DIAM=auto"
set /p CELL_MODEL_IN="Cell model (cyto, cyto2, cyto3, nuclei) [cyto3]: "
if not "%CELL_MODEL_IN%"=="" set "CELL_MODEL=%CELL_MODEL_IN%"

echo.
echo ==== NUCLEI SEGMENTATION PARAMETERS ====
set /p NUC_CH="Nuclei channel (R/G/B or 0/1/2; blank if grayscale) [B]: "
if "%NUC_CH%"=="" set "NUC_CH=B"
set /p NUC_DIAM="Nuclei diameter in pixels or 'auto' [auto]: "
if "%NUC_DIAM%"=="" set "NUC_DIAM=auto"
set /p NUC_MODEL_IN="Nuclei model (nuclei, cyto, cyto2, cyto3) [nuclei]: "
if not "%NUC_MODEL_IN%"=="" set "NUC_MODEL=%NUC_MODEL_IN%"

echo.
echo Optional cell L/W filter (keeps cells with L/W > cutoff)
set /p USE_LWR_FILTER="Enable L/W filter? (1=yes, 0=no) [0]: "
if "%USE_LWR_FILTER%"=="" set "USE_LWR_FILTER=0"
if "%USE_LWR_FILTER%"=="1" (
  set /p LWR_CUTOFF="Enter cutoff (e.g., 1.3): "
) else (
  set "LWR_CUTOFF=0.0"
)

if not exist "%OUTPUT_DIR%" mkdir "%OUTPUT_DIR%"
set "LOG_FILE=%OUTPUT_DIR%\run_log.txt"

echo.                                   >> "%LOG_FILE%"
echo ===== Started %DATE% %TIME% =====   >> "%LOG_FILE%"
echo INPUT_DIR=%INPUT_DIR%               >> "%LOG_FILE%"
echo OUTPUT_DIR=%OUTPUT_DIR%             >> "%LOG_FILE%"
echo CELL: ch=%CELL_CH% diam=%CELL_DIAM% model=%CELL_MODEL% >> "%LOG_FILE%"
echo NUC : ch=%NUC_CH%  diam=%NUC_DIAM%  model=%NUC_MODEL%  >> "%LOG_FILE%"
echo LWR : use=%USE_LWR_FILTER% cutoff=%LWR_CUTOFF%         >> "%LOG_FILE%"
echo.                                   >> "%LOG_FILE%"

"%PYTHON_EXE%" "%QUANT_SCRIPT%" "%INPUT_DIR%" "%OUTPUT_DIR%" ^
  --cell_channel "%CELL_CH%" --cell_diameter "%CELL_DIAM%" --cell_model "%CELL_MODEL%" ^
  --nuc_channel "%NUC_CH%" --nuc_diameter "%NUC_DIAM%" --nuc_model "%NUC_MODEL%" ^
  --lwr_filter %USE_LWR_FILTER% --lwr_cutoff %LWR_CUTOFF% ^
  --save_overlays 1 >> "%LOG_FILE%" 2>&1

if errorlevel 1 (
  echo [ERROR] Quantification script returned an error. See "%LOG_FILE%"
) else (
  echo Completed. See "%OUTPUT_DIR%"
)

echo ===== Finished %DATE% %TIME% =====  >> "%LOG_FILE%"
endlocal
pause
