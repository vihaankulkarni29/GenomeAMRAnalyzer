@echo off
setlocal

REM Simple user-friendly runner for GenomeAMRAnalyzer on Windows
REM - Uses the production config at config\snakemake_config.yaml
REM - Ensures Python is available and creates logs directory

echo === GenomeAMRAnalyzer ===
echo Running with config: config\snakemake_config.yaml

where python >nul 2>&1
if errorlevel 1 (
  echo ERROR: Python not found in PATH.
  echo Please install Python 3.9+ and re-run.
  exit /b 1
)

REM Ensure logs directory exists
if not exist logs (
  mkdir logs >nul 2>&1
)

REM Optional: remind user to set NCBI email
for /f "tokens=2 delims=:" %%E in ('findstr /b /c:"  email:" config\snakemake_config.yaml') do set "NCBI_EMAIL=%%E"
set NCBI_EMAIL=%NCBI_EMAIL:"=%
if "%NCBI_EMAIL%"=="  \"your.email@institution.edu\"" (
  echo.
  echo WARNING: NCBI email is still the placeholder. Edit config\snakemake_config.yaml (ncbi.email).
)

echo Starting pipeline...
python run_pipeline.py --config config\snakemake_config.yaml 1> logs\user_pipeline.out 2> logs\user_pipeline.err
set EXITCODE=%ERRORLEVEL%

if %EXITCODE% NEQ 0 (
  echo Pipeline failed. See logs\user_pipeline.err for details.
  exit /b %EXITCODE%
)

echo Pipeline completed successfully.
echo Outputs should be in the directories configured under 'directories' in the config file.
exit /b 0
