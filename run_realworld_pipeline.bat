@echo off
setlocal enabledelayedexpansion

REM ============================================
REM Real-World AMR Pipeline Runner
REM ============================================
echo ==========================================
echo AMR GENOMICS RESEARCH PIPELINE  
echo ==========================================
echo.
echo RESEARCH OBJECTIVE: 
echo "How are antimicrobial resistance genomes developing mutations?"
echo.
echo CRITICAL CHECK: RGI (Resistance Gene Identifier) Required
echo ==========================================

REM Check for RGI first - it's essential for research functionality
rgi --version >nul 2>&1
if errorlevel 1 (
    echo.
    echo ❌ CRITICAL: RGI not found - Pipeline will have limited research value!
    echo.
    echo Without RGI:
    echo   ❌ No resistance gene detection
    echo   ❌ No protein extraction  
    echo   ❌ No mutation analysis
    echo   ❌ No research insights
    echo.
    echo RECOMMENDATION: Install RGI first for full functionality
    echo Run: setup_rgi.bat for installation guide
    echo.
    set /p continue="Continue anyway? (y/N): "
    if /i not "!continue!"=="y" (
        echo.
        echo Please install RGI and try again.
        echo See: RGI_INSTALLATION_GUIDE.md
        pause
        exit /b 1
    )
    echo.
    echo ⚠️  WARNING: Continuing with limited functionality...
) else (
    echo ✅ RGI detected - Full research functionality available!
)

echo.
echo Starting Real-World AMR Pipeline...

REM Check if config file exists
if not exist "accession_pipeline.yaml" (
    echo ERROR: accession_pipeline.yaml not found in current directory.
    echo Please ensure you're running from the GenomeAMRAnalyzer root directory.
    pause
    exit /b 1
)

REM Check required files (use fixed names from config)
set "ACC_FILE=accession_list.txt"
set "GENE_FILE=gene_list.txt"

if not exist "%ACC_FILE%" (
    echo ERROR: %ACC_FILE% not found.
    pause
    exit /b 1
)

if not exist "%GENE_FILE%" (
    echo ERROR: %GENE_FILE% not found.
    pause
    exit /b 1
)

echo Found config and input files. Proceeding...

REM ============================================
REM Set PYTHONPATH and Environment
REM ============================================
set "PYTHONPATH=%CD%;%PYTHONPATH%"
echo PYTHONPATH set to: %PYTHONPATH%

REM Set NCBI email from config or prompt
set "NCBI_EMAIL=vihaankulkarni29@gmail.com"
echo Using NCBI email: %NCBI_EMAIL%

REM ============================================
REM Validate inputs before pipeline execution
REM ============================================
echo Validating pipeline inputs...

REM Check: Ensure all accessions are valid, non-empty, and deduplicated before download
python scripts\validate_inputs.py

if errorlevel 1 (
    echo Input validation failed. Please check your files.
    pause
    exit /b 1
)

echo ✅ Input validation completed successfully

REM ============================================
REM Auto-convert accessions using robust converter
REM ============================================
set "ACC_BAK=%ACC_FILE%.backup"

echo Converting accessions to assembly format...
copy /y "%ACC_FILE%" "%ACC_BAK%" >nul

echo Running robust accession converter...
python scripts\robust_accession_converter.py "%ACC_FILE%" --email "%NCBI_EMAIL%" --log-level INFO
if errorlevel 1 (
    echo Conversion had issues but continuing...
    echo Check logs\accession_conversion.log for details
) else (
    echo Conversion completed successfully.
)
    echo     try: >> temp_convert.py
    echo         search_params = urlencode({'db': 'assembly', 'term': f'{acc}[RefSeq Accession] OR {acc}[GenBank Accession]', 'retmode': 'xml'}) >> temp_convert.py
    echo         search_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?{search_params}' >> temp_convert.py
    echo         with urlopen(search_url) as response: >> temp_convert.py
    echo             search_xml = response.read().decode() >> temp_convert.py
    echo         search_root = ET.fromstring(search_xml) >> temp_convert.py
    echo         id_list = search_root.find('.//IdList') >> temp_convert.py
    echo         if id_list is None or len(id_list) == 0: >> temp_convert.py
    echo             print(f'Warning: No assembly found for {acc}', file=sys.stderr) >> temp_convert.py
    echo             return acc >> temp_convert.py
    echo         uid = id_list[0].text >> temp_convert.py
    echo         summary_params = urlencode({'db': 'assembly', 'id': uid, 'retmode': 'xml'}) >> temp_convert.py
    echo         summary_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?{summary_params}' >> temp_convert.py
    echo         with urlopen(summary_url) as response: >> temp_convert.py
    echo             summary_xml = response.read().decode() >> temp_convert.py
    echo         summary_root = ET.fromstring(summary_xml) >> temp_convert.py
    echo         for doc_elem in summary_root.findall('.//DocumentSummary') + summary_root.findall('.//DocSum'): >> temp_convert.py
    echo             acc_elem = doc_elem.find('.//AssemblyAccession') >> temp_convert.py
    echo             if acc_elem is not None and acc_elem.text: >> temp_convert.py
    echo                 return acc_elem.text.strip() >> temp_convert.py
    echo             for item in doc_elem.findall('.//Item[@Name="AssemblyAccession"]'): >> temp_convert.py
    echo                 if item.text: >> temp_convert.py
    echo                     return item.text.strip() >> temp_convert.py
    echo         print(f'Warning: Could not extract assembly accession for {acc}', file=sys.stderr) >> temp_convert.py
    echo         return acc >> temp_convert.py
    echo     except Exception as e: >> temp_convert.py
    echo         print(f'Error converting {acc}: {e}', file=sys.stderr) >> temp_convert.py
    echo         return acc >> temp_convert.py
    echo. >> temp_convert.py
    echo with open(sys.argv[1], 'r') as f: >> temp_convert.py
    echo     lines = f.readlines() >> temp_convert.py
    echo. >> temp_convert.py
    echo converted_lines = [] >> temp_convert.py
    echo for line in lines: >> temp_convert.py
    echo     acc = line.strip() >> temp_convert.py
    echo     if acc: >> temp_convert.py
    echo         converted = convert_accession(acc) >> temp_convert.py
    echo         converted_lines.append(converted + '\n') >> temp_convert.py
    echo         if converted != acc: >> temp_convert.py
    echo             print(f'Converted: {acc} -^> {converted}') >> temp_convert.py
    echo         time.sleep(0.1) >> temp_convert.py
    echo     else: >> temp_convert.py
    echo         converted_lines.append(line) >> temp_convert.py
    echo. >> temp_convert.py
    echo with open(sys.argv[1], 'w') as f: >> temp_convert.py
    echo     f.writelines(converted_lines) >> temp_convert.py
    
    python temp_convert.py "%ACC_FILE%"
    del temp_convert.py
    
    if errorlevel 1 (
        echo Conversion had issues but continuing...
    ) else (
        echo Conversion completed successfully.
    )
) else (
    echo accession_list.txt contains only assembly accessions (GCF_/GCA_). No conversion needed.
)

REM ============================================
REM Run the orchestrator as a module (disable resume by default)
REM ============================================
set "AMR_DISABLE_RESUME=1"
set "ORCH=src.priority3.pipeline.orchestrator"

REM If resume is disabled, remove any stale harvest checkpoint before running
if defined AMR_DISABLE_RESUME (
    if exist "genome_data\harvest_checkpoint.json" del /q "genome_data\harvest_checkpoint.json" >nul 2>&1
)

echo Running pipeline orchestrator...
python -m %ORCH% -c accession_pipeline.yaml
if errorlevel 1 (
    echo Pipeline execution failed. Check logs\pipeline.log for details.
    echo.
    echo Common issues and solutions:
    echo - RGI not installed: conda install -c bioconda rgi
    echo - Missing references: Check reference_dir in config
    echo - Network issues: Check internet connectivity for NCBI access
    echo.
    pause
    goto RESTORE_BACKUP
)

echo Pipeline completed successfully.

REM ============================================
REM Post-Pipeline Validation and Summary
REM ============================================
echo Validating pipeline results...

python scripts\validate_results.py
    # Comprehensive validation of pipeline outputs
    print('\\n=== Pipeline Results Summary ===')
    
    # Check genome downloads
    genome_dir = Path('genome_data')
    if genome_dir.exists():
        genome_files = list(genome_dir.glob('*.fasta')) + list(genome_dir.glob('*.fa')) + list(genome_dir.glob('*.fna'))
        print(f'Downloaded genomes: {len(genome_files)}')
    else:
        print('Downloaded genomes: 0 (directory not found)')
        
    # Check CARD results
    card_dir = Path('card_results')
    if card_dir.exists():
        card_files = list(card_dir.glob('*_rgi.txt'))
        print(f'CARD RGI results: {len(card_files)}')
    else:
        print('CARD RGI results: 0 (directory not found)')
        
    # Check protein extractions
    protein_dir = Path('extracted_proteins')
    if protein_dir.exists():
        protein_files = list(protein_dir.glob('*.faa'))
        print(f'Extracted proteins: {len(protein_files)}')
    else:
        print('Extracted proteins: 0 (directory not found)')
        
    # Check alignments
    align_dir = Path('alignments')
    if align_dir.exists():
        align_files = list(align_dir.glob('*.txt')) + list(align_dir.glob('*.aln'))
        print(f'Alignments: {len(align_files)}')
    else:
        print('Alignments: 0 (directory not found)')
        
    # Check mutations
    subscan_dir = Path('subscan_results')
    if subscan_dir.exists():
        subscan_files = list(subscan_dir.glob('*.json'))
        print(f'Mutation analyses: {len(subscan_files)}')
    else:
        print('Mutation analyses: 0 (directory not found)')
        
    # Check co-occurrence
    cooccur_dir = Path('cooccurrence_results')
    if cooccur_dir.exists():
        cooccur_files = list(cooccur_dir.glob('*.csv'))
        print(f'Co-occurrence analyses: {len(cooccur_files)}')
    else:
        print('Co-occurrence analyses: 0 (directory not found)')
        
    # Check reports
    report_dir = Path('report/html_reports')
    if report_dir.exists():
        report_files = list(report_dir.glob('*.html'))
        print(f'HTML reports: {len(report_files)}')
        if report_files:
            latest_report = max(report_files, key=os.path.getctime)
            print(f'Latest report: {latest_report.name}')
    else:
        print('HTML reports: 0 (directory not found)')
        
    # Check for common issues
    print('\\n=== Issue Detection ===')
    issues_found = False
    
    if not genome_dir.exists() or len(list(genome_dir.glob('*.fasta'))) == 0:
        print('❌ No genomes downloaded - check NCBI access and accession validity')
        issues_found = True
        
    if not card_dir.exists() or len(list(card_dir.glob('*_rgi.txt'))) == 0:
        print('⚠️  No CARD results - RGI may not be installed')
        print('   Solution: conda install -c bioconda rgi')
        
    if not protein_dir.exists() or len(list(protein_dir.glob('*.faa'))) == 0:
        print('⚠️  No proteins extracted - depends on CARD results')
        
    if not issues_found:
        print('✅ No critical issues detected')
        
    print('\\n=== Next Steps ===')

REM ============================================
REM Post-run validation and reporting
REM ============================================

REM Post-run validation: count downloaded genomes (multiple extensions)
set "GENOME_DIR=genome_data"
set "GENOME_COUNT=0"
if exist "%GENOME_DIR%\" (
    for %%E in (fasta fa fna gz) do (
        for /r "%GENOME_DIR%" %%F in (*.%%E) do set /a GENOME_COUNT+=1
    )
)

echo Downloaded genomes detected: %GENOME_COUNT%

REM Warn if no genomes found and offer to clear checkpoints and re-run
if "%GENOME_COUNT%"=="0" (
    echo.
    echo WARNING: No genomes were downloaded. This might be due to:
    echo - Stale checkpoints from previous runs
    echo - Network issues or API rate limiting
    echo - Accession conversion problems
    echo.
    set /p CLEAR_RETRY="Clear all checkpoints and re-run once? (Y/N): "
    if /i "!CLEAR_RETRY!"=="Y" (
        echo Clearing checkpoints and cache directories...
        call :CLEAR_CHECKPOINTS
        echo Re-running pipeline once after clearing checkpoints...
        REM ensure harvest checkpoint is removed
        if exist "genome_data\harvest_checkpoint.json" del /q "genome_data\harvest_checkpoint.json" >nul 2>&1
        python -m %ORCH% -c accession_pipeline.yaml --no-resume
        if errorlevel 1 (
            python -m %ORCH% -c accession_pipeline.yaml
        )
        REM recount
        set "GENOME_COUNT=0"
        if exist "%GENOME_DIR%\" (
            for %%E in (fasta fa fna gz) do (
                for /r "%GENOME_DIR%" %%F in (*.%%E) do set /a GENOME_COUNT+=1
            )
        )
        echo Downloaded genomes detected after re-run: %GENOME_COUNT%
    )
)

REM Auto-open the latest report
echo Checking for generated reports...
set "LATEST_REPORT="
for /f "delims=" %%F in ('dir /b /o:-d "report\html_reports\*.html" 2^>nul') do (
    if not defined LATEST_REPORT set "LATEST_REPORT=report\html_reports\%%F"
)
if not defined LATEST_REPORT (
    for /f "delims=" %%F in ('dir /b /o:-d "report\*.html" 2^>nul') do (
        if not defined LATEST_REPORT set "LATEST_REPORT=report\%%F"
    )
)

if defined LATEST_REPORT (
    echo Opening latest report: %LATEST_REPORT%
    start "" "%LATEST_REPORT%"
) else (
    echo No HTML reports found in report\ directories.
)

:RESTORE_BACKUP
REM Restore original accession file if backup exists
if exist "%ACC_BAK%" (
    echo Restoring original accession file...
    copy /y "%ACC_BAK%" "%ACC_FILE%" >nul
    del /q "%ACC_BAK%" >nul
)

echo.
echo Pipeline run complete. Check logs\pipeline.log for detailed execution information.
pause
exit /b 0

:CLEAR_CHECKPOINTS
for %%D in (checkpoints .checkpoints .cache cache) do (
    if exist "%%D" (
        echo Removing %%D ...
        rmdir /s /q "%%D"
    )
)
if exist "genome_data\harvest_checkpoint.json" del /q "genome_data\harvest_checkpoint.json" >nul 2>&1
exit /b 0