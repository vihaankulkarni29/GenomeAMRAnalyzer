# Output Directory

This directory will contain all pipeline results after execution.

## Expected Files:
- `analysis_report.html` - Main analysis report
- `resistance_genes.csv` - Detected resistance genes
- `analysis_summary.json` - Machine-readable summary
- `logs/` - Processing logs
- `proteins/` - Extracted protein sequences

## Note:
This directory is automatically created and populated by the Docker pipeline.
Results are preserved between runs unless manually deleted.