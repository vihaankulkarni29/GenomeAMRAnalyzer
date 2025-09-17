# Hardcore Test Suite

This suite validates URL-based genome discovery and download reliability, correctness, and performance.

## What’s covered
- URL parsing for multiple NCBI formats
- E-utilities integration (mocked) and metadata filtering
- Downloader concurrency, retries, and FASTA validation
- Full workflow orchestration and artifact generation
- Performance characteristics (parallel speed)

## How to run

Windows PowerShell (from project root):

```powershell
python -m pip install -r requirements.txt
pytest -q
```

## Notes
- All network calls are mocked using `aioresponses` so tests are fast and deterministic.
- Performance test uses synthetic payloads to measure concurrency behavior.
- End-to-end test writes outputs to a temporary directory.
