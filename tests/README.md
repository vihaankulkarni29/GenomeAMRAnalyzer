# GenomeAMRAnalyzer Testing Strategy

This document outlines the testing strategy for the GenomeAMRAnalyzer pipeline, ensuring code quality, correctness, and stability.

## Frameworks
- **Test Runner:** `pytest`
- **Mocking Library:** `pytest-mock`
- **Assertion Library:** Native `pytest` assertions

## Test Suite Structure
The `tests/` directory contains unit and integration tests for the core modules of the pipeline.

- **Unit Tests:** Each module in `src/` has a corresponding `test_*.py` file. These tests focus on a single unit of code in isolation, with external dependencies (like subprocess calls or file I/O) mocked out.
- **Integration Tests:** Tests for the main pipeline orchestrator (`test_production_pipeline_orchestrator.py`) verify that the modules interact correctly.

## How to Run Tests
To run the full test suite, navigate to the root of the repository and execute the following command from within the activated Docker or Conda environment:

```bash
pytest -v
```
This command will automatically discover and run all test files.
