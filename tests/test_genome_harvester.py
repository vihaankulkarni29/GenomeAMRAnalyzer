"""
Test Suite for Genome Harvester
-------------------------------
Comprehensive testing for validation, downloading, and error handling.
"""

import pytest
import tempfile
import asyncio
import sys
import os
from pathlib import Path
from unittest.mock import Mock, patch, AsyncMock
import aiohttp

# Add src directory to path
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), 'src', 'core'))

from genome_harvester import GenomeHarvester, AccessionStatus


class TestGenomeHarvester:
    
    @pytest.fixture
    def temp_dir(self):
        """Create temporary directory for testing"""
        with tempfile.TemporaryDirectory() as tmp_dir:
            yield tmp_dir
    
    @pytest.fixture
    def harvester(self, temp_dir):
        """Create test harvester instance"""
        return GenomeHarvester(
            output_dir=temp_dir,
            email="test@example.com",
            max_concurrent=2,
            retry_attempts=2
        )
    
    @pytest.fixture
    def valid_accession_file(self, temp_dir):
        """Create test accession file with valid entries"""
        file_path = Path(temp_dir) / "test_accessions.txt"
        with open(file_path, 'w') as f:
            f.write("GCF_000005825.2\n")
            f.write("GCA_000001405.29\n")
            f.write("NZ_CP012345.1\n")
            f.write("NC_000913.3\n")
        return str(file_path)
    
    @pytest.fixture
    def mixed_accession_file(self, temp_dir):
        """Create test accession file with valid and invalid entries"""
        file_path = Path(temp_dir) / "mixed_accessions.txt"
        with open(file_path, 'w') as f:
            f.write("GCF_000005825.2\n")  # Valid
            f.write("INVALID_123\n")       # Invalid format
            f.write("GCA_000001405.29\n")  # Valid
            f.write("GCF_000005825.2\n")   # Duplicate
            f.write("\n")                  # Empty line
            f.write("NC_000913.3\n")       # Valid
        return str(file_path)
    
    def test_validation_valid_accessions(self, harvester, valid_accession_file):
        """Test validation with all valid accessions"""
        valid_accessions = harvester.validate_accession_list(valid_accession_file)
        
        assert len(valid_accessions) == 4
        assert "GCF_000005825.2" in valid_accessions
        assert "GCA_000001405.29" in valid_accessions
        assert "NZ_CP012345.1" in valid_accessions
        assert "NC_000913.3" in valid_accessions
        
        # Check status tracking
        for acc in valid_accessions:
            assert harvester.status_log[acc].is_valid
            assert not harvester.status_log[acc].is_duplicate
    
    def test_validation_mixed_accessions(self, harvester, mixed_accession_file):
        """Test validation with mixed valid/invalid accessions"""
        valid_accessions = harvester.validate_accession_list(mixed_accession_file)
        
        # Should have 3 valid unique accessions
        assert len(valid_accessions) == 3
        assert "GCF_000005825.2" in valid_accessions
        assert "GCA_000001405.29" in valid_accessions
        assert "NC_000913.3" in valid_accessions
        
        # Check invalid accession
        assert "INVALID_123" in harvester.status_log
        assert not harvester.status_log["INVALID_123"].is_valid
        
        # Check duplicate detection
        duplicate_entries = [s for s in harvester.status_log.values() if s.is_duplicate]
        assert len(duplicate_entries) >= 1
    
    def test_validation_missing_file(self, harvester):
        """Test validation with missing file"""
        with pytest.raises(FileNotFoundError):
            harvester.validate_accession_list("nonexistent_file.txt")
    
    def test_url_construction(self, harvester):
        """Test URL construction for different accession types"""
        
        # Assembly accessions
        gcf_url = harvester._construct_download_url("GCF_000005825.2")
        assert "ftp.ncbi.nlm.nih.gov" in gcf_url
        assert "GCF_000005825_2" in gcf_url
        
        gca_url = harvester._construct_download_url("GCA_000001405.29")
        assert "ftp.ncbi.nlm.nih.gov" in gca_url
        assert "GCA_000001405_29" in gca_url
        
        # Individual sequences
        nc_url = harvester._construct_download_url("NC_000913.3")
        assert "eutils.ncbi.nlm.nih.gov" in nc_url
        assert "NC_000913.3" in nc_url
        
        nz_url = harvester._construct_download_url("NZ_CP012345.1")
        assert "eutils.ncbi.nlm.nih.gov" in nz_url
        assert "NZ_CP012345.1" in nz_url
    
    @pytest.mark.asyncio
    async def test_download_success(self, harvester):
        """Test successful download scenario"""
        
        # Mock successful response
        mock_response = AsyncMock()
        mock_response.status = 200
        mock_response.read.return_value = b">test\nATCG\n"
        
        mock_session = AsyncMock()
        mock_session.get.return_value.__aenter__.return_value = mock_response
        
        # Initialize status
        harvester.status_log["TEST_001"] = AccessionStatus(accession="TEST_001", is_valid=True)
        
        # Test download
        result = await harvester._download_single_genome(mock_session, "TEST_001")
        
        assert result is True
        assert harvester.status_log["TEST_001"].download_success
        assert harvester.status_log["TEST_001"].file_path is not None
    
    @pytest.mark.asyncio
    async def test_download_failure(self, harvester):
        """Test download failure scenarios"""
        
        # Mock failed response
        mock_response = AsyncMock()
        mock_response.status = 404
        mock_response.reason = "Not Found"
        
        mock_session = AsyncMock()
        mock_session.get.return_value.__aenter__.return_value = mock_response
        
        # Initialize status
        harvester.status_log["TEST_FAIL"] = AccessionStatus(accession="TEST_FAIL", is_valid=True)
        
        # Test download
        result = await harvester._download_single_genome(mock_session, "TEST_FAIL")
        
        assert result is False
        assert not harvester.status_log["TEST_FAIL"].download_success
        assert harvester.status_log["TEST_FAIL"].error_message is not None
    
    def test_status_report_generation(self, harvester, temp_dir):
        """Test status report generation"""
        
        # Add some test statuses
        harvester.status_log["SUCCESS_001"] = AccessionStatus(
            accession="SUCCESS_001",
            is_valid=True,
            download_success=True,
            file_size_bytes=1024,
            download_time_seconds=1.5
        )
        
        harvester.status_log["FAILED_001"] = AccessionStatus(
            accession="FAILED_001",
            is_valid=True,
            download_attempted=True,
            download_success=False,
            error_message="Network error"
        )
        
        harvester.status_log["INVALID_001"] = AccessionStatus(
            accession="INVALID_001",
            is_valid=False,
            error_message="Invalid format"
        )
        
        # Generate report
        report_file = Path(temp_dir) / "test_report.txt"
        report = harvester.generate_status_report(str(report_file))
        
        assert report_file.exists()
        assert "SUCCESS_001" in report
        assert "FAILED_001" in report
        assert "INVALID_001" in report
        assert "Network error" in report
        assert "Invalid format" in report


def test_integration_workflow(temp_dir):
    """Integration test for complete workflow"""
    
    # Create test accession file
    accession_file = Path(temp_dir) / "integration_test.txt"
    with open(accession_file, 'w') as f:
        f.write("GCF_000005825.2\n")
        f.write("INVALID_FORMAT\n")
        f.write("GCA_000001405.29\n")
    
    # Initialize harvester
    harvester = GenomeHarvester(
        output_dir=temp_dir,
        email="test@example.com",
        max_concurrent=1,
        retry_attempts=1
    )
    
    # Test validation
    valid_accessions = harvester.validate_accession_list(str(accession_file))
    assert len(valid_accessions) == 2
    
    # Check status tracking
    assert len(harvester.status_log) == 3  # 2 valid + 1 invalid
    invalid_count = sum(1 for s in harvester.status_log.values() if not s.is_valid)
    assert invalid_count == 1
    
    # Generate report
    report = harvester.generate_status_report()
    assert "INVALID_FORMAT" in report
    assert "Invalid accession format" in report


if __name__ == "__main__":
    pytest.main([__file__, "-v"])