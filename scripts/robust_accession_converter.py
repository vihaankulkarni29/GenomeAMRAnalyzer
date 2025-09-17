#!/usr/bin/env python3
"""
Robust Accession Converter for AMR Genomics Pipeline

Converts any NCBI accession (NZ_*, NC_*, CP_*, etc.) to assembly accessions (GCF_*, GCA_*)
with comprehensive error handling, deduplication, and validation for research-grade reliability.

Author: AMR Pipeline Engineering Team
Usage: python robust_accession_converter.py input_file.txt [output_file.txt]
"""

import sys
import re
import time
import logging
import argparse
from pathlib import Path
from urllib.request import urlopen, Request
from urllib.parse import urlencode
from urllib.error import HTTPError, URLError
import xml.etree.ElementTree as ET
from typing import List, Dict, Set, Tuple, Optional
import json

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/accession_conversion.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

class AccessionConverter:
    """Production-grade accession converter with comprehensive error handling."""
    
    def __init__(self, email: str = "amr-pipeline@example.com", api_key: str = "", rate_limit: float = 0.34):
        """
        Initialize converter with NCBI E-utilities configuration.
        
        Args:
            email: Required NCBI email for E-utilities
            api_key: Optional NCBI API key for higher rate limits
            rate_limit: Seconds between requests (0.34 = ~3 requests/second)
        """
        self.email = email
        self.api_key = api_key
        self.rate_limit = rate_limit
        self.conversion_cache = {}
        self.failed_accessions = set()
        self.stats = {
            'total': 0,
            'converted': 0,
            'already_assembly': 0,
            'failed': 0,
            'duplicates_removed': 0
        }
        
    def validate_accession(self, accession: str) -> bool:
        """Validate accession format."""
        if not accession:
            return False
        
        # Common NCBI accession patterns
        patterns = [
            r'^GCF_\d+\.\d+$',  # Assembly RefSeq
            r'^GCA_\d+\.\d+$',  # Assembly GenBank
            r'^NZ_[A-Z]{2}\d+\.\d+$',  # RefSeq contig/scaffold
            r'^NC_\d+\.\d+$',  # RefSeq chromosome
            r'^CP\d+\.\d+$',   # Complete genome
            r'^[A-Z]{1,2}\d+\.\d+$',  # GenBank
        ]
        
        return any(re.match(pattern, accession) for pattern in patterns)
    
    def is_assembly_accession(self, accession: str) -> bool:
        """Check if accession is already an assembly accession."""
        return accession.startswith(('GCF_', 'GCA_'))
    
    def make_eutils_request(self, url: str, max_retries: int = 3) -> str:
        """Make E-utilities request with retry logic and proper headers."""
        for attempt in range(max_retries):
            try:
                headers = {
                    'User-Agent': f'AMR-Pipeline/1.0 ({self.email})',
                    'Accept': 'application/xml'
                }
                
                request = Request(url, headers=headers)
                with urlopen(request, timeout=30) as response:
                    return response.read().decode('utf-8')
                    
            except (HTTPError, URLError) as e:
                if attempt < max_retries - 1:
                    wait_time = 2 ** attempt  # Exponential backoff
                    logger.warning(f"Request failed (attempt {attempt + 1}): {e}. Retrying in {wait_time}s...")
                    time.sleep(wait_time)
                else:
                    logger.error(f"Request failed after {max_retries} attempts: {e}")
                    raise
            except Exception as e:
                logger.error(f"Unexpected error in request: {e}")
                raise
                
    def convert_accession(self, accession: str) -> Tuple[str, str]:
        """
        Convert single accession to assembly format.
        
        Returns:
            Tuple of (converted_accession, status_message)
        """
        accession = accession.strip()
        
        if not self.validate_accession(accession):
            return accession, "INVALID_FORMAT"
            
        if self.is_assembly_accession(accession):
            self.stats['already_assembly'] += 1
            return accession, "ALREADY_ASSEMBLY"
            
        if accession in self.conversion_cache:
            return self.conversion_cache[accession], "CACHED"
            
        if accession in self.failed_accessions:
            return accession, "PREVIOUSLY_FAILED"
            
        try:
            # Step 1: Search for assembly using the accession
            search_params = {
                'db': 'assembly',
                'term': f'"{accession}"[RefSeq Accession] OR "{accession}"[GenBank Accession]',
                'retmode': 'xml',
                'email': self.email
            }
            
            if self.api_key:
                search_params['api_key'] = self.api_key
                
            search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?{urlencode(search_params)}"
            
            logger.debug(f"Searching for assembly: {accession}")
            search_xml = self.make_eutils_request(search_url)
            search_root = ET.fromstring(search_xml)
            
            # Check for errors
            error_list = search_root.find('.//ErrorList')
            if error_list is not None and len(error_list) > 0:
                error_msg = error_list[0].text if error_list[0].text else "Unknown error"
                logger.warning(f"NCBI search error for {accession}: {error_msg}")
                self.failed_accessions.add(accession)
                return accession, f"NCBI_ERROR: {error_msg}"
            
            id_list = search_root.find('.//IdList')
            if id_list is None or len(id_list) == 0:
                logger.warning(f"No assembly found for {accession}")
                self.failed_accessions.add(accession)
                return accession, "NO_ASSEMBLY_FOUND"
                
            # Step 2: Get assembly summary
            uid = id_list[0].text
            summary_params = {
                'db': 'assembly',
                'id': uid,
                'retmode': 'xml',
                'email': self.email
            }
            
            if self.api_key:
                summary_params['api_key'] = self.api_key
                
            summary_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?{urlencode(summary_params)}"
            
            logger.debug(f"Getting assembly summary for UID: {uid}")
            summary_xml = self.make_eutils_request(summary_url)
            summary_root = ET.fromstring(summary_xml)
            
            # Parse assembly accession
            assembly_acc = self.parse_assembly_accession(summary_root)
            if assembly_acc:
                self.conversion_cache[accession] = assembly_acc
                self.stats['converted'] += 1
                logger.info(f"Converted: {accession} → {assembly_acc}")
                return assembly_acc, "CONVERTED"
            else:
                logger.warning(f"Could not extract assembly accession for {accession}")
                self.failed_accessions.add(accession)
                return accession, "PARSE_ERROR"
                
        except Exception as e:
            logger.error(f"Error converting {accession}: {e}")
            self.failed_accessions.add(accession)
            return accession, f"EXCEPTION: {str(e)}"
        finally:
            # Rate limiting
            time.sleep(self.rate_limit)
            
    def parse_assembly_accession(self, summary_root: ET.Element) -> Optional[str]:
        """Extract assembly accession from eSummary XML."""
        # Try multiple XML structures
        for doc_elem in summary_root.findall('.//DocumentSummary') + summary_root.findall('.//DocSum'):
            # Direct child
            acc_elem = doc_elem.find('.//AssemblyAccession')
            if acc_elem is not None and acc_elem.text:
                return acc_elem.text.strip()
                
            # Item with Name attribute
            for item in doc_elem.findall('.//Item[@Name="AssemblyAccession"]'):
                if item.text:
                    return item.text.strip()
                    
        return None
        
    def deduplicate_accessions(self, accessions: List[str]) -> List[str]:
        """Remove duplicates while preserving order."""
        seen = set()
        result = []
        for acc in accessions:
            if acc and acc not in seen:
                seen.add(acc)
                result.append(acc)
            elif acc in seen:
                self.stats['duplicates_removed'] += 1
                
        return result
        
    def convert_file(self, input_file: Path, output_file: Path = None) -> Dict:
        """
        Convert accessions from input file.
        
        Returns:
            Dictionary with conversion results and statistics
        """
        if output_file is None:
            output_file = input_file
            
        logger.info(f"Converting accessions from {input_file}")
        
        # Read input file
        try:
            with open(input_file, 'r', encoding='utf-8') as f:
                lines = [line.strip() for line in f if line.strip()]
        except Exception as e:
            logger.error(f"Failed to read input file {input_file}: {e}")
            return {'success': False, 'error': str(e)}
            
        if not lines:
            logger.warning("Input file is empty")
            return {'success': False, 'error': 'Input file is empty'}
            
        # Deduplicate
        original_count = len(lines)
        lines = self.deduplicate_accessions(lines)
        self.stats['total'] = len(lines)
        
        logger.info(f"Processing {len(lines)} unique accessions (removed {original_count - len(lines)} duplicates)")
        
        # Convert accessions
        converted_lines = []
        conversion_log = []
        
        for i, accession in enumerate(lines, 1):
            logger.info(f"Processing {i}/{len(lines)}: {accession}")
            
            converted_acc, status = self.convert_accession(accession)
            converted_lines.append(converted_acc)
            
            conversion_log.append({
                'original': accession,
                'converted': converted_acc,
                'status': status
            })
            
            # Progress indicator
            if i % 10 == 0:
                logger.info(f"Progress: {i}/{len(lines)} ({i/len(lines)*100:.1f}%)")
                
        # Write output file
        try:
            with open(output_file, 'w', encoding='utf-8') as f:
                for acc in converted_lines:
                    f.write(f"{acc}\n")
        except Exception as e:
            logger.error(f"Failed to write output file {output_file}: {e}")
            return {'success': False, 'error': str(e)}
            
        # Save detailed conversion log
        log_file = output_file.parent / f"{output_file.stem}_conversion_log.json"
        try:
            with open(log_file, 'w', encoding='utf-8') as f:
                json.dump({
                    'conversion_log': conversion_log,
                    'statistics': self.stats,
                    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
                }, f, indent=2)
        except Exception as e:
            logger.warning(f"Failed to write conversion log: {e}")
            
        self.stats['failed'] = len(self.failed_accessions)
        
        logger.info("Conversion completed!")
        logger.info(f"Statistics: {self.stats}")
        
        return {
            'success': True,
            'statistics': self.stats,
            'conversion_log': conversion_log,
            'log_file': str(log_file)
        }


def main():
    """Command-line interface for accession converter."""
    parser = argparse.ArgumentParser(
        description="Convert NCBI accessions to assembly format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python robust_accession_converter.py accession_list.txt
  python robust_accession_converter.py input.txt output.txt --email user@example.com
  python robust_accession_converter.py input.txt --api-key YOUR_KEY --rate-limit 0.1
        """
    )
    
    parser.add_argument('input_file', help='Input file with accessions (one per line)')
    parser.add_argument('output_file', nargs='?', help='Output file (default: overwrite input)')
    parser.add_argument('--email', default='amr-pipeline@example.com', 
                       help='NCBI email (required for E-utilities)')
    parser.add_argument('--api-key', default='', help='NCBI API key (optional, for higher rate limits)')
    parser.add_argument('--rate-limit', type=float, default=0.34, 
                       help='Seconds between requests (default: 0.34)')
    parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], 
                       default='INFO', help='Logging level')
    
    args = parser.parse_args()
    
    # Set log level
    logging.getLogger().setLevel(getattr(logging, args.log_level))
    
    # Create logs directory
    Path('logs').mkdir(exist_ok=True)
    
    # Initialize converter
    converter = AccessionConverter(
        email=args.email,
        api_key=args.api_key,
        rate_limit=args.rate_limit
    )
    
    # Convert file
    input_path = Path(args.input_file)
    output_path = Path(args.output_file) if args.output_file else input_path
    
    result = converter.convert_file(input_path, output_path)
    
    if result['success']:
        print(f"\n✅ Conversion completed successfully!")
        print(f"📊 Statistics: {result['statistics']}")
        print(f"📝 Detailed log saved to: {result['log_file']}")
        sys.exit(0)
    else:
        print(f"\n❌ Conversion failed: {result['error']}")
        sys.exit(1)


if __name__ == "__main__":
    main()