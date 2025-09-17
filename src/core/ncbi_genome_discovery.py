"""
NCBI URL Parser and Genome Discovery Engine
===========================================
Professional URL parser that extracts search terms from NCBI URLs and 
discovers genome accessions using E-utilities.
"""

import re
import logging
import asyncio
import aiohttp
import xml.etree.ElementTree as ET
from typing import List, Dict, Optional, Tuple
from urllib.parse import parse_qs, urlparse, unquote_plus
from dataclasses import dataclass
import time


@dataclass
class GenomeMetadata:
    """Store metadata for each discovered genome"""
    accession: str
    title: str = ""
    organism: str = ""
    length: int = 0
    pub_date: str = ""
    molecule_type: str = ""
    topology: str = ""


class NCBIUrlParser:
    """
    Parse NCBI search URLs and extract search terms
    """
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        self.logger = logger or self._setup_logger()
    
    def _setup_logger(self) -> logging.Logger:
        """Setup logging for URL parser"""
        logger = logging.getLogger("NCBIUrlParser")
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(levelname)s %(name)s: %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        return logger
    
    def parse_search_url(self, ncbi_url: str) -> str:
        """
        Extract search term from NCBI URL
        
        Args:
            ncbi_url: NCBI search URL
            
        Returns:
            Cleaned search term for E-utilities
            
        Raises:
            ValueError: If URL format is invalid or unsupported
        """
        self.logger.info(f"Parsing NCBI URL: {ncbi_url}")
        
        try:
            parsed = urlparse(ncbi_url)
            
            # Handle different NCBI URL formats
            if 'nuccore' in parsed.path or 'nucleotide' in parsed.path:
                search_term = self._extract_nuccore_term(ncbi_url, parsed)
            elif 'genome' in parsed.path:
                search_term = self._extract_genome_term(ncbi_url, parsed)
            else:
                raise ValueError(f"Unsupported NCBI URL type: {parsed.path}")
            
            # Clean and validate search term
            cleaned_term = self._clean_search_term(search_term)
            self.logger.info(f"Extracted search term: {cleaned_term}")
            return cleaned_term
            
        except Exception as e:
            self.logger.error(f"Failed to parse NCBI URL: {e}")
            raise ValueError(f"Invalid NCBI URL format: {e}")
    
    def _extract_nuccore_term(self, url: str, parsed) -> str:
        """Extract search term from nuccore/nucleotide URLs"""
        query_params = parse_qs(parsed.query)
        
        # Method 1: Direct term parameter
        if 'term' in query_params:
            return query_params['term'][0]
        
        # Method 2: Extract from URL string (for complex URLs)
        if 'term=' in url:
            term_start = url.find('term=') + 5
            term_end = url.find('&', term_start)
            if term_end == -1:
                term_end = len(url)
            return unquote_plus(url[term_start:term_end])
        
        raise ValueError("No search term found in nuccore URL")
    
    def _extract_genome_term(self, url: str, parsed) -> str:
        """Extract search term from genome URLs"""
        # Genome URLs often have different parameter structures
        query_params = parse_qs(parsed.query)
        
        if 'term' in query_params:
            return query_params['term'][0]
        elif 'q' in query_params:
            return query_params['q'][0]
        
        raise ValueError("No search term found in genome URL")
    
    def _clean_search_term(self, term: str) -> str:
        """Clean and validate search term"""
        # URL decode if needed
        cleaned = unquote_plus(term)
        
        # Remove extra whitespace
        cleaned = ' '.join(cleaned.split())
        
        # Basic validation
        if len(cleaned.strip()) < 3:
            raise ValueError("Search term too short")
        
        return cleaned


class NCBIGenomeDiscovery:
    """
    Discover genome accessions using NCBI E-utilities
    """
    
    def __init__(self, 
                 email: str,
                 api_key: Optional[str] = None,
                 max_results: int = 10,
                 logger: Optional[logging.Logger] = None):
        
        self.email = email
        self.api_key = api_key
        self.max_results = max_results
        self.logger = logger or self._setup_logger()
        
        # E-utilities endpoints
        self.esearch_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        self.esummary_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        self.efetch_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        
        # Results storage
        self.discovered_genomes: Dict[str, GenomeMetadata] = {}
    
    def _setup_logger(self) -> logging.Logger:
        """Setup logging for genome discovery"""
        logger = logging.getLogger("NCBIGenomeDiscovery")
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(levelname)s %(name)s: %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        return logger
    
    async def discover_genomes(self, search_term: str) -> List[GenomeMetadata]:
        """
        Discover genomes using NCBI E-utilities
        
        Args:
            search_term: Search query for NCBI
            
        Returns:
            List of GenomeMetadata objects
        """
        self.logger.info(f"Discovering genomes for: {search_term}")
        
        try:
            # Step 1: Search for genome IDs
            genome_ids = await self._search_genomes(search_term)
            
            if not genome_ids:
                self.logger.warning("No genomes found for search term")
                return []
            
            # Step 2: Get metadata for discovered genomes
            genomes = await self._get_genome_metadata(genome_ids)
            
            # Step 3: Filter and validate results
            filtered_genomes = self._filter_genomes(genomes)
            
            self.logger.info(f"Discovered {len(filtered_genomes)} quality genomes")
            return filtered_genomes
            
        except Exception as e:
            self.logger.error(f"Genome discovery failed: {e}")
            raise
    
    async def _search_genomes(self, search_term: str) -> List[str]:
        """Search NCBI for genome IDs"""
        params = {
            'db': 'nuccore',
            'term': search_term,
            'retmax': str(self.max_results * 2),  # Get extra to allow filtering
            'retmode': 'xml',
            'tool': 'GenomeHarvester',
            'email': self.email
        }
        
        if self.api_key:
            params['api_key'] = self.api_key
        
        async with aiohttp.ClientSession() as session:
            async with session.get(self.esearch_base, params=params) as response:
                if response.status != 200:
                    raise Exception(f"NCBI search failed: HTTP {response.status}")
                
                xml_content = await response.text()
                root = ET.fromstring(xml_content)
                
                # Extract genome IDs
                id_list = []
                for id_elem in root.findall('.//Id'):
                    if id_elem.text:
                        id_list.append(id_elem.text)
                
                self.logger.info(f"Found {len(id_list)} genome IDs")
                return id_list
    
    async def _get_genome_metadata(self, genome_ids: List[str]) -> List[GenomeMetadata]:
        """Get detailed metadata for genome IDs"""
        if not genome_ids:
            return []
        
        # Process in batches to avoid URL length limits
        batch_size = 200
        all_genomes = []
        
        for i in range(0, len(genome_ids), batch_size):
            batch_ids = genome_ids[i:i + batch_size]
            batch_genomes = await self._get_batch_metadata(batch_ids)
            all_genomes.extend(batch_genomes)
        
        return all_genomes
    
    async def _get_batch_metadata(self, genome_ids: List[str]) -> List[GenomeMetadata]:
        """Get metadata for a batch of genome IDs"""
        params = {
            'db': 'nuccore',
            'id': ','.join(genome_ids),
            'retmode': 'xml',
            'tool': 'GenomeHarvester',
            'email': self.email
        }
        
        if self.api_key:
            params['api_key'] = self.api_key
        
        genomes = []
        
        async with aiohttp.ClientSession() as session:
            async with session.get(self.esummary_base, params=params) as response:
                if response.status != 200:
                    self.logger.warning(f"Batch metadata request failed: HTTP {response.status}")
                    return []
                
                xml_content = await response.text()
                root = ET.fromstring(xml_content)
                
                # Parse each genome summary
                for doc_sum in root.findall('.//DocSum'):
                    genome = self._parse_genome_summary(doc_sum)
                    if genome:
                        genomes.append(genome)
                        self.discovered_genomes[genome.accession] = genome
        
        return genomes
    
    def _parse_genome_summary(self, doc_sum) -> Optional[GenomeMetadata]:
        """Parse individual genome summary from XML"""
        try:
            accession = ""
            title = ""
            organism = ""
            length = 0
            pub_date = ""
            molecule_type = ""
            topology = ""
            
            # Extract all relevant fields
            for item in doc_sum.findall('.//Item'):
                name = item.get('Name', '')
                text = item.text or ""
                
                if name == 'AccessionVersion':
                    accession = text
                elif name == 'Title':
                    title = text
                elif name == 'Organism':
                    organism = text
                elif name == 'Length':
                    try:
                        length = int(text)
                    except ValueError:
                        length = 0
                elif name == 'PubDate':
                    pub_date = text
                elif name == 'MolType':
                    molecule_type = text
                elif name == 'Topology':
                    topology = text
            
            if accession:
                return GenomeMetadata(
                    accession=accession,
                    title=title,
                    organism=organism,
                    length=length,
                    pub_date=pub_date,
                    molecule_type=molecule_type,
                    topology=topology
                )
            
        except Exception as e:
            self.logger.warning(f"Failed to parse genome summary: {e}")
        
        return None
    
    def _filter_genomes(self, genomes: List[GenomeMetadata]) -> List[GenomeMetadata]:
        """Filter genomes based on quality criteria"""
        filtered = []
        
        for genome in genomes:
            # Size filtering
            if genome.length < 1_000_000:  # Less than 1MB
                self.logger.debug(f"Filtered {genome.accession}: too small ({genome.length} bp)")
                continue
            
            if genome.length > 20_000_000:  # More than 20MB
                self.logger.debug(f"Filtered {genome.accession}: too large ({genome.length} bp)")
                continue
            
            # Prefer complete genomes
            is_complete = any(term in genome.title.lower() for term in [
                'complete genome', 'complete sequence', 'chromosome'
            ])
            
            if is_complete or len(filtered) < self.max_results // 2:
                # Accept complete genomes or if we need more results
                filtered.append(genome)
            
            # Stop when we have enough results
            if len(filtered) >= self.max_results:
                break
        
        return filtered
    
    def get_accession_list(self) -> List[str]:
        """Get list of discovered accessions"""
        return [genome.accession for genome in self.discovered_genomes.values()]
    
    def generate_discovery_report(self) -> str:
        """Generate detailed discovery report"""
        genomes = list(self.discovered_genomes.values())
        
        report = f"""
NCBI Genome Discovery Report
===========================
Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}

Summary:
--------
Total genomes discovered: {len(genomes)}
Average genome size: {sum(g.length for g in genomes) / len(genomes) / 1_000_000:.2f} MB
Size range: {min(g.length for g in genomes) / 1_000_000:.2f} - {max(g.length for g in genomes) / 1_000_000:.2f} MB

Discovered Genomes:
------------------
"""
        
        for genome in genomes:
            report += f"{genome.accession}\t{genome.organism[:30]}\t{genome.length / 1_000_000:.2f} MB\t{genome.title[:50]}...\n"
        
        return report


class URLBasedGenomeDiscovery:
    """
    Complete URL-based genome discovery workflow
    """
    
    def __init__(self, 
                 email: str,
                 api_key: Optional[str] = None,
                 max_results: int = 10,
                 logger: Optional[logging.Logger] = None):
        
        self.url_parser = NCBIUrlParser(logger)
        self.genome_discovery = NCBIGenomeDiscovery(email, api_key, max_results, logger)
        self.logger = logger or self._setup_logger()
    
    def _setup_logger(self) -> logging.Logger:
        """Setup main workflow logger"""
        logger = logging.getLogger("URLBasedGenomeDiscovery")
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(levelname)s %(name)s: %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        return logger
    
    async def discover_from_url(self, ncbi_url: str) -> Tuple[List[GenomeMetadata], str]:
        """
        Complete workflow: URL → search term → genome discovery
        
        Args:
            ncbi_url: NCBI search URL
            
        Returns:
            Tuple of (discovered genomes, discovery report)
        """
        self.logger.info("Starting URL-based genome discovery")
        
        try:
            # Step 1: Parse URL and extract search term
            search_term = self.url_parser.parse_search_url(ncbi_url)
            
            # Step 2: Discover genomes using search term
            genomes = await self.genome_discovery.discover_genomes(search_term)
            
            # Step 3: Generate discovery report
            report = self.genome_discovery.generate_discovery_report()
            
            self.logger.info(f"Discovery complete: {len(genomes)} genomes ready for download")
            return genomes, report
            
        except Exception as e:
            self.logger.error(f"URL-based discovery failed: {e}")
            raise


# Example usage and testing
async def main():
    """Test the URL-based genome discovery"""
    
    # Test URL
    test_url = "https://www.ncbi.nlm.nih.gov/nuccore/?term=(Escherichia+coli+and+macrolide+resistance)+AND+%22Escherichia+coli%22%5Bporgn%3A__txid562%5D+and+complete+genome"
    
    # Create discovery engine
    discovery = URLBasedGenomeDiscovery(
        email="vihaankulkarni29@gmail.com",
        api_key="ef7622c2e716fa317fe04d24c42904211107",
        max_results=10
    )
    
    try:
        # Run discovery
        genomes, report = await discovery.discover_from_url(test_url)
        
        print(f"\n🎉 Discovered {len(genomes)} genomes!")
        print("\nFirst few genomes:")
        for i, genome in enumerate(genomes[:5]):
            print(f"{i+1}. {genome.accession} - {genome.organism} ({genome.length/1_000_000:.2f} MB)")
        
        print(f"\nAccession list for download:")
        accessions = discovery.genome_discovery.get_accession_list()
        for acc in accessions:
            print(f"  {acc}")
        
    except Exception as e:
        print(f"❌ Discovery failed: {e}")


if __name__ == "__main__":
    asyncio.run(main())