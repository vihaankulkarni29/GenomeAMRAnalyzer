import sys; import os; sys.path.append(os.path.join(os.path.dirname(__file__), '\\\\\\\\\src')); from configuration_manager import config_manager
#!/usr/bin/env python3
"""
EnhancedNCBIProteinHarvester - Scientifically validated metadata collection
Collects comprehensive taxonomic, clinical, and biological metadata for AMR research

Author: MetaDataHarvester Pipeline
Version: 2.0 - Species-Specific Enhanced
"""

import os
import sys
import logging
import time
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from dataclasses import dataclass, asdict
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import pandas as pd
import json
import re
from urllib.parse import urlparse
import requests


@dataclass
class EnhancedProteinMetadata:
    """Comprehensive metadata for AMR research"""
    # Basic identifiers
    accession: str
    protein_name: str
    description: str
    sequence_length: int

    # Enhanced taxonomic data
    organism_name: str
    genus: str
    species: str
    strain: str
    taxonomy_id: str
    subspecies: Optional[str] = None

    # Clinical metadata
    mic_data: Dict[str, float] = None  # antibiotic: MIC value
    resistance_phenotype: str = "unknown"  # resistant/intermediate/susceptible
    pubmed_studies: List[str] = None
    isolation_source: str = "unknown"
    geographic_location: str = "unknown"
    collection_date: str = "unknown"

    # Protein-specific metadata
    protein_family: str = "unknown"  # RND, ABC, MFS
    functional_annotation: List[str] = None  # GO terms, domains
    reference_mapping: Dict = None  # species-specific reference info

    # Quality control
    data_completeness: float = 0.0
    validation_warnings: List[str] = None


class EnhancedNCBIProteinHarvester:
    """
    Enhanced NCBI data harvester with comprehensive metadata collection
    """

    def __init__(self, email: str, output_dir: str = "enhanced_protein_data"):
        self.email = email
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)

        # Setup logging
        self.logger = logging.getLogger(__name__)
        self._setup_logging()

        # Initialize NCBI
        Entrez.email = email

        # Cache for API calls
        self.cache_dir = self.output_dir / "cache"
        self.cache_dir.mkdir(exist_ok=True)

    def _setup_logging(self):
        """Setup comprehensive logging"""
        log_file = self.output_dir / "harvester.log"
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )

    def harvest_enhanced_metadata(self, accession_list: List[str]) -> List[EnhancedProteinMetadata]:
        """
        Harvest comprehensive metadata for a list of protein accessions

        Args:
            accession_list: List of NCBI protein accession numbers

        Returns:
            List of enhanced metadata objects
        """
        self.logger.info(f"Starting enhanced metadata harvest for {len(accession_list)} accessions")

        enhanced_metadata = []

        for i, accession in enumerate(accession_list):
            self.logger.info(f"Processing {accession} ({i+1}/{len(accession_list)})")

            try:
                metadata = self._extract_enhanced_metadata(accession)
                enhanced_metadata.append(metadata)

                # Rate limiting for NCBI API
                time.sleep(0.5)

            except Exception as e:
                self.logger.error(f"Failed to process {accession}: {e}")
                # Create minimal metadata for failed accessions
                metadata = EnhancedProteinMetadata(
                    accession=accession,
                    protein_name="unknown",
                    description="Failed to retrieve",
                    sequence_length=0,
                    organism_name="unknown",
                    genus="unknown",
                    species="unknown",
                    strain="unknown",
                    taxonomy_id="unknown",
                    validation_warnings=[f"Metadata extraction failed: {str(e)}"]
                )
                enhanced_metadata.append(metadata)

        self.logger.info(f"Completed metadata harvest. Processed {len(enhanced_metadata)} accessions")
        return enhanced_metadata

    def _extract_enhanced_metadata(self, accession: str) -> EnhancedProteinMetadata:
        """Extract comprehensive metadata from NCBI"""
        # Get protein record
        protein_record = self._fetch_protein_record(accession)

        # Parse basic information
        basic_info = self._parse_basic_info(protein_record, accession)

        # Extract taxonomic information
        taxonomic_info = self._extract_taxonomic_info(protein_record)

        # Search for clinical data
        clinical_info = self._extract_clinical_info(accession, protein_record)

        # Determine protein family and reference mapping
        protein_info = self._classify_protein_family(protein_record, taxonomic_info)

        # Calculate data completeness
        completeness = self._calculate_data_completeness(basic_info, taxonomic_info, clinical_info)

        # Generate validation warnings
        warnings = self._generate_validation_warnings(basic_info, taxonomic_info, clinical_info)

        # Create enhanced metadata object
        metadata = EnhancedProteinMetadata(
            accession=accession,
            protein_name=basic_info['protein_name'],
            description=basic_info['description'],
            sequence_length=basic_info['sequence_length'],
            organism_name=taxonomic_info['organism_name'],
            genus=taxonomic_info['genus'],
            species=taxonomic_info['species'],
            strain=taxonomic_info['strain'],
            taxonomy_id=taxonomic_info['taxonomy_id'],
            subspecies=taxonomic_info.get('subspecies'),
            mic_data=clinical_info.get('mic_data'),
            resistance_phenotype=clinical_info.get('resistance_phenotype', 'unknown'),
            pubmed_studies=clinical_info.get('pubmed_studies', []),
            isolation_source=clinical_info.get('isolation_source', 'unknown'),
            geographic_location=clinical_info.get('geographic_location', 'unknown'),
            collection_date=clinical_info.get('collection_date', 'unknown'),
            protein_family=protein_info['protein_family'],
            functional_annotation=protein_info.get('functional_annotation', []),
            reference_mapping=protein_info.get('reference_mapping'),
            data_completeness=completeness,
            validation_warnings=warnings
        )

        return metadata

    def _fetch_protein_record(self, accession: str) -> str:
        """Fetch protein record from NCBI"""
        cache_file = self.cache_dir / f"{accession}.gb"

        if cache_file.exists():
            with open(cache_file, 'r') as f:
                return f.read()

        try:
            handle = Entrez.efetch(db="protein", id=accession, rettype="gb", retmode="text")
            record = handle.read()
            handle.close()

            # Cache the record
            with open(cache_file, 'w') as f:
                f.write(record)

            return record

        except Exception as e:
            self.logger.warning(f"Failed to fetch {accession}: {e}")
            raise

    def _parse_basic_info(self, protein_record: str, accession: str) -> Dict:
        """Parse basic protein information"""
        lines = protein_record.split('\n')

        # Extract protein name
        protein_name = "unknown"
        for line in lines:
            if line.startswith('PRODUCT'):
                protein_name = line.replace('PRODUCT', '').strip()
                break

        # Extract description
        description = "unknown"
        for line in lines:
            if line.startswith('DEFINITION'):
                description = line.replace('DEFINITION', '').strip()
                break

        # Get sequence length (this is approximate from the record)
        sequence_length = 0
        in_sequence = False
        sequence_lines = []

        for line in lines:
            if line.startswith('ORIGIN'):
                in_sequence = True
                continue
            elif in_sequence and line.strip() and not line.startswith('//'):
                # Remove numbers and spaces
                seq_part = re.sub(r'[0-9\s]', '', line)
                sequence_lines.append(seq_part)

        if sequence_lines:
            sequence_length = len(''.join(sequence_lines))

        return {
            'protein_name': protein_name,
            'description': description,
            'sequence_length': sequence_length
        }

    def _extract_taxonomic_info(self, protein_record: str) -> Dict:
        """Extract detailed taxonomic information"""
        lines = protein_record.split('\n')

        organism_name = "unknown"
        taxonomy_id = "unknown"

        for line in lines:
            if line.startswith('  ORGANISM'):
                organism_name = line.replace('  ORGANISM', '').strip()
            elif 'db_xref="taxon:' in line:
                tax_match = re.search(r'db_xref="taxon:(\d+)"', line)
                if tax_match:
                    taxonomy_id = tax_match.group(1)

        # Parse organism name into components
        genus, species, strain, subspecies = self._parse_organism_name(organism_name)

        return {
            'organism_name': organism_name,
            'genus': genus,
            'species': species,
            'strain': strain,
            'taxonomy_id': taxonomy_id,
            'subspecies': subspecies
        }

    def _parse_organism_name(self, organism_name: str) -> Tuple[str, str, str, Optional[str]]:
        """Parse organism name into taxonomic components"""
        # Examples:
        # "Escherichia coli str. K-12 substr. MG1655"
        # "Pseudomonas aeruginosa PAO1"
        # "Klebsiella pneumoniae subsp. pneumoniae"

        parts = organism_name.split()

        genus = parts[0] if len(parts) > 0 else "unknown"
        species = parts[1] if len(parts) > 1 else "unknown"

        strain = ""
        subspecies = None

        # Look for strain information
        if "str." in organism_name:
            str_index = organism_name.find("str.")
            strain = organism_name[str_index:].strip()
        elif "subsp." in organism_name:
            subsp_index = organism_name.find("subsp.")
            subspecies_part = organism_name[subsp_index:].strip()
            subspecies = subspecies_part.replace("subsp.", "").strip()

        return genus, species, strain, subspecies

    def _extract_clinical_info(self, accession: str, protein_record: str) -> Dict:
        """Extract clinical and phenotypic information"""
        clinical_info = {
            'mic_data': {},
            'resistance_phenotype': 'unknown',
            'pubmed_studies': [],
            'isolation_source': 'unknown',
            'geographic_location': 'unknown',
            'collection_date': 'unknown'
        }

        # Search for associated BioSample/BioProject data
        biosample_data = self._get_biosample_data(accession)
        if biosample_data:
            clinical_info.update(biosample_data)

        # Search for PubMed studies
        pubmed_ids = self._find_associated_pubmed_studies(accession, protein_record)
        clinical_info['pubmed_studies'] = pubmed_ids

        return clinical_info

    def _get_biosample_data(self, accession: str) -> Optional[Dict]:
        """Get clinical data from BioSample records"""
        try:
            # Find BioSample links
            handle = Entrez.elink(dbfrom="protein", db="biosample", id=accession)
            link_results = Entrez.read(handle)
            handle.close()

            if link_results and link_results[0]['LinkSetDb']:
                biosample_ids = [
                    link['Id'] for link in link_results[0]['LinkSetDb'][0]['Link']
                ]

                if biosample_ids:
                    # Fetch first BioSample record
                    handle = Entrez.efetch(db="biosample", id=biosample_ids[0], rettype="xml")
                    biosample_xml = handle.read()
                    handle.close()

                    return self._parse_biosample_xml(biosample_xml)

        except Exception as e:
            self.logger.debug(f"Could not fetch BioSample data for {accession}: {e}")

        return None

    def _parse_biosample_xml(self, xml_data: str) -> Dict:
        """Parse BioSample XML for clinical information"""
        # This is a simplified parser - in practice, you'd use proper XML parsing
        clinical_info = {}

        # Look for MIC data
        if "MIC" in xml_data.upper():
            # Extract MIC values (simplified)
            mic_matches = re.findall(r'MIC[^:]*:?\s*([^<\n]+)', xml_data, re.IGNORECASE)
            if mic_matches:
                clinical_info['mic_data'] = {'extracted': mic_matches[0]}

        # Look for isolation source
        if "isolation source" in xml_data.lower():
            source_match = re.search(r'isolation source[^:]*:?\s*([^<\n]+)', xml_data, re.IGNORECASE)
            if source_match:
                clinical_info['isolation_source'] = source_match.group(1).strip()

        # Look for geographic location
        if "geographic location" in xml_data.lower():
            location_match = re.search(r'geographic location[^:]*:?\s*([^<\n]+)', xml_data, re.IGNORECASE)
            if location_match:
                clinical_info['geographic_location'] = location_match.group(1).strip()

        return clinical_info

    def _find_associated_pubmed_studies(self, accession: str, protein_record: str) -> List[str]:
        """Find PubMed studies associated with the protein"""
        pubmed_ids = []

        try:
            # Search PubMed for studies mentioning this accession
            query = f"{accession}[Accession] OR {accession}[Protein]"
            handle = Entrez.esearch(db="pubmed", term=query, retmax=5)
            search_results = Entrez.read(handle)
            handle.close()

            pubmed_ids = search_results.get('IdList', [])

        except Exception as e:
            self.logger.debug(f"Could not find PubMed studies for {accession}: {e}")

        return pubmed_ids

    def _classify_protein_family(self, protein_record: str, taxonomic_info: Dict) -> Dict:
        """Classify protein family and determine reference mapping"""
        description = protein_record.lower()

        protein_family = "unknown"
        functional_annotation = []

        # Classify efflux pump family
        if any(term in description for term in ['rnd', 'resistance-nodulation-division']):
            protein_family = "RND"
            functional_annotation.append("GO:0043190")  # efflux
        elif any(term in description for term in ['abc', 'atp-binding cassette']):
            protein_family = "ABC"
            functional_annotation.append("GO:0043190")
        elif any(term in description for term in ['mfs', 'major facilitator']):
            protein_family = "MFS"
            functional_annotation.append("GO:0043190")

        # Specific protein identification
        protein_name = "unknown"
        if "acra" in description or "acr a" in description:
            protein_name = "AcrA"
        elif "acrb" in description or "acr b" in description:
            protein_name = "AcrB"
        elif "tolc" in description:
            protein_name = "TolC"

        # Determine reference mapping
        reference_mapping = self._determine_reference_mapping(
            taxonomic_info['genus'], protein_name, protein_family
        )

        return {
            'protein_family': protein_family,
            'protein_name': protein_name,
            'functional_annotation': functional_annotation,
            'reference_mapping': reference_mapping
        }

    def _determine_reference_mapping(self, genus: str, protein_name: str, protein_family: str) -> Dict:
        """Determine which reference sequence to use for alignment"""
        # Priority 1: Exact genus match
        reference_file = f"references/{genus}/{protein_name}_reference.fasta"

        # Priority 2: Phylogenetically close genus
        if genus == "Escherichia":
            close_genera = ["Shigella", "Salmonella", "Klebsiella"]
        elif genus == "Pseudomonas":
            close_genera = ["Acinetobacter", "Burkholderia"]
        else:
            close_genera = ["Escherichia"]  # Default fallback

        for close_genus in close_genera:
            alt_reference = f"references/{close_genus}/{protein_name}_reference.fasta"
            if Path(alt_reference).exists():
                reference_file = alt_reference
                break

        return {
            'genus': genus,
            'protein_name': protein_name,
            'reference_file': reference_file,
            'selection_method': 'genus_specific' if Path(reference_file).exists() else 'phylogenetic_fallback'
        }

    def _calculate_data_completeness(self, basic_info: Dict, taxonomic_info: Dict,
                                   clinical_info: Dict) -> float:
        """Calculate data completeness score"""
        score = 0.0
        max_score = 10.0

        # Basic information (2 points)
        if basic_info['protein_name'] != 'unknown':
            score += 1
        if basic_info['sequence_length'] > 0:
            score += 1

        # Taxonomic information (3 points)
        if taxonomic_info['genus'] != 'unknown':
            score += 1
        if taxonomic_info['species'] != 'unknown':
            score += 1
        if taxonomic_info['taxonomy_id'] != 'unknown':
            score += 1

        # Clinical information (3 points)
        if clinical_info.get('isolation_source') != 'unknown':
            score += 1
        if clinical_info.get('geographic_location') != 'unknown':
            score += 1
        if clinical_info.get('pubmed_studies'):
            score += 1

        # Protein classification (2 points)
        if clinical_info.get('protein_family') != 'unknown':
            score += 1
        if clinical_info.get('reference_mapping'):
            score += 1

        return score / max_score

    def _generate_validation_warnings(self, basic_info: Dict, taxonomic_info: Dict,
                                    clinical_info: Dict) -> List[str]:
        """Generate validation warnings"""
        warnings = []

        if basic_info['protein_name'] == 'unknown':
            warnings.append("Protein name not identified")

        if taxonomic_info['genus'] == 'unknown':
            warnings.append("Genus not determined")

        if not clinical_info.get('pubmed_studies'):
            warnings.append("No associated PubMed studies found")

        if clinical_info.get('resistance_phenotype') == 'unknown':
            warnings.append("Resistance phenotype not available")

        return warnings

    def save_enhanced_metadata(self, metadata_list: List[EnhancedProteinMetadata],
                              output_file: str = "enhanced_metadata.csv") -> None:
        """Save enhanced metadata to CSV"""
        output_path = self.output_dir / output_file

        # Convert to dictionary format
        data = []
        for metadata in metadata_list:
            row = asdict(metadata)
            # Convert lists/dicts to JSON strings for CSV storage
            if row['mic_data']:
                row['mic_data'] = json.dumps(row['mic_data'])
            if row['pubmed_studies']:
                row['pubmed_studies'] = json.dumps(row['pubmed_studies'])
            if row['functional_annotation']:
                row['functional_annotation'] = json.dumps(row['functional_annotation'])
            if row['reference_mapping']:
                row['reference_mapping'] = json.dumps(row['reference_mapping'])
            if row['validation_warnings']:
                row['validation_warnings'] = json.dumps(row['validation_warnings'])

            data.append(row)

        df = pd.DataFrame(data)
        df.to_csv(output_path, index=False)

        self.logger.info(f"Saved enhanced metadata for {len(metadata_list)} proteins to {output_path}")

    def generate_metadata_report(self, metadata_list: List[EnhancedProteinMetadata]) -> None:
        """Generate comprehensive metadata report"""
        report_file = self.output_dir / "metadata_report.txt"

        with open(report_file, 'w') as f:
            f.write("ENHANCED NCBI PROTEIN HARVESTER REPORT\n")
            f.write("=" * 50 + "\n\n")

            f.write(f"Total proteins processed: {len(metadata_list)}\n\n")

            # Taxonomic distribution
            genera = [m.genus for m in metadata_list if m.genus != 'unknown']
            genus_counts = pd.Series(genera).value_counts()

            f.write("TAXONOMIC DISTRIBUTION:\n")
            f.write("-" * 25 + "\n")
            for genus, count in genus_counts.head(10).items():
                f.write(f"{genus}: {count}\n")
            f.write("\n")

            # Protein family distribution
            families = [m.protein_family for m in metadata_list if m.protein_family != 'unknown']
            family_counts = pd.Series(families).value_counts()

            f.write("PROTEIN FAMILY DISTRIBUTION:\n")
            f.write("-" * 30 + "\n")
            for family, count in family_counts.items():
                f.write(f"{family}: {count}\n")
            f.write("\n")

            # Data completeness
            completeness_scores = [m.data_completeness for m in metadata_list]
            avg_completeness = sum(completeness_scores) / len(completeness_scores)

            f.write("DATA COMPLETENESS:\n")
            f.write("-" * 20 + "\n")
            f.write(f"Average completeness: {avg_completeness:.1f}\n")
            f.write(f"Min completeness: {min(completeness_scores):.1f}\n")
            f.write(f"Max completeness: {max(completeness_scores):.1f}\n")
            f.write("\n")

            # Clinical data availability
            with_pubmed = sum(1 for m in metadata_list if m.pubmed_studies)
            with_mic = sum(1 for m in metadata_list if m.mic_data)

            f.write("CLINICAL DATA AVAILABILITY:\n")
            f.write("-" * 30 + "\n")
            f.write(f"Proteins with PubMed studies: {with_pubmed}\n")
            f.write(f"Proteins with MIC data: {with_mic}\n")

        self.logger.info(f"Generated metadata report: {report_file}")


def main():
    """Main entry point"""
    import argparse

    parser = argparse.ArgumentParser(
        description="Enhanced NCBI protein harvester with comprehensive metadata collection"
    )

    parser.add_argument(
        "--accession-file",
        required=True,
        help="File containing protein accession numbers (one per line)"
    )

    parser.add_argument(
        "--email",
        required=True,
        help="NCBI email for API access"
    )

    parser.add_argument(
        "--output-dir",
        default="enhanced_protein_data",
        help="Output directory for results"
    )

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    # Read accession list
    with open(args.accession_file, 'r') as f:
        accessions = [line.strip() for line in f if line.strip()]

    # Initialize harvester
    harvester = EnhancedNCBIProteinHarvester(args.email, args.output_dir)

    # Harvest enhanced metadata
    metadata_list = harvester.harvest_enhanced_metadata(accessions)

    # Save results
    harvester.save_enhanced_metadata(metadata_list)
    harvester.generate_metadata_report(metadata_list)

    print(f"\nEnhanced metadata harvest complete!")
    print(f"Processed {len(metadata_list)} proteins")
    print(f"Results saved to: {args.output_dir}")


if __name__ == "__main__":
    main()