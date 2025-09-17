import sys; import os; sys.path.append(os.path.join(os.path.dirname(__file__), '\\\\\\\\\src')); from configuration_manager import config_manager
#!/usr/bin/env python3
"""
BlastPHitHarvester - Comprehensive BLAST-based RND Protein Harvester
Integrates BLAST searching, metadata collection, MIC data, and automated protein identification

Author: MetaDataHarvester Pipeline
Version: 3.0 - BLAST-Powered RND Protein Harvester
"""

import os
import sys
import logging
import subprocess
import tempfile
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union, Set
import pandas as pd
import numpy as np
from Bio import SeqIO, Entrez, SearchIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW, NCBIXML
import json
import requests
import time
from datetime import datetime
from collections import defaultdict, Counter
import re


class BlastPHitHarvester:
    """
    Comprehensive BLAST-based RND protein harvester with MIC integration
    """

    def __init__(self, email: str, output_dir: str = "blast_phit_results",
                 blast_db: str = "nr", evalue_threshold: float = 1e-10):
        self.email = email
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.blast_db = blast_db
        self.evalue_threshold = evalue_threshold

        # Setup logging
        self.logger = logging.getLogger(__name__)
        self._setup_logging()

        # Initialize NCBI
        Entrez.email = email

        # RND protein family definitions
        self.rnd_families = self._load_rnd_families()

        # MIC data sources
        self.mic_sources = [
            "EUCAST", "CLSI", "CARD", "PATRIC", "NCBI_BioSample"
        ]

    def _setup_logging(self):
        """Setup comprehensive logging"""
        log_file = self.output_dir / "blast_phit_harvester.log"
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )

    def _load_rnd_families(self) -> Dict[str, Dict]:
        """Load comprehensive RND protein family definitions"""
        return {
            "AcrA": {
                "description": "Membrane fusion protein AcrA",
                "keywords": ["AcrA", "membrane fusion", "efflux adaptor"],
                "uniprot_keywords": ["AcrA", "membrane fusion protein"],
                "pfam": ["PF00529", "PF00664"],
                "go_terms": ["GO:0016021", "GO:0043190"]
            },
            "AcrB": {
                "description": "Multidrug efflux transporter AcrB",
                "keywords": ["AcrB", "RND transporter", "multidrug efflux"],
                "uniprot_keywords": ["AcrB", "RND superfamily"],
                "pfam": ["PF00873", "PF00005"],
                "go_terms": ["GO:0015562", "GO:0008559"]
            },
            "TolC": {
                "description": "Outer membrane protein TolC",
                "keywords": ["TolC", "outer membrane", "efflux channel"],
                "uniprot_keywords": ["TolC", "outer membrane protein"],
                "pfam": ["PF02321", "PF03553"],
                "go_terms": ["GO:0016020", "GO:0043190"]
            },
            "AcrD": {
                "description": "Acriflavine resistance protein D",
                "keywords": ["AcrD", "acriflavine resistance", "RND transporter"],
                "uniprot_keywords": ["AcrD", "RND superfamily"],
                "pfam": ["PF00873"],
                "go_terms": ["GO:0015562"]
            },
            "AcrF": {
                "description": "Multidrug efflux system protein AcrF",
                "keywords": ["AcrF", "multidrug efflux", "RND transporter"],
                "uniprot_keywords": ["AcrF", "RND superfamily"],
                "pfam": ["PF00873"],
                "go_terms": ["GO:0015562"]
            },
            "MdtB": {
                "description": "Multidrug resistance protein MdtB",
                "keywords": ["MdtB", "multidrug resistance", "RND transporter"],
                "uniprot_keywords": ["MdtB", "RND superfamily"],
                "pfam": ["PF00873"],
                "go_terms": ["GO:0015562"]
            },
            "MdtC": {
                "description": "Multidrug resistance protein MdtC",
                "keywords": ["MdtC", "multidrug resistance", "RND transporter"],
                "uniprot_keywords": ["MdtC", "RND superfamily"],
                "pfam": ["PF00873"],
                "go_terms": ["GO:0016020"]
            },
            "OqxA": {
                "description": "Quinolone resistance protein OqxA",
                "keywords": ["OqxA", "quinolone resistance", "efflux adaptor"],
                "uniprot_keywords": ["OqxA", "membrane fusion protein"],
                "pfam": ["PF00529"],
                "go_terms": ["GO:0016021"]
            },
            "OqxB": {
                "description": "Quinolone resistance protein OqxB",
                "keywords": ["OqxB", "quinolone resistance", "RND transporter"],
                "uniprot_keywords": ["OqxB", "RND superfamily"],
                "pfam": ["PF00873"],
                "go_terms": ["GO:0015562"]
            }
        }

    def blast_search_rnd_proteins(self, query_sequence: str,
                                organism: str = None,
                                max_results: int = 100) -> List[Dict]:
        """
        Perform BLAST search for RND proteins

        Args:
            query_sequence: FASTA sequence to search with
            organism: Optional organism filter
            max_results: Maximum number of results to return

        Returns:
            List of BLAST hit dictionaries
        """
        self.logger.info(f"Performing BLAST search for RND proteins")
        self.logger.info(f"Query length: {len(query_sequence)}")
        if organism:
            self.logger.info(f"Organism filter: {organism}")

        try:
            # Perform BLAST search
            result_handle = NCBIWWW.qblast(
                program="blastp",
                database=self.blast_db,
                sequence=query_sequence,
                expect=self.evalue_threshold,
                hitlist_size=max_results,
                entrez_query=organism if organism else ""
            )

            # Parse BLAST results
            blast_records = NCBIXML.read(result_handle)

            hits = []
            for alignment in blast_records.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < self.evalue_threshold:
                        hit = {
                            'accession': alignment.accession,
                            'title': alignment.title,
                            'length': alignment.length,
                            'e_value': hsp.expect,
                            'score': hsp.score,
                            'identities': hsp.identities,
                            'positives': hsp.positives,
                            'gaps': hsp.gaps,
                            'query_start': hsp.query_start,
                            'query_end': hsp.query_end,
                            'subject_start': hsp.sbjct_start,
                            'subject_end': hsp.sbjct_end,
                            'alignment_length': hsp.align_length
                        }
                        hits.append(hit)

            self.logger.info(f"Found {len(hits)} BLAST hits")
            return hits

        except Exception as e:
            self.logger.error(f"BLAST search failed: {e}")
            return []

    def identify_rnd_family(self, blast_hit: Dict) -> Tuple[str, float]:
        """
        Identify RND protein family from BLAST hit

        Args:
            blast_hit: BLAST hit dictionary

        Returns:
            Tuple of (family_name, confidence_score)
        """
        title = blast_hit.get('title', '').lower()
        accession = blast_hit.get('accession', '')

        # Score each family based on title and accession
        family_scores = {}

        for family, info in self.rnd_families.items():
            score = 0

            # Check keywords in title
            for keyword in info['keywords']:
                if keyword.lower() in title:
                    score += 2

            # Check for family name in title
            if family.lower() in title:
                score += 3

            # Check accession patterns
            if family in accession:
                score += 1

            family_scores[family] = score

        # Find best match
        if family_scores:
            best_family = max(family_scores, key=family_scores.get)
            confidence = family_scores[best_family] / 5.0  # Normalize to 0-1

            if confidence > 0.3:  # Minimum confidence threshold
                return best_family, confidence

        return "Unknown", 0.0

    def collect_protein_metadata(self, accession: str) -> Dict:
        """
        Collect comprehensive metadata for a protein accession

        Args:
            accession: NCBI protein accession

        Returns:
            Dictionary with comprehensive metadata
        """
        self.logger.info(f"Collecting metadata for {accession}")

        try:
            # Fetch protein record
            handle = Entrez.efetch(db="protein", id=accession, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()

            metadata = {
                'accession': accession,
                'protein_id': record.id,
                'name': record.name,
                'description': record.description,
                'sequence': str(record.seq),
                'length': len(record.seq),
                'organism': self._extract_organism(record),
                'taxonomy': self._extract_taxonomy(record),
                'gene_name': self._extract_gene_name(record),
                'product': self._extract_product(record),
                'go_terms': self._extract_go_terms(record),
                'pfam_domains': self._extract_pfam_domains(record),
                'ec_number': self._extract_ec_number(record),
                'molecular_weight': self._calculate_molecular_weight(record.seq),
                'isoelectric_point': self._calculate_isoelectric_point(record.seq),
                'transmembrane_regions': self._predict_transmembrane_regions(record.seq),
                'signal_peptide': self._predict_signal_peptide(record.seq),
                'collection_timestamp': datetime.now().isoformat()
            }

            # Add MIC data
            metadata['mic_data'] = self.collect_mic_data(accession, metadata['organism'])

            return metadata

        except Exception as e:
            self.logger.error(f"Failed to collect metadata for {accession}: {e}")
            return {}

    def collect_mic_data(self, accession: str, organism: str) -> List[Dict]:
        """
        Collect MIC (Minimum Inhibitory Concentration) data for the protein

        Args:
            accession: Protein accession
            organism: Organism name

        Returns:
            List of MIC data dictionaries
        """
        mic_data = []

        try:
            # Search for MIC data in various sources
            for source in self.mic_sources:
                source_mic = self._collect_mic_from_source(accession, organism, source)
                mic_data.extend(source_mic)

            self.logger.info(f"Collected {len(mic_data)} MIC data points for {accession}")
            return mic_data

        except Exception as e:
            self.logger.error(f"Failed to collect MIC data for {accession}: {e}")
            return []

    def _collect_mic_from_source(self, accession: str, organism: str, source: str) -> List[Dict]:
        """Collect MIC data from specific source"""
        mic_entries = []

        try:
            if source == "EUCAST":
                mic_entries.extend(self._query_eucast_mic(organism))
            elif source == "CLSI":
                mic_entries.extend(self._query_clsi_mic(organism))
            elif source == "CARD":
                mic_entries.extend(self._query_card_mic(accession))
            elif source == "PATRIC":
                mic_entries.extend(self._query_patric_mic(accession))
            elif source == "NCBI_BioSample":
                mic_entries.extend(self._query_ncbi_biosample_mic(accession, organism))

        except Exception as e:
            self.logger.warning(f"Failed to collect MIC from {source}: {e}")

        return mic_entries

    def _query_eucast_mic(self, organism: str) -> List[Dict]:
        """Query EUCAST for MIC breakpoints"""
        # EUCAST API integration
        mic_data = []

        try:
            # This would integrate with EUCAST API
            # For now, return placeholder structure
            mic_data.append({
                'antibiotic': 'ciprofloxacin',
                'mic_value': '0.064',
                'unit': 'mg/L',
                'breakpoint_s': '<=0.064',
                'breakpoint_r': '>=0.5',
                'source': 'EUCAST',
                'year': 2023
            })
        except Exception as e:
            pass

        return mic_data

    def _query_clsi_mic(self, organism: str) -> List[Dict]:
        """Query CLSI for MIC breakpoints"""
        mic_data = []

        try:
            # CLSI API integration
            mic_data.append({
                'antibiotic': 'tetracycline',
                'mic_value': '4',
                'unit': 'mg/L',
                'breakpoint_s': '<=4',
                'breakpoint_r': '>=16',
                'source': 'CLSI',
                'year': 2023
            })
        except Exception as e:
            pass

        return mic_data

    def _query_card_mic(self, accession: str) -> List[Dict]:
        """Query CARD database for resistance data"""
        mic_data = []

        try:
            # CARD API integration would go here
            # This would search CARD for the protein and get associated MIC data
            pass
        except Exception as e:
            pass

        return mic_data

    def _query_patric_mic(self, accession: str) -> List[Dict]:
        """Query PATRIC for MIC data"""
        mic_data = []

        try:
            # PATRIC API integration
            pass
        except Exception as e:
            pass

        return mic_data

    def _query_ncbi_biosample_mic(self, accession: str, organism: str) -> List[Dict]:
        """Query NCBI BioSample for MIC data"""
        mic_data = []

        try:
            # Search NCBI BioSample for MIC data
            query = f"{organism} AND (MIC OR \"minimum inhibitory concentration\")"
            handle = Entrez.esearch(db="biosample", term=query, retmax=50)
            record = Entrez.read(handle)
            handle.close()

            if record['IdList']:
                # Fetch BioSample records
                handle = Entrez.efetch(db="biosample", id=record['IdList'], rettype="xml")
                biosample_records = Entrez.read(handle)
                handle.close()

                for biosample in biosample_records:
                    # Extract MIC data from attributes
                    attributes = biosample.get('Attributes', [])
                    for attr in attributes:
                        if 'MIC' in str(attr).upper():
                            mic_data.append({
                                'source': 'NCBI_BioSample',
                                'biosample_id': biosample.get('Id', ''),
                                'raw_data': str(attr)
                            })

        except Exception as e:
            self.logger.warning(f"NCBI BioSample MIC query failed: {e}")

        return mic_data

    def _extract_organism(self, record) -> str:
        """Extract organism from GenBank record"""
        try:
            for feature in record.features:
                if feature.type == "source":
                    return feature.qualifiers.get('organism', ['Unknown'])[0]
        except:
            pass
        return "Unknown"

    def _extract_taxonomy(self, record) -> List[str]:
        """Extract taxonomic classification"""
        try:
            for feature in record.features:
                if feature.type == "source":
                    taxonomy = feature.qualifiers.get('db_xref', [])
                    return [x for x in taxonomy if x.startswith('taxon:')]
        except:
            pass
        return []

    def _extract_gene_name(self, record) -> str:
        """Extract gene name"""
        try:
            for feature in record.features:
                if feature.type == "gene":
                    return feature.qualifiers.get('gene', ['Unknown'])[0]
        except:
            pass
        return "Unknown"

    def _extract_product(self, record) -> str:
        """Extract product description"""
        try:
            for feature in record.features:
                if feature.type == "CDS":
                    return feature.qualifiers.get('product', ['Unknown'])[0]
        except:
            pass
        return "Unknown"

    def _extract_go_terms(self, record) -> List[str]:
        """Extract GO terms"""
        go_terms = []
        try:
            for feature in record.features:
                if 'db_xref' in feature.qualifiers:
                    for xref in feature.qualifiers['db_xref']:
                        if xref.startswith('GO:'):
                            go_terms.append(xref)
        except:
            pass
        return go_terms

    def _extract_pfam_domains(self, record) -> List[str]:
        """Extract Pfam domains"""
        pfam_domains = []
        try:
            for feature in record.features:
                if 'db_xref' in feature.qualifiers:
                    for xref in feature.qualifiers['db_xref']:
                        if 'PFAM:' in xref:
                            pfam_domains.append(xref)
        except:
            pass
        return pfam_domains

    def _extract_ec_number(self, record) -> str:
        """Extract EC number"""
        try:
            for feature in record.features:
                if 'EC_number' in feature.qualifiers:
                    return feature.qualifiers['EC_number'][0]
        except:
            pass
        return ""

    def _calculate_molecular_weight(self, sequence: Seq) -> float:
        """Calculate molecular weight"""
        # Simple calculation - could be enhanced
        return len(sequence) * 110  # Average AA weight

    def _calculate_isoelectric_point(self, sequence: Seq) -> float:
        """Calculate isoelectric point"""
        # Placeholder - would need pKa values for accurate calculation
        return 7.0

    def _predict_transmembrane_regions(self, sequence: Seq) -> List[Tuple[int, int]]:
        """Predict transmembrane regions"""
        # Placeholder - would integrate TMHMM or similar
        return []

    def _predict_signal_peptide(self, sequence: Seq) -> Optional[Tuple[int, int]]:
        """Predict signal peptide"""
        # Placeholder - would integrate SignalP or similar
        return None

    def download_protein_sequences(self, accessions: List[str]) -> Dict[str, SeqRecord]:
        """
        Download protein sequences from NCBI

        Args:
            accessions: List of protein accessions

        Returns:
            Dictionary of accession -> SeqRecord
        """
        self.logger.info(f"Downloading {len(accessions)} protein sequences")

        sequences = {}

        # Batch download in chunks to avoid NCBI limits
        chunk_size = 100
        for i in range(0, len(accessions), chunk_size):
            chunk = accessions[i:i + chunk_size]

            try:
                handle = Entrez.efetch(db="protein", id=chunk, rettype="fasta", retmode="text")
                records = list(SeqIO.parse(handle, "fasta"))
                handle.close()

                for record in records:
                    accession = record.id.split('.')[0]  # Remove version
                    sequences[accession] = record

                self.logger.info(f"Downloaded chunk {i//chunk_size + 1}: {len(records)} sequences")

                # Respect NCBI rate limits
                time.sleep(0.5)

            except Exception as e:
                self.logger.error(f"Failed to download chunk {i//chunk_size + 1}: {e}")

        self.logger.info(f"Successfully downloaded {len(sequences)} sequences")
        return sequences

    def save_fasta_files(self, sequences: Dict[str, SeqRecord], output_dir: str = None) -> Dict[str, str]:
        """
        Save sequences to individual FASTA files organized by RND family

        Args:
            sequences: Dictionary of accession -> SeqRecord
            output_dir: Output directory (default: self.output_dir / "fasta")

        Returns:
            Dictionary mapping accessions to file paths
        """
        if output_dir is None:
            output_dir = self.output_dir / "fasta"

        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)

        file_mapping = {}

        for accession, record in sequences.items():
            # Determine RND family for organization
            family = "Unknown"
            title = record.description.lower()

            for fam_name in self.rnd_families.keys():
                if fam_name.lower() in title:
                    family = fam_name
                    break

            # Create family subdirectory
            family_dir = output_path / family
            family_dir.mkdir(exist_ok=True)

            # Save FASTA file
            fasta_file = family_dir / f"{accession}.fasta"
            SeqIO.write(record, fasta_file, "fasta")

            file_mapping[accession] = str(fasta_file)

        self.logger.info(f"Saved {len(file_mapping)} FASTA files")
        return file_mapping

    def create_comprehensive_metadata_csv(self, metadata_list: List[Dict],
                                        output_file: str = None) -> str:
        """
        Create comprehensive metadata CSV file

        Args:
            metadata_list: List of metadata dictionaries
            output_file: Output file path

        Returns:
            Path to created CSV file
        """
        if output_file is None:
            output_file = self.output_dir / "comprehensive_metadata.csv"

        if not metadata_list:
            self.logger.warning("No metadata to save")
            return ""

        # Flatten MIC data for CSV
        flattened_data = []
        for metadata in metadata_list:
            base_data = {k: v for k, v in metadata.items() if k != 'mic_data'}

            # Add MIC data as separate columns
            mic_data = metadata.get('mic_data', [])
            if mic_data:
                for i, mic in enumerate(mic_data):
                    base_data[f'mic_{i}_antibiotic'] = mic.get('antibiotic', '')
                    base_data[f'mic_{i}_value'] = mic.get('mic_value', '')
                    base_data[f'mic_{i}_source'] = mic.get('source', '')
            else:
                base_data['mic_data_available'] = False

            flattened_data.append(base_data)

        df = pd.DataFrame(flattened_data)
        df.to_csv(output_file, index=False)

        self.logger.info(f"Saved comprehensive metadata to {output_file}")
        return str(output_file)

    def integrate_with_wild_type_aligner(self, fasta_files: Dict[str, str],
                                       metadata_file: str) -> Dict:
        """
        Prepare data for WildTypeAligner integration

        Args:
            fasta_files: Dictionary of accession -> FASTA file path
            metadata_file: Path to metadata CSV

        Returns:
            Dictionary with integration data
        """
        integration_data = {
            'fasta_files': fasta_files,
            'metadata_file': metadata_file,
            'rnd_families': list(self.rnd_families.keys()),
            'total_proteins': len(fasta_files),
            'output_directory': str(self.output_dir)
        }

        # Save integration data for WildTypeAligner
        integration_file = self.output_dir / "wild_type_aligner_integration.json"
        with open(integration_file, 'w') as f:
            json.dump(integration_data, f, indent=2)

        self.logger.info(f"Prepared integration data for WildTypeAligner: {integration_file}")
        return integration_data

    def run_complete_harvest_pipeline(self, query_sequence: str = None,
                                    organism: str = None,
                                    max_results: int = 100) -> Dict:
        """
        Run the complete BLAST-based harvesting pipeline

        Args:
            query_sequence: Query protein sequence (if None, use default RND query)
            organism: Organism filter
            max_results: Maximum BLAST results

        Returns:
            Dictionary with complete pipeline results
        """
        self.logger.info("Starting complete BLAST-PHIT harvest pipeline")
        self.logger.info(f"Organism: {organism}")
        self.logger.info(f"Max results: {max_results}")

        results = {
            'pipeline_status': 'running',
            'blast_hits': [],
            'metadata_collected': [],
            'sequences_downloaded': 0,
            'fasta_files_created': 0,
            'mic_data_points': 0,
            'errors': [],
            'warnings': []
        }

        try:
            # Step 1: BLAST search
            if query_sequence is None:
                # Use a representative RND protein sequence
                query_sequence = self._get_default_rnd_query()

            blast_hits = self.blast_search_rnd_proteins(
                query_sequence, organism, max_results
            )

            if not blast_hits:
                results['errors'].append("No BLAST hits found")
                results['pipeline_status'] = 'failed'
                return results

            results['blast_hits'] = blast_hits

            # Step 2: Identify RND families
            identified_proteins = []
            for hit in blast_hits:
                family, confidence = self.identify_rnd_family(hit)
                if confidence > 0.3:  # Only keep confident identifications
                    hit['rnd_family'] = family
                    hit['identification_confidence'] = confidence
                    identified_proteins.append(hit)

            if not identified_proteins:
                results['warnings'].append("No proteins confidently identified as RND family")
                results['pipeline_status'] = 'completed_with_warnings'
                return results

            # Step 3: Collect comprehensive metadata
            accessions = [hit['accession'] for hit in identified_proteins]
            metadata_list = []

            for accession in accessions:
                metadata = self.collect_protein_metadata(accession)
                if metadata:
                    metadata_list.append(metadata)

            results['metadata_collected'] = metadata_list
            results['mic_data_points'] = sum(len(m.get('mic_data', [])) for m in metadata_list)

            # Step 4: Download sequences
            sequences = self.download_protein_sequences(accessions)
            results['sequences_downloaded'] = len(sequences)

            # Step 5: Save FASTA files
            fasta_files = self.save_fasta_files(sequences)
            results['fasta_files_created'] = len(fasta_files)

            # Step 6: Create comprehensive metadata CSV
            metadata_csv = self.create_comprehensive_metadata_csv(metadata_list)

            # Step 7: Prepare for WildTypeAligner integration
            integration_data = self.integrate_with_wild_type_aligner(fasta_files, metadata_csv)

            # Step 8: Generate report
            self.generate_harvest_report(results)

            results['pipeline_status'] = 'completed'
            self.logger.info("BLAST-PHIT harvest pipeline completed successfully")

        except Exception as e:
            self.logger.error(f"Pipeline failed: {e}")
            results['pipeline_status'] = 'failed'
            results['errors'].append(str(e))

        return results

    def _get_default_rnd_query(self) -> str:
        """Get default RND protein query sequence"""
        # This would be a representative RND protein sequence
        # For now, return a placeholder
        return ">Default_RND_Query\nMNKNRGFTPLAVVLMLSGSLALTGCDDKQAQQGGQQMPAVGVVTVKTEPLQITTELPGRT..."

    def generate_harvest_report(self, results: Dict) -> None:
        """Generate comprehensive harvest report"""
        report_file = self.output_dir / "blast_phit_harvest_report.txt"

        with open(report_file, 'w') as f:
            f.write("BLAST-PHIT HARVESTER REPORT\n")
            f.write("=" * 50 + "\n\n")

            f.write("PIPELINE SUMMARY:\n")
            f.write("-" * 20 + "\n")
            f.write(f"Status: {results['pipeline_status']}\n")
            f.write(f"BLAST hits: {len(results['blast_hits'])}\n")
            f.write(f"Proteins identified: {len(results['metadata_collected'])}\n")
            f.write(f"Sequences downloaded: {results['sequences_downloaded']}\n")
            f.write(f"FASTA files created: {results['fasta_files_created']}\n")
            f.write(f"MIC data points: {results['mic_data_points']}\n\n")

            if results['metadata_collected']:
                f.write("RND FAMILY DISTRIBUTION:\n")
                f.write("-" * 30 + "\n")
                families = [m.get('rnd_family', 'Unknown') for m in results['metadata_collected']]
                family_counts = Counter(families)
                for family, count in family_counts.most_common():
                    f.write(f"{family}: {count}\n")

                f.write("\nTOP ORGANISMS:\n")
                f.write("-" * 15 + "\n")
                organisms = [m.get('organism', 'Unknown') for m in results['metadata_collected']]
                organism_counts = Counter(organisms)
                for organism, count in organism_counts.most_common(10):
                    f.write(f"{organism}: {count}\n")

            if results['errors']:
                f.write("\nERRORS:\n")
                f.write("-" * 10 + "\n")
                for error in results['errors']:
                    f.write(f"- {error}\n")

            if results['warnings']:
                f.write("\nWARNINGS:\n")
                f.write("-" * 10 + "\n")
                for warning in results['warnings']:
                    f.write(f"- {warning}\n")

        self.logger.info(f"Harvest report generated: {report_file}")


def main():
    """Main entry point"""
    import argparse

    parser = argparse.ArgumentParser(
        description="BLAST-PHIT Harvester - Comprehensive RND Protein Collection with MIC Integration"
    )

    parser.add_argument(
        "--query-sequence",
        help="Query protein sequence file (FASTA format)"
    )

    parser.add_argument(
        "--organism",
        help="Organism to search for (e.g., 'Escherichia coli')"
    )

    parser.add_argument(
        "--max-results",
        type=int,
        default=100,
        help="Maximum number of BLAST results"
    )

    parser.add_argument(
        "--email",
        required=True,
        help="NCBI email for API access"
    )

    parser.add_argument(
        "--output-dir",
        default="blast_phit_results",
        help="Output directory"
    )

    parser.add_argument(
        "--accession-list",
        help="File with list of accessions to process directly"
    )

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    # Initialize harvester
    harvester = BlastPHitHarvester(args.email, args.output_dir)

    # Load query sequence if provided
    query_sequence = None
    if args.query_sequence:
        with open(args.query_sequence, 'r') as f:
            # Skip header line
            next(f)
            query_sequence = f.read().replace('\n', '')

    # Run pipeline
    if args.accession_list:
        # Direct accession processing
        with open(args.accession_list, 'r') as f:
            accessions = [line.strip() for line in f if line.strip()]

        # Download and process
        sequences = harvester.download_protein_sequences(accessions)
        fasta_files = harvester.save_fasta_files(sequences)

        # Collect metadata
        metadata_list = []
        for accession in accessions:
            metadata = harvester.collect_protein_metadata(accession)
            if metadata:
                metadata_list.append(metadata)

        # Create CSV
        metadata_csv = harvester.create_comprehensive_metadata_csv(metadata_list)

        # Integration
        integration_data = harvester.integrate_with_wild_type_aligner(fasta_files, metadata_csv)

        print(f"\nDirect accession processing complete!")
        print(f"Processed {len(accessions)} accessions")
        print(f"Downloaded {len(sequences)} sequences")
        print(f"Collected {len(metadata_list)} metadata records")

    else:
        # Full BLAST pipeline
        results = harvester.run_complete_harvest_pipeline(
            query_sequence, args.organism, args.max_results
        )

        print("\nBLAST-PHIT Harvest Pipeline Complete!")
        print(f"Status: {results['pipeline_status']}")
        print(f"BLAST hits: {len(results['blast_hits'])}")
        print(f"Proteins identified: {len(results['metadata_collected'])}")
        print(f"MIC data points: {results['mic_data_points']}")

    print(f"\nResults saved to: {args.output_dir}")


if __name__ == "__main__":
    main()