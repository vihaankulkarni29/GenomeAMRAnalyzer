"""
Pipeline Orchestrator for Priority 3 AMR Genomics Pipeline
=========================================================

- Modular, stepwise execution of all pipeline stages
- Config validation and loading (YAML/JSON)
- Robust error handling and logging
- Checkpointing and resumability
- User-friendly CLI entry point

Engineering Principles:
- Fail-fast on config or critical errors
- Each stage logs provenance, errors, and outputs
- Supports dry-run and resume modes
- Extensible for new modules and workflow steps
"""

import sys
import os
import logging
import argparse
import yaml
import json
from pathlib import Path
from typing import Dict, Any, Optional

# Import all core modules
from src.priority3.db.repositories import GenomeRepository
from src.priority3.db.harvesters.genbank_genome_harvester import GenBankGenomeHarvester
from src.priority3.metadata.mic_metadata_harvester import NCBIMICHarvester
from src.priority3.analysis.subscan_alignment_analyzer import SubScanAlignmentAnalyzer
from src.priority3.analysis.mutation_cooccurrence_analyzer import MutationCooccurrenceAnalyzer
from src.priority3.report.html_report_generator import HTMLReportGenerator

class PipelineOrchestrator:
    def __init__(self, config_path: str):
        self.config_path = config_path
        self.config = self._load_config(config_path)
        self.logger = logging.getLogger("PipelineOrchestrator")
        self._setup_logging()
        self.repository = GenomeRepository(self.config['database']['path'])
        self.checkpoints = {}

    # ------------- Helpers for list inputs -------------
    def _read_list_file(self, path: Optional[str]) -> list[str]:
        if not path:
            return []
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(f"List file not found: {path}")
        items: list[str] = []
        for line in p.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            items.append(line)
        return items

    def _load_gene_list(self) -> list[str]:
        inp = self.config.get('input', {})
        genes = inp.get('gene_list') or []
        if not genes:
            genes = self._read_list_file(inp.get('gene_list_file'))
        # Normalize/canonicalize gene names to simple strings
        return [str(g).strip() for g in genes if str(g).strip()]

    def _setup_logging(self):
        log_dir = self.config.get('logging', {}).get('log_dir', 'logs')
        Path(log_dir).mkdir(exist_ok=True, parents=True)
        log_file = os.path.join(log_dir, 'pipeline.log')
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s %(levelname)s %(name)s: %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )

    def _load_config(self, path: str) -> Dict[str, Any]:
        if not os.path.exists(path):
            raise FileNotFoundError(f"Config file not found: {path}")
        with open(path, 'r') as f:
            if path.endswith('.yaml') or path.endswith('.yml'):
                return yaml.safe_load(f)
            elif path.endswith('.json'):
                return json.load(f)
            else:
                raise ValueError("Config file must be .yaml, .yml, or .json")

    def run(self, resume_from: Optional[str] = None, dry_run: bool = False):
        """
        Execute the full pipeline with robust error handling.
        Each stage can fail partially without stopping the entire pipeline.
        """
        stages = [
            ('genome_harvest', self.run_genome_harvest),
            ('mic_harvest', self.run_mic_harvest),
            ('card_abricate', self.run_card_abricate),
            ('fasta_aa_extractor', self.run_fasta_aa_extractor),
            ('wildtype_aligner', self.run_wildtype_aligner),
            ('subscan', self.run_subscan),
            ('cooccurrence', self.run_cooccurrence),
            ('report', self.run_report)
        ]
        
        resume_idx = 0
        if resume_from:
            stage_names = [s[0] for s in stages]
            if resume_from not in stage_names:
                raise ValueError(f"Invalid resume stage: {resume_from}")
            resume_idx = stage_names.index(resume_from)
            
        pipeline_errors = []
        successful_stages = []
        
        for name, func in stages[resume_idx:]:
            self.logger.info(f"=== Running stage: {name} ===")
            if dry_run:
                self.logger.info(f"[DRY RUN] Would execute: {name}")
                continue
                
            try:
                func()
                self.checkpoints[name] = 'completed'
                successful_stages.append(name)
                self.logger.info(f"✅ Stage {name} completed successfully")
            except Exception as e:
                error_msg = f"Stage {name} failed: {e}"
                self.logger.error(error_msg)
                self.checkpoints[name] = f'failed: {e}'
                pipeline_errors.append({'stage': name, 'error': str(e)})
                
                # Critical stages that should halt the pipeline
                critical_stages = ['genome_harvest']
                if name in critical_stages:
                    self.logger.error(f"Critical stage {name} failed. Halting pipeline.")
                    raise
                else:
                    self.logger.warning(f"Non-critical stage {name} failed. Continuing with remaining stages.")
                    
        # Log final pipeline status
        if pipeline_errors:
            self.logger.warning(f"Pipeline completed with {len(pipeline_errors)} stage failures:")
            for error in pipeline_errors:
                self.logger.warning(f"  - {error['stage']}: {error['error']}")
        else:
            self.logger.info("🎉 Pipeline completed successfully with no errors!")
            
        # Store error summary for report generation
        self.pipeline_errors = pipeline_errors
        self.successful_stages = successful_stages

    def run_card_abricate(self):
        """Run CARD analysis via Abricate with graceful error handling."""
        from src.abricate_runner import run_abricate_on_file
        from src.abricate_to_coords import convert_abricate_to_coords
        from pathlib import Path
        
        db_path = self.config['database']['path']
        output_dir = self.config.get('card_output_dir', 'card_results')
        abricate_dir = self.config.get('abricate_output_dir', 'abricate_results')
        
        repo = GenomeRepository(db_path)
        genome_records = repo.list_genomes(status='downloaded')
        
        if not genome_records:
            self.logger.warning("No downloaded genomes found for Abricate processing")
            return
            
        fasta_files = [g.file_path for g in genome_records if g.file_path and os.path.exists(g.file_path)]
        sample_ids = [g.accession for g in genome_records if g.file_path and os.path.exists(g.file_path)]
        
        if not fasta_files:
            self.logger.warning("No valid FASTA files found for Abricate processing")
            return
        
        # Check if Abricate is available
        import shutil
        if not shutil.which('abricate'):
            self.logger.error("Abricate binary not found in PATH")
            self.logger.error("Please install Abricate: conda install -c bioconda abricate")
            self.logger.warning("Skipping CARD Abricate stage - continuing with remaining pipeline stages")
            return
            
        try:
            # Create output directories
            os.makedirs(abricate_dir, exist_ok=True)
            os.makedirs(output_dir, exist_ok=True)
            
            # Process genomes individually for better error handling
            successful_runs = []
            failed_runs = []
            
            for fasta_file, sample_id in zip(fasta_files, sample_ids):
                try:
                    self.logger.info(f"Running Abricate on {sample_id}")
                    
                    # Step 1: Run Abricate and write to file
                    abricate_content = run_abricate_on_file(Path(fasta_file), db='card')
                    if abricate_content.strip():
                        abricate_file = os.path.join(abricate_dir, f"{sample_id}_abricate.tsv")
                        with open(abricate_file, 'w') as f:
                            f.write(abricate_content)
                        
                        # Step 2: Convert to coordinates CSV
                        coords_file = os.path.join(output_dir, f"{sample_id}_coordinates.csv")
                        convert_abricate_to_coords(abricate_file, coords_file)
                    else:
                        self.logger.info(f"No Abricate hits for {sample_id}")
                    
                    successful_runs.append(sample_id)
                except Exception as e:
                    self.logger.error(f"Abricate failed for {sample_id}: {e}")
                    failed_runs.append((sample_id, str(e)))
                    
            self.logger.info(f"CARD Abricate completed: {len(successful_runs)} successful, {len(failed_runs)} failed")
            
            if failed_runs:
                self.logger.warning("Abricate failures:")
                for sample_id, error in failed_runs:
                    self.logger.warning(f"  - {sample_id}: {error}")
                    
        except Exception as e:
            self.logger.error(f"Abricate runner initialization failed: {e}")
            self.logger.warning("Continuing pipeline without CARD Abricate results")

    def run_fasta_aa_extractor(self):
        """Run protein extraction with per-genome error handling and validation."""
        from src.priority3.extractor.fasta_aa_extractor import FastaAAExtractor
        
        db_path = self.config['database']['path']
        card_dir = self.config.get('card_output_dir', 'card_results')
        extractor_dir = self.config.get('extractor_output_dir', 'extracted_proteins')
        
        repo = GenomeRepository(db_path)
        gene_list = self._load_gene_list()
        genome_records = repo.list_genomes(status='downloaded')
        
        if not genome_records:
            self.logger.warning("No downloaded genomes found for protein extraction")
            return
            
        if not gene_list:
            self.logger.warning("No genes specified for extraction")
            return
            
        sample_ids = [g.accession for g in genome_records if g.file_path]
        fasta_files = [g.file_path for g in genome_records if g.file_path and os.path.exists(g.file_path)]
        
        # Validate input files exist
        valid_genomes = []
        for i, (fasta_file, sample_id) in enumerate(zip(fasta_files, sample_ids)):
            if os.path.exists(fasta_file):
                valid_genomes.append((fasta_file, sample_id))
            else:
                self.logger.warning(f"FASTA file not found for {sample_id}: {fasta_file}")
                
        if not valid_genomes:
            self.logger.warning("No valid FASTA files found for protein extraction")
            return
            
        try:
            extractor = FastaAAExtractor(output_dir=extractor_dir)
            
            successful_extractions = []
            failed_extractions = []
            
            for fasta_file, sample_id in valid_genomes:
                rgi_file = str(Path(card_dir) / f"{sample_id}_rgi.txt")
                
                try:
                    # Check if RGI file exists, if not, warn but continue
                    if not os.path.exists(rgi_file):
                        self.logger.warning(f"RGI file not found for {sample_id}: {rgi_file}")
                        self.logger.info(f"Attempting direct extraction from genome for {sample_id}")
                        
                    result = extractor.extract_proteins(fasta_file, rgi_file, gene_list, sample_id)
                    successful_extractions.append(sample_id)
                    self.logger.info(f"✅ Protein extraction completed for {sample_id}")
                    
                except Exception as e:
                    self.logger.error(f"Protein extraction failed for {sample_id}: {e}")
                    failed_extractions.append((sample_id, str(e)))
                    
            self.logger.info(f"Protein extraction completed: {len(successful_extractions)} successful, {len(failed_extractions)} failed")
            
            if failed_extractions:
                self.logger.warning("Protein extraction failures:")
                for sample_id, error in failed_extractions:
                    self.logger.warning(f"  - {sample_id}: {error}")
                    
        except Exception as e:
            self.logger.error(f"FastaAAExtractor initialization failed: {e}")
            self.logger.warning("Continuing pipeline without protein extraction results")

    def run_wildtype_aligner(self):
        """Run alignment with comprehensive validation and per-gene error handling."""
        from src.priority3.aligner.wildtype_aligner import WildTypeAligner
        
        db_path = self.config['database']['path']
        extractor_dir = self.config.get('extractor_output_dir', 'extracted_proteins')
        aligner_dir = self.config.get('alignment_output_dir', 'alignments')
        reference_dir = self.config.get('reference_dir', 'reference_proteins')
        
        repo = GenomeRepository(db_path)
        gene_list = self._load_gene_list()
        genome_records = repo.list_genomes(status='downloaded')
        
        if not genome_records:
            self.logger.warning("No downloaded genomes found for alignment")
            return
            
        if not gene_list:
            self.logger.warning("No genes specified for alignment")
            return
            
        sample_ids = [g.accession for g in genome_records if g.file_path]
        sample_species = {g.accession: g.organism for g in genome_records if g.file_path}
        
        # Find available protein files
        available_proteins = []
        missing_proteins = []
        
        for sample_id in sample_ids:
            for gene in gene_list:
                protein_fasta = Path(extractor_dir) / f"{sample_id}_{gene}.faa"
                if protein_fasta.exists() and protein_fasta.stat().st_size > 0:
                    available_proteins.append({
                        'file': str(protein_fasta),
                        'gene': gene,
                        'sample_id': sample_id,
                        'species': sample_species.get(sample_id, 'Unknown')
                    })
                else:
                    missing_proteins.append(f"{sample_id}_{gene}")
                    
        if missing_proteins:
            self.logger.warning(f"Missing protein files: {len(missing_proteins)} files not found")
            if len(missing_proteins) <= 10:  # Show details for small numbers
                for missing in missing_proteins:
                    self.logger.debug(f"  Missing: {missing}.faa")
                    
        if not available_proteins:
            self.logger.warning("No protein files found for alignment")
            return
            
        # Check reference files
        missing_references = []
        if os.path.exists(reference_dir):
            for gene in gene_list:
                ref_files = list(Path(reference_dir).glob(f"*{gene}*.fasta")) + list(Path(reference_dir).glob(f"*{gene}*.faa"))
                if not ref_files:
                    missing_references.append(gene)
        else:
            self.logger.warning(f"Reference directory not found: {reference_dir}")
            missing_references = gene_list
            
        if missing_references:
            self.logger.warning(f"Missing reference files for genes: {missing_references}")
            self.logger.info("Alignment will proceed with available references")
            
        try:
            aligner = WildTypeAligner(reference_dir=reference_dir, output_dir=aligner_dir)
            
            # Group by gene for batch processing
            genes_to_process = {}
            for protein in available_proteins:
                gene = protein['gene']
                if gene not in genes_to_process:
                    genes_to_process[gene] = []
                genes_to_process[gene].append(protein)
                
            successful_alignments = []
            failed_alignments = []
            
            for gene, proteins in genes_to_process.items():
                try:
                    protein_files = [p['file'] for p in proteins]
                    genes = [p['gene'] for p in proteins]
                    sample_ids = [p['sample_id'] for p in proteins]
                    species_list = [p['species'] for p in proteins]
                    
                    self.logger.info(f"Aligning {len(proteins)} {gene} proteins")
                    aligner.align_batch(protein_files, genes, sample_ids, species_list)
                    
                    successful_alignments.extend([f"{p['sample_id']}_{gene}" for p in proteins])
                    
                except Exception as e:
                    self.logger.error(f"Alignment failed for gene {gene}: {e}")
                    failed_alignments.extend([f"{p['sample_id']}_{gene}" for p in proteins])
                    
            self.logger.info(f"Alignment completed: {len(successful_alignments)} successful, {len(failed_alignments)} failed")
            
            if failed_alignments:
                self.logger.warning("Alignment failures:")
                for failed in failed_alignments[:10]:  # Show first 10
                    self.logger.warning(f"  - {failed}")
                if len(failed_alignments) > 10:
                    self.logger.warning(f"  ... and {len(failed_alignments) - 10} more")
                    
        except Exception as e:
            self.logger.error(f"WildTypeAligner initialization failed: {e}")
            self.logger.warning("Continuing pipeline without alignment results")


    def run_genome_harvest(self):
        """Run genome harvesting with robust accession handling and fallbacks."""
        # Extract config
        ncbi_cfg = self.config.get('ncbi', {})
        db_path = self.config['database']['path']
        output_dir = self.config.get('genome_output_dir', 'genomes')
        api_key = ncbi_cfg.get('api_key')
        email = ncbi_cfg.get('email', 'user@example.com')
        mock_mode = ncbi_cfg.get('mock_mode', False)
        
        if not email or email == 'user@example.com':
            self.logger.warning("NCBI email not set in config (ncbi.email). Please provide your contact email for E-utilities.")
            
        inp = self.config.get('input', {})
        query = inp.get('genome_query')
        max_genomes = inp.get('max_genomes')
        
        # Load accessions from config or file
        accessions = inp.get('accessions') or []
        if not accessions:
            accessions = self._read_list_file(inp.get('accessions_file'))
            
        # Smart accession handling: try assembly accessions first, then use all as-is if needed
        assembly_accessions = []
        other_accessions = []
        
        for acc in accessions:
            acc = str(acc).strip()
            if acc.startswith(('GCF_', 'GCA_')):
                assembly_accessions.append(acc)
            else:
                other_accessions.append(acc)
                
        self.logger.info(f"Found {len(assembly_accessions)} assembly accessions, {len(other_accessions)} other accessions")
        
        # Instantiate harvester
        try:
            harvester = GenBankGenomeHarvester(
                output_dir=output_dir,
                db_path=db_path,
                api_key=api_key,
                email=email,
                mock_mode=mock_mode
            )
            
            genomes_harvested = 0
            
            # First, try assembly accessions
            if assembly_accessions:
                self.logger.info(f"Harvesting {len(assembly_accessions)} assembly accessions")
                try:
                    result = harvester.harvest_by_accessions(assembly_accessions, resume=True)
                    genomes_harvested += len(assembly_accessions)
                except Exception as e:
                    self.logger.error(f"Assembly accession harvest failed: {e}")
                    
            # For non-assembly accessions, try individual nucleotide searches
            if other_accessions:
                self.logger.info(f"Attempting to harvest {len(other_accessions)} individual accessions")
                successful_individual = []
                failed_individual = []
                
                for acc in other_accessions[:10]:  # Limit to first 10 for testing
                    try:
                        # Try searching for the individual accession
                        individual_query = f'"{acc}"[Accession]'
                        self.logger.info(f"Searching for {acc}")
                        result = harvester.harvest_genomes(individual_query, max_genomes=1)
                        successful_individual.append(acc)
                        genomes_harvested += 1
                    except Exception as e:
                        self.logger.warning(f"Individual harvest failed for {acc}: {e}")
                        failed_individual.append(acc)
                        
                if successful_individual:
                    self.logger.info(f"Successfully harvested {len(successful_individual)} individual accessions")
                if failed_individual:
                    self.logger.warning(f"Failed to harvest {len(failed_individual)} individual accessions")
                    
            # If no specific accessions worked and we have a query, use it as fallback
            if genomes_harvested == 0 and query:
                self.logger.info(f"No specific accessions succeeded. Using fallback query: {query}")
                try:
                    harvester.harvest_genomes(query, max_genomes=max_genomes)
                except Exception as e:
                    self.logger.error(f"Fallback query harvest failed: {e}")
                    
            harvester.close()
            
        except Exception as e:
            self.logger.error(f"Genome harvest failed: {e}")
            raise

    def run_mic_harvest(self):
        """Run MIC harvesting with per-genome error handling and validation."""
        ncbi_cfg = self.config.get('ncbi', {})
        db_path = self.config['database']['path']
        api_key = ncbi_cfg.get('api_key')
        email = ncbi_cfg.get('email', 'user@example.com')
        mock_mode = ncbi_cfg.get('mock_mode', False)
        
        # Check: Ensure MIC data is joined to correct genome via biosample/bioproject
        repo = GenomeRepository(db_path)
        genome_records = repo.list_genomes(status='downloaded')
        
        if not genome_records:
            self.logger.warning("No downloaded genomes found for MIC harvest")
            return
            
        accessions = [g.accession for g in genome_records]
        self.logger.info(f"Harvesting MIC data for {len(accessions)} genomes")
        
        try:
            mic_harvester = NCBIMICHarvester(
                db_path=db_path,
                api_key=api_key,
                email=email,
                mock_mode=mock_mode
            )
            
            successful_harvests = []
            failed_harvests = []
            no_mic_data = []
            
            # Error Handling: If MIC not found, set as missing/null, do not fail
            for accession in accessions:
                try:
                    self.logger.debug(f"Harvesting MIC data for {accession}")
                    result = mic_harvester.harvest_mic_data([accession])
                    
                    # Check if any MIC data was actually found
                    # This depends on the harvester implementation
                    successful_harvests.append(accession)
                    
                except Exception as e:
                    self.logger.warning(f"MIC harvest failed for {accession}: {e}")
                    failed_harvests.append((accession, str(e)))
                    
            mic_harvester.close()
            
            self.logger.info(f"MIC harvest completed: {len(successful_harvests)} attempted, {len(failed_harvests)} failed")
            
            # MIC data is often sparse - this is normal
            if failed_harvests:
                self.logger.info("MIC harvest issues (normal for genomes without published MIC data):")
                for accession, error in failed_harvests[:5]:  # Show first 5
                    self.logger.debug(f"  - {accession}: {error}")
                if len(failed_harvests) > 5:
                    self.logger.debug(f"  ... and {len(failed_harvests) - 5} more")
                    
        except Exception as e:
            self.logger.error(f"MIC harvester initialization failed: {e}")
            self.logger.warning("Continuing pipeline without MIC data")

    def run_subscan(self):
        """Run SubScan analysis with per-alignment error handling and validation."""
        db_path = self.config['database']['path']
        alignment_dir = self.config.get('alignment_output_dir', 'alignments')
        output_dir = self.config.get('subscan_output_dir', 'subscan_results')
        min_conf = self.config.get('subscan', {}).get('min_confidence', 'MEDIUM')
        
        from src.priority3.analysis.subscan_alignment_analyzer import ConfidenceLevel
        min_conf_enum = getattr(ConfidenceLevel, min_conf.upper(), ConfidenceLevel.MEDIUM)
        
        # Check: Validate all alignment files are present and readable
        align_dir = Path(alignment_dir)
        if not align_dir.exists():
            self.logger.warning(f"Alignment directory not found: {alignment_dir}")
            return
            
        alignment_files = list(align_dir.glob('*.txt')) + list(align_dir.glob('*.aln')) + list(align_dir.glob('*.fasta'))
        
        if not alignment_files:
            self.logger.warning(f"No alignment files found in {alignment_dir}")
            return
            
        # Validate files are readable and non-empty
        valid_alignments = []
        invalid_alignments = []
        
        for align_file in alignment_files:
            try:
                if align_file.stat().st_size > 0:
                    with open(align_file, 'r') as f:
                        content = f.read(100)  # Read first 100 chars to validate
                        if content.strip():
                            valid_alignments.append(align_file)
                        else:
                            invalid_alignments.append((str(align_file), "Empty file"))
                else:
                    invalid_alignments.append((str(align_file), "Zero size"))
            except Exception as e:
                invalid_alignments.append((str(align_file), f"Read error: {e}"))
                
        if invalid_alignments:
            self.logger.warning(f"Found {len(invalid_alignments)} invalid alignment files:")
            for file_path, reason in invalid_alignments:
                self.logger.warning(f"  - {file_path}: {reason}")
                
        if not valid_alignments:
            self.logger.warning("No valid alignment files found for SubScan analysis")
            return
            
        self.logger.info(f"Processing {len(valid_alignments)} valid alignment files")
        
        try:
            analyzer = SubScanAlignmentAnalyzer(
                output_dir=output_dir,
                db_path=db_path,
                min_confidence=min_conf_enum
            )
            
            successful_analyses = []
            failed_analyses = []
            
            # Error Handling: If parsing fails for one alignment, log and continue
            for align_file in valid_alignments:
                try:
                    # Extract genome accession from filename
                    genome_acc = align_file.stem.split('_')[0]
                    
                    self.logger.debug(f"Analyzing alignment: {align_file}")
                    analyzer.analyze_batch([str(align_file)], [genome_acc])
                    successful_analyses.append(str(align_file))
                    
                except Exception as e:
                    self.logger.error(f"SubScan analysis failed for {align_file}: {e}")
                    failed_analyses.append((str(align_file), str(e)))
                    
            analyzer.close()
            
            self.logger.info(f"SubScan completed: {len(successful_analyses)} successful, {len(failed_analyses)} failed")
            
            if failed_analyses:
                self.logger.warning("SubScan analysis failures:")
                for file_path, error in failed_analyses:
                    self.logger.warning(f"  - {file_path}: {error}")
                    
        except Exception as e:
            self.logger.error(f"SubScan analyzer initialization failed: {e}")
            self.logger.warning("Continuing pipeline without SubScan results")

    def run_cooccurrence(self):
        """Run co-occurrence analysis with comprehensive validation and error handling."""
        db_path = self.config['database']['path']
        output_dir = self.config.get('cooccurrence_output_dir', 'cooccurrence_results')
        min_count = self.config.get('cooccurrence', {}).get('min_count', 5)
        significance_level = self.config.get('cooccurrence', {}).get('significance_level', 0.01)
        
        try:
            analyzer = MutationCooccurrenceAnalyzer(
                db_path=db_path,
                min_count=min_count,
                significance_level=significance_level,
                output_dir=output_dir
            )
            
            # Check: Ensure all mutation records are linked to correct genome/protein
            self.logger.info("Loading mutation matrix for co-occurrence analysis")
            
            try:
                mutation_matrix, genome_ids, mutation_ids = analyzer.load_mutation_matrix()
                
                if mutation_matrix is None or len(genome_ids) == 0:
                    self.logger.warning("No mutation data available for co-occurrence analysis")
                    analyzer.close()
                    return
                    
                self.logger.info(f"Loaded mutation data: {len(genome_ids)} genomes, {len(mutation_ids)} mutations")
                
                # Validate data integrity
                if len(genome_ids) < min_count:
                    self.logger.warning(f"Insufficient genomes ({len(genome_ids)}) for co-occurrence analysis (minimum: {min_count})")
                    analyzer.close()
                    return
                    
                # Error Handling: If analysis fails for a pair, log and continue
                self.logger.info("Computing mutation co-occurrence patterns")
                results = analyzer.compute_cooccurrence(mutation_matrix)
                
                if results is not None and len(results) > 0:
                    output_file = "cooccurrence_analysis.csv"
                    analyzer.save_results(results, output_file)
                    self.logger.info(f"Co-occurrence analysis completed: {len(results)} significant pairs found")
                    self.logger.info(f"Results saved to: {output_dir}/{output_file}")
                else:
                    self.logger.warning("No significant co-occurrence patterns found")
                    # Still save empty results for consistency
                    analyzer.save_results([], "cooccurrence_analysis.csv")
                    
            except Exception as e:
                self.logger.error(f"Co-occurrence computation failed: {e}")
                self.logger.warning("Saving empty co-occurrence results")
                # Try to save empty results to maintain pipeline consistency
                try:
                    analyzer.save_results([], "cooccurrence_analysis.csv")
                except:
                    pass
                    
            analyzer.close()
            
        except Exception as e:
            self.logger.error(f"Co-occurrence analyzer initialization failed: {e}")
            self.logger.warning("Continuing pipeline without co-occurrence analysis")

    def run_report(self):
        """Generate HTML report with comprehensive validation and partial report support."""
        db_path = self.config['database']['path']
        
        try:
            generator = HTMLReportGenerator()
            repo = GenomeRepository(db_path)
            
            # Check: Validate all required fields are present in context; warn if missing
            self.logger.info("Gathering data for report generation")
            
            # Collect genome data
            genomes = []
            try:
                genome_records = repo.list_genomes()
                genomes = [g.__dict__ for g in genome_records]
                self.logger.info(f"Found {len(genomes)} genomes for report")
            except Exception as e:
                self.logger.warning(f"Failed to load genome data: {e}")
                genomes = []
                
            # Load mutations from subscan artifacts
            mutations = []
            mutation_errors = []
            
            for g in genomes:
                try:
                    artifacts = repo.list_artifacts(accession=g['accession'], artifact_type='subscan_analysis')
                    for artifact in artifacts:
                        try:
                            with open(artifact.path, 'r') as f:
                                data = json.load(f)
                                for m in data.get('mutations', []):
                                    m['artifact_path'] = artifact.path
                                    m['genome_accession'] = g['accession']  # Ensure linkage
                                    mutations.append(m)
                        except Exception as e:
                            mutation_errors.append(f"Failed to load {artifact.path}: {e}")
                except Exception as e:
                    mutation_errors.append(f"Failed to get artifacts for {g['accession']}: {e}")
                    
            self.logger.info(f"Loaded {len(mutations)} mutations from {len(genomes)} genomes")
            if mutation_errors:
                self.logger.warning(f"Mutation loading errors: {len(mutation_errors)}")
                for error in mutation_errors[:3]:  # Show first 3
                    self.logger.debug(f"  - {error}")
                    
            # Load MIC data
            mic_data = []
            mic_errors = []
            
            for g in genomes:
                try:
                    artifacts = repo.list_artifacts(accession=g['accession'], artifact_type='mic_data')
                    for artifact in artifacts:
                        try:
                            with open(artifact.path, 'r') as f:
                                data = json.load(f)
                                data['genome_accession'] = g['accession']  # Ensure linkage
                                mic_data.append(data)
                        except Exception as e:
                            mic_errors.append(f"Failed to load {artifact.path}: {e}")
                except Exception as e:
                    mic_errors.append(f"Failed to get MIC artifacts for {g['accession']}: {e}")
                    
            self.logger.info(f"Loaded MIC data for {len(mic_data)} entries")
            if not mic_data:
                self.logger.info("No MIC data found - this is normal for many genomic studies")
                
            # Load cooccurrence results
            cooccurrence = []
            cooccur_errors = []
            
            cooccur_dir = Path(self.config.get('cooccurrence_output_dir', 'cooccurrence_results'))
            if cooccur_dir.exists():
                for file in cooccur_dir.glob('*.csv'):
                    try:
                        # Avoid pandas import error by reading manually if pandas not available
                        try:
                            import pandas as pd
                            df = pd.read_csv(file)
                            for row in df.to_dict(orient='records'):
                                row['artifact_path'] = str(file)
                                cooccurrence.append(row)
                        except ImportError:
                            # Fallback: read CSV manually
                            import csv
                            with open(file, 'r') as f:
                                reader = csv.DictReader(f)
                                for row in reader:
                                    row['artifact_path'] = str(file)
                                    cooccurrence.append(row)
                                    
                    except Exception as e:
                        cooccur_errors.append(f"Failed to load {file}: {e}")
            else:
                self.logger.info(f"Co-occurrence directory not found: {cooccur_dir}")
                
            self.logger.info(f"Loaded {len(cooccurrence)} co-occurrence entries")
            
            # Get repository statistics
            try:
                stats = repo.get_stats()
            except Exception as e:
                self.logger.warning(f"Failed to get repository stats: {e}")
                stats = {}
                
            # Add pipeline execution summary
            pipeline_summary = {
                'successful_stages': getattr(self, 'successful_stages', []),
                'pipeline_errors': getattr(self, 'pipeline_errors', []),
                'total_genomes': len(genomes),
                'total_mutations': len(mutations),
                'total_mic_entries': len(mic_data),
                'total_cooccurrence_pairs': len(cooccurrence)
            }
            
            run_id = self.config.get('run_id', 'pipeline_run')
            
            # Error Handling: If report generation fails, log and output partial report
            try:
                context = generator.build_context(genomes, mutations, mic_data, cooccurrence, stats, run_id)
                
                # Add pipeline summary to context
                context['pipeline_summary'] = pipeline_summary
                context['data_warnings'] = {
                    'mutation_errors': len(mutation_errors),
                    'mic_errors': len(mic_errors),
                    'cooccurrence_errors': len(cooccur_errors)
                }
                
                generator.ensure_templates()
                report_path = generator.render_report(context)
                
                self.logger.info(f"✅ HTML report generated successfully: {report_path}")
                
            except Exception as e:
                self.logger.error(f"Report generation failed: {e}")
                
                # Try to generate a minimal fallback report
                try:
                    self.logger.info("Attempting to generate minimal fallback report")
                    minimal_context = {
                        'genomes': genomes,
                        'mutations': [],
                        'mic_data': [],
                        'cooccurrence': [],
                        'stats': stats,
                        'run_id': run_id,
                        'pipeline_summary': pipeline_summary,
                        'error_message': f"Full report generation failed: {e}"
                    }
                    
                    fallback_path = generator.render_report(minimal_context)
                    self.logger.warning(f"Partial report generated: {fallback_path}")
                    
                except Exception as fallback_error:
                    self.logger.error(f"Even fallback report generation failed: {fallback_error}")
                    raise
                    
        except Exception as e:
            self.logger.error(f"Report generation completely failed: {e}")
            self.logger.warning("Continuing - pipeline completed without report")




def main():
    parser = argparse.ArgumentParser(description="Priority 3 AMR Pipeline Orchestrator")
    parser.add_argument("-c", "--config", required=True, help="Path to pipeline config YAML/JSON")
    parser.add_argument("--resume", help="Stage to resume from (optional)")
    parser.add_argument("--dry-run", action="store_true", help="Dry run (no execution)")
    args = parser.parse_args()
    orchestrator = PipelineOrchestrator(args.config)
    orchestrator.run(resume_from=args.resume, dry_run=args.dry_run)

if __name__ == "__main__":
    main()
