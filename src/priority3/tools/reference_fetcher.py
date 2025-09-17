import os
import subprocess
import sys
from pathlib import Path
import sys; import os; sys.path.append(os.path.join(os.path.dirname(__file__), '\\\\\')); from .configuration_manager import config_manager

# Point to the actual SEPI location
SEPI_PATH = Path(__file__).parent.parent.parent.parent / "MetaDataHarvester" / "sepi2.0"


def ensure_sepi():
    """Ensure SEPI is available and accessible."""
    if not SEPI_PATH.exists():
        print("SEPI not found at expected location. Attempting to clone...")
        # Try to clone to the correct location
        subprocess.run([
            'git', 'clone', 'https://github.com/vihaankulkarni29/sepi2.0', str(SEPI_PATH)
        ], check=True)
    
    sepi_script = SEPI_PATH / 'sepi.py'
    if not sepi_script.exists():
        raise FileNotFoundError(f"SEPI script not found at {sepi_script}")


def fetch_reference_protein(genus, species, gene, output_fasta):
    """
    Fetches the reference protein using SEPI 2.0 for the given species and gene.
    
    Args:
        genus: Genus name (used for consistency, but SEPI uses full organism name)
        species: Full species name (e.g., "Escherichia coli")
        gene: Gene name (e.g., config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0])
        output_fasta: Output FASTA file path
        
    Returns:
        str: Path to the output FASTA file if successful
        
    Raises:
        RuntimeError: If SEPI fails
        FileNotFoundError: If SEPI doesn't produce output
    """
    ensure_sepi()
    sepi_script = SEPI_PATH / 'sepi.py'
    
    # SEPI 2.0 uses --organism and --proteins parameters
    output_base = Path(output_fasta).stem  # Remove .faa extension for SEPI base name
    
    cmd = [
        sys.executable, str(sepi_script),
        '--organism', species,  # Use full species name
        '--proteins', gene,     # Single gene/protein name
        '--output', output_base,
        '--assembly_level', 'complete_genome',  # Prefer complete genomes
        '--email', 'user@example.com'  # Required by NCBI
    ]
    
    # Set working directory to where we want the output
    work_dir = Path(output_fasta).parent
    
    try:
        result = subprocess.run(
            cmd, 
            capture_output=True, 
            text=True, 
            timeout=120,  # Longer timeout for network operations
            cwd=str(work_dir)
        )
        
        if result.returncode != 0:
            raise RuntimeError(f'SEPI failed (return code {result.returncode}): {result.stderr}')
        
        # SEPI 2.0 creates files with specific naming convention
        # It creates a zip file containing the FASTA files
        zip_file = work_dir / f"{output_base}.zip"
        
        if zip_file.exists():
            # Extract the zip file to find the FASTA
            import zipfile
            with zipfile.ZipFile(zip_file, 'r') as zf:
                zip_contents = zf.namelist()
                
                # Look for the gene-specific FASTA file
                gene_fasta = None
                for file_name in zip_contents:
                    if file_name.lower().startswith(gene.lower()) and file_name.endswith('.fasta'):
                        gene_fasta = file_name
                        break
                
                if gene_fasta:
                    # Extract the specific file
                    zf.extract(gene_fasta, work_dir)
                    extracted_file = work_dir / gene_fasta
                    
                    # Move to the expected output location
                    import shutil
                    shutil.move(str(extracted_file), str(output_fasta))
                    
                    return str(output_fasta)
        
        # If zip method didn't work, try other file patterns
        potential_files = [
            work_dir / f"{gene}_{output_base}.fasta",
            work_dir / f"{output_base}_{gene}.fasta",
            work_dir / f"{gene}.fasta"
        ]
        
        created_file = None
        for potential_file in potential_files:
            if potential_file.exists():
                created_file = potential_file
                break
        
        if created_file is None:
            # List files in work directory for debugging
            files_in_dir = list(work_dir.glob("*"))
            raise FileNotFoundError(
                f"SEPI completed but no FASTA file found. "
                f"Expected gene: {gene}, output: {output_fasta}, "
                f"Files in directory: {files_in_dir}"
            )
        
        # If the created file is not at the expected location, move it
        if created_file != Path(output_fasta):
            import shutil
            shutil.move(str(created_file), str(output_fasta))
        
        return str(output_fasta)
        
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"SEPI timed out after 120 seconds")
    except Exception as e:
        raise RuntimeError(f"SEPI execution error: {e}")
