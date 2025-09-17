import sys; import os; sys.path.append(os.path.join(os.path.dirname(__file__), '\\\src')); from configuration_manager import config_manager
#!/usr/bin/env python3
"""
Example: FastaAAExtractor Integration Workflow
==============================================

This example demonstrates how to use the FastaAAExtractor integration
to extract protein sequences from genome coordinates.
"""

import sys
import os
import tempfile
import shutil
from pathlib import Path

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))

from fasta_aa_extractor_integration import FastaAAExtractorIntegrator

def create_sample_files():
    """Create sample input files for demonstration"""
    
    # Create temporary directory
    temp_dir = Path(tempfile.mkdtemp(prefix="fasta_extractor_example_"))
    
    # Create sample CARD coordinates file
    coords_file = temp_dir / "sample_coordinates.csv"
    with open(coords_file, 'w') as f:
        f.write("genome_id,gene_name,start,end,strand,contig\n")
        f.write("GCF_000005825.2,mdtF,1245,2678,+,chromosome\n")
        f.write("GCF_000005825.2,acrA,2890,3825,+,chromosome\n")
        f.write("GCF_000005825.2,acrB,3840,6951,+,chromosome\n")
        f.write("GCF_000006945.2,mdtF,1180,2613,+,chromosome\n")
        f.write("GCF_000006945.2,acrA,2825,3760,+,chromosome\n")
        f.write("GCF_000006945.2,tolC,4100,5500,+,chromosome\n")
    
    # Create sample genome files (minimal FASTA)
    genomes_dir = temp_dir / "genomes"
    genomes_dir.mkdir()
    
    # Sample genome 1
    genome1_file = genomes_dir / "GCF_000005825.2_genomic.fna"
    with open(genome1_file, 'w') as f:
        f.write(">GCF_000005825.2_chromosome\n")
        f.write("ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTGAAAAAACCATTAGCGGCCAGGATGCTTTACCCAATATCAGCGATGCCGAACGTATTTTTGCCGAACTTTTGACGGGACTCGCCGCCGCCCAGCCGGGGTTCCCGCTGGCGCAATTGAAAACTTTCGTCGATCAGGAATTTGCCCAAATAAAACATGTCCTGCATGGCATTAGTTTGTTGGGGCAGTGCCCGGATAGCATCAACGCTGCGCTGATTTGCCGTGGCGAGAAAATGTCGATCGCCATTATGGCCGGCGTATTAGAAGCGCGCGGTCACAACGTTACTGTTATCGATCCGGTCGTTCTGAATGGCGGTTTCCGGGGCTGGCCGCATCATGGCGGCTTGGGTTGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTG" * 50 + "\n")  # Repeat to make it long enough
    
    # Sample genome 2
    genome2_file = genomes_dir / "GCF_000006945.2_genomic.fna"
    with open(genome2_file, 'w') as f:
        f.write(">GCF_000006945.2_chromosome\n")
        f.write("ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTGAAAAAACCATTAGCGGCCAGGATGCTTTACCCAATATCAGCGATGCCGAACGTATTTTTGCCGAACTTTTGACGGGACTCGCCGCCGCCCAGCCGGGGTTCCCGCTGGCGCAATTGAAAACTTTCGTCGATCAGGAATTTGCCCAAATAAAACATGTCCTGCATGGCATTAGTTTGTTGGGGCAGTGCCCGGATAGCATCAACGCTGCGCTGATTTGCCGTGGCGAGAAAATGTCGATCGCCATTATGGCCGGCGTATTAGAAGCGCGCGGTCACAACGTTACTGTTATCGATCCGGTCGTTCTGAATGGCGGTTTCCGGGGCTGGCCGCATCATGGCGGCTTGGGTTGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTG" * 50 + "\n")
    
    # Create sample reference files
    references_dir = temp_dir / "references"
    references_dir.mkdir()
    
    ref_sequences = {
        "mdtF": "MKRLAPPITTTITITTQVTVRADAYRKHRKPPPDSAGFFFFFFFFDQRVTVTACEPTFGKLSSANAMQNFLRCCDILESNARGRGAVATCLCCTKITLGLAIDEKTISGRMLYPISDPERISFELDRPAPAGFPLAQEKRHYLIMFPKKCSHGILFL",
        config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]: "MARAPLKTLLSLLQFRSLLLLGLLLLPSVREEEERERERMLRTNSLTLRKRTKRTKRTN",
        config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]: "MPLKTLLSLLQFRSLLLLGLLLLPSVREEEERERERMLRTNSLTLRKRTKRTKRTNERE", 
        config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2]: "MKTLLSLLQFRSLLLLGLLLLPSVREEEERERERMLRTNSLTLRKRTKRTKRTNERERE"
    }
    
    for gene, sequence in ref_sequences.items():
        ref_file = references_dir / f"{gene}_reference.faa"
        with open(ref_file, 'w') as f:
            f.write(f">{gene}_MG1655_reference\n")
            f.write(f"{sequence}\n")
    
    return temp_dir, coords_file, genomes_dir, references_dir

def main():
    """Run FastaAAExtractor integration example"""
    
    print("GenomeAMRAnalyzer - FastaAAExtractor Integration Example")
    print("=" * 65)
    
    # Create sample data
    print("\n1. Creating sample data files...")
    temp_dir, coords_file, genomes_dir, references_dir = create_sample_files()
    print(f"   Sample data created in: {temp_dir}")
    
    # Initialize output directory
    output_dir = temp_dir / "extraction_output"
    
    try:
        # Initialize integrator
        print(f"\n2. Initializing FastaAAExtractor integrator...")
        integrator = FastaAAExtractorIntegrator(
            output_dir=str(output_dir),
            log_level="INFO"
        )
        
        # Load CARD coordinates
        print(f"\n3. Loading CARD coordinates from {coords_file.name}...")
        integrator.load_card_coordinates(
            coordinates_file=coords_file,
            file_format="csv"
        )
        
        # Display coordinate summary
        total_coords = sum(len(coords) for coords in integrator.gene_coordinates.values())
        print(f"   Loaded coordinates for {len(integrator.gene_coordinates)} genomes")
        print(f"   Total gene coordinates: {total_coords}")
        
        print(f"\n   Genes found:")
        all_genes = set()
        for coords in integrator.gene_coordinates.values():
            all_genes.update(coord.gene_name for coord in coords)
        for gene in sorted(all_genes):
            print(f"     - {gene}")
        
        # Extract proteins (using internal method for this example)
        print(f"\n4. Extracting protein sequences...")
        print(f"   Using internal BioPython method...")
        
        output_files = integrator.extract_proteins(
            genomes_dir=genomes_dir,
            gene_list=["mdtF", config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1], config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2]],  # Filter to specific genes
            use_external_tool=False  # Use internal method
        )
        
        # Display extraction results
        print(f"\n5. Extraction Results:")
        print(f"   Genomes processed: {integrator.stats['total_genomes']}")
        print(f"   Successful extractions: {integrator.stats['successful_extractions']}")
        print(f"   Failed extractions: {integrator.stats['failed_extractions']}")
        print(f"   Sequences generated: {integrator.stats['sequences_generated']}")
        
        # Show output files
        print(f"\n6. Output Files Generated:")
        for file_type, file_path in output_files.items():
            if isinstance(file_path, list):
                print(f"   {file_type}: {len(file_path)} files")
                for fp in file_path[:3]:
                    file_name = Path(fp).name
                    print(f"     - {file_name}")
                if len(file_path) > 3:
                    print(f"     - ... and {len(file_path) - 3} more")
            else:
                file_name = Path(file_path).name
                print(f"   {file_type}: {file_name}")
        
        # Prepare for WildTypeAligner
        print(f"\n7. Preparing sequences for WildTypeAligner...")
        aligner_files = integrator.prepare_for_wild_type_aligner(
            reference_dir=references_dir,
            output_structure="by_gene"
        )
        
        print(f"   Prepared {len(aligner_files)} gene datasets for alignment:")
        for gene_name, files in aligner_files.items():
            print(f"     {gene_name}:")
            if isinstance(files, dict):
                query_file = files.get('query_sequences', '')
                ref_file = files.get('reference_sequence', '')
                output_dir_path = files.get('output_directory', '')
                
                if query_file:
                    print(f"       - Query sequences: {Path(query_file).name}")
                if ref_file:
                    print(f"       - Reference sequence: {Path(ref_file).name}")
                if output_dir_path:
                    print(f"       - Output directory: {Path(output_dir_path).name}")
            else:
                print(f"       - Files: {files}")
        
        # Show sample extracted sequence
        print(f"\n8. Sample Extracted Sequences:")
        if 'combined_fasta' in output_files:
            combined_file = Path(output_files['combined_fasta'])
            if combined_file.exists():
                with open(combined_file, 'r') as f:
                    lines = f.readlines()
                    if len(lines) >= 4:  # At least one complete sequence
                        print(f"   First sequence from {combined_file.name}:")
                        print(f"     Header: {lines[0].strip()}")
                        print(f"     Sequence: {lines[1][:60]}..." if len(lines[1]) > 60 else f"     Sequence: {lines[1].strip()}")
        
        print(f"\n9. Integration Example Complete!")
        print(f"   Next steps:")
        print(f"   - Run WildTypeAligner on prepared sequences")
        print(f"   - Use SubScan for mutation detection")
        print(f"   - Analyze results with Co-occurrence Analyzer")
        
        # Display file locations
        print(f"\n   All output files are in: {output_dir}")
        print(f"   WildTypeAligner input files are in: {output_dir}/wild_type_aligner_input/")
        
    except Exception as e:
        print(f"\nError during processing: {e}")
        import traceback
        traceback.print_exc()
        
    finally:
        # Cleanup
        print(f"\n10. Cleaning up temporary files...")
        try:
            shutil.rmtree(temp_dir)
            print(f"    Temporary directory cleaned up successfully")
        except Exception as e:
            print(f"    Warning: Could not clean up {temp_dir}: {e}")

if __name__ == "__main__":
    main()