"""
GenomeAMRAnalyzer: Production-grade antimicrobial resistance analysis pipeline
Enterprise-ready setup for bioinformatics workflows and bacterial genome analysis
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read the README file for long description
this_directory = Path(__file__).parent
try:
    long_description = (this_directory / "README.md").read_text(encoding='utf-8')
except FileNotFoundError:
    long_description = "Production-grade antimicrobial resistance gene analysis pipeline for bacterial genomes"

# Read requirements with error handling
requirements = []
try:
    with open('requirements.txt', 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#') and not line.startswith('-'):
                # Clean package names for setup.py compatibility
                package = line.split('>=')[0].split('==')[0].split('<')[0].split('[')[0]
                if package:
                    requirements.append(line)  # Keep full version specification
except FileNotFoundError:
    # Fallback requirements if file missing
    requirements = ["biopython>=1.80", "pandas>=1.5.0", "numpy>=1.24.0"]

setup(
    name="genomeamranalyzer",
    version="2.0.0",
    author="GenomeAMRAnalyzer Development Team",
    author_email="your.email@institution.edu",
    description="Production-grade antimicrobial resistance gene analysis pipeline for bacterial genomes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/GenomeAMRAnalyzer",
    project_urls={
        "Bug Reports": "https://github.com/yourusername/GenomeAMRAnalyzer/issues",
        "Source": "https://github.com/yourusername/GenomeAMRAnalyzer",
        "Documentation": "https://github.com/yourusername/GenomeAMRAnalyzer/wiki",
        "Changelog": "https://github.com/yourusername/GenomeAMRAnalyzer/blob/main/CHANGELOG.md",
    },
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Healthcare Industry",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Operating System :: OS Independent",
        "Environment :: Console",
        "Environment :: Web Environment",
        "Natural Language :: English",
    ],
    python_requires=">=3.9",
    install_requires=requirements,
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=4.0.0",
            "black>=23.0.0",
            "flake8>=6.0.0",
            "mypy>=1.0.0",
            "isort>=5.12.0",
        ],
        "performance": [
            "numba>=0.58.0",
            "dask>=2023.1.0",
            "joblib>=1.3.0",
        ],
        "visualization": [
            "matplotlib>=3.6.0",
            "seaborn>=0.12.0",
            "plotly>=5.0.0",
            "bokeh>=3.0.0",
            "holoviews>=1.17.0",
        ],
        "export": [
            "weasyprint>=59.0",
            "reportlab>=4.0.0",
        ],
        "docs": [
            "sphinx>=6.0",
            "sphinx-rtd-theme>=1.3",
            "myst-parser>=2.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "genomeamr=run_pipeline:main",
            "genomeamr-download=src.simple_genome_downloader:main",
            # RGI-based card runner removed; Abricate flow is handled via run_pipeline and url_to_card_pipeline
            "genomeamr-extract=src.fasta_aa_extractor_integration:main",
            "genomeamr-align=src.production_wildtype_aligner:main",
            "genomeamr-subscan=src.production_subscan_analyzer:main",
            "genomeamr-cooccur=src.production_cooccurrence_analyzer:main",
            "genomeamr-report=src.html_report_generator:main",
        ],
    },
    include_package_data=True,
    package_data={
        "": ["config/*.yaml", "config/*.json", "examples/*", "docs/*", "*.md"],
    },
    zip_safe=False,
    keywords=[
        "antimicrobial resistance", "AMR", "genomics", "bioinformatics", 
        "CARD", "RGI", "efflux pumps", "resistance genes", "cooccurrence",
        "network analysis", "bacterial genomes", "clinical microbiology"
    ],
    platforms=["any"],
    license="MIT",
)