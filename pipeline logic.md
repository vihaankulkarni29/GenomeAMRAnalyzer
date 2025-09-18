# GenomeAMRAnalyzer: Logical Workflow

This document explains the core logic behind each major module of the GenomeAMRAnalyzer pipeline, focusing on the "why" and "how" without deep technical jargon.

---

## 📥 User Input

The analysis begins when a user provides three key pieces of information:

1.  **Genome Source:** A URL from an NCBI search (e.g., "Escherichia coli AND Erythromycin resistance") or a simple text file listing specific genome accession numbers.
2.  **Gene List:** A list of genes the researcher is interested in studying for potential mutations (e.g., `acrB`, `acrA`, `mdtH`).
3.  **User Email:** An email address is required for programmatic access to NCBI's databases.

---

## 🧬 Pipeline Modules

### 1. NCBI Genome Extractor

**Purpose:** This is the starting block. Its job is to fetch the raw materials for the analysis—the genome sequences. Think of it as a specialized librarian for the vast NCBI database.

It uses the provided **NCBI search URL** or **accession list** to identify and download all the relevant genome files. Critically, it also gathers high-quality **metadata** for each genome, such as collection date, geographical origin, and the study it belongs to. This context is vital for later interpretation.

* **Input:** NCBI URL or accession list, user's email.
* **Output:** Genome sequence files (fastas) and a detailed metadata file.
* **For more details:** [NCBI Genome Extractor Repository](https://github.com/vihaankulkarni29/ncbi_genome_extractor)

*(Flow: User Input -> **NCBI Genome Extractor** -> CARD RGI Integration)*

---

### 2. CARD RGI Integration

**Purpose:** This module's primary job in this pipeline is to act as a precise gene mapper. While it can identify all AMR genes, we use it for a more clever purpose: to get the exact location (**coordinates**) of the genes the user wants to analyze.

The next step needs to know *where* in the massive genome file to find `acrB`, for example. CARD RGI scans the genome and creates a map, telling us, "The `acrB` gene starts at position X and ends at position Y." This map is essential for the next module.

* **Input:** Genome sequence files from the previous step.
* **Output:** Coordinate files that map the locations of all AMR genes in each genome.

*(Flow: NCBI Genome Extractor -> **CARD RGI Integration** -> FastaAAExtractor)*

---

### 3. FastaAAExtractor

**Purpose:** Genes are recipes for making proteins. This module acts as the chef, using the recipe (gene) to create the final dish (protein).

It takes the **genome sequences**, the **gene map (coordinates)** from CARD RGI, and the user's specific **gene list**. For each genome, it goes to the exact coordinates for a gene like `acrB`, "cuts out" that genetic sequence, and translates it into its corresponding protein sequence (an amino acid chain). Each resulting protein file is carefully named after the genome it came from to keep perfect records.

* **Input:** Genome files, coordinate files, and the user's gene list.
* **Output:** Protein sequence files (protein fastas) for each requested gene from each genome.

*(Flow: CARD RGI Integration -> **FastaAAExtractor** -> WildTypeAligner)*

---

### 4. WildTypeAligner

**Purpose:** This is the **heartbeat of the pipeline**. It's where we find out if the user's suspected genes are actually mutated. It works by comparing the extracted proteins against a "normal," non-resistant **reference protein**.

Imagine you have a sentence that might have a typo. To find it, you'd compare it letter-by-letter to the original, correct sentence. This module does the same for protein sequences. To do this, it needs a proper reference (e.g., the AcrB protein from a standard, non-resistant *E. coli* K-12 strain).

Instead of storing a huge library of all possible references, this module uses a smart tool called **SEPI** to fetch the perfect reference protein on-the-fly. This keeps the pipeline lightweight and dynamic. The final output is a detailed, side-by-side alignment that highlights every single difference.

* **Input:** Protein sequence files from FastaAAExtractor.
* **Output:** Standardized alignment files that are ready for analysis.
* **Helper Tool:** [SEPI 2.0 Repository](https://github.com/vihaankulkarni29/sepi2.0)

*(Flow: FastaAAExtractor -> **WildTypeAligner** -> SubScan)*

---

### 5. SubScan

**Purpose:** The alignment files from the previous step are detailed but visually cluttered and hard for a human to read. SubScan is an interpreter that transforms this complex data into a simple, clean summary table.

It scans the alignment files (which use symbols like `|`, `.`, `:`) and extracts only the critical information: **Which** protein has a mutation, **what** is the exact change (e.g., Alanine changed to Glycine), and **where** in the protein it occurred. The module provides this final data without making any scientific conclusions—that's the researcher's job.

* **Input:** Alignment files from WildTypeAligner.
* **Output:** A master CSV spreadsheet detailing every mutation found.
* **For more details:** [SubScan Repository](https://github.com/vihaankulkarni29/SubScan)

*(Flow: WildTypeAligner -> **SubScan** -> Co-Occurrence Analyzer)*

---

### 6. Co-Occurrence Analyzer

**Purpose:** This module tackles a more advanced question: do these mutations tend to happen together? It was inspired by research into protein complexes where multiple proteins must work as a team.

It analyzes the master mutation table from SubScan. For each genome, it checks if mutations appear in several of the genes from the user's list. It identifies patterns: In this genome, are all 3 target genes mutated? Or are only genes A and B mutated, but not C? This analysis helps researchers investigate if mutations in different genes are functionally linked.

* **Input:** The master CSV file from SubScan.
* **Output:** A CSV file summarizing the co-occurrence patterns of mutations.

*(Flow: SubScan -> **Co-Occurrence Analyzer** -> HTMLReportGenerator)*

---

### 7. HTML Report Generator

**Purpose:** This is the final step, focused on user-friendly presentation. It gathers all the key results from the entire pipeline and displays them in a single, interactive, and easy-to-understand web page (an HTML file).

It combines the mutation data from **SubScan**, the co-occurrence patterns from the **Co-Occurrence Analyzer**, and quantitative resistance data (MIC values) fetched by a helper module. The report allows a researcher to see the whole story at a glance: what genes are mutated, where the mutations are, whether they appear together, and how resistant the bacterium is.

* **Input:** Data from SubScan, Co-Occurrence Analyzer, and MICHarvester.
* **Output:** A comprehensive, visual HTML report for the user.

*(Flow: SubScan & Co-Occurrence Analyzer -> **HTML Report Generator** -> Pipeline Complete)*