# GenomeAMRAnalyzer: Project Investigator Testing Guide

## 🎯 Quick Start for Project Investigators

**Time Required**: 5-10 minutes for basic test  
**Technical Knowledge**: None required  
**Prerequisites**: Windows computer with internet connection

---

## 📋 Option 1: Instant Test (Recommended for First Time)

### Step 1: Download and Setup
1. **Download the pipeline** from GitHub or receive it from your team
2. **Open Command Prompt** (search "cmd" in Windows start menu)
3. **Navigate to the pipeline folder**:
   ```cmd
   cd C:\path\to\GenomeAMRAnalyzer
   ```

### Step 2: Run Basic Test
**Copy and paste this exact command** (replace with your email):
```cmd
python genomeamr_auto.py --accessions config/genes_default.txt --genes config/genes_default.txt --email your.email@institution.edu
```

**What this does**: Tests the pipeline with built-in example data

### Step 3: Expected Results
- ✅ Pipeline runs for 2-5 minutes
- ✅ Shows progress messages
- ✅ Generates HTML report in `reports/pipeline_report.html`
- ✅ No error messages

**If successful**: Open `reports/pipeline_report.html` in your web browser to see the analysis report.

---

## 📋 Option 2: Test with Real NCBI Data

### Step 1: Choose Your Research Interest
Pick one of these ready-to-use examples:

#### **Example A: E. coli Antibiotic Resistance**
```cmd
python genomeamr_auto.py --url "https://www.ncbi.nlm.nih.gov/nuccore/?term=escherichia+coli+antibiotic+resistance" --genes config/genes_default.txt --email your.email@institution.edu
```

#### **Example B: Erythromycin Resistance**
```cmd
python genomeamr_auto.py --url "https://www.ncbi.nlm.nih.gov/nuccore/?term=erythromycin+resistance" --genes config/genes_erythromycin.txt --email your.email@institution.edu
```

#### **Example C: Custom Gene List**
1. Create a text file called `my_genes.txt` with genes of interest (one per line):
   ```
   acrB
   acrA
   tolC
   mexA
   mexB
   ```
2. Run:
   ```cmd
   python genomeamr_auto.py --url "https://www.ncbi.nlm.nih.gov/nuccore/?term=your+search+terms" --genes my_genes.txt --email your.email@institution.edu
   ```

### Step 2: Monitor Progress
The pipeline will show:
- **[STEP 1]** Downloading genomes from NCBI
- **[STEP 2]** Running CARD analysis for resistance genes
- **[STEP 3]** Extracting protein sequences
- **[STEP 4]** Running alignment analysis
- **[STEP 5]** Analyzing mutations
- **[STEP 6]** Running co-occurrence analysis
- **[STEP 7]** Generating HTML report

**Total time**: 5-30 minutes depending on genome count

---

## 🖱️ Option 3: One-Click Testing (Easiest)

### Windows Users:
1. **Double-click** `quick_start.bat`
2. **Follow the prompts**:
   - Enter your email when asked
   - Provide NCBI URL or choose built-in example
   - Optionally upload a gene list file
3. **Wait for completion**
4. **View results** in automatically opened report

### Mac/Linux Users:
1. **Double-click** `quick_start.sh` (or run in terminal)
2. **Follow the same prompts as above**

---

## 📊 What to Look For in Results

### Successful Test Indicators:
✅ **No error messages during execution**  
✅ **HTML report generated** (`reports/pipeline_report.html`)  
✅ **Report opens in web browser**  
✅ **Report shows**:
   - Genome count processed
   - Resistance genes found
   - Mutation analysis results
   - Co-occurrence patterns
   - Professional formatting

### Report Contents to Verify:
1. **Executive Summary**: Key findings overview
2. **Genome Analysis**: Number of genomes processed successfully
3. **Resistance Gene Detection**: Genes found with confidence scores
4. **Mutation Analysis**: Specific mutations identified
5. **Co-occurrence Patterns**: Gene interaction networks
6. **Clinical Significance**: Interpretation of results

---

## 🔍 Testing Different Scenarios

### Test 1: Small Dataset (Quick Test)
```cmd
python genomeamr_auto.py --accessions config/genes_default.txt --genes config/genes_default.txt --email test@example.com --max-genomes 5
```
**Expected time**: 2-3 minutes

### Test 2: Medium Dataset (Realistic Test)
```cmd
python genomeamr_auto.py --url "https://www.ncbi.nlm.nih.gov/nuccore/?term=escherichia+coli+resistance" --genes config/genes_default.txt --email test@example.com --max-genomes 20
```
**Expected time**: 10-15 minutes

### Test 3: Custom Research Question
```cmd
python genomeamr_auto.py --url "YOUR_NCBI_SEARCH_URL" --genes YOUR_GENE_LIST.txt --email your.email@institution.edu
```

---

## 🚨 Troubleshooting Common Issues

### Problem: "Command not found" or "Python not recognized"
**Solution**: Install Python first
1. Visit [python.org](https://python.org)
2. Download Python 3.8 or newer
3. During installation, check "Add Python to PATH"
4. Restart command prompt and try again

### Problem: "Gene list file not found"
**Solution**: Check file path
```cmd
# Use quotes for paths with spaces
python genomeamr_auto.py --genes "C:\Users\YourName\Documents\gene_list.txt" --email test@example.com
```

### Problem: Pipeline stops with network error
**Solution**: Check internet connection and try again
- The pipeline downloads data from NCBI
- Temporary network issues can be resolved by rerunning

### Problem: "Permission denied" errors
**Solution**: Run as administrator
1. Right-click Command Prompt
2. Select "Run as administrator"
3. Navigate to pipeline folder and run again

### Problem: Long processing time
**Expected behavior**: 
- 1-5 genomes: 2-5 minutes
- 10-20 genomes: 10-20 minutes  
- 50+ genomes: 30-60 minutes

---

## 📈 Validating Scientific Results

### Check 1: Genome Processing Success Rate
- **Good**: 80-100% success rate
- **Acceptable**: 60-80% success rate
- **Investigate**: <60% success rate

### Check 2: Resistance Gene Detection
- **Realistic**: 1-10 resistance genes per genome for typical bacteria
- **High**: 10+ genes (multi-drug resistant strains)
- **Suspicious**: 0 genes (check gene list relevance)

### Check 3: Mutation Analysis
- **Normal**: 0-5 mutations per gene
- **High resistance**: 5+ mutations per gene
- **Quality check**: Mutations have confidence scores >0.7

### Check 4: Report Quality
- **Professional appearance**: Clean formatting, clear charts
- **Complete data**: All sections populated with results
- **Export capability**: Can download data as CSV/PDF

---

## 🎯 Test Success Criteria

### Minimal Success (Pipeline Works):
- [x] Command executes without errors
- [x] HTML report generated
- [x] Report contains data (not empty)

### Scientific Success (Results are Meaningful):
- [x] Appropriate number of genomes processed
- [x] Resistance genes detected match expectations
- [x] Mutations identified with confidence scores
- [x] Co-occurrence patterns make biological sense

### User Experience Success (Easy to Use):
- [x] Clear progress messages during execution
- [x] Intuitive report layout
- [x] Actionable error messages if issues occur
- [x] Results interpretable without technical expertise

---

## 📞 Getting Help

### If Tests Fail:
1. **Copy the exact error message**
2. **Note which step failed** (Step 1-7)
3. **Include the command you ran**
4. **Contact your technical team with this information**

### For Scientific Questions:
- **Gene selection**: Ask your team about relevant genes for your research
- **Result interpretation**: Consult with bioinformatics experts
- **Clinical significance**: Discuss with clinicians familiar with AMR

### For Technical Issues:
- **Installation problems**: Contact IT support
- **Performance issues**: Check computer specifications
- **Network issues**: Verify internet connectivity

---

## 🌟 Next Steps After Successful Testing

### For Research Use:
1. **Define your research question** clearly
2. **Prepare your gene list** based on research focus
3. **Identify relevant NCBI searches** for your organisms of interest
4. **Plan analysis schedule** based on processing times observed

### For Production Use:
1. **Validate results** against known datasets
2. **Establish quality control** procedures
3. **Train team members** on pipeline usage
4. **Document your specific workflows** for reproducibility

---

**Remember**: This pipeline is designed to be foolproof. If something doesn't work as described, it's likely a setup issue, not a fundamental problem with your approach. Don't hesitate to ask for technical support!