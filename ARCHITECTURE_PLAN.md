# URL-Based AMR Genome Analysis Pipeline: Complete Architecture Plan

## 🏗️ **System Architecture Overview**

### **Core Philosophy**
Transform from static accession lists to dynamic NCBI query-driven genome discovery, maintaining robust error handling and scientific reproducibility.

---

## 📋 **Phase 1: Foundation & URL Processing (Days 1-2)**

### **1.1 URL Parser & Query Engine**
```
Input:  NCBI Search URL (e.g., E. coli + macrolide resistance)
Output: Parsed search terms and parameters
```

**Components:**
- **URLParser**: Extract search terms from complex NCBI URLs
- **QueryValidator**: Ensure search terms are scientifically valid
- **ParameterExtractor**: Handle filters (organism, genome type, etc.)

**Deliverables:**
- `url_parser.py` - Robust URL parsing with multiple format support
- `query_validator.py` - Scientific term validation
- Unit tests for all URL formats

---

### **1.2 NCBI E-utilities Integration**
```
Input:  Parsed search terms
Output: List of genome accessions with metadata
```

**Components:**
- **ESearchClient**: Query NCBI databases
- **ESummaryClient**: Extract metadata (organism, size, title)
- **ResultFilter**: Apply quality filters (genome size, completeness)

**Deliverables:**
- `ncbi_client.py` - Professional E-utilities wrapper
- Metadata extraction and validation
- Rate limiting and API key support

---

## 📋 **Phase 2: Genome Acquisition (Days 2-3)**

### **2.1 Parallel Download Engine**
```
Input:  Filtered accession list (max 10 for testing)
Output: Downloaded FASTA files with validation
```

**Components:**
- **AsyncDownloader**: Parallel downloads with semaphore control
- **FileValidator**: FASTA format and content validation
- **RetryHandler**: Exponential backoff for failed downloads

**Deliverables:**
- `download_engine.py` - Production-grade downloader
- Comprehensive error handling and logging
- Progress tracking and status reporting

---

### **2.2 Quality Control System**
```
Input:  Raw downloaded genomes
Output: Validated, quality-controlled FASTA files
```

**Components:**
- **GenomeValidator**: Size, format, completeness checks
- **QualityReporter**: Detailed QC metrics
- **FileOrganizer**: Structured output directory management

---

## 📋 **Phase 3: AMR Gene Discovery (Days 3-4)**

### **3.1 CARD RGI Integration**
```
Input:  Quality-controlled genomes
Output: AMR gene coordinates and classifications
```

**Components:**
- **RGIRunner**: Robust RGI execution with error handling
- **CoordinateParser**: Extract gene positions from RGI JSON
- **ResultValidator**: Ensure RGI outputs are complete

**Deliverables:**
- `rgi_integration.py` - Professional RGI wrapper
- JSON parsing and coordinate extraction
- Fallback strategies for RGI failures

---

### **3.2 Gene Coordinate Extraction**
```
Input:  RGI results + target gene list
Output: Precise coordinates for user-specified genes
```

**Components:**
- **CoordinateExtractor**: Parse RGI JSON for target genes
- **GeneValidator**: Validate against user gene list
- **MissingGeneHandler**: Report and handle missing genes

---

## 📋 **Phase 4: Gene Sequence Extraction (Days 4-5)**

### **4.1 FastaAAExtractor Enhancement**
```
Input:  Genome FASTA + gene coordinates
Output: Extracted gene/protein sequences
```

**Components:**
- **SequenceExtractor**: Precise coordinate-based extraction
- **TranslationEngine**: DNA to protein conversion
- **SequenceValidator**: Quality checks on extracted sequences

**Deliverables:**
- Enhanced `fasta_aa_extractor.py`
- Support for both DNA and protein extraction
- Comprehensive sequence validation

---

### **4.2 Multi-Genome Processing**
```
Input:  Multiple genomes + coordinates
Output: Organized gene sequences by genome and gene
```

**Components:**
- **BatchProcessor**: Handle multiple genomes efficiently
- **OutputOrganizer**: Structure results by gene/genome
- **SequenceDeduplicator**: Handle identical sequences

---

## 📋 **Phase 5: Integration & Testing (Days 5-6)**

### **5.1 Pipeline Orchestration**
```
Input:  NCBI URL + target genes
Output: Complete AMR analysis results
```

**Components:**
- **PipelineRunner**: End-to-end workflow execution
- **StatusTracker**: Real-time progress monitoring
- **ErrorRecovery**: Graceful failure handling

**Deliverables:**
- `amr_pipeline.py` - Complete workflow orchestrator
- Configuration-driven execution
- Comprehensive logging and reporting

---

### **5.2 Testing Framework**
```
Input:  Test URLs and gene lists
Output: Validated pipeline functionality
```

**Components:**
- **UnitTests**: Individual component testing
- **IntegrationTests**: End-to-end workflow testing
- **PerformanceTests**: Speed and resource usage validation

---

## 📋 **Phase 6: Configuration & Documentation (Day 6)**

### **6.1 Configuration Updates**
Update `snakemake_config.yaml` to support URL-based workflow:

```yaml
# URL-based genome discovery
ncbi_search:
  max_genomes: 10
  quality_filters:
    min_size: 1000000
    max_size: 20000000
  search_url: ""  # User-provided NCBI URL

# Enhanced RGI settings
rgi:
  timeout_seconds: 300
  fallback_enabled: true
  coordinate_validation: true
```

### **6.2 Documentation & Examples**
- Complete README with URL-based workflow
- Example URLs for different organisms
- Troubleshooting guide
- API documentation

---

## 🎯 **Implementation Strategy**

### **Day 1: Foundation**
1. ✅ URL parsing and validation
2. ✅ NCBI E-utilities client
3. ✅ Basic testing framework

### **Day 2: Download Engine**
1. ✅ Parallel download implementation
2. ✅ Quality control system
3. ✅ Error handling and retries

### **Day 3: AMR Integration**
1. ✅ RGI wrapper and execution
2. ✅ Coordinate parsing and validation
3. ✅ Gene extraction logic

### **Day 4: Sequence Processing**
1. ✅ Enhanced FastaAAExtractor
2. ✅ Multi-genome batch processing
3. ✅ Output organization

### **Day 5: Integration**
1. ✅ Complete pipeline orchestration
2. ✅ End-to-end testing
3. ✅ Performance optimization

### **Day 6: Polish**
1. ✅ Configuration updates
2. ✅ Documentation completion
3. ✅ Final validation

---

## 🔧 **Technical Decisions**

### **Dependencies**
- **asyncio + aiohttp**: Parallel downloads
- **xml.etree.ElementTree**: NCBI XML parsing
- **pathlib**: Modern file handling
- **dataclasses**: Clean data structures
- **logging**: Comprehensive logging

### **Architecture Patterns**
- **Async/Await**: Non-blocking I/O operations
- **Dependency Injection**: Testable components
- **Error Boundaries**: Isolated failure handling
- **Configuration-Driven**: Flexible execution

### **Quality Assurance**
- **Type Hints**: Full type coverage
- **Docstrings**: Comprehensive documentation
- **Unit Tests**: >90% code coverage
- **Integration Tests**: Real-world scenarios

---

## 🎯 **Success Metrics**

1. **Functionality**: Successfully process NCBI URLs with 3000+ results
2. **Performance**: Download and process 10 genomes in <5 minutes
3. **Reliability**: Handle network failures gracefully
4. **Usability**: Single command execution from URL to results
5. **Maintainability**: Clean, documented, testable code

---

## 🚀 **Next Steps**

Ready to begin implementation? The plan provides:
- Clear component boundaries
- Incremental development approach
- Built-in testing and validation
- Production-ready architecture

**Shall we start with Phase 1: URL Parser & NCBI Integration?**