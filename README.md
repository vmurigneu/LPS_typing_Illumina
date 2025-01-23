# LPS_typing_Illumina
Bioinformatics pipeline for Pasteurella multocida LPS typing using Illumina sequencing data

  - [Overall pipeline](#Overall-pipeline)
  - [User guide](#Step-by-step-user-guide)
  - [Example data](#Example-data)
  - [Optional parameters](#Optional-parameters)
  - [Output files](#structure-of-the-output-folders)
  - [Acknowledgements/citations/credits](#acknowledgements--citations--credits)
    
## Overall pipeline 

### 1. Read trimming

The raw Illumina reads are trimmed using [fastp](https://github.com/OpenGene/fastp) v0.24.0.   

### 2. Illumina reads quality metrics 

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is used to compute Illumina read metrics for each barcode on the trimmed reads.   

### 3. Genome assembly using Shovill

The Illumina paired-end reads are assembled using the software [Shovill](https://github.com/tseemann/shovill) v1.1.0. Shovill is a pipeline which uses the SPAdes genome assembler at its core.    

### 4. 	Assembly quality assessment with QUAST

The software [QUAST](https://quast.sourceforge.net/quast.html) v5.2.0 is used to compute genome assembly metrics on the polished assemblies.  

### 5. Assembly quality assessment with CheckM

The software [CheckM](https://github.com/Ecogenomics/CheckM) v1.2.2 (command [lineage_wf](https://github.com/Ecogenomics/CheckM/wiki/Workflows#lineage-specific-workflow)) is used to compute genome assembly completeness and contamination, based on the presence or absence of marker genes. 

### 6. Kraken2/Bracken taxonomy classification

Illumina reads are used as input to the taxonomy classifier [Kraken2](https://github.com/DerrickWood/kraken2) v2.1.3 followed by [Bracken](https://github.com/jenniferlu717/Bracken) v3.0 to estimate abundance of species within a sample. The default Kraken2 database is the PlusPF which contains Standard plus RefSeq protozoa & fungi (see [details](https://benlangmead.github.io/aws-indexes/k2)). The database was downloaded from https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20240605.tar.gz. 

### 7. LPS typing using Kaptive

The LPS type of the sample is obtained using the software [Kaptive](https://kaptive.readthedocs.io/en/latest/) v3. The genomes assemblies are used as input to this tool. The 9 LPS database is used but can be modified (config parameter "kaptive_db_9lps").  

### 8. 	Variant calling using Snippy

- The reads are mapped to the reference LPS type sequence identified by Kaptive using the sequence alignment program BWA-mem as part of the [snippy](https://github.com/tseemann/snippy) v4.6.0 pipeline.     
- Then Snippy calls variants in the reads as compared to the reference LPS type sequence. It will find both substitutions (snps) and insertions/deletions (indels). 
- Then Snippy uses [SnpEff](https://pcingola.github.io/SnpEff/) to annotate and predict the effects of the variants on genes and proteins (such as amino acid changes).  
- Finally, a bash command is used to extract the variants predicted to have a high impact on the protein (frameshift and stop_gained variants).    

### 9. 	MLST typing

The software [mlst](https://github.com/tseemann/mlst) is used to scan the genome assemblies against the  PubMLST typing scheme "pmultocida_2".    

