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

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is used to compute Illumina read metrics for each barcode on the trimmed reads. [MultiQC](https://github.com/MultiQC/MultiQC) is used to produce a report containing the FastQC results for all the samples.     

### 3. Genome assembly using Shovill

The Illumina paired-end reads are assembled using the software [Shovill](https://github.com/tseemann/shovill) v1.1.0. Shovill is a pipeline which uses the SPAdes genome assembler at its core. The pipeline uses a default minimum contig length of 200 and a default minimum contig coverage of 10.     

### 4. 	Assembly quality assessment with QUAST

The software [QUAST](https://quast.sourceforge.net/quast.html) v5.2.0 is used to compute genome assembly metrics on the polished assemblies.  

### 5. Assembly quality assessment with CheckM

The software [CheckM](https://github.com/Ecogenomics/CheckM) v1.2.2 (command [lineage_wf](https://github.com/Ecogenomics/CheckM/wiki/Workflows#lineage-specific-workflow)) is used to compute genome assembly completeness and contamination, based on the presence or absence of marker genes. 

### 6. Kraken2/Bracken taxonomy classification

Illumina reads are used as input to the taxonomy classifier [Kraken2](https://github.com/DerrickWood/kraken2) v2.1.3 followed by [Bracken](https://github.com/jenniferlu717/Bracken) v3.0 to estimate abundance of species within a sample. The default Kraken2 database is the PlusPF which contains the Standard database (RefSeq archaea, bacteria, viral, plasmid, human, UniVec_Core) plus RefSeq protozoa and fungi, see [details](https://benlangmead.github.io/aws-indexes/k2). The database was downloaded from https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20240605.tar.gz.   

### 7. LPS typing using Kaptive

The LPS type of the sample is obtained using the software [Kaptive](https://kaptive.readthedocs.io/en/latest/) v3. The genomes assemblies are used as input to this tool. The 9 LPS database is used but can be modified (config parameter "kaptive_db_9lps").  

### 8. 	Variant calling using Snippy

- The reads are mapped to the reference LPS type sequence identified by Kaptive using the sequence alignment program BWA-mem as part of the [snippy](https://github.com/tseemann/snippy) v4.6.0 pipeline.     
- Then Snippy calls variants in the reads as compared to the reference LPS type sequence. It will find both substitutions (snps) and insertions/deletions (indels). 
- Then Snippy uses [SnpEff](https://pcingola.github.io/SnpEff/) to annotate and predict the effects of the variants on genes and proteins (such as amino acid changes).  
- Finally, a bash command is used to extract the variants predicted to have a high impact on the protein (frameshift and stop_gained variants).    

### 9. 	MLST typing

The software [mlst](https://github.com/tseemann/mlst) is used to scan the genome assemblies against the  PubMLST typing scheme "pmultocida_2" by default (RIRDC). The typing scheme can be modified by specifying the parameter --mlst_scheme (e.g. --mlst_scheme "pmultocida"). 

### 10.  LPS subtype report

The pipeline generates a subtype report file (10_Illumina_subtype_report.tsv) summarising the variants found in the subtype database. To be reported, the variant identified by snippy must be present in the subtype database with the following conditions:  
- the variant must be identified at the same position in the same reference sequence and
- both the reference allele and the alternate allele must be matching their corresponding allele from the variant in the database.  

### 11. 	Genome annotation using Bakta

The software [Bakta](https://github.com/oschwengers/bakta) is used to annotate the genome assemblies. The default database is v6.0 from 2025-02-24, https://zenodo.org/records/14916843.    

### 12. 	Antimicrobial Resistance genes

The software [AMRFinderPlus](https://github.com/ncbi/amr) is used to identify AMR genes in the genome assemblies. The default database is the version 2025-03-25.1 that was downloaded using the command amrfinder_update.  

## Step by step user guide

Some files required to use the pipeline are provided to the user (see sections 1a, 1b and 2 below). Additional files must be created/modified by the user (see sections 1c, 3 and 4 below). 

**1) Clone the Github pipeline repository**

Navigate to a folder to your scratch space where you would like to run the pipeline (e.g. $raw_dir below) and clone the pipeline repository to import the required files: 
```
raw_dir=/scratch/project_mnt/SXXX/PIPELINE
cd $raw_dir
git clone https://github.com/vmurigneu/LPS_typing_Illumina.git
```

It will create a repository called "LPS_typing_Illumina". The following three files can be found in the pipeline repository:
- **a) Nextflow configuration file (nextflow.config)**  

When a Nexflow pipeline script is launched, Nextflow looks for a file named **nextflow.config** in the current directory. The configuration file defines default parameters values for the pipeline and cluster settings such as the executor (e.g. "slurm", "local") and queues to be used (https://www.nextflow.io/docs/latest/config.html).  

The pipeline uses separated Singularity containers for all processes. Nextflow will automatically pull the singularity images required to run the pipeline and cache those images in the singularity directory in the pipeline work directory by default or in the singularity.cacheDir specified in the nextflow.config file ([see documentation](https://www.nextflow.io/docs/latest/singularity.html)). Ensure that you have sufficient space in your assigned singularity directory as images can be large.   

An example configuration file can be found [here](https://github.com/vmurigneu/LPS_typing_Illumina/blob/main/nextflow.config). 

- **b) Nextflow main script (main.nf)**

The main.nf script contains the pipeline code and is generally not user-modifiable. 

- **c) Nextflow execution bash script (nextflow.sh)**
This is the bash script used to launch the workflow on the HPC. The template slurm script provided can be used to launch the pipeline on UQ HPC Bunya and is available [here](https://github.com/vmurigneu/LPS_typing_Illumina/blob/main/nextflow.sh). This file should be modified by the user to provide the path to the samplesheet file, Illumina data files etc (see section "Step by step user guide" below). 

**2) Database files for Kraken, Kaptive, CheckM, Bakta and AMRFinderPlus**

Copy the databases folder from the RDM to the cloned pipeline repository on the scratch space (named "dir" below):
```
dir=/scratch/project_mnt/SXXX/PIPELINE/LPS_typing_Illumina
cp -r /QRISdata/Q2313/Valentine/PIPELINES/databases ${dir}
```

**3) Prepare the samplesheet file (csv)**

- The raw Illumina fastq files must be copied in a directory (parameter "--fqdir").
```
fastq=/scratch/project_mnt/SXXX/PIPELINE/LPS_typing_pipeline_Illumina/fastq
mkdir $fastq
cp /path/to/fastq/files/ $fastq
```
  
- The user must specify the path to the raw Illumina files in the samplesheet file. The samplesheet file is a comma-separated values files that defines the names of the samples with their corresponding input fastq files. The header line should match the header line in the examples below. The samplesheet can be saved in a folder named samplesheet e.g. 
```
mkdir /scratch/project_mnt/SXXX/PIPELINE/LPS_typing_Illumina/samplesheet
vim /scratch/project_mnt/SXXX/PIPELINE/LPS_typing_Illumina/samplesheet/samples.csv
```

The samplesheet contains one line for each sample with the following information: the sample identifier (column "sample_id") and the path to the corresponding Illumina paired-end reads file (columns "short_fastq_1" and "short_fastq_2"). File paths are given in relation to the workflow base directory, they are not absolute paths. 
```
sample_id,short_fastq_1,short_fastq_2
PM3034,fastq/1_22VH7WLT3_ATGTCGTATT-TTCTTGCTGG_L002_R1.fastq.gz,fastq/1_22VH7WLT3_ATGTCGTATT-TTCTTGCTGG_L002_R2.fastq.gz
PM3065,fastq/3_22VH7WLT3_GCAATATTCA-GGCGCCAATT_L002_R1.fastq.gz,fastq/3_22VH7WLT3_GCAATATTCA-GGCGCCAATT_L002_R2.fastq.gz
```

**4) Run the pipeline**

The pipeline will be launched on the HPC Bunya using the bash script nextflow.sh. The command to start the pipeline is:  
`nextflow main.nf --samplesheet /path/to/samples.csv --fqdir /path/to/fastq/directory/ --outdir /path/to/outdir/ --slurm_account 'account' `

```
--samplesheet: path to the samplesheet file
--outdir: path to the output directory to be created
--fqdir: path to the directory containing the Illumina fastq files
--slurm_account: name of the Bunya account (default='a_uqds') 
```

Note: To run the assembly and assembly metrics steps only (skip LPS typing and variant calling):  
`nextflow main.nf --samplesheet /path/to/samples.csv --fqdir /path/to/fastq/directory/ --outdir /path/to/outdir/ --slurm_account 'account' --skip_kaptive3 --skip_snippy`

Once the nextflow.sh file is ready, the user can submit the pipeline on Bunya using the command:
```
sbatch nextflow.sh
```

## Optional parameters

Some parameters can be added to the command line in order to include or skip some steps and modify some parameters:

1. Read trimming:
* `--skip_fastp`: skip the read trimming step (default=false). Not recommended. 

2. FastQC reads quality metrics:
* `--skip_fastqc`: skip the FastQC step (default=false)
* `--skip_summary_fastqc`: skip the summary FastQC step using MultiQC (default=false)

3. Genome assembly:
* `--skip_assembly`: skip the assembly step (default=false). Note: it is not recommended to skip assembly as many steps in the downstream processing depends on the assembly results.   
* `--shovill_threads`: number of threads for the assembly (default=4)
* `--shovill_args`: Shovill optional parameters (default="--minlen 200 --mincov 10")
* `--genome_size`: estimated genome size (default="2.3M")

4. Assembly quality assessment with QUAST:
* `--skip_quast`: skip the QUAST step (default=false)
* `--quast_threads`: number of threads for QUAST (default=2)

5. Assembly quality assessment with CheckM:
* `--skip_checkm`: skip the CheckM step (default=false)
* `--checkm_db`: path to the CheckM database folder (default="../../../databases/CheckM-1.2.2)

6. Kraken2/Bracken taxonomy classification:
* `--skip_kraken`: skip the Kraken2/Bracken classification step (default=false)
* `--kraken_db`: path to the Kraken2 database folder (default="../../../databases/k2_pluspf_20240605")

7. LPS typing using Kaptive:
* `--skip_kaptive3`: skip the Kaptive typing step (default=false). note: it will automatically skip the variant calling step.
* `--kaptive_db_9lps`: path to the Kaptive database file (default=""../../../databases/kaptive3_LPS_db_v1/9lps.gbk")

8. Variant calling using Snippy:
* `--skip_snippy`: skip the variant calling Snippy pipeline (default=false)
* `--snippy_threads`: number of threads for the Snippy pipeline (default=6)
* `--reference_LPS`: path to the file summarising the reference LPS sequence files (default="../../../databases/LPS/reference_LPS.txt")

9. MLST typing:
* `--skip_mlst`: skip the MLST typing step (default=false)
* `--mlst_scheme`: MLST typing scheme (default="pmultocida_2")

10. Report:
* `--subtype_db`: path to the subtype database file (default="../../../databases/LPS/LPS_subtype_database_v1.txt")

11. Genome annotation using Bakta:
* `--skip_bakta`: skip the genome annotation step (default=false)
* `--bakta_threads`: number of threads for the Bakta step (default=8)  
* `--bakta_db`: path to the Bakta database files (default="../../../databases/bakta/db")
* `--bakta_args`: Bakta optional parameters (default="--proteins ../../../databases/LPS/NC_002663_LPS.gb")  

12. AMR genes identification using AMRFinderPlus:
* `--skip_amrfinder`: skip the AMR genes identification step (default=false)
* `--amrfinder_db`: path to the AMRFinderPlus database files (default="../../../databases/amrfinderplus/2025-03-25.1")
* `--amrfinder_args`: AMRFinderPlus optional parameters (default="")
  
## Structure of the output folders

The pipeline will create several folders corresponding to the different steps of the pipeline. 
The main output folder (`--outdir`) will contain a folder per sample (the folder is named as in the column sample_id in the samplesheet file).

Each sample folder will contain the following folders:
* **1_trimming:** Paired-end trimmed fastq files (sample_id_R1_trimmed.fastq.gz and sample_id_R2_trimmed.fastq.gz).
* **2_fastqc:** FastQC quality control results for the paired-end reads:
    * FastQC report in html format (sample_id_R1_trimmed_fastqc.html and sample_id_R2_trimmed_fastqc.html)
    * FastQC zipped results folder (sample_id_R1_trimmed_fastqc.zip and sample_id_R2_trimmed_fastqc.zip)
* **3_assembly:** Shovill assembly output files, see [details](https://github.com/tseemann/shovill?tab=readme-ov-file#output-files).
    * Final genome assembly in fasta format (sample_id_contigs.fa)
    * SPADES assembly graph (sample_id_contigs.gfa)
* **4_quast:** QUAST output report file (sample_id_report.tsv).
* **5_checkm:** CheckM output file (sample_id_checkm_lineage_wf_results.tsv).  
* **6_kraken:**  Kraken2/Bracken taxonomy classification results, see output files format details [here](https://github.com/DerrickWood/kraken2/wiki/Manual#output-formats) and [here](https://ccb.jhu.edu/software/bracken/index.shtml?t=manual#format)
  * Kraken2 classification report (sample_id_kraken2_report.txt)  
  * Kraken2 classification assignments for a read (sample_id_kraken2.tsv.gz)  
  * Bracken species abundance results (sample_id_bracken_species.tsv)   
  * Bracken results in Kraken style report format (sample_id_bracken_report.txt)  
* **7_kaptive_v3:** Kaptive output files, see [details](https://kaptive.readthedocs.io/en/latest/Outputs.html)
    * LPS type results (sample_id_kaptive_results.tsv)
    * LPS sequence in fasta format (sample_id_kaptive_results.fna)
* **8_snippy:** Mapping files and variant calling results from Snippy, see [details](https://github.com/tseemann/snippy?tab=readme-ov-file#output-files):
    * BWA mapping file in bam format (sample_id_snps.bam_mapped.bam and .bai index). 
    * Unfiltered variants from Freebayes in VCF format (sample_id_clair_snps.raw.vcf) 
    * Filtered variants from Freebayes in VCF format (sample_id_snps.filt.vcf)
    * Summary of variants in tabular format (sample_id_snps.tab)
    * Summary of high impact variants (frameshift_variant and stop_gained) in tabular format (sample_id_snps.high_impact.tab)
* **9_mlst:** MLST typing output file (sample_id_mlst_pmultocida_rirdc.csv) 
* **10_report:** Summary of results for all samples
    * MultiQC report in html format (2_Illumina_multiqc_report.html) and general statistics in tabular format (2_Illumina_multiqc_general_stats.txt)
    * Shovill assembly statistics: assembly coverage, number of contigs, assembly size (3_Illumina_shovill_stats.tsv) 
    * QUAST combined report file (4_Illumina_quast_report.tsv)  
    * Checkm results (5_Illumina_checkm_lineage_wf_results.tsv)  
    * Kraken/Bracken taxonomy results:  
        - Abundance of P. multocida reads: 6_Illumina_bracken_pasteurella_multocida_species_abundance.tsv
        - Information about the most abundant species identified: 6_Illumina_bracken_most_abundant_species.tsv    
    * Kaptive results (7_Illumina_kaptive_results.tsv)  
    * Snippy variants results:  
        - all variants: 8_Illumina_snippy_snps.tsv 
        - only variants predicted to have a high impact on the protein: 8_Illumina_snippy_snps.high_impact.tsv 
    * MLST results (9_Illumina_mlst.csv)  
    * Subtype results summarising the variants found in the subtype database (10_Illumina_subtype_report.tsv). To be reported, the variant identified by snippy must be present in the subtype database with the following conditions:
       - the variant must be identified at the same position in the reference sequence and
       - both the reference allele and the alternate allele must be matching their corresponding allele from the variant in the database.
    * AMRFinderPlus results (12_Illumina_amrfinder.tsv)
* **11_bakta:** Bakta genome annotation output files. The output files are described [here](https://github.com/oschwengers/bakta?tab=readme-ov-file#output) and include:  
    * Annotations & sequences in (multi) GenBank format (sample_id_bakta.gbff)  
    * Inference metrics (score, evalue, coverage, identity) for annotated accessions as TSV (sample_id_bakta.inference.tsv)  
    * Annotation summary in text format (sample_id_bakta.txt)    
* **12_amrfinder:** AMRFinderPlus output file (sample_id_amrfinder.tsv). The output in tab-delimited format is described [here](https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus#output-format). 
## Running the workflow in assembly mode for other organisms

The default parameters are suited for Pasteurella multocida. The LPS typing and variant calling are specific to Pasteurella multocida. Here are the paraneters to use the workflow to assemble another species:  
* `--genome_size`: estimated genome size (default="2.3M")
* `--mlst_scheme`: MLST typing scheme (default="pmultocida_2")
* `--skip_kaptive3`: skip the Kaptive typing step (default=false). note: it will automatically skip the variant calling step. 
