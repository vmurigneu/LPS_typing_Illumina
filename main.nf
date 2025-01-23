#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================
        Pasteurella multocida LPS analysis pipeline
========================================================================================
 #### Documentation
 #https://github.com/vmurigneu/LPS_typing_Illumina
 #### Authors
 Valentine Murigneux <v.murigneux@uq.edu.au>
========================================================================================
*/

def helpMessage() {
	log.info"""
	=========================================
	Pasteurella multocida LPS analysis pipeline v${workflow.manifest.version}
	=========================================
	Usage:
	nextflow main.nf --fqdir /path/to/fastq/directory/ --outdir /path/to/outdir/

	Required arguments:
		--fqdir					Path to the directory containing the Illumina fastq files
		--outdir				Path to the output directory to be created
		--samplesheet				Path to the samplesheet file
    
	Optional parameters:
		--threads				Number of threads (default=4)


    """.stripIndent()
}

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

process fastp {
	cpus "${params.threads}"
	tag "${sample}"
	label "cpu"
	publishDir "$params.outdir/$sample/1_trimming",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/1_trimming",  mode: 'copy', pattern: '*trimmed.fastq.gz', saveAs: { filename -> "${sample}_$filename" }
	input:
		tuple val(sample), path(reads1), path(reads2)
	output:
		tuple val(sample), path(reads1), path(reads2), path("R1_trimmed.fastq.gz"), path("R2_trimmed.fastq.gz"),  emit: trimmed_fastq
		path("fastp.log")
		path("*fastq.gz")
	script:
	"""
	fastp -i ${reads1} -I ${reads2} -o R1_trimmed.fastq.gz -O R2_trimmed.fastq.gz
	cp .command.log fastp.log
	"""
}

process fastqc {
        cpus "${params.threads}"
        tag "${sample}"
        label "cpu"
        publishDir "$params.outdir/$sample/2_fastqc",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
        publishDir "$params.outdir/$sample/2_fastqc",  mode: 'copy', pattern: '*fastqc.zip', saveAs: { filename -> "${sample}_$filename" }
        publishDir "$params.outdir/$sample/2_fastqc",  mode: 'copy', pattern: '*fastqc.html', saveAs: { filename -> "${sample}_$filename" }
	input:
                tuple val(sample), path(reads1), path(reads2), path(reads1_trimmed), path(reads2_trimmed)
        output:
                tuple val(sample), path(reads1), path(reads2), path(reads1_trimmed), path(reads2_trimmed), emit: reads_qc
                path("fastqc.log")
		path("*fastqc.zip")
                path("*fastqc.html")
        when:
        !params.skip_fastqc
        script:
        """
        fastqc -o \$PWD ${reads1_trimmed} 
	fastqc -o \$PWD ${reads2_trimmed} 
        cp .command.log fastqc.log
        """
}

process shovill {
        cpus "${params.shovill_threads}"
        tag "${sample}"
        label "cpu"
        label "high_memory"
	publishDir "$params.outdir/$sample/3_assembly",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
        publishDir "$params.outdir/$sample/3_assembly",  mode: 'copy', pattern: '*fa', saveAs: { filename -> "${sample}_$filename" }
        input:
                tuple val(sample), path(reads1), path(reads2)
        output:
                tuple val(sample), path("contigs.fa"), emit: assembly_out
                path("shovill.log")
		path("*fa")
        script:
        """
        shovill --outdir \$PWD --R1 ${reads1} --R2 ${reads2} --gsize ${params.genome_size} --force --cpus ${params.shovill_threads} 
        cp .command.log shovill.log
        """
}

process quast {
        cpus "${params.threads}"
        tag "${sample}"
        label "cpu"
	publishDir "$params.outdir/$sample/4_quast",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/4_quast",  mode: 'copy', pattern: '*tsv', saveAs: { filename -> "${sample}_$filename" }
	input:
                tuple val(sample), path(assembly)
        output:
		tuple val(sample), path("report.tsv"), emit: quast_results
                path("quast.log")
        when:
        !params.skip_quast
        script:
        """
	quast.py ${assembly} --threads ${params.threads} -o \$PWD
        cp .command.log quast.log
        """
}

process checkm {
        cpus "${params.threads}"
        tag "${sample}"
        label "cpu"
        label "high_memory"
        publishDir "$params.outdir/$sample/5_checkm",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
        publishDir "$params.outdir/$sample/5_checkm",  mode: 'copy', pattern: '*tsv', saveAs: { filename -> "${sample}_$filename" }
        input:
                tuple val(sample), path(assembly)
        output:
                tuple val(sample), path("checkm_lineage_wf_results.tsv"),  emit: checkm_results
                path("checkm.log")
        when:
        !params.skip_checkm
        script:
        """
        export CHECKM_DATA_PATH=${params.checkm_db}
        checkm data setRoot ${params.checkm_db}
        checkm lineage_wf --reduced_tree `dirname ${assembly}` \$PWD --threads ${params.threads} --pplacer_threads ${params.threads} --tab_table -f checkm_lineage_wf_results.tsv -x fa
        cp .command.log checkm.log
        """
}

process kraken {
        cpus "${params.threads}"
        tag "${sample}"
        label "cpu"
        label "high_memory"
	publishDir "$params.outdir/$sample/6_kraken",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/6_kraken",  mode: 'copy', pattern: '*txt', saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/6_kraken",  mode: 'copy', pattern: '*tsv.gz', saveAs: { filename -> "${sample}_$filename" }
        input:
                tuple val(sample), path(reads1), path(reads2), path(reads1_trimmed), path(reads2_trimmed)
        output:
                tuple val(sample), path(reads1_trimmed), path(reads2_trimmed), path("kraken2_report.txt"), path("kraken2.tsv.gz"),  emit: kraken_results
                path("kraken.log")
        when:
        !params.skip_kraken
        script:
        """
	kraken2 --gzip-compressed --db ${params.kraken_db} --report kraken2_report.txt --paired ${reads1_trimmed} ${reads2_trimmed} > kraken2.tsv
        gzip kraken2.tsv
	cp .command.log kraken.log
        """
}

process bracken {
        cpus "${params.bracken_threads}"
        tag "${sample}"
        label "cpu"
        label "high_memory"
        publishDir "$params.outdir/$sample/6_kraken",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
        publishDir "$params.outdir/$sample/6_kraken",  mode: 'copy', pattern: '*txt', saveAs: { filename -> "${sample}_$filename" }
        publishDir "$params.outdir/$sample/6_kraken",  mode: 'copy', pattern: '*tsv', saveAs: { filename -> "${sample}_$filename" }
        input:
                tuple val(sample), path(reads1_trimmed), path(reads2_trimmed), path(kraken2_report), path(kraken_tsv)
        output:
                tuple val(sample), path(reads1_trimmed), path(reads2_trimmed), path("bracken_report.txt"), path("bracken_species.tsv"),  emit: bracken_results
                path("bracken.log")
        when:
        !params.skip_kraken
        script:
        """
        bracken -d ${params.kraken_db} -i ${kraken2_report} -o bracken_species.tsv  -w bracken_report.txt -r 100 -l S -t ${params.bracken_threads}
        cp .command.log bracken.log
        """
}

process kaptive3 {
        cpus "${params.threads}"
        tag "${sample}"
        label "cpu"
        publishDir "$params.outdir/$sample/7_kaptive_v3",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
        publishDir "$params.outdir/$sample/7_kaptive_v3",  mode: 'copy', pattern: '*tsv', saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/7_kaptive_v3",  mode: 'copy', pattern: '*fna', saveAs: { filename -> "${sample}_$filename" }
        input:
                tuple val(sample), path(assembly)
        output:
                tuple val(sample), path("kaptive_results.tsv"),  emit: kaptive_results
		path("*fna")
                path("kaptive_v3.log")
        when:
        !params.skip_kaptive3
        script:
        """
	kaptive assembly ${params.kaptive_db_9lps} ${assembly} -f \$PWD -o kaptive_results.tsv
        cp .command.log kaptive_v3.log
        """
}

process snippy {
        cpus "${params.snippy_threads}"
        tag "${sample}"
        label "cpu"
        publishDir "$params.outdir/$sample/8_snippy",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
        publishDir "$params.outdir/$sample/8_snippy",  mode: 'copy', pattern: '*tab', saveAs: { filename -> "${sample}_$filename" }
        publishDir "$params.outdir/$sample/8_snippy",  mode: 'copy', pattern: 'snps*', saveAs: { filename -> "${sample}_$filename" }
	input:
                tuple val(sample), path(reads1), path(reads2), path(reads1_trimmed), path(reads2_trimmed), path(kaptive_report)
	output:
                tuple val(sample), path("snps.tab"), path("snps.high_impact.tab"), path("snps.raw.vcf"), path("snps.filt.vcf"),  path("snps.bam"), path("snps.bam.bai"),  emit: snippy_results
		path("snippy.log")
        when:
        !params.skip_snippy
        shell:
        '''
	locus=`tail -1 !{kaptive_report} | cut -f3`
	ref_gb=`grep ${locus:0:2} !{params.reference_LPS} | cut -f2`
	snippy --cpus !{params.snippy_threads} --force --outdir \$PWD --ref $ref_gb --R1 !{reads1_trimmed} --R2 !{reads2_trimmed}
        egrep "^CHROM|frameshift_variant|stop_gained" snps.tab > snps.high_impact.tab
	cp .command.log snippy.log
        '''
}

process mlst {
        cpus "${params.threads}"
        tag "${sample}"
        label "cpu"
        publishDir "$params.outdir/$sample/9_mlst",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
        publishDir "$params.outdir/$sample/9_mlst",  mode: 'copy', pattern: '*csv', saveAs: { filename -> "${sample}_$filename" }
        input:
                tuple val(sample), path(assembly)
        output:
                tuple val(sample), path("mlst_pmultocida_rirdc.csv"),  emit: mlst_results
                path("mlst.log")
        when:
        !params.skip_mlst
        script:
        """
	mlst --scheme pmultocida_2 ${assembly} --quiet --csv --threads ${params.threads} > mlst_pmultocida_rirdc.csv
        cp .command.log mlst.log
        """
}

workflow {
	Channel.fromPath( "${params.samplesheet}", checkIfExists:true )
	.splitCsv(header:true, sep:',')
	.map { row -> tuple(row.sample_id, file(row.short_fastq_1, checkIfExists: true), file(row.short_fastq_2, checkIfExists: true)) }
	.set { ch_samplesheet_illumina }
	ch_samplesheet_illumina.view()
	fastp(ch_samplesheet_illumina)
	if (!params.skip_fastqc) {
		fastqc(fastp.out.trimmed_fastq)
	}
	shovill(ch_samplesheet_illumina)
	if (!params.skip_quast) {
		quast(shovill.out.assembly_out)
	}
	if (!params.skip_kraken) {
		kraken(fastp.out.trimmed_fastq)
		bracken(kraken.out.kraken_results)
	}
	if (!params.skip_checkm) {
		checkm(shovill.out.assembly_out)
	}
	if (!params.skip_kaptive3) {
		kaptive3(shovill.out.assembly_out)
	}
	if (!params.skip_snippy) {
		snippy(fastp.out.trimmed_fastq.join(kaptive3.out.kaptive_results))
	}
	if (!params.skip_mlst) {
		mlst(shovill.out.assembly_out)
	}
}
