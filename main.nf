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
        publishDir "$params.outdir/$sample/2_fastqc",  mode: 'copy', pattern: '*fastqc.zip'
        publishDir "$params.outdir/$sample/2_fastqc",  mode: 'copy', pattern: '*fastqc.html', saveAs: { filename -> "${sample}_$filename" }
	input:
                tuple val(sample), path(reads1), path(reads2), path(reads1_trimmed), path(reads2_trimmed)
        output:
                tuple val(sample), path(reads1), path(reads2), path(reads1_trimmed), path(reads2_trimmed), emit: reads_qc
		path("fastqc.log")
		path("*fastqc.zip"), emit: fastqc_zip
                path("*fastqc.html")
        when:
        !params.skip_fastqc
        script:
        """
        fastqc -o \$PWD ${reads1_trimmed} 
	fastqc -o \$PWD ${reads2_trimmed} 
        mv R1_trimmed_fastqc.zip ${sample}_R1_trimmed_fastqc.zip
	mv R2_trimmed_fastqc.zip ${sample}_R2_trimmed_fastqc.zip
	cp .command.log fastqc.log
        """
}

process summary_fastqc {
	publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*html'
	publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*txt'
	input:
		path(fastqc_files)
	output:
		path("2_multiqc_report.html"), emit: fastqc_summary
		path("2_multiqc_general_stats.txt"), emit: fastqc_stats
	when:
	!params.skip_summary_fastqc
	script:
	"""
	multiqc --fn_as_s_name .
	cp .command.log summary_fastqc.log
	mv multiqc_report.html 2_multiqc_report.html
	mv multiqc_data/multiqc_general_stats.txt 2_multiqc_general_stats.txt
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
	publishDir "$params.outdir/$sample/4_quast",  mode: 'copy', pattern: '*tsv'
	input:
                tuple val(sample), path(assembly)
        output:
		path("*report.tsv"), emit: quast_results
                path("quast.log")
        when:
        !params.skip_quast
        script:
        """
	quast.py ${assembly} --threads ${params.threads} -o \$PWD
	sed "s/contigs\$/${sample}/" report.tsv > ${sample}_report.tsv
        rm transposed_report.tsv report.tsv
	cp .command.log quast.log
        """
}

process summary_quast {
	publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*tsv'
	input:
		path(quast_files)
	output:
		path("4_quast_report.tsv"), emit: quast_summary	
	script:
	"""
	for file in `ls *report.tsv`; do cut -f2 \$file > \$file.tmp.txt; cut -f1 \$file > rownames.txt; done
	paste rownames.txt *tmp.txt > 4_quast_report.tsv
	"""
}

process checkm {
        cpus "${params.threads}"
        tag "${sample}"
        label "cpu"
        label "high_memory"
        publishDir "$params.outdir/$sample/5_checkm",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/5_checkm",  mode: 'copy', pattern: '*tsv'
        input:
                tuple val(sample), path(assembly)
        output:
                path("checkm.log")
		path("*checkm_lineage_wf_results.tsv"), emit: checkm_results
        when:
        !params.skip_checkm
        script:
        """
        export CHECKM_DATA_PATH=${params.checkm_db}
        checkm data setRoot ${params.checkm_db}
        checkm lineage_wf --reduced_tree `dirname ${assembly}` \$PWD --threads ${params.threads} --pplacer_threads ${params.threads} --tab_table -f checkm_lineage_wf_results.tsv -x fa
        mv checkm_lineage_wf_results.tsv ${sample}_checkm_lineage_wf_results.tsv
	cp .command.log checkm.log
        """
}

process summary_checkm {
	publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*tsv'
	input:
		path(checkm_files)
	output:
		path("5_checkm_lineage_wf_results.tsv"), emit: checkm_summary
	script:
	"""
	echo -e  sampleID\\\tMarker_lineage\\\tNbGenomes\\\tNbMarkers\\\tNbMarkerSets\\\t0\\\t1\\\t2\\\t3\\\t4\\\t5+\\\tCompleteness\\\tContamination\\\tStrain_heterogeneity > header_checkm
	for file in `ls *checkm_lineage_wf_results.tsv`; do fileName=\$(basename \$file); sample=\${fileName%%_checkm_lineage_wf_results.tsv}; grep -v Bin \$file | sed s/^contigs/\${sample}/  >> 5_checkm_lineage_wf_results.tsv.tmp; done
	cat header_checkm 5_checkm_lineage_wf_results.tsv.tmp > 5_checkm_lineage_wf_results.tsv
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
        publishDir "$params.outdir/$sample/6_kraken",  mode: 'copy', pattern: '*tsv'
        input:
                tuple val(sample), path(reads1_trimmed), path(reads2_trimmed), path(kraken2_report), path(kraken_tsv)
        output:
                path("bracken.log")
		path("*bracken_species.tsv"), emit: bracken_results
        when:
        !params.skip_kraken
        script:
        """
        bracken -d ${params.kraken_db} -i ${kraken2_report} -o bracken_species.tsv  -w bracken_report.txt -r 100 -l S -t ${params.bracken_threads}
	mv bracken_species.tsv ${sample}_bracken_species.tsv
        cp .command.log bracken.log
        """
}

process summary_bracken {
        publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*tsv'
        input:
                path(bracken_files)
        output:
                tuple path("6_bracken_pasteurella_multocida_species_abundance.tsv") ,path("6_bracken_most_abundant_species.tsv"), emit: bracken_summary
        script:
        """
        echo -e sampleID\\\tname\\\ttaxonomy_id\\\ttaxonomy_lvl\\\tkraken_assigned_reads\\\tadded_reads\\\tnew_est_reads\\\tfraction_total_reads > header_bracken
        for file in `ls *_bracken_species.tsv`; do fileName=\$(basename \$file); sample=\${fileName%%_bracken_species.tsv}; grep Pasteurella \$file | grep multocida | sed s/^/\${sample}\\\t/  >> 6_bracken_pasteurella_multocida_species_abundance.tsv.tmp; done
        cat header_bracken 6_bracken_pasteurella_multocida_species_abundance.tsv.tmp > 6_bracken_pasteurella_multocida_species_abundance.tsv
	for file in `ls *_bracken_species.tsv`; do fileName=\$(basename \$file); sample=\${fileName%%_bracken_species.tsv}; grep -v taxonomy_id \$file | sort -t\$'\t' -k7gr | head -1 | sed s/^/\${sample}\\\t/  >> 6_bracken_most_abundant_species.tsv.tmp; done
	cat header_bracken 6_bracken_most_abundant_species.tsv.tmp > 6_bracken_most_abundant_species.tsv
        """
}

process kaptive3 {
        cpus "${params.threads}"
        tag "${sample}"
        label "cpu"
        publishDir "$params.outdir/$sample/7_kaptive_v3",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/7_kaptive_v3",  mode: 'copy', pattern: '*fna'
	publishDir "$params.outdir/$sample/7_kaptive_v3",  mode: 'copy', pattern: '*tsv'
        input:
                tuple val(sample), path(assembly)
        output:
		tuple val(sample), path("*kaptive_results.tsv"), emit: kaptive_results
		path("*kaptive_results.tsv"),  emit: kaptive_tsv
		path("*fna")
                path("kaptive_v3.log")
        when:
        !params.skip_kaptive3
        script:
        """
	kaptive assembly ${params.kaptive_db_9lps} ${assembly} -f \$PWD -o kaptive_results.tsv
	mv kaptive_results.tsv ${sample}_kaptive_results.tsv
	sed s/contigs/${sample}/ contigs_kaptive_results.fna > ${sample}_contigs_kaptive_results.fna
	rm contigs_kaptive_results.fna
	cp .command.log kaptive_v3.log
        """
}

process summary_kaptive {
        publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*tsv'
	input:
		path(kaptive_files)
	output:
		path("7_kaptive_results.tsv"), emit: kaptive_summary
	script:
	"""
	echo -e sampleID\\\tBest match locus\\\tBest match type\\\tMatch confidence\\\tProblems\\\tIdentity\\\tCoverage\\\tLength discrepancy\\\tExpected genes in locus\\\tExpected genes in locus, details\\\tMissing expected genes\\\tOther genes in locus\\\tOther genes in locus, details\\\tExpected genes outside locus\\\tExpected genes outside locus, details\\\tOther genes outside locus\\\tOther genes outside locus, details\\\tTruncated genes, details\\\tExtra genes, details >  header_kaptive3
	for file in `ls *_kaptive_results.tsv`; do fileName=\$(basename \$file); sample=\${fileName%%_kaptive_results.tsv}; grep -v Assembly \$file | sed s/contigs/\${sample}/  >> 7_kaptive_results.tsv.tmp; done
	cat header_kaptive3 7_kaptive_results.tsv.tmp > 7_kaptive_results.tsv
	"""
}

process snippy {
        cpus "${params.snippy_threads}"
        tag "${sample}"
        label "cpu"
        publishDir "$params.outdir/$sample/8_snippy",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
        publishDir "$params.outdir/$sample/8_snippy",  mode: 'copy', pattern: 'snps*', saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/8_snippy",  mode: 'copy', pattern: '*tab'
	input:
                tuple val(sample), path(reads1), path(reads2), path(reads1_trimmed), path(reads2_trimmed), path(kaptive_report)
	output:
                tuple val(sample), path("*snps.tab"), path("*snps.high_impact.tab"), path("snps.raw.vcf"), path("snps.filt.vcf"),  path("snps.bam"), path("snps.bam.bai"),  emit: snippy_results
		path("snippy.log")
		tuple path("*snps.tab"), path("*snps.high_impact.tab"), emit: snippy_impact_tab
        when:
        !params.skip_snippy
        shell:
        '''
	locus=`tail -1 !{kaptive_report} | cut -f3`
	ref_gb=`grep ${locus:0:2} !{params.reference_LPS} | cut -f2`
	snippy --cpus !{params.snippy_threads} --force --outdir \$PWD --ref $ref_gb --R1 !{reads1_trimmed} --R2 !{reads2_trimmed}
        egrep "^CHROM|frameshift_variant|stop_gained" snps.tab > snps.high_impact.tab
	mv snps.high_impact.tab !{sample}_snps.high_impact.tab
	mv snps.tab !{sample}_snps.tab
	cp .command.log snippy.log
        '''
}

process report {
	publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*tsv'	
	input:
		path(snippy_files)
	output:
		tuple path("8_snippy_snps.tsv"), path("8_snippy_snps.high_impact.tsv"), path("10_genotype_report.tsv"), emit: genotype_report	
	script:
	"""
	echo -e sampleID\\\tCHROM\\\tPOS\\\tTYPE\\\tREF\\\tALT\\\tEVIDENCE\\\tFTYPE\\\tSTRAND\\\tNT_POS\\\tAA_POS\\\tEFFECT\\\tLOCUS_TAG\\\tGENE\\\tPRODUCT > header_snippy
	for file in `ls *_snps.high_impact.tab`; do fileName=\$(basename \$file); sample=\${fileName%%_snps.high_impact.tab}; grep -v EVIDENCE \$file | sed s/^/\${sample}\\\t/  >> 8_snippy_snps.high_impact.tsv.tmp; done
	cat header_snippy 8_snippy_snps.high_impact.tsv.tmp > 8_snippy_snps.high_impact.tsv
	for file in `ls *_snps.tab`; do fileName=\$(basename \$file); sample=\${fileName%%_snps.tab}; grep -v EVIDENCE \$file | sed s/^/\${sample}\\\t/  >> 8_snippy_snps.tsv.tmp; done
	cat header_snippy 8_snippy_snps.tsv.tmp > 8_snippy_snps.tsv
	touch 10_genotype_report.tsv
	while IFS=\$'\t' read sample chrom pos type ref alt evidence ftype strand nt_pos aa_pos effect locus_tag gene product; do
		while IFS=\$'\t' read db_LPStype db_genotype db_isolate db_chrom db_pos db_type db_ref db_alt db_gene; do 
			if [[ \$chrom == \$db_chrom && \$pos == \$db_pos && \$ref == \$db_ref && \$alt == \$db_alt ]]; then
				if [[ \$sample != "sampleID" ]]; then
					echo "sample" \$sample": found genotype" \$db_genotype "with" \$db_type "(similar to isolate" \$db_isolate")" >> 10_genotype_report.tsv
				fi
			fi
		done < ${params.genotype_db}
	done < 8_snippy_snps.tsv
	"""
}

process mlst {
        cpus "${params.threads}"
        tag "${sample}"
        label "cpu"
        publishDir "$params.outdir/$sample/9_mlst",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/9_mlst",  mode: 'copy', pattern: '*csv'
        input:
                tuple val(sample), path(assembly)
        output:
                path("mlst.log")
		path("*_mlst.csv"), emit: mlst_results
        when:
        !params.skip_mlst
        script:
        """
	mlst --scheme ${params.mlst_scheme} ${assembly} --quiet --csv --threads ${params.threads} > mlst.csv
        mv mlst.csv ${sample}_mlst.csv
	cp .command.log mlst.log
        """
}

process summary_mlst {
	publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*csv'
	input:
		path(mlst_files)
	output:
		path("9_mlst.csv"), emit: mlst_summary
	script:
	"""
	for file in `ls *_mlst.csv`; do fileName=\$(basename \$file); sample=\${fileName%%_mlst.csv};  sed s/^/\${sample}_/ \$file >> 9_mlst.csv; done
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
		if (!params.skip_summary_fastqc) {
			summary_fastqc(fastqc.out.fastqc_zip.collect())
		}
	}
	shovill(ch_samplesheet_illumina)
	if (!params.skip_quast) {
		quast(shovill.out.assembly_out)
		summary_quast(quast.out.quast_results.collect())
	}
	if (!params.skip_kraken) {
		kraken(fastp.out.trimmed_fastq)
		bracken(kraken.out.kraken_results)
		summary_bracken(bracken.out.bracken_results.collect())
	}
	if (!params.skip_checkm) {
		checkm(shovill.out.assembly_out)
		summary_checkm(checkm.out.checkm_results.collect())
	}
	if (!params.skip_kaptive3) {
		kaptive3(shovill.out.assembly_out)
		summary_kaptive(kaptive3.out.kaptive_tsv.collect())
		if (!params.skip_snippy) {
			snippy(fastp.out.trimmed_fastq.join(kaptive3.out.kaptive_results))
			report(snippy.out.snippy_impact_tab.collect())
		}
	}
	if (!params.skip_mlst) {
		mlst(shovill.out.assembly_out)
		summary_mlst(mlst.out.mlst_results.collect())
	}
}
