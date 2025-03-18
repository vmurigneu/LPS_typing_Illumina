#!/bin/bash
#
#SBATCH --time=30:00:00
#SBATCH --job-name=Illumina_LPS_pipeline
#SBATCH --output=./s%j_job.pipeline_Illumina_L6_ONT_test.out
#SBATCH --error=./s%j_job.pipeline_Illumina_L6_ONT_test.error
#SBATCH --account=a_uqds
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1

module load nextflow/23.04.2

#Directory containing the nextflow.config file and the main.nf script
dir=/scratch/project/uqds/valentine/LPS/PIPELINE
cd ${dir}

#Samplesheet file
samplesheet=${dir}/samplesheet/samples_L6_test.csv

#Directory that will be created to contain the output files
out_dir=${dir}/results_test

#Bunya Slurm account 
slurm_account='a_uqds'

#Directory containing the Illumina fastq files
fqdir=${dir}/fastq
#cp /QRISdata/Q2313/AGRF-raw-Illumina/AGRF_CAGRF24030272-2_22VH7WLT3/*fastq.gz ${fqdir}

#Run the pipeline (-resume is used to restart the pipeline if it fails)
nextflow main.nf --outdir ${out_dir} --fqdir ${fqdir} --samplesheet ${samplesheet} -resume

# Collect main results in one file for all samples
cd ${out_dir}
#5_checkm
echo -e  sampleID\\tMarker_lineage\\tNbGenomes\\tNbMarkers\\tNbMarkerSets\\t0\\t1\\t2\\t3\\t4\\t5+\\tCompleteness\\tContamination\\tStrain_heterogeneity > header_checkm
for file in `ls ${out_dir}/*/5_checkm/*checkm_lineage_wf_results.tsv`; do fileName=$(basename $file); sample=${fileName%%_checkm_lineage_wf_results.tsv}; grep -v Bin $file | sed s/^/${sample}_/  >> 5_checkm_lineage_wf_results.tsv.tmp; done
cat header_checkm 5_checkm_lineage_wf_results.tsv.tmp > 5_checkm_lineage_wf_results.tsv

#6_kraken
echo -e sampleID\\tname\\ttaxonomy_id\\ttaxonomy_lvl\\tkraken_assigned_reads\\tadded_reads\\tnew_est_reads\\tfraction_total_reads > header_bracken
for file in `ls ${out_dir}/*/6_kraken/*_bracken_species.tsv`; do fileName=$(basename $file); sample=${fileName%%_bracken_species.tsv}; grep Pasteurella $file | grep multocida | sed s/^/${sample}\\t/  >> 6_bracken_species.tsv.tmp; done
cat header_bracken 6_bracken_species.tsv.tmp > 6_bracken_species.tsv

#7_kaptive
echo -e sampleID\\tBest match locus\\tBest match type\\tMatch confidence\\tProblems\\tIdentity\\tCoverage\\tLength discrepancy\\tExpected genes in locus\\tExpected genes in locus, details\\tMissing expected genes\\tOther genes in locus\\tOther genes in locus, details\\tExpected genes outside locus\\tExpected genes outside locus, details\\tOther genes outside locus\\tOther genes outside locus, details\\tTruncated genes, details\\tExtra genes, details >  header_kaptive3
for file in `ls ${out_dir}/*/7_kaptive_v3/*_kaptive_results.tsv`; do fileName=$(basename $file); sample=${fileName%%_kaptive_results.tsv}; grep -v Assembly $file | sed s/^/${sample}_/  >> 7_kaptive_v3_9lps_results.tsv.tmp; done
cat header_kaptive3 7_kaptive_v3_9lps_results.tsv.tmp > 7_kaptive_v3_9lps_results.tsv

#8_snippy
echo -e sampleID\\tCHROM\\tPOS\\tTYPE\\tREF\\tALT\\tEVIDENCE\\tFTYPE\\tSTRAND\\tNT_POS\\tAA_POS\\tEFFECT\\tLOCUS_TAG\\tGENE\\tPRODUCT > header_snippy
for file in `ls ${out_dir}/*/8_snippy/*_snps.high_impact.tab`; do fileName=$(basename $file); sample=${fileName%%_snps.high_impact.tab}; grep -v EVIDENCE $file | sed s/^/${sample}\\t/  >> 8_snippy_snps.high_impact.tsv.tmp; done
cat header_snippy 8_snippy_snps.high_impact.tsv.tmp > 8_snippy_snps.high_impact.tsv

#9_mlst
for file in `ls ${out_dir}/*/9_mlst/*_mlst_pmultocida_rirdc.csv`; do fileName=$(basename $file); sample=${fileName%%_mlst_pmultocida_rirdc.csv};  sed s/^/${sample}_/ $file >> 9_mlst_pmultocida_rirdc.csv; done

#clean folder $out_dir
rm $out_dir/*tmp  $out_dir/header*


