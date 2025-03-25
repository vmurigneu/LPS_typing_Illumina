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

