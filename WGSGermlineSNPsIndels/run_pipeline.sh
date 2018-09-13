#!/bin/bash
#SBATCH -p long
#SBATCH -N 1
#SBATCH --tasks-per-node=16
#SBATCH -t 3-00:00:00
#SBATCH --mem=60000
#SBATCH -o slurm-wgs-gatk-%j.out
#SBATCH -e slurm-wgs-gatk-%j.err

source /etc/profile.d/profile.modules.sh
module load java-jre/1.8.0_77
module load gatk/4.0.6.0
module load bwa/0.7.15
module load picard/2.10.10
module load python/2.7.14
module load python/3.4.5
module load proj/4.9.3/gcc gdal/2.2.4/gcc gcc/6.1.0 tiff/4.0.9/gcc pcre/3.38/gcc curl/7.40/gcc xz/5.2.2/gcc R/3.5-GA/gcc
module load samtools/1.7
module load htslib/1.3.1
module load java-jre/1.8.0_77 qualimap/2.2.1

java -Dconfig.file=/lustre/jlorsal/genomics/pipelines/wgs/gatk4/FDA-Challenge/pipeline/backend-slurm.conf \
-jar /lustre/amunoz/tools/cromwell-30.2.jar \
run /lustre/jlorsal/genomics/pipelines/wgs/gatk4/FDA-Challenge/pipeline/WholeGenomeSequencingGATK4.wdl \
-i /lustre/jlorsal/genomics/pipelines/wgs/gatk4/FDA-Challenge/pipeline/TruthChallengeFDA.inputs.json \
-m metadata.json