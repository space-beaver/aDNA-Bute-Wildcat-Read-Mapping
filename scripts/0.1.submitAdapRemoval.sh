#!/bin/sh
# Grid Engine options
#$ -N runAdapRemoval
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l h_vmem=16G
#$ -pe sharedmem 4
#$ -t 1-8
#$ -l rl9=true
#$ -o e_adap_files
#$ -e o_adap_files

#script to trim, qc and merge illumina PE reads from ancient DNA data
#to run on UoE Eddie compute cluster as an array 
#written by Melissa M. Marr mmarr3@ed.ac.uk 20/07/2024 

#Initialise the environment modules
. /etc/profile.d/modules.sh
source /exports/applications/apps/community/roslin/conda/4.9.1/etc/profile.d/conda.sh
conda activate adapter_removal 

#define dirs
target_dir=/exports/eddie/scratch/mmarr3/fsil_fastqs
out_dir=/exports/eddie/scratch/mmarr3/fsil_trim
sample_list="/exports/eddie/scratch/mmarr3/fsil_fastqs/bute_cat_samples.txt"

#filelist 
base=`sed -n "$SGE_TASK_ID"p $sample_list | awk '{print $1}'`
R1=`sed -n "$SGE_TASK_ID"p $sample_list | awk '{print $2}'`
R2=`sed -n "$SGE_TASK_ID"p $sample_list | awk '{print $3}'`


#process
echo Processing sample: ${base} on $HOSTNAME
echo Processing $R1
echo Processing $R2

#trim 
AdapterRemoval --file1 $target_dir/${base}_R1_001.fastq.gz --file2 $target_dir/${base}_R2_001.fastq.gz \
--threads 4 \
--gzip --mm 3 \
--minlength 30  \
--basename $out_dir/${base} \
--trimns \
--trimqualities \
--collapse \
--output1 $out_dir/${base}_R1.only.fastq.gz \
--output2  $out_dir/${base}_R2.only.fastq.gz \
--outputcollapsedtruncated $out_dir/${base}_truncated.fastq.gz  \
--outputcollapsed $out_dir/${base}_collapsed.fastq.gz \
--discarded $out_dir/${base}_discarded.fastq.gz \
--settings $out_dir/${base}.settings \
--singleton $out_dir/${base}_singleton.fastq.gz \
--adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
--adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT 

