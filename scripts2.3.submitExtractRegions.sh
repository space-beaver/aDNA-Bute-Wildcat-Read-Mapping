#!/bin/sh
# SGE options (lines prefixed with #$)
#$ -N runExtractRegions.sh
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l h_vmem=8G
#$ -pe sharedmem 4 
#$ -e e_extbowtiw
#$ -o o_extbowtie 


#Initialise the environment modules
. /etc/profile.d/modules.sh
module load igmm/apps/samtools/1.16.1

#initialise conda 
source /exports/applications/apps/community/roslin/conda/4.9.1/etc/profile.d/conda.sh
conda activate mapstats

#define dirs
target_dir=/exports/eddie/scratch/mmarr3/fsil_map_bowtie

#extract
samtools view -b $target_dir/FS091.merged.bowtie.RG.unique.bam CM001378.3 CM001379.3 CM001380.3 CM001381.3 CM001382.3 CM001383.3 CM001384.3 CM001385.3 CM001386.3 CM001387.3 CM001388.3 CM001389.3 CM001390.3 CM001391.3 CM001392.3 CM001393.3 CM001394.3 CM001395.3 > $target_dir/FS091_final.bowtie.bam 

# Sort and Index the new BAM file
samtools sort $target_dir/FS091_final.bowtie.bam > $target_dir/FS091_final.sort.bowtie.bam
samtools index  $target_dir/FS091_final.sort.bowtie.bam

#get stats 
samtools flagstat $target_dir/FS091_final.sort.bowtie.bam  > $target_dir/FS091_final.sort.bowtie.flagstat
bedtools genomecov -ibam $target_dir/FS091_final.sort.bowtie.bam  > $target_dir/FS091_final.sort.bowtie.cov.txt
mosdepth -n -t 4 -Q30 $target_dir/FS091_final.sort.bowtie $target_dir/FS091_final.sort.bowtie.bam
