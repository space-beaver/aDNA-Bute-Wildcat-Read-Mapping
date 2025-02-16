#!/bin/sh
# SGE options (lines prefixed with #$)
#$ -N runMergeBamsBT.sh
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l h_vmem=8G
#$ -pe sharedmem 4 
#$ -e e_mergebtie
#$ -o o_mergebtie

#Jobscript to merge bam files  - need to make into an array at some point. 

#Initialise the environment modules
. /etc/profile.d/modules.sh
module load igmm/apps/samtools/1.16.1
module load roslin/bwa/2.1.0


#initialise conda 
source /exports/applications/apps/community/roslin/conda/4.9.1/etc/profile.d/conda.sh
conda activate mapstats

#define dirs
target_dir=/exports/eddie/scratch/mmarr3/fsil_map_bowtie 
picard=/exports/cmvm/eddie/eb/groups/ogden_grp/software/picard/picard.jar

#merge
java -Xmx10g -jar $picard MergeSamFiles I=$target_dir/091SepRun.bowtie.sort.RG.bam I=$target_dir/091x.bowtie.sort.RG.bam I=$target_dir/M01.bowtie.sort.RG.bam I=$target_dir/M02.bowtie.sort.RG.bam I=$target_dir/M03.bowtie.sort.RG.bam I=$target_dir/M04.bowtie.sort.RG.bam I=$target_dir/MM091fs.bowtie.sort.RG.bam I=$target_dir/MM091fs.bowtie.sort.RG.bam O=$target_dir/FS091.merged.bowtie.RG.bam TMP_DIR=$tmp_dir

#sort 
samtools sort $target_dir/FS091.merged.bowtie.RG.bam > $target_dir/FS091.merged.bowtie.RG.sorted.bam
samtools index $target_dir/FS091.merged.bowtie.RG.sorted.bam

#mark dups
java -Xmx4g -jar /exports/cmvm/eddie/eb/groups/ogden_grp/software/picard/picard.jar MarkDuplicates \
     I=$target_dir/FS091.merged.bowtie.RG.bam  \
     O=$target_dir/FS091.merged.bowtie.RG.dedup.bam \
     M=$target_dir/FS091.bowtie.metrics.txt \
     REMOVE_DUPLICATES=false \
     VALIDATION_STRINGENCY=SILENT
	   TMP_DIR=tmp

# Create BAM with only duplicate reads
samtools view -b -f 1024 $target_dir/FS091.merged.bowtie.RG.dedup.bam > $target_dir/FS091.merged.bowtie.RG.dups.bam 
samtools index $target_dir/FS091.merged.bowtie.RG.dups.bam 

# Create BAM with duplicates removed
samtools view -b -F 1024 $target_dir/FS091.merged.bowtie.RG.dedup.bam > $target_dir/FS091.merged.bowtie.RG.unique.bam
samtools index $target_dir/FS091.merged.bowtie.RG.unique.bam

#get stats 
samtools flagstat $target_dir/FS091.merged.bowtie.RG.unique.bam  > $target_dir/FS091.bowtie.merged.unique.flagstat
bedtools genomecov -ibam $target_dir/FS091.merged.bowtie.RG.unique.bam  > $target_dir/FS091.bowtie.merged.unique.cov.txt
mosdepth -n -t 4 -Q30 $target_dir/FS091.bowtie.merged $target_dir/FS091.merged.bowtie.RG.unique.bam
