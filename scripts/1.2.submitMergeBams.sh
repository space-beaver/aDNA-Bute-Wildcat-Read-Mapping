#!/bin/sh
# SGE options (lines prefixed with #$)
#$ -N runMergeBams.sh
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l h_vmem=8G
#$ -pe sharedmem 4 
#$ -e e_merge
#$ -o o_merge

#Jobscript to merge bam files  - need to make into an array at some point. 

#Initialise the environment modules
. /etc/profile.d/modules.sh
module load igmm/apps/samtools/1.16.1
module load roslin/bwa/2.1.0


#initialise conda 
source /exports/applications/apps/community/roslin/conda/4.9.1/etc/profile.d/conda.sh
conda activate mapstats

#define dirs
target_dir=/exports/eddie/scratch/mmarr3/fsil_map_bwa
picard=/exports/cmvm/eddie/eb/groups/ogden_grp/software/picard/picard.jar

#merge
java -Xmx10g -jar $picard MergeSamFiles I=$target_dir/091SepRun.map.sort.RG.bam I=$target_dir/091x.map.sort.RG.bam I=$target_dir/M01.map.sort.RG.bam I=$target_dir/M02.map.sort.RG.bam I=$target_dir/M03.map.sort.RG.bam I=$target_dir/M04.map.sort.RG.bam I=$target_dir/MM091fs.map.sort.RG.bam I=$target_dir/MM091fs.map.sort.RG.bam O=$target_dir/FS091.merged.map.RG.bam TMP_DIR=$tmp_dir

#sort 
samtools sort $target_dir/FS091.merged.map.RG.bam > $target_dir/FS091.merged.map.RG.sorted.bam
samtools index $target_dir/FS091.merged.map.RG.sorted.bam

#mark dups
java -Xmx4g -jar /exports/cmvm/eddie/eb/groups/ogden_grp/software/picard/picard.jar MarkDuplicates \
     I=$target_dir/FS091.merged.map.RG.bam  \
     O=$target_dir/FS091.merged.map.RG.dedup.bam \
     M=$target_dir/FS091.metrics.txt \
     REMOVE_DUPLICATES=false \
     VALIDATION_STRINGENCY=SILENT
	   TMP_DIR=tmp

# Create BAM with only duplicate reads
samtools view -b -f 1024 $target_dir/FS091.merged.map.RG.dedup.bam > $target_dir/FS091.merged.map.RG.dups.bam 
samtools index $target_dir/FS091.merged.map.RG.dups.bam 

# Create BAM with duplicates removed
samtools view -b -F 1024 $target_dir/FS091.merged.map.RG.dedup.bam > $target_dir/FS091.merged.map.RG.unique.bam
samtools index $target_dir/FS091.merged.map.RG.unique.bam

#get stats 
samtools flagstat $target_dir/FS091.merged.map.RG.unique.bam  > $target_dir/FS091.merged.unique.flagstat
bedtools genomecov -ibam $target_dir/FS091.merged.map.RG.unique.bam  > $target_dir/FS091.merged.unique.cov.txt
mosdepth -n -t 4 -Q30 $targer_dir/FS091.merged $target_dir/FS091.merged.map.RG.unique.bam
