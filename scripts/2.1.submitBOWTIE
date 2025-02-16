#!/bin/sh
# Grid Engine options
#$ -N runBOWTIEaln.sh
#$ -hold_jid runAdapRemoval.sh
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l h_vmem=16G
#$ -pe sharedmem 4
#$ -t 1-8
#$ -o o_BOWTIE
#$ -e e_BOWTIE

##script to map illumina PE reads from ancient DNA data to a reference genome
#to run on UoE Eddie compute cluster as an array 
#written by Melissa M. Marr mmarr3@ed.ac.uk 20/07/2024 
#assumes no merging of bam files 

#Initialise the environment modules
. /etc/profile.d/modules.sh
module load igmm/apps/samtools/1.16.1
module load igmm/apps/bowtie/2.5.3 

#initialise conda 
source /exports/applications/apps/community/roslin/conda/4.9.1/etc/profile.d/conda.sh
conda activate mapstats 

#define dirs
refgen=/exports/eddie/scratch/mmarr3/felCat_9.2_X
sample_list="/exports/eddie/scratch/mmarr3/fsil_fastqs/bute_cat_samples.txt"
in_dir=/exports/eddie/scratch/mmarr3/fsil_trim
out_dir=/exports/eddie/scratch/mmarr3/fsil_map_bowtie
fastqs=/exports/eddie/scratch/mmarr3/fsil_fastqs
picard=/exports/cmvm/eddie/eb/groups/ogden_grp/software/picard/picard.jar

#get sample lists
base=`sed -n "$SGE_TASK_ID"p $sample_list | awk '{print $1}'`

#process
echo Processing sample: ${base} on $HOSTNAME

#align 
bowtie2 -x $refgen/felcat9.2_X -U $in_dir/${base}_collapsed.fastq.gz --local -N 1 --mp 4 | samtools view -bF 4 -q 30 | samtools sort -o $out_dir/${base}.bowtie.sort.bam

#index
samtools index -@ 4 $out_dir/${base}.bowtie.sort.bam

#add read group info 
file=`sed -n "$SGE_TASK_ID"p $sample_list | awk '{print $2}'` 
flow_cell=`sed -n "$SGE_TASK_ID"p $sample_list | awk '{print $4}'` 

echo $file 
echo $flow_cell 


java -Xmx4g -jar /exports/cmvm/eddie/eb/groups/ogden_grp/software/picard/picard.jar AddOrReplaceReadGroups \
     I=$out_dir/${base}.bowtie.sort.bam \
     O=$out_dir/${base}.bowtie.sort.RG.bam \
     RGID=$base \
     RGPL=illumina \
     RGLB=lib1 \
     RGPU=$flow_cell \
     RGSM=$base \
     VALIDATION_STRINGENCY=SILENT \
     SORT_ORDER=coordinate \
     TMP_DIR=tmp


#mark dups
java -Xmx4g -jar /exports/cmvm/eddie/eb/groups/ogden_grp/software/picard/picard.jar MarkDuplicates \
     I=$out_dir/${base}.bowtie.sort.RG.bam \
     O=$out_dir/${base}.bowtie.sort.RG.dedup.bam \
     M=$out_dir/${base}.metrics.txt \
     REMOVE_DUPLICATES=false \
     VALIDATION_STRINGENCY=SILENT
	   TMP_DIR=tmp

# Create BAM with only duplicate reads
samtools view -b -f 1024 $out_dir/${base}.bowtie.sort.RG.dedup.bam > $out_dir/${base}.bowtie.sort.RG.dups.bam
samtools index $out_dir/${base}.bowtie.sort.RG.dups.bam

# Create BAM with duplicates removed
samtools view -b -F 1024 $out_dir/${base}.bowtie.sort.RG.dedup.bam > $out_dir/${base}.bowtie.sort.RG.unique.bam
samtools index $out_dir/${base}.bowtie.sort.RG.unique.bam

#get stats 
samtools flagstat $out_dir/${base}.bowtie.sort.RG.dups.bam  > $out_dir/${base}.bowtie.dups.flagstat.txt
samtools flagstat $out_dir/${base}.bowtie.sort.RG.unique.bam  > $out_dir/${base}.bowtie.unique.flagstat.txt
bedtools genomecov -ibam $out_dir/${base}.bowtie.sort.RG.unique.bam  > $out_dir/${base}.bowtie.unique.cov.txt
mosdepth -n -t 4 -Q30 $out_dir/${base} $out_dir/${base}.bowtie.sort.RG.unique.bam

#multiqc $in_dir -o $out_dir 
