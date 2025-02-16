#!/bin/sh
# SGE options (lines prefixed with #$)
#$ -N runExtractRegions.sh
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l h_vmem=8G
#$ -pe sharedmem 4 
#$ -e e_extbwa
#$ -o o_extbwa


#Initialise the environment modules
. /etc/profile.d/modules.sh
module load igmm/apps/samtools/1.16.1

#initialise conda 
source /exports/applications/apps/community/roslin/conda/4.9.1/etc/profile.d/conda.sh
conda activate mapstats

#define dirs
target_dir=/exports/eddie/scratch/mmarr3/fsil_map_bwa

# Generate list of chromosome names from CM001378.3 to CM001395.3
chromosomes=""
for i in $(seq 1378 1395); do
    chromosomes+="CM$(printf "%06d" $i).3 "
done

# Extract specified chromosomes without listing them explicitly
samtools view -b "$target_dir/FS091.merged.map.RG.unique.bam" $chromosomes > "$target_dir/FS091_final.bam"

# Sort and index the new BAM file
samtools sort "$target_dir/FS091_final.bam" -o "$target_dir/FS091_final.sort.bam"
samtools index "$target_dir/FS091_final.sort.bam"

# Get mapping statistics
samtools flagstat "$target_dir/FS091_final.sort.bam" > "$target_dir/FS091_final.sort.flagstat"
bedtools genomecov -ibam "$target_dir/FS091_final.sort.bam" > "$target_dir/FS091_final.sort.cov.txt"
mosdepth -n -t 4 -Q30 "$target_dir/FS091_final.sort" "$target_dir/FS091_final.sort.bam"
