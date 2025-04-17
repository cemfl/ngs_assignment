#!/bin/bash
set -euo pipefail # exits when errors are found

#=========================================================
# This script runs my entire ngs pipeline.
# Run using 'bash ngs_pipelin.sh' in the terminal
# with these files in the folder:

# NexteraPE-PE.fa 
# environment.yml
# annovar.latest.tar.gz 
# ngs_pipeline.sh

# ======================================================

# CONDA ENVIRONMENT
# The conda environment for the project is set up below.

# ======================================================

if [[ ! -d "$HOME/miniconda3/envs/bioinfo" ]]; then
  conda env create -f environment.yml
else
  echo "Environment already installed, intialising environment...."
fi

eval "$(conda shell.bash hook)"
conda activate bioinfo2

#=======================================================

#DIRECTORIES

#Directories for the pipline are created below.

#========================================================

mkdir -p ngs_pipeline/data/{untrimmed_fastq,trimmed_fastq,reference,aligned,beds} \
         ngs_pipeline/logs \
         ngs_pipeline/meta \
         ngs_pipeline/results/{depth_coverage,fastqc,flagstats,idxstats,insert_size,marked_duplicates} \

#==============================================================

#VARIABLE-ASSIGNMENT

#Variables for file directories and filenames are set out below.
#This is for easy adjustment of the paths and files used,
#as well as an easy to read and adjust pipeline.

#==============================================================

# Sets the number of threads for the pipeline tools, where applicable.
threads=4

# Data directory path variables
data_dir="ngs_pipeline/data/"
untrimmed_dir="${data_dir}untrimmed_fastq/"
trimmed_dir="${data_dir}trimmed_fastq/"
aligned_dir="${data_dir}aligned/"
reference_dir="${data_dir}reference/"
bed_dir="${data_dir}beds/"

# Results directory path variables
results_dir="ngs_pipeline/results/"
dep_cov="${results_dir}depth_coverage/"
fastqc="${results_dir}fastqc/"
flagstats="${results_dir}flagstats/"
idxstats="${results_dir}idxstats/"
insert_size="${results_dir}insert_size/"
mark_dup="${results_dir}marked_duplicates/"

# Sets the variables for the data files that will be used.

# Sequence files
rep1url="https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz"
rep2url="https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz"
filename=$(basename "$rep1url")
sample_id=$(echo "$filename" | cut -d'.' -f1) 
rep1="${untrimmed_dir}${sample_id}.R1.fastq.gz"
rep2="${untrimmed_dir}${sample_id}.R2.fastq.gz"
trimmed_R_1P="${trimmed_dir}${sample_id}_trimmed_R_1P.fastq"
trimmed_R_2P="${trimmed_dir}${sample_id}_trimmed_R_2P.fastq"

# Bedfiles
bedfilename="annotation"
bedfile="${bed_dir}${bedfilename}.bed"
bedfile_stripped="${bed_dir}${bedfilename}_stripped.bed"

# Reference files
ref_zip="${reference_dir}hg19.fa.gz"
ref_unzip="${reference_dir}hg19.fa"

# Depth_coverage
dep_cov_prefix="${dep_cov}${sample_id}_coverage"
dep_cov_file="${dep_cov}${sample_id}_coverage.regions.bed.gz"

# SAM and BAM files
samfile="${aligned_dir}${sample_id}.sam"
bamfile="${aligned_dir}${sample_id}.bam"
sorted_bamfile="${aligned_dir}${sample_id}_sorted.bam"
marked_bamfile="${aligned_dir}${sample_id}_sorted_marked.bam"
filtered_bamfile="${aligned_dir}${sample_id}_sorted_filtered.bam"

# VCF files
vcfFile_unzip="${results_dir}${sample_id}.vcf"
vcfFile_zip="${vcfFile_unzip}.gz"
vcfFile_filtered="${results_dir}${sample_id}_filtered.vcf"
vcfFile_filtered_zip="${vcfFile_filtered}.gz"
vcfFile_annotated="${results_dir}${sample_id}_filtered_annotation.vcf"
vcfFile_annotated_zip="${vcfFile_annotated}.gz"

# Annovar files
avinput="${results_dir}${sample_id}_filtered_annotation.avinput"

#==============================================================

# FILES - DOWNLOAD AND RENAME

#==============================================================


# Downloads the data files you are using to the directory they will be used from.
if [[ ! -f $rep2 ]]; then
  wget -P "$untrimmed_dir" "$rep1url" "$rep2url"
else
  echo "Sequences already downloaded - skipping...."
fi

# Ensures file end in .gz and not .qz or .cz etc. 
if [[ ! -f $rep2 ]]; then
  for f in ngs_pipeline/data/untrimmed_fastq/*.fastq.*; do
      newname="${f%%.fastq.*}.fastq.gz"
      mv "$f" "$newname"
      echo "Renamed: \"$f\" → \"$newname\""
  done
else
  echo "File extensions are correct, skipping...."
fi

# Downloads the bed file into the directory it will be used from.
if [[ ! -f $bedfile ]]; then
  wget -P $bed_dir "https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed"
else
  echo "BED file already downloaded - skipping...."
fi

# Downloads the reference file with chr1, chr2 etc. to fit with the bed file.
# It downloads into the directory it will be used from.
if [[ ! -f "$ref_zip" ]]; then
  wget -P "$reference_dir" http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
else
  echo "REFERENCE file already downloaded - skipping...."
fi

# Downloads hg19 data from SNPEFF
if [[ ! -d "$HOME/.snpEff/data/hg19" ]]; then
  snpEff download hg19
else
  echo "SNPEGG HG19 alredy downloaded, skipping..."
fi

# Annovar must be registered and downloaded to the same folder
# this script is run from, from here:
# http://www.openbioinformatics.org/annovar/annovar_download_form.php

#===============================================================

# PIPELINE

# The pipeline begins here. Large processes have been skipped
# if the files they produce already exist. This is to allow for
# faster re-running of the script after adjusting other parts.
# If you wish to adjust these parts you will have to delete the
# files already produced.

#===============================================================
# QUALITY ASSESSMENT - FASTQC
#===============================================================

# Comment on code
if [[ ! -f "ngs_pipeline/results/fastqc/NGS0001.R1_fastqc.html" ]]; then
  fastqc --threads $threads \
  "${untrimmed_dir}"*.fastq.gz \
  -o "$fastqc"
else
  echo "fastqc already performed, skipping fastqc"
fi

#==============================================================
# TRIMMOMATIC
#==============================================================

# Comment on code

if [[ ! -f "${trimmed_dir}${sample_id}_trimmed_R_1P.fastq" ]]; then
  trimmomatic PE \
  -threads "$threads" \
  -phred33 \
  "$rep1" "$rep2" \
  -baseout "${trimmed_dir}${sample_id}_trimmed_R.fastq" \
  ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 \
  TRAILING:25 MINLEN:50
else
  echo "TRIMMING already performed, skipping TRIMMING"
fi

#===============================================================
# REFERENCE DOWNLOAD AND UNZIP
#===============================================================
# Downloads reference if it doesn't already exist
if [[ ! -f "$ref_zip" ]]; then
  wget -O "$ref_zip" \
  ftp://ftp.ensembl.org/pub/grch37/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
else
  echo "Reference already exists, download skipped"
fi

# Unzips reference if the unzipped referenced doesn't exist.
if [[ ! -f "$ref_unzip" ]]; then
  zcat "$ref_zip" > "$ref_unzip"
else
  echo "Reference already unzipped, unzipping skipped"
fi

#================================================================
# REFERENCE INDEXING
#================================================================

# Indexes unzipped reference file with BWA, if one of the outputs doesn't exist.

if [[ ! -f "${ref_unzip}.amb" || \
      ! -f "${ref_unzip}.ann" || \
      ! -f "${ref_unzip}.bwt" || \
      ! -f "${ref_unzip}.pac" || \
      ! -f "${ref_unzip}.sa" ]]; then
  bwa index "$ref_unzip"
else
  echo "Reference already indexed, indexing skipped"
fi

#==========================================================
# SEQUENCE MAPPING
#==========================================================

# Runs BWA MEM if the SAM file doesn't exist # Around 60 mins
# Data is produced in the form of a Sequence Alignment Map, hence
# SAM file.

# Read group data can be taken from the first line of one of the reps,
# e.g. zcat NGS0001.R1.fastq.gz | head -n 1

if [[ ! -f $samfile ]]; then
  bwa mem \
  -t 2 \
  -v 1 \
  -R "@RG\tID:11V6WR1.111.D1375ACXX.1\tSM:NGS0001\tPL:ILLUMINA\tLB:NGS0001-library\tDT:2025-04-15\tPU:11V6WR1.111.D1375ACXX.1" \
  -I 250,50 \
  "$ref_unzip" \
  "$trimmed_R_1P" \
  "$trimmed_R_2P" \
  > "$samfile"
else
  echo "SAM file already exists, skipping...."
fi

# note: only use up to half of all of your machines threads (-t).

#==========================================================
# SAMFILE COMPRESSION TO BAMFILE
#==========================================================

# This step compresses the samfile into a bam file. Another form
# of compressed SAM files are CRAM. GATK tools (like Picard used further
# along in the pipeline require BAM and it's compression is lossless)

if [[ ! -f $bamfile ]]; then
  samtools view -@ "$threads" -h -b "$samfile" > "$bamfile"
else
  echo "BAMFILE already exists, SAM to BAM skipped"
fi

#=========================================================

# This code sorts the bamfile produced in the previous step.
# It organises the reads by the position on the genome to speed up
# further steps. It is required for indexing and other tools like variant callers.

if [[ ! -f "$sorted_bamfile" ]]; then
  samtools sort -@ "$threads" "$bamfile" > "$sorted_bamfile"
else
  echo "SORTED BAMFILE already exists - skipping...."
fi

#=========================================================

if [[ ! -f "${sorted_bamfile}.bai" ]]; then
  samtools index "$sorted_bamfile" # This will generate a .bai index file
else
  echo "INDEX already exists - skipping...."
fi

#=========================================================

# this uses Picard to mark duplicated reads
if [[ ! -f "$marked_bamfile" ]]; then
  picard MarkDuplicates I="$sorted_bamfile" O="$marked_bamfile" M="${mark_dup}marked_dup_metrics.txt"
else
  echo "BAM file already marked, skipping...."
fi
#=========================================================

# Index the marked bamfile
if [[ ! -f "${marked_bamfile}.bai" ]]; then
  samtools index "$marked_bamfile"
else
  echo "Marked BAM.BAI file already exists, skipping...."
fi
#=========================================================

# Filter the reads based on various things 
if [[ ! -f "$filtered_bamfile" ]]; then
  samtools view -F 1796  -q 20 -o "$filtered_bamfile" "$marked_bamfile"
else
  echo " Filtered BAM file already exists, skipping....."
fi

# Index filtered bamfile
if [[ ! -f "${filtered_bamfile}.bai" ]]; then
  samtools index "$filtered_bamfile"
else
  echo "Filtered BAM.BAI already exists, skipping..."
fi

# Perform flagstats
if [[ ! -f "${flagstats}bam_filtered_flagstat.txt" ]]; then
  samtools flagstat "$bamfile" > "${flagstats}bam_flagstat.txt"
  samtools flagstat "$sorted_bamfile" > "${flagstats}bam_sorted_flagstat.txt"
  samtools flagstat "$marked_bamfile" > "${flagstats}bam_marked_flagstat.txt"
  samtools flagstat "$filtered_bamfile" > "${flagstats}bam_filtered_flagstat.txt"
else
  echo "FLAGSTATS already run, skipping..."
fi

# Performs idxstats

if [[ ! -f "${idxstats}${sample_id}_sorted_filtered_idxstats.txt" ]]; then
  samtools idxstats "$filtered_bamfile" > "${idxstats}${sample_id}_sorted_filtered_idxstats.txt"
else
  echo "IDXSTATS file already exists, skipping...."
fi

# Performs insert size metrics

if [[ ! -f "${insert_size}metrics.txt" ]]; then
  picard CollectInsertSizeMetrics \
      I="$filtered_bamfile" \
      O="${insert_size}metrics.txt" \
      H="${insert_size}histogram.pdf" \
      M=0.5
else
  echo "Insert_size metrics exist, skipping...."
fi

# Cuts bedfile for more efficient processing
if [[ ! -f "$bedfile_stripped" ]]; then
  cut -f1,2,3 "$bedfile" > "$bedfile_stripped"
else
  echo "BEDFILE already stripped, skipping...."
fi

# Performs Depth of Coverage with mosdepth
if [[ ! -f $dep_cov_file ]]; then
  mosdepth --threads "$threads" --by "$bedfile_stripped" "$dep_cov_prefix" "$filtered_bamfile"
else
  echo "MOSDEPTH output already exists, skipping...."
fi

#=======================================================

# FREEBAYES - Make sure the code for the alternative
# variant caller, BCFTOOLS (below), is commented out.
# Comment this out nad uncomment BCFTOOLS to use it.

#======================================================

Uses FreeBayes to call varients using the bedfile to
restrict analysis of the region.
if [[ ! -f "$vcfFile_unzip" ]]; then
  freebayes \
    --bam "$filtered_bamfile" \
    --fasta-reference "$ref_unzip" \
    --targets "$bedfile" \
    --vcf "$vcfFile_unzip"
else
  echo
fi

bgzip "$vcfFile_unzip"
tabix -p vcf "$vcfFile_zip" # creates .vcf.gz.tbi

vcffilter -f "QUAL > 30 & QUAL / AO > 20 & SAF > 1 & SAR > 1 & RPR > 2 & RPL > 2" \
  "$vcfFile_zip" > "$vcfFile_filtered"

bgzip "$vcfFile_filtered"
tabix -p vcf "$vcfFile_filtered_zip"

#=================================================
# ALTERNATIVE TO FREEBAYES - BCFTOOLS
#=================================================

# if [[ ! -f $vcfFile_zip ]]; then
#   bcftools mpileup \
#     -f $ref_unzip \
#     -a AD,DP \
#     --threads $threads \
#     -Ou $filtered_bamfile |
#   bcftools call \
#     -mv \
#     --threads $threads \
#     -Oz \
#     -o $vcfFile_zip
# else
#   echo "BCFTOOLS output already exists, skipping...."
# fi

# tabix -f -p vcf $vcfFile_zip #index the compressed vcf

# bcftools filter -i "QUAL > 10 && FORMAT/DP > 20 && FORMAT/AD[0:1] > 8 && FORMAT/AD[0:1]/FORMAT/DP > 0.3" \
#  -o $vcfFile_filtered_zip $vcfFile_zip

# tabix -f -p vcf $vcfFile_filtered_zip

#===============================================
# BEDTOOLS INTERSECT
#===============================================

if [[ ! -f "$vcfFile_annotated_zip" ]]; then
  zcat "$vcfFile_filtered_zip" | bedtools intersect \
    -header \
    -wa \
    -a stdin \
    -b "$bedfile" \
    > "$vcfFile_annotated"

  bgzip "$vcfFile_annotated"
  tabix -p vcf "$vcfFile_annotated_zip"
else
  echo "BEDTOOLS intersect already run, skipping...."
fi

#===============================================================

# ANNOTATION USING ANNOVAR

#===============================================================

# Install Annovar
if [[ ! -d "./annovar" ]]; then
  tar -zxvf annovar.latest.tar.gz
else
  echo "Annovar already installed, skipping"
fi

annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/
annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro humandb/

annovar/convert2annovar.pl -format vcf4 "$vcfFile_annotated_zip" > "$avinput"

annovar/table_annovar.pl $avinput humandb/ \
  -buildver hg19 \
  -out "${results_dir}${sample_id}_filtered_annotation" \
  -remove \
  -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro \
  -operation g,g,f,f,f \
  -otherinfo \
  -nastring . \
  -csvout

###############################################################
#SNPEFF
###############################################################
  
if [[ ! -f "ngs_pipeline/results/NGS0001_snpeff.vcf" ]]; then
  snpEff -v hg19 "ngs_pipeline/results/NGS0001_filtered.vcf.gz" \
  > "ngs_pipeline/results/NGS0001_snpeff.vcf"
else
 echo "SNPEFF file already exists, skipping..."
fi

###############################################################
# VARIANT PRIORITISATION
################################################################

if [[ ! -f "ngs_pipeline/results/prioritized.vcf.gz" ]]; then
  bcftools view -i 'INFO/ANN ~ "exonic"' \
  ngs_pipeline/results/NGS0001_snpeff.vcf \
  -Oz -o ngs_pipeline/results/prioritized.vcf.gz
else
  echo "Variant prioritisation file already exists, skipping..."
fi
