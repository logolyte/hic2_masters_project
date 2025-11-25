# # Variables
# INDEX_DIR=STAR_index  # STAR genome index directory
# BAM_INPUT=bam_files/Pool1_mESC.bam       # BAM file (mapped or unmapped)
# OUT_DIR=star_counts             # output directory
#
# STAR \
#   --genomeDir $INDEX_DIR \
#   --readFilesIn $BAM_INPUT \
#   --readFilesType SAM SE \
#   --readFilesCommand "samtools view -F 0x100" \
#   --soloType CB_UMI_Simple \
#   --soloCBwhitelist None \
#   --soloUMIlen 12 \
#   --soloInputSAMattrBarcodeSeq CR UR \
#   --soloInputSAMattrBarcodeQual CY UY \
#   --soloFeatures Gene GeneFull SJ Velocyto \
#   --soloOutFileNames Solo.out/ \
#   --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
#   --outSAMtype BAM SortedByCoordinate \
#   --soloCellFilter EmptyDrops_CR \
#   --runThreadN 16

#!/bin/bash

# -----------------------------
# Variables (CHANGE AS NEEDED)
# -----------------------------
INDEX_DIR="star_index"     # STAR genome index
WHITELIST="10x_index/v3-whitelist-february-2018.txt.gz"

FASTQ_DIR="fastq_files/Pool1_mESC"
OUT_DIR="starsolo_Pool1_mESC"

THREADS=16

# -----------------------------
# Collect FASTQ files
# -----------------------------
# All R1 from all lanes/subdirectories
R1_FILES=$(ls ${FASTQ_DIR}/**/*_R1_001.fastq.gz 2>/dev/null | tr '\n' ',' | sed 's/,$//')

# All R2 from all lanes/subdirectories
R2_FILES=$(ls ${FASTQ_DIR}/**/*_R2_001.fastq.gz 2>/dev/null | tr '\n' ',' | sed 's/,$//')

echo "Found R1 files:"
echo $R1_FILES
echo "Found R2 files:"
echo $R2_FILES

# -----------------------------
# Run STARsolo
# -----------------------------
STAR \
  --genomeDir $INDEX_DIR \
  --readFilesIn $R2_FILES $R1_FILES \
  --readFilesCommand zcat \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist $WHITELIST \
  --soloUMIlen 12 \
  --soloFeatures Gene GeneFull SJ Velocyto \
  --soloOutFileNames Solo.out/ \
  --soloCellFilter EmptyDrops_CR \
  --outFileNamePrefix ${OUT_DIR}/ \
  --runThreadN $THREADS \
  --clipAdapterType CellRanger4 \
  --outFilterScoreMin 30 \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts\
  --soloUMIfiltering MultiGeneUMI_CR\
  --soloUMIdedup 1MM_CR

#   --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
#   --outSAMtype BAM SortedByCoordinate \

