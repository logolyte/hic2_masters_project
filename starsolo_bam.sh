# Variables
INDEX_DIR=gencode_star_index  # STAR genome index directory
BAM_INPUT=bam_files/Repro_Day2_BFP.bam       # BAM file (mapped or unmapped)
OUT_DIR=day2_bfp/             # output directory

STAR \
  --genomeDir $INDEX_DIR \
  --readFilesIn $BAM_INPUT \
  --readFilesType SAM SE \
  --readFilesCommand "samtools view -F 0x100" \
  --soloType CB_UMI_Simple \
  --soloUMIlen 12 \
  --soloInputSAMattrBarcodeSeq CR UR \
  --soloInputSAMattrBarcodeQual CY UY \
  --soloFeatures Gene GeneFull SJ Velocyto \
  --soloOutFileNames "./" \
  --outFileNamePrefix $OUT_DIR \
  --outSAMtype None \
  --soloCellFilter EmptyDrops_CR \
  --soloCBwhitelist "10x_index/v3-whitelist-february-2018.txt" \
  --runThreadN 16 \
  --clipAdapterType CellRanger4 \
  --outFilterScoreMin 30 \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts\
  --soloUMIfiltering MultiGeneUMI_CR\
  --soloUMIdedup 1MM_CR
