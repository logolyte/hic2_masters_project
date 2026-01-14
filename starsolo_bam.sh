# Variables
base=$1
INDEX_DIR="/scratch/lema/m26_losu/mm10 gencode reference/star_index"  # STAR genome index directory
BAM_INPUT="/scratch/lema/m26_losu/bam_files/$base.bam"       # BAM file (mapped or unmapped)
out_dir="/scratch/lema/m26_losu/star_counts_gencode_mm10/$base/"              # output directory

STAR \
  --genomeDir "$INDEX_DIR" \
  --readFilesIn "$BAM_INPUT" \
  --readFilesType SAM SE \
  --readFilesCommand "samtools view -F 0x100" \
  --soloType CB_UMI_Simple \
  --soloUMIlen 12 \
  --soloInputSAMattrBarcodeSeq CR UR \
  --soloInputSAMattrBarcodeQual CY UY \
  --soloFeatures Gene GeneFull SJ Velocyto \
  --soloMultiMappers EM \
  --soloOutFileNames "./" \
  --outFileNamePrefix "$out_dir" \
  --outSAMtype None \
  --soloCellFilter EmptyDrops_CR \
  --soloCBwhitelist "/scratch/lema/m26_losu/10x-v3-whitelist-february-2018.txt" \
  --runThreadN 32 \
  --clipAdapterType CellRanger4 \
  --outFilterScoreMin 30 \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts\
  --soloUMIfiltering MultiGeneUMI_CR\
  --soloUMIdedup 1MM_CR
