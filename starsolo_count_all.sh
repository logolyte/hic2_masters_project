# Variables
INDEX_DIR=gencode_star_index  # STAR genome index directory

for bam in /scratch/lema/m26_losu/bam_files/*.bam; do
    base=$(basename "$bam" .bam)
    BAM_INPUT=$bam       # BAM file (mapped or unmapped)
    out_dir="/scratch/lema/m26_losu/star_counts_gencode_new/$base/"             # output directory
    echo "Processing $bam..."
    STAR \
      --genomeDir $INDEX_DIR \
      --readFilesIn $bam \
      --readFilesType SAM SE \
      --readFilesCommand "samtools view -F 0x100" \
      --soloType CB_UMI_Simple \
      --soloUMIlen 12 \
      --soloInputSAMattrBarcodeSeq CR UR \
      --soloInputSAMattrBarcodeQual CY UY \
      --soloFeatures Gene GeneFull SJ Velocyto \
      --soloOutFileNames "./" \
      --outFileNamePrefix $out_dir \
      --outSAMtype None \
      --soloCellFilter EmptyDrops_CR \
      --soloCBwhitelist "10x_index/v3-whitelist-february-2018.txt" \
      --runThreadN 16 \
      --clipAdapterType CellRanger4 \
      --outFilterScoreMin 30 \
      --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts\
      --soloUMIfiltering MultiGeneUMI_CR\
      --soloUMIdedup 1MM_CR
    echo "Finished $bam"
done
