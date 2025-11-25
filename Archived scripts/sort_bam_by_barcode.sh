for bam in bam_files/*.bam; do
    base=$(basename "$bam" .bam)
    echo "Processing $bam..."
    samtools sort -l 1 -m 1500M -t CB -O BAM -@ 16 \
        -o "sorted_bam/cellsorted_$base.bam" \
        $bam
    echo "Finished $bam"
done
