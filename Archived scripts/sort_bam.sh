for bam in bam_files/*.bam; do
    base=$(basename "$bam" .bam)
    echo "Processing $bam..."
    samtools sort -m 1500M -@ 16 \
        -o "sorted_bam/sorted_$base.bam" \
        $bam
    echo "Finished $bam"
done
