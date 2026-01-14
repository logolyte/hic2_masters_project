STAR --runMode genomeGenerate \
    --runThreadN 32 \
     --genomeDir "/scratch/lema/m26_losu/mm10 gencode reference/star_index/" \
     --genomeFastaFiles "/scratch/lema/m26_losu/mm10 gencode reference/GRCm38.primary_assembly.genome.fa" \
     --sjdbGTFfile "/scratch/lema/m26_losu/mm10 gencode reference/gencode.vM10.primary_assembly.annotation.gtf" \
     --genomeSAsparseD 3 \
