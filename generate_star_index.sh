STAR --runMode genomeGenerate \
    --runThreadN 32 \
     --genomeDir "RNA Sequencing Data/mm10 gencode reference/star_index/" \
     --genomeFastaFiles "RNA Sequencing Data/mm10 gencode reference/GRCm38.primary_assembly.genome.fa" \
     --sjdbGTFfile "RNA Sequencing Data/mm10 gencode reference/gencode.vM10.primary_assembly.annotation.gtf" \
     --genomeSAsparseD 3 \
