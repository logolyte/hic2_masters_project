# Studying role of Hic2 in cellular reprogramming with computational models

Environments:
- `velocity`: Used for almost all notebooks and scripts.
- `unitvelo`: Used for the "Unitvelo.ipynb" notebook.

How to run:
1. Download [single-cell RNA sequencing data](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13029?key=fc47d3b3-2fed-4241-92ef-9b071a873d19) and put all BAM files in the directory called "RNA Sequencing Data/bam_files".
2. Put the file 220701 Pool1-Pool3.rds in the directory "RNA Sequencing Data". Convert it to a h5ad file with the script "convert_rds.r" and process it with the notebook "Analysis of Seurat data.ipynb".
3. Put narrowPeak files from ChIP-seq in the directory called "ChIP Sequencing Data".
4. Create a Starsolo index in the directory "RNA Sequencing Data/mm10 gencode reference". This can be done by putting an mm10 gtf and fa file with the primary assembly from [Gencode](https://www.gencodegenes.org/mouse/release_M10.html) in the directory and running "generate_star_index.sh". Also download the [10x v3 whitelist](https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/3M-february-2018.txt.gz), name it "10x-v3-whitelist-february-2018.txt" and put it in "RNA Sequencing Data".
5. Create a count matrix from the single-cell RNA sequencing reads with "starsolo_count_all.sh". Convert the resulting mtx files into an adata file with "Convert mtx files to Anndata.ipynb". This will create the Anndata object "RNA Sequencing Data/splice_counts.h5ad".
6. Run the velocity analysis with "SDEvelo.ipynb". The results will be saved to "RNA Sequencing Data/adata_sdevelo.h5ad". Results can be compared to those given with UniTVelo with "Unitvelo.ipynb" to those given by scVelo with "Velocity analysis scVelo.ipynb".
7. Run Cellrank analysis with "Cellrank.ipynb".
8. Find HIC2 and KLF4 targets with the notebook "Analyze ChIPSeq data.ipynb".
9. Perform differential expression analyses with "Differential expression.ipynb".
10. Perform gene set enrichment analyses and transcription factor estimations with ULM with "Compare macrostates.ipynb".
11. Identify driver genes with "Driver genes.ipynb".
12. Define and fit stochastic differential equation model of the reprogramming network with "Gene regulatory network.ipynb". For a motivation for why AUCell and ULM were chosen to score gene sets and transcription factor activity respectively, see the notebooks "Scoring gene sets.ipynb" and "Scoring transcription factor activities.ipynb".