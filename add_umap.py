import scanpy
import anndata
from matplotlib import pyplot
import hdf5plugin
import numpy
import scvelo
import seaborn
import pandas

# Read Seurat object

seurat_path = "RNA Sequencing Data/seurat_object.h5ad"
adata_seurat = scanpy.read_h5ad(seurat_path)

seq_file_path = "RNA Sequencing Data/splice_counts.h5ad"
adata = scanpy.read_h5ad(seq_file_path)

# Add umap from Seurat object
sample_dataframes = []
for sample in adata.obs["sample"].unique():
    print(sample)
    sample_star = adata[adata.obs["sample"]==sample]
    sample_seurat = adata_seurat[adata_seurat.obs["sample"]==sample]
    barcodes_star = set(sample_star.obs["barcode"])
    barcodes_seurat = set(sample_seurat.obs["barcode"])
    barcodes_both = barcodes_star & barcodes_seurat
    print(f"Intersection: {len(barcodes_both)}")
    print(f"Only star: {len(barcodes_star - barcodes_seurat)}")
    print(f"Only Seurat: {len(barcodes_seurat - barcodes_star)}")
    print("Average counts per cell")
    print(f"Average any:                 {sample_star.layers["any"].sum() / len(sample_star)}")
    print(f"Spliced+unspliced+ambiguous: {(sample_star.layers["spliced"].sum() + sample_star.layers["unspliced"].sum() + sample_star.layers["ambiguous"].sum()) / len(sample_star)}")
    print(f"Cellranger total:            {sample_seurat.obs["nCount_RNA"].sum()/ len(sample_seurat)}")

    # Filter out barcodes not in Seurat object
    sample_combined = sample_star[sample_star.obs["barcode"].isin(barcodes_seurat)].copy()
    
    # Add umap from Seurat object
    umap_map = pandas.DataFrame(
        sample_seurat.obsm['X_umap'], 
        index=sample_seurat.obs['barcode']
    )
    sample_combined.obsm['X_umap'] = umap_map.loc[sample_combined.obs["barcode"]].values

    sample_dataframes.append(sample_combined)
    print()

# Concatenate
adata_combined = sample_dataframes[0].concatenate(
    sample_dataframes[1:],
    index_unique=None,
    batch_key=None
)

print(adata_combined)

# Save combined object
adata_combined.write_h5ad(
    "RNA Sequencing Data/splice_counts_umap.h5ad",
    compression=hdf5plugin.FILTERS["zstd"]
)