#%% 
import scanpy
import hdf5plugin
import pandas
import numpy
import anndata
import scipy

#%% Define groups

base_names = [
    "MEF",
    "Pool1_mESC",
    "Pool2_mESC",
    "Repro_Day2_BFP",
    "Repro_Day2_Hic2",
    "Repro_Day4_BFP",
    "Repro_Day4_Hic2",
    "Repro_Day6_BFP",
    "Repro_Day6_Hic2",
    "Repro_Day9_BFP",
    "Repro_Day9_Hic2",
    "Repro_Day12_BFP",
    "Repro_Day12_Hic2",
]

groups = [
    "MEF",
    "mESC",
    "mESC",
    "control",
    "Hic2",
    "control",
    "Hic2",
    "control",
    "Hic2",
    "control",
    "Hic2",
    "control",
    "Hic2",
]

days = [
    0,
    0,
    0,
    2,
    2,
    4,
    4,
    6,
    6,
    9,
    9,
    12,
    12,
]

#%% Function to convert mtx to anndata

def mtx_to_anndata(path):
    """Generate an anndata object from the STAR aligner output folder.
    Based on https://github.com/alexdobin/STAR/issues/774#issuecomment-850477636."""
    path=path
    # Load Read Counts
    X = scanpy.read_mtx(path+'Velocyto/raw/spliced.mtx')

    # Transpose counts matrix to have Cells as rows and Genes as cols as expected by AnnData objects
    X = X.X.transpose()

    # Load the 3 matrices containing Spliced, Unspliced and Ambigous reads
    mtxU = numpy.loadtxt(path+'Velocyto/raw/unspliced.mtx', skiprows=3, delimiter=' ')
    mtxS = numpy.loadtxt(path+'Velocyto/raw/spliced.mtx', skiprows=3, delimiter=' ')
    mtxA = numpy.loadtxt(path+'Velocyto/raw/ambiguous.mtx', skiprows=3, delimiter=' ')
    mtxG = numpy.loadtxt(path+'Gene/raw/matrix.mtx', skiprows=3, delimiter=' ')

    # Extract sparse matrix shape informations from the third row
    shapeU = numpy.loadtxt(path+'Velocyto/raw/unspliced.mtx', skiprows=2, max_rows = 1 ,delimiter=' ')[0:2].astype(int)
    shapeS = numpy.loadtxt(path+'Velocyto/raw/spliced.mtx', skiprows=2, max_rows = 1 ,delimiter=' ')[0:2].astype(int)
    shapeA = numpy.loadtxt(path+'Velocyto/raw/ambiguous.mtx', skiprows=2, max_rows = 1 ,delimiter=' ')[0:2].astype(int)
    shapeG = numpy.loadtxt(path+'Gene/raw/matrix.mtx', skiprows=2, max_rows = 1 ,delimiter=' ')[0:2].astype(int)

    # Read the sparse matrix with csr_matrix((data, (row_ind, col_ind)), shape=(M, N))
    # Subract -1 to rows and cols index because csr_matrix expects a 0 based index
    # Traspose counts matrix to have Cells as rows and Genes as cols as expected by AnnData objects

    spliced = scipy.sparse.csr_matrix((mtxS[:,2], (mtxS[:,0]-1, mtxS[:,1]-1)), shape = shapeS).transpose()
    unspliced = scipy.sparse.csr_matrix((mtxU[:,2], (mtxU[:,0]-1, mtxU[:,1]-1)), shape = shapeU).transpose()
    ambiguous = scipy.sparse.csr_matrix((mtxA[:,2], (mtxA[:,0]-1, mtxA[:,1]-1)), shape = shapeA).transpose()
    any = scipy.sparse.csr_matrix((mtxG[:,2], (mtxG[:,0]-1, mtxG[:,1]-1)), shape = shapeG).transpose()

    # Load Genes and Cells identifiers
    obs = pandas.read_csv(path+'Velocyto/raw/barcodes.tsv',
                  header = None, index_col = 0)

    # Remove index column name to make it compliant with the anndata format
    obs.index.name = None
    obs["barcode"] = obs.index

    var = pandas.read_csv(path+'Velocyto/raw/features.tsv', sep='\t', usecols=[0, 1], header=None,
                                    names = ["gene_ids", "gene_names"], index_col = 0)
  
    # Build AnnData object to be used with ScanPy and ScVelo
    adata = anndata.AnnData(X = X, obs = obs, var = var,
                                                 layers = {'spliced': spliced, 'unspliced': unspliced,
                                                           'ambiguous': ambiguous, "any": any})
    adata.var_names_make_unique()

    # Subset Cells based on STAR filtering
    selected_barcodes = pandas.read_csv(path+'Gene/filtered/barcodes.tsv', header = None)
    adata = adata[selected_barcodes[0]]

    return adata.copy()

#%% Read mtx files

adata_files = []

path = "RNA Sequencing Data/raw_reads/star_counts_gencode"

for index, base_name in enumerate(base_names):
    star_path = f"{path}/{base_name}/"
    group = groups[index]
    day = days[index]
    print(f"Loading file: {base_name}")
    
    adata = mtx_to_anndata(star_path)
    adata.obs["sample"] = base_name
    adata.obs["group"] = groups[index]
    adata.obs["day"] = days[index]
    adata_files.append(adata)


adata_combined = adata_files[0].concatenate(
    adata_files[1:],
    join="outer",
    batch_key="sample",
    batch_categories=base_names,
    index_unique="-"
)

# Choose gene names as variable names
adata_combined.var["ensemble_ids"] = adata_combined.var_names
adata_combined.var_names = adata_combined.var["gene_names"].astype(str)
if not adata_combined.var_names.is_unique:
    # Make duplicates unique by adding suffix
    adata_combined.var_names_make_unique()
# Rename gene name
adata_combined.var['gene_symbol'] = adata_combined.var['gene_names']
adata_combined.var = adata_combined.var.drop(columns='gene_names')

print("Finished reading files")

# %% Save adata file

adata_combined.write_h5ad(
    "RNA Sequencing Data/splice_counts.h5ad",
    compression=hdf5plugin.FILTERS["zstd"]
)

