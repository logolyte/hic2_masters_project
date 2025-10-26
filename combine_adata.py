#%% 
import scanpy
import anndata
import pandas

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
    12,
    12,
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

#%% Read adata files

adata_files = []

path = "RNA Sequencing Data/raw_reads/h5ad_files/"
filename = "counts_unfiltered/adata.h5ad"

for index, base_name in enumerate(base_names):
    full_path = f"{path}/{base_name}/{filename}"
    group = groups[index]
    day = days[index]
    print(f"Loading file: {base_name}")
    adata = scanpy.read_h5ad(full_path)
    adata.obs["sample"] = base_name
    adata.obs["group"] = groups[index]
    adata.obs["day"] = days[index]
    adata_files.append(adata)

adata_combined = adata_files[0].concatenate(
    adata_files[1:],
    join="outer",
    batch_key="sample",
    batch_categories=base_names
)

# Choose gene names as variable names
adata_combined.var["ensembl_id"] = adata_combined.var_names
adata_combined.var_names = adata_combined.var["gene_name"].astype(str)
if not adata_combined.var_names.is_unique:
    # Make duplicates unique by adding suffix
    adata_combined.var_names_make_unique()

# # Add var layer back
# # grab all var DataFrames from our dictionary
# all_var = [x.var for x in adata_files]
# # concatenate them
# all_var = pandas.concat(all_var, join="outer")
# # remove duplicates
# all_var = all_var[~all_var.duplicated()]
# adata_combined.var = all_var.loc[adata.var_names]

print("Finished reading files")

#%% Check output

print(adata_combined.obs)
print(adata_combined.X.shape)
print(adata_combined.layers)
print(adata_combined.var)

#print(adata_combined.X)
#%% Process data



# print(adata_combined.var)

# print(adata_combined.var_names)
print(adata_combined[:, "Krtdap"].la)

# %%
