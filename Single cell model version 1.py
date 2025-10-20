#%% 

import scanpy
import anndata
from matplotlib import pyplot
import numpy

#%% Read input file

scanpy.settings.set_figure_params(dpi=100, facecolor="white")

seq_file_path = "RNA Sequencing Data/220701Pool1-Pool3.h5ad"
adata = scanpy.read_h5ad(seq_file_path)

#%% Analyze file

print(adata)
print(adata.obs.columns)
print("Identities: ", numpy.sort(adata.obs["orig.ident"].unique()))
print("Days: ", numpy.sort(adata.obs["Day"].unique()))
# Value counts
print(adata.obs["Day"].value_counts())

#%%

# Check how orig.ident relates to day
for identity in numpy.sort(adata.obs['orig.ident'].unique()):
    print(f"Identity: {identity}")
    day_data = adata.obs.loc[adata.obs['orig.ident']==identity, 'Day']
    print("Value counts:")
    print(day_data.value_counts())
    print()

#%% Map identities to day and group

print("Identities (Days):")
print(adata.obs["Day"].value_counts())

day_map={
    0: None,
    1: 2,
    2: 2,
    3: 4,
    4: 4,
    5: 6,
    6: 6,
    7: 9,
    8: 9,
    9: 12,
    10: 12,
    11: None,
}

group_map={
    0: "MEFs",
    1: "control",
    2: "Hic2",
    3: "control",
    4: "Hic2",
    5: "control",
    6: "Hic2",
    7: "control",
    8: "Hic2",
    9: "control",
    10: "Hic2",
    11: "ESCs",
}

adata.obs["exp_day"] = adata.obs["Day"].map(day_map)
adata.obs["group"] = adata.obs["Day"].map(group_map)

experiment_groups = {"control", "Hic2"}
reference_groups = {"MEFs", "ESCs"}

print("Experiment groups:")
print(f"{adata.obs[["exp_day", "group"]].value_counts(sort=False)}\n")

print("Reference groups:")
print(f"{adata.obs.query("group in @reference_groups")["group"].value_counts(sort=False)}\n")

#%% Quality control

# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^Hb[^(p)]")

scanpy.pp.calculate_qc_metrics(
    adata, qc_vars=["ribo", "hb"], inplace=True, log1p=True
)

scanpy.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts"],
    jitter=0.4,
    multi_panel=True,
)

scanpy.pl.scatter(adata, "total_counts", "n_genes_by_counts")

#%% Normalization if necessary

adata.var.loc["Krt2"]

#%% Clustering

scanpy.tl.pca(adata)
scanpy.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
scanpy.pl.pca(
    adata,
    color=["exp_day", "exp_day", "group", "group"],
    dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
    ncols=2,
    size=2,
)
