import scanpy
import anndata
import matplotlib
from matplotlib import pyplot
import hdf5plugin
import numpy
import scvelo
import seaborn
import pandas
import warnings
import cellrank
import gseapy
import decoupler
import IPython
import pygam
import textwrap

def plot_group_composition(adata, group_key, sample_key="sample", normalize=False):
    if f"{sample_key}_colors" in adata.uns.keys():
        sample_color_map = {}
        for sample, color in zip(adata.obs["sample"].cat.categories, adata.uns["sample_colors"]):
            sample_color_map[sample] = color
    else:
        sample_color_map = None

    # Calculate the frequency of each sample within each cluster
    composition = pandas.crosstab(adata.obs["macrostate"], adata.obs["sample"])
    if normalize:
        composition = composition.divide(composition.sum(axis=1), axis=0)

    # Plot the stacked bar chart
    composition.plot(kind="bar", stacked=True, figsize=(10, 6), color=sample_color_map)
    axis = pyplot.gca()
    axis.legend(title="Sample", bbox_to_anchor=(1.05, 1), loc="upper left")
    axis.grid(False)
    axis.set_xlabel(group_key.capitalize())
    if normalize:
        ylabel = "Cells"
    else:
        ylabel = "Proportion"
    axis.set_ylabel(ylabel)
    figure = pyplot.gcf()
    figure.tight_layout()
    pyplot.show()

# Function for plotting trends independent of trajectory

def plot_grouped_gene_trend(adata, genes, group_key="group", time_key="latent_time", layer="X", groups=None, columns=4, figure_width=14, plot_height=4, n_splines=6):
    """
    Plots smoothed gene expression using generalized additive models (GAMs).
    This provides true smooth local regression with confidence intervals.
    """

    number_of_genes = len(genes)
    number_of_columns = min(columns, number_of_genes)
    
    # Calculate rows
    rows = (number_of_genes + number_of_columns - 1) // number_of_columns
    
    figure, axes = pyplot.subplots(rows, number_of_columns)
    figure.set_figwidth(figure_width)
    figure.set_figheight(plot_height * rows)

    # Handle single gene case (axes is not an array) and 1D array case
    if number_of_genes == 1:
        axes_flat = [axes]
    else:
        axes_flat = axes.flatten()

    # Define colors for groups
    unique_groups = adata.obs[group_key].unique()
    if f"{group_key}_colors" in adata.uns_keys():
        palette = adata.uns[f"{group_key}_colors"]
        color_map = dict(zip(unique_groups, palette))
    else:
        palette = seaborn.color_palette("tab10", n_colors=len(unique_groups))
        color_map = dict(zip(unique_groups, palette))


    for gene_index, gene in enumerate(genes):
        axis = axes_flat[gene_index]
        
        # Fetch expression data
        if layer == "raw":
            expr_data = adata.raw.to_adata()[:, gene].X
        elif layer == "X":
            expr_data = adata[:, gene].X
        elif layer == "obs":
            expr_data = adata.obs[gene]
        else:
            expr_data = adata[:, gene].layers[layer]
        if hasattr(expr_data, "toarray"):
            expr_data = expr_data.toarray()
        if hasattr(expr_data, "flatten"):
            expr_data = expr_data.flatten()
        
        # Create DataFrame
        df = pandas.DataFrame({
            "Latent Time": adata.obs[time_key],
            "Expression": expr_data,
            "Group": adata.obs[group_key]
        })
        
        # Determine which groups to plot
        if groups is None:
            # If no specific list, plot all groups found in the data
            plot_groups = df["Group"].unique()
        else:
            plot_groups = groups

        # Loop through each group and fit a GAM
        for group in plot_groups:
            group_data = df[df["Group"] == group]
            color = color_map[group]
            
            # Skip empty groups
            if len(group_data) < 10:
                continue

            # Prepare X and y for pygam
            X = group_data["Latent Time"].values.reshape(-1, 1)
            y = group_data["Expression"].values
            
            # Fit LinearGAM using a spline term s(0)
            gam = pygam.LinearGAM(pygam.s(0, n_splines=n_splines)).fit(X, y)
            
            # Greate grid
            X_grid = numpy.linspace(X.min(), X.max(), 500)
            
            # Predict values and confidence intervals
            y_pred = gam.predict(X_grid)
            confidence_intervals = gam.confidence_intervals(X_grid, width=0.95)
            
            # Plot line
            axis.plot(X_grid, y_pred, label=group, color=color, linewidth=3)
            
            # Plot confidence interval
            axis.fill_between(
                X_grid, 
                confidence_intervals[:, 0], 
                confidence_intervals[:, 1], 
                color=color, 
                alpha=0.2
            )

        axis.set_title(f"{gene} expression")
        axis.set_xlabel("Latent Time")
        
        # Add legend to the first plot
        if gene_index == 0:
            axis.legend()

    # Clean up empty subplots
    for i in range(number_of_genes, len(axes_flat)):
        axes_flat[i].axis('off')

    pyplot.tight_layout()
    pyplot.show()

def compare_gsea(adata, gene_set, group_key="macrostate", layer="matrix", n_top_terms=5, significance_cutoff=0.05, sort_key="NES", remove_parenthesis=False,
                 wrap_width=30, font_size=10, min_dot_size=5, top_from="all", significance_metric="FDR q-val"):

    if sort_key == "FDR q-val":
        sort_ascending = True
        signed_key = False
    elif sort_key == "NES":
        sort_ascending = False
        signed_key =True

    # extract unique groups from the categorical column
    groups = adata.obs[group_key].cat.categories

    scanpy.tl.rank_genes_groups(
        adata,
        groupby=group_key,
        method="wilcoxon",
        key_added="macrostate_rank_genes",
        layer=layer,
        use_raw=False
    )

    # compute gene set enrichment for each individual group
    gsea_results = []
    for group in groups:
        rank_df = scanpy.get.rank_genes_groups_df(
            adata,
            key=f"{group_key}_rank_genes",
            group=group
        )
        rank_df.index = rank_df["names"]
        rank_df = rank_df.drop(columns="names")

        pre_res = gseapy.prerank(
            rnk=rank_df[["scores"]], 
            gene_sets=gene_set,
            threads=32,
            min_size=5,
            max_size=1000,
        )
        results_dataframe = pre_res.res2d.copy()
        results_dataframe["group"] = group
        gsea_results.append(results_dataframe)

    # concatenate all enrichment results into a single long-form dataframe
    gsea_dataframe = pandas.concat(gsea_results, ignore_index=True)

    # identify the top significant gene sets across all groups
    top_gene_sets = []
    if top_from == "all":
        for group in groups:
            significant_results = gsea_dataframe[(gsea_dataframe["group"] == group) & (gsea_dataframe[significance_metric] < significance_cutoff)]
            if signed_key:
                top_terms = significant_results.sort_values(by=sort_key, ascending=sort_ascending).head(n_top_terms)["Term"].tolist()
            else:
                top_terms = significant_results.query("NES > 0").sort_values(by=sort_key, ascending=sort_ascending).head(n_top_terms)["Term"].tolist()
            top_gene_sets.extend(top_terms)
    else:
        group = top_from
        significant_results = gsea_dataframe[(gsea_dataframe["group"] == group) & (gsea_dataframe[significance_metric] <= significance_cutoff)]
        if signed_key:
            top_terms = significant_results.sort_values(by=sort_key, ascending=sort_ascending).head(n_top_terms)["Term"].tolist()
        else:
            top_terms = significant_results.query("NES > 0").sort_values(by=sort_key, ascending=sort_ascending).head(n_top_terms)["Term"].tolist()
        top_gene_sets.extend(top_terms)

    # retain unique gene sets while preserving their initial order
    unique_top_gene_sets = list(dict.fromkeys(top_gene_sets))

    # filter the main dataframe to include only the selected top gene sets
    plot_dataframe = gsea_dataframe[gsea_dataframe["Term"].isin(unique_top_gene_sets)].copy()

    # sort dataframe
    plot_dataframe.index = plot_dataframe["Term"]
    plot_dataframe = plot_dataframe.loc[unique_top_gene_sets]

    # remove parentheses using vectorized string operations
    if remove_parenthesis:
        plot_dataframe["Term"] = plot_dataframe["Term"].str.split("(").str[0].str.strip()

    # wrap long text strings to prevent horizontal axes compression
    plot_dataframe["Term"] = plot_dataframe["Term"].apply(lambda text: textwrap.fill(text, width=wrap_width))

    # convert false discovery rate to negative log10 scale for dot size scaling
    plot_dataframe[significance_metric] = plot_dataframe[significance_metric].replace(0.0, 1e-10)
    
    # strictly cast the column to float to prevent numpy object array errors
    plot_dataframe[significance_metric] = plot_dataframe[significance_metric].astype(float)
    
    plot_dataframe["-log10(FDR)"] = -numpy.log10(plot_dataframe[significance_metric])
    
    # convert -0 to +0
    plot_dataframe["-log10(FDR)"] = numpy.abs(plot_dataframe["-log10(FDR)"])

    # determine the maximum absolute enrichment score to strictly center the color map at zero
    max_absolute_score = plot_dataframe["NES"].abs().max()

    # dynamically adjust figure size based on the number of terms and groups
    number_of_terms = len(plot_dataframe["Term"].unique())
    number_of_groups = len(groups)
    figure_height = max(6.0, number_of_terms * 0.6)
    figure_width = max(8.0, number_of_groups * 1.2)

    # generate the dotplot visualization
    matplotlib.pyplot.figure(figsize=(figure_width, figure_height))
    axis = seaborn.scatterplot(
        data=plot_dataframe,
        x="group",
        y="Term",
        hue="NES",
        size="-log10(FDR)",
        sizes=(min_dot_size, 300),
        palette="coolwarm",
        hue_norm=(-max_absolute_score, max_absolute_score)
    )

    # format the axes and reposition the legend outside the plot area
    axis.set_xlabel(group_key.capitalize(), fontsize=font_size)
    axis.set_xticklabels(axis.get_xticklabels(), rotation=45, rotation_mode="anchor", ha="right")
    axis.set_ylabel("Gene set", fontsize=font_size)
    axis.grid(False)
    axis.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.0)
    axis.tick_params(axis="both", labelsize=font_size)

    figure = axis.get_figure()
    figure.tight_layout()
    pyplot.show()

# Score TF activity in all macrostates

def plot_tf_activity(adata, network, group_key="macrostate", layer="matrix", n_top=5, top_from="all", significance_cutoff=0.05, min_dot_size=5, font_size=12,
    activity_cutoff=-numpy.inf, reverse = False, pathway_key="TF"):

    groups = adata.obs[group_key].cat.categories

    scanpy.tl.rank_genes_groups(
        adata,
        groupby=group_key,
        method="wilcoxon",
        key_added="macrostate_rank_genes",
        layer=layer,
        use_raw=False
    )

    tf_activity_dfs = []
    for group in groups:
        rank_df = scanpy.get.rank_genes_groups_df(
            adata,
            key=f"{group_key}_rank_genes",
            group=group
        )
        rank_df.index = rank_df["names"]
        rank_df = rank_df.drop(columns="names")

        data = rank_df[["scores"]].T.dropna(axis=1)
        tf_acts, tf_padj = decoupler.mt.ulm(data=data, net=network)
        activities_df = pandas.DataFrame(data={pathway_key: tf_acts.columns, "Activity": tf_acts.transpose()["scores"], "Adjusted p-value": tf_padj.transpose()["scores"],
                                               "group": group})
        tf_activity_dfs.append(activities_df)
    tf_activities = pandas.concat(tf_activity_dfs, ignore_index=True)

    top_gene_sets = []
    if top_from == "all":
        for group in groups:
            if reverse:
                significant_results = tf_activities[(tf_activities["group"] == group) & (tf_activities["Adjusted p-value"] < significance_cutoff) &\
                                                (tf_activities["Activity"] < activity_cutoff)]
            else:
                significant_results = tf_activities[(tf_activities["group"] == group) & (tf_activities["Adjusted p-value"] < significance_cutoff) &\
                                                (tf_activities["Activity"] > activity_cutoff)]
            top_terms = significant_results.sort_values(by="Activity", ascending=reverse).head(n_top)[pathway_key].tolist()
            top_gene_sets.extend(top_terms)
    else:
        group = top_from
        if reverse:
                significant_results = tf_activities[(tf_activities["group"] == group) & (tf_activities["Adjusted p-value"] < significance_cutoff) &\
                                            (tf_activities["Activity"] < activity_cutoff)]
        else:
            significant_results = tf_activities[(tf_activities["group"] == group) & (tf_activities["Adjusted p-value"] < significance_cutoff) &\
                                            (tf_activities["Activity"] > activity_cutoff)]
        top_terms = significant_results.sort_values(by="Activity", ascending=reverse).head(n_top)["TF"].tolist()
        top_gene_sets.extend(top_terms)

    # retain unique gene sets while preserving their initial order
    unique_top_gene_sets = list(dict.fromkeys(top_gene_sets))

    # filter the main dataframe to include only the selected top gene sets
    plot_dataframe = tf_activities[tf_activities[pathway_key].isin(unique_top_gene_sets)].copy()

    # sort dataframe
    plot_dataframe.index = plot_dataframe[pathway_key]
    plot_dataframe = plot_dataframe.loc[unique_top_gene_sets]

    # convert false discovery rate to negative log10 scale for dot size scaling
    plot_dataframe["Adjusted p-value"] = plot_dataframe["Adjusted p-value"].replace(0.0, 1e-10)
    
    # strictly cast the column to float to prevent numpy object array errors
    plot_dataframe["Adjusted p-value"] = plot_dataframe["Adjusted p-value"].astype(float)
    
    plot_dataframe["-log10(padj)"] = -numpy.log10(plot_dataframe["Adjusted p-value"])
    
    # convert -0 to +0
    plot_dataframe["-log10(padj)"] = numpy.abs(plot_dataframe["-log10(padj)"])

    # determine the maximum absolute enrichment score to strictly center the color map at zero
    max_absolute_score = plot_dataframe["Activity"].abs().max()

    # dynamically adjust figure size based on the number of terms and groups
    number_of_terms = len(plot_dataframe[pathway_key].unique())
    number_of_groups = len(groups)
    figure_height = max(6.0, number_of_terms * 0.6)
    figure_width = max(8.0, number_of_groups * 1.2)

    # generate the dotplot visualization
    matplotlib.pyplot.figure(figsize=(figure_width, figure_height))
    axis = seaborn.scatterplot(
        data=plot_dataframe,
        x="group",
        y=pathway_key,
        hue="Activity",
        size="-log10(padj)",
        sizes=(min_dot_size, 300),
        palette="coolwarm",
        hue_norm=(-max_absolute_score, max_absolute_score)
    )

    # format the axes and reposition the legend outside the plot area
    axis.set_xlabel(group_key.capitalize(), fontsize=font_size)
    axis.set_xticklabels(axis.get_xticklabels(), rotation=45, rotation_mode="anchor", ha="right")
    axis.set_ylabel(pathway_key, fontsize=font_size)
    axis.grid(False)
    axis.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.0)
    axis.tick_params(axis="both", labelsize=font_size)

    figure = axis.get_figure()
    figure.tight_layout()
    pyplot.show()