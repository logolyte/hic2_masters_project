"""
single_cell_analysis

A small library of utility functions for single-cell RNA-seq and trajectory analysis.
This module includes plotting helpers, gene set enrichment comparisons, transcription factor
activity scoring, and correlation utilities for AnnData objects.

Author: Love Sundin
Date: 2026-06-08
"""

import scanpy
import anndata
import matplotlib
from matplotlib import pyplot
import mpl_toolkits
import hdf5plugin
import numpy
import seaborn
import pandas
import warnings
import cellrank
import gseapy
import decoupler
import pygam
import textwrap
import scipy
import mpl_toolkits.axes_grid1

# Suppress non-critical warnings
warnings.simplefilter("ignore", category=UserWarning)
warnings.simplefilter("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=matplotlib.MatplotlibDeprecationWarning)

def plot_group_composition(adata, group_key, sample_key="sample", normalize=False, show=True, rotation=45):
    """
    Plot sample composition across groups.

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data matrix containing observations and sample metadata.
    group_key : str
        Column name in ``adata.obs`` containing the grouping variable to compare.
    sample_key : str, default "sample"
        Column name in ``adata.obs`` containing sample labels.
    normalize : bool, default False
        If True, show normalized proportions within each group.
    show : bool, default True
        If True, display the figure immediately. Otherwise return the axis object.
    rotation : int, default 45
        Rotation angle for x-axis tick labels.

    Returns
    -------
    matplotlib.axes.Axes or None
        Axis object when ``show`` is False, otherwise None.
    """
    # Retrieve pre-computed or auto-generate color palette for samples
    if f"{sample_key}_colors" in adata.uns.keys():
        sample_color_map = {}
        for sample, color in zip(adata.obs[sample_key].cat.categories, adata.uns[f"{sample_key}_colors"]):
            sample_color_map[sample] = color
    else:
        sample_color_map = None

    # Calculate sample frequency by group via contingency table
    composition = pandas.crosstab(adata.obs[group_key], adata.obs[sample_key])
    if normalize:
        composition = composition.divide(composition.sum(axis=1), axis=0)

    # Plot stacked bar chart
    composition.plot(kind="bar", stacked=True, figsize=(10, 6), color=sample_color_map)
    axis = pyplot.gca()
    axis.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    legend = axis.get_legend()
    legend.set_title(sample_key.capitalize(), {"size": matplotlib.rcParams.get("font.size")})
    axis.grid(False)
    axis.set_xlabel(group_key.capitalize())
    axis.set_xticklabels(axis.get_xticklabels(), rotation=rotation, rotation_mode="anchor", ha="right")
    if normalize:
        ylabel = "Proportion"
    else:
        ylabel = "Cells"
    axis.set_ylabel(ylabel)
    figure = pyplot.gcf()
    figure.tight_layout()
    if show:
        pyplot.show()
        return
    else:
        return axis

# Plot trends for grouped genes without trajectory-specific plotting
def plot_grouped_gene_trend(adata, genes, group_key="group", time_key="latent_time", layer=None, groups="all", columns=4, figure_width=14, plot_height=4, n_splines=10,
    obsm=False, show=True):
    """Plot smoothed gene expression trends by group.

    This function fits a generalized additive model (GAM) for each gene in ``genes`` and for each
    group defined by ``group_key``. Confidence intervals are shown around the fitted curves.

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data matrix containing expression values and cell annotations.
    genes : list[str]
        List of gene names to plot.
    group_key : str, default "group"
        Column name in ``adata.obs`` used to define the groups for trend comparison.
    time_key : str, default "latent_time"
        Column name in ``adata.obs`` containing pseudotime or latent time values.
    layer : str or None, default None
        AnnData layer key to use for expression values. Supported values include ``obs``, ``raw``, ``X``,
        or a custom layer name. If ``None`` and the gene is present in ``adata.obs``, that observation
        field will be used.
    groups : list[str] or "all", default "all"
        Specific groups to plot. If set to ``"all"``, all categories in ``adata.obs[group_key]`` are used.
    columns : int, default 4
        Number of columns in the subplot layout.
    figure_width : float, default 14
        Width of the overall figure in inches.
    plot_height : float, default 4
        Height of each row in the subplot grid in inches.
    n_splines : int, default 10
        Number of spline basis functions used by the GAM.
    obsm : bool, default False
        If True, read gene expression from ``adata.obsm[layer]`` instead of layers or ``X``.
    show : bool, default True
        If True, display the plot. When False, return the axes list.

    Returns
    -------
    list[matplotlib.axes.Axes] or None
        A list of axes objects when ``show`` is False; otherwise None.
    """

    number_of_genes = len(genes)
    number_of_columns = min(columns, number_of_genes)
    
    # Compute subplot row count
    rows = (number_of_genes + number_of_columns - 1) // number_of_columns
    
    figure, axes = pyplot.subplots(rows, number_of_columns)
    figure.set_figwidth(figure_width)
    figure.set_figheight(plot_height * rows)

    # Normalize axes to a flat list so indexing is consistent
    if number_of_genes == 1:
        axes_flat = [axes]
    else:
        axes_flat = axes.flatten()

    # Choose color mapping for groups (use stored palette when available)
    unique_groups = adata.obs[group_key].cat.categories
    # Use stored color mapping for consistency; fall back to default palette if unavailable
    if f"{group_key}_colors" in adata.uns_keys():
        palette = adata.uns[f"{group_key}_colors"]
        color_map = dict(zip(unique_groups, palette))
    else:
        palette = seaborn.color_palette("tab10", n_colors=len(unique_groups))
        color_map = dict(zip(unique_groups, palette))


    for gene_index, gene in enumerate(genes):
        axis = axes_flat[gene_index]
        
        # Fetch expression values from the requested source
        # Support multiple storage formats with fallback order: obsm -> obs -> raw -> X -> layers
        # Convert sparse or matrix objects to 1D numpy arrays for downstream modeling
        if obsm:
            expr_data = adata.obsm[layer][gene]
        else:
            if (layer == None and gene in adata.obs_keys()) or layer == "obs":
                expr_data = adata.obs[gene]
            elif layer == "raw":
                expr_data = adata.raw.to_adata()[:, gene].X
            elif layer == "X" or layer == None:
                expr_data = adata[:, gene].X
            else:
                expr_data = adata[:, gene].layers[layer]
        if hasattr(expr_data, "toarray"):
            expr_data = expr_data.toarray()
        if hasattr(expr_data, "flatten"):
            expr_data = expr_data.flatten()
        
        # Assemble a tidy DataFrame for modeling and plotting
        df = pandas.DataFrame({
            "Latent Time": adata.obs[time_key],
            "Expression": expr_data,
            "Group": adata.obs[group_key]
        })
        
        # Determine which groups to plot
        if groups is None:
            # Default to all groups when no selection is provided
            plot_groups = df["Group"].unique()
        else:
            plot_groups = groups

        # Fit a GAM for each selected group to smooth expression over latent time
        for group in plot_groups:
            group_data = df[df["Group"] == group]
            color = color_map[group]
            
            # Skip groups with insufficient cells for reliable GAM fitting
            if len(group_data) < 10:
                continue

            # Prepare predictor and response arrays for GAM fitting
            X = group_data["Latent Time"].values.reshape(-1, 1)
            y = group_data["Expression"].values
            
            # Fit a spline-based GAM with smoothing splines on latent time
            gam = pygam.LinearGAM(pygam.s(0, n_splines=n_splines)).fit(X, y)
            
            # Create a dense time grid for prediction and CI estimation
            X_grid = numpy.linspace(X.min(), X.max(), 500)
            
            # Compute predictions and 95% confidence intervals from the GAM
            y_pred = gam.predict(X_grid)
            confidence_intervals = gam.confidence_intervals(X_grid, width=0.95)
            
            # Plot the fitted trend line and shaded 95% CI
            axis.plot(X_grid, y_pred, label=group, color=color, linewidth=3)
            
            # Shade the confidence interval
            axis.fill_between(
                X_grid, 
                confidence_intervals[:, 0], 
                confidence_intervals[:, 1], 
                color=color, 
                alpha=0.2
            )

        axis.set_title(f"{gene}")
        axis.set_xlabel("Latent Time")
        
        # Add legend to the first plot
        if gene_index == 0:
            axis.legend()

    # Turn off axes for any unused subplot cells
    for i in range(number_of_genes, len(axes_flat)):
        axes_flat[i].axis('off')

    pyplot.tight_layout()
    if show:
        pyplot.show()
    else:
        return axes_flat

def compare_gsea(data, gene_set, group_key="macrostate", layer=None, n_top_terms=5, significance_cutoff=0.05, sort_key="NES", remove_parenthesis=False,
                 wrap_width=30, font_size=14, min_dot_size=5, top_from="all", significance_metric="FDR q-val", enrichment_cutoff=None, reverse=False,
                 terms="all", permutation_num=1000, show=False, xlabel="", ylabel="", row_height=0.4, column_width=0.4):
    """Compare gene set enrichment analysis results across groups.

    This utility performs preranked GSEA using either an AnnData object or a precomputed
    dictionary of ranked gene lists. Results are plotted as a bubble chart with NES values
    and FDR-based point sizes.

    Parameters
    ----------
    data : anndata.AnnData or dict
        Input data. When AnnData is provided, differential ranking is computed from ``adata.obs[group_key]``.
        When a dictionary is provided, each value should be a DataFrame with ranking results.
    gene_set : str or dict
        Gene set database or dictionary accepted by :func:`gseapy.prerank`.
    group_key : str, default "macrostate"
        Observation column used for group ranking within AnnData.
    layer : str or None, default None
        Optional layer to use for rank gene computation when ``data`` is AnnData.
    n_top_terms : int, default 5
        Number of top gene sets to display per group.
    significance_cutoff : float, default 0.05
        Maximum adjusted p-value to consider a gene set significant.
    sort_key : str, default "NES"
        Column used to rank terms within each group.
    remove_parenthesis : bool, default False
        If True, strip trailing parenthetical text from gene set names.
    wrap_width : int, default 30
        Maximum text width for wrapped gene set names.
    font_size : int, default 14
        Font size used for axis labels and ticks.
    min_dot_size : int, default 5
        Minimum point size for the scatterplot.
    top_from : str, default "all"
        If ``"all"``, top terms are selected across all groups; otherwise only from the specified group.
    significance_metric : str, default "FDR q-val"
        Column name used to define significance for plotting.
    enrichment_cutoff : float or None, default None
        Optional cutoff on NES values for term selection.
    reverse : bool, default False
        If True, select terms in ascending rather than descending order.
    terms : list[str] or "all", default "all"
        Specific gene set names to include in the final plot.
    permutation_num : int, default 1000
        Number of permutations used for prerank GSEA.
    show : bool, default False
        If True, display the figure immediately. Otherwise return the axis object.
    xlabel : str, default ""
        Label for the x-axis.
    ylabel : str, default ""
        Label for the y-axis.
    row_height : float, default 0.4
        Height of each term row in the figure grid.
    column_width : float, default 0.4
        Width of each group column in the figure grid.

    Returns
    -------
    matplotlib.axes.Axes or None
        The axis object when ``show`` is False, otherwise None.
    """

    if type(data) == anndata.AnnData:
        adata = data
        # Extract group categories from the annotation column
        groups = adata.obs[group_key].cat.categories

        scanpy.tl.rank_genes_groups(
            adata,
            groupby=group_key,
            method="wilcoxon",
            key_added=f"{group_key}_rank_genes",
            layer=layer,
            use_raw=False
        )
    elif type(data) == dict:
        gsea_results = []
        groups = data.keys()

    # Compute gene set enrichment for each group
    gsea_results = []
    for group in groups:
        if type(data) == anndata.AnnData:
            rank_df = scanpy.get.rank_genes_groups_df(
                adata,
                key=f"{group_key}_rank_genes",
                group=group
            )
        elif type(data) == dict:
            rank_df = data[group]
        if "names" in rank_df.columns:
            rank_df.index = rank_df["names"]
            rank_df = rank_df.drop(columns="names")

        pre_res = gseapy.prerank(
            rnk=rank_df[["scores"]], 
            gene_sets=gene_set,
            threads=32,
            min_size=5,
            max_size=1000,
            permutation_num=permutation_num
        )
        results_dataframe = pre_res.res2d.copy()
        results_dataframe["group"] = group
        gsea_results.append(results_dataframe)

    # Concatenate all enrichment results into a single long-form DataFrame
    gsea_dataframe = pandas.concat(gsea_results, ignore_index=True)

    # Filter to user-specified term names if provided
    if not terms == "all":
        gsea_dataframe = gsea_dataframe[gsea_dataframe["Term"].isin(terms)]

    # Select top significant gene sets based on filtering criteria
    # Filters: group membership, significance threshold, and optional enrichment direction
    top_gene_sets = []
    if top_from == "all":
        group_set = groups
    else:
        group_set = [top_from]
    for group in group_set:
            # Apply significance and enrichment cutoff filters
            gene_set_filter = (gsea_dataframe["group"] == group) & (gsea_dataframe[significance_metric] <= significance_cutoff)
            if not enrichment_cutoff is None:
                if reverse:
                    gene_set_filter &= (gsea_dataframe["NES"] < enrichment_cutoff)
                else:
                    gene_set_filter &= (gsea_dataframe["NES"] > enrichment_cutoff)
            significant_results = gsea_dataframe[gene_set_filter]
            top_terms = significant_results.sort_values(by=sort_key, ascending = reverse).head(n_top_terms)["Term"].tolist()
            top_gene_sets.extend(top_terms)

    # Retain unique gene sets while preserving their initial order
    unique_top_gene_sets = list(dict.fromkeys(top_gene_sets))[::-1]

    # Filter the DataFrame to include only the selected top gene sets
    plot_dataframe = gsea_dataframe[gsea_dataframe["Term"].isin(unique_top_gene_sets)].copy()

    # Reorder the results by the selected gene set order
    plot_dataframe.index = plot_dataframe["Term"]
    plot_dataframe = plot_dataframe.loc[unique_top_gene_sets]

    # Optionally remove parenthetical text from term labels
    if remove_parenthesis:
        plot_dataframe["Term"] = plot_dataframe["Term"].str.split("(").str[0].str.strip()

    # Wrap long term labels to reduce horizontal label overlap
    plot_dataframe["Term"] = plot_dataframe["Term"].apply(lambda text: textwrap.fill(text, width=wrap_width))

    # Convert the significance metric to a numeric scale for plotting
    plot_dataframe[significance_metric] = plot_dataframe[significance_metric].replace(0.0, 1/permutation_num)
    
    # Cast significance values to float to avoid dtype issues
    plot_dataframe[significance_metric] = plot_dataframe[significance_metric].astype(float)
    
    fdr_key = "-log$_{10}$(FDR)"
    plot_dataframe[fdr_key] = -numpy.log10(plot_dataframe[significance_metric])
    
    # Convert negative zero values to positive zero
    plot_dataframe[fdr_key] = numpy.abs(plot_dataframe[fdr_key])

    # Determine the maximum absolute enrichment score for color scaling
    max_absolute_score = plot_dataframe["NES"].abs().max()

    # Dynamically adjust the figure size for the number of terms and groups
    number_of_terms = len(plot_dataframe["Term"].unique())
    number_of_groups = len(groups)
    
    desired_axes_width = max(2.0, number_of_groups * column_width)
    desired_axes_height = max(2.0, number_of_terms * row_height)

    # Allocate safe margins around the plot
    left_margin = 8.0
    right_margin = 4.0
    bottom_margin = 3.0
    top_margin = 1.0

    figure_width = left_margin + desired_axes_width + right_margin
    figure_height = bottom_margin + desired_axes_height + top_margin

    # Create figure with a fixed grid layout sized to terms and groups
    figure = matplotlib.pyplot.figure(figsize=(figure_width, figure_height))

    # Define fixed size grid regions
    horizontal_regions = [
        mpl_toolkits.axes_grid1.Size.Fixed(left_margin),
        mpl_toolkits.axes_grid1.Size.Fixed(desired_axes_width),
        mpl_toolkits.axes_grid1.Size.Fixed(right_margin)
    ]

    vertical_regions = [
        mpl_toolkits.axes_grid1.Size.Fixed(bottom_margin),
        mpl_toolkits.axes_grid1.Size.Fixed(desired_axes_height),
        mpl_toolkits.axes_grid1.Size.Fixed(top_margin)
    ]

    divider = mpl_toolkits.axes_grid1.Divider(
        figure,
        pos=(0.0, 0.0, 1.0, 1.0),
        horizontal=horizontal_regions,
        vertical=vertical_regions
    )

    # Place the main plotting axis into the central layout cell
    axis = figure.add_axes(
        rect=(0.0, 0.0, 1.0, 1.0),
        axes_locator=divider.new_locator(nx=1, ny=1)
    )

    seaborn.scatterplot(
        data=plot_dataframe,
        x="group",
        y="Term",
        hue="NES",
        size=fdr_key,
        sizes=(min_dot_size, 300),
        palette="coolwarm",
        hue_norm=(-max_absolute_score, max_absolute_score),
        ax=axis
    )

    # Explicitly set categorical axis limits to prevent marker clipping when plotting
    axis.set_xlim(-0.5, number_of_groups - 0.5)
    axis.set_ylim(-0.5, number_of_terms - 0.5)

    # Format axes and move legend outside the plotting area for readability
    axis.set_xlabel(xlabel, fontsize=font_size)
    axis.set_xticklabels(axis.get_xticklabels(), rotation=45, rotation_mode="anchor", ha="right")
    axis.set_ylabel(ylabel, fontsize=font_size)
    axis.grid(False)
    axis.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.0)
    axis.tick_params(axis="both", labelsize=font_size)

    if show:
        matplotlib.pyplot.show()
    else:
        return axis

def plot_tf_activity(adata, network, group_key="macrostate", layer=None, n_top=5, top_from="all", significance_cutoff=0.05, min_dot_size=5, font_size=14,
    activity_cutoff=-numpy.inf, reverse = False, pathway_key="", show=False, xlabel=""):
    """Visualize transcription factor activity across groups.

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data matrix with expression and metadata.
    network : array-like or pandas.DataFrame
        Regulatory network used by :func:`decoupler.mt.ulm`.
    group_key : str, default "macrostate"
        Observation column used to define comparison groups.
    layer : str or None, default None
        Optional layer used for ranking genes in differential expression.
    n_top : int, default 5
        Number of top activities to display per group.
    top_from : str, default "all"
        Group used to select top activities. If ``"all"``, selections come from every group.
    significance_cutoff : float, default 0.05
        Adjusted p-value threshold for including activities.
    min_dot_size : int, default 5
        Minimum dot size for the scatterplot.
    font_size : int, default 14
        Font size for axis labels and tick labels.
    activity_cutoff : float, default -numpy.inf
        Activity threshold used to filter selected terms.
    reverse : bool, default False
        If True, select terms with lower activity scores preferentially.
    pathway_key : str, default ""
        Column name for the pathway or transcription factor labels.
    show : bool, default False
        If True, display the figure immediately. Otherwise return the axis object.
    xlabel : str, default ""
        Label for the x-axis.

    Returns
    -------
    matplotlib.axes.Axes or None
        Axis object when ``show`` is False; otherwise None.
    """

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
        group_set = groups
    else:
        group_set = [top_from]
    for group in group_set:
            gene_set_filter = (tf_activities["group"] == group) & (tf_activities["Adjusted p-value"] <= significance_cutoff)
            if not activity_cutoff is None:
                if reverse:
                    gene_set_filter &= (tf_activities["Activity"] < activity_cutoff)
                else:
                    gene_set_filter &= (tf_activities["Activity"] > activity_cutoff)
            significant_results = tf_activities[gene_set_filter]
            top_terms = significant_results.sort_values(by="Activity", ascending = reverse).head(n_top)[pathway_key].tolist()
            top_gene_sets.extend(top_terms)

    # Retain unique gene sets while preserving their initial order
    unique_top_gene_sets = list(dict.fromkeys(top_gene_sets))[::-1]

    # Filter the DataFrame to include only the selected top gene sets
    plot_dataframe = tf_activities[tf_activities[pathway_key].isin(unique_top_gene_sets)].copy()

    # Reorder the results by the selected gene set order
    plot_dataframe.index = plot_dataframe[pathway_key]
    plot_dataframe = plot_dataframe.loc[unique_top_gene_sets]

    # Convert adjusted p-values to a numeric scale for plotting
    plot_dataframe["Adjusted p-value"] = plot_dataframe["Adjusted p-value"].replace(0.0, 1e-10)
    
    # Cast adjusted p-values to float to avoid dtype issues
    plot_dataframe["Adjusted p-value"] = plot_dataframe["Adjusted p-value"].astype(float)

    padj_key = "-log$_{10}(p_{adj})$"
    plot_dataframe[padj_key] = -numpy.log10(plot_dataframe["Adjusted p-value"])
    
    # Convert negative zero values to positive zero
    plot_dataframe[padj_key] = numpy.abs(plot_dataframe[padj_key])

    # Determine the maximum absolute activity score for color scaling
    max_absolute_score = plot_dataframe["Activity"].abs().max()

    # Dynamically adjust figure size for the number of terms and groups
    number_of_terms = len(plot_dataframe[pathway_key].unique())
    number_of_groups = len(groups)
    figure_height = max(6.0, number_of_terms * 0.6)
    figure_width = max(5.0, number_of_groups * 0.7)

    # Create dotplot figure and axis
    matplotlib.pyplot.figure(figsize=(figure_width, figure_height))
    axis = seaborn.scatterplot(
        data=plot_dataframe,
        x="group",
        y=pathway_key,
        hue="Activity",
        size=padj_key,
        sizes=(min_dot_size, 300),
        palette="coolwarm",
        hue_norm=(-max_absolute_score, max_absolute_score)
    )

    # Format axes and move legend outside the plotting area for readability
    axis.set_xlabel(xlabel, fontsize=font_size)
    axis.set_xticklabels(axis.get_xticklabels(), rotation=45, rotation_mode="anchor", ha="right")
    axis.set_ylabel(pathway_key, fontsize=font_size)
    axis.grid(False)
    axis.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.0)
    axis.tick_params(axis="both", labelsize=font_size)

    # Explicitly set categorical axis limits to prevent marker clipping when plotting
    axis.set_xlim(-0.5, number_of_groups - 0.5)
    axis.set_ylim(-0.5, number_of_terms - 0.5)

    figure = axis.get_figure()
    figure.tight_layout()

    if show:
        pyplot.show()
    else:
        return axis

def calculate_correlation(adata, variable_1, variable_2, group_key="group", groups="all", layers=["Ms", "Ms"], print_results=True, method="spearman"):
    """Calculate correlation between two variables across groups.

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data matrix containing variables and layers.
    variable_1 : str
        Name of the first variable to correlate.
    variable_2 : str
        Name of the second variable to correlate.
    group_key : str, default "group"
        Observation column used to split the data into groups.
    groups : list[str] or "all", default "all"
        Specific groups to evaluate. If ``"all"``, all categories in ``adata.obs[group_key]`` are used.
    layers : list[str], default ["Ms", "Ms"]
        Layer names for each variable when they are not found in ``adata.obs``.
    print_results : bool, default True
        If True, print summary statistics for each group and total correlation.
    method : str, default "spearman"
        Correlation method, either ``"spearman"`` or ``"pearson"``.

    Returns
    -------
    dict or None
        Dictionary of correlation results when ``print_results`` is False; otherwise None.
    """
    if print_results:
        print(f"Correlation between {variable_1} and {variable_2}")
    # Select correlation function based on method parameter
    if method == "spearman":
        correlation_method = scipy.stats.spearmanr
    elif method == "pearson":
        correlation_method = scipy.stats.pearsonr
    results_dict = {}
    # Expand "all" to the full list of group categories if specified
    if groups == "all":
        groups = adata.obs[group_key].cat.categories
    # Compute correlation for each group
    for group in groups:
        indices = adata.obs[group_key] == group
        arrays = []
        # Extract expression values for both variables from obs or specified layers
        for variable, layer in zip([variable_1, variable_2], layers):
            if variable in adata.obs.columns:
                # Fetch from observation metadata (e.g., pre-computed scores)
                array = adata[indices].obs[variable].to_numpy()
            else:
                # Fetch from expression matrix, handling sparse and dense formats
                if layer == "raw":
                    array = adata.raw.to_adata()[indices, variable].X
                    if hasattr(array, "toarray"):
                        array = array.toarray()
                    else:
                        array = numpy.array(array)
                elif layer == "X":
                    array = adata[indices, variable].X
                    if hasattr(array, "toarray"):
                        array = array.toarray()
                    else:
                        array = numpy.array(array)
                else:
                    array = numpy.array(adata[indices, variable].layers[layer])
            array = array.flatten()
            arrays.append(array)
        # Compute correlation and handle scalar vs array return values
        correlation_result = correlation_method(arrays[0], arrays[1])
        if isinstance(correlation_result.statistic, numpy.ndarray):
            statistic = correlation_result.statistic[0]
        else:
            statistic = correlation_result.statistic
        if isinstance(correlation_result.pvalue, numpy.ndarray):
            pvalue = correlation_result.pvalue[0]
        else:
            pvalue = correlation_result.pvalue
        if print_results:
            print(f"{group}: {statistic:.3g}, p-value: {pvalue:.3g}")
        else:
            results_dict[group] = {"statistic": statistic, "pvalue": pvalue}
    
    # Compute overall correlation across all selected groups
    if len(groups) > 1:
        indices = adata.obs[group_key].isin(groups)
        arrays = []
        for variable, layer in zip([variable_1, variable_2], layers):
            if variable in adata.obs.columns:
                array = adata[indices].obs[variable].to_numpy()
            else:
                if layer == "raw":
                    array = adata.raw.to_adata()[indices, variable].X
                    if hasattr(array, "toarray"):
                        array = array.toarray()
                    else:
                        array = numpy.array(array)
                elif layer == "X":
                    array = adata[indices, variable].X
                    if hasattr(array, "toarray"):
                        array = array.toarray()
                    else:
                        array = numpy.array(array)
                else:
                    array = numpy.array(adata[indices, variable].layers[layer])
            array = array.flatten()
            arrays.append(array)
        correlation_result = correlation_method(arrays[0], arrays[1])
        if isinstance(correlation_result.statistic, numpy.ndarray):
            statistic = correlation_result.statistic[0]
        else:
            statistic = correlation_result.statistic
        if isinstance(correlation_result.pvalue, numpy.ndarray):
            pvalue = correlation_result.pvalue[0]
        else:
            pvalue = correlation_result.pvalue
        if print_results:
            print(f"Total: {statistic:.3g}, p-value: {pvalue:.3g}")
        else:
            results_dict["total"] = {"statistic": statistic, "pvalue": pvalue}
    if not print_results:
        return results_dict

def plot_gene_correlation(adata, x_genes, y_genes, layer="Ms", annotation_style="stars", figsize=(10, 8), title="", group_key=None, groups=None):
    """Plot pairwise gene correlation as a heatmap.

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data matrix containing gene expression data.
    x_genes : list[str]
        List of genes plotted on the x-axis.
    y_genes : list[str]
        List of genes plotted on the y-axis.
    layer : str, default "Ms"
        AnnData layer containing expression values to use for correlations.
    annotation_style : str, default "stars"
        Annotation style for significance labels: ``"stars"`` or ``"p-value"``.
    figsize : tuple, default (10, 8)
        Size of the output figure in inches.
    title : str, default ""
        Plot title.
    group_key : str or None, default None
        Observation column used to subset the data before plotting.
    groups : list[str] or str or None, default None
        Specific group(s) to include when subsetting via ``group_key``.

    Returns
    -------
    None
        The function displays the heatmap directly.
    """
    
    # Subset AnnData when specific groups are requested
    if group_key is not None and groups is not None:
        if isinstance(groups, str):
            groups = [groups]
        adata_subset = adata[adata.obs[group_key].isin(groups)]
    else:
        adata_subset = adata
    
    # Filter genes to those present in the subsetted dataset
    valid_x_genes = [gene for gene in x_genes if gene in adata_subset.var_names]
    valid_y_genes = [gene for gene in y_genes if gene in adata_subset.var_names]
    
    # Extract expression matrices for the selected gene sets
    x_data = pandas.DataFrame(
        data=adata_subset[:, valid_x_genes].layers[layer], 
        index=adata_subset.obs.index, 
        columns=valid_x_genes
    )
    y_data = pandas.DataFrame(
        data=adata_subset[:, valid_y_genes].layers[layer], 
        index=adata_subset.obs.index, 
        columns=valid_y_genes
    )
    
    # Initialize result matrices for correlation and significance
    correlation_matrix = pandas.DataFrame(index=valid_y_genes, columns=valid_x_genes, dtype=float)
    pvalue_matrix = pandas.DataFrame(index=valid_y_genes, columns=valid_x_genes, dtype=float)
    
    # Calculate Spearman correlation for each gene pair
    for y_gene in valid_y_genes:
        for x_gene in valid_x_genes:
            correlation, pvalue = scipy.stats.spearmanr(y_data[y_gene], x_data[x_gene])
            correlation_matrix.loc[y_gene, x_gene] = correlation
            pvalue_matrix.loc[y_gene, x_gene] = pvalue
            
    # Build annotation matrix for significance symbols or p-values based on user preference
    annotation_matrix = pandas.DataFrame(index=valid_y_genes, columns=valid_x_genes, dtype=str)
    
    for y_gene in valid_y_genes:
        for x_gene in valid_x_genes:
            p = pvalue_matrix.loc[y_gene, x_gene]
            # Assign symbol (*, **, ***) or p-value text based on significance threshold
            if annotation_style == "stars":
                if p < 0.001:
                    annotation = "***"
                elif p < 0.01:
                    annotation = "**"
                elif p < 0.05:
                    annotation = "*"
                else:
                    annotation = ""
            elif annotation_style == "p-value":
                annotation = f"{p:.1e}"
            else:
                annotation = ""
                
            annotation_matrix.loc[y_gene, x_gene] = annotation
            
    # Create figure and main axis for the heatmap
    figure, axis = matplotlib.pyplot.subplots(figsize=figsize)

    # Add a colorbar axis whose height matches the heatmap height
    divider = mpl_toolkits.axes_grid1.make_axes_locatable(axis)
    colorbar_axis = divider.append_axes("right", size="5%", pad=0.1)

    # Plot the correlation heatmap with significance annotations and colorbar
    seaborn.heatmap(
        data=correlation_matrix,
        annot=annotation_matrix,
        fmt="",
        cmap="coolwarm",
        vmin=-1,
        vmax=1,
        center=0,
        square=True,
        xticklabels=True,
        yticklabels=True,
        ax=axis,
        cbar_ax=colorbar_axis,
        cbar_kws={"label": "Spearman Correlation"}
    )

    # Set axis labels
    axis.set_xlabel("Target Genes")
    axis.set_ylabel("Source Genes")

    # Rotate x-axis labels and normalize y-axis orientation for readability
    axis.set_xticklabels(axis.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    axis.set_yticklabels(axis.get_yticklabels(), rotation=0)

    axis.set_title(title)

    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.show()


def plot_trends_by_trajectory(
    adata, 
    tfs, 
    time_key="latent_time", 
    model=None, 
    group_key=None, 
    groups="all",
    lineages="all",
    layer=None, 
    n_knots=10, 
    smoothing_penalty=10.0,
    n_columns=4
):
    """Visualize transcription factor trends along cell fate trajectories.

    This function constructs an AnnData object containing the selected transcription factor
    expression values and applies a CellRank GAM model to plot gene trends for one or more
    lineages.

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data matrix with expression values and lineage metadata.
    tfs : list[str]
        List of transcription factors or variables to plot.
    time_key : str, default "latent_time"
        Observation column containing pseudotime or latent time values.
    model : cellrank.models.GAMR or None, default None
        Precomputed CellRank model. If None, a new model is fit on the selected data.
    group_key : str or None, default None
        Observation column used to filter cells by group.
    groups : list[str] or "all", default "all"
        Specific groups to include when filtering. Applied only if ``group_key`` is provided.
    lineages : list[str] or "all", default "all"
        Lineages to plot. If ``"all"``, all categories in ``adata.obs["term_states_fwd"]`` are used.
    layer : str or None, default None
        Data source for transcription factor values. Supported values include ``None``, ``"obsm"``,
        ``"obs"``, ``"X"``, ``"raw"``, or a custom layer name.
    n_knots : int, default 10
        Number of knots for the GAM basis.
    smoothing_penalty : float, default 10.0
        Penalty applied to the GAM smoothing term.
    n_columns : int, default 4
        Number of columns in the gene trends plot layout.

    Returns
    -------
    None
        The function displays the trend plot directly.
    """

    if group_key is None or groups == "all":
        adata_sub = adata
    else:
        group_mask = adata.obs[group_key].isin(groups)
        adata_sub = adata[group_mask].copy()

    term_state_colors = adata.uns["term_states_fwd_colors"].copy()

    # Extract transcription factor or activity values from the specified layer
    if layer == None:
        # Default: fetch from obs (metadata) or X (expression matrix)
        data = {}
        index = adata_sub.obs_names
        for var in tfs:
            if var in adata.obs_keys():
                # Use pre-computed metadata (e.g., scores)
                data[var] = adata_sub.obs[var]
            else:
                # Use expression counts, handling sparse and dense matrices
                data[var] = adata_sub[:, var].X
                if hasattr(data[var], "toarray"):
                    data[var] = data[var].toarray()
                if hasattr(data[var], "flatten"):
                    data[var] = data[var].flatten()
        activity_matrix = pandas.DataFrame(data = data, index = index)

    elif layer == "obsm":
        activity_matrix = pandas.DataFrame(
            adata_sub.obsm["score_ulm"], 
            index=adata_sub.obs_names
        )
    elif layer == "obs":
        # Fetch from observation metadata
        activity_matrix = pandas.DataFrame(
            data={var: adata_sub.obs[var] for var in tfs}, 
            index=adata_sub.obs_names
        )
    elif layer == "X":
        # Fetch from main expression matrix, converting sparse to dense
        X_data = adata_sub[:, tfs].X
        if hasattr(X_data, "toarray"):
            X_data = X_data.toarray()
        activity_matrix = pandas.DataFrame(
            X_data, 
            columns=tfs, 
            index=adata_sub.obs_names
        )
    elif layer == "raw":
        # Fetch from raw (unprocessed) expression matrix
        X_data = adata_sub.raw.to_adata()[:, tfs].X
        if hasattr(X_data, "toarray"):
            X_data = X_data.toarray()
        activity_matrix = pandas.DataFrame(
            X_data, 
            columns=tfs, 
            index=adata_sub.obs_names
        )
    else:
        # Fetch from a custom named layer
        layer_data = adata_sub[:, tfs].layers[layer]
        if hasattr(layer_data, "toarray"):
            layer_data = layer_data.toarray()
        activity_matrix = pandas.DataFrame(
            layer_data, 
            columns=tfs, 
            index=adata_sub.obs_names
        )
    
    # Create a new AnnData object with TF/activity values and inherited metadata
    tf_adata = scanpy.AnnData(X=activity_matrix)
    tf_adata.obs = adata_sub.obs.copy()
    tf_adata.obsm = adata_sub.obsm.copy()
    tf_adata.uns = adata_sub.uns.copy()
    tf_adata.uns["term_states_fwd_colors"] = term_state_colors
    
    # Restore original categories to prevent missing lineages when subsetting drops states
    tf_adata.obs["term_states_fwd"] = pandas.Categorical(
        adata_sub.obs["term_states_fwd"],
        categories=adata.obs["term_states_fwd"].cat.categories
    )
    
    # Reconstruct the lineage object using the correct forward fate probabilities key
    if "lineages_fwd" in tf_adata.obsm:
        tf_adata.obsm["lineages_fwd"] = cellrank.Lineage.from_adata(
            tf_adata, 
            kind="fate_probs"
        )

    # Fit a GAM model if one is not provided
    if model is None:
        model = cellrank.models.GAMR(
            tf_adata, 
            n_knots=n_knots, 
            smoothing_penalty=smoothing_penalty
        )

    if lineages == "all":
        lineages = adata.obs["term_states_fwd"].cat.categories
        
    cellrank.pl.gene_trends(
        tf_adata,
        model=model,
        genes=tfs,
        data_key="X",
        same_plot=True,
        lineages=lineages,
        ncols=n_columns,
        time_key=time_key,
        hide_cells=True,
        weight_threshold=(1e-3, 1e-3),
    )

def plot_differential_expression(differential_expression_dataframe, genes, groups="all"):
    """Plot differential expression log fold changes with significance annotations.

    Parameters
    ----------
    differential_expression_dataframe : pandas.DataFrame
        DataFrame containing log fold changes and adjusted p-values for multiple groups.
    genes : list[str]
        Gene names to include in the heatmap.
    groups : list[str] or "all", default "all"
        Groups to include in the plot. If ``"all"``, groups are inferred from column names.

    Returns
    -------
    None
        The function displays the heatmap directly.
    """

    if groups == "all":
        groups = differential_expression_dataframe.columns[differential_expression_dataframe.columns.str.contains("LFC")].str.split(" ").str[:-1].str.join(sep=" ")

    # Select and order genes for plotting based on Reprogramming LFC
    plot_data = differential_expression_dataframe.loc[genes].sort_values("Reprogramming LFC")

    # Separate log fold change and adjusted p-value columns for plotting
    lfc_columns = [f"{group} LFC" for group in groups]
    padj_columns = [f"{group} p adjusted" for group in groups]

    lfc_matrix = plot_data[lfc_columns]
    padj_matrix = plot_data[padj_columns]

    # Set clean column labels for plotting
    lfc_matrix.columns = groups

    # Convert adjusted p-values to a float array and replace NaN with 1
    padj_array = padj_matrix.to_numpy(dtype=float)
    numpy.nan_to_num(padj_array, nan=1, copy=False)

    # Build annotation strings by progressively adding asterisks for increasing significance
    # This creates overlapping "*" symbols for multi-asterisk notation (e.g., " * *" for p < 0.01)
    significance_annotations = numpy.full(padj_array.shape, "", dtype=numpy.dtypes.StrDType)
    significance_annotations[padj_array < 0.05] += "*"
    significance_annotations[padj_array < 0.01] += " *"
    significance_annotations[padj_array < 0.001] += " *"

    # Initialize the figure canvas for the heatmap
    pyplot.figure(figsize=(8, 6))

    # Plot the heatmap using a divergent colormap centered at zero
    seaborn.heatmap(
        data=lfc_matrix,
        annot=significance_annotations,
        fmt="s",
        cmap="coolwarm",
        center=0,
        linewidths=0.5,
        cbar_kws={"label": "Log2 Fold Change"}
    )

    # Label axes and finalize layout
    pyplot.xlabel("Condition")
    pyplot.ylabel("Gene")
    pyplot.tight_layout()
    pyplot.show()