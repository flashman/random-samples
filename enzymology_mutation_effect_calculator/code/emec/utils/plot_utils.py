from typing import Any, Dict, Optional

from Bio.Align import MultipleSeqAlignment
import matplotlib.pylab as plt
import numpy as np
import pandas as pd
import seaborn as sns
from statsmodels.base.model import Results

from . import alignment_utils as au
from . import statsmodels_utils as su


DEFAULT_H = 10
DEFAULT_W = 10
DEFAULT_SIZE = (DEFAULT_W, DEFAULT_H)
DEFAULT_COEF_PALETTE = "coolwarm_r"
DEFAULT_COUNT_PALETTE = "gist_yarg"


def plot_model_coefficients(
    results: Results, p_threshold=None, saveas: str = None
) -> plt.Figure:
    """
    Create a bar chart of a fitted model's feature coefficients.

    The method assumes that the provided dataframe is indexed by feature name.

    Args:
        results: the fitted Results object.
        p_threshold: coefficient significance threshold.  If none, all coefficients are returned.
        saveas: optionally, save the figure to the given file path.

    Returns:
        a barplot Figure.
    """
    df = su.get_model_coefficients(results, p_threshold=p_threshold)

    df.index.name = "feature"
    df = df.reset_index()

    df.sort_values("coef", ascending=False, inplace=True)

    h = len(df)
    w = DEFAULT_W
    fig, ax = plt.subplots(1, 1, figsize=(w, h))

    err = df["std err"].values
    pallete = _get_colors(
        df.coef, keys=df.feature, palette=DEFAULT_COEF_PALETTE, center=0
    )
    sns.barplot(
        y="feature",
        x="coef",
        data=df,
        orient="h",
        ax=ax,
        ci=None,
        xerr=err,
        palette=pallete,
    )

    ax.xaxis.set_ticks_position("top")

    if saveas:
        fig.savefig(saveas)

    return fig


def plot_model_coefficient_interactions(
    results: Results,
    drop_na: Optional[bool] = False,
    p_threshold: Optional[float] = None,
    saveas: Optional[str] = None,
) -> plt.Figure:
    """
    Plot a heatmap of interaction term coefficients for a fitted model.

    First order coefficients are also plotted along the diagonal for comparison.

    Args:
        results: the fitted Results object.
        drop_na: if True, drop parameters lacking value estimates.
        p_threshold: coefficient significance threshold.  If none, all coefficients are returned.
        saveas: optionally, save the figure to the given file path.

    Returns:
        a heatmap Figure.
    """
    # Get matrix of first and second order terms.
    df = su.get_model_coefficients_matrix(
        results, drop_na=drop_na, p_threshold=p_threshold
    )

    # Make it big...
    width = df.shape[0] / 2
    height = df.shape[0] / 2

    # Plot matrix as heatmap.
    fig, ax = plt.subplots(1, 1, figsize=(width, height))
    sns.heatmap(
        df,
        linewidths=0.01,
        linecolor="grey",
        square=True,
        cmap=DEFAULT_COEF_PALETTE,
        cbar_kws={"label": "Coefficient"},
        ax=ax,
        center=0,
    )

    # Set missing value color.
    ax.set_facecolor("white")

    # Label axes.
    ax.set_xlabel("Feature 1")
    ax.set_ylabel("Feature 2")

    if saveas:
        fig.savefig(saveas)

    return fig


def plot_model_fit(results: Results, saveas: str = None) -> plt.Figure:
    """
    Plot the in-sample prediction vs that actual training data for a `statsmodels` Results object.

    Args:
        results: the fitted Results object.
        saveas: optionally, save the figure to the given file path.

    Returns:
        a lmplot Figure.
    """
    sns.set(style="whitegrid")

    X = results.model.data.exog
    y = results.model.data.endog
    y_pred = results.predict(X)

    df = pd.DataFrame({"actual": y, "in-sample prediction": y_pred})

    g = sns.lmplot("actual", "in-sample prediction", df, height=DEFAULT_H)

    if saveas:
        g.savefig(saveas)

    return g.fig


def plot_model_coefficient_support(
    results: Results, drop_na=False, p_threshold=None, saveas: str = None
) -> plt.Figure:
    """
    Plot the performance distributions
    """
    df = su.get_model_coefficient_support(
        results, p_threshold=p_threshold, drop_na=drop_na
    )
    cdf = su.get_model_coefficients(results, p_threshold=p_threshold, drop_na=drop_na)

    # Initialize the figure.
    h = df.feature.nunique() / 2
    w = DEFAULT_W
    fig, ax = plt.subplots(figsize=(w, h))

    # Plot the per-mutation value distributions.
    pallete = _get_colors(
        cdf.coef, keys=cdf.index, palette=DEFAULT_COEF_PALETTE, center=0
    )
    sns.boxplot(
        x="value", y="feature", data=df, whis=[0, 100], width=0.6, palette=pallete
    )
    sns.stripplot(x="value", y="feature", data=df, linewidth=0, size=4, color=".3")

    # Tweak the visual presentation
    ax.xaxis.grid(True)
    sns.despine(trim=True, left=True)
    ax.set(ylabel="Mutation Name", xlabel=results.model.endog_names)

    if saveas:
        fig.savefig(saveas)

    return fig


def plot_library_coverage(
    msa: MultipleSeqAlignment, mincount: int = 1, saveas: str = None
) -> plt.Figure:
    """
    Plot the position-by-amino acid coverage for the given multi-sequence alignment.

    Args:
        msa: multiple sequence alignment, with first record assumed to be the reference.
        mincount: only include positions that contain at least mincount mutations.
        saveas: optionally, save the figure to the given file path.

    Returns:
        a heatmap Figure.
    """
    df = au.msa_to_coverage_df(msa, mincount).T

    width = df.shape[1] / 2
    height = df.shape[0] / 2

    fig, ax = plt.subplots(1, 1, figsize=(width, height))

    ax = sns.heatmap(
        df,
        linewidths=0.01,
        linecolor="grey",
        square=True,
        annot=True,
        fmt="d",
        cbar_kws={"label": "Count"},
        cmap=DEFAULT_COUNT_PALETTE,
        ax=ax,
        vmin=0,
    )

    if saveas:
        fig.savefig(saveas)

    return fig


def plot_mutation_co_occurence(
    msa: MultipleSeqAlignment, mincount: int = 1, saveas: str = None
) -> plt.Figure:
    """
    Plot mutation co-occurrence between sequences in the multiple sequence alignment.

    Args:
        msa: multiple sequence alignment, with first record assumed to be the reference.
        mincount: include mutations that appear at least this many times.
        saveas: optionally, save the figure to the given file path.

    Returns:
        a mutation correlation Figure.
    """
    # Get composition df.
    comp_df = au.msa_to_mutation_composition_df(msa)

    # Drop low frequency mutations.
    drop_cols = comp_df.columns[comp_df.sum(axis=0) < mincount]
    comp_df.drop(drop_cols, axis=1, inplace=True)

    # Get co-occurrence count.
    cocount_df = comp_df.T.dot(comp_df)

    # Plot co-occurrence.
    width = cocount_df.shape[1] / 1.5
    height = cocount_df.shape[0] / 2

    fig, ax = plt.subplots(1, 1, figsize=(width, height))
    sns.heatmap(
        cocount_df,
        linewidths=0.01,
        linecolor="grey",
        square=True,
        annot=True,
        fmt="d",
        cbar_kws={"label": "Count"},
        cmap=DEFAULT_COUNT_PALETTE,
        ax=ax,
        vmin=0,
    )

    ax.set_xlabel("Mutation Name")
    ax.set_ylabel("Mutation Name")

    if saveas:
        fig.savefig(saveas)

    return fig


def plot_sequence_distance(msa: MultipleSeqAlignment, saveas: str = None) -> plt.Figure:
    """
    Plot the pairwise hamming distance between all sequences in the multiple sequence alignment.

    Args:
        msa: multiple sequence alignment, with first record assumed to be the reference.
        saveas: optionally, save the figure to the given file path.

    Returns:
        a sequence similarity Figure.
    """
    sim_df = au.msa_to_hamming_distance(msa)

    g = sns.clustermap(sim_df, figsize=DEFAULT_SIZE)

    g.ax_heatmap.set_xlabel("Sequence ID")
    g.ax_heatmap.set_ylabel("Sequence ID")

    if saveas:
        g.savefig(saveas)

    return g.fig


def _get_colors(
    values: np.array, palette: str = None, keys: np.array = None, center: float = None
) -> Dict[Any, np.array]:
    """
    Returns colors associated with the input values using the specified palette.

    Args:
        values: an array of values to map to colors
        palette: the name of the color palette to use. If None, the current default is used.
        keys: an array of value keys.  If provided, the color map is returned as a dict of
            key-color pairs.
        center: if provided, the palette is centered at this value.  Note that this can reduce the
            resulting color range.

    Returns:
        either a list or dict of colors corresponding to the input values.
    """
    # Normalize the values to range [0, 1]
    if center is not None:
        # Locate the center value at 0.5
        normalized = (values - center) / max(abs(values - center)) / 2 + 0.5
    else:
        # Use full color range.
        normalized = (values - min(values)) / (max(values) - min(values))
    # Convert to indices
    indices = np.ceil(normalized * (len(values) - 1)).astype(np.int32)
    # Use the indices to get the colors
    palette = sns.color_palette(palette, len(values))
    colors = np.array(palette).take(indices, axis=0)
    # Construct color map.
    if keys is not None:
        return dict(zip(keys, colors))
    else:
        return colors
