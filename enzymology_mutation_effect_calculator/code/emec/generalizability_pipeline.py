import logging
from typing import List, Optional, Tuple

import matplotlib.pylab as plt
import numpy as np
import pandas as pd
from pydantic import BaseModel
import seaborn as sns
from sklearn.model_selection import GroupShuffleSplit

from emec.pipeline import fit_elastic_net_cv
from emec.utils.pipeline_utils import get_degree
from emec.utils.plot_utils import DEFAULT_H, DEFAULT_W, DEFAULT_SIZE

GENERALIZABILITY_PARAMS_KEY = "generalizability_analysis_parameters"


class GeneralizabilityParams(BaseModel):
    """Parameters for running the model generalizability pipeline.

    Attributes:
        min_train_frac: minimum train/test splitting fraction.
        max_train_frac: maximum train/test splitting fraction.
        n_steps: number of train/test splitting fractions, including the min and max fractions.
        n_split: number of random splits to perform for each train/test splitting fraction.
    """

    min_train_frac: float = 0.2
    max_train_frac: float = 0.8
    n_steps: int = 4
    n_splits: int = 8


def run_pipeline(
    results, params: "GeneralizabilityParams", model_fit_params: dict
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Run pipeline to evaluate how model fit and coefficients change with the the size of the training
    dataset.

    The pipeline performs the following steps:

    1. Select train/test split ratios based on provided parameters: min_train_frac, max_train_frac,
    n_steps.
    2. For each split ratio, use SKLearn's GroupShuffleSplit class to generate `n_splits` random
    train/test splits of the full dataset, in proportion to the split ration, taking are not to
    split up strains of the same id.
    3. For each training set, fit the ElasticNetCV model. Note that this fitting step includes
    hyper-parameter tuning based on the user-provided pipeline parameters.
    4. Score (using R2) the fitted model on both the training set and the test set.

    In addition to in-sample and out-of-sample models scores, for each fitted model, we also store
    the optimized hyper-parameters and the full set of model coefficients.

    Args:
        results: the fitted Statsmodels Results object, returned by the mutation effect evaluation
            pipeline.
        params: parameters to control the model generalizability pipeline. See
            GeneralizabilityParams for parameter details.
        model_fit_params: a dictionary of modeling parameters.  Theses should be the same as

    Returns:
        a dataframe of model fit results for each iteration.
        a dataframe of model coefficients for each iteration.
        a dataframe of training feature counts for each iteration.

    """
    # Get full training data from results.
    # Note that SKLearn handles the intercept differently than Statsmodels.
    coef_cols = results.model.data.orig_exog.columns
    X = results.model.data.orig_exog.iloc[:, 1:]
    y = results.model.data.orig_endog

    # A safe place to store trial results.
    model_fits = []
    model_coefs = []
    feature_counts = []

    # Fit the ElasticNetCV model using random sub-samples of the full training dataset.
    for i, train_frac in enumerate(
        np.linspace(params.min_train_frac, params.max_train_frac, params.n_steps)
    ):
        logging.warning(
            f"Fitting {params.n_splits} models using train/test splitting fraction {train_frac}..."
        )

        # Perform n_split random train/test splits of the full dataset.
        # But don't split up samples with the same strain id to avoid information leakage.
        gss = GroupShuffleSplit(n_splits=params.n_splits, test_size=1 - train_frac)
        for j, (train_inds, test_inds) in enumerate(gss.split(X, y, X.index)):
            logging.warning(".")

            # Prepare train and test data.
            X_train = X.iloc[train_inds]
            X_test = X.iloc[test_inds]
            y_train = y.iloc[train_inds]
            y_test = y.iloc[test_inds]

            # Fit the model.
            encv = fit_elastic_net_cv(X_train, y_train, model_fit_params)

            # Store model fit trial result.
            mf = {
                "group_id": i,
                "split_id": j,
                "train_frac": train_frac,
                "train_size": len(X_train),
                "total_size": len(X),
                "alpha": encv.alpha_,
                "l1_ratio": encv.l1_ratio_,
                "train_R2": encv.score(X_train, y_train),
                "test_R2": encv.score(X_test, y_test),
            }

            # Compute other library statistics such as observations per feature.
            feature_degrees = X_train.columns.to_series().apply(get_degree)
            for d in feature_degrees.unique():
                deg_feat_cnts = X_train.loc[:, feature_degrees == d].sum(axis=0)
                mf[f"train_coverage_deg{d}_mean"] = deg_feat_cnts.mean()
                mf[f"train_coverage_deg{d}_std"] = deg_feat_cnts.std()

            model_fits.append(mf)

            # Store model coefficient trial result.
            # Manually add the intercept for easier comparison with Statsmodel results.
            model_coefs.append(np.concatenate(([encv.intercept_], encv.coef_)))

            # Store model fit coefficient counts.
            feature_counts.append(X_train.sum(axis=0))

        # Prepare model fit summary table.
        model_fit_df = pd.DataFrame(model_fits)
        # Training size can differ slightly from split to split. Use mean across groups.
        model_fit_df.train_size = model_fit_df.groupby("group_id").train_size.transform(
            "mean"
        )

        # Prepare model coefficients summary table
        model_coef_df = pd.DataFrame(
            model_coefs, columns=coef_cols, index=model_fit_df.index
        )
        # Attach group_id, split_id, and train_size for easy summary later.
        model_coef_df = pd.concat(
            [model_fit_df[["group_id", "split_id", "train_size"]], model_coef_df],
            axis=1,
        )

        # Prepare model coefficients counts table
        feature_counts_df = pd.DataFrame(feature_counts, index=model_fit_df.index)
        # Attach group_id, split_id, and train_size for easy summary later.
        feature_counts_df = pd.concat(
            [model_fit_df[["group_id", "split_id", "train_size"]], feature_counts_df],
            axis=1,
        )

    return model_fit_df, model_coef_df, feature_counts_df


# TODO(flash): Find a better home for this plotting methods.


def plot_model_fit_vs_library_size(
    model_fit_df: pd.DataFrame, saveas: Optional[str] = None
) -> plt.Figure:
    """
    Plot how in-sample and out-of-sample model fit (R-squared) vary with library size.

    Args:
        model_fit_df: a dataframe of model fit results returned by run_pipeline().
        saveas: optionally, save the figure to the given file path.

    Returns:
        the model fit Figure.
    """
    # Melt train/test fit columns.
    plot_df = model_fit_df.melt(
        "train_size", ["train_R2", "test_R2"], "split", "R-squared"
    )
    plot_df["split"] = plot_df["split"].str.split("_", expand=True)[0]

    # Initialize the figure.
    fig, ax = plt.subplots(1, 1, figsize=DEFAULT_SIZE)

    # Plot R-squared by training size, and color by train/test.
    sns.scatterplot(
        x="train_size", y="R-squared", hue="split", data=plot_df, legend=None, ax=ax
    )
    sns.lineplot(x="train_size", y="R-squared", hue="split", data=plot_df, ax=ax)

    # Tweak the axes.
    ax.set_xlabel("Library Size")

    if saveas:
        fig.savefig(saveas)

    return fig


def plot_model_hyperparameters_vs_library_size(
    model_fit_df: pd.DataFrame, saveas: Optional[str] = None
) -> plt.Figure:
    """
    Plot how model hyperparameters, alpha and l1_ratio, vary with library size.

    Args:
        model_fit_df: a dataframe of model fit results returned by run_pipeline().
        saveas: optionally, save the figure to the given file path.

    Returns:
        the hyperparameters Figure.
    """
    # Initialize the figure.
    fig, axes = plt.subplots(1, 2, figsize=(DEFAULT_W * 2, DEFAULT_H))

    # Plot alpha by training size.
    sns.scatterplot(x="train_size", y="alpha", data=model_fit_df, ax=axes[0])
    sns.lineplot(x="train_size", y="alpha", data=model_fit_df, ax=axes[0])

    # Plot l1_ratio by training size.
    sns.scatterplot(x="train_size", y="l1_ratio", data=model_fit_df, ax=axes[1])
    sns.lineplot(x="train_size", y="l1_ratio", data=model_fit_df, ax=axes[1])

    # Tweak the axes.
    axes[1].set_ylim(-0.1, 1.1)
    axes[0].set_xlabel("Library Size")
    axes[1].set_xlabel("Library Size")

    if saveas:
        fig.savefig(saveas)

    return fig


def plot_model_coefficients_vs_library_size(
    model_coefs_df: pd.DataFrame,
    coefficients: Optional[List[str]] = None,
    drop_zeros: bool = True,
    saveas: Optional[str] = None,
) -> plt.Figure:
    """
    Plot how model coefficients vary with library size.

    Args:
        model_coefs_df: a dataframe of model coefficients returned by run_pipeline().
        coefficients: optionally, only display coefficients in this list.
        drop_zeros: optionally, remove coefficients whose value is always zero.
        saveas: optionally, save the figure to the given file path.

    Returns:
        the coefficient Figure.
    """
    # Melt coefficient values and apply optionally filters.
    COEF_COL_START = 3
    plot_df = model_coefs_df.melt(
        "train_size", model_coefs_df.columns[COEF_COL_START:], "coef"
    )
    if coefficients is not None:
        plot_df = plot_df[plot_df.coef.isin(coefficients)]
    if drop_zeros:
        non_zero_coefs = plot_df[plot_df.value != 0].coef.unique()
        plot_df = plot_df[plot_df.coef.isin(non_zero_coefs)]

    # Initialize the figure.
    H = plot_df.coef.nunique() * plot_df.train_size.nunique() / 4
    fig, ax = plt.subplots(1, 1, figsize=(DEFAULT_W, H))

    # Plot coefficient values by training size.
    sns.barplot(y="coef", x="value", data=plot_df, hue="train_size", orient="h", ax=ax)

    # Tweak the axes.
    [ax.axhline(y + 0.5, linestyle="--", color="black") for y in ax.get_yticks()[:-1]]
    ax.set_ylabel("Mutation Name")

    if saveas:
        fig.savefig(saveas)

    return fig


def plot_feature_counts_vs_library_size(
    feature_counts_df: pd.DataFrame, saveas: Optional[str] = None
) -> plt.Figure:
    """
    Plot how training data mutation feature counts vary with library size.

    Args:
        feature_counts_df: a dataframe of training data feature counts returned by run_pipeline().
        saveas: optionally, save the figure to the given file path.

    Returns:
        the feature counts Figure.
    """
    # Melt counts and add feature degree.
    COEF_COL_START = 3
    plot_df = feature_counts_df.melt(
        id_vars="train_size", value_vars=feature_counts_df.columns[COEF_COL_START:]
    )
    plot_df["degree"] = plot_df.variable.apply(get_degree)

    # Plot feature counts by training size.
    fig, ax = plt.subplots(1, 1, figsize=DEFAULT_SIZE)
    sns.boxenplot(x="train_size", y="value", hue="degree", data=plot_df, ax=ax)

    # Tweak the axes.
    ax.set_xlabel("Library Size")
    ax.set_ylabel("Mutation Feature Count")
    ax.legend_.set_title("Mutation Feature Degree")

    if saveas:
        fig.savefig(saveas)

    return fig


def plot_feature_counts_vs_model_coefficients(
    model_coefs_df: pd.DataFrame,
    feature_counts_df: pd.DataFrame,
    coefficients: Optional[List[str]] = None,
    drop_zeros: bool = True,
    saveas: Optional[str] = None,
) -> plt.Figure:
    """
    Plot how model coefficients vary with mutation feature counts as library size is varied.

    Args:
        model_coefs_df: a dataframe of model coefficients returned by run_pipeline().
        feature_counts_df: a dataframe of training data feature counts returned by run_pipeline().
        coefficients: optionally, only display coefficients in this list.
        drop_zeros: optionally, remove coefficients whose value is always zero.
        saveas: optionally, save the figure to the given file path.

    Returns:
        the Figure.
    """
    # Melt feature counts and coefficient values.
    DATA_COL_START = 3
    model_coefs_melted = model_coefs_df.melt(
        id_vars=["group_id", "split_id", "train_size"],
        value_vars=model_coefs_df.columns[DATA_COL_START:],
        var_name="coef_name",
        value_name="coef_value",
    )
    feature_counts_melted = feature_counts_df.melt(
        id_vars=["group_id", "split_id", "train_size"],
        value_vars=feature_counts_df.columns[DATA_COL_START:],
        var_name="coef_name",
        value_name="coef_count",
    )

    # Merge melted data for plotting.
    plot_df = pd.merge(
        model_coefs_melted,
        feature_counts_melted,
        on=["group_id", "split_id", "train_size", "coef_name"],
    )

    # Feature counts can differ slightly from split to split. Use group mean.
    plot_df = (
        plot_df.groupby(["group_id", "train_size", "coef_name"]).mean().reset_index()
    )

    # Add feature degree column.
    plot_df["degree"] = plot_df["coef_name"].apply(get_degree)

    # Apply user-provided data filters.
    if coefficients is not None:
        plot_df = plot_df[plot_df["coef_name"].isin(coefficients)]
    if drop_zeros:
        non_zero_coefs = plot_df[plot_df["coef_value"] != 0]["coef_name"].unique()
        plot_df = plot_df[plot_df["coef_name"].isin(non_zero_coefs)]

    # Prep feature labels data for plotting. Labels go next to right-most data points.
    labels_df = plot_df[plot_df["group_id"] == max(plot_df["group_id"])]

    # Initialize the figure.
    fig, ax = plt.subplots(1, 1, figsize=(16, 12))

    # Plot feature counts by feature coefficient as line and points.
    sns.lineplot(
        data=plot_df,
        x="coef_count",
        y="coef_value",
        hue="coef_name",
        style="degree",
        legend=False,
        ax=ax,
    )
    ax = sns.scatterplot(
        data=plot_df,
        x="coef_count",
        y="coef_value",
        size="train_size",
        style="degree",
        legend="full",
        ax=ax,
    )

    # Add labels to plot.
    for _, row in labels_df.iterrows():
        ax.text(row["coef_count"] + 1, row["coef_value"], row["coef_name"])

    # Tweak the axes.
    ax.set_xlabel("Mutation Feature Count (Mean)")
    ax.set_ylabel("Mutation Feature Coefficient Value (Mean)")

    # Tweak the legend.
    handles, labels = ax.get_legend_handles_labels()
    labels = _replace(labels, "train_size", "Library Size")
    labels = _replace(labels, "degree", "Feature Degree")
    ax.legend(handles, labels, loc="upper left", bbox_to_anchor=(1, 1))

    if saveas:
        fig.savefig(saveas)

    return fig


def _replace(lst, old, new):
    """Helper method replaces all instances of an old value with a new value in the given list"""
    return [new if v == old else v for v in lst]
