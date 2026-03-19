from typing import Optional

import numpy as np
import pandas as pd
from statsmodels.base.model import Results

from emec.utils.alignment_utils import sort_mutations


def get_model_summary(results: Results, props: Optional[dict] = None) -> pd.DataFrame:
    """
    Get summary statistics for a `statsmodels` Results object.

    The summary dataframe includes model properties such as:
        Model, Method, Date, Time, No. Observations, R-squared, F-statistic, Log-Likelihood, etc.

    Args:
        results: the fitted Results object
        props: any additional model properties to include with the model summary

    Returns:
        a dataframe containing the model summary.
    """
    # Load data.
    df = pd.DataFrame(results.summary().tables[0].data)

    # Stack cols 1-2 and 3-4
    df1 = df.loc[:, 0:1].copy()
    df2 = df.loc[:, 2:3].copy()
    df2.columns = [0, 1]
    df = pd.concat([df1, df2], axis=0, ignore_index=True)

    # Fix values
    df.replace(r"^[ ]+$", np.nan, regex=True, inplace=True)
    df.dropna(inplace=True)

    # Fix index and column headers.
    df.set_index(0, inplace=True)
    df.index = df.index.str.strip().str.strip(":")
    df.columns = ["value"]
    df.index.name = None

    # Attach any additional parameters to the dataframe
    if props:
        for k, v in props.items():
            df.loc[k, "value"] = v

    return df


def get_model_coefficients(
    results: Results, drop_na: bool = True, p_threshold: Optional[float] = None
) -> pd.DataFrame:
    """
    Get fitted model parameter coefficients for a `statsmodels` Results object.

    Args:
        results: the fitted Results object
        drop_na: if True, drop parameters lacking value estimates.
        p_threshold: coefficient significance threshold.  If none, all coefficients are returned.

    Returns:
        a dataframe containing model parameter coefficients.
    """
    # Load data.
    df = pd.DataFrame(results.summary().tables[1].data)

    # Fix index and columns headers
    df.set_index(0, inplace=True)
    df = df.rename(columns=df.iloc[0]).drop(df.index[0])
    df.index.name = "feature"
    df.columns = df.columns.str.strip()
    df.index = df.index.str.strip()

    # Fix nans and data types.
    df.replace(r".*nan", np.nan, regex=True, inplace=True)
    df = df.astype(float)

    # Add feature observation counts.
    X_df = results.model.data.orig_exog
    df["count"] = X_df.sum(axis=0)

    # Apply filters.
    if drop_na:
        df.dropna(inplace=True)

    if p_threshold:
        df = df[df["P>|t|"] < p_threshold]

    return df


def get_model_coefficients_matrix(
    results: Results,
    drop_na: Optional[bool] = False,
    p_threshold: Optional[float] = None,
):
    """
    Get 1st and 2nd order model coefficients for a `statsmodels` Results object, as a
    symmetric matrix.

    Here the 2nd order interaction terms appear in the off-diagonal positions, while 1st order
    terms appear along the diagonal.  This assumes that the 1st order features are themselves a
    one-hot encoding.

    The method additionally assumes that second order feature names are space delimited
    i.e. the feature corresponding to the second order interaction betwen "feature_i" and
    "feature_j" is called "feature_i feature_j".

    Args:
        results: the fitted Results object
        drop_na: if True, drop parameters lacking coefficient value estimates.
        p_threshold: filter coefficient by significance threshold.  If none, all coefficients are
            returned.

    Returns:
        a dataframe containing 1st and 2nd order model coefficients.
    """
    # Get the model coefficient dataframe the results, with optional filtering.
    cdf = get_model_coefficients(results, drop_na=drop_na, p_threshold=p_threshold)

    # Isolate first and second order coefficients.
    # Assumes that second order feature names are space delimited e.g. "feature_i feature_j"
    first_order_coeffs = cdf.index.str.match(r"^\w\d+\w$")
    second_order_coeffs = cdf.index.str.match(r"^\w\d+\w \w\d+\w$")

    # Prepare first order coefficients for the matrix construction.
    fo_df = cdf[first_order_coeffs].reset_index()
    fo_df["feature 1"] = fo_df.feature
    fo_df["feature 2"] = fo_df.feature

    # Prepare second order coefficients for matrix construction, provided they exist.
    if sum(second_order_coeffs) > 0:
        # Assumes that second order features take the form "feature1 feature2"
        so_df = cdf[second_order_coeffs].reset_index()
        so_df[["feature 1", "feature 2"]] = so_df.feature.str.split(expand=True)

        # Lets make our matrix symmetric.
        sso_df = cdf[second_order_coeffs].reset_index()
        sso_df[["feature 2", "feature 1"]] = sso_df.feature.str.split(expand=True)

        # Build the coefficients matrix.
        df = pd.concat([fo_df, so_df, sso_df], axis=0)
    else:
        # Build the coefficients matrix.
        df = pd.concat([fo_df], axis=0)

    matrix = df.pivot(index="feature 1", columns="feature 2", values="coef")

    # Order rows and columns based on original mutation ordering.
    ordered = sort_mutations(matrix.index)
    matrix = matrix.loc[ordered, ordered]

    return matrix


def get_model_coefficient_support(
    results: Results, drop_na: bool = False, p_threshold: Optional[float] = None
) -> pd.DataFrame:
    """
    Get supporting data for the coefficients of a fitted model.

    The dataframe consists of the following columns:
        feature: the name of the coefficient
        id: the id of the supporting datum
        value: the endogenous variable value
        name: the endogenous variable name

    Args:
        results: the fitted statsmodels Results object
        drop_na: if True, drop parameters lacking coefficient value estimates.
        p_threshold: filter coefficient by significance threshold.  If none, all coefficients are
            returned.

    Returns:
        a dataframe of supporting data for each feature.
    """

    X = results.model.data.orig_exog
    y = results.model.data.orig_endog
    cdf = get_model_coefficients(results, drop_na=drop_na, p_threshold=p_threshold)

    dfs = []
    for c in cdf.index:
        selected = (X.index[X[c] > 0]).drop_duplicates()
        values = y.loc[selected]
        df = pd.DataFrame(
            {"feature": c, "id": values.index, "value": values, "name": y.name}
        )
        dfs.append(df)

    df = pd.concat(dfs, axis=0, ignore_index=True)

    return df


def get_model_predictions_df(results: Results) -> pd.DataFrame:
    """
    Get a dataframe of actual and predicted values for the input data of a fitted model.

    The dataframe consists of the following columns:
        id: the id of the input datum
        actual: the actual value of the endogenous variable
        prediction: the value of the endogenous variable predicted by the model
        name: the endogenous variable name

    Args:
        results: the fitted statsmodels Results object.

    Returns:
        a dataframe of actual and predicted values.
    """
    X = results.model.data.exog
    y = results.model.data.endog
    y_pred = results.predict(X)

    df = pd.DataFrame(
        {
            "id": results.model.data.row_labels,
            "actual": y,
            "prediction": y_pred,
            "name": results.model.data.ynames,
        }
    )

    return df
