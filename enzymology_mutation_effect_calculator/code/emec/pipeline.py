import copy
import json
import logging

from Bio.Align import MultipleSeqAlignment
import pandas as pd
from sklearn.linear_model import ElasticNetCV
from sklearn.model_selection import KFold, GroupKFold
from sklearn.preprocessing import PolynomialFeatures
import statsmodels.api as sm
from statsmodels.base.elastic_net import RegularizedResults


from emec.utils import alignment_utils
from emec.utils import supervised_pca
from emec.utils import pipeline_utils


PIPELINE_PARAMS_KEY = "pipeline_parameters"
PIPELINE_PARAMS = {
    "precondition": False,
    "degree": 1,
    "min_observations": 1,
    "drop_redundant_features": True,
    "group_k_fold": {"n_splits": 5},
    "elastic_net_cv": {
        "l1_ratio": [0.1, 0.5, 0.7, 0.9, 0.99, 1.0],
        "selection": "random",
    },
}


def run_pipeline(
    msa: MultipleSeqAlignment, y: pd.Series, params: dict
) -> RegularizedResults:
    """
    Run the mutation effect evaluation pipeline.

    The pipeline consists of three steps:

        1. Featurize the multiple sequence alignment
           a. Optionally add polynomial terms.
           b. Optionally drop features with low observation counts.
           c. Optionally drop redundant features.

        2. Use sklearn's built-in cross validation capability to determine the optimal
        regularization parameters for ElasticNet regression.  To avoid genetic information leakage
        between train and test sets, we use the GroupKFold splitting class with strain ID as group
        identifier.

        3. Use the tuned alpha and L1 weight to fit Statsmodel's version of ElasticNet
        regression. The benefit of the Statsmodel version is that it provides confidence estimates
        of for model coefficient.

    Args:
        msa: a multiple sequence alignment where the first sequence is the reference to which all
            other sequences will be compared.
        y: n_sequence-length performance data.
        params: a dictionary of model parameters.  See PIPELINE_PARAMS for valid parameter schema.

    Returns:
        reg: the statsmodel result object.
        (alpha, L1_wt): the tuned alpha and L1_weight ElasticNet parameters.
        best-fit ElasticNetCV model.

    """
    # Apply user params to defaults.
    parameters = copy.deepcopy(PIPELINE_PARAMS)
    parameters.update(params)

    # Featurize the multiple sequence alignment by converting to the mutation composition format.
    features_df = alignment_utils.msa_to_mutation_composition_df(msa)

    # Drop null values from input data.
    y = y.dropna()

    # Drop any rows not represented in both datasets.
    joined_df = features_df.join(y, how="inner")
    X = joined_df[features_df.columns]
    y = joined_df[y.name]

    # Add higher order terms, but not powers.
    pf = PolynomialFeatures(
        degree=parameters["degree"], include_bias=False, interaction_only=True
    )
    X_t = pf.fit_transform(X)
    feature_names = pf.get_feature_names(input_features=X.columns)
    X = pd.DataFrame(X_t, columns=feature_names, index=X.index)

    # Drop all features with insufficient observations so as not to dwell on useless features.
    drop_cols = X.columns[X.sum(axis=0) < parameters["min_observations"]]
    if len(drop_cols) > 0:
        X.drop(drop_cols, axis=1, inplace=True)
        logging.warning(
            f"Dropped {len(drop_cols)} features with fewer than {parameters['min_observations']} "
            f"observations:\n {drop_cols.values}"
        )

    # Drop all features that are redundant with lower-order features.
    # First identify all such features.
    dedup_mapping = alignment_utils.get_covariate_mutations(X, min_covariates=2)
    n_duped_features = sum(len(cs) for cs in dedup_mapping.values())
    if n_duped_features > 0:
        # Always log presence of redundant features.
        logging.warning(
            f"Identified {n_duped_features} redundant features:\n"
            f"{json.dumps(dedup_mapping, indent=2)}"
        )

        # Optionally drop redundant features..
        if parameters["drop_redundant_features"]:
            # Identify redundant columns with higher degree than the first observed.
            # Here we make use of the assumption that features are ordered by degree.
            drop_cols = []
            for col, cvs in dedup_mapping.items():
                for c in cvs:
                    if pipeline_utils.get_degree(col) < pipeline_utils.get_degree(c):
                        drop_cols.append(c)

            # Drop the redundant columns.
            X.drop(drop_cols, axis=1, inplace=True)
            logging.warning(
                f"Dropped {len(drop_cols)} redundant features:\n {drop_cols}"
            )

    # Save for later.
    y_tune = y.copy()

    # Remove noisy features with preconditioning prior to hyper-parameter tuning.
    if parameters["precondition"]:
        spca = supervised_pca.SupervisedPCARegressor()
        spca.fit(X.values, y_tune)
        y_tune = pd.Series(spca.predict(X.values), index=y_tune.index, name=y_tune.name)

    # Perform hyper-parameter tuning to determine optimal ElasticNet regularization via cross
    # validation.  To that end, we fit ElasticNetCV to the data.
    encv = fit_elastic_net_cv(X, y_tune, parameters)

    # Store tuned parameters for later.
    fitted_parameters = {"alpha": encv.alpha_, "L1_wt": encv.l1_ratio_}

    # Fit Statsmodel ElasticNet using original raw data and tuned parameter.
    exog = sm.add_constant(X)
    m = sm.OLS(y, exog)
    reg = m.fit_regularized(refit=True, maxiter=1000, **fitted_parameters)

    # Return the results, tuned parameters, and fitted sklearn ElasticNet model (for comparison)
    return reg, fitted_parameters, encv


def fit_elastic_net_cv(X: pd.DataFrame, y: pd.Series, parameters: dict) -> ElasticNetCV:
    """ "
    Helper method to fit the ElasticNetCV model to the given input data using provided pipeline
    parameters.
    """
    # Manually shuffle inputs for cross validation. (The preffered splitter doesn't have a shuffle
    # option so we do it manually.)
    shuffled_idx = X.index.drop_duplicates().to_series().sample(frac=1)
    X_cv = X.loc[shuffled_idx]
    y_cv = y.loc[shuffled_idx]
    groups = X_cv.index.values

    # Initialize cross validation splitter.
    if "k_fold" in parameters:
        # Use deprecated KFold splitter.
        # TODO(flash): Remove this option once the other approach is fully vetted.
        logging.warning(
            "KFold cross validation is deprecated.  Use GroupKFold instead."
        )
        cv = KFold(**parameters["k_fold"]).split(X_cv, y_cv)
    elif "group_k_fold" in parameters:
        # Use GroupKFold splitter.
        cv = GroupKFold(**parameters["group_k_fold"]).split(X_cv, y_cv, groups)

    # Fit ElasticNet using cross validation for parameter tuning.
    encv = ElasticNetCV(cv=cv, **parameters["elastic_net_cv"])
    encv.fit(X_cv, y_cv)

    return encv
