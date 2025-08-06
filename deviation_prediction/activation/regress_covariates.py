import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.formula.api import ols

def regress_covariates(data, covariates, covariate_names=None, categorical_covariates=None):
    """
    regress_covariates

    :parameters
    - data: shape (n_samples, n_features)
    - covariates: shape (n_samples, n_covariates)
    - covariate_names: list of str
    - categorical_covariates: list of str

    return:
    - resid_data: shape (n_samples, n_features)
    """

    if isinstance(covariates, np.ndarray):
        cov_df = pd.DataFrame(covariates)
        if covariate_names is None:
            covariate_names = [f'cov_{i}' for i in range(covariates.shape[1])]
        cov_df.columns = covariate_names
    else:
        cov_df = covariates.copy()
        covariate_names = cov_df.columns.tolist()

    if isinstance(data, np.ndarray):
        data_df = pd.DataFrame(data, columns=[f'feature_{i}' for i in range(data.shape[1])])
    else:
        data_df = data.copy()

    formula_parts = []
    for name in covariate_names:
        if categorical_covariates and name in categorical_covariates:
            formula_parts.append(f'C({name})')
        else:
            formula_parts.append(name)
    formula_base = ' + '.join(formula_parts)
    full_formula = f'data_col ~ {formula_base}'

    resid_data = np.zeros_like(data, dtype=float)

    for col_idx in range(data_df.shape[1]):
        df_temp = cov_df.copy()
        df_temp['data_col'] = data_df.iloc[:, col_idx]

        try:
            model = ols(full_formula, data=df_temp).fit()
            resid_data[:, col_idx] = model.resid
        except Exception as e:
            print(f"Feature {col_idx} regression failed: {e}")
            resid_data[:, col_idx] = np.nan

    return resid_data