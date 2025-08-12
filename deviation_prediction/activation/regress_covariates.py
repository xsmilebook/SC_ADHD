import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.formula.api import ols

def regress_covariates(train_data, train_covariates,
                       test_data=None, test_covariates=None,
                       covariate_names=None, categorical_covariates=None):
    """
    regress_covariates

    :parameters
    - train_data: shape (n_train_samples, n_features)
    - train_covariates: shape (n_train_samples, n_covariates)
    - test_data: shape (n_test_samples, n_features), optional
    - test_covariates: shape (n_test_samples, n_covariates), optional
    - covariate_names: list of str
    - categorical_covariates: list of str

    return:
    - train_resid_data: shape (n_train_samples, n_features)
    - test_resid_data: shape (n_test_samples, n_features) if test_data is provided
    - models: list of fitted OLS models (one per feature)
    """

    if isinstance(train_covariates, np.ndarray):
        train_cov_df = pd.DataFrame(train_covariates)
        if covariate_names is None:
            covariate_names = [f'cov_{i}' for i in range(train_covariates.shape[1])]
        train_cov_df.columns = covariate_names
    else:
        train_cov_df = train_covariates.copy()
        covariate_names = train_cov_df.columns.tolist()

    if isinstance(train_data, np.ndarray):
        train_data_df = pd.DataFrame(train_data, columns=[f'feature_{i}' for i in range(train_data.shape[1])])
    else:
        train_data_df = train_data.copy()

    if test_data is not None and test_covariates is not None:
        if isinstance(test_covariates, np.ndarray):
            test_cov_df = pd.DataFrame(test_covariates, columns=covariate_names)
        else:
            test_cov_df = test_covariates.copy()

        if isinstance(test_data, np.ndarray):
            test_data_df = pd.DataFrame(test_data, columns=[f'feature_{i}' for i in range(test_data.shape[1])])
        else:
            test_data_df = test_data.copy()
    else:
        test_cov_df = None
        test_data_df = None

    formula_parts = []
    for name in covariate_names:
        if categorical_covariates and name in categorical_covariates:
            formula_parts.append(f'C({name})')
        else:
            formula_parts.append(name)
    formula_base = ' + '.join(formula_parts)
    full_formula = f'data_col ~ {formula_base}'


    train_resid_data = np.zeros_like(train_data, dtype=float)
    test_resid_data = np.zeros_like(test_data, dtype=float) if test_data is not None else None

    models = []

    for col_idx in range(train_data_df.shape[1]):
        df_temp_train = train_cov_df.copy()
        df_temp_train['data_col'] = train_data_df.iloc[:, col_idx]

        try:
            model = ols(full_formula, data=df_temp_train).fit()
            models.append(model)

            train_resid_data[:, col_idx] = model.resid

            if test_data is not None and test_covariates is not None:
                df_temp_test = test_cov_df.copy()
                # residual = real_value - fitting_value
                df_temp_test['data_col'] = test_data_df.iloc[:, col_idx]
                test_pred = model.predict(df_temp_test)
                test_resid_data[:, col_idx] = test_data_df.iloc[:, col_idx] - test_pred

        except Exception as e:
            print(f"Feature {col_idx} regression failed: {e}")
            train_resid_data[:, col_idx] = np.nan
            if test_resid_data is not None:
                test_resid_data[:, col_idx] = np.nan
            models.append(None)

    if test_data is not None:
        return train_resid_data, test_resid_data, models
    else:
        return train_resid_data, models