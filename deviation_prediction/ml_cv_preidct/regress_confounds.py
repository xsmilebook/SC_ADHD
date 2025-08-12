import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.formula.api import ols

def regress_confounds(train_data, train_confounds,
                       test_data=None, test_confounds=None,
                       confound_names=None, categorical_confounds=None):
    """
    regress_confounds

    :parameters
    - train_data: shape (n_train_samples, n_features)
    - train_confounds: shape (n_train_samples, n_confounds)
    - test_data: shape (n_test_samples, n_features), optional
    - test_confounds: shape (n_test_samples, n_confounds), optional
    - confound_names: list of str
    - categorical_confounds: list of str

    return:
    - train_resid_data: shape (n_train_samples, n_features)
    - test_resid_data: shape (n_test_samples, n_features) if test_data is provided
    - models: list of fitted OLS models (one per feature)
    """

    if isinstance(train_confounds, np.ndarray):
        train_cov_df = pd.DataFrame(train_confounds)
        if confound_names is None:
            confound_names = [f'cov_{i}' for i in range(train_confounds.shape[1])]
        train_cov_df.columns = confound_names
    else:
        train_cov_df = train_confounds.copy()
        confound_names = train_cov_df.columns.tolist()

    if isinstance(train_data, np.ndarray):
        train_data_df = pd.DataFrame(train_data, columns=[f'feature_{i}' for i in range(train_data.shape[1])])
    else:
        train_data_df = train_data.copy()

    if test_data is not None and test_confounds is not None:
        if isinstance(test_confounds, np.ndarray):
            test_cov_df = pd.DataFrame(test_confounds, columns=confound_names)
        else:
            test_cov_df = test_confounds.copy()

        if isinstance(test_data, np.ndarray):
            test_data_df = pd.DataFrame(test_data, columns=[f'feature_{i}' for i in range(test_data.shape[1])])
        else:
            test_data_df = test_data.copy()
    else:
        test_cov_df = None
        test_data_df = None

    formula_parts = []
    for name in confound_names:
        if categorical_confounds and name in categorical_confounds:
            formula_parts.append(f'C({name})')
        else:
            formula_parts.append(name)
    formula_base = ' + '.join(formula_parts)
    full_formula = f'data_col ~ {formula_base}'

    train_resid_df = pd.DataFrame(index=train_data_df.index, columns=train_data_df.columns)
    test_resid_df = None
    if test_data_df is not None:
        test_resid_df = pd.DataFrame(index=test_data_df.index, columns=test_data_df.columns)

    # train_resid_data = np.zeros_like(train_data, dtype=float)
    # test_resid_data = np.zeros_like(test_data, dtype=float) if test_data is not None else None

    models = []

    for col_name in train_data_df.columns:
        df_temp_train = train_cov_df.copy()
        df_temp_train['data_col'] = train_data_df[col_name]

        try:
            model = ols(full_formula, data=df_temp_train).fit()
            models.append(model)

            train_resid_df[col_name] = model.resid

            if test_data_df is not None:
                test_pred = model.predict(test_cov_df)
                test_resid_df[col_name] = test_data_df[col_name] - test_pred

        except Exception as e:
            print(f"特征 '{col_name}' 回归失败: {e}")
            train_resid_df[col_name] = np.nan
            if test_resid_df is not None:
                test_resid_df[col_name] = np.nan
            models.append(None)

    if test_data is not None:
        return train_resid_df, test_resid_df, models
    else:
        return train_resid_df, models