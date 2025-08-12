import numpy as np
import pandas as pd
from sklearn.model_selection import KFold, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from sklearn.feature_selection import SelectKBest, f_regression
from sklearn.pipeline import Pipeline
from scipy.stats import pearsonr
from regress_confounds import regress_confounds
import warnings

warnings.filterwarnings("ignore")
def run_svr_nested_cv(
    data_df,
    feature_cols,
    target_col,
    outer_folds=10,
    inner_folds=5,
    n_runs=1,
    random_state_base=42,
    confounds_df=None,
    permutation_flag=False,
    verbose=True
):
    """
    nested_cv for svr：外层评估性能，内层进行特征选择 + 超参数调优。

    Parameters:
    -----------
    data_df : pd.DataFrame
        包含特征和目标的数据
    feature_cols : list of str
        所有候选特征列名
    target_col : str
        目标变量
    outer_folds : int
        外层CV折数（性能评估）
    inner_folds : int
        内层CV折数（调参）
    n_runs : int
        重复嵌套CV的次数（用于稳定结果）
    random_state_base : int
        随机种子基础
    permutation_flag : bool
        是否使用打乱的标签（用于置换检验）
    verbose : bool

    Returns:
    --------
    dict : 包含性能、参数、特征等信息
    """

    valid_df = data_df[feature_cols + [target_col]].dropna()
    X = valid_df[feature_cols]
    y = valid_df[target_col]

    if len(X) == 0:
        raise ValueError(f"No valid data for {target_col}")

    all_outer_corrs = []
    best_params_list = []
    selected_k_list = []
    selected_features_list = []

    for run in range(n_runs):
        rs = random_state_base + run
        if verbose:
            mode = "Permuted" if permutation_flag else "Original"
            print(f"[{mode}] {target_col} - Run {run+1}/{n_runs}")

        # outer CV loop
        outer_kf = KFold(n_splits=outer_folds, shuffle=True, random_state=rs)
        fold_corrs = []

        for fold, (train_idx, test_idx) in enumerate(outer_kf.split(X)):
            X_train_outer, X_test_outer = X.iloc[train_idx], X.iloc[test_idx]
            y_train_outer, y_test_outer = y.iloc[train_idx], y.iloc[test_idx]

            confounds_train = None
            confounds_test = None
            confounds_var_list = []
            category_confounds_var_list = []
            if(confounds_df is not None):
                confounds_var_list = confounds_df.columns
                category_confounds_var_list = confounds_df.select_dtypes(include=['category']).columns.tolist()
                confounds_train = confounds_df.iloc[train_idx]
                confounds_test = confounds_df.iloc[test_idx]

            X_train_outer, X_test_outer, reg_resid_models = regress_confounds(
                train_data=X_train_outer,
                train_confounds=confounds_train,
                test_data=X_test_outer,
                test_confounds=confounds_test,
                confound_names=confounds_var_list,
                categorical_confounds=category_confounds_var_list
            )

            if permutation_flag:
                y_train_outer = y_train_outer.sample(frac=1, random_state=rs + fold).reset_index(drop=True)

            # -----------------------------
            # Inner CV loop：超参数调优 + 特征选择
            # -----------------------------
            pipe = Pipeline([
                ('scaler', StandardScaler()),
                ('selector', SelectKBest(f_regression)),
                ('svr', SVR())
            ])

            param_grid = {
                'selector__k': list(range(1, len(feature_cols) + 1)),  # k=1~10
                'svr__C': [0.1, 0.6, 1],
                'svr__gamma': ['scale', 'auto'],
                'svr__kernel': ['rbf']
            }

            inner_kf = KFold(n_splits=inner_folds, shuffle=True, random_state=rs)
            grid = GridSearchCV(
                pipe,
                param_grid,
                cv=inner_kf,
                scoring='r2',
                n_jobs=-1,
                verbose=0
            )
            grid.fit(X_train_outer, y_train_outer)

            best_model = grid.best_estimator_
            best_params_list.append(grid.best_params_)
            selected_k_list.append(grid.best_params_['selector__k'])

            # get selected features
            mask = best_model.named_steps['selector'].get_support()
            selected_features = X_train_outer.columns[mask].tolist()
            selected_features_list.append(selected_features)

            # Outer CV Loop
            y_pred_outer = best_model.predict(X_test_outer)
            corr, _ = pearsonr(y_test_outer, y_pred_outer)
            fold_corrs.append(corr)

        run_mean_corr = np.mean(fold_corrs)
        all_outer_corrs.append(run_mean_corr)

        if verbose and n_runs > 1:
            print(f"  Run {run+1} mean corr: {run_mean_corr:.4f}")

    # 汇总结果
    mean_correlation = float(np.mean(all_outer_corrs))
    std_correlation = float(np.std(all_outer_corrs))

    result = {
        'target': target_col,
        'mean_correlation': mean_correlation,
        'std_correlation': std_correlation,
        'all_correlations': all_outer_corrs,
        'permutation_flag': permutation_flag,
        'n_runs': n_runs,
        'outer_folds': outer_folds,
        'inner_folds': inner_folds,
        'best_params': best_params_list,
        'selected_k': selected_k_list,
        'selected_features': selected_features_list
    }

    if verbose:
        print(f"[{target_col}] Nested CV Done. Mean corr: {mean_correlation:.4f} ± {std_correlation:.4f}")

    return result