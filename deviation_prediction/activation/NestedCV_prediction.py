import os
import pandas as pd
from regress_covariates import regress_covariates
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold, GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LinearRegression
from sklearn.svm import SVR
from scipy.stats import pearsonr
import numpy as np

df = pd.read_csv('D:/code/SC_ADHD/datasets/ABCD/task-fMRI/deviation/SCdeviation_nbackActivation_deviation_ADHD_TDtest.csv', low_memory=False)

# 1. filter the ADHD subset
df = df[df["if_TD"] == "ADHD"].reset_index(drop=True)
feature_columns = [
    "SC.109_deviationZ", "SC.111_deviationZ", "SC.116_deviationZ", "SC.118_deviationZ",
    "SC.119_deviationZ", "SC.120_deviationZ", "SC.18_deviationZ", "SC.57_deviationZ",
    "SC.59_deviationZ", "SC.60_deviationZ", "SC.65_deviationZ", "SC.70_deviationZ",
    "SC.79_deviationZ", "SC.84_deviationZ", "SC.8_deviationZ", "SC.93_deviationZ",
    "SC.98_deviationZ"
]
predict_columns = [
    "nbackFP.A_deviationZ",
    "nbackFP.B_deviationZ",
    "nbackSM.A_deviationZ",
    "nbackSM.B_deviationZ",
    "nbackVA.B_deviationZ"
]
covariates_list = ["age", "sex", "if_TD"]
df = df[feature_columns + predict_columns + covariates_list]

# 2. regress confounds
df['sex'] = df['sex'].astype('category')
covariates_list = ["age", "sex"]
covariates_data = df[covariates_list]
category_cov_list = ["sex"]

df[predict_columns] = regress_covariates(df[predict_columns].copy(), covariates_data, covariates_list, category_cov_list)

# 3. select features
temp_df = df[feature_columns + predict_columns]
correlation_matrix = temp_df.corr()
selected_feature_col_dict = {}
for col in predict_columns:
    correlations = correlation_matrix.loc[feature_columns, col].abs()

    top_5 = correlations.nlargest(5)
    # print(top_5)
    selected_feature_col_dict[col] = top_5.index.tolist()

# 4. normalize
# df_features = df[feature_columns]
#
# # Initialize and fit the scaler
# scaler = StandardScaler()
# scaler.fit(df_features)
#
# # Transform the features
# df_scaled_features = pd.DataFrame(scaler.transform(df_features), columns=feature_columns)
#
# # Combine scaled features with target variables
# df_prepared = pd.concat([df_scaled_features, df[predict_columns]], axis=1)


# 5. Nested_CV
outer_scores = {col: [] for col in predict_columns}
outer_kf = KFold(n_splits=10, shuffle=True, random_state=42)

model_name = "SVR"

param_grid = {
    'svr__C': [0.1, 0.5, 1.0, 2.0],
    'svr__gamma': ['scale', 'auto', 0.01, 0.001],
    'svr__kernel': ['linear', 'rbf']
}

# outer CV
for fold_idx, (train_idx, test_idx) in enumerate(outer_kf.split(df)):
    print(f"\nOuter Fold {fold_idx + 1}/10")

    # split datasets
    X_train_outer = df.iloc[train_idx][feature_columns]
    X_test_outer  = df.iloc[test_idx][feature_columns]
    y_train_outer = df.iloc[train_idx][predict_columns]
    y_test_outer  = df.iloc[test_idx][predict_columns]

    # inner CV(parameter optimize)
    inner_kf = KFold(n_splits=10, shuffle=True, random_state=42)

    for predict_column in predict_columns:
        print(f"  Tuning for {predict_column}...")

        # pipeline: normalize + model
        model = Pipeline([
            ('scaler', StandardScaler()),
            ('svr', SVR())
        ])

        # grid search, inner CV
        grid_search = GridSearchCV(
            estimator=model,
            param_grid=param_grid,
            cv=inner_kf,
            scoring='neg_mean_squared_error',
            n_jobs=-1
        )

        # fitting and select the best parameters
        grid_search.fit(X_train_outer, y_train_outer[predict_column])

        # test in outer test set
        y_pred = grid_search.predict(X_test_outer)

        # compute the correlation
        corr, _ = pearsonr(y_test_outer[predict_column], y_pred)
        outer_scores[predict_column].append(corr)

        print(f"    Best params: {grid_search.best_params_}")
        print(f"    Test corr: {corr:.4f}")


print("\n" + "="*50)
print("Final Nested CV Results (10×10)")
print("="*50)
for col in predict_columns:
    mean_corr = np.mean(outer_scores[col])
    std_corr  = np.std(outer_scores[col])
    print(f"{col}: {mean_corr:.4f} ± {std_corr:.4f}")

# 6. save results


