import os
import pandas as pd
from sklearn.inspection import permutation_importance
from regress_covariates import regress_covariates
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LinearRegression
from sklearn.svm import SVR
from scipy.stats import pearsonr
import numpy as np
import matplotlib.pyplot as plt

data_df = pd.read_csv(
    'D:/code/SC_ADHD/datasets/ABCD/task-fMRI/deviation/SCdeviation_nbackActivation_deviation_ADHD_TDtest.csv',
    low_memory=False)

# 1. filter the ADHD subset
data_df = data_df[data_df["if_TD"] == "ADHD"].reset_index(drop=True)
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
data_df = data_df[feature_columns + predict_columns + covariates_list]

# 2. regress confounds
data_df['sex'] = data_df['sex'].astype('category')
covariates_list = ["age", "sex"]
covariates_data = data_df[covariates_list]
category_cov_list = ["sex"]

data_df[predict_columns] = regress_covariates(data_df[predict_columns].copy(), covariates_data, covariates_list,
                                              category_cov_list)

# 3. select features
temp_df = data_df[feature_columns + predict_columns]
correlation_matrix = temp_df.corr()
selected_feature_col_dict = {}
for col in predict_columns:
    correlations = correlation_matrix.loc[feature_columns, col].abs()
    top_5 = correlations.nlargest(5)
    selected_feature_col_dict[col] = top_5.index.tolist()

model_name = "LR"
# model = SVR(kernel='rbf', C=0.6)
model = LinearRegression()

original_correlations = {col: [] for col in predict_columns}
permuted_correlations = {col: [] for col in predict_columns}

# 101 Times, 10-fold cross-validation
print("Running 101 runs of 10-fold cross-validation for original correlations...")
n_original_runs = 101

all_original_correlations = {col: [] for col in predict_columns}

for run in range(n_original_runs):
    if run % 10 == 0:
        print(f"Original run {run}/{n_original_runs}")

    kf = KFold(n_splits=10, shuffle=True, random_state=run)
    splits = []

    for train_index, test_index in kf.split(data_df):
        X_train = data_df.iloc[train_index][feature_columns]
        X_test = data_df.iloc[test_index][feature_columns]
        y_train = data_df.iloc[train_index][predict_columns]
        y_test = data_df.iloc[test_index][predict_columns]
        splits.append((X_train, X_test, y_train, y_test))

    run_correlations = {col: [] for col in predict_columns}

    for predict_column in predict_columns:
        correlations = []
        selected_features = selected_feature_col_dict[predict_column]
        for X_train, X_test, y_train, y_test in splits:
            X_train_sel = X_train[selected_features]
            X_test_sel = X_test[selected_features]

            scaler = StandardScaler()
            X_train_scaled = scaler.fit_transform(X_train_sel)
            X_test_scaled = scaler.transform(X_test_sel)

            model.fit(X_train_scaled, y_train[predict_column])
            y_pred = model.predict(X_test_scaled)
            corr = np.corrcoef(y_test[predict_column], y_pred)[0, 1]
            correlations.append(corr)
        run_correlations[predict_column] = np.mean(correlations)

    for col in predict_columns:
        all_original_correlations[col].append(run_correlations[col])

# 计算101次运行的平均相关性作为原始相关性
final_original_correlations = {}
for col in predict_columns:
    final_original_correlations[col] = np.mean(all_original_correlations[col])
    original_correlations[col] = all_original_correlations[col]  # 保存所有101次的结果用于后续分析

print("Final original correlations (101 runs average):")
for col in predict_columns:
    print(f"  {col}: {final_original_correlations[col]:.4f}")

# 置换检验 1000 次
n_permutations = 1000
print(f"\nRunning {n_permutations} permutation tests...")

for perm_run in range(n_permutations):
    if perm_run % 100 == 0:
        print(f"Permutation run {perm_run}/{n_permutations}")

    # 使用固定的10-fold分割（与原始运行保持一致，使用相同的random_state）
    kf_perm = KFold(n_splits=10, shuffle=True, random_state=perm_run % 101)  # 循环使用原始的random_state
    perm_splits = []

    for train_index, test_index in kf_perm.split(data_df):
        X_train = data_df.iloc[train_index][feature_columns]
        X_test = data_df.iloc[test_index][feature_columns]
        y_train = data_df.iloc[train_index][predict_columns]
        y_test = data_df.iloc[test_index][predict_columns]

        # 打乱训练集的标签（即 y_train）
        y_train_shuffled = y_train.sample(frac=1, random_state=perm_run).reset_index(drop=True)
        perm_splits.append((X_train, X_test, y_train_shuffled, y_test))

    for predict_column in predict_columns:
        correlations = []
        selected_features = selected_feature_col_dict[predict_column]
        for X_train, X_test, y_train, y_test in perm_splits:
            X_train_sel = X_train[selected_features]
            X_test_sel = X_test[selected_features]

            scaler = StandardScaler()
            X_train_scaled = scaler.fit_transform(X_train_sel)
            X_test_scaled = scaler.transform(X_test_sel)

            model.fit(X_train_scaled, y_train[predict_column])
            y_pred = model.predict(X_test_scaled)
            corr = np.corrcoef(y_test[predict_column], y_pred)[0, 1]
            correlations.append(corr)
        permuted_correlations[predict_column].append(np.mean(correlations))

# 计算 p 值：有多少次置换结果 >= 原始结果（使用101次的平均值）
p_values = {}
for col in predict_columns:
    orig = final_original_correlations[col]
    perm = np.array(permuted_correlations[col])
    p_value = np.mean(perm >= orig)
    p_values[col] = p_value

print("\nPermutation test results:")
for col in predict_columns:
    orig = final_original_correlations[col]
    p_val = p_values[col]
    significance = "**" if p_val < 0.05 else ""
    print(f"  {col}: original corr = {orig:.4f}, p-value = {p_val:.4f} {significance}")

# 可视化分布
output_dir = os.path.join("D:/code/SC_ADHD/datasets/prediction/activation_deviation", model_name)
os.makedirs(output_dir, exist_ok=True)

for col in predict_columns:
    plt.figure(figsize=(8, 5))
    plt.hist(permuted_correlations[col], bins=50, alpha=0.7, color='gray', label='Permuted correlations')
    plt.axvline(final_original_correlations[col], color='red', linestyle='--', linewidth=2,
                label='Original correlation')
    plt.title(f'{col} - Permutation Test')
    plt.xlabel('Correlation')
    plt.ylabel('Frequency')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/permutation_test_{col}.png", dpi=300)
    plt.close()

# 保存详细结果到CSV
results_df = pd.DataFrame({
    'target': predict_columns,
    'original_correlation_mean': [final_original_correlations[col] for col in predict_columns],
    'original_correlation_std': [np.std(original_correlations[col]) for col in predict_columns],
    'p_value': [p_values[col] for col in predict_columns],
    'significant': [p_values[col] < 0.05 for col in predict_columns]
})

results_df.to_csv(f"{output_dir}/permutation_test_results.csv", index=False)
print(f"\nResults saved to {output_dir}/permutation_test_results.csv")