import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.svm import SVR
from sklearn.model_selection import cross_val_score, KFold, LeaveOneOut, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import pearsonr, spearmanr
import warnings

warnings.filterwarnings('ignore')

# 读取数据
data_df = pd.read_csv(r'/datasets/deviation/PKU6/SCdata_deviationZ_medication_PKU6.csv')
columns = ["SC.109_deviationZ", "SC.120_deviationZ", "SC.59_deviationZ", "SC.65_deviationZ", "SC.79_deviationZ",
           "SC.84_deviationZ"]
predict_label = ["IA_delta", "HI_delta"]

data_clean = data_df[columns + predict_label].dropna()

X = data_clean[columns].values
y_IA = data_clean[predict_label[0]].values
y_HI = data_clean[predict_label[1]].values

def global_standardize(data):
    global_mean = data.mean()
    global_std = data.std()
    return (data - global_mean) / global_std

X_scaled = global_standardize(X)



# 超参数网格
param_grid = {
    'C': [0.1, 1, 10, 100],
    'gamma': ['scale', 'auto', 0.001, 0.01, 0.1, 1]
}

# 内部交叉验证调参 + 外部评估
def perform_nested_cv(X, y, target_name):
    print(f"\n=== Nested Cross-Validation for {target_name} ===")

    outer_cv = KFold(n_splits=10, shuffle=True, random_state=42)
    inner_cv = KFold(n_splits=5, shuffle=True, random_state=42)

    mse_scores = []
    r2_scores = []

    for train_idx, test_idx in outer_cv.split(X):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        # 内部 GridSearchCV
        grid_search = GridSearchCV(
            estimator=SVR(kernel='rbf'),
            param_grid=param_grid,
            cv=inner_cv,
            scoring='neg_mean_squared_error',
            n_jobs=-1
        )
        grid_search.fit(X_train, y_train)

        best_model = grid_search.best_estimator_

        y_pred = best_model.predict(X_test)
        mse = mean_squared_error(y_test, y_pred)
        r2 = r2_score(y_test, y_pred)

        mse_scores.append(mse)
        r2_scores.append(r2)

    mse_scores = np.array(mse_scores)
    r2_scores = np.array(r2_scores)

    print(f"MSE scores: {mse_scores}")
    print(f"R² scores: {r2_scores}")
    print(f"Average MSE: {mse_scores.mean():.4f} (+/- {mse_scores.std() * 2:.4f})")
    print(f"Average R²: {r2_scores.mean():.4f} (+/- {r2_scores.std() * 2:.4f})")

    return mse_scores, r2_scores


def perform_loocv(X, y, target_name):
    print(f"\n=== LOOCV for {target_name} ===")

    loo = LeaveOneOut()
    y_true_loo = []
    y_pred_loo = []

    inner_cv = KFold(n_splits=5, shuffle=True, random_state=42)

    for train_index, test_index in loo.split(X):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]

        # 内部 GridSearchCV
        grid_search = GridSearchCV(
            estimator=SVR(kernel='rbf'),
            param_grid=param_grid,
            cv=inner_cv,
            scoring='neg_mean_squared_error',
            n_jobs=-1
        )
        grid_search.fit(X_train, y_train)

        best_model = grid_search.best_estimator_
        y_pred = best_model.predict(X_test)

        y_true_loo.append(y_test[0])
        y_pred_loo.append(y_pred[0])

    y_true_loo = np.array(y_true_loo)
    y_pred_loo = np.array(y_pred_loo)

    mse = mean_squared_error(y_true_loo, y_pred_loo)
    r2 = r2_score(y_true_loo, y_pred_loo)

    pearson_corr, pearson_p = pearsonr(y_true_loo, y_pred_loo)
    spearman_corr, spearman_p = spearmanr(y_true_loo, y_pred_loo)

    print(f"MSE: {mse:.4f}")
    print(f"R²: {r2:.4f}")
    print(f"Pearson correlation: {pearson_corr:.4f} (p-value: {pearson_p:.4f})")
    print(f"Spearman correlation: {spearman_corr:.4f} (p-value: {spearman_p:.4f})")

    plt.figure(figsize=(8, 6))
    plt.scatter(y_true_loo, y_pred_loo, alpha=0.6)
    plt.plot([y_true_loo.min(), y_true_loo.max()], [y_true_loo.min(), y_true_loo.max()], 'r--', lw=2)
    plt.xlabel('Actual Values')
    plt.ylabel('Predicted Values')
    plt.title(f'LOOCV: Actual vs Predicted for {target_name}\n'
              f'Pearson r = {pearson_corr:.3f}, Spearman ρ = {spearman_corr:.3f}')
    plt.grid(True, alpha=0.3)
    plt.show()

    print(f"\n--- Correlation Analysis Report for {target_name} ---")
    print(f"Number of samples: {len(y_true_loo)}")
    print(f"Pearson correlation coefficient (r): {pearson_corr:.4f}")
    print(f"Pearson p-value: {pearson_p:.4f}")

    if pearson_p < 0.05:
        print("Pearson correlation is statistically significant (p < 0.05)")
    else:
        print("Pearson correlation is NOT statistically significant (p >= 0.05)")

    print(f"\nSpearman correlation coefficient (ρ): {spearman_corr:.4f}")
    print(f"Spearman p-value: {spearman_p:.4f}")

    if spearman_p < 0.05:
        print("Spearman correlation is statistically significant (p < 0.05)")
    else:
        print("Spearman correlation is NOT statistically significant (p >= 0.05)")

    return mse, r2, y_true_loo, y_pred_loo


results = {}

for i, target in enumerate(predict_label):
    y_target = data_clean[target].values

    print(f"\n{'=' * 50}")
    print(f"Predicting {target}")
    print(f"{'=' * 50}")

    # 嵌套交叉验证
    mse_nested, r2_nested = perform_nested_cv(X_scaled, y_target, target)

    # LOOCV
    mse_loocv, r2_loocv, y_true, y_pred = perform_loocv(X_scaled, y_target, target)

    results[target] = {
        'nested_cv': {'mse': mse_nested, 'r2': r2_nested},
        'loocv': {'mse': mse_loocv, 'r2': r2_loocv, 'y_true': y_true, 'y_pred': y_pred}
    }

print(f"\n{'=' * 60}")
print("SUMMARY RESULTS")
print(f"{'=' * 60}")

for target in predict_label:
    print(f"\n{target}:")
    print(
        f"  Nested CV - MSE: {results[target]['nested_cv']['mse'].mean():.4f}, R²: {results[target]['nested_cv']['r2'].mean():.4f}")
    print(f"  LOOCV     - MSE: {results[target]['loocv']['mse']:.4f}, R²: {results[target]['loocv']['r2']:.4f}")


# 特征重要性分析（使用最优参数）
def feature_importance_analysis(X, y, feature_names, target_name):
    print(f"\n=== Feature Importance for {target_name} ===")

    # 基线模型（使用最优参数）
    inner_cv = KFold(n_splits=5, shuffle=True, random_state=42)
    grid_search = GridSearchCV(
        estimator=SVR(kernel='rbf'),
        param_grid=param_grid,
        cv=inner_cv,
        scoring='r2',
        n_jobs=-1
    )
    grid_search.fit(X, y)
    best_model = grid_search.best_estimator_

    baseline_r2 = cross_val_score(best_model, X, y, cv=KFold(n_splits=10, shuffle=True, random_state=42), scoring='r2').mean()
    print(f"Baseline R² (with best params): {baseline_r2:.4f}")

    feature_importance = []
    for i, feature_name in enumerate(feature_names):
        X_permuted = X.copy()
        np.random.shuffle(X_permuted[:, i])

        permuted_r2 = cross_val_score(best_model, X_permuted, y, cv=KFold(n_splits=10, shuffle=True, random_state=42), scoring='r2').mean()

        importance = baseline_r2 - permuted_r2
        feature_importance.append(importance)
        print(f"{feature_name}: {importance:.4f}")

    plt.figure(figsize=(10, 6))
    indices = np.argsort(feature_importance)[::-1]
    plt.bar(range(len(feature_importance)), np.array(feature_importance)[indices])
    plt.xticks(range(len(feature_importance)), np.array(feature_names)[indices], rotation=45, ha='right')
    plt.xlabel('Features')
    plt.ylabel('Importance (R² decrease)')
    plt.title(f'Feature Importance for {target_name}')
    plt.tight_layout()
    plt.show()

    return feature_importance


for i, target in enumerate(predict_label):
    y_target = data_clean[target].values
    importance = feature_importance_analysis(X_scaled, y_target, columns, target)