import numpy as np
import pandas as pd
from sklearn.model_selection import LeaveOneOut, KFold, cross_val_score
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, r2_score
from scipy import stats
import warnings

warnings.filterwarnings('ignore')


def robust_prediction_evaluation(X, y, feature_names, target_name, n_bootstrap=100):
    """
    稳健的预测能力评估方案
    """
    print(f"\n{'=' * 60}")
    print(f"ROBUST PREDICTION EVALUATION FOR {target_name}")
    print(f"{'=' * 60}")
    print(f"Sample size: {len(y)}")

    # 1. 使用Ridge回归（适合小数据集）
    model = Ridge(alpha=1.0)

    # 2. LOOCV评估
    print("\n1. Leave-One-Out Cross Validation (LOOCV):")
    loo = LeaveOneOut()
    y_true_loo, y_pred_loo = [], []

    for train_idx, test_idx in loo.split(X):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)

        y_true_loo.append(y_test[0])
        y_pred_loo.append(y_pred[0])

    y_true_loo = np.array(y_true_loo)
    y_pred_loo = np.array(y_pred_loo)

    loo_mse = mean_squared_error(y_true_loo, y_pred_loo)
    loo_r2 = r2_score(y_true_loo, y_pred_loo)

    # 相关性分析
    pearson_r, pearson_p = stats.pearsonr(y_true_loo, y_pred_loo)
    spearman_r, spearman_p = stats.spearmanr(y_true_loo, y_pred_loo)

    print(f"   LOOCV R²: {loo_r2:.4f}")
    print(f"   LOOCV MSE: {loo_mse:.4f}")
    print(f"   Pearson r: {pearson_r:.4f} (p={pearson_p:.4f})")
    print(f"   Spearman ρ: {spearman_r:.4f} (p={spearman_p:.4f})")

    # 3. Bootstrap评估（提供置信区间）
    print(f"\n2. Bootstrap Evaluation ({n_bootstrap} iterations):")
    bootstrap_r2 = []
    bootstrap_corr = []

    np.random.seed(42)
    for i in range(n_bootstrap):
        # 有放回采样
        indices = np.random.choice(len(X), size=len(X), replace=True)
        X_boot = X[indices]
        y_boot = y[indices]

        # 训练测试分割
        test_indices = np.setdiff1d(np.arange(len(X)), indices)
        if len(test_indices) > 0:
            X_train, X_test = X_boot, X[test_indices]
            y_train, y_test = y_boot, y[test_indices]

            model.fit(X_train, y_train)
            y_pred = model.predict(X_test)

            if len(y_test) > 1:
                r2 = r2_score(y_test, y_pred)
                corr, _ = stats.pearsonr(y_test, y_pred)
                bootstrap_r2.append(r2)
                bootstrap_corr.append(corr)

    if bootstrap_r2:
        print(f"   Bootstrap R²: {np.mean(bootstrap_r2):.4f} "
              f"[{np.percentile(bootstrap_r2, 2.5):.4f}, "
              f"{np.percentile(bootstrap_r2, 97.5):.4f}]")
        print(f"   Bootstrap Correlation: {np.mean(bootstrap_corr):.4f} "
              f"[{np.percentile(bootstrap_corr, 2.5):.4f}, "
              f"{np.percentile(bootstrap_corr, 97.5):.4f}]")

    # 4. 特征重要性分析
    print(f"\n3. Feature Importance (Ridge Coefficients):")
    model.fit(X, y)
    coefficients = model.coef_

    # 按重要性排序
    feature_importance = list(zip(feature_names, coefficients))
    feature_importance.sort(key=lambda x: abs(x[1]), reverse=True)

    for i, (feature, coef) in enumerate(feature_importance):
        print(f"   {i + 1}. {feature}: {coef:.4f}")

    # 5. 预测稳定性分析
    print(f"\n4. Prediction Stability Analysis:")
    abs_errors = np.abs(y_true_loo - y_pred_loo)
    print(f"   Mean Absolute Error: {np.mean(abs_errors):.4f}")
    print(f"   Median Absolute Error: {np.median(abs_errors):.4f}")
    print(f"   Max Absolute Error: {np.max(abs_errors):.4f}")

    # 6. 统计显著性判断
    print(f"\n5. Statistical Significance:")
    if pearson_p < 0.05:
        print(f"   ✓ Pearson correlation is statistically significant (p < 0.05)")
    else:
        print(f"   ✗ Pearson correlation is NOT statistically significant (p ≥ 0.05)")

    if loo_r2 > 0:
        print(f"   ✓ Model explains {loo_r2 * 100:.1f}% of variance")
    else:
        print(f"   ✗ Model performs worse than simple mean prediction")

    # 7. 实用性评估
    print(f"\n6. Practical Utility Assessment:")
    effect_size = abs(pearson_r)
    if effect_size >= 0.5:
        utility = "Strong"
    elif effect_size >= 0.3:
        utility = "Moderate"
    elif effect_size >= 0.1:
        utility = "Weak"
    else:
        utility = "Negligible"

    print(f"   Effect size: {utility} (r = {pearson_r:.3f})")

    if pearson_p < 0.05 and effect_size >= 0.3:
        print(f"   ✓ Prediction has practical utility")
    elif pearson_p < 0.05:
        print(f"   ⚠ Prediction is statistically significant but weak")
    else:
        print(f"   ✗ Prediction lacks practical utility")

    return {
        'loo_r2': loo_r2,
        'loo_mse': loo_mse,
        'pearson_r': pearson_r,
        'pearson_p': pearson_p,
        'spearman_r': spearman_r,
        'spearman_p': spearman_p,
        'y_true': y_true_loo,
        'y_pred': y_pred_loo,
        'coefficients': dict(zip(feature_names, coefficients))
    }


results = {}

data_df = pd.read_csv(r'/datasets/deviation/PKU6/SCdata_deviationZ_medication_PKU6.csv')
columns = ["SC.109_deviationZ", "SC.120_deviationZ", "SC.59_deviationZ", "SC.65_deviationZ", "SC.79_deviationZ",
           "SC.84_deviationZ"]
# columns = [item for item in data_df.columns if item.endswith('deviationZ')]
predict_label = ["IA_delta", "HI_delta"]

data_clean = data_df[columns + predict_label].dropna()
X = data_clean[columns].values
y_IA = data_clean[predict_label[0]].values
y_HI = data_clean[predict_label[1]].values


for target in predict_label:
    y_target = data_clean[target].values

    # 标准化（小数据集建议标准化）
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    result = robust_prediction_evaluation(
        X_scaled, y_target, columns, target, n_bootstrap=100
    )

    results[target] = result

# 汇总报告
print(f"\n{'=' * 60}")
print("FINAL SUMMARY REPORT")
print(f"{'=' * 60}")

for target in predict_label:
    result = results[target]
    print(f"\n{target}:")
    print(f"  Prediction Performance:")
    print(f"    R² = {result['loo_r2']:.4f}")
    print(f"    Correlation = {result['pearson_r']:.4f} (p = {result['pearson_p']:.4f})")

    # 判断预测能力
    if result['pearson_p'] < 0.05 and abs(result['pearson_r']) >= 0.3:
        print(f"  ✓ Prediction capability: STRONG")
    elif result['pearson_p'] < 0.05 and abs(result['pearson_r']) >= 0.1:
        print(f"  ⚠ Prediction capability: WEAK but significant")
    else:
        print(f"  ✗ Prediction capability: INSUFFICIENT")