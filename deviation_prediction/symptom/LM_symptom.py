import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import cross_val_score, KFold, LeaveOneOut
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

# 使用sklearn的标准标准化方法
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)


# 调试函数：检查预测是否颠倒
def debug_predictions(y_true, y_pred, method_name):
    """调试预测结果"""
    print(f"\n--- {method_name} 调试信息 ---")
    print(f"真实值范围: [{y_true.min():.3f}, {y_true.max():.3f}]")
    print(f"预测值范围: [{y_pred.min():.3f}, {y_pred.max():.3f}]")
    print(f"真实值均值: {y_true.mean():.3f}")
    print(f"预测值均值: {y_pred.mean():.3f}")

    # 检查是否系统性负相关
    correlation_check = np.corrcoef(y_true, y_pred)[0, 1]
    print(f"相关系数检查: {correlation_check:.4f}")

    # 简单基线：总是预测均值
    baseline_pred = np.full_like(y_true, y_true.mean())
    baseline_r2 = r2_score(y_true, baseline_pred)
    print(f"基线R²(总是预测均值): {baseline_r2:.4f}")


def perform_10fold_cv(X, y, target_name):
    print(f"\n=== 10-Fold CV for {target_name} ===")

    kf = KFold(n_splits=10, shuffle=True, random_state=42)

    y_true_all = []
    y_pred_all = []

    # 手动进行交叉验证
    for train_index, test_index in kf.split(X):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]

        lr_model = LinearRegression()
        lr_model.fit(X_train, y_train)
        y_pred = lr_model.predict(X_test)

        y_true_all.extend(y_test)
        y_pred_all.extend(y_pred)

    y_true_all = np.array(y_true_all)
    y_pred_all = np.array(y_pred_all)

    # 调试检查
    debug_predictions(y_true_all, y_pred_all, "10-Fold CV")

    # 计算指标（确保正确顺序）
    mse = mean_squared_error(y_true_all, y_pred_all)
    r2 = r2_score(y_true_all, y_pred_all)
    pearson_corr, pearson_p = pearsonr(y_true_all, y_pred_all)  # 正确：真实值, 预测值
    spearman_corr, spearman_p = spearmanr(y_true_all, y_pred_all)

    print(f"\n结果指标:")
    print(f"MSE: {mse:.4f}")
    print(f"R²: {r2:.4f}")
    print(f"Pearson correlation (r): {pearson_corr:.4f} (p-value: {pearson_p:.4f})")
    print(f"Spearman correlation (ρ): {spearman_corr:.4f} (p-value: {spearman_p:.4f})")

    # 绘制预测结果
    plt.figure(figsize=(8, 6))
    plt.scatter(y_true_all, y_pred_all, alpha=0.6)
    plt.plot([y_true_all.min(), y_true_all.max()], [y_true_all.min(), y_true_all.max()], 'r--', lw=2)
    plt.xlabel('Actual Values')
    plt.ylabel('Predicted Values')
    plt.title(f'10-Fold CV: Actual vs Predicted for {target_name}\n'
              f'r = {pearson_corr:.3f}, R² = {r2:.3f}')
    plt.grid(True, alpha=0.3)
    plt.show()

    return mse, r2, pearson_corr, spearman_corr, y_true_all, y_pred_all


def perform_loocv(X, y, target_name):
    print(f"\n=== LOOCV for {target_name} ===")

    loo = LeaveOneOut()
    y_true_loo = []
    y_pred_loo = []

    for train_index, test_index in loo.split(X):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]

        lr_model = LinearRegression()
        lr_model.fit(X_train, y_train)
        y_pred = lr_model.predict(X_test)

        y_true_loo.append(y_test[0])
        y_pred_loo.append(y_pred[0])

    y_true_loo = np.array(y_true_loo)
    y_pred_loo = np.array(y_pred_loo)

    # 调试检查
    debug_predictions(y_true_loo, y_pred_loo, "LOOCV")

    # 计算指标
    mse = mean_squared_error(y_true_loo, y_pred_loo)
    r2 = r2_score(y_true_loo, y_pred_loo)
    pearson_corr, pearson_p = pearsonr(y_true_loo, y_pred_loo)  # 正确顺序
    spearman_corr, spearman_p = spearmanr(y_true_loo, y_pred_loo)

    print(f"\n结果指标:")
    print(f"MSE: {mse:.4f}")
    print(f"R²: {r2:.4f}")
    print(f"Pearson correlation (r): {pearson_corr:.4f} (p-value: {pearson_p:.4f})")
    print(f"Spearman correlation (ρ): {spearman_corr:.4f} (p-value: {spearman_p:.4f})")

    plt.figure(figsize=(8, 6))
    plt.scatter(y_true_loo, y_pred_loo, alpha=0.6)
    plt.plot([y_true_loo.min(), y_true_loo.max()], [y_true_loo.min(), y_true_loo.max()], 'r--', lw=2)
    plt.xlabel('Actual Values')
    plt.ylabel('Predicted Values')
    plt.title(f'LOOCV: Actual vs Predicted for {target_name}\n'
              f'r = {pearson_corr:.3f}, R² = {r2:.3f}')
    plt.grid(True, alpha=0.3)
    plt.show()

    return mse, r2, pearson_corr, spearman_corr, y_true_loo, y_pred_loo


# 数据探索函数
def explore_basic_statistics(data_clean, columns, predict_label):
    print("=== 基本统计信息 ===")

    for target in predict_label:
        y = data_clean[target].values
        print(f"\n{target}:")
        print(f"  样本数: {len(y)}")
        print(f"  均值: {y.mean():.4f}")
        print(f"  标准差: {y.std():.4f}")
        print(f"  范围: [{y.min():.4f}, {y.max():.4f}]")

    print(f"\n特征统计:")
    for col in columns:
        x = data_clean[col].values
        print(f"{col}: 均值={x.mean():.4f}, std={x.std():.4f}")


# 主程序执行
print("开始数据分析...")

# 基本统计信息
explore_basic_statistics(data_clean, columns, predict_label)

results = {}

for i, target in enumerate(predict_label):
    y_target = data_clean[target].values

    print(f"\n{'=' * 50}")
    print(f"Predicting {target}")
    print(f"{'=' * 50}")

    # 10-Fold CV
    try:
        mse_10fold, r2_10fold, r_10fold, rho_10fold, y_true_10fold, y_pred_10fold = perform_10fold_cv(X_scaled,
                                                                                                      y_target, target)
    except Exception as e:
        print(f"10-Fold CV 出错: {e}")
        continue

    # LOOCV
    try:
        mse_loocv, r2_loocv, r_loocv, rho_loocv, y_true_loocv, y_pred_loocv = perform_loocv(X_scaled, y_target, target)
    except Exception as e:
        print(f"LOOCV 出错: {e}")
        continue

    results[target] = {
        '10fold': {'mse': mse_10fold, 'r2': r2_10fold, 'r': r_10fold, 'rho': rho_10fold,
                   'y_true': y_true_10fold, 'y_pred': y_pred_10fold},
        'loocv': {'mse': mse_loocv, 'r2': r2_loocv, 'r': r_loocv, 'rho': rho_loocv,
                  'y_true': y_true_loocv, 'y_pred': y_pred_loocv}
    }

print(f"\n{'=' * 60}")
print("SUMMARY RESULTS")
print(f"{'=' * 60}")

for target in predict_label:
    if target in results:
        print(f"\n{target}:")
        print(f"  10-Fold CV - MSE: {results[target]['10fold']['mse']:.4f}, "
              f"R²: {results[target]['10fold']['r2']:.4f}, "
              f"r: {results[target]['10fold']['r']:.4f}")
        print(f"  LOOCV      - MSE: {results[target]['loocv']['mse']:.4f}, "
              f"R²: {results[target]['loocv']['r2']:.4f}, "
              f"r: {results[target]['loocv']['r']:.4f}")


# 添加特征-目标相关性分析
def analyze_feature_target_correlations(data_clean, columns, predict_label):
    print("=== 特征与目标变量相关性分析 ===")

    for target in predict_label:
        print(f"\n{target}:")
        y = data_clean[target].values

        print("特征相关性:")
        correlations = []
        for col in columns:
            x = data_clean[col].values
            corr, p_val = pearsonr(x, y)
            correlations.append(corr)
            significance = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else ""
            print(f"  {col}: r={corr:.3f} (p={p_val:.3f}) {significance}")

        # 计算平均相关性
        mean_abs_corr = np.mean(np.abs(correlations))
        max_corr = max(np.abs(correlations))
        print(f"  平均|相关性|: {mean_abs_corr:.3f}")
        print(f"  最大|相关性|: {max_corr:.3f}")

        if mean_abs_corr < 0.1:
            print("  ⚠️  警告: 特征与目标变量相关性极弱!")


# 在主程序中调用
analyze_feature_target_correlations(data_clean, columns, predict_label)