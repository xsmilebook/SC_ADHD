import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge, ElasticNet, RidgeCV, ElasticNetCV
from sklearn.model_selection import KFold, LeaveOneOut, cross_val_predict
from scipy.stats import pearsonr
import warnings
warnings.filterwarnings("ignore")

# -------------------------------
# 1. 加载数据
# -------------------------------
file_path = 'D:/code/SC_ADHD/datasets/ABCD/task-fMRI/deviation/SCdeviation_nbackActivation_deviation_ADHD_TDtest.csv'
data_df = pd.read_csv(file_path, low_memory=False)

# 特征列和预测列
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

# 提取特征 X 和目标 Y
X = data_df[feature_columns].values
Y = data_df[predict_columns].values

# 检查缺失值并删除含缺失的样本（或可改为插补）
print(f"原始样本数: {len(data_df)}")
print(f"缺失值统计（X）: {np.isnan(X).sum(axis=0)}")
print(f"缺失值统计（Y）: {np.isnan(Y).sum(axis=0)}")

# 删除任意位置有 NaN 的行
mask = ~(np.isnan(X).any(axis=1) | np.isnan(Y).any(axis=1))
X = X[mask]
Y = Y[mask]
print(f"清理后样本数: {X.shape[0]}")

# 标准化
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# -------------------------------
# 2. 定义交叉验证策略
# -------------------------------
kfold = KFold(n_splits=5, shuffle=True, random_state=42)
loo = LeaveOneOut()

# -------------------------------
# 3. 定义模型（自动调参）
# -------------------------------
# Ridge with built-in CV for alpha
ridge_model = RidgeCV(alphas=np.logspace(-3, 5, 50), cv=5)

# ElasticNet with built-in CV
elastic_model = ElasticNetCV(alphas=np.logspace(-5, 1, 50),
                             l1_ratio=[0.1, 0.5, 0.7, 0.9, 0.95, 0.99, 1.0],
                             cv=5, max_iter=1000, random_state=42)

models = {
    'Ridge': ridge_model,
    'ElasticNet': elastic_model
}

# 存储结果
results = []

# -------------------------------
# 4. 对每个模型和每个目标变量进行评估
# -------------------------------
for model_name, model in models.items():
    print(f"\n=== 正在运行 {model_name} ===")

    for i, target_name in enumerate(predict_columns):
        y = Y[:, i]

        # ---- K-Fold CV (5-fold) ----
        try:
            y_pred_kfold = cross_val_predict(model, X_scaled, y, cv=kfold)
            r_kfold, p_kfold = pearsonr(y, y_pred_kfold)
        except Exception as e:
            r_kfold = np.nan
            print(f"{model_name} - {target_name} k-fold CV 错误: {e}")

        # ---- LOOCV ----
        # 注意：LOOCV 计算较慢，但样本量2000尚可接受
        try:
            y_pred_loo = cross_val_predict(model, X_scaled, y, cv=loo)
            r_loo, p_loo = pearsonr(y, y_pred_loo)
        except Exception as e:
            r_loo = np.nan
            print(f"{model_name} - {target_name} LOOCV 错误: {e}")

        # 保存结果
        results.append({
            'Model': model_name,
            'Target': target_name,
            'K-fold_r': r_kfold,
            'LOOCV_r': r_loo
        })

        print(f"{target_name}: k-fold r = {r_kfold:.3f}, LOOCV r = {r_loo:.3f}")

# -------------------------------
# 5. 结果整理为 DataFrame
# -------------------------------
results_df = pd.DataFrame(results)
print("\n" + "="*60)
print("最终结果汇总 (Pearson r)")
print("="*60)
print(results_df)

# 保存结果到 CSV（可选）
results_df.to_csv('prediction_correlations_kfold_loocv.csv', index=False)
print("\n结果已保存至 'prediction_correlations_kfold_loocv.csv'")

# -------------------------------
# 6. 可视化：热图或条形图
# -------------------------------
# 转为宽表格式（每个模型一列）
pivot_kfold = results_df.pivot(index='Target', columns='Model', values='K-fold_r')
pivot_loocv = results_df.pivot(index='Target', columns='Model', values='LOOCV_r')

fig, ax = plt.subplots(1, 2, figsize=(12, 5))

# K-fold 图
pivot_kfold.plot(kind='bar', ax=ax[0], title='K-fold Pearson r', ylim=[0, 0.6])
ax[0].set_ylabel('Pearson Correlation (r)')
ax[0].legend()
ax[0].tick_params(axis='x', rotation=45)

# LOOCV 图
pivot_loocv.plot(kind='bar', ax=ax[1], title='LOOCV Pearson r', ylim=[0, 0.6])
ax[1].set_ylabel('Pearson Correlation (r)')
ax[1].legend()
ax[1].tick_params(axis='x', rotation=45)

plt.tight_layout()
plt.savefig('prediction_correlations.png', dpi=300, bbox_inches='tight')
plt.show()