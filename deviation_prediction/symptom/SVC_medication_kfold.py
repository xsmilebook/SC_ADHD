import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, roc_auc_score, roc_curve
import seaborn as sns
import warnings

warnings.filterwarnings('ignore')

# 读取数据
data_df = pd.read_csv(r'/datasets/deviation/PKU6/SCdata_deviationZ_medication_PKU6.csv')
columns = ["SC.109_deviationZ", "SC.120_deviationZ", "SC.59_deviationZ", "SC.65_deviationZ", "SC.79_deviationZ",
           "SC.84_deviationZ"]
predict_label = ["medication"]  # 假设标签列名为 medication

# 数据清洗
data_clean = data_df[columns + predict_label].dropna()

# 提取特征和标签
X = data_clean[columns].values
y = data_clean[predict_label[0]].values  # 字符串标签

# 标签编码为整数
label_encoder = LabelEncoder()
y_encoded = label_encoder.fit_transform(y)

# 特征标准化
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# 超参数网格
param_grid = {
    'C': [0.1, 1, 10, 100],
    'gamma': ['scale', 'auto', 0.001, 0.01, 0.1, 1],
    'kernel': ['rbf']
}


# 计算敏感性和特异性
def calculate_sensitivity_specificity(y_true, y_pred):
    """
    计算二分类的敏感性和特异性
    """
    cm = confusion_matrix(y_true, y_pred)

    # 对于二分类，假设:
    # cm[0,0] = TN, cm[0,1] = FP
    # cm[1,0] = FN, cm[1,1] = TP
    tn, fp, fn, tp = cm.ravel()

    sensitivity = tp / (tp + fn)  # 召回率
    specificity = tn / (tn + fp)  # 特异性

    return sensitivity, specificity


# 内外层都是10-fold分层交叉验证
def perform_nested_stratified_cv(X, y, target_name):
    print(f"\n=== Nested Stratified 10-Fold Cross-Validation for {target_name} ===")

    # 外层使用分层10-fold
    outer_cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
    # 内层也使用分层10-fold
    inner_cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

    accuracy_scores = []
    sensitivity_scores = []
    specificity_scores = []
    y_true_all = []
    y_pred_all = []

    for fold, (train_idx, test_idx) in enumerate(outer_cv.split(X, y)):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        # 内部 GridSearchCV 使用分层交叉验证
        grid_search = GridSearchCV(
            estimator=SVC(probability=True),
            param_grid=param_grid,
            cv=inner_cv,  # 内层分层10-fold
            scoring='accuracy',
            n_jobs=-1
        )
        grid_search.fit(X_train, y_train)

        best_model = grid_search.best_estimator_

        y_pred = best_model.predict(X_test)

        acc = accuracy_score(y_test, y_pred)

        # 计算敏感性和特异性（仅适用于二分类）
        if len(np.unique(y)) == 2:
            sens, spec = calculate_sensitivity_specificity(y_test, y_pred)
            sensitivity_scores.append(sens)
            specificity_scores.append(spec)
            print(f"Fold {fold + 1}: Accuracy = {acc:.4f}, Sensitivity = {sens:.4f}, Specificity = {spec:.4f}")
        else:
            print(f"Fold {fold + 1}: Accuracy = {acc:.4f}")

        accuracy_scores.append(acc)

        # 收集所有预测结果用于最终评估
        y_true_all.extend(y_test)
        y_pred_all.extend(y_pred)

    accuracy_scores = np.array(accuracy_scores)
    y_true_all = np.array(y_true_all)
    y_pred_all = np.array(y_pred_all)

    print(f"\nStratified 10-Fold CV Results:")
    print(f"Average Accuracy: {accuracy_scores.mean():.4f} (+/- {accuracy_scores.std() * 2:.4f})")

    if len(np.unique(y)) == 2:
        sensitivity_scores = np.array(sensitivity_scores)
        specificity_scores = np.array(specificity_scores)
        print(f"Average Sensitivity: {sensitivity_scores.mean():.4f} (+/- {sensitivity_scores.std() * 2:.4f})")
        print(f"Average Specificity: {specificity_scores.mean():.4f} (+/- {specificity_scores.std() * 2:.4f})")

    # 输出整体分类报告
    print("\nOverall Classification Report:")
    print(classification_report(y_true_all, y_pred_all, target_names=label_encoder.classes_))

    # 计算总体的敏感性和特异性
    if len(np.unique(y)) == 2:
        overall_sens, overall_spec = calculate_sensitivity_specificity(y_true_all, y_pred_all)
        print(f"\nOverall Sensitivity: {overall_sens:.4f}")
        print(f"Overall Specificity: {overall_spec:.4f}")

    return accuracy_scores, sensitivity_scores, specificity_scores, y_true_all, y_pred_all


# 主流程
results = {}

for i, target in enumerate(predict_label):
    y_target = y_encoded

    print(f"\n{'=' * 60}")
    print(f"Predicting {target}")
    print(f"{'=' * 60}")

    # 分层10-fold嵌套交叉验证
    acc_scores, sens_scores, spec_scores, y_true, y_pred = perform_nested_stratified_cv(X_scaled, y_target, target)

    results[target] = {
        'cv': {
            'acc': acc_scores,
            'sens': sens_scores,
            'spec': spec_scores,
            'y_true': y_true,
            'y_pred': y_pred
        }
    }

# 总结输出
print(f"\n{'=' * 60}")
print("SUMMARY RESULTS")
print(f"{'=' * 60}")

for target in predict_label:
    print(f"\n{target}:")
    print(
        f"  Accuracy:    {results[target]['cv']['acc'].mean():.4f} (+/- {results[target]['cv']['acc'].std() * 2:.4f})")

    if len(results[target]['cv']['sens']) > 0:  # 二分类情况
        print(
            f"  Sensitivity: {results[target]['cv']['sens'].mean():.4f} (+/- {results[target]['cv']['sens'].std() * 2:.4f})")
        print(
            f"  Specificity: {results[target]['cv']['spec'].mean():.4f} (+/- {results[target]['cv']['spec'].std() * 2:.4f})")

    # 绘制混淆矩阵
    y_true_final = results[target]['cv']['y_true']
    y_pred_final = results[target]['cv']['y_pred']

    cm = confusion_matrix(y_true_final, y_pred_final)
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt="d", cmap="Blues",
                xticklabels=label_encoder.classes_, yticklabels=label_encoder.classes_)
    plt.xlabel('Predicted')
    plt.ylabel('Actual')
    plt.title(f'Confusion Matrix for {target}\nAccuracy: {results[target]["cv"]["acc"].mean():.4f}')
    plt.tight_layout()
    plt.show()

    # 如果是二分类，绘制ROC曲线并计算AUC
    if len(label_encoder.classes_) == 2:
        # 重新训练模型获取概率（使用分层CV）
        outer_cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
        inner_cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

        y_true_roc = []
        y_proba_roc = []

        for train_idx, test_idx in outer_cv.split(X_scaled, y_target):
            X_train, X_test = X_scaled[train_idx], X_scaled[test_idx]
            y_train, y_test = y_target[train_idx], y_target[test_idx]

            grid_search = GridSearchCV(
                estimator=SVC(probability=True),
                param_grid=param_grid,
                cv=inner_cv,
                scoring='accuracy',
                n_jobs=-1
            )
            grid_search.fit(X_train, y_train)

            best_model = grid_search.best_estimator_
            y_proba = best_model.predict_proba(X_test)[:, 1]

            y_true_roc.extend(y_test)
            y_proba_roc.extend(y_proba)

        y_true_roc = np.array(y_true_roc)
        y_proba_roc = np.array(y_proba_roc)

        fpr, tpr, _ = roc_curve(y_true_roc, y_proba_roc)
        auc = roc_auc_score(y_true_roc, y_proba_roc)

        plt.figure(figsize=(8, 6))
        plt.plot(fpr, tpr, label=f'ROC curve (AUC = {auc:.3f})')
        plt.plot([0, 1], [0, 1], 'k--')
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(f'ROC Curve for {target}')
        plt.legend(loc='lower right')
        plt.grid(True)
        plt.tight_layout()
        plt.show()

        print(f"  AUC Score: {auc:.4f}")

# 显示类别分布
print(f"\n{'=' * 60}")
print("CLASS DISTRIBUTION")
print(f"{'=' * 60}")
unique, counts = np.unique(y_encoded, return_counts=True)
for cls, count in zip(label_encoder.classes_, counts):
    print(f"{cls}: {count} samples ({count / len(y_encoded) * 100:.1f}%)")

# 显示详细的性能指标表格
print(f"\n{'=' * 60}")
print("DETAILED PERFORMANCE METRICS")
print(f"{'=' * 60}")

for target in predict_label:
    print(f"\n{target}:")
    print(f"  Metric        Value     Std Dev    95% CI")
    print(f"  {'-' * 40}")
    acc_mean = results[target]['cv']['acc'].mean()
    acc_std = results[target]['cv']['acc'].std()
    print(
        f"  Accuracy      {acc_mean:.4f}    {acc_std:.4f}    [{acc_mean - 1.96 * acc_std:.4f}, {acc_mean + 1.96 * acc_std:.4f}]")

    if len(results[target]['cv']['sens']) > 0:
        sens_mean = results[target]['cv']['sens'].mean()
        sens_std = results[target]['cv']['sens'].std()
        spec_mean = results[target]['cv']['spec'].mean()
        spec_std = results[target]['cv']['spec'].std()

        print(
            f"  Sensitivity   {sens_mean:.4f}    {sens_std:.4f}    [{sens_mean - 1.96 * sens_std:.4f}, {sens_mean + 1.96 * sens_std:.4f}]")
        print(
            f"  Specificity   {spec_mean:.4f}    {spec_std:.4f}    [{spec_mean - 1.96 * spec_std:.4f}, {spec_mean + 1.96 * spec_std:.4f}]")