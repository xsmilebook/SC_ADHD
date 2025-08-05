import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.svm import SVC
from sklearn.model_selection import cross_val_score, KFold, LeaveOneOut, GridSearchCV
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, roc_auc_score, roc_curve
import seaborn as sns
import warnings

warnings.filterwarnings('ignore')

# 读取数据
data_df = pd.read_csv(r'/datasets/deviation/PKU6/SCdata_deviationZ_medication_PKU6.csv')
columns = ["SC.109_deviationZ", "SC.120_deviationZ", "SC.59_deviationZ", "SC.65_deviationZ", "SC.79_deviationZ",
           "SC.84_deviationZ"]
predict_label = ["medication"]  # 假设你有一个字符串标签列，如 "ClassA" 和 "ClassB"

# 数据清洗
data_clean = data_df[columns + predict_label].dropna()

# 提取特征和标签
X = data_clean[columns].values
y = data_clean[predict_label[0]].values  # 字符串标签

# 标签编码为整数（SVC 需要数值型标签）
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

# 内部交叉验证调参 + 外部评估
def perform_nested_cv(X, y, target_name):
    print(f"\n=== Nested Cross-Validation for {target_name} ===")

    outer_cv = KFold(n_splits=10, shuffle=True, random_state=42)
    inner_cv = KFold(n_splits=5, shuffle=True, random_state=42)

    accuracy_scores = []

    for train_idx, test_idx in outer_cv.split(X):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        # 内部 GridSearchCV
        grid_search = GridSearchCV(
            estimator=SVC(probability=True),  # probability=True 用于 ROC 曲线
            param_grid=param_grid,
            cv=inner_cv,
            scoring='accuracy',
            n_jobs=-1
        )
        grid_search.fit(X_train, y_train)

        best_model = grid_search.best_estimator_

        y_pred = best_model.predict(X_test)
        acc = accuracy_score(y_test, y_pred)

        accuracy_scores.append(acc)

    accuracy_scores = np.array(accuracy_scores)

    print(f"Accuracy scores: {accuracy_scores}")
    print(f"Average Accuracy: {accuracy_scores.mean():.4f} (+/- {accuracy_scores.std() * 2:.4f})")

    return accuracy_scores


# LOOCV
def perform_loocv(X, y, target_name):
    print(f"\n=== LOOCV for {target_name} ===")

    loo = LeaveOneOut()
    y_true_loo = []
    y_pred_loo = []
    y_proba_loo = []

    inner_cv = KFold(n_splits=5, shuffle=True, random_state=42)

    for train_index, test_index in loo.split(X):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]

        # 内部 GridSearchCV
        grid_search = GridSearchCV(
            estimator=SVC(probability=True),
            param_grid=param_grid,
            cv=inner_cv,
            scoring='accuracy',
            n_jobs=-1
        )
        grid_search.fit(X_train, y_train)

        best_model = grid_search.best_estimator_
        y_pred = best_model.predict(X_test)
        y_proba = best_model.predict_proba(X_test)[:, 1]  # 假设正类为 1

        y_true_loo.append(y_test[0])
        y_pred_loo.append(y_pred[0])
        y_proba_loo.append(y_proba[0])

    y_true_loo = np.array(y_true_loo)
    y_pred_loo = np.array(y_pred_loo)
    y_proba_loo = np.array(y_proba_loo)

    acc = accuracy_score(y_true_loo, y_pred_loo)

    print(f"Accuracy: {acc:.4f}")
    print(classification_report(y_true_loo, y_pred_loo, target_names=label_encoder.classes_))

    # 混淆矩阵
    cm = confusion_matrix(y_true_loo, y_pred_loo)
    sns.heatmap(cm, annot=True, fmt="d", cmap="Blues",
                xticklabels=label_encoder.classes_, yticklabels=label_encoder.classes_)
    plt.xlabel('Predicted')
    plt.ylabel('Actual')
    plt.title(f'Confusion Matrix for {target_name}')
    plt.show()

    # ROC Curve
    fpr, tpr, _ = roc_curve(y_true_loo, y_proba_loo)
    auc = roc_auc_score(y_true_loo, y_proba_loo)

    plt.figure()
    plt.plot(fpr, tpr, label=f'ROC curve (AUC = {auc:.2f})')
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f'ROC Curve for {target_name}')
    plt.legend(loc='lower right')
    plt.grid(True)
    plt.show()

    return acc, y_true_loo, y_pred_loo


# 主流程
results = {}

for i, target in enumerate(predict_label):
    y_target = y_encoded

    print(f"\n{'=' * 50}")
    print(f"Predicting {target}")
    print(f"{'=' * 50}")

    # 嵌套交叉验证
    acc_nested = perform_nested_cv(X_scaled, y_target, target)

    # LOOCV
    acc_loocv, y_true, y_pred = perform_loocv(X_scaled, y_target, target)

    results[target] = {
        'nested_cv': {'acc': acc_nested},
        'loocv': {'acc': acc_loocv, 'y_true': y_true, 'y_pred': y_pred}
    }

# 总结输出
print(f"\n{'=' * 60}")
print("SUMMARY RESULTS")
print(f"{'=' * 60}")

for target in predict_label:
    print(f"\n{target}:")
    print(f"  Nested CV - Accuracy: {results[target]['nested_cv']['acc'].mean():.4f}")
    print(f"  LOOCV     - Accuracy: {results[target]['loocv']['acc']:.4f}")