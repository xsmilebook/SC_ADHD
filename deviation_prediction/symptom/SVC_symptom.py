import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.svm import SVC
from sklearn.model_selection import cross_val_score, KFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, confusion_matrix, ConfusionMatrixDisplay
import warnings


warnings.filterwarnings('ignore')


data_df = pd.read_csv(r'/datasets/deviation/PKU6/SCdata_deviationZ_medication_PKU6.csv')
# columns = ["SC.109_deviationZ", "SC.120_deviationZ", "SC.59_deviationZ", "SC.65_deviationZ", "SC.79_deviationZ", "SC.84_deviationZ"]
columns = [item for item in data_df.columns if item.endswith('deviationZ')]
predict_label = ["IA_delta", "HI_delta"]


data_clean = data_df[columns + predict_label].dropna()


X = data_clean[columns].values

scaler_X = StandardScaler()
X_scaled = scaler_X.fit_transform(X)

svc_model = SVC(kernel='rbf', C=1.0, gamma='scale')

# 10折交叉验证函数（分类）
def perform_10fold_cv_classification(X, y, target_name):
    print(f"\n=== 10-Fold CV for {target_name} ===")

    # 10折交叉验证
    kf = KFold(n_splits=10, shuffle=True, random_state=42)

    # 计算准确率
    acc_scores = cross_val_score(svc_model, X, y, cv=kf, scoring='accuracy')

    print(f"Accuracy scores: {acc_scores}")
    print(f"Average Accuracy: {acc_scores.mean():.4f} (+/- {acc_scores.std() * 2:.4f})")

    return acc_scores


# LOOCV函数（分类）
def perform_loocv_classification(X, y, target_name):
    print(f"\n=== LOOCV for {target_name} ===")

    from sklearn.model_selection import LeaveOneOut

    loo = LeaveOneOut()
    y_true_loo = []
    y_pred_loo = []

    for train_index, test_index in loo.split(X):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]

        svc_model.fit(X_train, y_train)
        y_pred = svc_model.predict(X_test)

        y_true_loo.append(y_test[0])
        y_pred_loo.append(y_pred[0])

    y_true_loo = np.array(y_true_loo)
    y_pred_loo = np.array(y_pred_loo)

    acc = accuracy_score(y_true_loo, y_pred_loo)
    print(f"LOOCV Accuracy: {acc:.4f}")

    # 混淆矩阵
    cm = confusion_matrix(y_true_loo, y_pred_loo)
    disp = ConfusionMatrixDisplay(confusion_matrix=cm)
    disp.plot(cmap='Blues')
    plt.title(f'LOOCV Confusion Matrix for {target_name}')
    plt.show()

    return acc, y_true_loo, y_pred_loo


results = {}

for i, target in enumerate(predict_label):
    y_reg = data_clean[target].values
    y_class = (y_reg > 0).astype(int)

    print(f"\n{'=' * 50}")
    print(f"Classifying {target} (positive vs negative)")
    print(f"{'=' * 50}")

    # 10折交叉验证
    acc_10fold = perform_10fold_cv_classification(X_scaled, y_class, target)

    # LOOCV
    acc_loocv, y_true, y_pred = perform_loocv_classification(X_scaled, y_class, target)

    # 保存结果
    results[target] = {
        '10fold': {'acc': acc_10fold},
        'loocv': {'acc': acc_loocv, 'y_true': y_true, 'y_pred': y_pred}
    }

# 汇总结果
print(f"\n{'=' * 60}")
print("SUMMARY RESULTS (Classification)")
print(f"{'=' * 60}")

for target in predict_label:
    print(f"\n{target}:")
    print(f"  10-Fold CV - Accuracy: {results[target]['10fold']['acc'].mean():.4f}")
    print(f"  LOOCV      - Accuracy: {results[target]['loocv']['acc']:.4f}")