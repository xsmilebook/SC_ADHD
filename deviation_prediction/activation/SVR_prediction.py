import os
import pandas as pd
from sklearn.inspection import permutation_importance
from regress_covariates import regress_covariates
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold, GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LinearRegression
from sklearn.svm import SVR
from scipy.stats import pearsonr, false_discovery_control
import numpy as np
import matplotlib.pyplot as plt

data_df = pd.read_csv('D:/code/SC_ADHD/datasets/ABCD/task-fMRI/deviation/SCdeviation_nbackActivation_deviation_ADHD_TDtest.csv', low_memory=False)

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

data_df[predict_columns] = regress_covariates(data_df[predict_columns].copy(), covariates_data, covariates_list, category_cov_list)

# 3. permutation
permutation_flag = False

# 4. select features
temp_df = data_df[feature_columns + predict_columns]
correlation_matrix = temp_df.corr()
selected_feature_col_dict = {}
for col in predict_columns:
    correlations = correlation_matrix.loc[feature_columns, col].abs()

    top_5 = correlations.nlargest(5)
    # print(top_5)
    selected_feature_col_dict[col] = top_5.index.tolist()

model_name = "SVR"
model = SVR(kernel='rbf', C=0.6)

all_runs_avg_correlations = {col: [] for col in predict_columns}
all_importances = []

for run in range(10):
    kf = KFold(n_splits=10, shuffle=True, random_state=run)
    splits = []

    for train_index, test_index in kf.split(data_df):
        X_train = data_df.iloc[train_index][feature_columns]
        X_test = data_df.iloc[test_index][feature_columns]
        y_train = data_df.iloc[train_index][predict_columns]
        y_test = data_df.iloc[test_index][predict_columns]
        splits.append((X_train, X_test, y_train, y_test))

    model_predictions = {}

    for fold_index, (X_train, X_test, y_train, y_test) in enumerate(splits):
        model_predictions[fold_index] = {}
        for predict_column in predict_columns:
            model = SVR(kernel='rbf', C=0.6)

            selected_features = selected_feature_col_dict[predict_column]
            X_train_sel = X_train[selected_features]
            X_test_sel = X_test[selected_features]

            scaler = StandardScaler()
            X_train_scaled = scaler.fit_transform(X_train_sel)
            X_test_scaled = scaler.transform(X_test_sel)

            model.fit(X_train_scaled, y_train[predict_column])

            predictions = model.predict(X_test_scaled)
            model_predictions[fold_index][predict_column] = predictions

            perm_result = permutation_importance(
                model, X_test_scaled, y_test[predict_column],
                n_repeats=5,
                random_state=run,
                scoring='neg_mean_squared_error'
            )

            importances = perm_result.importances_mean

            imp_dict = {
                feat: imp for feat, imp in zip(selected_features, importances)
            }
            imp_dict['run'] = run
            imp_dict['fold'] = fold_index
            imp_dict['target'] = predict_column

            all_importances.append(imp_dict)

    average_correlations = {}
    for predict_column in predict_columns:
        correlations = []
        for fold_index in range(len(splits)):
            y_true = splits[fold_index][3][predict_column]
            y_pred = model_predictions[fold_index][predict_column]
            corr = np.corrcoef(y_true, y_pred)[0, 1]
            correlations.append(corr)
        average_correlations[predict_column] = np.mean(correlations)

    print(f"Run {run + 1} average correlations:")
    for col in predict_columns:
        print(f"  {col}: {average_correlations[col]:.4f}")
        all_runs_avg_correlations[col].append(average_correlations[col])

print("\nFinal average correlations over 10 runs:")
for col in predict_columns:
    avg_over_runs = np.mean(all_runs_avg_correlations[col])
    print(f"  {col}: {avg_over_runs:.4f}")

importance_df = pd.DataFrame(all_importances)

output_dir = os.path.join("D:/code/SC_ADHD/datasets/prediction/activation_deviation", model_name)
os.makedirs(output_dir, exist_ok=True)

for target in predict_columns:
    target_df = importance_df[importance_df['target'] == target]
    target_df.to_csv(f"{output_dir}/feature_importance_{target}.csv", index=False)

importance_df.to_csv(f"{output_dir}/all_feature_importance.csv", index=False)


for target in predict_columns:
    target_data = importance_df[importance_df['target'] == target]
    selected_features = selected_feature_col_dict[target]

    avg_imp = target_data[selected_features].mean(axis=0).sort_values(ascending=False)

    plt.figure(figsize=(10, 6))
    avg_imp.plot(kind='bar', color='skyblue', edgecolor='black')
    plt.title(f'Permutation Importance - {target}')
    plt.ylabel('Importance (MSE reduction)')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/importance_{target}.png", dpi=300)
    plt.close()




