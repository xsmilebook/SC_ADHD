import os
import pickle
import pandas as pd

results_folder = r"D:\code\SC_ADHD\datasets\results\task_act_prediction"
task_id_list = ["nback", "sst", "MIDallvn", "MIDalrvn", "MIDrpvnf"]
predict_columns = [
    "FP.A_h",
    "FP.B_h",
    "SM.A_h",
    "SM.B_h",
    "VA.B_h"
]
permutation_flag = False

summary_per_condition = []


for task_id in task_id_list:
    for region_base in predict_columns:
        full_column = task_id + region_base

        if permutation_flag:
            tag = "permutation"
        else:
            tag = "formal"

        model_list = ["linear", "rbf"]
        for model in model_list:
            data_folder = os.path.join(results_folder, task_id, full_column, tag, model)

            # 检查路径是否存在
            if not os.path.exists(data_folder):
                print(f"Warning: Path not found: {data_folder}")
                continue

            corr_list = []  # 收集该条件下的所有 run 的 corr
            run_count = 0

            # 遍历 Time_0, Time_1, ..., 子文件夹
            for folder_name in os.listdir(data_folder):
                folder_path = os.path.join(data_folder, folder_name)
                if not os.path.isdir(folder_path):
                    continue

                result_file = os.path.join(folder_path, "result.pkl")
                if not os.path.exists(result_file):
                    print(f"Missing result file: {result_file}")
                    continue

                try:
                    with open(result_file, 'rb') as f:
                        result = pickle.load(f)
                        corr = result['mean_correlation']  # 假设这是 outer CV 平均 corr
                        corr_list.append(corr)
                        run_count += 1
                except Exception as e:
                    print(f"Error reading {result_file}: {e}")

            # 如果没有有效运行，跳过
            if len(corr_list) == 0:
                print(f"No valid runs for {task_id}, {full_column}, {model}")
                continue

            # ✅ 在 101 次运行上取平均，作为该条件下的预测性能
            avg_corr = sum(corr_list) / len(corr_list)
            summary_per_condition.append({
                'task_id': task_id,
                'region': region_base,
                'model': model,
                'num_runs': run_count,
                'mean_correlation': round(avg_corr, 4)
            })

            print(f"{task_id:8} | {region_base:7} | {model:6} | N={run_count:2} | avg corr = {avg_corr:.4f}")

# 转为 DataFrame
df_summary = pd.DataFrame(summary_per_condition)

# 保存每个条件的平均性能
output_csv = os.path.join(results_folder, "summary_per_condition.csv")
df_summary.to_csv(output_csv, index=False)
print(f"\n✅ 每个 (任务, 脑区, 模型) 的平均预测性能已保存至:\n   {output_csv}")

# ✅ 计算所有任务、所有脑区、所有模型的 **总体平均预测 corr**
overall_mean = df_summary['mean_correlation'].mean()
overall_std = df_summary['mean_correlation'].std()

print(f"\n🎯 总体平均预测相关系数（所有 {len(df_summary)} 个条件）: {overall_mean:.4f}")
print(f"🎯 标准差（反映不同条件间稳定性）: {overall_std:.4f}")

# （可选）按任务和模型分组的平均
print("\n📈 按任务和模型分组的平均 corr:")
grouped = df_summary.groupby(['task_id', 'model'])['mean_correlation'].agg(['mean', 'std', 'count']).round(4)
print(grouped)

grouped_csv = os.path.join(results_folder, "grouped_summary_by_task_model.csv")
grouped.to_csv(grouped_csv)
print(f"📈 分组结果已保存至:\n   {grouped_csv}")