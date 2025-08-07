# -*- coding: utf-8 -*-
import os
import scipy.io as sio
import numpy as np
import pandas as pd

def calculate_significance(base_result_folder, perm_result_folder_base, num_permutations=1000, metric='Mean_Corr'):
    """
    读取原始预测和置换检验结果，计算显著性 p-value。

    Parameters:
    ----------
    base_result_folder : str
        原始预测（未置换）结果的根目录。
        例如: '/ibmgpfs/cuizaixu_lab/xuhaoshu/SC_ADHD/datasets/prediction/activation/nbackFP.A_deviationZ/LR'
        假设 Res_NFold.mat 在这个目录下。
    perm_result_folder_base : str
        置换检验结果的根目录基础路径。
        例如: '/ibmgpfs/cuizaixu_lab/xuhaoshu/SC_ADHD/datasets/prediction/activation/nbackFP.A_deviationZ/LR_Permutation'
        假设里面有 Time_0, Time_1, ..., Time_999 等子目录。
    num_permutations : int, optional
        计划或实际运行的置换检验次数 (默认 1000)。
    metric : str, optional
        要用于计算显著性的性能指标名称 (默认 'Mean_Corr')。

    Returns:
    -------
    dict
        包含 'observed_metric', 'perm_metrics', 'p_value' 的字典。
    """

    # 1. 读取原始预测（非置换）的性能指标
    original_result_file = os.path.join(base_result_folder, 'Res_NFold.mat')
    if not os.path.exists(original_result_file):
        raise FileNotFoundError(f"原始结果文件未找到: {original_result_file}")

    try:
        original_data = sio.loadmat(original_result_file)
        observed_metric = float(original_data[metric]) # 确保是标量
    except Exception as e:
        raise RuntimeError(f"读取原始结果文件失败 {original_result_file}: {e}")

    print(f"原始预测 {metric}: {observed_metric}")

    # 2. 读取 1000 次置换检验的性能指标
    perm_metrics = []
    missing_files = []
    error_files = []

    for i in range(num_permutations):
        perm_folder = os.path.join(perm_result_folder_base, f'Time_{i}')
        perm_result_file = os.path.join(perm_folder, 'Res_NFold.mat')

        if not os.path.exists(perm_result_file):
            missing_files.append(perm_result_file)
            # 可以选择跳过或报错
            # print(f"警告: 置换结果文件未找到: {perm_result_file}")
            continue # 跳过这个缺失的文件

        try:
            perm_data = sio.loadmat(perm_result_file)
            perm_metric = float(perm_data[metric]) # 确保是标量
            perm_metrics.append(perm_metric)
        except Exception as e:
            error_files.append((perm_result_file, str(e)))
            # print(f"警告: 读取置换结果文件失败 {perm_result_file}: {e}")
            continue # 跳过这个读取失败的文件

    if missing_files:
        print(f"警告: 未找到 {len(missing_files)} 个置换结果文件。")
        # 可以选择根据情况处理，例如如果缺失太多则报错
        # if len(missing_files) > num_permutations * 0.1: # 超过10%缺失
        #     raise RuntimeError(f"缺失过多置换结果文件 ({len(missing_files)} / {num_permutations})")

    if error_files:
        print(f"警告: 读取 {len(error_files)} 个置换结果文件时出错。")

    if not perm_metrics:
        raise RuntimeError("没有成功读取到任何置换检验结果。")

    perm_metrics = np.array(perm_metrics)
    num_valid_perm = len(perm_metrics)
    print(f"成功读取 {num_valid_perm} 个有效的置换检验结果。")

    # 3. 计算 p-value (双侧检验)
    # 对于双侧检验，我们关心的是绝对值 |metric|
    # 或者更常见的是，如果原始metric是相关性（越高越好），我们看有多少置换的metric >= 原始metric
    # 这里假设 metric 越高越好 (如 Corr), 使用单侧检验 (>=)
    # 如果是 MAE (越低越好), 则逻辑需要调整 (<=)

    if metric.lower().endswith('corr'): # 假设 Corr 类指标越高越好
        # 计算 >= observed_metric 的置换次数
        num_better_or_equal = np.sum(perm_metrics >= observed_metric)
        # 使用 +1 平滑避免 p=0
        p_value = (num_better_or_equal + 1) / (num_valid_perm + 1)
    elif 'mae' in metric.lower(): # 假设 MAE 类指标越低越好
         # 计算 <= observed_metric 的置换次数
        num_worse_or_equal = np.sum(perm_metrics <= observed_metric)
        p_value = (num_worse_or_equal + 1) / (num_valid_perm + 1)
    else:
        # 默认假设越高越好
        print(f"警告: 无法确定指标 '{metric}' 的方向，默认按越高越好处理。")
        num_better_or_equal = np.sum(perm_metrics >= observed_metric)
        p_value = (num_better_or_equal + 1) / (num_valid_perm + 1)


    print(f"观测到的 {metric}: {observed_metric}")
    print(f"置换检验 {metric} 的均值: {np.mean(perm_metrics):.6f}")
    print(f"置换检验 {metric} 的标准差: {np.std(perm_metrics):.6f}")
    print(f"计算得到的 p-value (单侧, >=): {p_value:.6f}")

    return {
        'observed_metric': observed_metric,
        'perm_metrics': perm_metrics, # 返回数组，方便后续分析
        'p_value': p_value,
        'num_permutations_run': num_valid_perm, # 实际用于计算的有效置换次数
        'metric_name': metric
    }


if __name__ == "__main__":


    targetStr_List = [
        "nbackFP.A_deviationZ",
        "nbackFP.B_deviationZ",
        "nbackSM.A_deviationZ",
        "nbackSM.B_deviationZ",
        "nbackVA.B_deviationZ"
    ]

    base_path = '/ibmgpfs/cuizaixu_lab/xuhaoshu/SC_ADHD/datasets/prediction/activation'

    results_summary = []

    for targetStr in targetStr_List:
        print(f"\n--- 处理目标: {targetStr} ---")

        # 原始预测结果文件夹 (你之前运行未置换预测的目录)
        original_result_folder = os.path.join(base_path, targetStr, 'LR', "Time_0")
        # 置换检验结果文件夹 (你之前运行 1000 次置换的目录)
        # 注意：你需要确保这个路径是正确的，并且里面包含了 Time_0 到 Time_999
        # 如果你是在同一个 'LR' 文件夹下运行的置换，可能需要调整结构，
        # 例如创建一个子文件夹 'Permutation_Results' 来存放置换结果
        # 假设你为置换检验创建了一个独立的文件夹
        perm_result_folder_base = os.path.join(base_path, targetStr, 'LR_Permutation') 

        try:
            # --- 执行显著性计算 ---
            sig_result = calculate_significance(
                base_result_folder=original_result_folder,
                perm_result_folder_base=perm_result_folder_base,
                num_permutations=1000, # 确保与你运行的次数一致
                metric='Mean_Corr'     # 根据需要更改指标
            )

            # --- 保存或打印结果 ---
            summary_entry = {
                'target': targetStr,
                'observed_Mean_Corr': sig_result['observed_metric'],
                'p_value': sig_result['p_value'],
                'num_permutations_used': sig_result['num_permutations_run'], # 实际用到的数量
                # 'perm_metrics': sig_result['perm_metrics'] # 如果需要保存所有置换结果可以取消注释
            }
            results_summary.append(summary_entry)

            print(f"目标 '{targetStr}' 的显著性检验完成。")
            print("-" * 20)

        except (FileNotFoundError, RuntimeError) as e:
            print(f"处理目标 '{targetStr}' 时出错: {e}")
            # 可以选择跳过或停止
            continue

    # --- 输出所有目标的汇总结果 ---
    if results_summary:
        print("\n=== 所有目标的显著性检验结果汇总 ===")
        summary_df = pd.DataFrame(results_summary)
        print(summary_df.to_string(index=False))

        # 可选：保存到 CSV 文件
        output_csv_path = os.path.join(base_path, 'LR_Significance_Results.csv')
        summary_df.to_csv(output_csv_path, index=False)
        print(f"\n汇总结果已保存到: {output_csv_path}")
    else:
        print("\n没有成功处理任何目标。")
