import pickle

import pandas
import pandas as pd
import os
import argparse
from cv_predict import run_svr_nested_cv
import random

def parse_args():
    parser = argparse.ArgumentParser(description="Run SVR nested CV for task fMRI activation prediction.")
    parser.add_argument('--permutation', action='store_true', help='Enable permutation mode (shuffle y)')
    parser.add_argument('--time_id', type=str, default=None, help='Time ID for result folder (e.g., Time_0)')

    return parser.parse_args()

def main():
    args = parse_args()

    demo_table = pd.read_csv("D:/code/SC_ADHD/datasets/ABCD/task-fMRI/activation/demo_table.csv")
    result_folder = r"D:\code\SC_ADHD\datasets\results\task_act_prediction"
    permutation_flag = args.permutation
    time_id = args.time_id
    # task_id_list = ["nback", "sst", "MIDallvn", "MIDalrvn", "MIDrpvnf"]
    task_name = "nback"
    task_id = "nback"

    print(f"Starting job: task={task_id}, permutation={permutation_flag}, time_id={time_id}")

    feature_columns = [
        "SC.109_deviationZ", "SC.111_deviationZ", "SC.116_deviationZ", "SC.118_deviationZ",
        "SC.119_deviationZ", "SC.120_deviationZ", "SC.18_deviationZ", "SC.57_deviationZ",
        "SC.59_deviationZ", "SC.60_deviationZ", "SC.65_deviationZ", "SC.70_deviationZ",
        "SC.79_deviationZ", "SC.84_deviationZ", "SC.8_deviationZ", "SC.93_deviationZ",
        "SC.98_deviationZ"
    ]

    predict_columns = [
        "FP.A_h",
        "FP.B_h",
        "SM.A_h",
        "SM.B_h",
        "VA.B_h"
    ]

    predict_columns = [task_id + region_id for region_id in predict_columns]


    demo_table['sex'] = demo_table['sex'].astype('category')
    demo_table['if_TD'] = demo_table['if_TD'].astype('category')
    confounds_list = ["age", "sex", "meanFD"+task_name, "if_TD"]
    confounds_data = demo_table[confounds_list]

    demo_table = demo_table[feature_columns + predict_columns]



    for target_column in predict_columns:
        task_result_folder = os.path.join(result_folder, task_id, target_column)
        if permutation_flag:
            task_result_folder = os.path.join(task_result_folder, "permutation", "linear", "Time_" + time_id)
        else:
            task_result_folder = os.path.join(task_result_folder, "formal", "linear", "Time_" + time_id)
        os.makedirs(task_result_folder, exist_ok=True)

        result = run_svr_nested_cv(
            data_df=demo_table,
            feature_cols=feature_columns,
            target_col=target_column,
            outer_folds=5,
            inner_folds=5,
            n_runs=1,
            random_state_base=abs(hash(time_id)) % 100000 + 1,
            confounds_df=confounds_data,
            permutation_flag=permutation_flag,
            verbose=True,
        )

        result_file = os.path.join(task_result_folder, "result.pkl")
        with open(result_file, 'wb') as f:
            pickle.dump(result, f)

        summary = {
            'target': result['target'],
            'mean_correlation': result['mean_correlation'],
            'std_correlation': result['std_correlation'],
            'n_runs': result['n_runs'],
            'permutation_flag': result['permutation_flag'],
            'time_id': args.time_id,
            'task_id': task_id
        }
        pd.DataFrame([summary]).to_csv(os.path.join(task_result_folder, "summary.csv"), index=False)


if __name__ == "__main__":
    main()


