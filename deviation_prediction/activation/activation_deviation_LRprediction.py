#coding: utf-8
import scipy.io as sio
import numpy as np
import pandas as pd
import os
import sys
from LR_Kfold_MultiTimes import LR_KFold_RandomCV_MultiTimes
sys.path.append('/ibmgpfs/cuizaixu_lab/xuhaoshu/SC_ADHD/code/deviation_prediction/activation')

feature_columns = [
    "SC.109_deviationZ", "SC.111_deviationZ", "SC.116_deviationZ", "SC.118_deviationZ",
    "SC.119_deviationZ", "SC.120_deviationZ", "SC.18_deviationZ", "SC.57_deviationZ",
    "SC.59_deviationZ", "SC.60_deviationZ", "SC.65_deviationZ", "SC.70_deviationZ",
    "SC.79_deviationZ", "SC.84_deviationZ", "SC.8_deviationZ", "SC.93_deviationZ",
    "SC.98_deviationZ"
]

targetStr_List = [
    "nbackFP.A_deviationZ",
    "nbackFP.B_deviationZ",
    "nbackSM.A_deviationZ",
    "nbackSM.B_deviationZ",
    "nbackVA.B_deviationZ"
]

for targetStr in targetStr_List:
    ResultsFolder = '/ibmgpfs/cuizaixu_lab/xuhaoshu/SC_ADHD/datasets/prediction/activation/'+ targetStr
    # Import data
    # 1. atlas loading
    datapath = '/ibmgpfs/cuizaixu_lab/xuhaoshu/SC_ADHD/datasets/ABCD/task_activation/deviation/SCdeviation_nbackActivation_deviation_ADHD_TDtest.csv'
    data_df = pd.read_csv(datapath, low_memory=False)

    temp_df = data_df[feature_columns + [targetStr]]
    correlation_matrix = temp_df.corr()
    selected_feature_cols = []
    correlations = correlation_matrix.loc[feature_columns, targetStr].abs()

    top_5 = correlations.nlargest(5)
    # print(top_5)
    selected_feature_cols = top_5.index.tolist()
    SubjectsData = data_df[selected_feature_cols].values

    # 2. subject label: cognition score
    dimention = targetStr
    label = data_df[dimention]
    y_label = np.array(label)
    Predict_label = y_label

    # 3. covariates
    Covariates = data_df[['age', 'sex']].values
    Covariates = Covariates[:, 0:].astype(float)

    # Range of parameters
    FoldQuantity = 10
    Parallel_Quantity = 1
    # CVtimes = 101

    # Predict
    ResultantFolder = ResultsFolder + '/LR'
    # LR_KFold_RandomCV_MultiTimes(SubjectsData, Predict_label, Covariates, FoldQuantity, ComponentNumber_Range, CVtimes, ResultantFolder, Parallel_Quantity, 0)

    # Permutation
    LR_KFold_RandomCV_MultiTimes(SubjectsData, Predict_label, Covariates, FoldQuantity, 10, ResultantFolder, Parallel_Quantity, False, '')


