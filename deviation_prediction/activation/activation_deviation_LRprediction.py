#coding: utf-8
import scipy.io as sio
import numpy as np
import pandas as pd
import os
import sys
from LR_Kfold_MultiTimes import LR_KFold_RandomCV_MultiTimes

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
    SubjectsData = data_df

    # 2. subject label: cognition score
    dimention = targetStr
    label = data_df[dimention]
    y_label = np.array(label)
    Predict_label = y_label

    # 3. covariates
    Covariates = data_df["age", "sex"].values
    Covariates = Covariates[:, 1:].astype(float)

    # Range of parameters
    ComponentNumber_Range = np.arange(10) + 1
    FoldQuantity = 5
    Parallel_Quantity = 1
    # CVtimes = 101

    # Predict
    ResultantFolder = ResultsFolder + '/LR'
    # LR_KFold_RandomCV_MultiTimes(SubjectsData, Predict_label, Covariates, FoldQuantity, ComponentNumber_Range, CVtimes, ResultantFolder, Parallel_Quantity, 0)

    # Permutation
    LR_KFold_RandomCV_MultiTimes(SubjectsData, Predict_label, Covariates, FoldQuantity, 1000, ResultantFolder, Parallel_Quantity, 1, 'all.q')


