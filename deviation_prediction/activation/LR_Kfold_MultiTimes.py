# -*- coding: utf-8 -*-
import os
import scipy.io as sio
import numpy as np
import pandas as pd
from sklearn import linear_model
from sklearn.preprocessing import StandardScaler
from joblib import Parallel, delayed
import statsmodels.formula.api as sm

# CodesPath = '/ibmgpfs/cuizaixu_lab/xuhaoshu/SC_ADHD/code/deviation_prediction/activation'

def LR_KFold_RandomCV_MultiTimes(Subjects_Data, Subjects_Score, Covariates, Fold_Quantity,
                                    CVRepeatTimes, ResultantFolder, Parallel_Quantity, Permutation_Flag,
                                    RandIndex_File_List=''):
    if not os.path.exists(ResultantFolder):
        os.makedirs(ResultantFolder)
    Subjects_Data_Mat = {'Subjects_Data': Subjects_Data}
    Subjects_Data_Mat_Path = ResultantFolder + '/Subjects_Data.mat'
    sio.savemat(Subjects_Data_Mat_Path, Subjects_Data_Mat)
    Finish_File = []
    Corr_MTimes = np.zeros(CVRepeatTimes)
    MAE_MTimes = np.zeros(CVRepeatTimes)
    for i in np.arange(CVRepeatTimes):
        ResultantFolder_TimeI = ResultantFolder + '/Time_' + str(i)
        if not os.path.exists(ResultantFolder_TimeI):
            os.makedirs(ResultantFolder_TimeI)
        if not os.path.exists(ResultantFolder_TimeI + '/Res_NFold.mat'):
            if RandIndex_File_List != '':
                RandIndex_File = RandIndex_File_List[i]
            else:
                RandIndex_File = ''
            Configuration_Mat = {'Subjects_Data_Mat_Path': Subjects_Data_Mat_Path, 'Subjects_Score': Subjects_Score,
                                 'Covariates': Covariates, 'Fold_Quantity': Fold_Quantity,
                                 'CVRepeatTimes': CVRepeatTimes,
                                 'ResultantFolder_TimeI': ResultantFolder_TimeI, 'Parallel_Quantity': Parallel_Quantity,
                                 'Permutation_Flag': Permutation_Flag, 'RandIndex_File': RandIndex_File};
            sio.savemat(ResultantFolder_TimeI + '/Configuration.mat', Configuration_Mat)
            system_cmd = 'python3 -c ' + '\'import sys;\
                sys.path.append("' + os.path.dirname(os.path.abspath(__file__)) + '");\
                from LR_CZ_Random_RegressCovariates import LR_KFold_RandomCV_MultiTimes_Sub; \
                import os;\
                import scipy.io as sio;\
                Configuration = sio.loadmat("' + ResultantFolder_TimeI + '/Configuration.mat");\
                Subjects_Data_Mat_Path = Configuration["Subjects_Data_Mat_Path"];\
                Subjects_Score = Configuration["Subjects_Score"];\
                Covariates = Configuration["Covariates"];\
                Fold_Quantity = Configuration["Fold_Quantity"];\
                # ComponentNumber_Range = Configuration["ComponentNumber_Range"];\
                ResultantFolder_TimeI = Configuration["ResultantFolder_TimeI"];\
                Permutation_Flag = Configuration["Permutation_Flag"];\
                RandIndex_File = Configuration["RandIndex_File"];\
                Parallel_Quantity = Configuration["Parallel_Quantity"];\
                LR_KFold_RandomCV_MultiTimes_Sub(Subjects_Data_Mat_Path[0], Subjects_Score[0], Covariates, Fold_Quantity[0][0], ResultantFolder_TimeI[0], Parallel_Quantity[0][0], Permutation_Flag[0][0], RandIndex_File)\' ';
            system_cmd = system_cmd + ' > "' + ResultantFolder_TimeI + '/Time_' + str(i) + '.log" 2>&1\n'
            Finish_File.append(ResultantFolder_TimeI + '/Res_NFold.mat')
            script = open(ResultantFolder_TimeI + '/script.sh', 'w')
            # # Submit jobs
            script.write('#!/bin/bash\n');
            script.write('#SBATCH --job-name=LINca' + str(i) + '\n')
            script.write('#SBATCH --nodes=1\n')
            script.write('#SBATCH --ntasks=1\n')
            script.write('#SBATCH --cpus-per-task=2\n')
            script.write('#SBATCH --mem-per-cpu 5G\n')
            script.write('#SBATCH -p q_fat_c\n')
            script.write('#SBATCH -q high_c\n')
            script.write('#SBATCH -o ' + ResultantFolder_TimeI + '/job.%j.out\n')
            script.write('#SBATCH -e ' + ResultantFolder_TimeI + '/job.%j.error.txt\n')
            script.write(system_cmd)
            script.close()
            os.system('chmod +x ' + ResultantFolder_TimeI + '/script.sh')
            os.system('sbatch ' + ResultantFolder_TimeI + '/script.sh')

            # Option = ' -V -o "' + ResultantFolder_TimeI + '/RandomCV_' + str(i) + '.o" -e "' + ResultantFolder_TimeI + '/RandomCV_' + str(i) + '.e" ';
            # os.system(' -l h_vmem=10G,s_vmem=10G -q ' + Queue + ' -N RandomCV_' + str(i) + Option + ResultantFolder_TimeI + '/script.sh')

def LR_KFold_RandomCV_MultiTimes_Sub(Subjects_Data_Mat_Path, Subjects_Score, Covariates, Fold_Quantity,
                                        # ComponentNumber_Range,
                                        ResultantFolder, Parallel_Quantity,
                                        Permutation_Flag, RandIndex_File=''):
    data = sio.loadmat(Subjects_Data_Mat_Path)
    Subjects_Data = data['Subjects_Data']
    LR_KFold_RandomCV(Subjects_Data, Subjects_Score, Covariates, Fold_Quantity,
                         # ComponentNumber_Range,
                         ResultantFolder, Parallel_Quantity, Permutation_Flag, RandIndex_File)



def LR_KFold_RandomCV(Subjects_Data, Subjects_Score, Covariates, Fold_Quantity,
                         ResultantFolder, Parallel_Quantity, Permutation_Flag,
                         RandIndex_File=''):
    if not os.path.exists(ResultantFolder):
        os.makedirs(ResultantFolder)
    # remove outliers
    mean_score = np.mean(Subjects_Score)
    std_score = np.std(Subjects_Score)
    inlier_idx = np.where(np.abs(Subjects_Score - mean_score) <= 3 * std_score)[0]
    Subjects_Score = Subjects_Score[inlier_idx]
    Subjects_Data = Subjects_Data[inlier_idx, :]
    Covariates = Covariates[inlier_idx, :]
    Subjects_Quantity = len(Subjects_Score)
    EachFold_Size = int(np.fix(np.divide(Subjects_Quantity, Fold_Quantity)))
    Remain = np.mod(Subjects_Quantity, Fold_Quantity)
    if len(RandIndex_File) == 0:
        RandIndex = np.arange(Subjects_Quantity)
        np.random.shuffle(RandIndex)
    else:
        tmpData = sio.loadmat(RandIndex_File[0])
        RandIndex = tmpData['RandIndex'][0]
    RandIndex_Mat = {'RandIndex': RandIndex}
    sio.savemat(ResultantFolder + '/RandIndex.mat', RandIndex_Mat)
    Fold_Corr = []
    Fold_MAE = []
    Fold_Weight = []
    Features_Quantity = np.shape(Subjects_Data)[1]
    for j in np.arange(Fold_Quantity):
        Fold_J_Index = RandIndex[EachFold_Size * j + np.arange(EachFold_Size)]
        if Remain > j:
            Fold_J_Index = np.insert(Fold_J_Index, len(Fold_J_Index), RandIndex[EachFold_Size * Fold_Quantity + j])
        Subjects_Data_test = Subjects_Data[Fold_J_Index, :]
        Subjects_Score_test = Subjects_Score[Fold_J_Index]
        Covariates_test = Covariates[Fold_J_Index, :]
        Subjects_Data_train = np.delete(Subjects_Data, Fold_J_Index, axis=0)
        Subjects_Score_train = np.delete(Subjects_Score, Fold_J_Index)
        Covariates_train = np.delete(Covariates, Fold_J_Index, axis=0)
        Covariates_Quantity = np.shape(Covariates)[1]
        # Controlling covariates from brain data
        df = {};
        df_test = {}
        for k in np.arange(Covariates_Quantity):
            df['Covariate_' + str(k)] = Covariates_train[:, k]
            df_test['Covariate_' + str(k)] = Covariates_test[:, k]
        # Construct formula
        Formula = 'Data ~ Covariate_0'
        for k in np.arange(Covariates_Quantity - 1) + 1:
            Formula = Formula + ' + Covariate_' + str(k)

        # Regress covariates from each brain features
        for k in np.arange(Features_Quantity):
            df['Data'] = Subjects_Data_train[:, k]
            df_test['Data'] = Subjects_Data_test[:, k]
            df_pd = pd.DataFrame(df)
            LinModel_Res = sm.ols(formula=Formula, data=df_pd).fit()
            Subjects_Data_train[:, k] = LinModel_Res.resid
            y_test_pred = LinModel_Res.predict(df_test)
            Subjects_Data_test[:, k] = Subjects_Data_test[:, k] - y_test_pred

        if Permutation_Flag:
            Subjects_Index_Random = np.arange(len(Subjects_Score_train))
            np.random.shuffle(Subjects_Index_Random)
            Subjects_Score_train = Subjects_Score_train[Subjects_Index_Random]
            if j == 0:
                PermutationIndex = {'Fold_0': Subjects_Index_Random}
            else:
                PermutationIndex['Fold_' + str(j)] = Subjects_Index_Random

        scaler = StandardScaler()
        Subjects_Data_train = scaler.fit_transform(Subjects_Data_train)
        Subjects_Data_test = scaler.transform(Subjects_Data_test)


        clf = linear_model.LinearRegression()
        clf.fit(Subjects_Data_train, Subjects_Score_train)
        Fold_J_Score = clf.predict(Subjects_Data_test)

        if Fold_J_Score.ndim > 1:
             Fold_J_Score = Fold_J_Score.ravel()

        print("********shape(fold, test)**********", Fold_J_Score.shape, Subjects_Score_test.shape)
        Fold_J_Corr = np.corrcoef(Fold_J_Score, Subjects_Score_test)
        Fold_J_Corr = Fold_J_Corr[0, 1]
        Fold_Corr.append(Fold_J_Corr)
        Fold_J_MAE = np.mean(np.abs(np.subtract(Fold_J_Score, Subjects_Score_test)))
        Fold_MAE.append(Fold_J_MAE)

        Weight = clf.coef_

        Fold_J_result = {'Index': Fold_J_Index, 'Test_Score': Subjects_Score_test, 'Predict_Score': Fold_J_Score,
                         'Corr': Fold_J_Corr, 'MAE': Fold_J_MAE, 'w_Brain': Weight}
        Fold_J_FileName = 'Fold_' + str(j) + '_Score.mat'
        ResultantFile = os.path.join(ResultantFolder, Fold_J_FileName)
        sio.savemat(ResultantFile, Fold_J_result)

    Fold_Corr = [0 if np.isnan(x) else x for x in Fold_Corr]
    Mean_Corr = np.mean(Fold_Corr)
    Mean_MAE = np.mean(Fold_MAE)

    Res_NFold = {'Mean_Corr': Mean_Corr, 'Mean_MAE': Mean_MAE};
    ResultantFile = os.path.join(ResultantFolder, 'Res_NFold.mat')
    sio.savemat(ResultantFile, Res_NFold)
    if Permutation_Flag:
        sio.savemat(ResultantFolder + '/PermutationIndex.mat', PermutationIndex)
    return (Mean_Corr, Mean_MAE)

