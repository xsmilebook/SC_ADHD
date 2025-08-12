#!/bin/bash

#SBATCH --job-name=rbf_nback
#SBATCH --array=1-101          
#SBATCH --output=/ibmgpfs/cuizaixu_lab/xuhaoshu/SC_ADHD/datasets/logs/CV_%A_%a.out  
#SBATCH --error=/ibmgpfs/cuizaixu_lab/xuhaoshu/SC_ADHD/datasets/logs/CV_%A_%a.err   
#SBATCH --mem=8G                 
#SBATCH --cpus-per-task=1        
#SBATCH --partition=q_cn

#=======================================================================
# Script Execution
#=======================================================================


echo "================================================================"
echo "Starting Slurm Job: $SLURM_JOB_NAME"
echo "Slurm Job ID: $SLURM_JOB_ID"
echo "Slurm Array Job ID: $SLURM_ARRAY_JOB_ID"
echo "Slurm Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "================================================================"


source /GPFS/cuizaixu_lab_permanent/xuhaoshu/miniconda3/bin/activate ML

PYTHON_SCRIPT="predict_task_act.py"

TIME_ID=$SLURM_ARRAY_TASK_ID

python $PYTHON_SCRIPT --time_id $TIME_ID
# python $PYTHON_SCRIPT --permutation --time_id $TIME_ID

echo "================================================================"
echo "Job finished with exit code $?."
echo "================================================================"