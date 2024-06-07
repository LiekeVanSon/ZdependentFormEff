#!/bin/bash

export DISBATCH_KVSSTCP_HOST=10.128.146.28:37939 PYTHONPATH=/mnt/home/carriero/projects/disBatch/beta/disBatch:${PYTHONPATH}

if [[ $1 == '--mon' ]]
then
    exec /mnt/home/lvanson/venvs/Jupyter_venv/bin/python3 /mnt/home/carriero/projects/disBatch/beta/disBatch/disbatchc/dbMon.py /mnt/home/lvanson/ZdependentFormEff/code/Single_StarTasks_disBatch_240607145120_473
elif [[ $1 == '--engine' ]]
then
    exec /mnt/home/lvanson/venvs/Jupyter_venv/bin/python3 /mnt/home/carriero/projects/disBatch/beta/disBatch/disBatch "$@"
else
    exec /mnt/home/lvanson/venvs/Jupyter_venv/bin/python3 /mnt/home/carriero/projects/disBatch/beta/disBatch/disBatch --context /mnt/home/lvanson/ZdependentFormEff/code/Single_StarTasks_disBatch_240607145120_473_dbUtil.sh "$@" < /dev/null 1> /mnt/home/lvanson/ZdependentFormEff/code/Single_StarTasks_disBatch_240607145120_473_${BASHPID-$$}_context_launch.log
fi
