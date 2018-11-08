#!/bin/bash
#SBATCH -J script_art633
#SBATCH --job-name=test_emb_arr
#SBATCH --output=res_emb_arr.txt
#
#SBATCH -c 10
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks-per-core=1
#SBATCH --time=00:03:00
#SBATCH --mail-user=fernandezperez.inaki@gmail.com
#SBATCH --array=1-256
#SBATCH --mem-per-cpu=100

srun ./commandeALaCon.sh $SLURM_ARRAY_TASK_ID

