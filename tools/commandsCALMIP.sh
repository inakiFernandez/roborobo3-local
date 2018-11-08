#!/bin/bash
#SBATCH -J script_art633
#SBATCH -N 4
#SBATCH -n 128
#SBATCH --ntasks-per-node=36
#SBATCH --ntasks-per-core=1
#SBATCH --time=00:03:00
#SBATCH --mail-user=fernandezperez.inaki@gmail.com
#SBATCH --array=1-256
srun ./commandeALaCon.sh $SLURM_ARRAY_TASK_ID

