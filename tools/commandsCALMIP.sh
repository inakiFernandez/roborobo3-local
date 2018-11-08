#!/bin/bash
#SBATCH -J calmiptestInaki
#SBATCH --job-name=test_emb_arr
#SBATCH --output=res_emb_arr.txt
#
#SBATCH -N 4
#SBATCH -n 20
#SBATCH --time=00:03:00
#SBATCH --mail-user=fernandezperez.inaki@gmail.com
#SBATCH --mail-type=ALL


parallel -j80 ./commandeALaCon.sh $SLURM_ARRAY_TASK_ID

