#!/bin/bash
#SBATCH -J calmiptestInaki
#SBATCH --job-name=test_emb_arr
#SBATCH --output=res_emb_arr.txt
#
#SBATCH -N 5
#SBATCH -n 20
#SBATCH --mail-user=fernandezperez.inaki@gmail.com
#SBATCH --mail-type=ALL


parallel -j80 -a ./config/evolvability.properties.parallel
#$SLURM_JOB_ID

