#!/bin/bash
#SBATCH -J calmiptestInaki
#SBATCH --job-name=test_emb_arr
#SBATCH --output=res_emb_arr.txt
#
#SBATCH -N 15
#SBATCH -n 20
#SBATCH --mail-user=fernandezperez.inaki@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --time=24:00:00

parallel -j300 -a ./config/evolvability.properties.parallel

#$SLURM_JOB_ID

