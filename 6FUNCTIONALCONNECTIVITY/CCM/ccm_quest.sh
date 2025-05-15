#!/bin/bash

#SBATCH --account=p32032
#SBATCH --partition=long
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=64
#SBATCH --time=168:00:00
#SBATCH --mem=243g
#SBATCH --job-name=diya_recalculate

module load matlab/r2022b
matlab -singleCompThread -batch "ccm_qscript"
