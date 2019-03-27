#!/bin/bash -l
# created: 
# author: peterson
#SBATCH -J variant_101
#SBATCH --constraint="snb|hsw"
#SBATCH -o variant101_output
#SBATCH -e variant101_error
#SBATCH -p serial
#SBATCH -n 1
#SBATCH -t 40:00:00
#SBATCH --mem-per-cpu=130000
#SBATCH --mail-type=END
#SBATCH --mail-user=isaac.peterson@helsinki.fi

# example run commands
module load geo-env 