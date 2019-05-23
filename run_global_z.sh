#!/bin/bash -l
# created: 
# author: peterson
#SBATCH -J variant_block
#SBATCH --constraint="snb|hsw"
#SBATCH -o variant_block_output
#SBATCH -e variant_block_error
#SBATCH -p serial
#SBATCH -n 1
#SBATCH -t 40:00:00
#SBATCH --mem-per-cpu=130000
#SBATCH --mail-type=END
#SBATCH --mail-user=isaac.peterson@helsinki.fi
module load geo-env
/bin/sh  z_variants/z_variant_1.sh
/bin/sh  z_variants/z_variant_2.sh
/bin/sh  z_variants/z_variant_3.sh
/bin/sh  z_variants/z_variant_4.sh
/bin/sh  z_variants/z_variant_5.sh
/bin/sh  z_variants/z_variant_6.sh
