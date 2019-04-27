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
module load geo-env
/bin/sh  z_variants/z_variant_1.sh
/bin/sh  z_variants/z_variant_2.sh
/bin/sh  z_variants/z_variant_3.sh
/bin/sh  z_variants/z_variant_4.sh
/bin/sh  z_variants/z_variant_5.sh
/bin/sh  z_variants/z_variant_6.sh
/bin/sh  z_variants/z_variant_7.sh
/bin/sh  z_variants/z_variant_8.sh
/bin/sh  z_variants/z_variant_9.sh
/bin/sh  z_variants/z_variant_10.sh
/bin/sh  z_variants/z_variant_11.sh
/bin/sh  z_variants/z_variant_12.sh
/bin/sh  z_variants/z_variant_13.sh
/bin/sh  z_variants/z_variant_14.sh
/bin/sh  z_variants/z_variant_15.sh
/bin/sh  z_variants/z_variant_16.sh
/bin/sh  z_variants/z_variant_17.sh
/bin/sh  z_variants/z_variant_18.sh
/bin/sh  z_variants/z_variant_19.sh
/bin/sh  z_variants/z_variant_20.sh
