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
/bin/sh zsetup/amphibians/01_amp_caz.sh
/bin/sh zsetup/amphibians/02_amp_abf.sh
/bin/sh zsetup/amphibians/03_amp_caz_wgt.sh
/bin/sh zsetup/amphibians/04_amp_abf_wgt.sh
/bin/sh zsetup/amphibians/05_amp_caz_wgt_con.sh
/bin/sh zsetup/amphibians/06_amp_abf_wgt_con.sh
/bin/sh zsetup/amphibians/07_amp_caz_wgt_con_hm3.sh
/bin/sh zsetup/amphibians/08_amp_abf_wgt_con_hm3.sh
/bin/sh zsetup/birds/09_bir_caz.sh
/bin/sh zsetup/birds/10_bir_abf.sh
/bin/sh zsetup/birds/11_bir_caz_wgt.sh
/bin/sh zsetup/birds/12_bir_abf_wgt.sh
/bin/sh zsetup/birds/13_bir_caz_wgt_con.sh
/bin/sh zsetup/birds/14_bir_abf_wgt_con.sh
/bin/sh zsetup/birds/15_bir_caz_wgt_con_hm3.sh
/bin/sh zsetup/birds/16_bir_abf_wgt_con_hm3.sh
/bin/sh zsetup/freshwater_fish/17_frf_caz.sh
/bin/sh zsetup/freshwater_fish/18_frf_abf.sh
/bin/sh zsetup/freshwater_fish/19_frf_caz_wgt.sh
/bin/sh zsetup/freshwater_fish/20_frf_abf_wgt.sh
/bin/sh zsetup/freshwater_fish/21_frf_caz_wgt_con.sh
/bin/sh zsetup/freshwater_fish/22_frf_abf_wgt_con.sh
/bin/sh zsetup/freshwater_fish/23_frf_caz_wgt_con_hm3.sh
/bin/sh zsetup/freshwater_fish/24_frf_abf_wgt_con_hm3.sh
/bin/sh zsetup/mammals/25_mam_caz.sh
/bin/sh zsetup/mammals/26_mam_abf.sh
/bin/sh zsetup/mammals/27_mam_caz_wgt.sh
/bin/sh zsetup/mammals/28_mam_abf_wgt.sh
/bin/sh zsetup/mammals/29_mam_caz_wgt_con.sh
/bin/sh zsetup/mammals/30_mam_abf_wgt_con.sh
/bin/sh zsetup/mammals/31_mam_caz_wgt_con_hm3.sh
/bin/sh zsetup/mammals/32_mam_abf_wgt_con_hm3.sh
/bin/sh zsetup/plants/33_pla_caz.sh
/bin/sh zsetup/plants/34_pla_abf.sh
/bin/sh zsetup/plants/35_pla_caz_wgt.sh
/bin/sh zsetup/plants/36_pla_abf_wgt.sh
/bin/sh zsetup/plants/37_pla_caz_wgt_con.sh
/bin/sh zsetup/plants/38_pla_abf_wgt_con.sh
/bin/sh zsetup/plants/39_pla_caz_wgt_con_hm3.sh
/bin/sh zsetup/plants/40_pla_abf_wgt_con_hm3.sh
/bin/sh zsetup/reptiles/41_rep_caz.sh
/bin/sh zsetup/reptiles/42_rep_abf.sh
/bin/sh zsetup/reptiles/43_rep_caz_wgt.sh
/bin/sh zsetup/reptiles/44_rep_abf_wgt.sh
/bin/sh zsetup/reptiles/45_rep_caz_wgt_con.sh
/bin/sh zsetup/reptiles/46_rep_abf_wgt_con.sh
/bin/sh zsetup/reptiles/47_rep_caz_wgt_con_hm3.sh
/bin/sh zsetup/reptiles/48_rep_abf_wgt_con_hm3.sh
/bin/sh zsetup/taxa_all/49_caz.sh
/bin/sh zsetup/taxa_all/50_abf.sh
/bin/sh zsetup/taxa_all/51_caz_wgt.sh
/bin/sh zsetup/taxa_all/52_abf_wgt.sh
/bin/sh zsetup/taxa_all/53_caz_wgt_con.sh
/bin/sh zsetup/taxa_all/54_abf_wgt_con.sh
/bin/sh zsetup/taxa_all/55_caz_wgt_con_hm3.sh
/bin/sh zsetup/taxa_all/56_abf_wgt_con_hm3.sh
