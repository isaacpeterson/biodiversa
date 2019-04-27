#!/bin/sh
cd "$(dirname $0)"
zig4 -r 24_frf_abf_wgt_con_hm3/24_frf_abf_wgt_con_hm3.dat 24_frf_abf_wgt_con_hm3/24_frf_abf_wgt_con_hm3.spp 24_frf_abf_wgt_con_hm3/24_frf_abf_wgt_con_hm3_out/24_frf_abf_wgt_con_hm3.txt 0.0 0 1.0 0 --grid-output-formats=compressed-tif --image-output-formats=png