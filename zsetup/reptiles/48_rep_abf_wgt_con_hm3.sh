#!/bin/sh
cd "$(dirname $0)"
zig4 -r 48_rep_abf_wgt_con_hm3/48_rep_abf_wgt_con_hm3.dat 48_rep_abf_wgt_con_hm3/48_rep_abf_wgt_con_hm3.spp 48_rep_abf_wgt_con_hm3/48_rep_abf_wgt_con_hm3_out/48_rep_abf_wgt_con_hm3.txt 0.0 0 1.0 0 --grid-output-formats=compressed-tif --image-output-formats=png
