#!/bin/sh
cd "$(dirname $0)"
zig4 -r 47_rep_caz_wgt_con_hm3/47_rep_caz_wgt_con_hm3.dat 47_rep_caz_wgt_con_hm3/47_rep_caz_wgt_con_hm3.spp 47_rep_caz_wgt_con_hm3/47_rep_caz_wgt_con_hm3_out/47_rep_caz_wgt_con_hm3.txt 0.0 0 1.0 0 --grid-output-formats=compressed-tif --image-output-formats=png
