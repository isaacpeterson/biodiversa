#!/bin/sh
cd "$(dirname $0)"
zig4 -r 16_bir_abf_wgt_con_hm3/16_bir_abf_wgt_con_hm3.dat 16_bir_abf_wgt_con_hm3/16_bir_abf_wgt_con_hm3.spp 16_bir_abf_wgt_con_hm3/16_bir_abf_wgt_con_hm3_out/16_bir_abf_wgt_con_hm3.txt 0.0 0 1.0 0 --grid-output-formats=compressed-tif --image-output-formats=png
