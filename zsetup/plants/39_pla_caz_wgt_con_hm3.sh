#!/bin/sh
cd "$(dirname $0)"
zig4 -r 39_pla_caz_wgt_con_hm3/39_pla_caz_wgt_con_hm3.dat 39_pla_caz_wgt_con_hm3/39_pla_caz_wgt_con_hm3.spp 39_pla_caz_wgt_con_hm3/39_pla_caz_wgt_con_hm3_out/39_pla_caz_wgt_con_hm3.txt 0.0 0 1.0 0 --grid-output-formats=compressed-tif --image-output-formats=png
