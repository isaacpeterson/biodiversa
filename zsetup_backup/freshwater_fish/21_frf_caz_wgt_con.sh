#!/bin/sh
cd "$(dirname $0)"
zig4 -r 21_frf_caz_wgt_con/21_frf_caz_wgt_con.dat 21_frf_caz_wgt_con/21_frf_caz_wgt_con.spp 21_frf_caz_wgt_con/21_frf_caz_wgt_con_out/21_frf_caz_wgt_con.txt 0.0 0 1.0 0 --grid-output-formats=compressed-tif --image-output-formats=png