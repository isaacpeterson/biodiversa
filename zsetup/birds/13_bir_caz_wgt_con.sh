#!/bin/sh
cd "$(dirname $0)"
zig4 -r 13_bir_caz_wgt_con/13_bir_caz_wgt_con.dat 13_bir_caz_wgt_con/13_bir_caz_wgt_con.spp 13_bir_caz_wgt_con/13_bir_caz_wgt_con_out/13_bir_caz_wgt_con.txt 0.0 0 1.0 0 --grid-output-formats=compressed-tif --image-output-formats=png
