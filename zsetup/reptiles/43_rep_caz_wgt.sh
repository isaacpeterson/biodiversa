#!/bin/sh
cd "$(dirname $0)"
zig4 -r 43_rep_caz_wgt/43_rep_caz_wgt.dat 43_rep_caz_wgt/43_rep_caz_wgt.spp 43_rep_caz_wgt/43_rep_caz_wgt_out/43_rep_caz_wgt.txt 0.0 0 1.0 0 --grid-output-formats=compressed-tif --image-output-formats=png
