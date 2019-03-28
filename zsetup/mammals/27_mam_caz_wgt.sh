#!/bin/sh
cd "$(dirname $0)"
zig4 -r 27_mam_caz_wgt/27_mam_caz_wgt.dat 27_mam_caz_wgt/27_mam_caz_wgt.spp 27_mam_caz_wgt/27_mam_caz_wgt_out/27_mam_caz_wgt.txt 0.0 0 1.0 0 --grid-output-formats=compressed-tif --image-output-formats=png
