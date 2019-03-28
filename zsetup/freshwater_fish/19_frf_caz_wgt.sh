#!/bin/sh
cd "$(dirname $0)"
zig4 -r 19_frf_caz_wgt/19_frf_caz_wgt.dat 19_frf_caz_wgt/19_frf_caz_wgt.spp 19_frf_caz_wgt/19_frf_caz_wgt_out/19_frf_caz_wgt.txt 0.0 0 1.0 0 --grid-output-formats=compressed-tif --image-output-formats=png
