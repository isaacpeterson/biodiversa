#!/bin/sh
cd "$(dirname $0)"
zig4 -r 03_amp_caz_wgt/03_amp_caz_wgt.dat 03_amp_caz_wgt/03_amp_caz_wgt.spp 03_amp_caz_wgt/03_amp_caz_wgt_out/03_amp_caz_wgt.txt 0.0 0 1.0 0 --grid-output-formats=compressed-tif --image-output-formats=png