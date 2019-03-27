#!/bin/sh
zig4 -r 01_amp_caz/01_amp_caz.dat 01_amp_caz/01_amp_caz.spp 01_amp_caz/01_amp_caz_out/01_amp_caz.txt 0.0 0 1.0 0 --grid-output-formats=compressed-tif --image-output-formats=png
