#!/bin/sh
zig4 -r 04_amp_abf_wgt/04_amp_abf_wgt.dat 04_amp_abf_wgt/04_amp_abf_wgt.spp 04_amp_abf_wgt/04_amp_abf_wgt_out/04_amp_abf_wgt.txt 0.0 0 1.0 0 --grid-output-formats=compressed-tif --image-output-formats=png
