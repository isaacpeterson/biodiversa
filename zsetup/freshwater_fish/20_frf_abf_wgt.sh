#!/bin/sh
cd "$(dirname $0)"
zig4 -r 20_frf_abf_wgt/20_frf_abf_wgt.dat 20_frf_abf_wgt/20_frf_abf_wgt.spp 20_frf_abf_wgt/20_frf_abf_wgt_out/20_frf_abf_wgt.txt 0.0 0 1.0 0 --grid-output-formats=compressed-tif --image-output-formats=png
