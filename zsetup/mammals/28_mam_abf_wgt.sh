#!/bin/sh
cd "$(dirname $0)"
zig4 -r 28_mam_abf_wgt/28_mam_abf_wgt.dat 28_mam_abf_wgt/28_mam_abf_wgt.spp 28_mam_abf_wgt/28_mam_abf_wgt_out/28_mam_abf_wgt.txt 0.0 0 1.0 0 --grid-output-formats=compressed-tif --image-output-formats=png
