#!/bin/sh
cd "$(dirname $0)"
zig4 -r 44_rep_abf_wgt/44_rep_abf_wgt.dat 44_rep_abf_wgt/44_rep_abf_wgt.spp 44_rep_abf_wgt/44_rep_abf_wgt_out/44_rep_abf_wgt.txt 0.0 0 1.0 0 --grid-output-formats=compressed-tif --image-output-formats=png
