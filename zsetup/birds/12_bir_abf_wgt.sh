#!/bin/sh
cd "$(dirname $0)"
zig4 -r 12_bir_abf_wgt/12_bir_abf_wgt.dat 12_bir_abf_wgt/12_bir_abf_wgt.spp 12_bir_abf_wgt/12_bir_abf_wgt_out/12_bir_abf_wgt.txt 0.0 0 1.0 0 --grid-output-formats=compressed-tif --image-output-formats=png
