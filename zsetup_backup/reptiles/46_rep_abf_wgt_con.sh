#!/bin/sh
cd "$(dirname $0)"
zig4 -r 46_rep_abf_wgt_con/46_rep_abf_wgt_con.dat 46_rep_abf_wgt_con/46_rep_abf_wgt_con.spp 46_rep_abf_wgt_con/46_rep_abf_wgt_con_out/46_rep_abf_wgt_con.txt 0.0 0 1.0 0 --grid-output-formats=compressed-tif --image-output-formats=png