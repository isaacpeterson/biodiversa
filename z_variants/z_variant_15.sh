#!/bin/sh
cd "$(dirname $0)"
zig4 -r z_variant_15/z_variant_15.dat z_variant_15/z_variant_15.spp z_variant_15/out/z_variant_15_out.txt 0.0 0 1.0 0 --grid-output-formats=compressed-tif --image-output-formats=png
