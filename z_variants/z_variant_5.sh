#!/bin/sh
cd "$(dirname $0)"
zig4 -r z_variant_5/z_variant_5.dat z_variant_5/z_variant_5.spp z_variant_5/out/z_variant_5_out.txt 0.0 0 1.0 0 --grid-output-formats=compressed-tif --image-output-formats=png
