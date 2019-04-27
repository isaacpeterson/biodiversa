#!/bin/sh
cd "$(dirname $0)"
zig4 -r z_variant_19/z_variant_19.dat z_variant_19/z_variant_19.spp z_variant_19/out/z_variant_19_out.txt 0.0 0 1.0 0 --grid-output-formats=compressed-tif --image-output-formats=png
