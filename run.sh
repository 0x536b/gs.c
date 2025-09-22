#!/bin/bash

set -ex

echo "== delete files.."
rm -rf render_gs_c.ppm
rm -rf render_gs_c.png
rm -rf render

echo "== compiling and running"
gcc gs.c -o render -lm
./render

echo "== converting to png"
magick render_gs_c.ppm render_gs_c.png