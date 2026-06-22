#!/bin/bash
# Build the standalone BMW Wilson-flow + w0 tools (no MPI, no gwu-qcd).
cd "$(dirname "$0")"
gcc -O3 -std=c99 -fopenmp -o wilson_flow wilson_flow.c -lm
gcc -O3 -o w0_scale w0_scale.c -lm
echo "built wilson_flow, w0_scale"
