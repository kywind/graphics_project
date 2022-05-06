#!/usr/bin/env bash

# If project not ready, generate cmake file.
if [[ ! -d build ]]; then
    echo "good"
else
    rm -rf build
fi
mkdir -p build
cd build
cmake ..
make -j
cd ..

# Run all testcases. 
# You can comment some lines to disable the run of specific examples.

mkdir -p output
bin/FinalProject testcases/demo.txt output/demo.bmp path
# bin/FinalProject testcases/scene01_path_tracing.txt output/scene01_path_tracing.bmp path
# bin/FinalProject testcases/scene02_path_tracing.txt output/scene02_path_tracing.bmp path
# bin/FinalProject testcases/scene03_ray_tracing.txt output/scene03_ray_tracing.bmp ray