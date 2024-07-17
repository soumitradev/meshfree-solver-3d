#!/bin/sh
rm ./build/meshfree_solver_test.out
mkdir -p build
apptainer exec --nv --env LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/usr/local/cuda-12/compat:/legion/bindings/regent --env CPATH=\$CPATH:/usr/local/cuda-12/targets/x86_64-linux/include/ --bind ./src:/src --bind ./data:/data --bind ./build:/build ../regent_env/v24.03.0/v24.03.0.sif /legion/language/regent.py /src/meshfree_solver.rg