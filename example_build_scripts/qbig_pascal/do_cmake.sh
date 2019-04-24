#!/bin/bash
peram_gen_src_dir=..

CXX=/opt/openmpi-2.0.2a1-with-pmi/bin/mpicxx \
CC=/opt/openmpi-2.0.2a1-with-pmi/bin/mpicc \
CXXFLAGS="-fopenmp -O3 -mtune=broadwell -march=broadwell -g" \
CFLAGS="-fopenmp -O3 -mtune=broadwell -march=broadwell -g" \
LDFLAGS="-fopenmp" \
cmake \
-DTMLQCD_SRC=/hadron/bartek/code/bleeding_edge/tmLQCD.quda_work \
-DTMLQCD_BUILD=/hadron/bartek/build/bleeding_edge/pascal/tmLQCD.quda_work.quda_develop-dynamic_clover \
-DLIME_HOME=/qbigwork2/bartek/libs/lime \
-DLEMON_HOME=/qbigwork2/bartek/libs/lemon/broadwell \
-DQUDA_HOME=/qbigwork2/bartek/libs/bleeding_edge/pascal/quda_develop-dynamic_clover \
-DSUPPORT_QUDA_DIRECT=ON \
${peram_gen_src_dir}

