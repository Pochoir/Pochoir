#!/bin/bash

set -x

./run_3D_perf.sh ./3dfd_loop >& 3dfd_loop_perf.dat
./run_3D_perf.sh ./3dfd_seq >& 3dfd_seq_perf.dat
./run_3D_perf.sh ./3dfd_sim >& 3dfd_sim_perf.dat

./run_4D_perf.sh ./heat_4D_NP_loop >& heat_4D_NP_loop_perf.dat
./run_4D_perf.sh ./heat_4D_NP_seq >& heat_4D_NP_seq_perf.dat
./run_4D_perf.sh ./heat_4D_NP_sim >& heat_4D_NP_sim_perf.dat

set +x

