#!/bin/bash

set -x
##echo "./run_3dfd_span.sh ./3dfd_loop >& 3dfd_loop_span.dat"
##./run_3dfd_span.sh ./3dfd_loop >& 3dfd_loop_span.dat
##
##echo "./run_3dfd_span.sh ./3dfd_seq >& 3dfd_seq_span.dat"
##./run_3dfd_span.sh ./3dfd_seq >& 3dfd_seq_span.dat
##
##echo "./run_3dfd_span.sh ./3dfd_sim >& 3dfd_sim_span.dat"
##./run_3dfd_span.sh ./3dfd_sim >& 3dfd_sim_span.dat
##
##echo "./run_4D_span.sh ./heat_4D_NP_loop >& heat_4D_NP_loop_span.dat"
##./run_4D_span.sh ./heat_4D_NP_loop >& heat_4D_NP_loop_span.dat

##echo "./run_4D_span.sh ./heat_4D_NP_seq >& heat_4D_NP_seq_span.dat"
##./run_4D_span.sh ./heat_4D_NP_seq >& heat_4D_NP_seq_span.dat
##
##echo "./run_4D_span.sh ./heat_4D_NP_sim >& heat_4D_NP_sim_span.dat"
##./run_4D_span.sh ./heat_4D_NP_sim >& heat_4D_NP_sim_span.dat

echo "./run_heat_span.sh ./heat_sim >& heat_sim_span.dat"
./run_heat_span.sh ./heat_sim >& heat_sim_span.dat

echo "./run_heat_span.sh ./heat_duo >& heat_duo_span.dat"
./run_heat_span.sh ./heat_duo >& heat_duo_span.dat
##
##echo "./run_heat_span.sh ./heat_2D_NP_sim >& heat_2D_NP_sim_span.dat"
##./run_heat_span.sh ./heat_2D_NP_sim >& heat_2D_NP_sim_span.dat
##
##echo "./run_heat_span.sh ./heat_2D_P_loop >& heat_2D_P_loop_span.dat"
##./run_heat_span.sh ./heat_2D_P_loop >& heat_2D_P_loop_span.dat
##
##echo "./run_heat_span.sh ./heat_2D_P_seq >& heat_2D_P_seq_span.dat"
##./run_heat_span.sh ./heat_2D_P_seq >& heat_2D_P_seq_span.dat
##
##echo "./run_heat_span.sh ./heat_2D_P_sim >& heat_2D_P_sim_span.dat"
##./run_heat_span.sh ./heat_2D_P_sim >& heat_2D_P_sim_span.dat

set +x
