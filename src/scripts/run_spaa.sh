#!/bin/bash

set -x
./run_2D_spaa.sh ./tb_heat_2D_NP_spaa >& heat_2D_NP_spaa.dat
./run_2D_spaa.sh ./tb_heat_2D_P_spaa >& heat_2D_P_spaa.dat
./run_2D_spaa.sh ./tb_life_spaa >& life_spaa.dat
./run_rna_spaa.sh ./rna_spaa >& rna_spaa.dat
./run_psa_spaa.sh ./psa_spaa >& psa_spaa.dat
./run_psa_spaa.sh ./lcs_spaa >& lcs_spaa.dat
./run_3D_spaa.sh ./tb_3dfd_spaa >& 3dfd_spaa.dat
./run_4D_spaa.sh ./tb_heat_4D_NP_spaa >& heat_4D_NP_spaa.dat
set +x

