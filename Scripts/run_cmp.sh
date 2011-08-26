#!/bin/sh

set -x

./run_2D_spaa.sh ./tb_heat_2D_P_raw_spaa >& heat_2D_P_raw.dat
./run_2D_spaa.sh ./tb_heat_2D_P_macro_spaa >& heat_2D_P_macro.dat

./run_2D_spaa.sh ./tb_heat_2D_P_c_spaa >& heat_2D_P_c.dat
./run_2D_spaa.sh ./tb_heat_2D_P_opt_spaa >& heat_2D_P_opt.dat
./run_2D_spaa.sh ./tb_heat_2D_P_pt_spaa >& heat_2D_P_pt.dat
set +x
