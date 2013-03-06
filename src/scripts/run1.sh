#!/bin/bash
#set -x
nsize_low=200
nsize_high=6400
tstep=1000
for ((size=${nsize_low}; size <= ${nsize_high}; size += ${size})) do
	echo "tb_heat_2D_P_macro $size $tstep"
	./tb_heat_2D_P_macro $size $tstep
	echo "tb_heat_2D_P_pointer $size $tstep"
	./tb_heat_2D_P_pointer $size $tstep
	echo "tb_heat_2D_P_opt_pointer $size $tstep"
	./tb_heat_2D_P_opt_pointer $size $tstep
done
#set +x

