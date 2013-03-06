#!/bin/bash

step1D_low=2000
step1D_high=32000
step1D_step=2000
step2D_low=500
step2D_high=500
step2D_step=200
step3D_low=100
step3D_high=1000
step3D_step=100
step4D_low=10
step4D_high=200
step4D_step=10
size1D_low=50000
size1D_high=400000
size1D_step=50000
size2D_low=200
size2D_high=3200
size2D_step=200
size3D_low=100
size3D_high=1000
size3D_step=100
size4D_low=20
size4D_high=400
size4D_step=20
size_rna_low=100
size_rna_high=1200
size_rna_step=100

set -x
#for ((size = ${size1D_low}; size <= ${size1D_high}; size += ${size1D_step})) do
#	echo "psa_struct_pochoir -r $size $size -i -d"
#	./psa_struct_pochoir -r $size $size -i -d
#done

#for ((step = ${step2D_low}; step <= ${step2D_high}; step += ${step})) do
#	for ((size = ${size2D_low}; size <= ${size2D_high}; size += ${size})) do
#		echo "heat_2D_NP_serial_loop $size $step"
#		./heat_2D_NP_sl $size $step
#	done
#done
#
#echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
#
#for ((step = ${step2D_low}; step <= ${step2D_high}; step += ${step})) do
#	for ((size = ${size2D_low}; size <= ${size2D_high}; size += ${size})) do
#		echo "heat_2D_NP_parallel_loop $size $step"
#		./heat_2D_NP_pl $size $step
#	done
#done
#
#echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
#
#for ((step = ${step2D_low}; step <= ${step2D_high}; step += ${step})) do
#	for ((size = ${size2D_low}; size <= ${size2D_high}; size += ${size})) do
#		echo "heat_2D_NP_pochoir $size $step"
#		./heat_2D_NP_pochoir $size $step
#	done
#done
#
#echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
#
#for ((step = ${step2D_low}; step <= ${step2D_high}; step += ${step})) do
#	for ((size = ${size2D_low}; size <= ${size2D_high}; size += ${size})) do
#		echo "heat_2D_P_serial_loop $size $step"
#		./heat_2D_P_sl $size $step
#	done
#done
#
#echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
#
#for ((step = ${step2D_low}; step <= ${step2D_high}; step += ${step})) do
#	for ((size = ${size2D_low}; size <= ${size2D_high}; size += ${size})) do
#		echo "heat_2D_P_parallel_loop $size $step"
#		./heat_2D_P_pl $size $step
#	done
#done
#
#echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
#
#for ((step = ${step2D_low}; step <= ${step2D_high}; step += ${step})) do
#	for ((size = ${size2D_low}; size <= ${size2D_high}; size += ${size})) do
#		echo "heat_2D_P_pochoir $size $step"
#		./heat_2D_P_pochoir $size $step
#	done
#done
#
#echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
#
#for ((step = ${step2D_low}; step <= ${step2D_high}; step += ${step})) do
#	for ((size = ${size2D_low}; size <= ${size2D_high}; size += ${size})) do
#		echo "life_serial_loop $size $step"
#		./life_sl $size $step
#	done
#done
#
#echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
#
#for ((step = ${step2D_low}; step <= ${step2D_high}; step += ${step})) do
#	for ((size = ${size2D_low}; size <= ${size2D_high}; size += ${size})) do
#		echo "life_parallel_loop $size $step"
#		./life_pl $size $step
#	done
#done
#
#echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
#
#for ((step = ${step2D_low}; step <= ${step2D_high}; step += ${step})) do
#	for ((size = ${size2D_low}; size <= ${size2D_high}; size += ${size})) do
#		echo "life_pochoir $size $step"
#		./life_pochoir $size $step
#	done
#done
#
#echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
#
#for ((step = ${step2D_low}; step <= ${step2D_high}; step += ${step})) do
#	for ((size = ${size2D_low}; size <= ${size2D_high}; size += ${size})) do
#		echo "life_pochoir_bit_trick $size $step"
#		./life_pochoir_bt $size $step
#	done
#done
#
#echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
#
#	for ((size = ${size1D_low}; size <= ${size1D_high}; size += ${size})) do
#		echo "lcs $size $step"
#		./lcs_pochoir -r $size $size -i -d
#	done
#
#echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
#
#	for ((size = ${size2D_low}; size <= ${size2D_high}; size += ${size})) do
#		echo "rna serial loop and pochoir $size $step"
#		./rna_pochoir_serial -r $size -i 
#	done
#
#echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
#
#	for ((size = ${size2D_low}; size <= ${size2D_high}; size += ${size})) do
#		echo "rna parallel loop and pochoir $size $step"
#		./rna_pochoir_parallel -r $size -i 
#	done
#
#echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

for ((step = ${step3D_low}; step <= ${step3D_high}; step += ${step3D_step})) do
	for ((size = ${size3D_low}; size <= ${size3D_high}; size += ${size3D_step})) do
		echo "tb_3dfd_pochoir $size $size $size $step"
		./tb_3dfd_pochoir $size $size $size $step
	done
done

set +x
