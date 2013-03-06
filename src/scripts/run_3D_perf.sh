#!/bin/bash
nsize_low=100
nsize_high=1000
#tstep=1000
set -x
for ((size=${nsize_low}; size <= ${nsize_high}; size += ${size})) do
	echo "perf -- $1 $size $size $size $size"
	./perf_3D.sh $1 $size $size $size $size
done
set +x

