#!/bin/bash
nsize_low=10
nsize_high=160
tstep=1000
set -x
for ((size=${nsize_low}; size <= ${nsize_high}; size += ${size})) do
	echo "perf -- $1 $size $tstep "
	./perf_4D.sh $1 $size $tstep
done
set +x

