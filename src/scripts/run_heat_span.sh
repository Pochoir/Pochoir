#!/bin/bash
nsize_low=100
nsize_high=6400
tstep=1000
set -x
for ((size=${nsize_low}; size <= ${nsize_high}; size += ${size})) do
	echo "cilkview -- $1 $size $tstep"
	cilkview $1 $size $tstep
done
set +x

