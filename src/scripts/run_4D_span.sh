#!/bin/bash
#set -x
nsize_low=10
nsize_high=160
tstep=1000
for ((size=${nsize_low}; size <= ${nsize_high}; size += ${size})) do
	echo "cilkview -- $1 $size $tstep"
	cilkview $1 $size $tstep
done
#set +x

