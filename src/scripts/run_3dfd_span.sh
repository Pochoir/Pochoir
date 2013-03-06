#!/bin/bash
#set -x
nsize_low=100
nsize_high=1000
tstep=1000
for ((size=${nsize_low}; size <= ${nsize_high}; size += ${size})) do
	echo "cilkview -- $1 $size $size $size $size"
	cilkview $1 $size $size $size $size
done
#set +x

