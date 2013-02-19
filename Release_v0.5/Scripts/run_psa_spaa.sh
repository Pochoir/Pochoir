#!/bin/sh

size_low=100
size_up=1280000

for ((size = ${size_low}; size <= ${size_up}; size += ${size})) do
	echo "$1 -r $size $size"
	$1 -r $size $size
done
