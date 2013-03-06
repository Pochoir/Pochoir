#!/bin/sh

size_low=0
size_up=1000

for ((size = ${size_low}; size <= ${size_up}; size += 200)) do
	echo "$1 $size $size $size $size" 
	$1 $size $size $size $size
done
