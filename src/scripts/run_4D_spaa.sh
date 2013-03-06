#!/bin/sh

size_low=10
size_up=160

for ((size = ${size_low}; size <= ${size_up}; size += $size)) do
	echo "$1 $size $size" 
	$1 $size $size
done
