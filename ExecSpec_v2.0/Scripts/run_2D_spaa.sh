#!/bin/sh

size_low=100
size_up=12800

for ((size = ${size_low}; size <= ${size_up}; size += ${size})) do
	echo "$1 $size 1000"
	$1 $size 1000
done
