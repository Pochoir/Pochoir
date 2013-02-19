#!/bin/bash

inFiles="`ls *.hs`"

for i in ${inFiles}; do
    prefix=${i:0:1}
    rest=${i:1}
    dest="${prefix}Gen${rest}"
    echo "prefix = $prefix, rest = $rest, dest = $dest"
    echo "cp $i $dest"
    cp $i $dest
done

