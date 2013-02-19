#!/bin/bash

srcFile=$1
dstFile="$1.orig"

set -x
cp $srcFile $dstFile
sed -f sed.sh $dstFile >& $srcFile
rm $dstFile
set +x
