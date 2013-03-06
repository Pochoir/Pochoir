#!/bin/bash

set -x

./run_cmp.sh
./run_perf.sh
./run_scal.sh

set +x
