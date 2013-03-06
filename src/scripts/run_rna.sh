#!/bin/bash

set -x
perf stat -e L1-dcache-loads -e L1-dcache-load-misses -e L1-dcache-stores -e L1-dcache-store-misses -e L1-dcache-prefetches -e L1-dcache-prefetch-misses -e L1-icache-load -e L1-icache-load-misses -e L1-icache-prefetches -e L1-icache-prefetch-misses -e LLC-loads -e LLC-load-misses -e LLC-stores -e LLC-store-misses -e LLC-prefetches -e LLC-prefetch-misses -e cache-references -e cache-misses -e branches -e branch-misses -e instructions -- ./rna -r 300 >& rna_pochoir.perf

perf stat -e L1-dcache-loads -e L1-dcache-load-misses -e L1-dcache-stores -e L1-dcache-store-misses -e L1-dcache-prefetches -e L1-dcache-prefetch-misses -e L1-icache-load -e L1-icache-load-misses -e L1-icache-prefetches -e L1-icache-prefetch-misses -e LLC-loads -e LLC-load-misses -e LLC-stores -e LLC-store-misses -e LLC-prefetches -e LLC-prefetch-misses -e cache-references -e cache-misses -e branches -e branch-misses -e instructions -- ./rna -r 300 -i >& rna_iter.perf

cilkview ./rna -r 300 >& rna_pochoir.span

cilkview ./rna -r 300 -i >& rna_iter.span
set +x
