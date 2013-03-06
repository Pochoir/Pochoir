#!/bin/bash

set -x
cilkview ./heat_duo 2000 2000 >& duo.scal

cilkview ./heat_sim 2000 2000 >& sim.scal
set +x

