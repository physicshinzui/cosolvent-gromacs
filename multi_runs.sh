#!/bin/bash
set -eu

for i in $(seq 1 5); do
    jsub -v id=$i eq2prod.sh -N run_$i
    sleep 0.1
done
