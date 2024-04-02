#!/bin/bash

NSITES="40 80"
NSPECIES="26 52"
NREPS=1000

parallel --jobs 128 --keep-order --plus --progress "
    mkdir -p job-{0#} && cp scenario.qmd job-{0#} &&
    quarto render job-{0#}/scenario.qmd -P nsites:{1} -P nspecies:{2} -P job:{#} &&
    rm job-{0#}/scenario.qmd" ::: $NSITES ::: $NSPECIES ::: $(seq $NREPS)
