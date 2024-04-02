#!/bin/bash

TAXA="col sma veg"
DISTRIBUTIONS="negbin poisson yule"

parallel --jobs 12 --progress "
    mkdir -p {1}-{2} && cp rtg-2300.qmd {1}-{2} &&
    quarto render {1}-{2}/rtg-2300.qmd -P taxon:{1} -P distribution:{2} &&
    rm {1}-{2}/rtg-2300.qmd" ::: $TAXA ::: $DISTRIBUTIONS
