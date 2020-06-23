#!/bin/bash

R CMD BATCH 1-figtab-prev.R &
R CMD BATCH 2-figtab-box-plots.R &
R CMD BATCH 3-fig-box-plots.R &
R CMD BATCH 4-fig-results-effectmod-qpcr.R & 
R CMD BATCH 5-fig-sens.R &
R CMD BATCH 6-fig-pscores.R & 
