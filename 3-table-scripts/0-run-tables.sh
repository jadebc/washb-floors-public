#!/bin/bash

R CMD BATCH 1a-table-characteristics-bd.R &
R CMD BATCH 1b-table-characteristics-ke.R &
R CMD BATCH 1c-table-characteristics-full.R &
R CMD BATCH 2a-table-prev-cq-qpcr.R &
R CMD BATCH 2b-table-prev-mean-kk.R &
R CMD BATCH 3a-table-characteristics-pos-bd.R &
R CMD BATCH 3b-table-characteristics-pos-ke.R &
R CMD BATCH 3c-table-characteristics-pos-full.R &
R CMD BATCH 4-table-prim-vs-sens-adj-unadj-pr.R &
R CMD BATCH 5-table-pr-tmle.R & 
R CMD BATCH 6-table-pr-glm-kk.R & 
R CMD BATCH 7-table-fecr-kk.R &
R CMD BATCH 8-table-fecr-qpcr-tmle.R & 
R CMD BATCH 9-table-evalues.R & 
R CMD BATCH 10-table-kk-intensity.R & 