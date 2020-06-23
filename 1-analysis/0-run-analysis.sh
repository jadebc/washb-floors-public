#!/bin/bash

R CMD BATCH 1-washb-floors-mde.R &
R CMD BATCH 2a-washb-prev-mean-bd.R & 
R CMD BATCH 2b-washb-prev-mean-ke.R &
R CMD BATCH 3-washb-floors-pfloors.R &
R CMD BATCH 4-washb-analysis.R &
R CMD BATCH 5-washb-analysis-intensity.R &
R CMD BATCH 6-washb-analysis-effectmod.R &
R CMD BATCH 7-washb-analysis-pos.R &
R CMD BATCH 8-e-values.R & 
R CMD BATCH 9-demographics.R & 