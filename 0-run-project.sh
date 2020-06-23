#!/bin/bash

cd 1-analysis
bash 0-run-analysis.sh 

cd ..
cd 2-figure-scripts
bash 0-run-figures.sh  

cd .. 
cd 3-table-scripts
bash 0-run-tables.sh  
