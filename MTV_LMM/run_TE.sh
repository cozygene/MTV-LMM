#!/bin/bash
set -ex



Rscript "Time_explainability/TE_step_1.R" $2 $3
chmod +x Time_explainability/TE_step_2.sh
Time_explainability/TE_step_2.sh $1 
Rscript "Time_explainability/TE_step_3.R" $2 $3



