#!/bin/bash
set -ex

#mkdir -p ../results/allo_HSCT



# time # time a script
Rscript "Time_explainability/TE_step_1.R"
Time_explainability/TE_step_2.sh $(cat init_example.txt)
Rscript "Time_explainability/TE_step_3.R"




wait
