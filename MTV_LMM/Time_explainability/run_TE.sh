#!/bin/bash
set -ex

#mkdir -p ../results/allo_HSCT



# time # time a script
Rscript "TE_step_1.R"
./TE_step_2.sh $(cat init.txt)
Rscript "TE_step_3.R"




wait
