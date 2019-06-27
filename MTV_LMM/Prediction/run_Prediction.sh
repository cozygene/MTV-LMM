#!/bin/bash
set -ex

#mkdir -p ../results/allo_HSCT



# time # time a script
#Rscript "/Prediction/Prediction_step_1.R"
/Prediction/Prediction_step_2.sh $(cat init.txt)
Rscript "/Prediction/Prediction_step_3.R"




wait
