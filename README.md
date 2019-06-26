# MTV-LMM

MTV-LMM is a linear mixed model method for the prediction of the microbial community temporal dynamics based on the community composition at previous time stamps. MTV-LMM thus allows identifying auto-regressive OTUs in time series microbiome datasets, which can be taken to further analysis of the microbiome trajectory over time.

Software Requirements:
-----------------------

- R version 3.4.0 and above 
- [GCTA_1.26 and above](https://cnsgenomics.com/software/gcta/#Download)


Content
-----------------------
The MTV-LMM package contains:
- otu_table_example.csv : A time-series OTU table (7244 OTUs x 1101 Time points, 40 individuals) - example dataset.
- metadata_example.csv : A meta-data file - example dataset. 
- src.R : Source code for MTV_LMM.
- Main_MTV_LMM_TE.R : Calculating Time_explainability for the example dataset 'otu_table_example.csv'.  
- Main_MTV_LMM_Prediction.R : Generating the predictions for the example dataset 'otu_table_example.csv'.
- run_mgrm.sh : Shell code, running GCTA software as part of the Time_explainability calculation.
- run_predictions.sh : Shell code, running GCTA software  as part of the prediction generation.


Identifying auto-regressive OTUs - calculating 'Time_explainability'
--------------------------
#Step 1 : Data pre-processing, creating the the kinship matrices, creating the fixed effects files and the target files
(relative abundance per OTU) - step 1 in 'Main_MTV_LMM_TE.R' 

#Step 2 : Applying LMM with reml using GCTA software -  using 'run_mgrm.sh' 

#Step 3 : Calculating the 'Time-explainability for each OTU - step 3 in 'Main_MTV_LMM_TE.R' 



Predicting auto-regressive OTUs
--------------------------
#Step 1 : Data pre-processing, creating the the kinship matrices, creating the fixed effects files and the target files (relative abundance per OTU)- step 1 in 'Main_MTV_LMM_Prediction.R'

#Step 2 : Iteratively applying LMM with reml using GCTA software - using 'run_predictions.sh' (GCTA)

#Step 3 : Generating predictions using BLUP - step 3 in 'Main_MTV_LMM_Prediction.R'

Contact: liashenhav@gmail.com
