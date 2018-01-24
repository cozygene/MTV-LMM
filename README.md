#MTV_LMM


Software Requirements:
-----------------------

- R version 3.4.0
- GCTA_1.26 and above


Content
-----------------------
The MTV-LMM package contains:
- 'otu_table_example.csv': A simulated time-series OTU table (7244 OTUs x 1101 Time points, 40 individuals).
- 'metadata_example.csv' :A simulated meta-data file. 
- 'src.R' : Source code for MTV_LMM.
- 'Main_MTV_LMM_TE.R': Calculating Time_explainability for the simulated data set 'otu_table_example.csv'.  
- 'Main_MTV_LMM_Prediction.R' - Generating the predictions for the simulated data set 'otu_table_example.csv'.
- 'run_mgrm.sh': Shell code, running GCTA software as part of the Time_explainability calculation.
- 'run_predictions.sh': Shell code, running GCTA software  as part of the prediction generation.


Running MTV_LMM - Time_explainability
--------------------------
#Step 1 : Data pre-processing, creating the the kinship matrices, creating the fixed 		  effects files and the target files (relative abundance per OTU). 
	- Run step 1 in 'Main_MTV_LMM_TE.R' 

#Step 2 : Applying LMM with reml using GCTA software. 
	- Run step 2 using 'run_mgrm.sh' 

#Step 3 : Calculating the 'Time-explainability for each OTU.
	- Run step 3 in 'Main_MTV_LMM_TE.R' 



Running MTV_LMM - Time_explainability
--------------------------
#Step 1 : Data pre-processing, creating the the kinship matrices, creating the fixed        
	  effects files and the target files (relative abundance per OTU). 
	- Run step 1 in 'Main_MTV_LMM_Prediction.R'

#Step 2 : Iteratively applying LMM with reml using GCTA software.  
	- Run step 2 using 'run_predictions.sh' (GCTA)

#Step 3 : Generating predictions using BLUP. 
	- Run step 3 in 'Main_MTV_LMM_Prediction.R'

Contact: liashenhav@gmail.com
