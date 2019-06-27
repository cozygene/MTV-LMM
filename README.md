# MTV-LMM

MTV-LMM (Microbial Temporal Variability Linear Mixed Model) is a scalable and computationally efficient tool to identify autoregressive taxa (i.e., taxa whose abundance can be predicted based on microbial composition at previous time-points), quantify their temporal effects (i.e., 'Time_explainability') and predict their trajectory in future time points. For more details see [Modeling the temporal dynamics of the gut microbial community in adults and infants, Shenhav et al. 2019](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006960)


Support
-----------------------

For support using MTV-LMM, please email: liashenhav@gmail.com


Software Requirements and dependencies
-----------------------

MTV-LMM is written R. In addition to R 3.4.0 (and higher), it has the following dependencies:

- "dplyr", "vegan" (R packages)

- [GCTA_1.26 and higher. We recommend using the latest beta version](https://cnsgenomics.com/software/gcta/#Download)


Input format
-----------------------
The input to MTV-LMM is composed of two .csv files :

count table  - A matrix of temporal samples by taxa, across multiple hosts. The first row contains the sample ids. The first column includes taxa ids. Then every consecutive column contains read counts for each sample. Note that this order must be respected (see example below).

metadata -  The first row contains the headers ('sample_id', 'ind_id', 'Id_num', 'ind_time', 'Sampling_day'). The first column contains the sample ids. The second column contains the subject ids; the third column is the index of each subject (between 1 -  number of subjects in the data). The fourth column is the time index where the scale is the experiment's sampling rate (e.g., days, weeks, months). The fifth column is the sampling day (similar to the fourth column, but on the scale of days). Note that these names must be respected  (see examples below).


Output format
-----------------------

The output is a matrix of autoregressive taxa and their temporal effects.  
Time_explainability = the estimate of variance explained by the microbial community composition at previous time points, 
SD_Time_explainability = the standard deviation of the time_explainability, 
Ind_effect = the estimate of variance explained by the individual at previous time points, 
SD_ind_effect  = the standard deviation of the Ind_effect, 
Fixed_effect = the effect of  the previous time points of the focal OTU (returns NA if fixed_effect_flag = 0), 
logL0 =  log-likelihood under the null (no termpral effect), 
logL = log-likelihood under the alternative, 
OTU_index = the index of the focal taxa in the count matrix, 
p_value_adjusted = the FDR adjusted p-value of the log-ratio test.


Usage instructions
---------------------------

1. Clone this repository ('MTV-LMM') and save it on your computer.
2. Save your input data (metadata and count table) in the directory 'Data_files'.
3. Open your terminal and navigave to the cloned repository. 
4. Run chmod +x run_TE.sh ; chmod +x Time_explainability/TE_step_2.sh (execute permissions). 
6. Run the file 'run_TE.sh' from the main directory after inserting the following arguments as input:


| ARGUMENT | LOCATION |DESCRIPTION |
| ------------- | ------------- |------------- |
| dir_path  |  TE_step_1.R, TE_step_3.R |The path in which you saved the main directory  (e.g., "~/Dropbox/MTV-LMM/") |
| gcta path  |  init.txt |The path in which you saved GCTA  (e.g., "/Users/liatshenhav/Downloads/gcta_1.91.3beta_mac/bin/gcta64") |
| otu_table file name  |  init.txt |The full name of your taxa count matrix, including file type (e.g., 'otu_table_example.csv')  |
| count_matrix   |  init.txt |The full name of your metadata file, including file type (e.g., 'metadata_example.csv')  |
| fixed_effect_flag  | init.txt  | If you wish to use the previous time points of the focal OTU as a fixed effect, re, otherwise = 0 |
| train-set proportion  |  init.txt |The proportion of the data used for model training |



Example
---------------------------

To run MTV-LMM on example data (using multiple sinks) do:


1. Clone this repository ('MTV-LMM') and save it on your computer.
2. Set 'dir_path' in files : E_step_1.R, TE_step_3.R to the path in which you saved the main directory  (e.g., "~/Dropbox/MTV-LMM/") .
3. Set 'gcta path' in the init.txt file to the path in which you saved GCTA  (e.g., "/Users/liatshenhav/Downloads/gcta_1.91.3beta_mac/bin/gcta64")
4. Open your terminal and navigave to the cloned repository. 
5. Run chmod +x run_TE_example.sh ; chmod +x Time_explainability/TE_step_2.sh (execute permissions).
6. Run the file 'run_TE_example.sh' from the main directory.


Input - 

metadata

| sample_id | ind_id |Id_num | ind_time | Sampling_day|
| ------------- | ------------- |------------- |-------------|-------------|
| E000823.1.8  |  E000823 | 1 | 1.8| 54 |
| E000823.2.6  |  E000823 | 1 | 2.6| 78 |
| E000823.4.0   |  E000823 | 1| 4 | 120 |
| E000823.5.0  |  E000823 | 1 | 5 | 150 |



count matrix (first 4 rows and columns):

| | E000823.1.8 |E000823.2.6 | E000823.4.0| E000823.5.0|
| ------------- | ------------- |------------- |------------- |------------- |
| taxa_1  |  1 | 1 | 2| 0 |
| taxa_2  |  0 | 0 | 0|5 |
| taxa_3  |  0 | 0 | 200|0 |
| taxa_4  |  4 | 5 | 0|0 |



Output - 


| Time_explainability | SD_Time_explainability |Ind_effect | SD_ind_effect | intercept | Fixed_effect| logL0| logL| OTU_index| p_value_adjusted| 
| ------------- | ------------- |------------- |-------------|-------------|-------------|-------------|-------------|-------------|-------------|
| 0.469181  |  0.046004 | 0.031580 | 0.022000|  0.000022 | 0.870526| 8136.467| 8145.550| 5890 | 3.542700e-05|
| 0.419181  |  0.033438 | 0.010985 | 0.011545| 0.003023 | 0.531058| 4590.989| 4610.062| 7005 | 5.055050e-09|

