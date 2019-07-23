# MTV-LMM

MTV-LMM (Microbial Temporal Variability Linear Mixed Model) is a scalable and computationally efficient tool to identify autoregressive taxa (i.e., taxa whose abundance can be predicted based on microbial composition at previous time-points), quantify their temporal effects (i.e., 'Time_explainability') and predict their trajectory in future time points. For more details see [Modeling the temporal dynamics of the gut microbial community in adults and infants, Shenhav et al. 2019](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006960)


Support
-----------------------

For support using MTV-LMM, please email: liashenhav@gmail.com


Software Requirements and dependencies
-----------------------

MTV-LMM is written R. In addition to R 3.4.0 (and higher), it has the following dependencies:

- "dplyr", "vegan" (R packages)

-  [GCTA_1.26 and higher. We recommend using the latest beta version](https://cnsgenomics.com/software/gcta/#Download). 
After downloading GCTA move it into your home directory using

```
mv ~/Downloads/gcta_[version] ~/. 
```



Input format
-----------------------
The input to MTV-LMM is composed of two .csv files :

count table  - A matrix of temporal samples by taxa, across multiple hosts. The first row contains the sample ids. The first column includes taxa ids. Then every consecutive column contains read counts for each sample. Note that this order must be respected (see example below).

metadata -  The first row contains the headers ('sample_id', 'ind_id', 'Id_num', 'ind_time', 'Sampling_day'). The first column contains the sample ids. The second column contains the subject ids; the third column is the index of each subject (between 1 -  number of subjects in the data). The fourth column is the time index where the scale is the experiment's sampling rate (e.g., days, weeks, months). The fifth column is the sampling day (similar to the fourth column, but on the scale of days). Note that these names must be respected  (see examples below).


Output format
-----------------------

The output is a matrix of taxa and their temporal effects. Taxa is considered to be 'autoregressive' if the Time_explainability component is significant (p_value_adjusted <= 0.05). Only qc taxa are presented (i.e., taxa that are present in at least 10% of the temporal samples).   

| VALUE  |DESCRIPTION |
| ------------- | ------------- |
| Time_explainability    | estimate of variance explained by the microbial community composition at previous time points|
| SD_Time_explainability   | standard deviation of the time_explainability  |
| Ind_effect   | estimate of variance explained by the host at previous time points |
| SD_ind_effect   | standard deviation of the host effect |
| logL0   | log-likelihood under the null hypothesis: no termpral effect |
| logL   | log-likelihood under the alternative |
| taxa_index   | index of the focal taxa in the count table|
| p_value_adjusted   | FDR adjusted p-value of the log-ratio test (null hypothesis: no termpral effect)|





Usage instructions
---------------------------

1. Clone this repository ('MTV-LMM') using the terminal.
```
git clone https://github.com/cozygene/MTV-LMM.git
```
2. Move the cloned directory into your home directory and change permissions by executing the following commands
```
mv MTV-LMM ~/. 
cd ~/MTV-LMM/MTV_LMM
chmod +x run_TE.sh
chmod +x run_Prediction.sh
```
3. Save your input data (metadata and count table) in directory 'Data_files'.

4.  Execute

```
./run_TE.sh [your_gcta_path] [count_table] [metadata_file] 
```
where

| ARGUMENT  |DESCRIPTION |
| ------------- | ------------- |
| your_gcta_path    |The path in which the file 'gcta64' is saved  (e.g., "~/gcta_1.91.3beta_mac/bin/gcta64"). In windows and ios 'gcta64' file is bin directory (in the GCTA folder). In Linux 'gcta64' is in the main directory|
| count_table   |The full name of your taxa count table, including file type (e.g., otu_table_example.csv)  |
| metadata_file   |The full name of your metadata file, including file type (e.g., metadata_example.csv)  |




Example
---------------------------

To run MTV-LMM on example data:

1. Follow steps 1 and 2 in the usage instructions above.
2.  Execute

```
./run_TE.sh [your_gcta_path] otu_table_example.csv metadata_example.csv
```


Input - 

metadata (first few rows):

| sample_id | ind_id |Id_num | ind_time | Sampling_day|
| ------------- | ------------- |------------- |-------------|-------------|
| E000823.1.8  |  E000823 | 1 | 1.8| 54 |
| E000823.2.6  |  E000823 | 1 | 2.6| 78 |
| E000823.4.0   |  E000823 | 1| 4 | 120 |
| E000823.5.0  |  E000823 | 1 | 5 | 150 |
| ...  |  ... | ... | ... | ... |
| E001958.2.0  |  E001958 | 2	 | 2 | 60 |
| E001958.2.9  |  E001958 | 2	 | 2.9 | 87 |
| E001958.4.2  |  E001958 | 2	 | 4.2 | 126 |


count table (first 4 rows and columns):

| | E000823.1.8 |E000823.2.6 | E000823.4.0| E000823.5.0|
| ------------- | ------------- |------------- |------------- |------------- |
| taxa_1  |  1 | 1 | 2| 0 |
| taxa_2  |  0 | 0 | 0|5 |
| taxa_3  |  0 | 0 | 200|0 |
| taxa_4  |  4 | 5 | 0|0 |



Output (first 2 rows) - 


| Time_explainability | SD_Time_explainability |Ind_effect | SD_ind_effect | logL0| logL| taxa_index| p_value_adjusted| 
| ------------- | ------------- |------------- |-------------|-------------|-------------|-------------|-------------|
| 0.477828  |  0.07661 | 1.00E-06 | 0.014884| 3740.364| 3776.3| 6 | 0|
| 0.130703  |  0.066163 | 0.045816 | 0.047176| 5234.124| 5238.217| 10 | 0.002760124|



