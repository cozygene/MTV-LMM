# rm(list = ls())
# gc()


dir_path = paste('')
if (dir_path == '') {
  message("Please tell me where to find the code by setting an environment variable 'dir_path' ", dir_path)
}

setwd(dir_path)
source("src.R")
init = read.table("init.txt")

setwd(paste(dir_path, "Data_files", sep = ""))


#Set the arguments of your data
count_matrix = as.character(init$V1[2])
metadata_file = as.character(init$V1[3])
fixed_effect_flag = 0


###STEP 1 

OTU_table = read.csv(count_matrix, header = T, fill = T, row.names = 1)
meta_data = read.csv(metadata_file, header = TRUE)
names(meta_data) = c("sample_id",  "ind_id" , "Id_num", "ind_time", "Sampling_day") # This is the meta-data format

# Extract only those samples in common between the two tables
common.sample.ids <- intersect(meta_data$sample_id, colnames(OTU_table))
OTU_table <- OTU_table[,common.sample.ids]
meta_data$sample_id <- meta_data$sample_id
# Double-check that the mapping file and otu table
# had overlapping samples
if(length(common.sample.ids) <= 1) {
   message <- paste(sprintf('Error: there are %d sample ids in common '),
                     'between the metadata file and data table')
   stop(message)
 }


#extract the number of time points per individual from the meta data
ind_data = ind_num_func(meta_data = meta_data, X_start = 1)

dir.create("Processed_data_Files")
Results_3_bins = Data_Preproc_func(OTU_table = OTU_table, meta_data = meta_data, Data_set_name = "Example",
                                   path_to_save = paste(dir_path, "Data_files/Processed_data_Files", sep = ""),
                                   ids = ind_data$ind_num, ind = c(1:dim(OTU_table)[1]), 
                                   individuals =  c(1:ind_data$ind_num))

#Create 3 data OTU matrices without the last 'X_start' time points - binned matrix (3 bins), count matrix and a relative abundance matrix
Data_list = create_data_matrices(Data_set = Results_3_bins, X_start = 1, ids = c(1:ind_data$ind_num), 
                                 OTU_names = Results_3_bins$All_OTUs) ##change to 3 individuals

#prevalent OTUs - according to the user's threshold
ind = qc_OTUs(Data = Data_list$x_train, threshold = 0.1, create_file_flag = T, 
              path_to_save = paste(dir_path, "Data_files", sep = "")) 

dat = Data_list$x_train[,ind] 
scaled.dat <- scale(dat)
scaled.dat[is.na(scaled.dat)] = 0


###STEP 3 - Calculating 'Time_explainability'###

  if(fixed_effect_flag == 1){
    config_flag = "2_GRM_FF"
    TE_path = "Time_explainability_FF"
  }else{
    config_flag = "2_GRM"
    TE_path = "Time_explainability"
  }



  if(dir.exists(file.path(paste0(dir_path, "Data_files/", TE_path)))){

    TE_Results = Calculate_TE_func(path_hsq_files = paste0(dir_path,"Data_files/", TE_path),
                                   Data = Results_3_bins$binned_data_all, All_individuals = c(1:ind_data$ind_num),
                                   num_time_points = sum(ind_data$t_points),Fixed_effect_flag = fixed_effect_flag ,ind_1 = NULL,
                                   config = config_flag)
    
    colnames(TE_Results) = c("Time_explainability", "SD_Time_explainability", "Ind_effect","SD_ind_effect", "p_value", 
                             "intercept", "Fixed_effect",
                             "logL0", "logL", "num_taxa", "taxa_index")
    

    TE_Results = TE_Results[order(TE_Results$OTU_index),]
    p_adjust = p.adjust(p = TE_Results$p_value, method = "BH", n = length(TE_Results$p_value))


    TE_Results$p_value_adjusted = p_adjust
    
    TE_Results = TE_Results[,c("Time_explainability", "SD_Time_explainability", "Ind_effect","SD_ind_effect", 
                               "logL0", "logL","taxa_index", "p_value_adjusted)]


    print(TE_Results)

  }else{

    print("You skipped step 2. Run the script run_mgrm.sh (GCTA)")
  }

setwd(dir_path)
dir.create("Results")
setwd(paste(dir_path, "Results", sep = ""))
write.csv(TE_Results, file = "TE_Results.csv")

