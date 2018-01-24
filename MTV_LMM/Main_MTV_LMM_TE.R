x<-c("vegan", "dplyr")
lapply(x, require, character.only = TRUE)
library("vegan")
library("dplyr")

print("Change dir_path")
dir_path = paste("./MTV_LMM/")
setwd(dir_path)
source("src.R")

setwd(paste(dir_path, "Data_files", sep = ""))

###STEP 1 
  
  OTU_table = read.csv("otu_table_example.csv", header = T, fill = T, row.names = 1)
  meta_data = read.csv("metadata_example.csv", header = TRUE)
  names(meta_data) = c("sample_id",  "ind_id" , "Id_num", "ind_time", "Sampling_day") # This is the meta-data format
  
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
  
  #Create the kinship matrix
  Kinship_test = cosine(scaled.dat)
  
  dir.create("Kinship_mat_files")
  setwd("./Kinship_mat_files")
  
  #GRM file for the OTU matrix W
  GRM_1 = create_GRM(data_OTUs = scaled.dat, 
                     Kinship_mat = Kinship_test ,
                     start_t = 1,
                     end_t = dim(Data_list$x_train)[1], 
                     id =  c(1:dim(Data_list$x_train)[1]),
                     save_flag = 1, 
                     save_path = paste(dir_path, "Data_files/Kinship_mat_files", sep = ""),
                     GRM_title = "Kinship_3_bins")
  
  #Permute the time points within each individual
  block_permut = Block_permutation(blocksizes = ind_data$t_points, t_max = NA)
  
  dat_shuff = Data_list$x_train[block_permut,ind]
  scaled.dat_shuff <- scale(dat_shuff)
  scaled.dat_shuff[is.na(scaled.dat_shuff)] = 0
  
  #Create the shuffled kinship matrix
  Kinship_test_shuff = cosine(scaled.dat_shuff)
  
  #GRM file for the shuffled matrix H
  GRM_2 = create_GRM(data_OTUs = scaled.dat_shuff, 
                     Kinship_mat = Kinship_test_shuff ,
                     start_t = 1,
                     end_t = dim(Data_list$x_train)[1], 
                     id =  c(1:dim(Data_list$x_train)[1]),
                     save_flag = 1, 
                     save_path = paste(dir_path, "Data_files/Kinship_mat_files", sep = ""),
                     GRM_title = "Kinship_3_bins_shuff")
  
  create_mgrm_file(t_min = 1, t_max = sum(ind_data$t_points), num_GRM = 2, prediction_flag = 0, saving_path = paste(dir_path, "Data_files", sep = ""))
  
  ###Craete Fixed effects files - can take a few minutes
  
  if(dir.exists(file.path(paste(dir_path, "Data_files/Fixed_effect_per_OTU", sep = ""))) == FALSE){
    create_Fixed_effect_files_per_OTU(Data = Data_list$x_train_rel, ind = ind, path_to_save = paste(dir_path,"Data_files/", sep = "")) }
  
  ###Craete 'phenotype' files - relative abundance over time without the first 'X_start' time points
  if(dir.exists(file.path(paste(dir_path, "Data_files/Phen_files", sep = ""))) == FALSE){
    create_pheno_files_per_OTU(Data = Results_3_bins$rel_data_all, ind = ind, ids = c(1:ind_data$ind_num), 
                              X_start = 1, path_to_save = paste(dir_path,"Data_files/", sep = ""), time_series_len = sum(ind_data$t_points))}
###STEP 2 - Use GCTA by running the script run_mgrm.sh###

  #run_mgrm.sh

###STEP 3 - Calculating 'Time_explainability'###

  if(dir.exists(file.path(paste(dir_path, "Data_files/Time_explainability_1_time_point", sep = "")))){
  
    TE_Results = Calculate_TE_func(path_hsq_files = paste(dir_path,"Data_files/Time_explainability_1_time_point", sep = ""), 
                              Data = Results_3_bins$binned_data_all, All_individuals = c(1:ind_data$ind_num),
                              num_time_points = sum(ind_data$t_points),fixed_effect_flag = 0,var_comp_num = 2,ind_1 = NULL)
                        
    colnames(TE_Results) = c("TE", "SE_TE", "Ind_effect","SE_ind_effect", "p_value", "LL0", "LL", "num_otus", "OTU_index")
    
    TE_Results = TE_Results[order(TE_Results$OTU_index),]
    p_adjust = p.adjust(p = TE_Results$p_value, method = "BH", n = length(TE_Results$p_value))
    
    
    TE_Results$p_value_adjusted = rep(NA, length(TE_Results$p_value))
    
    for(k in 1:length(ind)){
      TE_Results$p_value_adjusted[TE_Results$OTU_index == ind[k]] = p_adjust[k]
  
    }
    print(TE_Results)
    
  }else{
    
    print("You skipped step 2. Run the script run_mgrm.sh (GCTA)")
  }



