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
  
    
  #Creating the mgrm files for each kinship matrix
  create_mgrm_file(t_min = round(sum(ind_data$t_points)*0.67), t_max = sum(ind_data$t_points), num_GRM = 2, prediction_flag = 1, saving_path = paste(dir_path, "Data_files", sep = ""))
  create_t_file(t_min = round(sum(ind_data$t_points)*0.67), t_max = sum(ind_data$t_points), saving_path = paste(dir_path, "Data_files", sep = ""))
  #extract the number of time points per individual from the meta data
    
    
  t_min = round(sum(ind_data$t_points)*0.67)
  t_max = sum(ind_data$t_points)
    
  dir.create("Kinship_mat_files")
  setwd("./Kinship_mat_files")
    
  for(t in t_min:t_max){
      
      start_t = 1
      end_t = t
    
      dat = Data_list$x_train_rel[c(start_t:end_t),ind] 
      scaled.dat <- scale(dat)
      scaled.dat[is.na(scaled.dat)] = 0
      
      #Create the kinship matrix
      Kinship_test = cosine(scaled.dat)
      
      #GRM file for the OTU matrix W
      GRM_1 = create_GRM(data_OTUs = scaled.dat, 
                         Kinship_mat = Kinship_test ,
                         start_t = 1,
                         end_t = end_t, 
                         id =  c(1:end_t),
                         save_flag = 1, 
                         save_path = paste(dir_path, "Data_files/Kinship_mat_files", sep = ""),
                         GRM_title = paste("Kinship_3_bins_", t, sep = ""))
      
      #extract the number of time points per individual from the meta data
      ind_data = ind_num_func(meta_data = meta_data, X_start = 1)
      
      #Permute the time points within each individual
      
      t_stop_ind = max(which(cumsum(ind_data$t_points) <= end_t))
      
      block_permut = Block_permutation(blocksizes = ind_data$t_points[c(1:(t_stop_ind + 1))], t_max = (t-1))
      
      dat_shuff = Data_list$x_train[block_permut,ind]
      scaled.dat_shuff <- scale(dat_shuff)
      scaled.dat_shuff[is.na(scaled.dat_shuff)] = 0
      
      #Create the shuffled kinship matrix
      Kinship_test_shuff = cosine(scaled.dat_shuff)
      
      #GRM file for the shuffled matrix H
      GRM_2 = create_GRM(data_OTUs = scaled.dat_shuff, 
                         Kinship_mat = Kinship_test_shuff ,
                         start_t = 1,
                         end_t = end_t, 
                         id =  c(1:end_t),
                         save_flag = 1, 
                         save_path = paste(dir_path, "Data_files/Kinship_mat_files", sep = ""),
                         GRM_title = paste("Kinship_3_bins_shuff_", t, sep = ""))
    }
  
  
  ###Craete Fixed effects files - can take a few minutes
    
  if(dir.exists(file.path(paste(dir_path, "Data_files/Fixed_effect_per_OTU", sep = ""))) == FALSE){
     create_Fixed_effect_files_per_OTU(Data = Data_list$x_train_rel, ind = ind, path_to_save = paste(dir_path,"Data_files/", sep = ""))}  
  
  ###Craete 'phenotype' files - relative abundance over time without the first 'X_start' time points
  
  if(dir.exists(file.path(paste(dir_path, "Data_files/Phen_files", sep = ""))) == FALSE){
     create_pheno_files_per_OTU(Data = Results_3_bins$rel_data_all, ind = ind, ids = c(1:ind_data$ind_num), 
                             X_start = 1, path_to_save = paste(dir_path,"Data_files/", sep = ""), time_series_len = sum(ind_data$t_points))}
  
###STEP 2 - Use GCTA by running the script run_predictions.sh###

  #run_predictions.sh

###STEP 3 - Prediction###

  if(dir.exists(file.path(paste(dir_path, "Data_files/Prediction/Prediction_1_time_point_*", sep = "")))){
    
    setwd(paste(dir_path, "Data_files/Prediction", sep = ""))
    save_path = paste(dir_path, "Data_files/Prediction", sep = "")
    
    Prediction_Results_1 = list()
    TH = 0.67
    
    for(k in 1:length(ind)){
      
      Prediction_Results <- c()
      
      Prediction_Results <- try(Prediction_function(X_train_data = Data_list$x_train, x_train_rel = Data_list$x_train_rel, Data = Results_3_bins$rel_data_all, 
                                                      index = ind, OTU = ind[k], X_start = 1, T_start = round(sum(ind_data$t_points)*TH), END = 5, 
                                                      All_individuals = c(1:ind_data$ind_num),
                                                      path_hsq_files = paste(dir_path, "Data_files/Prediction/Prediction_1_time_point_", sep = ""),
                                                      path_fixed_effect = paste(dir_path, "", sep = ""),
                                                      norm_flag = 0, Fixed_effect_flag = 0,
                                                      fixed_effect_file = "prev_t_1_times__", plot_flag = 0), silent=FALSE)
      
      if ('try-error' %in% class(Prediction_Results))
      {
        
        print(paste("OTU",ind[k], "was not evaluated"))
        OTU = ind[k]
        Prediction_Results_1[[k]] = NA
        setwd(save_path)
        save(Prediction_Results_1, file = "Prediction_Results_1.RData")
        next
      }
      
      
      Prediction_Results_1[[k]] = Prediction_Results
      OTU = ind[k]
      print(paste("OTU = ",ind[k] ,";","Prediction cor = ",round(Prediction_Results_1[[k]]$cor, 4)))
      setwd(save_path)
      save(Prediction_Results_1, file = "Prediction_Results_1.RData")
      
      
    }
    
    LMM_1 = Prediction_Rsq_calc(ind = ind, Prediction_results = Prediction_Results_1) 
    R_2_LMM_1 = LMM_1$cor_vec[ind]
    print("Prediction_R^2")
    print(R_2_LMM_1)
    
  }else{
    
    print("You skipped step 2. Run the script run_predictions.sh (GCTA)")
  }



