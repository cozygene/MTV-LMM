rm(list = ls())
gc()


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
# fixed_effect_flag can take the values 1 or 0
# if fixed_effect_flag = 1 -> previous time points of the focal OTU will be used as fixed effects
fixed_effect_flag = as.numeric(as.character(init$V1[4]))
TH = as.numeric(as.character(init$V1[5]))


###STEP 1 

OTU_table = read.csv(count_matrix, header = T, fill = T, row.names = 1)
meta_data = read.csv(metadata_file, header = TRUE)
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


###STEP 3 - Prediction###

if(fixed_effect_flag == 1){
  config_flag = "2_GRM_FF"
  pred_path = "Prediction_FF_"
}else{
  config_flag = "2_GRM"
  pred_path = "Prediction_"
}

if(dir.exists(file.path(paste(dir_path, "Data_files/Prediction/", sep = "")))){
  
  setwd(paste(dir_path, "Data_files/Prediction", sep = ""))
  save_path = paste(dir_path, "Data_files/Prediction", sep = "")

  
  Prediction_Results_1 = list()
  
  for(k in 1:length(ind)){
    
    Prediction_Results <- c()
    
  
    Prediction_Results <-  try(Prediction_function(X_train_data = Data_list$x_train, 
                            Data = Results_3_bins$rel_data_all, 
                            x_train_rel = Data_list$x_train_rel, 
                            index = ind, OTU = ind[k], 
                            X_start = 1, T_start = round(sum(ind_data$t_points)*TH), 
                            All_individuals = c(1:ind_data$ind_num),
                            path_hsq = paste0(dir_path, "Data_files/Prediction/",pred_path),
                            path_fixed_effect =  paste0(dir_path, "Data_files/Fixed_effect_per_OTU/"), 
                            norm_flag = 1, 
                            Fixed_effect_flag = fixed_effect_flag,
                            fixed_effect_file = "prev_t_1_times__", 
                            plot_flag = 0, 
                            config = config_flag), silent=FALSE)
    
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
    setwd(save_path)
    save(Prediction_Results_1, file = "Prediction_Results_1.RData")
    
    
  }
  

  
}else{
  
  print("You skipped step 2. Run the script run_predictions.sh (GCTA)")
}
