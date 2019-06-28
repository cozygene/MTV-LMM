rm(list = ls())
gc()


dir_path = "~/MTV-LMM/MTV_LMM/"
setwd(dir_path)
source("src.R")

init = read.table("init.txt")


setwd(paste(dir_path, "Data_files", sep = ""))


#Set the arguments of your data
count_matrix = as.character(init$V1[2])
metadata_file = as.character(init$V1[3])
fixed_effect_flag = 0
TH = as.numeric(as.character(init$V1[4]))



OTU_table = read.csv(count_matrix, header = T, fill = T, row.names = 1)
meta_data = read.csv(metadata_file, header = TRUE)
names(meta_data) = c("sample_id",  "ind_id" , "Id_num", "ind_time", "Sampling_day") # This is the meta-data format

#Order the metadata file by time (per subject) 
ids = as.character(unique(meta_data$ind_id))
meta_data_ordered = c()

for(i in 1:length(ids)){
  
  tmp = meta_data[which(meta_data$ind_id == ids[i]),]
  tmp = tmp[order(tmp$ind_time),]
  meta_data_ordered = rbind(meta_data_ordered, tmp)
  
}

meta_data = meta_data_ordered

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
Results_3_bins_tmp = Data_Preproc_func(OTU_table = OTU_table, meta_data = meta_data, Data_set_name = "Example",
                                   path_to_save = paste(dir_path, "Data_files/Processed_data_Files", sep = ""),
                                   ids = ind_data$ind_num, ind = c(1:dim(OTU_table)[1]), 
                                   individuals =  c(1:ind_data$ind_num))

#Create 3 data OTU matrices without the last 'X_start' time points - binned matrix (3 bins), count matrix and a relative abundance matrix
Data_list_tmp = create_data_matrices(Data_set = Results_3_bins_tmp, X_start = 1, ids = c(1:ind_data$ind_num), 
                                 OTU_names = Results_3_bins_tmp$All_OTUs) ##change to 3 individuals

#prevalent OTUs - according to the user's threshold
ind = qc_OTUs(Data = Data_list_tmp$x_train, threshold = 0.5, create_file_flag = T, 
              path_to_save = paste(dir_path, "Data_files", sep = "")) 

OTU_table_qc = matrix(NA, ncol = dim(OTU_table)[2], nrow = (length(ind) + 1))

for(j in 1:length(ind)){
  
  OTU_table_qc[j,] = as.numeric(OTU_table[ind[j],])
}

if((dim(OTU_table)[1] - length(ind)) > 1)
  other = as.numeric(apply(OTU_table[-ind,], 2, sum))
if((dim(OTU_table)[1] - length(ind)) == 1)
  other = sum(OTU_table[-ind,])
OTU_table_qc[(length(ind) + 1),] = other

colnames(OTU_table_qc) = colnames(OTU_table)
rownames(OTU_table_qc) = c(rownames(OTU_table)[ind], "Other")

ind = c(1:dim(OTU_table_qc)[1])


Results_3_bins = Data_Preproc_func(OTU_table = OTU_table_qc, meta_data = meta_data, Data_set_name = "Example",
                                   path_to_save = paste(dir_path, "Data_files/Processed_data_Files", sep = ""),
                                   ids = ind_data$ind_num, ind = c(1:dim(OTU_table_qc)[1]), 
                                   individuals =  c(1:ind_data$ind_num))

#Create 3 data OTU matrices without the last 'X_start' time points - binned matrix (3 bins), count matrix and a relative abundance matrix
Data_list = create_data_matrices(Data_set = Results_3_bins, X_start = 1, ids = c(1:ind_data$ind_num), 
                                 OTU_names = Results_3_bins$All_OTUs) ##change to 3 individuals


ind = qc_OTUs(Data = Data_list$x_train, threshold = 0.5, create_file_flag = T, 
              path_to_save = paste(dir_path, "Data_files", sep = "")) 

###STEP 3 - Prediction###


config_flag = "2_GRM"
pred_path = "Prediction_"

setwd(dir_path)
dir.create("Results")
save_path = paste(dir_path, "Results", sep = "")


if(dir.exists(file.path(paste(dir_path, "Data_files/Prediction/", sep = "")))){
  
  setwd(paste(dir_path, "Data_files/Prediction", sep = ""))  
  
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
                                                   config = config_flag, meta_data = meta_data), silent=FALSE)
    
    
    
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
    
    print(paste("Predicting taxa", k))
    
    
  }
  
  
  
}else{
  
  print("You skipped step 2. Run the script Prediction_step_2.sh ")
}
