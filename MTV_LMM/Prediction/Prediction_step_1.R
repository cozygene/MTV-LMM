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
fixed_effect_flag = 0
TH = as.numeric(as.character(init$V1[5]))


###STEP 1 

OTU_table = read.csv(count_matrix, header = T, fill = T, row.names = 1)
meta_data = read.csv(metadata_file, header = TRUE)
names(meta_data) = c("sample_id",  "ind_id" , "Id_num", "ind_time", "Sampling_day") # This is the meta-data format

# meta_data = meta_data[meta_data$Id_num < 21,]

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


#Creating the mgrm files for each kinship matrix
create_mgrm_file(t_min = round(sum(ind_data$t_points)*TH), t_max = sum(ind_data$t_points), num_GRM = 1, prediction_flag = 1, saving_path = paste(dir_path, "Data_files", sep = ""))
create_t_file(t_min = round(sum(ind_data$t_points)*TH), t_max = sum(ind_data$t_points), saving_path = paste(dir_path, "Data_files", sep = ""))
#extract the number of time points per individual from the meta data


###Craete Fixed effects files - can take a few minutes

if(dir.exists(file.path(paste(dir_path, "Data_files/Fixed_effect_per_OTU", sep = ""))) == FALSE){
  create_Fixed_effect_files_per_OTU(Data = Data_list$x_train_rel, ind = ind, path_to_save = paste(dir_path,"Data_files/", sep = ""))}  

###Craete 'phenotype' files - relative abundance over time without the first 'X_start' time points

if(dir.exists(file.path(paste(dir_path, "Data_files/Phen_files", sep = ""))) == FALSE){
  create_pheno_files_per_OTU(Data = Results_3_bins$rel_data_all, ind = ind, ids = c(1:ind_data$ind_num), 
                             X_start = 1, path_to_save = paste(dir_path,"Data_files/", sep = ""), time_series_len = sum(ind_data$t_points))}

t_min = round(sum(ind_data$t_points)*TH)
t_max = sum(ind_data$t_points)


setwd(paste(dir_path, "Data_files", sep = ""))
dir.create("Kinship_mat_files")
setwd(paste(dir_path, "Data_files/Kinship_mat_files", sep = ""))

for(t in t_min:t_max){
  
  start_t = 1
  end_t = t
  
  dat = Data_list$x_train[c(start_t:end_t),ind] 
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
  
  if(length(which(cumsum(ind_data$t_points) <= end_t)) > 0){
    
    t_stop_ind = max(which(cumsum(ind_data$t_points) <= end_t))
    block_permut = Block_permutation(blocksizes = ind_data$t_points[c(1:(t_stop_ind + 1))], t_max = (t-1))
  }
  
  else{
    index = c(1:t)
    block_permut = sample(index)
  }
  

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


