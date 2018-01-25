
#OTU table and meta-data preprocessing
Data_Preproc_func<- function(OTU_table, meta_data, Data_set_name, ids, path_to_save,
                             ind, individuals){
  
  
  setwd(path_to_save)
  
  #intersection between the metadata table and the OTU table
  col_names_all_data = meta_data$sample_id[ meta_data$sample_id %in% colnames(OTU_table)] 
  
  All_data_NEW = OTU_table[,  colnames(OTU_table) %in% col_names_all_data]
  New_data = matrix(NA, ncol = dim(All_data_NEW)[2], nrow  = dim(All_data_NEW)[1])
  
  #Re-ordering the OTU table by individuals according to the meta-data 
  for(i in 1:dim(All_data_NEW)[2]){
    
    New_data[,i] = All_data_NEW[,colnames(All_data_NEW) == col_names_all_data[i]]
    
  }
  
  colnames(New_data) = col_names_all_data
  rownames(New_data) = rownames(OTU_table)
  
  #New metadata
  New_meta_data = meta_data[ meta_data$sample_id %in% col_names_all_data, ]
  
  save(New_data, file = paste("New_data_",Data_set_name,".RData", sep = ""))
  save(New_meta_data, file = paste("New_meta_data_",Data_set_name,".RData", sep = ""))
  
  
  individual_len = c()
  individual_id = unique(New_meta_data$ind_id)
  
  individual_len_func <- function(c, meta_data){
    
    individual_len = length(meta_data$Sampling_day[ meta_data$Id_num == c])
    return(individual_len)
  }
  
  for(c in 1:ids){
    
    individual_len[c] = individual_len_func(c, meta_data = New_meta_data)
    
  }
  
  create_individual <- function(c, meta_data, OTU_table){
    
    individual_meta_data = meta_data[meta_data$Id_num == c,]
    
    sample_id = individual_meta_data$sample_id[individual_meta_data$sample_id%in% colnames(OTU_table)]
    sample_day = individual_meta_data$Sampling_day[ individual_meta_data$sample_id   %in% sample_id ]
    
    data_no_order = OTU_table[,colnames(OTU_table) %in% individual_meta_data$sample_id]
    
    individual_data = matrix(NA, nrow = dim(data_no_order)[1] , ncol = dim(data_no_order)[2])
    
    for(j in 1:dim(individual_data)[2]){
      
      individual_data[,j] = OTU_table[,colnames(OTU_table) == sample_id[j]]
      
    }
    
    colnames(individual_data) = sample_day
    rownames(individual_data) = rownames(OTU_table) 
    return(individual_data)
    
  }
  
  C = length(unique(meta_data$ind_id))
  Per_individual_Data = list()
  
  for(c in 1:C){
    
    Per_individual_Data[[c]] = create_individual(c, meta_data = New_meta_data, OTU_table = New_data)
    
  }
  
  save(Per_individual_Data, file = paste("Per_individual_Data_",Data_set_name,".RData", sep = ""))
  
  na_tofill = c()
  T_max = max(individual_len)
  C = length(individual_len) #number of unique individuals
  t_per_individual_Data  = array(NA, dim = c((dim(Per_individual_Data[[1]])[1]+2),T_max,C)) #creating matrix for each individual
  
  for(c in 1:C){
    
    na_tofill[c] = T_max - dim(Per_individual_Data[[c]])[2]
    
    if(na_tofill[c] > 0){
      
      MAT_fill = matrix(NA, ncol =  na_tofill[c], nrow = dim(Per_individual_Data[[c]])[1])
      
      
      t_per_individual_Data[1,,c] = rep(c, T_max) 
      t_per_individual_Data[2,,c] = c(as.numeric(meta_data$Sampling_day[meta_data$Id_num == c]),rep(NA, na_tofill[c]))
      t_per_individual_Data[3:(dim(Per_individual_Data[[c]])[1]+2),,c] = as.matrix(cbind(Per_individual_Data[[c]], MAT_fill), ncol = T_max)
    }
    
    if(na_tofill[c] == 0){
      
      t_per_individual_Data[1,,c] = rep(c, T_max) 
      t_per_individual_Data[2,,c] = c(as.numeric(meta_data$Sampling_day[meta_data$Id_num == c]),rep(NA, na_tofill[c]))
      t_per_individual_Data[3:(dim(Per_individual_Data[[c]])[1]+2),,c] = as.matrix(Per_individual_Data[[c]])
    }
    
  }
  
  max_time = c()
  for(c in 1:C){
    
    max_time[c] = max(na.omit(t_per_individual_Data[2,,c]))
    
  }
  
  Max_T_all = max(max_time)
  
  
  ##Creating a time counter - for each time stamp in each individual
  
  time_counter = array(0, dim = c(1,Max_T_all,C))
  
  for(c in 1:C){
    
    for(i in 1:dim(t_per_individual_Data[,,])[2]){
      
      if(t_per_individual_Data[2,i,c] %in% seq(1,Max_T_all,1))
        
        time_counter[,t_per_individual_Data[2,i,c],c] = time_counter[,t_per_individual_Data[2,i,c],c] + 1
    }
  }
  
  #problematic individuals with double time index
  counter = 0
  time_rep_individuals = c()
  for(c in 1:C){
    
    for(t in 1:Max_T_all){
      
      if(time_counter[,t,c] > 1){
        
        counter = counter+1
        
        time_rep_individuals[counter] = c
      }
    }
  }
  
  
  #fixing the problematic individuals by keeping only the first point 
  
  if(length(time_rep_individuals) > 0){
    
    for(i in time_rep_individuals){
      
      j = 1
      while( j < length(na.omit(t_per_individual_Data[2,,i])) ) {
        
        if(t_per_individual_Data[2,j,i] == t_per_individual_Data[2,j+1,i]) {
          print(j)
          
          t_per_individual_Data[,,i] = cbind(t_per_individual_Data[,-(j+1),i], rep(NA, dim(t_per_individual_Data[,,i])[1]))
        }
        
        j = j+1
      }
      
    }
    
  }
  
  save(t_per_individual_Data, file = paste("t_per_individual_Data_",Data_set_name,".RData", sep = ""))
  
  ####Binning starts
  
  #each individual's time length:
  individual_lenght = c()
  
  individual_len_func <- function(c, meta_data){
    
    individual_len = length(meta_data$Sampling_day[ meta_data$Id_num == c])
    return(individual_len)
  }
  
  for(c in individuals){
    
    individual_lenght[c] = individual_len_func(c, meta_data = meta_data)
    
  }
  
  
  individual_lenght = na.omit(individual_lenght)
  
  print('Creating a matrix for each individual')
  C = length(individuals)
  c = 1
  individuals_data_NEW = array(NA, dim = c(length(3:dim(t_per_individual_Data[,,c])[1]),max(individual_lenght),C)) 
  
  for(c in 1:C){
    
    individuals_data_NEW[,,c] = t_per_individual_Data[ 3:dim(t_per_individual_Data[,,individuals[c]])[1] ,,individuals[c]]   
    
  }
  
  individuals_data_bin = individuals_data_NEW
  
  temp_mat = c()
  temp_mat_2 = c()
  na_fill = c()
  
  
  active_otu<- function(n){
    
    
    if(n >= bining_TH_1 & n <=  bining_TH_2)
      return(1)
    if(n > bining_TH_2 & n <=  bining_TH_3)
      return(2)
    if(n > bining_TH_3)
      return(3)
  }
  
  Active_OTUs <- function(vec){
    
    
    bining_TH_1 = min(vec); bining_TH_2 =  quantile(vec, 0.25) ; bining_TH_3 = quantile(vec, 0.75)
    
    new_vec = c()
    for(i in 1:length(vec)){
      
      n =  vec[i]
      if(n >= bining_TH_1 & n <=  bining_TH_2)
        new_vec[i]  = 1
      if(n > bining_TH_2 & n <=  bining_TH_3)
        new_vec[i]  = 2
      if(n > bining_TH_3)
        new_vec[i]  = 3
      
    }  
    
    
    # table(new_vec)
    
    return(new_vec)
    
    
  }
  
  
  #creating the binning matrix
  
  print('Binning the OTU table')
  
  individuals_data_bin = array(NA, dim = c(dim(individuals_data_NEW)[1], dim(individuals_data_NEW)[2], dim(individuals_data_NEW)[3]))
  
  temp_mat = c()
  temp_mat_2 = c()
  na_fill = c()
  
  for(c in 1:length(individuals)){
    
    
    temp_mat = matrix(individuals_data_NEW[,1:individual_lenght[c],c], nrow = dim(individuals_data_NEW[,,c])[1], 
                      ncol = individual_lenght[c])
    
    
    rownames(temp_mat) = rownames(New_data)
    colnames(temp_mat) = c(1:individual_lenght[c])
    
    
    temp_mat_2 = t(apply(temp_mat, 1, Active_OTUs) )
    
    na_fill[c] = dim(individuals_data_NEW)[2] - dim(temp_mat_2)[2]
    
    if(na_fill[c] > 0){
      
      individuals_data_bin[,,c] = matrix(c(temp_mat_2, rep(NA, nrow(temp_mat_2)*na_fill[c])), 
                                         nrow = nrow(temp_mat_2), 
                                         ncol = dim(individuals_data_NEW)[2])
    }
    
    if(na_fill[c] == 0){
      
      individuals_data_bin[,,c] = matrix(c(temp_mat_2), 
                                         nrow = dim(individuals_data_NEW)[1], 
                                         ncol = dim(individuals_data_NEW)[2])
    }
    
  }
  
  All_OTUs = rownames(New_data)[ind]
  
  # print('Generating the OTUs indicies')
  
  OTU_ind = c()
  for(i in 1:length(All_OTUs)){
    
    OTU_ind[i] = which(rownames(New_data)%in%All_OTUs[i])
    
  }
  
  print('Creating the OTUs 3-bins matrix for each individual')
  C = length(individuals)
  individuals_data_bin_all  = array(NA, dim = c(length(OTU_ind),max(individual_lenght),C)) 
  
  for(c in 1:C){
    
    individuals_data_bin_all[,,c] = individuals_data_bin[OTU_ind,,c]
    
  }
  
  
  print('Creating the OTUs count matrix for each individual')
  C = length(individuals)
  individuals_data_count_all  = array(NA, dim = c(length(OTU_ind),max(individual_lenght),C)) #creating matrix for each individual with var OTUs
  
  for(c in 1:C){
    
    individuals_data_count_all[,,c] = individuals_data_NEW[OTU_ind,,c]
    
  }
  
  
  
  print('Creating the OTUs relative abundance matrix for each individual')
  individuals_data_rel_all = array(NA, dim = c(dim(individuals_data_NEW)[1], dim(individuals_data_NEW)[2], dim(individuals_data_NEW)[3]))
  
  
  for(c in 1:C){
    for(t in 1:length(na.omit(individuals_data_NEW[1,,c]))){
      temp_sum = sum(individuals_data_NEW[,t,c])
      individuals_data_rel_all[,t,c] = individuals_data_NEW[,t,c]/temp_sum
    }
    
  }
  
  
  ######
  
  time_data = matrix(t_per_individual_Data[2,,1:C], ncol = max(individual_lenght), byrow = T)
  time_data_new = time_data
  individual_len = c()
  for(c in 1:C){
    
    for(t in 1:max(individual_lenght)){
      
      
      if(is.na(time_data[c,t])){
        
        time_data_new[c,t] = NA
      }
      
    }
    individual_len[c] = length(na.omit(time_data_new[c,]))
  }
  
  
  time_points_num = length(t_per_individual_Data[2,,1])
  time_data_new = time_data_new[,1:time_points_num]
  
  
  #### creating the adjusted individual data:
  num_otus = dim(individuals_data_bin)[1]
  all_OTU_num = dim(individuals_data_NEW)[1]
  
  binned_observation_NEW_all =  array(NA, dim = c(all_OTU_num,time_points_num,C)) #creating matrix for each individual with var OTUs
  real_observation_NEW_all = array(NA, dim = c(all_OTU_num,time_points_num,C)) #creating matrix for each individual with all OTUs
  rel_abun_observation_NEW_all = array(NA, dim = c(all_OTU_num,time_points_num,C)) #creating matrix for each individual with all OTUs
  
  diff = c()
  
  for(c in 1:C){
    
    diff[c] = (time_points_num - individual_len[c])
    
    
    binned_observation_NEW_all[,,c] = matrix( c(individuals_data_bin_all[c(1:all_OTU_num),1:individual_len[c],c], rep(NA, all_OTU_num*diff[c])),
                                              ncol = time_points_num)
    
    real_observation_NEW_all[,,c] = matrix( c(individuals_data_count_all[c(1:all_OTU_num),1:individual_len[c],c], rep(NA, all_OTU_num*diff[c])),
                                            ncol = time_points_num)
    
    
    rel_abun_observation_NEW_all[,,c] = matrix( c(individuals_data_rel_all[c(1:all_OTU_num),1:individual_len[c],c], rep(NA, all_OTU_num*diff[c])),
                                                ncol = time_points_num)
    
  }
  
  
  Result = list(binned_data_all = binned_observation_NEW_all, 
                real_data_all = real_observation_NEW_all, 
                rel_data_all = rel_abun_observation_NEW_all,
                time = time_data_new,
                All_OTUs = All_OTUs)
  
  return(Result)
}

#normalize
normalize <- function (m, norm = c("l1", "l2", "none")) {
  norm = match.arg(norm)
  if (norm == "none") 
    return(m)
  norm_vec = switch(norm, l1 = 1/rowSums(m), l2 = 1/sqrt(rowSums(m^2)))
  norm_vec[is.infinite(norm_vec)] = 0
  if (inherits(m, "sparseMatrix")) 
    Diagonal(x = norm_vec) %*% m
  else m * norm_vec
}

#calculate the time series length
length_func <- function(data, individuals){
  i = 1
  len_cow = c()
  for(c in individuals){
    
    len_cow[i] = length(na.omit(data[1,,c]))
    i = i+1
  }
  
  return(len_cow)
}

#Create y
create_y_vec <- function(data, individuals, start, feature){
  
  max_length = max(length_func(data = data, individuals = individuals))
  y = na.omit(as.vector(data[feature,c(start:max_length), individuals]))
  return(y)
}

#Create OTU matrix
create_x_mat <- function(data, individuals,start ,num_features, feature){
  
  
  if(start == 1){
    
    obs = matrix(NA, ncol = num_features, nrow = (length(na.omit(data[feature,,individuals[1]]))-1))
    obs = t(data[,c(1:(length(na.omit(data[feature,,individuals[1]]))-1)),individuals[1]])
    
    for(c in individuals[-1]){
      
      obs_to_add = matrix(NA, ncol = num_features,  nrow= (length(na.omit(data[feature,,c]))-1))
      
      obs_to_add = t(data[,c(1:(length(na.omit(data[feature,,c]))-1)),c])
      
      obs = rbind(obs, obs_to_add)
      
    }
  }
  
  if(start == 2){
    
    # data = real_observation_NEW;
    # individuals = All_individuals[-cow_to_pred];
    # num_features = 141;
    # feature = 1;
    
    obs = matrix(NA, ncol = num_features*2, nrow = (length(na.omit(data[feature,,individuals[1]]))-2))
    
    
    obs = t(data[, c(2:(length(na.omit(data[feature,,individuals[1]]))-1))  ,individuals[1]])
    add_obs = t(data[, c(1:(length(na.omit(data[feature,,individuals[1]]))-2))  ,individuals[1]])
    obs = cbind(obs, add_obs)
    
    for(c in individuals[-1]){
      
      
      obs_to_add = matrix(NA, ncol = num_features*2, nrow = (length(na.omit(data[feature,,c]))-2))
      
      obs_to_add = t(data[, c(2:(length(na.omit(data[feature,,c]))-1))  ,c])
      add = t(data[, c(1:(length(na.omit(data[feature,,c]))-2))  ,c])
      obs_to_add = cbind(obs_to_add, add)
      
      obs = rbind(obs, obs_to_add)
      
    }
  }
  
  if(start == 3){
    
    # data = Results_3_bins_M3$rel_data_Core;
    # individuals = All_individuals; start = X_start+1; feature = OTU
    # num_features = dim(Results_3_bins_M3$real_data_all)[1]
    
    obs = matrix(NA, ncol = num_features*3, nrow = (length(na.omit(data[feature,,individuals[1]]))-3))
    
    
    obs = t(data[, c(3:(length(na.omit(data[feature,,individuals[1]]))-1))  ,individuals[1]])
    add_obs_1 = t(data[, c(2:(length(na.omit(data[feature,,individuals[1]]))-2))  ,individuals[1]])
    add_obs_2 = t(data[, c(1:(length(na.omit(data[feature,,individuals[1]]))-3))  ,individuals[1]])
    obs = cbind(obs, add_obs_1, add_obs_2)
    
    for(c in individuals[-1]){
      
      
      obs_to_add = matrix(NA, ncol = num_features*3, nrow = (length(na.omit(data[feature,,c]))-3))
      
      add_1 = t(data[, c(3:(length(na.omit(data[feature,,c]))-1))  ,c])
      add_2 = t(data[, c(2:(length(na.omit(data[feature,,c]))-2))  ,c])
      add_3 = t(data[, c(1:(length(na.omit(data[feature,,c]))-3))  ,c])
      
      obs_to_add = cbind(add_1, add_2, add_3)
      
      obs = rbind(obs, obs_to_add)
      
      
      
    }
  }
  
  return(obs)
  
}

### qc - prevalent OTUs
qc_OTUs <- function(Data, threshold, create_file_flag = 0, path_to_save = NULL){
  
  p_all = c()
  ind = c()
  for(j in 1:dim(Data)[2]){
    
    if(  (length(table(Data[,j])) >= 2 & sum(table(Data[,j])[-1]) > dim(Data)[1] * threshold) 
         ||  (sum(Data[,j]) == 2*length(Data[,j])) ||  (sum(Data[,j]) == 3*length(Data[,j])) ){
      
      if( length(table(Data[,j])) >= 2 )
        p_all[j] = sum((table(Data[,j])/dim(Data)[1])[-1])
      
      if((sum(Data[,j]) == 2*length(Data[,j])))
        p_all[j] = sum((table(Data[,j])/dim(Data)[1])[-1])
      
    }
    
    else
      p_all[j] = 0
    
  }
  ind = which(p_all > 0)
  
  if(create_file_flag == T){
    setwd(path_to_save)
    write.table(ind, file = "ind_NEW.txt", quote = F, row.names = F, col.names = F)
  }

  
  return(ind)
} 

###Creating the time lagged data sets 
create_data_matrices<-function(Data_set, X_start, ids, OTU_names){
  
  ##Binned data
  feature_OTU = sample(x = c(1:dim(Data_set$binned_data_all)[1]), size = 1)
  
  x_train = create_x_mat(data = Data_set$binned_data_all  , individuals = ids, start = X_start, 
                         num_features = dim(Data_set$binned_data_all)[1], feature = feature_OTU)
  colnames(x_train) = OTU_names
  
  ##Count data
  x_train_real = create_x_mat(data = Data_set$real_data_all, individuals = ids, start = X_start, 
                              num_features = dim(Data_set$real_data_all)[1], feature = feature_OTU)
  colnames(x_train_real) = OTU_names
  
  
  ##Relative abundance data
  x_train_rel = create_x_mat(data = Data_set$rel_data_all, individuals = ids, start = X_start,
                             num_features = dim(Data_set$real_data_all)[1], feature = feature_OTU)
  
  colnames(x_train_rel) = OTU_names
  
  
  data_list = list(x_train = x_train, x_train_real = x_train_real, x_train_rel = x_train_rel)
  return(data_list)
  
}

#extract the number of individuals from and number of time points per individual from the meta_data
ind_num_func <-function(meta_data, X_start){
  
  ind_num = length(unique(meta_data$Id_num))
  t_points = c()
  for(j in 1:ind_num){
    
    t_points[j] = length(which(meta_data$Id_num == j)) - X_start
  }
  
  result = list(ind_num = ind_num, t_points = t_points)
  return(result)
}

#Block permutation within each individual's time series
Block_permutation <- function(blocksizes, t_max){
  
  if(is.na(t_max)){
    
    index = c(0,cumsum(blocksizes))
    
    ind_permutation = list()
    for(j in 1:length(blocksizes)){
      
      ind_permutation[[j]] = sample(c((index[j]+1):(index[j+1])))
    }
    
  }else{
    index = c(0,cumsum(blocksizes))
    
    ind_permutation = list()
    for(j in 1:(length(blocksizes) - 1)){
      
      ind_permutation[[j]] = sample(c((index[j]+1):(index[j+1])))
    }
    
    tmp = unlist(ind_permutation)
    
    if(length(tmp) < (t_max+1)){
      
      if(length(c((index[length(blocksizes)]):(t_max))) == 1)
        ind_permutation[[length(blocksizes)]] = (t_max+1)
      else
      ind_permutation[[length(blocksizes)]] = sample(x = c((index[length(blocksizes)]):(t_max)))
    }
    
  }
  
  ind_permutation_new = unlist(ind_permutation)
  return(ind_permutation_new)
}

#Craete Fixed effects files
create_Fixed_effect_files_per_OTU <- function(Data, ind, path_to_save){
  
  setwd(path_to_save)
  
  dir.create("Fixed_effect_per_OTU") ## create a directory for output files
  setwd(paste(path_to_save, "Fixed_effect_per_OTU", sep = ""))
  
  for( j in 1:length(ind)){
    OTU = ind[j]
    prev_t = c()
    prev_t = data.frame(rep(0, dim(Data[,ind])[1]), c(1:dim(Data[,ind])[1]), Data[,OTU])
    write.table(prev_t, paste("prev_t_1_times__", OTU, ".txt", sep = ""), col.names = FALSE, row.names = FALSE)
  }
  
} 

#Craete 'phenotype' files (OTU relative abundance over time)
create_pheno_files_per_OTU <- function(Data, ind, ids, time_series_len, X_start, path_to_save){
  
  setwd(path_to_save)
  
  dir.create("Phen_files") ## create a directory for output files
  setwd(paste(path_to_save, "Phen_files", sep = ""))
  
  for(i in ind){
    
    OTU = i
    y_train_vec = as.numeric(create_y_vec(data = Data, individuals = ids, start = X_start+1, feature = OTU)) 
    phen_file = matrix(NA, ncol = 3, nrow = time_series_len)
    id = c(1:time_series_len)
    
    phen_file[,1] = rep(0,  time_series_len)
    phen_file[,2] = id
    phen_file[,3] = y_train_vec
    
    write.table(phen_file, paste("phen_file__",OTU,".txt", sep = ""), col.names = FALSE, row.names = FALSE)
  }
  
}

#cosine distance
cosine <- function(m) {
  m_normalized <- m / sqrt(rowSums(m ^ 2))
  tcrossprod(m_normalized)
} 

#Create GRM files in GCTA format
create_GRM <- function(data_OTUs, Kinship_mat, save_flag = 1, save_path, GRM_title, start_t = 1, 
                       end_t = dim(x_train)[1], id){
  
  
  print('Crearing the data structure - part A')
  diff = end_t - start_t + 1
  index_row = rep(c(start_t:end_t), times = c(1:(diff)) )
  index_col_list = list()
  for(j in start_t:end_t){
    
    index_col_list[[j]] = seq(from = start_t,to = j,by = 1)
    
  }
  
  index_col = unlist(index_col_list)
  
  n = dim(Kinship_mat)[1]
  Size = n*(n+1)/2
  
  GRM_data_new = matrix(NA, ncol = 4, nrow = Size)
  GRM_data_new[,1] = index_row
  GRM_data_new[,2] = index_col
  
  print('Crearing the data structure - part B')
  
  for(j in 1:Size){
    
    GRM_data_new[j,4] = Kinship_mat[(index_row[j]- start_t +1), (index_col[j]- start_t +1)]
    
  }
  
  GRM_data_new[,3] = rep(dim(data_OTUs)[2], Size)
  
  # setwd("/Users/liatshenhav/Dropbox/OTU/Data_liat/Phen_files")
  # id = read.table("phen_file1.txt")[,2]
  # id = c(1:n)
  
  GRM_OTUs.grm.id = data.frame(rep(0, length(id)), id)
  
  print('Saving the data structure')
  if(save_flag == 1){
    setwd(save_path)
    write.table(GRM_data_new, paste(GRM_title, ".grm", sep = ""), row.names = FALSE, col.names = FALSE)
    write.table(GRM_OTUs.grm.id, paste(GRM_title,".grm.id", sep = ""), row.names = FALSE, col.names = FALSE)
  }
  
  return(GRM_data_new)
}


#Create t_ind.txt file
create_t_file <- function(t_min, t_max, saving_path){
  setwd(saving_path)
  t_data = seq(from = t_min, to = t_max, by = 1)
  t_file = matrix(t_data, ncol = 1, nrow = length(t_data))
  write.table(t_file, file = "ind_t.txt", row.names = F, col.names = F)
}

#Create mgrm.txt file
create_mgrm_file <- function(t_min, t_max, num_GRM, saving_path, prediction_flag = 0){
  
  setwd(saving_path)
  if(prediction_flag == 1){
  
    dir.create("multi_grm_files")
    setwd(paste(saving_path,"/multi_grm_files", sep = ""))
    t_data = seq(from = t_min, to = t_max, by = 1)
    GRM = matrix(NA, ncol = num_GRM, nrow = length(t_data))
    
    for(k in 1:length(t_data)){
      
      tmp = c()
      for(i in 1:num_GRM){
        
        GRM[k,i] = paste("Data_files/GRM_files/GRM_3_bins__",i,"_",t_data[k], sep = "")
        tmp = rbind(tmp, GRM[k,i]) 
      }
      
      write.table(tmp, paste("multi_grm_",t_data[k], ".txt", sep = ""), row.names = F, col.names = F, quote = F)
      
    }
    
    
  }else{
    GRM = matrix(NA, ncol = num_GRM, nrow = 1)
    tmp = c()
    for(i in 1:num_GRM){
      
      GRM[1,i] = paste("GRM_3_bins_",i, sep = "")
      tmp = rbind(tmp, GRM[1,i]) 
    }
    
    write.table(tmp, "multi_grm.txt", row.names = F, col.names = F, quote = F)
    
    
  }
}

#Parse hsq files from GCTA
Create_hsq_data_Fixed_effects <- function(path_hsq_files, hsq_txt, num_otus, var_comp_num, 
                                          Refactor_flag, data=NULL, individuals = All_individuals, fixed_effect_flag = 1, ind_1){
  

  setwd(path_hsq_files)
  
  if(var_comp_num == 1 & fixed_effect_flag == 1)
    OTU_data = array(NA, dim = c(13, 3, num_otus))
  
  if(var_comp_num == 1 & fixed_effect_flag == 0)
    OTU_data = array(NA, dim = c(10, 3, num_otus))
  
  if(var_comp_num == 2)
    OTU_data = array(NA, dim = c(13, 3, num_otus))
  
  if(var_comp_num == "mixed")
    OTU_data = array(NA, dim = c(13, 3, num_otus))
  
  files <- list.files(path=path_hsq_files, pattern="*.hsq", full.names=T, recursive=FALSE)
  
  index = c()
  
  for(j in 1:length(files)){
    split = strsplit(x = as.character(files[j]), split = "__")
    split_2 = strsplit(x = as.character(split[[1]][2]), split = ".hsq")
    index = c(index, as.numeric(split_2[[1]]))
  }
  
  index = sort(index)
  
  k = 1
  for(i in na.omit(index)){
    
    temp =  as.matrix(read.table(paste(hsq_txt ,i,".hsq", sep = ""), sep = "\t", fill = T, header = T))
    
    if(dim(temp)[1] == 13)
      OTU_data[,,k] = as.matrix(read.table(paste(hsq_txt ,i,".hsq", sep = ""), sep = "\t", fill = T, header = T))
    else{
      
      temp = rbind(temp, matrix(rep(NA, 9), ncol = 3))
      colnames(temp) = c("Source", "Variance", "SE")
      OTU_data[,,k] = temp
      
    }
    k = k+1
  }
  
  index_1 = which(index %in% ind_1)
  
  hsq_time = c()
  p_val_hsq  = c()
  hsq_individual = c()
  intercept = c()
  fixed_effect = c()
  logL0 = c()
  logL = c()
  se_hsq_time = c()
  hsq_ind = c()
  se_hsq_ind = c()
  
  if(var_comp_num == 1 & fixed_effect_flag == 1){
    
    for(j in 1:num_otus){
      
      hsq_time[j] = as.numeric(OTU_data[4,2,j])
      se_hsq_time[j] = as.numeric(OTU_data[4,3,j])
      p_val_hsq[j] = as.numeric(OTU_data[9,2,j])
      intercept[j] = as.numeric(OTU_data[12,1,j])
      fixed_effect[j] = as.numeric(OTU_data[13,1,j])
      logL0[j] = as.numeric(OTU_data[6,2,j])
      logL[j] = as.numeric(OTU_data[5,2,j])
    }
  }
  
  if(var_comp_num == 1 & fixed_effect_flag == 0){
    
    for(j in 1:num_otus){
      
      hsq_time[j] = as.numeric(OTU_data[4,2,j])
      se_hsq_time[j] = as.numeric(OTU_data[4,3,j])
      p_val_hsq[j] = as.numeric(OTU_data[9,2,j])
      logL0[j] = as.numeric(OTU_data[6,2,j])
      logL[j] = as.numeric(OTU_data[5,2,j])
    }
  }
  
  
  if(var_comp_num == 2 & fixed_effect_flag == 1){
    
    for(j in 1:num_otus){
      
      hsq_time[j] = as.numeric(OTU_data[4,2,j])
      se_hsq_time[j] = as.numeric(OTU_data[4,3,j])
      p_val_hsq[j] = as.numeric(OTU_data[9,2,j])
      intercept[j] = as.numeric(OTU_data[12,1,j])
      fixed_effect[j] = as.numeric(OTU_data[13,1,j])
      logL0[j] = as.numeric(OTU_data[6,2,j])
      logL[j] = as.numeric(OTU_data[5,2,j])
    }
  }
  
  if(var_comp_num == 2 & fixed_effect_flag == 0){
    
    for(j in 1:num_otus){
      
      hsq_time[j] = as.numeric(OTU_data[5,2,j])
      se_hsq_time[j] = as.numeric(OTU_data[5,3,j])
      
      hsq_ind[j] = as.numeric(OTU_data[6,2,j])
      se_hsq_ind[j] = as.numeric(OTU_data[6,3,j])
      
      p_val_hsq[j] = as.numeric(OTU_data[12,2,j])
      logL0[j] = as.numeric(OTU_data[9,2,j])
      logL[j] = as.numeric(OTU_data[8,2,j])
    }
  }
  
  if(var_comp_num == "mixed" & fixed_effect_flag == 0){
    
    for(j in 1:num_otus){
      
      if(j %in% index_1){
        
        hsq_time[j] = as.numeric(OTU_data[4,2,j])
        se_hsq_time[j] = as.numeric(OTU_data[4,3,j])
        p_val_hsq[j] = as.numeric(OTU_data[9,2,j])
        logL0[j] = as.numeric(OTU_data[6,2,j])
        logL[j] = as.numeric(OTU_data[5,2,j])
      }
      
      else{
        
        hsq_time[j] = as.numeric(OTU_data[5,2,j])
        se_hsq_time[j] = as.numeric(OTU_data[5,3,j])
        
        hsq_ind[j] = as.numeric(OTU_data[6,2,j])
        se_hsq_ind[j] = as.numeric(OTU_data[6,3,j])
        
        p_val_hsq[j] = as.numeric(OTU_data[12,2,j])
        logL0[j] = as.numeric(OTU_data[9,2,j])
        logL[j] = as.numeric(OTU_data[8,2,j])
      }
      
    }
  }
  
  if(length(data) == 0){
    
    X_start = 1
    two_len = c()
    
    for(j in 1:num_otus){
      
      y_train = as.numeric(create_y_vec(data = data, individuals = individuals, start = X_start+1, feature = index[j]))
      
      two_len[j] = length(y_train[ y_train == 2])
      
      
    }
    
    if(var_comp_num == 1)
      hsq_data = data.frame( na.omit(two_len), hsq_time, se_hsq_time, p_val_hsq)
    if(var_comp_num == 2)
      hsq_data = data.frame( na.omit(two_len), hsq_time, se_hsq_time, hsq_individual, p_val_hsq, c(1:num_otus), index)
    
  }
  
  else{
    
    if(var_comp_num == 1 & fixed_effect_flag == 1)
      hsq_data = data.frame( hsq_time, se_hsq_time, p_val_hsq,intercept, fixed_effect, logL0, logL, c(1:num_otus), index)
    if(var_comp_num == 1 & fixed_effect_flag == 0)
      hsq_data = data.frame( hsq_time, se_hsq_time,  p_val_hsq, logL0, logL, c(1:num_otus), index)
    if(var_comp_num == 2 & fixed_effect_flag == 1)
      hsq_data = data.frame( hsq_time, se_hsq_time, p_val_hsq,intercept, fixed_effect, logL0, logL, c(1:num_otus), index)
    if(var_comp_num == 2 & fixed_effect_flag == 0)
      hsq_data = data.frame( hsq_time, se_hsq_time, hsq_ind, se_hsq_ind,  p_val_hsq, logL0, logL, c(1:num_otus), index)
    if(var_comp_num == "mixed" & fixed_effect_flag == 0)
      hsq_data = data.frame( hsq_time, se_hsq_time, hsq_ind, se_hsq_ind,  p_val_hsq, logL0, logL, c(1:num_otus), index)
    
  }
  
  # hsq_data = hsq_data[order(hsq_data$hsq_time),]
  
  
  
  return(hsq_data)
  
}

#Calculate Time explainability for each OTU
Calculate_TE_func <- function(path_hsq_files,
                             Data, All_individuals, num_time_points, 
                             fixed_effect_flag = 0, var_comp_num, ind_1){
  

  files <- list.files(path=path_hsq_files, pattern="*.hsq", full.names=T, recursive=FALSE)
  N = length(files)
  
  hsq_data_NEW = Create_hsq_data_Fixed_effects(path_hsq_files = path_hsq_files,
                                               hsq_txt = "GRM_OTUs_reml__", 
                                               num_otus = N, var_comp_num = var_comp_num, 
                                               data = Data, individuals = All_individuals, fixed_effect_flag = fixed_effect_flag,
                                               ind_1 = ind_1)
  
  
  
  return(hsq_data_NEW)
  
}

#Calculate iterative prediction
Prediction_function <- function(X_train_data, Data, 
                                index, OTU, X_start, T_start, END = 5, All_individuals,
                                path_hsq_files, path_fixed_effect = NULL, norm_flag = 0, Fixed_effect_flag = 0,
                                fixed_effect_file = "prev_t_1_times__", plot_flag = 0){
  
  
  
  path_hsq_files = paste(path_hsq_files ,OTU, sep = "")
  
  if(Fixed_effect_flag == 1){
    
    setwd(path_fixed_effect)
    
    fixed_effect_file = read.table(paste(fixed_effect_file ,OTU,".txt", sep = ""))
    fixed_effect_vec = x_train_rel[,OTU]
    
    
    setwd(path_hsq_files)
    files <- list.files(path=path_hsq_files, pattern="*.blp", full.names=T, recursive=FALSE)
    
    
    N = length(files)
    
    hsq_data_NEW_Fixed_effects = Create_hsq_data_Fixed_effects(path_hsq_files = path_hsq_files,
                                                               hsq_txt = "GRM_OTUs_reml__", 
                                                               num_otus = N, var_comp_num = 1, 
                                                               data = Data, individuals = 1)
    
    
  }
  
  y_train_vec = as.numeric(create_y_vec(data = Data, 
                                        individuals = All_individuals, start = X_start+1, feature = OTU))
  
  setwd(path_hsq_files)
  files <- list.files(path=path_hsq_files, pattern="*.blp", full.names=T, recursive=FALSE)
  
  hsq_txt = "GRM_OTUs_reml__"
  
  par_index = c()
  
  for(j in 1:length(files)){
    split = strsplit(x = as.character(files[j]), split = "__")
    split_2 = strsplit(x = as.character(split[[1]][2]), split = ".indi.blp")
    par_index = c(par_index, as.numeric(split_2[[1]]))
  }
  
  par_index = sort(par_index)
  
  par_index = par_index[c(which(par_index == T_start) : length(par_index))]
  
  g_new = c()
  g_train_vec = c()
  Pred = c()
  Inter_Pred = c()
  fixed_intercept = c()
  fixed_coeff = c()
  
  
  for(t in T_start:(dim(X_train_data)[1] - 1)){
    
    
    if(t %% 100 == 0)
      print(t)
    
    par_test_mat = as.matrix(read.table(paste(hsq_txt ,par_index[(t-T_start+1)],".indi.blp", sep = ""), 
                                        sep = "\t", fill = T, header = F))[,1:6]
    
    tail(par_test_mat)
    
    g_train = par_test_mat[,4]
    
    if(t >= T_start){
      
      if(norm_flag == 1){
        
        C_data_train = X_train_data[c(1:t),index]
        # dim(C_data_train)
        # C_data_train[1:10, 1:10]
        
        dat = C_data_train
        scaled.dat <- scale(dat)
        
        scaled.dat[is.na(scaled.dat)] = 0
        
        
        cosine <- function(m) {
          m_normalized <- m / sqrt(rowSums(m ^ 2))
          tcrossprod(m_normalized)
        }
        
        
        Kinship_test = cosine(scaled.dat[-c(1:END),])
        # Kinship_test[1:10, 1:10]
        
        Inverse_A = solve(Kinship_test)
        # Inverse_A[1:10, 1:10]
        
        
        t_W = t(scaled.dat[-c(1:END),])
        # dim(t_W)
        
        u_hat = (t_W %*% Inverse_A)  %*% (g_train[-c(1:END)])/dim(t_W)[1]
      }
      
      if(norm_flag == 0){
        
        C_data_train = X_train_data[c(1:t),index]
        # dim(C_data_train)
        # C_data_train[1:10, 1:10]
        
        cosine <- function(m) {
          m_normalized <- m / sqrt(rowSums(m ^ 2))
          tcrossprod(m_normalized)
        }
        
        
        Kinship_test = cosine(C_data_train)
        # Kinship_test[1:10, 1:10]
        
        Inverse_A = solve(Kinship_test)
        # Inverse_A[1:10, 1:10]
        
        C_data_train_scaled = normalize(C_data_train, norm = "l2")
        t_W = t(C_data_train_scaled)
        # dim(t_W)
        
        u_hat = (t_W %*% Inverse_A)  %*% (g_train)/dim(t_W)[1]
        
      }
      
    }
    
    
    if(Fixed_effect_flag == 1){
      
      
      if(norm_flag == 1){
        
        C_data_test = X_train_data[c(1:(t+1)),index]
        # dim(C_data_test)
        dat_test = C_data_test
        scaled.dat_test <- scale(dat_test)
        scaled.dat_test[is.na(scaled.dat_test)] = 0
        
        g_new[t+1] = scaled.dat_test[(t+1),] %*% u_hat
      }
      
      if(norm_flag == 0){
        
        C_data_test = X_train_data[c(1:(t+1)),ind]
        # dim(C_data_test)
        C_data_test_scaled = normalize(C_data_test, norm = "l2")
        g_new[t+1] = C_data_test_scaled[(t+1),] %*% u_hat
      }
      
      
      fixed_intercept[t+1] = hsq_data_NEW_Fixed_effects$intercept[t-T_start+1]
      fixed_coeff[t+1] = hsq_data_NEW_Fixed_effects$fixed_effect[t-T_start+1]
      
      
      
      Pred[t+1] = fixed_intercept[t+1] + fixed_coeff[t+1]*fixed_effect_vec[t+1] +  g_new[t+1]
      # Pred[t+1] = fixed_intercept[T_start+1] + fixed_coeff[T_start+1]*fixed_effect_vec[t+1] +  g_new[t+1]
      
    }
    
    
    if(Fixed_effect_flag == 0){
      
      
      if(norm_flag == 1){
        
        C_data_test = X_train_data[c(1:(t+1)),index]
        # dim(C_data_test)
        dat_test = C_data_test
        scaled.dat_test <- scale(dat_test)
        scaled.dat_test[is.na(scaled.dat_test)] = 0
        
        g_new[t+1] = scaled.dat_test[(t+1),] %*% u_hat
      }
      
      if(norm_flag == 0){
        
        C_data_test = X_train_data[c(1:(t+1)),ind]
        # dim(C_data_test)
        C_data_test_scaled = normalize(C_data_test, norm = "l2")
        g_new[t+1] = C_data_test_scaled[(t+1),] %*% u_hat
      }
      
      Pred[t+1] = g_new[t+1] #Only random effects
    }
    
    
  }
  
  y = y_train_vec[(T_start+1):length(y_train_vec)]
  y_hat = na.omit(Pred)
  
  if(plot_flag == 1){
    
    
    plot(y, type = "l", ylim = c(min(y, na.omit(y_hat)), max(y, na.omit(y_hat))))
    lines(y_hat, col = "red")
    
  }
  
  cor_test = cor.test(y, y_hat)
  
  Result = list(y = y , y_hat = y_hat,cor = cor_test$estimate^2, mse = sqrt(mean((y-y_hat)^2)) , mae = sqrt(mean(abs(y-y_hat))) )
  
  return(Result)
  
}
#Calculate Prediction R^2
Prediction_Rsq_calc <- function(ind, Prediction_results){
  
  cor_vec = rep(NA, max(ind))  
  mse_vec = rep(NA, max(ind))  
  mae_vec = rep(NA, max(ind))  
  
  for(k in 1:length(Prediction_results)){
    
    if(is.na(Prediction_results[[k]]) == FALSE){
      
      y = Prediction_results[[k]]$y
      y_hat = Prediction_results[[k]]$y_hat   
      cor_vec[ind[k]] = (cov(y, y_hat)^2)/(var(y)*var(y_hat))
      mse_vec[ind[k]] = Prediction_results[[k]]$mse
      mae_vec[ind[k]] = Prediction_results[[k]]$mae
      
      
    }
    
  }
  
  Result = list(cor_vec = cor_vec, mse_vec = mse_vec, mae_vec = mae_vec)
  return(Result)
  
}
