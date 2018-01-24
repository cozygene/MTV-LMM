OTU_id="Data_files/ind_NEW.txt"
time_points="Data_files/ind_t.txt"
mkdir Data_files/Prediction

while IFS= read -r line

do 
	
	number="$line"
	
	mkdir Data_files/Prediction/Prediction_1_time_point_$number
	mkdir Data_files/Prediction/Prediction_1_time_point_FF_$number

  	while read line
	
  	do

  	t="$line"	

	  	gzip Data_files/Kinship_mat_files/Kinship_3_bins_$t.grm
	  	gzip Data_files/Kinship_mat_files/Kinship_3_bins_shuff_$t.grm

		/home/cozygene/cozygene/software/GCTA_1.26/gcta64 --grm-gz Data_files/Kinship_mat_files/Kinship_3_bins_$t --make-grm --out GRM_3_bins_1_$t
		/home/cozygene/cozygene/software/GCTA_1.26/gcta64 --grm-gz Data_files/Kinship_mat_files/Kinship_3_bins_shuff_$t --make-grm --out GRM_3_bins_2_$t

		#Random effects only
		/home/cozygene/cozygene/software/GCTA_1.26/gcta64 --reml --mgrm Data_files/multi_grm_files/multi_grm_$t.txt --pheno Data_files/Phen_files/phen_file__$number.txt --reml-pred-rand --qcovar Data_files/Fixed_effect_per_OTU/prev_t_1_times__$number.txt --reml-est-fix --out Data_files/Prediction/Prediction_1_time_point_FF_$number/GRM_OTUs_reml__$t

		#Fixed effects and random effects
		/home/cozygene/cozygene/software/GCTA_1.26/gcta64 --reml --mgrm Data_files/multi_grm_files/multi_grm_$t.txt --pheno Data_files/Phen_files/phen_file__$number.txt --reml-pred-rand --out Data_files/Prediction/Prediction_1_time_point_$number/GRM_OTUs_reml__$t

  	done < "$time_points"

done < "$OTU_id"
exit 0
