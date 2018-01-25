mkdir Data_files/Time_explainability_1_time_point
mkdir Data_files/Time_explainability_1_time_point_FF
mkdir Data_files/GRM_files

cat Data_files/ind_NEW.txt | while read -r line

do
number="$line"

gzip Data_files/Kinship_mat_files/Kinship_3_bins.grm
gzip Data_files/Kinship_mat_files/Kinship_3_bins_shuff.grm

/home/cozygene/cozygene/software/GCTA_1.26/gcta64 --grm-gz Data_files/Kinship_mat_files/Kinship_3_bins --make-grm --out Data_files/GRM_files/GRM_3_bins_1
/home/cozygene/cozygene/software/GCTA_1.26/gcta64 --grm-gz Data_files/Kinship_mat_files/Kinship_3_bins_shuff --make-grm --out Data_files/GRM_files/GRM_3_bins_2

/home/cozygene/cozygene/software/GCTA_1.26/gcta64 --reml --mgrm Data_files/multi_grm.txt --pheno Data_files/Phen_files/phen_file__$number.txt  --reml-pred-rand --qcovar Data_files/Fixed_effect_per_OTU/prev_t_1_times__$number.txt --reml-est-fix --out Data_files/Time_explainability_1_time_point_FF/GRM_OTUs_reml__$number

/home/cozygene/cozygene/software/GCTA_1.26/gcta64 --reml --mgrm Data_files/multi_grm.txt --pheno Data_files/Phen_files/phen_file__$number.txt  --reml-pred-rand --out Data_files/Time_explainability_1_time_point/GRM_OTUs_reml__$number


done
exit 0
