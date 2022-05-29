#0;136;0c#!/bin/bash>
 
#### Rscript
 
Rscript=/software/R-4.1.0/bin/Rscript



MASTER_ROUTE=$1

path_bash_scripts=$2

path_Rscripts=$(echo "$path_bash_scripts""Rscripts/")
path_Dependencies=$(echo "$path_bash_scripts""Dependencies/")

indir=$3

mem=$4
pc=$5
queue=$6

output="$path_bash_scripts""MPRA_pipeline_bash_print.sh"

touch $output
echo -n "" > $output

echo "#!/bin/bash"  >> $output

fdr_Threshold=$(echo '0.05')
enhancer_empirical_log_pval_Threshold=$(echo '1.3')
enhancer_pval_Threshold=$(echo '1.3')
ASE_log_pval_Threshold=$(echo '1.3')
output_dir=$MASTER_ROUTE

# # CAREFUL!!!!

#rm -rf $output_dir
#mkdir -p $output_dir

 
  
Rscript_data_wrangling_normalization_metrics=$(echo "$path_Rscripts""251_MPRA_1_DW_PLUS_NOR.R")



type=$(echo "data_wrangling_normalization_metrics")
outfile_data_wrangling_normalization_metrics=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_data_wrangling_normalization_metrics
echo -n "" > $outfile_data_wrangling_normalization_metrics
name_data_wrangling_normalization_metrics=$(echo "$type""_job")



LONG_MATRIX=$(echo "$path_Dependencies""Correspondence_barcode_variant.txt")
Initial_Selection=$(echo "$path_Dependencies""ER_Labelling_Initial_Selection.rds") 


echo "bsub -G team151 -o $outfile_data_wrangling_normalization_metrics -M $mem  -J $name_data_wrangling_normalization_metrics -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_data_wrangling_normalization_metrics \\" >> $output
echo "--LONG_MATRIX $LONG_MATRIX \\" >> $output
echo "--Initial_Selection $Initial_Selection \\" >> $output
echo "--indir $indir \\" >> $output
echo "--type $type --out $output_dir\"" >> $output



#bash $output
#exit

Rscript_QC_graphs=$(echo "$path_Rscripts""251_MPRA_2_QC_GRAPHS.R")

type=$(echo "QC_graphs")
outfile_QC_graphs=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_QC_graphs
echo -n "" > $outfile_QC_graphs
name_QC_graphs=$(echo "$type""_job")


output_dir=$MASTER_ROUTE


#echo "bsub -G team151 -o $outfile_QC_graphs -M $mem  -J $name_QC_graphs -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "bsub -G team151 -o $outfile_QC_graphs -M $mem  -w\"done($name_data_wrangling_normalization_metrics)\" -J $name_QC_graphs -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_QC_graphs \\" >> $output
echo "--type $type --out $output_dir\"" >> $output




echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "#####################################################################-----> POST QC <-----###############################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

echo "#####################################################################-----> MPRAnalyze <-----###############################################################################"  >> $output

K562_replicates_QC_PASS=$(echo 'K562_6,K562_7,K562_14,K562_16,K562_17,K562_18,K562_19')
CHRF_replicates_QC_PASS=$(echo 'CHRF_R1,CHRF_R2,CHRF_R11,CHRF_R12,CHRF_R13,CHRF_R14,CHRF_R15')
HL60_replicates_QC_PASS=$(echo 'HL60_R5,HL60_R7,HL60_R9,HL60_R10,HL60_R12,HL60_R14,HL60_R15')
THP1_replicates_QC_PASS=$(echo 'THP1_R0plus,THP1_R1,THP1_R6,THP1_R7,THP1_R8,THP1_R9,THP1_R10')
Initial_Selection=$(echo "$path_Dependencies""ER_Labelling_Initial_Selection.rds") 

output_dir=$MASTER_ROUTE


Rscript_QC_pass=$(echo "$path_Rscripts""251_MPRA_3_QC_PASS.R")

type=$(echo "QC_pass")
outfile_QC_pass=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_QC_pass
echo -n "" > $outfile_QC_pass
name_QC_pass=$(echo "$type""_job")

  

echo "bsub -G team151 -o $outfile_QC_pass -M $mem -w\"done($name_data_wrangling_normalization_metrics)\"  -J $name_QC_pass -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_QC_pass -M $mem -J $name_QC_pass -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_QC_pass \\" >> $output
echo "--Initial_Selection $Initial_Selection \\" >> $output
echo "--K562_replicates_QC_PASS $K562_replicates_QC_PASS \\" >> $output
echo "--CHRF_replicates_QC_PASS $CHRF_replicates_QC_PASS \\" >> $output
echo "--HL60_replicates_QC_PASS $HL60_replicates_QC_PASS \\" >> $output
echo "--THP1_replicates_QC_PASS $THP1_replicates_QC_PASS \\" >> $output
echo "--type $type --out $output_dir\"" >> $output

Cell_Type_string=$(echo 'K562;''CHRF;''HL60;''THP1;')
#Cell_Type_string=$(echo 'HL60;''THP1;')




echo $Cell_Type_string

a=($(echo "$Cell_Type_string" | tr ";" '\n'))

declare -a arr

for i  in "${a[@]}"
    do
        Cell_Type=${i}
        echo "$Cell_Type"

	Rscript_QC_pass_MPRAnalyze=$(echo "$path_Rscripts""251_MPRA_3_5_MPRAnalyze_partII.R")

	type=$(echo "QC_pass_MPRAnalyze""_""$Cell_Type")
	outfile_QC_pass_MPRAnalyze=$(echo "$output_dir""outfile""_""$type"".out")
	touch $outfile_QC_pass_MPRAnalyze
	echo -n "" > $outfile_QC_pass_MPRAnalyze
	name_QC_pass_MPRAnalyze=$(echo "$type""_job")



	echo "bsub -G team151 -o $outfile_QC_pass_MPRAnalyze -M $mem -w\"done($name_QC_pass)\"  -J $name_QC_pass_MPRAnalyze -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
	echo "\"$Rscript $Rscript_QC_pass_MPRAnalyze \\" >> $output
	echo "--Cell_Type $Cell_Type \\" >> $output
	echo "--type $type --out $output_dir\"" >> $output

	QC_pass_MPRAnalyze_string=$(echo "&& done($name_QC_pass_MPRAnalyze)")
        echo "->>>$QC_pass_MPRAnalyze_string"
	arr[${#arr[@]}]="$QC_pass_MPRAnalyze_string"
	
done




 echo "#####################################################################-----> parameters_recalculation.R  <-----###############################################################################"  >> $output


Rscript_Parameter_recalculation=$(echo "$path_Rscripts""251_MPRA_4_QC_parameters_recalculation.R")

type=$(echo "Parameter_recalculation")
outfile_Parameter_recalculation=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_Parameter_recalculation
echo -n "" > $outfile_Parameter_recalculation
name_Parameter_recalculation=$(echo "$type""_job")

K562_replicates_QC_PASS=$(echo 'K562_6,K562_7,K562_14,K562_16,K562_17,K562_18,K562_19')
CHRF_replicates_QC_PASS=$(echo 'CHRF_R1,CHRF_R2,CHRF_R11,CHRF_R12,CHRF_R13,CHRF_R14,CHRF_R15')
HL60_replicates_QC_PASS=$(echo 'HL60_R5,HL60_R7,HL60_R9,HL60_R10,HL60_R12,HL60_R14,HL60_R15')
THP1_replicates_QC_PASS=$(echo 'THP1_R0plus,THP1_R1,THP1_R6,THP1_R7,THP1_R8,THP1_R9,THP1_R10')

output_dir=$MASTER_ROUTE

echo "bsub -G team151 -o $outfile_Parameter_recalculation -M $mem -w\"done($name_QC_pass)\"  -J $name_Parameter_recalculation -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_Parameter_recalculation -M $mem  -J $name_Parameter_recalculation -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_Parameter_recalculation \\" >> $output
echo "--K562_replicates_QC_PASS $K562_replicates_QC_PASS \\" >> $output
echo "--CHRF_replicates_QC_PASS $CHRF_replicates_QC_PASS \\" >> $output
echo "--HL60_replicates_QC_PASS $HL60_replicates_QC_PASS \\" >> $output
echo "--THP1_replicates_QC_PASS $THP1_replicates_QC_PASS \\" >> $output
echo "--type $type --out $output_dir\"" >> $output

first_half_done_string=$(echo "\"""done($name_Parameter_recalculation) ")
# #echo "$first_half_done_string"

done_string=$(echo ""${arr[@]}"""\"")
# #echo "$done_string"

complete_done_string=$(echo "$first_half_done_string""$done_string")

echo "$complete_done_string"



# echo "#########################################################################################################################################################################"  >> $output
# echo "#####################################################################-----> Volcanos overview <-----###############################################################################"  >> $output

Rscript_volcano_overview=$(echo "$path_Rscripts""251_MPRA_5_Volcanos_and_ACTIVE_tile_definition_v2.R")

type=$(echo "volcano_overview")
outfile_volcano_overview=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_volcano_overview
echo -n "" > $outfile_volcano_overview
name_volcano_overview=$(echo "$type""_job")

enhancer_result_K562=$(echo "$output_dir""MPRAnalyze_enhancer_results_K562.txt")
ASE_result_K562=$(echo "$output_dir""MPRAnalyze_ASE_results_ASE_K562.txt")

enhancer_result_CHRF=$(echo "$output_dir""MPRAnalyze_enhancer_results_CHRF.txt")
ASE_result_CHRF=$(echo "$output_dir""MPRAnalyze_ASE_results_ASE_CHRF.txt")

enhancer_result_HL60=$(echo "$output_dir""MPRAnalyze_enhancer_results_HL60.txt")
ASE_result_HL60=$(echo "$output_dir""MPRAnalyze_ASE_results_ASE_HL60.txt")

enhancer_result_THP1=$(echo "$output_dir""MPRAnalyze_enhancer_results_THP1.txt")
ASE_result_THP1=$(echo "$output_dir""MPRAnalyze_ASE_results_ASE_THP1.txt")


enhancer_logval_Threshold=$(echo "1.3")
FC_Threshold=$(echo "0.05")
Vockley_REF_Threshold=$(echo "0.95,1.05")




echo "bsub -G team151 -o $outfile_volcano_overview -M $mem -w$complete_done_string  -J $name_volcano_overview -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_volcano_overview -M $mem -J $name_volcano_overview -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_volcano_overview \\" >> $output
echo "--fdr_Threshold $fdr_Threshold \\" >> $output
echo "--FC_Threshold $FC_Threshold \\" >> $output
echo "--Vockley_REF_Threshold $Vockley_REF_Threshold \\" >> $output
echo "--enhancer_logval_Threshold $enhancer_logval_Threshold \\" >> $output
echo "--ASE_log_pval_Threshold $ASE_log_pval_Threshold \\" >> $output
echo "--enhancer_result_K562 $enhancer_result_K562 \\" >> $output
echo "--ASE_result_K562 $ASE_result_K562 \\" >> $output
echo "--enhancer_result_CHRF $enhancer_result_CHRF \\" >> $output
echo "--ASE_result_CHRF $ASE_result_CHRF \\" >> $output
echo "--enhancer_result_HL60 $enhancer_result_HL60 \\" >> $output
echo "--ASE_result_HL60 $ASE_result_HL60 \\" >> $output
echo "--enhancer_result_THP1 $enhancer_result_THP1 \\" >> $output
echo "--ASE_result_THP1 $ASE_result_THP1 \\" >> $output
echo "--type $type --out $output_dir\"" >> $output






echo "#####################################################################-----> KTP_collapse <-----###############################################################################"  >> $output

Rscript_KTP_collapse=$(echo "$path_Rscripts""251_MPRA_6_KTP_collapse.R")

type=$(echo "KTP_collapse")
outfile_KTP_collapse=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_KTP_collapse
echo -n "" > $outfile_KTP_collapse
name_KTP_collapse=$(echo "$type""_""job")


ACTIVE_TILES=$(echo "$output_dir""ACTIVE_TILES.txt")
MPRA_Real_tile_QC2_PASS=$(echo "$output_dir""MPRA_Real_Tile_QC2_PASS.rds")
CSQ_colors=$(echo "$path_Dependencies""df_CSQ_colors.rds")


#echo "bsub -G team151 -o $outfile_KTP_collapse -M $mem -J $name_KTP_collapse -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output	
echo "bsub -G team151 -o $outfile_KTP_collapse -M $mem -w\"done($name_volcano_overview)\" -J $name_KTP_collapse -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_KTP_collapse \\" >> $output
echo "--MPRA_Real_tile_QC2_PASS $MPRA_Real_tile_QC2_PASS \\" >> $output
echo "--ACTIVE_TILES $ACTIVE_TILES \\" >> $output
echo "--CSQ_colors $CSQ_colors \\" >> $output
echo "--type $type --out $output_dir\"" >> $output



echo "#####################################################################-----> ACTIVE_TILES_BARPLOT <-----###############################################################################"  >> $output

Rscript_ACTIVE_TILES_BARPLOT=$(echo "$path_Rscripts""251_MPRA_7_ACTIVE_TILES_BARPLOTS.R")

type=$(echo "ACTIVE_TILES_BARPLOT")
outfile_ACTIVE_TILES_BARPLOT=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_ACTIVE_TILES_BARPLOT
echo -n "" > $outfile_ACTIVE_TILES_BARPLOT
name_ACTIVE_TILES_BARPLOT=$(echo "$type""_""job")


CSQ_colors=$(echo "$path_Dependencies""df_CSQ_colors.rds")

 
echo "bsub -G team151 -o $outfile_ACTIVE_TILES_BARPLOT -M $mem -w\"done($name_KTP_collapse)\" -J $name_ACTIVE_TILES_BARPLOT -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_ACTIVE_TILES_BARPLOT -M $mem -J $name_ACTIVE_TILES_BARPLOT -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_ACTIVE_TILES_BARPLOT \\" >> $output
echo "--CSQ_colors $CSQ_colors \\" >> $output
echo "--type $type --out $output_dir\"" >> $output

echo "#####################################################################-----> CUMULATIVE_CURVES <-----###############################################################################"  >> $output

Rscript_CUMULATIVE_CURVES=$(echo "$path_Rscripts""251_MPRA_7_5_CUMMULATIVE_CURVES.R")

type=$(echo "CUMULATIVE_CURVES")
outfile_CUMULATIVE_CURVES=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_CUMULATIVE_CURVES
echo -n "" > $outfile_CUMULATIVE_CURVES
name_CUMULATIVE_CURVES=$(echo "$type""_""job")

ACTIVE_TILES=$(echo "$output_dir""ACTIVE_TILES_for_cummulative.txt")
CSQ_colors=$(echo "$path_Dependencies""df_CSQ_colors.rds")

 
echo "bsub -G team151 -o $outfile_CUMULATIVE_CURVES -M $mem -w\"done($name_KTP_collapse)\" -J $name_CUMULATIVE_CURVES -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_CUMULATIVE_CURVES -M $mem -J $name_CUMULATIVE_CURVES -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_CUMULATIVE_CURVES \\" >> $output
echo "--CSQ_colors $CSQ_colors \\" >> $output
echo "--ACTIVE_TILES $ACTIVE_TILES \\" >> $output
echo "--type $type --out $output_dir\"" >> $output


echo "#####################################################################-----> UPSETR <-----###############################################################################"  >> $output
echo "#####################################################################-----> Tier Printer <-----###############################################################################"  >> $output

Rscript_UPSETR=$(echo "$path_Rscripts""251_MPRA_8_UPSETR_PLOTS_v2.R")

type=$(echo "UPSETR")
outfile_UPSETR=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_UPSETR
echo -n "" > $outfile_UPSETR
name_UPSETR=$(echo "$type""_""job")


ACTIVE_TILES=$(echo "$output_dir""ACTIVE_TILES.txt")
MPRA_Real_tile_QC2_PASS=$(echo "$output_dir""MPRA_Real_Tile_QC2_PASS.rds")
CSQ_colors=$(echo "$path_Dependencies""df_CSQ_colors.rds")


echo "bsub -G team151 -o $outfile_UPSETR -M $mem -w\"done($name_KTP_collapse)\" -J $name_UPSETR -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_UPSETR -M $mem  -J $name_UPSETR -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_UPSETR \\" >> $output
echo "--CSQ_colors $CSQ_colors \\" >> $output
echo "--type $type --out $output_dir\"" >> $output

  

 

echo "#####################################################################-----> LINEAGE_PLOTS <-----###############################################################################"  >> $output

Rscript_LINEAGE_PLOTS=$(echo "$path_Rscripts""251_MPRA_9_Lineage_plots.R")

type=$(echo "LINEAGE_PLOTS")
outfile_LINEAGE_PLOTS=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_LINEAGE_PLOTS
echo -n "" > $outfile_LINEAGE_PLOTS
name_LINEAGE_PLOTS=$(echo "$type""_""job")


CSQ_colors=$(echo "$path_Dependencies""df_CSQ_colors.rds")
dB=$(echo "$path_Dependencies""ALL_db.tsv")
dB=$(echo '/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv')

TOME_correspondence=$(echo "$path_Dependencies""Correspondence_phenotype_TOME.txt")
finemap_prob_Threshold=$(echo '0.1')



echo "bsub -G team151 -o $outfile_LINEAGE_PLOTS -M $mem -w\"done($name_KTP_collapse)\" -J $name_LINEAGE_PLOTS -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_LINEAGE_PLOTS -M $mem -J $name_LINEAGE_PLOTS -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_LINEAGE_PLOTS \\" >> $output
echo "--CSQ_colors $CSQ_colors \\" >> $output
echo "--dB $dB \\" >> $output
echo "--TOME_correspondence $TOME_correspondence \\" >> $output
echo "--finemap_prob_Threshold $finemap_prob_Threshold \\" >> $output
echo "--type $type --out $output_dir\"" >> $output


echo "#####################################################################-----> EFFECT_SIZE_STATS_and_LM <-----###############################################################################"  >> $output
 
Rscript_EFFECT_SIZE_STATS_and_LM=$(echo "$path_Rscripts""251_MPRA_10_Effect_sizes.R")

type=$(echo "EFFECT_SIZE_STATS_and_LM")
outfile_EFFECT_SIZE_STATS_and_LM=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_EFFECT_SIZE_STATS_and_LM
echo -n "" > $outfile_EFFECT_SIZE_STATS_and_LM
name_EFFECT_SIZE_STATS_and_LM=$(echo "$type""_""job")


dB=$(echo "$path_Dependencies""ALL_db.tsv")
dB=$(echo '/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv')
finemap_prob_Threshold=$(echo '0.1')
KEY_collpase_Plus_Variant_lineage_CLASSIFICATION=$(echo "$output_dir""LINEAGE_plots/KEY_collpase_Plus_Variant_lineage_CLASSIFICATION.rds")
TOME_correspondence=$(echo "$path_Dependencies""Correspondence_phenotype_TOME.txt")
MPRA_Real_tile_QC2_PASS=$(echo "$output_dir""MPRA_Real_Tile_QC2_PASS.rds")


echo "bsub -G team151 -o $outfile_EFFECT_SIZE_STATS_and_LM -M $mem -w\"done($name_LINEAGE_PLOTS)\" -J $name_EFFECT_SIZE_STATS_and_LM -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "bsub -G team151 -o $outfile_EFFECT_SIZE_STATS_and_LM -M $mem -J $name_EFFECT_SIZE_STATS_and_LM -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_EFFECT_SIZE_STATS_and_LM \\" >> $output
echo "--KEY_collpase_Plus_Variant_lineage_CLASSIFICATION $KEY_collpase_Plus_Variant_lineage_CLASSIFICATION \\" >> $output
echo "--dB $dB \\" >> $output
echo "--TOME_correspondence $TOME_correspondence \\" >> $output
echo "--MPRA_Real_tile_QC2_PASS $MPRA_Real_tile_QC2_PASS \\" >> $output
echo "--finemap_prob_Threshold $finemap_prob_Threshold \\" >> $output
echo "--type $type --out $output_dir\"" >> $output

echo "#####################################################################-----> EFFECT_SIZE_PLOTS <-----###############################################################################"  >> $output

Rscript_EFFECT_SIZE_PLOTS=$(echo "$path_Rscripts""251_MPRA_11_Effect_sizes_selected_plots_v2.R")

type=$(echo "EFFECT_SIZE_PLOTS")
outfile_EFFECT_SIZE_PLOTS=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_EFFECT_SIZE_PLOTS
echo -n "" > $outfile_EFFECT_SIZE_PLOTS
name_EFFECT_SIZE_PLOTS=$(echo "$type""_""job")



echo "bsub -G team151 -o $outfile_EFFECT_SIZE_PLOTS -M $mem -w\"done($name_EFFECT_SIZE_STATS_and_LM)\" -J $name_EFFECT_SIZE_PLOTS -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_EFFECT_SIZE_PLOTS -M $mem -J $name_EFFECT_SIZE_PLOTS -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_EFFECT_SIZE_PLOTS \\" >> $output
echo "--type $type --out $output_dir\"" >> $output


echo "#####################################################################-----> VJ_directionality  <-----###############################################################################"  >> $output
 

Rscript_VJ_directionality=$(echo "$path_Rscripts""251_MPRA_12_VJ_check_directionality.R")

type=$(echo "VJ_directionality")
outfile_VJ_directionality=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_VJ_directionality
echo -n "" > $outfile_VJ_directionality
name_VJ_directionality=$(echo "$type""_job")

Sankaran_MPRA=$(echo "$path_Dependencies""MPRA_Sankaran.csv")
KEY_collpase_Plus_Variant_lineage_CLASSIFICATION=$(echo "$output_dir""LINEAGE_plots/KEY_collpase_Plus_Variant_lineage_CLASSIFICATION.rds")
MPRA_Real_tile_QC2_PASS=$(echo "$output_dir""MPRA_Real_Tile_QC2_PASS.rds")
output_dir=$MASTER_ROUTE

echo "bsub -G team151 -o $outfile_VJ_directionality -M $mem -w\"done($name_LINEAGE_PLOTS) && done($name_Parameter_recalculation)\"  -J $name_VJ_directionality -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_VJ_directionality -M $mem  -w\"done($name_LINEAGE_PLOTS)\" -J $name_VJ_directionality -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_VJ_directionality \\" >> $output
echo "--Sankaran_MPRA $Sankaran_MPRA \\" >> $output
echo "--MPRA_Real_tile_QC2_PASS $MPRA_Real_tile_QC2_PASS \\" >> $output
echo "--KEY_collpase_Plus_Variant_lineage_CLASSIFICATION $KEY_collpase_Plus_Variant_lineage_CLASSIFICATION \\" >> $output
echo "--type $type --out $output_dir\"" >> $output


echo "#####################################################################-----> POST_QC_PLOTS  <-----###############################################################################"  >> $output
 
Rscript_POST_QC_graphs=$(echo "$path_Rscripts""251_MPRA_13_POST_QC_PLOTS.R")

type=$(echo "POST_QC_graphs")
outfile_POST_QC_graphs=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_POST_QC_graphs
echo -n "" > $outfile_POST_QC_graphs
name_POST_QC_graphs=$(echo "$type""_job")


output_dir=$MASTER_ROUTE
MPRA_Real_tile_QC2_PASS=$(echo "$output_dir""MPRA_Real_Tile_QC2_PASS.rds")


echo "bsub -G team151 -o $outfile_POST_QC_graphs -M $mem -w\"done($name_LINEAGE_PLOTS) && done($name_Parameter_recalculation)\" -R\"select[model==Intel_Platinum]\" -J $name_POST_QC_graphs  -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_POST_QC_graphs -M $mem  -w\"done($name_LINEAGE_PLOTS)\" -J $name_POST_QC_graphs  -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_POST_QC_graphs \\" >> $output
echo "--MPRA_Real_tile_QC2_PASS $MPRA_Real_tile_QC2_PASS \\" >> $output
echo "--type $type --out $output_dir\"" >> $output


echo "#########################################################################################################################################################################"  >> $output
echo "#####################################################################-----> Per variant graphs <-----###############################################################################"  >> $output

Rscript_per_variant_plots=$(echo "$path_Rscripts""251_MPRA_14_Per_variant_plots.R")

type=$(echo "per_variant_plots")
outfile_per_variant_plots=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_per_variant_plots
echo -n "" > $outfile_per_variant_plots
name_per_variant_plots=$(echo "$type""_job")

MPRA_Real_tile_QC2_PASS=$(echo "$output_dir""MPRA_Real_Tile_QC2_PASS.rds")

echo "bsub -G team151 -o $outfile_per_variant_plots -M $mem -w\"done($name_LINEAGE_PLOTS) && done($name_Parameter_recalculation)\"  -J $name_per_variant_plots -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_per_variant_plots -M $mem -w\"done($name_LINEAGE_PLOTS)\" -J $name_per_variant_plots -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_per_variant_plots -M $mem -J $name_per_variant_plots -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_per_variant_plots \\" >> $output
echo "--MPRA_Real_tile_QC2_PASS $MPRA_Real_tile_QC2_PASS \\" >> $output
echo "--type $type --out $output_dir\"" >> $output


Rscript_SUSHI_plots_Per_variant=$(echo "$path_Rscripts""251_MPRA_15_SUSHI_plots_Per_variant.R")
 
type=$(echo "SUSHI_plots_Per_variant")
outfile_SUSHI_plots_Per_variant=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_SUSHI_plots_Per_variant
echo -n "" > $outfile_SUSHI_plots_Per_variant
name_SUSHI_plots_Per_variant=$(echo "$type""_job")

MPRA_Real_tile_QC2_PASS=$(echo "$output_dir""MPRA_Real_Tile_QC2_PASS.rds")
LONG_MATRIX=$(echo "$path_Dependencies""Correspondence_barcode_variant.txt")


echo "bsub -G team151 -o $outfile_SUSHI_plots_Per_variant -M $mem -w\"done($name_LINEAGE_PLOTS) && done($name_Parameter_recalculation)\"  -J $name_SUSHI_plots_Per_variant -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_SUSHI_plots_Per_variant -M $mem  -w\"done($name_LINEAGE_PLOTS)\" -J $name_SUSHI_plots_Per_variant -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_SUSHI_plots_Per_variant \\" >> $output
echo "--LONG_MATRIX $LONG_MATRIX \\" >> $output
echo "--MPRA_Real_tile_QC2_PASS $MPRA_Real_tile_QC2_PASS \\" >> $output
echo "--type $type --out $output_dir\"" >> $output

# echo "#########################################################################################################################################################################"  >> $output
# echo "#####################################################################-----> Per variant graphs <-----###############################################################################"  >> $output

# Rscript_dot_plot_variant=$(echo "$path_Rscripts""251_MPRA_16_Per_variant_dot_plot.R")

# type=$(echo "dot_plot_variant")
# outfile_dot_plot_variant=$(echo "$output_dir""outfile""_""$type"".out")
# touch $outfile_dot_plot_variant
# echo -n "" > $outfile_dot_plot_variant
# name_dot_plot_variant=$(echo "$type""_job")

# MPRA_Real_tile_QC2_PASS=$(echo "$output_dir""MPRA_Real_Tile_QC2_PASS.rds")

# echo "bsub -G team151 -o $outfile_dot_plot_variant -M $mem -w\"done($name_LINEAGE_PLOTS) && done($name_Parameter_recalculation)\"  -J $name_dot_plot_variant -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
# #echo "bsub -G team151 -o $outfile_dot_plot_variant -M $mem -w\"done($name_LINEAGE_PLOTS)\" -J $name_dot_plot_variant -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
# #echo "bsub -G team151 -o $outfile_dot_plot_variant -M $mem -J $name_dot_plot_variant -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
# echo "\"$Rscript $Rscript_dot_plot_variant \\" >> $output
# echo "--MPRA_Real_tile_QC2_PASS $MPRA_Real_tile_QC2_PASS \\" >> $output
# echo "--type $type --out $output_dir\"" >> $output

bash $output

exit
