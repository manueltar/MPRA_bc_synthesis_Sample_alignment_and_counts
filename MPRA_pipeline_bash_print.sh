#!/bin/bash
bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/outfile_data_wrangling_normalization_metrics.out -M 16000  -J data_wrangling_normalization_metrics_job -R"select[model==Intel_Platinum]" -R"select[mem>=16000] rusage[mem=16000] span[hosts=1]" -n4 -q normal -- \
"/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Rscripts/251_MPRA_1_DW_PLUS_NOR.R \
--LONG_MATRIX /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Dependencies/Correspondence_barcode_variant.txt \
--Initial_Selection /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Dependencies/ER_Labelling_Initial_Selection.rds \
--indir /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST/ \
--type data_wrangling_normalization_metrics --out /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/"
bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/outfile_QC_graphs.out -M 16000  -w"done(data_wrangling_normalization_metrics_job)" -J QC_graphs_job -R"select[model==Intel_Platinum]" -R"select[mem>=16000] rusage[mem=16000] span[hosts=1]" -n4 -q normal -- \
"/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Rscripts/251_MPRA_2_QC_GRAPHS.R \
--type QC_graphs --out /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/"
#########################################################################################################################################################################
#########################################################################################################################################################################
#####################################################################-----> POST QC <-----###############################################################################
#########################################################################################################################################################################
#########################################################################################################################################################################
#####################################################################-----> MPRAnalyze <-----###############################################################################
bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/outfile_QC_pass.out -M 16000 -w"done(data_wrangling_normalization_metrics_job)"  -J QC_pass_job -R"select[model==Intel_Platinum]" -R"select[mem>=16000] rusage[mem=16000] span[hosts=1]" -n4 -q normal -- \
"/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Rscripts/251_MPRA_3_QC_PASS.R \
--Initial_Selection /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Dependencies/ER_Labelling_Initial_Selection.rds \
--K562_replicates_QC_PASS K562_6,K562_7,K562_14,K562_16,K562_17,K562_18,K562_19 \
--CHRF_replicates_QC_PASS CHRF_R1,CHRF_R2,CHRF_R11,CHRF_R12,CHRF_R13,CHRF_R14,CHRF_R15 \
--HL60_replicates_QC_PASS HL60_R5,HL60_R7,HL60_R9,HL60_R10,HL60_R12,HL60_R14,HL60_R15 \
--THP1_replicates_QC_PASS THP1_R0plus,THP1_R1,THP1_R6,THP1_R7,THP1_R8,THP1_R9,THP1_R10 \
--type QC_pass --out /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/"
bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/outfile_QC_pass_MPRAnalyze_K562.out -M 16000 -w"done(QC_pass_job)"  -J QC_pass_MPRAnalyze_K562_job -R"select[model==Intel_Platinum]" -R"select[mem>=16000] rusage[mem=16000] span[hosts=1]" -n4 -q normal -- \
"/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Rscripts/251_MPRA_3_5_MPRAnalyze_partII.R \
--Cell_Type K562 \
--type QC_pass_MPRAnalyze_K562 --out /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/"
bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/outfile_QC_pass_MPRAnalyze_CHRF.out -M 16000 -w"done(QC_pass_job)"  -J QC_pass_MPRAnalyze_CHRF_job -R"select[model==Intel_Platinum]" -R"select[mem>=16000] rusage[mem=16000] span[hosts=1]" -n4 -q normal -- \
"/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Rscripts/251_MPRA_3_5_MPRAnalyze_partII.R \
--Cell_Type CHRF \
--type QC_pass_MPRAnalyze_CHRF --out /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/"
bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/outfile_QC_pass_MPRAnalyze_HL60.out -M 16000 -w"done(QC_pass_job)"  -J QC_pass_MPRAnalyze_HL60_job -R"select[model==Intel_Platinum]" -R"select[mem>=16000] rusage[mem=16000] span[hosts=1]" -n4 -q normal -- \
"/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Rscripts/251_MPRA_3_5_MPRAnalyze_partII.R \
--Cell_Type HL60 \
--type QC_pass_MPRAnalyze_HL60 --out /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/"
bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/outfile_QC_pass_MPRAnalyze_THP1.out -M 16000 -w"done(QC_pass_job)"  -J QC_pass_MPRAnalyze_THP1_job -R"select[model==Intel_Platinum]" -R"select[mem>=16000] rusage[mem=16000] span[hosts=1]" -n4 -q normal -- \
"/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Rscripts/251_MPRA_3_5_MPRAnalyze_partII.R \
--Cell_Type THP1 \
--type QC_pass_MPRAnalyze_THP1 --out /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/"
#####################################################################-----> parameters_recalculation.R  <-----###############################################################################
bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/outfile_Parameter_recalculation.out -M 16000 -w"done(QC_pass_job)"  -J Parameter_recalculation_job -R"select[model==Intel_Platinum]" -R"select[mem>=16000] rusage[mem=16000] span[hosts=1]" -n4 -q normal -- \
"/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Rscripts/251_MPRA_4_QC_parameters_recalculation.R \
--K562_replicates_QC_PASS K562_6,K562_7,K562_14,K562_16,K562_17,K562_18,K562_19 \
--CHRF_replicates_QC_PASS CHRF_R1,CHRF_R2,CHRF_R11,CHRF_R12,CHRF_R13,CHRF_R14,CHRF_R15 \
--HL60_replicates_QC_PASS HL60_R5,HL60_R7,HL60_R9,HL60_R10,HL60_R12,HL60_R14,HL60_R15 \
--THP1_replicates_QC_PASS THP1_R0plus,THP1_R1,THP1_R6,THP1_R7,THP1_R8,THP1_R9,THP1_R10 \
--type Parameter_recalculation --out /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/"
bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/outfile_volcano_overview.out -M 16000 -w"done(Parameter_recalculation_job) && done(QC_pass_MPRAnalyze_K562_job) && done(QC_pass_MPRAnalyze_CHRF_job) && done(QC_pass_MPRAnalyze_HL60_job) && done(QC_pass_MPRAnalyze_THP1_job)"  -J volcano_overview_job -R"select[model==Intel_Platinum]" -R"select[mem>=16000] rusage[mem=16000] span[hosts=1]" -n4 -q normal -- \
"/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Rscripts/251_MPRA_5_Volcanos_and_ACTIVE_tile_definition_v2.R \
--fdr_Threshold 0.05 \
--FC_Threshold 0.05 \
--Vockley_REF_Threshold 0.95,1.05 \
--enhancer_logval_Threshold 1.3 \
--ASE_log_pval_Threshold 1.3 \
--enhancer_result_K562 /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/MPRAnalyze_enhancer_results_K562.txt \
--ASE_result_K562 /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/MPRAnalyze_ASE_results_ASE_K562.txt \
--enhancer_result_CHRF /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/MPRAnalyze_enhancer_results_CHRF.txt \
--ASE_result_CHRF /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/MPRAnalyze_ASE_results_ASE_CHRF.txt \
--enhancer_result_HL60 /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/MPRAnalyze_enhancer_results_HL60.txt \
--ASE_result_HL60 /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/MPRAnalyze_ASE_results_ASE_HL60.txt \
--enhancer_result_THP1 /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/MPRAnalyze_enhancer_results_THP1.txt \
--ASE_result_THP1 /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/MPRAnalyze_ASE_results_ASE_THP1.txt \
--type volcano_overview --out /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/"
#####################################################################-----> KTP_collapse <-----###############################################################################
bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/outfile_KTP_collapse.out -M 16000 -w"done(volcano_overview_job)" -J KTP_collapse_job -R"select[model==Intel_Platinum]" -R"select[mem>=16000] rusage[mem=16000] span[hosts=1]" -n4 -q normal -- \
"/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Rscripts/251_MPRA_6_KTP_collapse.R \
--MPRA_Real_tile_QC2_PASS /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/MPRA_Real_Tile_QC2_PASS.rds \
--ACTIVE_TILES /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/ACTIVE_TILES.txt \
--CSQ_colors /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Dependencies/df_CSQ_colors.rds \
--type KTP_collapse --out /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/"
#####################################################################-----> ACTIVE_TILES_BARPLOT <-----###############################################################################
bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/outfile_ACTIVE_TILES_BARPLOT.out -M 16000 -w"done(KTP_collapse_job)" -J ACTIVE_TILES_BARPLOT_job -R"select[model==Intel_Platinum]" -R"select[mem>=16000] rusage[mem=16000] span[hosts=1]" -n4 -q normal -- \
"/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Rscripts/251_MPRA_7_ACTIVE_TILES_BARPLOTS.R \
--CSQ_colors /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Dependencies/df_CSQ_colors.rds \
--type ACTIVE_TILES_BARPLOT --out /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/"
#####################################################################-----> CUMULATIVE_CURVES <-----###############################################################################
bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/outfile_CUMULATIVE_CURVES.out -M 16000 -w"done(KTP_collapse_job)" -J CUMULATIVE_CURVES_job -R"select[model==Intel_Platinum]" -R"select[mem>=16000] rusage[mem=16000] span[hosts=1]" -n4 -q normal -- \
"/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Rscripts/251_MPRA_7_5_CUMMULATIVE_CURVES.R \
--CSQ_colors /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Dependencies/df_CSQ_colors.rds \
--ACTIVE_TILES /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/ACTIVE_TILES_for_cummulative.txt \
--type CUMULATIVE_CURVES --out /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/"
#####################################################################-----> UPSETR <-----###############################################################################
#####################################################################-----> Tier Printer <-----###############################################################################
bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/outfile_UPSETR.out -M 16000 -w"done(KTP_collapse_job)" -J UPSETR_job -R"select[model==Intel_Platinum]" -R"select[mem>=16000] rusage[mem=16000] span[hosts=1]" -n4 -q normal -- \
"/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Rscripts/251_MPRA_8_UPSETR_PLOTS_v2.R \
--CSQ_colors /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Dependencies/df_CSQ_colors.rds \
--type UPSETR --out /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/"
#####################################################################-----> LINEAGE_PLOTS <-----###############################################################################
bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/outfile_LINEAGE_PLOTS.out -M 16000 -w"done(KTP_collapse_job)" -J LINEAGE_PLOTS_job -R"select[model==Intel_Platinum]" -R"select[mem>=16000] rusage[mem=16000] span[hosts=1]" -n4 -q normal -- \
"/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Rscripts/251_MPRA_9_Lineage_plots.R \
--CSQ_colors /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Dependencies/df_CSQ_colors.rds \
--dB /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv \
--TOME_correspondence /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Dependencies/Correspondence_phenotype_TOME.txt \
--finemap_prob_Threshold 0.1 \
--type LINEAGE_PLOTS --out /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/"
#####################################################################-----> EFFECT_SIZE_STATS_and_LM <-----###############################################################################
bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/outfile_EFFECT_SIZE_STATS_and_LM.out -M 16000 -w"done(LINEAGE_PLOTS_job)" -J EFFECT_SIZE_STATS_and_LM_job -R"select[model==Intel_Platinum]" -R"select[mem>=16000] rusage[mem=16000] span[hosts=1]" -n4 -q normal -- \
bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/outfile_EFFECT_SIZE_STATS_and_LM.out -M 16000 -J EFFECT_SIZE_STATS_and_LM_job -R"select[model==Intel_Platinum]"  -R"select[mem>=16000] rusage[mem=16000] span[hosts=1]" -n4 -q normal -- \
"/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Rscripts/251_MPRA_10_Effect_sizes.R \
--KEY_collpase_Plus_Variant_lineage_CLASSIFICATION /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/LINEAGE_plots/KEY_collpase_Plus_Variant_lineage_CLASSIFICATION.rds \
--dB /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv \
--TOME_correspondence /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Dependencies/Correspondence_phenotype_TOME.txt \
--MPRA_Real_tile_QC2_PASS /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/MPRA_Real_Tile_QC2_PASS.rds \
--finemap_prob_Threshold 0.1 \
--type EFFECT_SIZE_STATS_and_LM --out /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/"
#####################################################################-----> EFFECT_SIZE_PLOTS <-----###############################################################################
bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/outfile_EFFECT_SIZE_PLOTS.out -M 16000 -w"done(EFFECT_SIZE_STATS_and_LM_job)" -J EFFECT_SIZE_PLOTS_job -R"select[model==Intel_Platinum]" -R"select[mem>=16000] rusage[mem=16000] span[hosts=1]" -n4 -q normal -- \
"/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Rscripts/251_MPRA_11_Effect_sizes_selected_plots_v2.R \
--type EFFECT_SIZE_PLOTS --out /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/"
#####################################################################-----> VJ_directionality  <-----###############################################################################
bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/outfile_VJ_directionality.out -M 16000 -w"done(LINEAGE_PLOTS_job) && done(Parameter_recalculation_job)"  -J VJ_directionality_job -R"select[model==Intel_Platinum]" -R"select[mem>=16000] rusage[mem=16000] span[hosts=1]" -n4 -q normal -- \
"/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Rscripts/251_MPRA_12_VJ_check_directionality.R \
--Sankaran_MPRA /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Dependencies/MPRA_Sankaran.csv \
--MPRA_Real_tile_QC2_PASS /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/MPRA_Real_Tile_QC2_PASS.rds \
--KEY_collpase_Plus_Variant_lineage_CLASSIFICATION /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/LINEAGE_plots/KEY_collpase_Plus_Variant_lineage_CLASSIFICATION.rds \
--type VJ_directionality --out /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/"
#####################################################################-----> POST_QC_PLOTS  <-----###############################################################################
bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/outfile_POST_QC_graphs.out -M 16000 -w"done(LINEAGE_PLOTS_job) && done(Parameter_recalculation_job)" -R"select[model==Intel_Platinum]" -J POST_QC_graphs_job  -R"select[model==Intel_Platinum]" -R"select[mem>=16000] rusage[mem=16000] span[hosts=1]" -n4 -q normal -- \
"/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Rscripts/251_MPRA_13_POST_QC_PLOTS.R \
--MPRA_Real_tile_QC2_PASS /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/MPRA_Real_Tile_QC2_PASS.rds \
--type POST_QC_graphs --out /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/"
#########################################################################################################################################################################
#####################################################################-----> Per variant graphs <-----###############################################################################
bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/outfile_per_variant_plots.out -M 16000 -w"done(LINEAGE_PLOTS_job) && done(Parameter_recalculation_job)"  -J per_variant_plots_job -R"select[model==Intel_Platinum]" -R"select[mem>=16000] rusage[mem=16000] span[hosts=1]" -n4 -q normal -- \
"/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Rscripts/251_MPRA_14_Per_variant_plots.R \
--MPRA_Real_tile_QC2_PASS /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/MPRA_Real_Tile_QC2_PASS.rds \
--type per_variant_plots --out /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/"
bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/outfile_SUSHI_plots_Per_variant.out -M 16000 -w"done(LINEAGE_PLOTS_job) && done(Parameter_recalculation_job)"  -J SUSHI_plots_Per_variant_job -R"select[model==Intel_Platinum]" -R"select[mem>=16000] rusage[mem=16000] span[hosts=1]" -n4 -q normal -- \
"/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Rscripts/251_MPRA_15_SUSHI_plots_Per_variant.R \
--LONG_MATRIX /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/Dependencies/Correspondence_barcode_variant.txt \
--MPRA_Real_tile_QC2_PASS /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/MPRA_Real_Tile_QC2_PASS.rds \
--type SUSHI_plots_Per_variant --out /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/"
