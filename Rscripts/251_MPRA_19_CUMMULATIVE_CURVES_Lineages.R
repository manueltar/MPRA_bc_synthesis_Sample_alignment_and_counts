
suppressMessages(library("plyr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("data.table", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("crayon", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("withr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("ggplot2", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("optparse", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("dplyr", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("withr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("backports", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("broom", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("rstudioapi", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("tzdb", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("cli", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("tidyverse", lib.loc="/nfs/team151/software/manuel_R_libs_4_1//"))
library("svglite", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("cowplot",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("digest",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("farver",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("labeling",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("reshape2",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")

library("ggeasy",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")

library("viridisLite",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")

library("RColorBrewer",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")

library("viridisLite",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
suppressMessages(library("cowplot", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("ggforce", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("ggnewscale", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("scales", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))

opt = NULL

options(warn=1)

Classification_of_variants_by_lineages=function(option_list)
{
  suppressMessages(library("nph", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ and transform out2 ----
  
  out2 = opt$out2
  
  cat("out2_\n")
  cat(sprintf(as.character(out2)))
  cat("\n")
  
  ### CUMMULATIVE_CLASSES #----
  
  CUMMULATIVE_CLASSES<-readRDS(file=opt$CUMMULATIVE_CLASSES)
  
  
  cat("CUMMULATIVE_CLASSES\n")
  cat(str(CUMMULATIVE_CLASSES))
  cat("\n")
  cat(str(unique(CUMMULATIVE_CLASSES$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(CUMMULATIVE_CLASSES$Lineages)))))
  cat("\n")
  cat(sprintf(as.character(summary(CUMMULATIVE_CLASSES$Lineages))))
  cat("\n")
  
  CUMMULATIVE_CLASSES$comparison_VAR<-gsub("^chr","",CUMMULATIVE_CLASSES$VAR)
  
  Condition_DEBUG <- 1
  
  if(Condition_DEBUG == 1)
  {
    cat("CUMMULATIVE_CLASSES\n")
    cat(str(CUMMULATIVE_CLASSES))
    cat("\n")
    cat(sprintf(as.character(names(summary(CUMMULATIVE_CLASSES$Lineages)))))
    cat("\n")
    cat(sprintf(as.character(summary(CUMMULATIVE_CLASSES$Lineages))))
    cat("\n")
  }
  
  CUMMULATIVE_CLASSES_restricted<-CUMMULATIVE_CLASSES[which(CUMMULATIVE_CLASSES$carried_variants == CUMMULATIVE_CLASSES$comparison_VAR &
                                                              CUMMULATIVE_CLASSES$Label_2 == "ASSAYED_VARIANT"),]
  
  if(Condition_DEBUG == 1)
  {
    cat("CUMMULATIVE_CLASSES_restricted_0\n")
    cat(str(CUMMULATIVE_CLASSES_restricted))
    cat("\n")
  }
  
  
 
  
  CUMMULATIVE_CLASSES_restricted<-droplevels(CUMMULATIVE_CLASSES_restricted[which(CUMMULATIVE_CLASSES_restricted$Cell_Type != "ALL_CT"),])
  
  CUMMULATIVE_CLASSES_restricted$enhancer_classif<-"NA"
  CUMMULATIVE_CLASSES_restricted$E_Plus_ASE_classif<-"NA"
  
  if(Condition_DEBUG == 1)
  {
    cat("CUMMULATIVE_CLASSES_restricted_1\n")
    cat(str(CUMMULATIVE_CLASSES_restricted))
    cat("\n")
    cat(str(unique(CUMMULATIVE_CLASSES_restricted$VAR)))
    cat("\n")
    cat(str(unique(CUMMULATIVE_CLASSES_restricted$KEY_Plus_carried_variants)))
    cat("\n")
  }
  
  #### READ finemap_prob_Threshold ----
  
  finemap_prob_Threshold = opt$finemap_prob_Threshold
  
  cat("finemap_prob_Threshold\n")
  cat(sprintf(as.character(finemap_prob_Threshold)))
  cat("\n")
  
  #### Read TOME_correspondence ----
  
  TOME_correspondence = read.table(opt$TOME_correspondence, sep="\t", stringsAsFactors = F, header = T)
  
  cat("TOME_correspondence_\n")
  str(TOME_correspondence)
  cat("\n")
  
  TOME_correspondence$Lineage[TOME_correspondence$phenotype == "wbc"]<-"ALL_wbc_lineage"
  
  indx.int<-c(which(colnames(TOME_correspondence) == "phenotype"),
              which(colnames(TOME_correspondence) == "Lineage"))
  
  
  TOME_correspondence_subset<-unique(TOME_correspondence[,indx.int])
  
  cat("TOME_correspondence_subset_\n")
  str(TOME_correspondence_subset)
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(TOME_correspondence_subset$Lineage))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(TOME_correspondence_subset$Lineage)))))
  cat("\n")
  
  ### Read dB file ----
  
  dB = as.data.frame(fread(file=opt$dB, sep="\t", stringsAsFactors = F, header = T), stringsAsFactors =F)
  
  cat("dB_\n")
  cat(str(dB))
  cat("\n")
  
  
  dB_subset<-dB[which(dB$VAR%in%CUMMULATIVE_CLASSES_restricted$VAR),]
  
  cat("dB_subset_\n")
  cat(str(dB_subset))
  cat("\n")
  
  dB_subset_thresholded<-dB_subset[which(dB_subset$finemap_prob >= finemap_prob_Threshold),]
  
  cat("dB_subset_thresholded_\n")
  cat(str(dB_subset_thresholded))
  cat("\n")
  
  ##### LOOP classify variants -----
  
  
  VARS<-unique(as.character(dB_subset_thresholded$VAR))
  
  
  cat("VARS_\n")
  cat(str(VARS))
  cat("\n")
  
  list_result<-list()
  
  for(i in 1:length(VARS))
  {
    VAR_sel<-VARS[i]
    
    cat("--------------------->\t")
    cat(sprintf(as.character(VAR_sel)))
    cat("\n")
    
    dB_subset_thresholded_sel<-dB_subset_thresholded[which(dB_subset_thresholded$VAR == VAR_sel),]
    
    # cat("dB_subset_thresholded_sel_\n")
    # cat(str(dB_subset_thresholded_sel))
    # cat("\n")
    
    CLASS_erythroid<-NULL
    
    TOME_correspondence_subset_erythroid<-TOME_correspondence_subset[which(TOME_correspondence_subset$Lineage == "erythroid_lineage"),]
    
    # cat("TOME_correspondence_subset_erythroid_\n")
    # cat(str(TOME_correspondence_subset_erythroid))
    # cat("\n")
    # cat(sprintf(as.character(TOME_correspondence_subset_erythroid$phenotype)))
    # cat("\n")
    
    TOME_correspondence_subset_erythroid_VAR_sel<-TOME_correspondence_subset_erythroid[which(TOME_correspondence_subset_erythroid$phenotype%in%dB_subset_thresholded_sel$phenotype),]
    
    # cat("TOME_correspondence_subset_erythroid_VAR_sel_\n")
    # cat(str(TOME_correspondence_subset_erythroid_VAR_sel))
    # cat("\n")
    
    if(dim(TOME_correspondence_subset_erythroid_VAR_sel)[1] >0)
    {
      CLASS_erythroid<-"erythroid_lineage"
      
    }else{
      
      CLASS_erythroid<-"OTHER"
    }
    
    CLASS_gran_mono<-NULL
    
    TOME_correspondence_subset_gran_mono<-TOME_correspondence_subset[which(TOME_correspondence_subset$Lineage == "gran_mono_lineage"),]
    
    # cat("TOME_correspondence_subset_gran_mono_\n")
    # cat(str(TOME_correspondence_subset_gran_mono))
    # cat("\n")
    # cat(sprintf(as.character(TOME_correspondence_subset_gran_mono$phenotype)))
    # cat("\n")
    
    TOME_correspondence_subset_gran_mono_VAR_sel<-TOME_correspondence_subset_gran_mono[which(TOME_correspondence_subset_gran_mono$phenotype%in%dB_subset_thresholded_sel$phenotype),]
    
    # cat("TOME_correspondence_subset_gran_mono_VAR_sel_\n")
    # cat(str(TOME_correspondence_subset_gran_mono_VAR_sel))
    # cat("\n")
    
    if(dim(TOME_correspondence_subset_gran_mono_VAR_sel)[1] >0)
    {
      CLASS_gran_mono<-"gran_mono_lineage"
      
    }else{
      
      CLASS_gran_mono<-"OTHER"
    }
    
    CLASS_mega<-NULL
    
    TOME_correspondence_subset_mega_lineage<-TOME_correspondence_subset[which(TOME_correspondence_subset$Lineage == "mega_lineage"),]
    
    # cat("TOME_correspondence_subset_mega_lineage_\n")
    # cat(str(TOME_correspondence_subset_mega_lineage))
    # cat("\n")
    # cat(sprintf(as.character(TOME_correspondence_subset_mega_lineage$phenotype)))
    # cat("\n")
    
    TOME_correspondence_subset_mega_lineage_VAR_sel<-TOME_correspondence_subset_mega_lineage[which(TOME_correspondence_subset_mega_lineage$phenotype%in%dB_subset_thresholded_sel$phenotype),]
    
    # cat("TOME_correspondence_subset_mega_lineage_VAR_sel_\n")
    # cat(str(TOME_correspondence_subset_mega_lineage_VAR_sel))
    # cat("\n")
    
    
    if(dim(TOME_correspondence_subset_mega_lineage_VAR_sel)[1] >0)
    {
      CLASS_mega<-"mega_lineage"
      
    }else{
      
      CLASS_mega<-"OTHER"
    }
    
    CLASS_lymph<-NULL
    
    TOME_correspondence_subset_lymph_lineage<-TOME_correspondence_subset[which(TOME_correspondence_subset$Lineage == "lymph_lineage"),]
    
    # cat("TOME_correspondence_subset_lymph_lineage_\n")
    # cat(str(TOME_correspondence_subset_lymph_lineage))
    # cat("\n")
    # cat(sprintf(as.character(TOME_correspondence_subset_lymph_lineage$phenotype)))
    # cat("\n")
    
    
    TOME_correspondence_subset_lymph_lineage_VAR_sel<-TOME_correspondence_subset_lymph_lineage[which(TOME_correspondence_subset_lymph_lineage$phenotype%in%dB_subset_thresholded_sel$phenotype),]
    
    # cat("TOME_correspondence_subset_lymph_lineage_VAR_sel_\n")
    # cat(str(TOME_correspondence_subset_lymph_lineage_VAR_sel))
    # cat("\n")
    
    if(dim(TOME_correspondence_subset_lymph_lineage_VAR_sel)[1] >0)
    {
      CLASS_lymph<-"lymph_lineage"
      
    }else{
      
      CLASS_lymph<-"OTHER"
    }
    
    CLASS_ALL_wbc<-NULL
    
    TOME_correspondence_subset_ALL_wbc_lineage<-TOME_correspondence_subset[which(TOME_correspondence_subset$Lineage == "ALL_wbc_lineage"),]
    
    # cat("TOME_correspondence_subset_ALL_wbc_lineage_\n")
    # cat(str(TOME_correspondence_subset_ALL_wbc_lineage))
    # cat("\n")
    # cat(sprintf(as.character(TOME_correspondence_subset_ALL_wbc_lineage$phenotype)))
    # cat("\n")
    
    
    TOME_correspondence_subset_ALL_wbc_lineage_VAR_sel<-TOME_correspondence_subset_ALL_wbc_lineage[which(TOME_correspondence_subset_ALL_wbc_lineage$phenotype%in%dB_subset_thresholded_sel$phenotype),]
    
    # cat("TOME_correspondence_subset_ALL_wbc_lineage_VAR_sel_\n")
    # cat(str(TOME_correspondence_subset_ALL_wbc_lineage_VAR_sel))
    # cat("\n")
    
    if(dim(TOME_correspondence_subset_ALL_wbc_lineage_VAR_sel)[1] >0)
    {
      CLASS_ALL_wbc<-"ALL_wbc_lineage"
      
    }else{
      
      CLASS_ALL_wbc<-"OTHER"
    }
    
    A.df<-as.data.frame(cbind(VAR_sel,CLASS_erythroid,CLASS_mega,CLASS_gran_mono,CLASS_lymph,CLASS_ALL_wbc), stringsAsFactors=F)
    colnames(A.df)<-c("VAR","CLASS_erythroid","CLASS_mega","CLASS_gran_mono","CLASS_lymph","CLASS_ALL_wbc")
    
    # cat("A.df_\n")
    # cat(str(A.df))
    # cat("\n")
    
    # quit(status = 1)
    
    list_result[[i]]<-A.df
    
    
    
    
  }#i VARS
  
  
  VAR_lineage_CLASSIF = unique(as.data.frame(data.table::rbindlist(list_result, fill = T)))
  
  cat("VAR_lineage_CLASSIF_0\n")
  cat(str(VAR_lineage_CLASSIF))
  cat("\n")
  
  VAR_lineage_CLASSIF$CLASS_erythroid<-factor(VAR_lineage_CLASSIF$CLASS_erythroid,
                                              levels=c("OTHER","erythroid_lineage"),ordered=T)
  
  VAR_lineage_CLASSIF$CLASS_mega<-factor(VAR_lineage_CLASSIF$CLASS_mega,
                                         levels=c("OTHER","mega_lineage"),ordered=T)
  
  VAR_lineage_CLASSIF$CLASS_gran_mono<-factor(VAR_lineage_CLASSIF$CLASS_gran_mono,
                                              levels=c("OTHER","gran_mono_lineage"),ordered=T)
  
  VAR_lineage_CLASSIF$CLASS_lymph<-factor(VAR_lineage_CLASSIF$CLASS_lymph,
                                          levels=c("OTHER","lymph_lineage"),ordered=T)
  
  VAR_lineage_CLASSIF$CLASS_ALL_wbc<-factor(VAR_lineage_CLASSIF$CLASS_ALL_wbc,
                                            levels=c("OTHER","ALL_wbc_lineage"),ordered=T)
  
  cat("VAR_lineage_CLASSIF_1\n")
  cat(str(VAR_lineage_CLASSIF))
  cat("\n")
  cat(sprintf(as.character(names(summary(VAR_lineage_CLASSIF$CLASS_erythroid)))))
  cat("\n")
  cat(sprintf(as.character(summary(VAR_lineage_CLASSIF$CLASS_erythroid))))
  cat("\n")
  cat(sprintf(as.character(names(summary(VAR_lineage_CLASSIF$CLASS_mega)))))
  cat("\n")
  cat(sprintf(as.character(summary(VAR_lineage_CLASSIF$CLASS_mega))))
  cat("\n")
  cat(sprintf(as.character(names(summary(VAR_lineage_CLASSIF$CLASS_gran_mono)))))
  cat("\n")
  cat(sprintf(as.character(summary(VAR_lineage_CLASSIF$CLASS_gran_mono))))
  cat("\n")
  cat(sprintf(as.character(names(summary(VAR_lineage_CLASSIF$CLASS_lymph)))))
  cat("\n")
  cat(sprintf(as.character(summary(VAR_lineage_CLASSIF$CLASS_lymph))))
  cat("\n")
  cat(sprintf(as.character(names(summary(VAR_lineage_CLASSIF$CLASS_ALL_wbc)))))
  cat("\n")
  cat(sprintf(as.character(summary(VAR_lineage_CLASSIF$CLASS_ALL_wbc))))
  cat("\n")
  
  
  
  
  ####  MERGES  ####
  
  CUMMULATIVE_CLASSES_restricted<-merge(CUMMULATIVE_CLASSES_restricted,
                      VAR_lineage_CLASSIF,
                      by="VAR",
                      all.x=T)
  
  cat("CUMMULATIVE_CLASSES_restricted_\n")
  cat(str(CUMMULATIVE_CLASSES_restricted))
  cat("\n")
  cat(sprintf(as.character(names(summary(CUMMULATIVE_CLASSES_restricted$CLASS_erythroid)))))
  cat("\n")
  cat(sprintf(as.character(summary(CUMMULATIVE_CLASSES_restricted$CLASS_erythroid))))
  cat("\n")
  cat(sprintf(as.character(names(summary(CUMMULATIVE_CLASSES_restricted$CLASS_mega)))))
  cat("\n")
  cat(sprintf(as.character(summary(CUMMULATIVE_CLASSES_restricted$CLASS_mega))))
  cat("\n")
  cat(sprintf(as.character(names(summary(CUMMULATIVE_CLASSES_restricted$CLASS_gran_mono)))))
  cat("\n")
  cat(sprintf(as.character(summary(CUMMULATIVE_CLASSES_restricted$CLASS_gran_mono))))
  cat("\n")
  cat(sprintf(as.character(names(summary(CUMMULATIVE_CLASSES_restricted$CLASS_lymph)))))
  cat("\n")
  cat(sprintf(as.character(summary(CUMMULATIVE_CLASSES_restricted$CLASS_lymph))))
  cat("\n")
  cat(sprintf(as.character(names(summary(CUMMULATIVE_CLASSES_restricted$CLASS_ALL_wbc)))))
  cat("\n")
  cat(sprintf(as.character(summary(CUMMULATIVE_CLASSES_restricted$CLASS_ALL_wbc))))
  cat("\n")
  
  
  ##### LINEAGE_plots  ----
  
  
  path5<-paste(out2,'LINEAGE_CUMULATIVE_CURVES','/', sep='')
  
  cat("path5\n")
  cat(sprintf(as.character(path5)))
  cat("\n")
  
  
  if (file.exists(path5)){
    
    
    
    
  } else {
    dir.create(file.path(path5))
    
  }
  
  
  ####    SAVE ----
  
  setwd(path5)
  
  filename_1<-paste("cummulative_classification_ONLY_ASSAYED_with_Lineage_classification",".rds", sep='')
  saveRDS(file=filename_1,CUMMULATIVE_CLASSES_restricted)
}

cummulative_CT_CLASS_mega_and_Log_Rank_test_E_Plus_ASE = function(option_list)
{
  suppressMessages(library("nph", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ and transform out2 ----
  
  out2 = opt$out2
  
  cat("out2_\n")
  cat(sprintf(as.character(out2)))
  cat("\n")
  
  setwd(out)
  
  
  filename<-paste('df_Cell_colors','.rds',sep='')
  
  df_Cell_colors<-readRDS(file=filename)
  
  cat("df_Cell_colors_0\n")
  cat(str(df_Cell_colors))
  cat("\n")
  
 
  ### CUMMULATIVE_CLASSES #----
  
  path5<-paste(out2,'LINEAGE_CUMULATIVE_CURVES','/', sep='')
  
  cat("path5\n")
  cat(sprintf(as.character(path5)))
  cat("\n")
  
  setwd(path5)
  
  filename_1<-paste("cummulative_classification_ONLY_ASSAYED_with_Lineage_classification",".rds", sep='')
  
  
  CUMMULATIVE_CLASSES_restricted<-readRDS(file=filename_1)
  
  
  cat("CUMMULATIVE_CLASSES_restricted\n")
  cat(str(CUMMULATIVE_CLASSES_restricted))
  cat("\n")
  cat(str(unique(CUMMULATIVE_CLASSES_restricted$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(CUMMULATIVE_CLASSES_restricted$CLASS_mega)))))
  cat("\n")
  cat(sprintf(as.character(summary(CUMMULATIVE_CLASSES_restricted$CLASS_mega))))
  cat("\n")
  
  Condition_DEBUG <- 1
  
  ########### Build df ---------------
  
  Cell_Type_levels<-levels(CUMMULATIVE_CLASSES_restricted$Cell_Type)
  n_Cell_Type_levels<-length(levels(CUMMULATIVE_CLASSES_restricted$Cell_Type))
  n_MPRA_TILES<-5
  VAR_vector<-as.character(unique(CUMMULATIVE_CLASSES_restricted$VAR))
  KEY_Plus_carried_variants_vector<-as.character(unique(CUMMULATIVE_CLASSES_restricted$KEY_Plus_carried_variants))
  
 
  
  n_VAR<-length(VAR_vector)
  
  if(Condition_DEBUG == 1)
  {
    cat("Parameters\n")
    cat(str(Cell_Type_levels))
    cat("\n")
    cat(str(n_Cell_Type_levels))
    cat("\n")
    cat(str(n_MPRA_TILES))
    cat("\n")
    cat(str(VAR_vector))
    cat("\n")
    cat(str(KEY_Plus_carried_variants_vector))
    cat("\n")
   
   
  }
  
 
  #### Build LogRank matrix ----
  
  
  
  Log_rank_df<-as.data.frame(cbind(rep(VAR_vector,n_MPRA_TILES),
                                   rep(KEY_Plus_carried_variants_vector,n_MPRA_TILES),
                                   c(rep(1,n_VAR),rep(2,n_VAR),rep(3,n_VAR),rep(4,n_VAR),rep(5,n_VAR)),
                                   rep("NA",n_MPRA_TILES*n_VAR)),
                             stringsAsFactors=F)
  
  colnames(Log_rank_df)<-c("VAR","KEY_Plus_carried_variants","TILE","ACTIVE")
  
  
  if(Condition_DEBUG == 1)
  {
    cat("Log_rank_df_0\n")
    cat(str(Log_rank_df))
    cat("\n")
    cat(str(unique(Log_rank_df$VAR)))
    cat("\n")
    cat(str(unique(Log_rank_df$KEY_Plus_carried_variants)))
    cat("\n")
  }
  
  
  
  CUMMULATIVE_CLASSES_restricted_subset<-unique(CUMMULATIVE_CLASSES_restricted[,c(which(colnames(CUMMULATIVE_CLASSES_restricted) == "VAR"),
                                                                                which(colnames(CUMMULATIVE_CLASSES_restricted) == "KEY_Plus_carried_variants"),
                                                                                which(colnames(CUMMULATIVE_CLASSES_restricted) == "CLASS_mega"),
                                                                                which(colnames(CUMMULATIVE_CLASSES_restricted) == "Cell_Type"),
                                                                                which(colnames(CUMMULATIVE_CLASSES_restricted) == "E_Plus_ASE_CLASS_TILES"))])
  levels_CLASS_mega<-levels(CUMMULATIVE_CLASSES_restricted$CLASS_mega)
  
  if(Condition_DEBUG == 1)
  {
    cat("CUMMULATIVE_CLASSES_restricted_subset_0\n")
    cat(str(CUMMULATIVE_CLASSES_restricted_subset))
    cat("\n")
    cat(str(unique(CUMMULATIVE_CLASSES_restricted_subset$VAR)))
    cat("\n")
    cat(str(unique(CUMMULATIVE_CLASSES_restricted_subset$KEY_Plus_carried_variants)))
    cat("\n")
  }
  
  Log_rank_df<-merge(Log_rank_df,
                     CUMMULATIVE_CLASSES_restricted_subset,
                     by=c("VAR","KEY_Plus_carried_variants"),
                     all=T)
  
  Log_rank_df$CLASS_mega<-Log_rank_df$CLASS_mega
  
  Log_rank_df$CLASS_mega<-factor(Log_rank_df$CLASS_mega,
                            levels=levels_CLASS_mega,
                            ordered=T)
  
  Log_rank_df$Cell_Type<-factor(Log_rank_df$Cell_Type,
                              levels=Cell_Type_levels,
                              ordered=T)
  
  Log_rank_df$TILE<-as.integer(Log_rank_df$TILE)
  
  if(Condition_DEBUG == 1)
  {
    cat("Log_rank_df_1\n")
    cat(str(Log_rank_df))
    cat("\n")
    cat(str(unique(Log_rank_df$VAR)))
    cat("\n")
    cat(str(unique(Log_rank_df$KEY_Plus_carried_variants)))
    cat("\n")
  }
  
  
  
  
  
  ################# Classification of TILES ----------------
  
  Log_rank_df$ACTIVE[which(Log_rank_df$TILE <= Log_rank_df$E_Plus_ASE_CLASS_TILES)]<-"1"
  Log_rank_df$ACTIVE[which(Log_rank_df$TILE > Log_rank_df$E_Plus_ASE_CLASS_TILES)]<-"0"
  
  Log_rank_df$ACTIVE<-as.numeric(Log_rank_df$ACTIVE)
  
  
  
  if(Condition_DEBUG == 1)
  {
    cat("Log_rank_df_2\n")
    cat(str(Log_rank_df))
    cat("\n")
    cat(str(unique(Log_rank_df$VAR)))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(Log_rank_df$ACTIVE))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(Log_rank_df$ACTIVE)))))
    cat("\n")
  }
  
  Log_rank_df<-Log_rank_df[order(Log_rank_df$KEY_Plus_carried_variants,Log_rank_df$Cell_Type, Log_rank_df$TILE),]
  
  
  setwd(out2)
  
  write.table(Log_rank_df, file="check2.tsv",sep="\t",quote=F,row.names = F)
  
  # ################################################################
  # quit(status = 1)
  
  
  ############## Freq Table ----------------------------
  
  ### TOTAL ---
  
  CUMMULATIVE_CLASSES_restricted.dt<-data.table(CUMMULATIVE_CLASSES_restricted, key=c("Cell_Type","CLASS_mega"))
  
  Freq_TOTAL<-as.data.frame(CUMMULATIVE_CLASSES_restricted.dt[,.(TOTAL=.N),by=key(CUMMULATIVE_CLASSES_restricted.dt)], stringsAsFactors=F)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_TOTAL_0\n")
    cat(str(Freq_TOTAL))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_TOTAL.dt<-data.table(Freq_TOTAL, key=c("CLASS_mega"))
  
  
  MAXED_df<-as.data.frame(Freq_TOTAL.dt[,.(MAX=max(TOTAL)),by=key(Freq_TOTAL.dt)], stringsAsFactors=F)
  
  
  
  if(Condition_DEBUG == 1)
  {
    cat("MAXED_df_0\n")
    cat(str(MAXED_df))
    cat("\n")
    #quit(status = 1)
  }
  
  ########## Active _tiles --------------------
  
  Log_rank_df.dt<-data.table(Log_rank_df, key=c("Cell_Type","CLASS_mega","TILE"))
  
  Freq_table<-as.data.frame(Log_rank_df.dt[,.(Freq=sum(ACTIVE)),by=key(Log_rank_df.dt)], stringsAsFactors=F)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_table_0\n")
    cat(str(Freq_table))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_table<-merge(MAXED_df,
                    Freq_table,
                    by=c("CLASS_mega"),
                    all=T)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_table_1\n")
    cat(str(Freq_table))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_table$Perc<-round(100*(Freq_table$Freq/Freq_table$MAX),2)
  
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_table_2\n")
    cat(str(Freq_table))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_table$CLASS_mega<-factor(Freq_table$CLASS_mega,
                              levels=levels_CLASS_mega,
                              ordered=T)
  
  Freq_table$Cell_Type<-factor(Freq_table$Cell_Type,
                                levels=Cell_Type_levels,
                                ordered=T)
  
  Freq_table$TILE<-as.integer(Freq_table$TILE)
  
  Freq_table<-Freq_table[order(Freq_table$CLASS_mega,Freq_table$Cell_Type, Freq_table$TILE),]
  
  
  setwd(out2)
  
  write.table(Freq_table, file="check3.tsv", row.names = F,quote = F,sep="\t")
 
  # ####################################################################
  # quit(status = 1)
  
  ############### Per CT LOOP ----
  
  
  CT_levels<-Cell_Type_levels
  
  list_ABSOLUTE_DEF<-list()
  list_ABSOLUTE_DEF_Freq<-list()
  
  Condition_DEBUG <- 1
  
  for(i in 1:length(CT_levels))
  {
    CT_sel<-CT_levels[i]
    
    cat("CT:---------------------------------------->\t")
    cat(sprintf(as.character(CT_sel)))
    cat("\n")
    
    Log_rank_df_sel<-Log_rank_df[which(Log_rank_df$Cell_Type == CT_sel),]
    
    if(Condition_DEBUG == 1)
    {
      cat("Log_rank_df_sel_0\n")
      cat(str(Log_rank_df_sel))
      cat("\n")
      cat(str(unique(Log_rank_df_sel$VAR)))
      cat("\n")
    }
    
  
    #### Log Rank test ----
    
    Condition_DEBUG <- 0
    
    list_DEF<-list()
    
    CLASS_mega_labels<-levels(Log_rank_df_sel$CLASS_mega)
    
    for(k in 1 :length(CLASS_mega_labels))
    {
      BAIT<-CLASS_mega_labels[k]
      
      cat("-----BAIT------>\t")
      cat(sprintf(as.character(BAIT)))
      cat("\t")
      
      list_PREY<-list()
      
      for(h in 1 :length(CLASS_mega_labels))
      {
        
        PREY<-CLASS_mega_labels[h]
        
        cat("------PREY----->\t")
        cat(sprintf(as.character(PREY)))
        cat("\n")
        
        Mini_comparison<-unique(c(BAIT,PREY))
        
        if(Condition_DEBUG == 1)
        {
          cat("Mini_comparison_0\n")
          cat(str(Mini_comparison))
          cat("\n")
        }
        
        Log_rank_df_sel_subset<-Log_rank_df_sel[which(Log_rank_df_sel$CLASS_mega%in%Mini_comparison),]
        
        if(Condition_DEBUG == 1)
        {
          cat("Log_rank_df_sel_subset_0\n")
          cat(str(Log_rank_df_sel_subset))
          cat("\n")
        }
        
        FLAG_diversity<-as.integer(length(unique(Log_rank_df_sel_subset$CLASS_mega)))
        
        if(Condition_DEBUG == 1)
        {
          cat("FLAG_diversity_0\n")
          cat(sprintf(as.character(FLAG_diversity)))
          cat("\n")
          
          cat(str(FLAG_diversity))
          cat("\n")
        }
        
        if(FLAG_diversity >1)
        {
          cat("Hello_world\n")
          
          Mini_comparison_c1<-as.character(Mini_comparison[1])
          Mini_comparison_c2<-as.character(Mini_comparison[2])
          
          FLAG_diversity_2<-length(unique(Log_rank_df_sel_subset$ACTIVE))
          
          if(Condition_DEBUG == 1)
          {
            cat("FLAG_diversity_2_0\n")
            cat(str(FLAG_diversity_2))
            cat("\n")
          }
          
          if(FLAG_diversity_2 >1)
          {
            if(Condition_DEBUG == 1)
            {
              cat("----------->\t")
              cat(sprintf(as.character(Mini_comparison_c1)))
              cat("\t")
              cat(sprintf(as.character(Mini_comparison_c2)))
              cat("\n")
              
              cat("Log_rank_df_sel_subset_0\n")
              cat(str(Log_rank_df_sel_subset))
              cat("\n")
              cat(str(unique(Log_rank_df_sel_subset$VAR)))
              cat("\n")
              cat(sprintf(as.character(names(summary(as.factor(Log_rank_df_sel_subset$ACTIVE))))))
              cat("\n")
              cat(sprintf(as.character(summary(as.factor(Log_rank_df_sel_subset$ACTIVE)))))
              cat("\n")
            }
            
            #, "less", "greater"
            
            STATS_test<-logrank.test(
              Log_rank_df_sel_subset$TILE,
              Log_rank_df_sel_subset$ACTIVE,
              Log_rank_df_sel_subset$CLASS_mega,
              alternative = c("two.sided"),
              rho = 0,
              gamma = 0,
              event_time_weights = NULL
            )
            
            # if(Condition_DEBUG == 1)
            # {
            #   cat("STATS_test_0\n")
            #   cat(str(STATS_test))
            #   cat("\n")
            # }
            
            
            p_value<-STATS_test$test$p
            minus_log_val<-round(-1*log10(p_value),2)
            
            if(Condition_DEBUG == 1)
            {
              cat("minus_log_val_0\n")
              cat(str(minus_log_val))
              cat("\n")
            }
            
            
            a.df<-as.data.frame(cbind(Mini_comparison_c1,Mini_comparison_c2,minus_log_val), stringsAsFactors=F)
            
            
            colnames(a.df)<-c("CLASS_mega_c1","CLASS_mega_c2","minus_logpval")
            a.df$minus_logpval<-as.numeric(a.df$minus_logpval)
            
            if(Condition_DEBUG == 1)
            {
              cat("a.df\n")
              cat(str(a.df))
              cat("\n")
            }
            
            list_PREY[[h]]<-a.df
            
            # ###############################################
            # quit(status = 1)
            
          }#FLAG_diversity_2 >1
        }#FLAG_diversity >1
      }# h in 1 :length(unique(Log_rank_df_sel$CLASS_mega)
      
      Condition_DEBUG <- 0
      
      if(length(list_PREY) >0)
      {
        
        PREY_df = unique(as.data.frame(data.table::rbindlist(list_PREY, fill=T), stringsAsFactors=F))
        
        if(Condition_DEBUG == 1)
        {
          cat("PREY_df_0\n")
          cat(str(PREY_df))
          cat("\n")
          #quit(status = 1)
        }
        
        list_DEF[[k]]<-PREY_df
        
      }# length(list_PREY) >0
      
      Condition_DEBUG <- 0
      
      
    }#k in 1 :length(unique(Log_rank_df_sel$CLASS_mega))
    
    
    Condition_DEBUG <- 0
    
    if(length(list_DEF) >0)
    {
      FINAL_df_per_CT = unique(as.data.frame(data.table::rbindlist(list_DEF, fill=T), stringsAsFactors=F))
      
      FINAL_df_per_CT$Cell_Type<-CT_sel
      
      if(Condition_DEBUG == 1)
      {
        cat("FINAL_df_per_CT_0\n")
        cat(str(FINAL_df_per_CT))
        cat("\n")
        #quit(status = 1)
      }
      
      list_ABSOLUTE_DEF[[i]]<-FINAL_df_per_CT
      
    }#length(list_DEF) >0
    
    Condition_DEBUG <- 0
    
    #### GRAPH -----
    
    Condition_DEBUG <- 1
    
    vector_colors_CLASS_mega<-c("gray","brown")
    
    Freq_table_sel<-Freq_table[which(Freq_table$Cell_Type == CT_sel),]
    
    if(Condition_DEBUG == 1)
    {
      cat("Freq_table_sel_0\n")
      cat(str(Freq_table_sel))
      cat("\n")
      #quit(status = 1)
    }
    
    
   
    
    path5<-paste(out2,'LINEAGE_CUMULATIVE_CURVES','/', sep='')
    
    if (file.exists(path5)){
      
      
      
      
    } else {
      dir.create(file.path(path5))
      
    }
    
    setwd(path5)
    
    v_parameter<-"NA" 
    
    indx_CLASS_mega<-which(colnames(Freq_table_sel) == "CLASS_mega")
    indx.yaxis<-which(colnames(Freq_table_sel) == "Perc")
    mycols <- vector_colors_CLASS_mega
    
    pdfname<-paste(CT_sel,"_CLASS_mega","_Cummulative_frequency_","E_Plus_ASE_activity",".pdf",sep='')
    pdf(file=pdfname, width=5, height=4, pointsize=12)
    
    par(mai=c(0.9,0.9,0.3,0.2))
    lab <- as.character(unique(Freq_table_sel[,indx_CLASS_mega]))
    
    cat("lab\n")
    cat(sprintf(as.character(lab)))
    cat("\n")
    
    
    plot(Freq_table_sel$TILE, Freq_table_sel[,indx.yaxis],
         ty="n", xlab="TILES with E_Plus_ASE activity",
         ylab="Cummulative % of variants",
         axes=F, cex.lab=1.2, cex.lab=1.3, ylim=c(0,100))
    
    
    
    # cat("Hello_world1\n")
    
    # points(Freq_table_sel$TILE, Freq_table_sel[,indx.yaxis], col="darkgrey", pch=19)
    
    # cat("Hello_world2\n")
    # 
    # 
    for (iteration_graph in 1:length(lab))  {
      
      cat("Hello_world3\t")
      cat(sprintf(as.character(lab[iteration_graph])))
      cat("\n")
      
      ind <- which(Freq_table_sel[,indx_CLASS_mega]==lab[iteration_graph])
      
      cat("------------------->ind\n")
      cat(str(ind))
      cat("\n")
      color_sel<-mycols[iteration_graph]
      
      cat(sprintf(as.character(color_sel)))
      cat("\n")
      
      points(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind],
             pch=21,lwd=5, col=color_sel)
      lines(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind],
            lty=1,lwd=5, col=color_sel)
      
      cat("END_3\n")
      
      # points(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind], pch=19, col=color_sel)
    }#iteration_graph in 1:length(lab)
    
    cat("END_4\n")
    
    
    legend("topright", legend=lab, fill=mycols, border=mycols, bty="n")
    axis(1, at=seq(0,100))
    axis(2, las=1)
    
    dev.off()
    
    
    
    list_ABSOLUTE_DEF_Freq[[i]]<-Freq_table_sel
    
    
    Condition_DEBUG <- 0
    
    # ############################################################### HERE HERE
    # quit(status = 1)
    
  }################################ i 1:length(CT_levels)
  
  
  
  
  
  
  
  
  
  Condition_DEBUG <- 1
  
  if(length(list_ABSOLUTE_DEF) >0)
  {
    FINAL_df_ABSOLUTE = unique(as.data.frame(data.table::rbindlist(list_ABSOLUTE_DEF, fill=T), stringsAsFactors=F))
    
   
    if(Condition_DEBUG == 1)
    {
      cat("FINAL_df_ABSOLUTE_0\n")
      cat(str(FINAL_df_ABSOLUTE))
      cat("\n")
      #quit(status = 1)
    }
    
    path5<-paste(out2,'LINEAGE_CUMULATIVE_CURVES','/', sep='')
    
    
    if (file.exists(path5)){
      
      
      
      
    } else {
      dir.create(file.path(path5))
      
    }
    
    setwd(path5)
   

   
    write.table(FINAL_df_ABSOLUTE, file="E_Plus_ASE_CLASS_mega_log_Rank_test_STAT.tsv", sep="\t", row.names = F, quote = F)
    
    
  }#length(list_ABSOLUTE_DEF) >0
  
  if(length(list_ABSOLUTE_DEF_Freq) >0)
  {
    FINAL_ABSOLUTE_DEF_Freq = unique(as.data.frame(data.table::rbindlist(list_ABSOLUTE_DEF_Freq, fill=T), stringsAsFactors=F))
    
    
    if(Condition_DEBUG == 1)
    {
      cat("FINAL_ABSOLUTE_DEF_Freq_0\n")
      cat(str(FINAL_ABSOLUTE_DEF_Freq))
      cat("\n")
      #quit(status = 1)
    }
    
    path5<-paste(out2,'LINEAGE_CUMULATIVE_CURVES','/', sep='')
    
    if (file.exists(path5)){
      
      
      
      
    } else {
      dir.create(file.path(path5))
      
    }
    
    setwd(path5)
    
    write.table(FINAL_ABSOLUTE_DEF_Freq, file="E_Plus_ASE_CLASS_mega_Cummulative_Freq_table.tsv", sep="\t", quote=F,row.names = F)
    
    
  }#length(list_ABSOLUTE_DEF_Freq) >0
  
  
 
}
 
cummulative_CT_CLASS_erythroid_and_Log_Rank_test_E_Plus_ASE = function(option_list)
{
  suppressMessages(library("nph", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ and transform out2 ----
  
  out2 = opt$out2
  
  cat("out2_\n")
  cat(sprintf(as.character(out2)))
  cat("\n")
  
  setwd(out)
  
  
  filename<-paste('df_Cell_colors','.rds',sep='')
  
  df_Cell_colors<-readRDS(file=filename)
  
  cat("df_Cell_colors_0\n")
  cat(str(df_Cell_colors))
  cat("\n")
  
  
  ### CUMMULATIVE_CLASSES #----
  
  path5<-paste(out2,'LINEAGE_CUMULATIVE_CURVES','/', sep='')
  
  cat("path5\n")
  cat(sprintf(as.character(path5)))
  cat("\n")
  
  setwd(path5)
  
  filename_1<-paste("cummulative_classification_ONLY_ASSAYED_with_Lineage_classification",".rds", sep='')
  
  
  CUMMULATIVE_CLASSES_restricted<-readRDS(file=filename_1)
  
  
  cat("CUMMULATIVE_CLASSES_restricted\n")
  cat(str(CUMMULATIVE_CLASSES_restricted))
  cat("\n")
  cat(str(unique(CUMMULATIVE_CLASSES_restricted$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(CUMMULATIVE_CLASSES_restricted$CLASS_erythroid)))))
  cat("\n")
  cat(sprintf(as.character(summary(CUMMULATIVE_CLASSES_restricted$CLASS_erythroid))))
  cat("\n")
  
  Condition_DEBUG <- 1
  
  ########### Build df ---------------
  
  Cell_Type_levels<-levels(CUMMULATIVE_CLASSES_restricted$Cell_Type)
  n_Cell_Type_levels<-length(levels(CUMMULATIVE_CLASSES_restricted$Cell_Type))
  n_MPRA_TILES<-5
  VAR_vector<-as.character(unique(CUMMULATIVE_CLASSES_restricted$VAR))
  KEY_Plus_carried_variants_vector<-as.character(unique(CUMMULATIVE_CLASSES_restricted$KEY_Plus_carried_variants))
  
  
  
  n_VAR<-length(VAR_vector)
  
  if(Condition_DEBUG == 1)
  {
    cat("Parameters\n")
    cat(str(Cell_Type_levels))
    cat("\n")
    cat(str(n_Cell_Type_levels))
    cat("\n")
    cat(str(n_MPRA_TILES))
    cat("\n")
    cat(str(VAR_vector))
    cat("\n")
    cat(str(KEY_Plus_carried_variants_vector))
    cat("\n")
    
    
  }
  
  
  #### Build LogRank matrix ----
  
  
  
  Log_rank_df<-as.data.frame(cbind(rep(VAR_vector,n_MPRA_TILES),
                                   rep(KEY_Plus_carried_variants_vector,n_MPRA_TILES),
                                   c(rep(1,n_VAR),rep(2,n_VAR),rep(3,n_VAR),rep(4,n_VAR),rep(5,n_VAR)),
                                   rep("NA",n_MPRA_TILES*n_VAR)),
                             stringsAsFactors=F)
  
  colnames(Log_rank_df)<-c("VAR","KEY_Plus_carried_variants","TILE","ACTIVE")
  
  
  if(Condition_DEBUG == 1)
  {
    cat("Log_rank_df_0\n")
    cat(str(Log_rank_df))
    cat("\n")
    cat(str(unique(Log_rank_df$VAR)))
    cat("\n")
    cat(str(unique(Log_rank_df$KEY_Plus_carried_variants)))
    cat("\n")
  }
  
  
  
  CUMMULATIVE_CLASSES_restricted_subset<-unique(CUMMULATIVE_CLASSES_restricted[,c(which(colnames(CUMMULATIVE_CLASSES_restricted) == "VAR"),
                                                                                  which(colnames(CUMMULATIVE_CLASSES_restricted) == "KEY_Plus_carried_variants"),
                                                                                  which(colnames(CUMMULATIVE_CLASSES_restricted) == "CLASS_erythroid"),
                                                                                  which(colnames(CUMMULATIVE_CLASSES_restricted) == "Cell_Type"),
                                                                                  which(colnames(CUMMULATIVE_CLASSES_restricted) == "E_Plus_ASE_CLASS_TILES"))])
  levels_CLASS_erythroid<-levels(CUMMULATIVE_CLASSES_restricted$CLASS_erythroid)
  
  if(Condition_DEBUG == 1)
  {
    cat("CUMMULATIVE_CLASSES_restricted_subset_0\n")
    cat(str(CUMMULATIVE_CLASSES_restricted_subset))
    cat("\n")
    cat(str(unique(CUMMULATIVE_CLASSES_restricted_subset$VAR)))
    cat("\n")
    cat(str(unique(CUMMULATIVE_CLASSES_restricted_subset$KEY_Plus_carried_variants)))
    cat("\n")
  }
  
  Log_rank_df<-merge(Log_rank_df,
                     CUMMULATIVE_CLASSES_restricted_subset,
                     by=c("VAR","KEY_Plus_carried_variants"),
                     all=T)
  
  Log_rank_df$CLASS_erythroid<-Log_rank_df$CLASS_erythroid
  
  Log_rank_df$CLASS_erythroid<-factor(Log_rank_df$CLASS_erythroid,
                                 levels=levels_CLASS_erythroid,
                                 ordered=T)
  
  Log_rank_df$Cell_Type<-factor(Log_rank_df$Cell_Type,
                                levels=Cell_Type_levels,
                                ordered=T)
  
  Log_rank_df$TILE<-as.integer(Log_rank_df$TILE)
  
  if(Condition_DEBUG == 1)
  {
    cat("Log_rank_df_1\n")
    cat(str(Log_rank_df))
    cat("\n")
    cat(str(unique(Log_rank_df$VAR)))
    cat("\n")
    cat(str(unique(Log_rank_df$KEY_Plus_carried_variants)))
    cat("\n")
  }
  
  
  
  
  
  ################# Classification of TILES ----------------
  
  Log_rank_df$ACTIVE[which(Log_rank_df$TILE <= Log_rank_df$E_Plus_ASE_CLASS_TILES)]<-"1"
  Log_rank_df$ACTIVE[which(Log_rank_df$TILE > Log_rank_df$E_Plus_ASE_CLASS_TILES)]<-"0"
  
  Log_rank_df$ACTIVE<-as.numeric(Log_rank_df$ACTIVE)
  
  
  
  if(Condition_DEBUG == 1)
  {
    cat("Log_rank_df_2\n")
    cat(str(Log_rank_df))
    cat("\n")
    cat(str(unique(Log_rank_df$VAR)))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(Log_rank_df$ACTIVE))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(Log_rank_df$ACTIVE)))))
    cat("\n")
  }
  
  Log_rank_df<-Log_rank_df[order(Log_rank_df$KEY_Plus_carried_variants,Log_rank_df$Cell_Type, Log_rank_df$TILE),]
  
  
  setwd(out2)
  
  write.table(Log_rank_df, file="check2.tsv",sep="\t",quote=F,row.names = F)
  
  # ################################################################
  # quit(status = 1)
  
  
  ############## Freq Table ----------------------------
  
  ### TOTAL ---
  
  CUMMULATIVE_CLASSES_restricted.dt<-data.table(CUMMULATIVE_CLASSES_restricted, key=c("Cell_Type","CLASS_erythroid"))
  
  Freq_TOTAL<-as.data.frame(CUMMULATIVE_CLASSES_restricted.dt[,.(TOTAL=.N),by=key(CUMMULATIVE_CLASSES_restricted.dt)], stringsAsFactors=F)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_TOTAL_0\n")
    cat(str(Freq_TOTAL))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_TOTAL.dt<-data.table(Freq_TOTAL, key=c("CLASS_erythroid"))
  
  
  MAXED_df<-as.data.frame(Freq_TOTAL.dt[,.(MAX=max(TOTAL)),by=key(Freq_TOTAL.dt)], stringsAsFactors=F)
  
  
  
  if(Condition_DEBUG == 1)
  {
    cat("MAXED_df_0\n")
    cat(str(MAXED_df))
    cat("\n")
    #quit(status = 1)
  }
  
  ########## Active _tiles --------------------
  
  Log_rank_df.dt<-data.table(Log_rank_df, key=c("Cell_Type","CLASS_erythroid","TILE"))
  
  Freq_table<-as.data.frame(Log_rank_df.dt[,.(Freq=sum(ACTIVE)),by=key(Log_rank_df.dt)], stringsAsFactors=F)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_table_0\n")
    cat(str(Freq_table))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_table<-merge(MAXED_df,
                    Freq_table,
                    by=c("CLASS_erythroid"),
                    all=T)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_table_1\n")
    cat(str(Freq_table))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_table$Perc<-round(100*(Freq_table$Freq/Freq_table$MAX),2)
  
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_table_2\n")
    cat(str(Freq_table))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_table$CLASS_erythroid<-factor(Freq_table$CLASS_erythroid,
                                levels=levels_CLASS_erythroid,
                                ordered=T)
  
  Freq_table$Cell_Type<-factor(Freq_table$Cell_Type,
                               levels=Cell_Type_levels,
                               ordered=T)
  
  Freq_table$TILE<-as.integer(Freq_table$TILE)
  
  Freq_table<-Freq_table[order(Freq_table$CLASS_erythroid,Freq_table$Cell_Type, Freq_table$TILE),]
  
  
  setwd(out2)
  
  write.table(Freq_table, file="check3.tsv", row.names = F,quote = F,sep="\t")
  
  # ####################################################################
  # quit(status = 1)
  
  ############### Per CT LOOP ----
  
  
  CT_levels<-Cell_Type_levels
  
  list_ABSOLUTE_DEF<-list()
  list_ABSOLUTE_DEF_Freq<-list()
  
  Condition_DEBUG <- 1
  
  for(i in 1:length(CT_levels))
  {
    CT_sel<-CT_levels[i]
    
    cat("CT:---------------------------------------->\t")
    cat(sprintf(as.character(CT_sel)))
    cat("\n")
    
    Log_rank_df_sel<-Log_rank_df[which(Log_rank_df$Cell_Type == CT_sel),]
    
    if(Condition_DEBUG == 1)
    {
      cat("Log_rank_df_sel_0\n")
      cat(str(Log_rank_df_sel))
      cat("\n")
      cat(str(unique(Log_rank_df_sel$VAR)))
      cat("\n")
    }
    
    
    #### Log Rank test ----
    
    Condition_DEBUG <- 0
    
    list_DEF<-list()
    
    CLASS_erythroid_labels<-levels(Log_rank_df_sel$CLASS_erythroid)
    
    for(k in 1 :length(CLASS_erythroid_labels))
    {
      BAIT<-CLASS_erythroid_labels[k]
      
      cat("-----BAIT------>\t")
      cat(sprintf(as.character(BAIT)))
      cat("\t")
      
      list_PREY<-list()
      
      for(h in 1 :length(CLASS_erythroid_labels))
      {
        
        PREY<-CLASS_erythroid_labels[h]
        
        cat("------PREY----->\t")
        cat(sprintf(as.character(PREY)))
        cat("\n")
        
        Mini_comparison<-unique(c(BAIT,PREY))
        
        if(Condition_DEBUG == 1)
        {
          cat("Mini_comparison_0\n")
          cat(str(Mini_comparison))
          cat("\n")
        }
        
        Log_rank_df_sel_subset<-Log_rank_df_sel[which(Log_rank_df_sel$CLASS_erythroid%in%Mini_comparison),]
        
        if(Condition_DEBUG == 1)
        {
          cat("Log_rank_df_sel_subset_0\n")
          cat(str(Log_rank_df_sel_subset))
          cat("\n")
        }
        
        FLAG_diversity<-as.integer(length(unique(Log_rank_df_sel_subset$CLASS_erythroid)))
        
        if(Condition_DEBUG == 1)
        {
          cat("FLAG_diversity_0\n")
          cat(sprintf(as.character(FLAG_diversity)))
          cat("\n")
          
          cat(str(FLAG_diversity))
          cat("\n")
        }
        
        if(FLAG_diversity >1)
        {
          cat("Hello_world\n")
          
          Mini_comparison_c1<-as.character(Mini_comparison[1])
          Mini_comparison_c2<-as.character(Mini_comparison[2])
          
          FLAG_diversity_2<-length(unique(Log_rank_df_sel_subset$ACTIVE))
          
          if(Condition_DEBUG == 1)
          {
            cat("FLAG_diversity_2_0\n")
            cat(str(FLAG_diversity_2))
            cat("\n")
          }
          
          if(FLAG_diversity_2 >1)
          {
            if(Condition_DEBUG == 1)
            {
              cat("----------->\t")
              cat(sprintf(as.character(Mini_comparison_c1)))
              cat("\t")
              cat(sprintf(as.character(Mini_comparison_c2)))
              cat("\n")
              
              cat("Log_rank_df_sel_subset_0\n")
              cat(str(Log_rank_df_sel_subset))
              cat("\n")
              cat(str(unique(Log_rank_df_sel_subset$VAR)))
              cat("\n")
              cat(sprintf(as.character(names(summary(as.factor(Log_rank_df_sel_subset$ACTIVE))))))
              cat("\n")
              cat(sprintf(as.character(summary(as.factor(Log_rank_df_sel_subset$ACTIVE)))))
              cat("\n")
            }
            
            #, "less", "greater"
            
            STATS_test<-logrank.test(
              Log_rank_df_sel_subset$TILE,
              Log_rank_df_sel_subset$ACTIVE,
              Log_rank_df_sel_subset$CLASS_erythroid,
              alternative = c("two.sided"),
              rho = 0,
              gamma = 0,
              event_time_weights = NULL
            )
            
            # if(Condition_DEBUG == 1)
            # {
            #   cat("STATS_test_0\n")
            #   cat(str(STATS_test))
            #   cat("\n")
            # }
            
            
            p_value<-STATS_test$test$p
            minus_log_val<-round(-1*log10(p_value),2)
            
            if(Condition_DEBUG == 1)
            {
              cat("minus_log_val_0\n")
              cat(str(minus_log_val))
              cat("\n")
            }
            
            
            a.df<-as.data.frame(cbind(Mini_comparison_c1,Mini_comparison_c2,minus_log_val), stringsAsFactors=F)
            
            
            colnames(a.df)<-c("CLASS_erythroid_c1","CLASS_erythroid_c2","minus_logpval")
            a.df$minus_logpval<-as.numeric(a.df$minus_logpval)
            
            if(Condition_DEBUG == 1)
            {
              cat("a.df\n")
              cat(str(a.df))
              cat("\n")
            }
            
            list_PREY[[h]]<-a.df
            
            # ###############################################
            # quit(status = 1)
            
          }#FLAG_diversity_2 >1
        }#FLAG_diversity >1
      }# h in 1 :length(unique(Log_rank_df_sel$CLASS_erythroid)
      
      Condition_DEBUG <- 0
      
      if(length(list_PREY) >0)
      {
        
        PREY_df = unique(as.data.frame(data.table::rbindlist(list_PREY, fill=T), stringsAsFactors=F))
        
        if(Condition_DEBUG == 1)
        {
          cat("PREY_df_0\n")
          cat(str(PREY_df))
          cat("\n")
          #quit(status = 1)
        }
        
        list_DEF[[k]]<-PREY_df
        
      }# length(list_PREY) >0
      
      Condition_DEBUG <- 0
      
      
    }#k in 1 :length(unique(Log_rank_df_sel$CLASS_erythroid))
    
    
    Condition_DEBUG <- 0
    
    if(length(list_DEF) >0)
    {
      FINAL_df_per_CT = unique(as.data.frame(data.table::rbindlist(list_DEF, fill=T), stringsAsFactors=F))
      
      FINAL_df_per_CT$Cell_Type<-CT_sel
      
      if(Condition_DEBUG == 1)
      {
        cat("FINAL_df_per_CT_0\n")
        cat(str(FINAL_df_per_CT))
        cat("\n")
        #quit(status = 1)
      }
      
      list_ABSOLUTE_DEF[[i]]<-FINAL_df_per_CT
      
    }#length(list_DEF) >0
    
    Condition_DEBUG <- 0
    
    #### GRAPH -----
    
    Condition_DEBUG <- 1
    
    vector_colors_CLASS_erythroid<-c("gray","firebrick1")
    
    Freq_table_sel<-Freq_table[which(Freq_table$Cell_Type == CT_sel),]
    
    if(Condition_DEBUG == 1)
    {
      cat("Freq_table_sel_0\n")
      cat(str(Freq_table_sel))
      cat("\n")
      #quit(status = 1)
    }
    
    
    
    
    path5<-paste(out2,'LINEAGE_CUMULATIVE_CURVES','/', sep='')
    
    if (file.exists(path5)){
      
      
      
      
    } else {
      dir.create(file.path(path5))
      
    }
    
    setwd(path5)
    
    v_parameter<-"NA" 
    
    indx_CLASS_erythroid<-which(colnames(Freq_table_sel) == "CLASS_erythroid")
    indx.yaxis<-which(colnames(Freq_table_sel) == "Perc")
    mycols <- vector_colors_CLASS_erythroid
    
    pdfname<-paste(CT_sel,"_CLASS_erythroid","_Cummulative_frequency_","E_Plus_ASE_activity",".pdf",sep='')
    pdf(file=pdfname, width=5, height=4, pointsize=12)
    
    par(mai=c(0.9,0.9,0.3,0.2))
    lab <- as.character(unique(Freq_table_sel[,indx_CLASS_erythroid]))
    
    cat("lab\n")
    cat(sprintf(as.character(lab)))
    cat("\n")
    
    
    plot(Freq_table_sel$TILE, Freq_table_sel[,indx.yaxis],
         ty="n", xlab="TILES with E_Plus_ASE activity",
         ylab="Cummulative % of variants",
         axes=F, cex.lab=1.2, cex.lab=1.3, ylim=c(0,100))
    
    
    
    # cat("Hello_world1\n")
    
    # points(Freq_table_sel$TILE, Freq_table_sel[,indx.yaxis], col="darkgrey", pch=19)
    
    # cat("Hello_world2\n")
    # 
    # 
    for (iteration_graph in 1:length(lab))  {
      
      cat("Hello_world3\t")
      cat(sprintf(as.character(lab[iteration_graph])))
      cat("\n")
      
      ind <- which(Freq_table_sel[,indx_CLASS_erythroid]==lab[iteration_graph])
      
      cat("------------------->ind\n")
      cat(str(ind))
      cat("\n")
      color_sel<-mycols[iteration_graph]
      
      cat(sprintf(as.character(color_sel)))
      cat("\n")
      
      points(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind],
             pch=21,lwd=5, col=color_sel)
      lines(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind],
            lty=1,lwd=5, col=color_sel)
      
      cat("END_3\n")
      
      # points(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind], pch=19, col=color_sel)
    }#iteration_graph in 1:length(lab)
    
    cat("END_4\n")
    
    
    legend("topright", legend=lab, fill=mycols, border=mycols, bty="n")
    axis(1, at=seq(0,100))
    axis(2, las=1)
    
    dev.off()
    
    
    
    list_ABSOLUTE_DEF_Freq[[i]]<-Freq_table_sel
    
    
    Condition_DEBUG <- 0
    
    # ############################################################### HERE HERE
    # quit(status = 1)
    
  }################################ i 1:length(CT_levels)
  
  
  
  
  
  
  
  
  
  Condition_DEBUG <- 1
  
  if(length(list_ABSOLUTE_DEF) >0)
  {
    FINAL_df_ABSOLUTE = unique(as.data.frame(data.table::rbindlist(list_ABSOLUTE_DEF, fill=T), stringsAsFactors=F))
    
    
    if(Condition_DEBUG == 1)
    {
      cat("FINAL_df_ABSOLUTE_0\n")
      cat(str(FINAL_df_ABSOLUTE))
      cat("\n")
      #quit(status = 1)
    }
    
    path5<-paste(out2,'LINEAGE_CUMULATIVE_CURVES','/', sep='')
    
    
    if (file.exists(path5)){
      
      
      
      
    } else {
      dir.create(file.path(path5))
      
    }
    
    setwd(path5)
    
    
    
    write.table(FINAL_df_ABSOLUTE, file="E_Plus_ASE_CLASS_erythroid_log_Rank_test_STAT.tsv", sep="\t", row.names = F, quote = F)
    
    
  }#length(list_ABSOLUTE_DEF) >0
  
  if(length(list_ABSOLUTE_DEF_Freq) >0)
  {
    FINAL_ABSOLUTE_DEF_Freq = unique(as.data.frame(data.table::rbindlist(list_ABSOLUTE_DEF_Freq, fill=T), stringsAsFactors=F))
    
    
    if(Condition_DEBUG == 1)
    {
      cat("FINAL_ABSOLUTE_DEF_Freq_0\n")
      cat(str(FINAL_ABSOLUTE_DEF_Freq))
      cat("\n")
      #quit(status = 1)
    }
    
    path5<-paste(out2,'LINEAGE_CUMULATIVE_CURVES','/', sep='')
    
    if (file.exists(path5)){
      
      
      
      
    } else {
      dir.create(file.path(path5))
      
    }
    
    setwd(path5)
    
    write.table(FINAL_ABSOLUTE_DEF_Freq, file="E_Plus_ASE_CLASS_erythroid_Cummulative_Freq_table.tsv", sep="\t", quote=F,row.names = F)
    
    
  }#length(list_ABSOLUTE_DEF_Freq) >0
  
  
  
}

cummulative_CT_CLASS_gran_mono_and_Log_Rank_test_E_Plus_ASE = function(option_list)
{
  suppressMessages(library("nph", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ and transform out2 ----
  
  out2 = opt$out2
  
  cat("out2_\n")
  cat(sprintf(as.character(out2)))
  cat("\n")
  
  setwd(out)
  
  
  filename<-paste('df_Cell_colors','.rds',sep='')
  
  df_Cell_colors<-readRDS(file=filename)
  
  cat("df_Cell_colors_0\n")
  cat(str(df_Cell_colors))
  cat("\n")
  
  
  ### CUMMULATIVE_CLASSES #----
  
  path5<-paste(out2,'LINEAGE_CUMULATIVE_CURVES','/', sep='')
  
  cat("path5\n")
  cat(sprintf(as.character(path5)))
  cat("\n")
  
  setwd(path5)
  
  filename_1<-paste("cummulative_classification_ONLY_ASSAYED_with_Lineage_classification",".rds", sep='')
  
  
  CUMMULATIVE_CLASSES_restricted<-readRDS(file=filename_1)
  
  
  cat("CUMMULATIVE_CLASSES_restricted\n")
  cat(str(CUMMULATIVE_CLASSES_restricted))
  cat("\n")
  cat(str(unique(CUMMULATIVE_CLASSES_restricted$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(CUMMULATIVE_CLASSES_restricted$CLASS_gran_mono)))))
  cat("\n")
  cat(sprintf(as.character(summary(CUMMULATIVE_CLASSES_restricted$CLASS_gran_mono))))
  cat("\n")
  
  Condition_DEBUG <- 1
  
  ########### Build df ---------------
  
  Cell_Type_levels<-levels(CUMMULATIVE_CLASSES_restricted$Cell_Type)
  n_Cell_Type_levels<-length(levels(CUMMULATIVE_CLASSES_restricted$Cell_Type))
  n_MPRA_TILES<-5
  VAR_vector<-as.character(unique(CUMMULATIVE_CLASSES_restricted$VAR))
  KEY_Plus_carried_variants_vector<-as.character(unique(CUMMULATIVE_CLASSES_restricted$KEY_Plus_carried_variants))
  
  
  
  n_VAR<-length(VAR_vector)
  
  if(Condition_DEBUG == 1)
  {
    cat("Parameters\n")
    cat(str(Cell_Type_levels))
    cat("\n")
    cat(str(n_Cell_Type_levels))
    cat("\n")
    cat(str(n_MPRA_TILES))
    cat("\n")
    cat(str(VAR_vector))
    cat("\n")
    cat(str(KEY_Plus_carried_variants_vector))
    cat("\n")
    
    
  }
  
  
  #### Build LogRank matrix ----
  
  
  
  Log_rank_df<-as.data.frame(cbind(rep(VAR_vector,n_MPRA_TILES),
                                   rep(KEY_Plus_carried_variants_vector,n_MPRA_TILES),
                                   c(rep(1,n_VAR),rep(2,n_VAR),rep(3,n_VAR),rep(4,n_VAR),rep(5,n_VAR)),
                                   rep("NA",n_MPRA_TILES*n_VAR)),
                             stringsAsFactors=F)
  
  colnames(Log_rank_df)<-c("VAR","KEY_Plus_carried_variants","TILE","ACTIVE")
  
  
  if(Condition_DEBUG == 1)
  {
    cat("Log_rank_df_0\n")
    cat(str(Log_rank_df))
    cat("\n")
    cat(str(unique(Log_rank_df$VAR)))
    cat("\n")
    cat(str(unique(Log_rank_df$KEY_Plus_carried_variants)))
    cat("\n")
  }
  
  
  
  CUMMULATIVE_CLASSES_restricted_subset<-unique(CUMMULATIVE_CLASSES_restricted[,c(which(colnames(CUMMULATIVE_CLASSES_restricted) == "VAR"),
                                                                                  which(colnames(CUMMULATIVE_CLASSES_restricted) == "KEY_Plus_carried_variants"),
                                                                                  which(colnames(CUMMULATIVE_CLASSES_restricted) == "CLASS_gran_mono"),
                                                                                  which(colnames(CUMMULATIVE_CLASSES_restricted) == "Cell_Type"),
                                                                                  which(colnames(CUMMULATIVE_CLASSES_restricted) == "E_Plus_ASE_CLASS_TILES"))])
  levels_CLASS_gran_mono<-levels(CUMMULATIVE_CLASSES_restricted$CLASS_gran_mono)
  
  if(Condition_DEBUG == 1)
  {
    cat("CUMMULATIVE_CLASSES_restricted_subset_0\n")
    cat(str(CUMMULATIVE_CLASSES_restricted_subset))
    cat("\n")
    cat(str(unique(CUMMULATIVE_CLASSES_restricted_subset$VAR)))
    cat("\n")
    cat(str(unique(CUMMULATIVE_CLASSES_restricted_subset$KEY_Plus_carried_variants)))
    cat("\n")
  }
  
  Log_rank_df<-merge(Log_rank_df,
                     CUMMULATIVE_CLASSES_restricted_subset,
                     by=c("VAR","KEY_Plus_carried_variants"),
                     all=T)
  
  Log_rank_df$CLASS_gran_mono<-Log_rank_df$CLASS_gran_mono
  
  Log_rank_df$CLASS_gran_mono<-factor(Log_rank_df$CLASS_gran_mono,
                                      levels=levels_CLASS_gran_mono,
                                      ordered=T)
  
  Log_rank_df$Cell_Type<-factor(Log_rank_df$Cell_Type,
                                levels=Cell_Type_levels,
                                ordered=T)
  
  Log_rank_df$TILE<-as.integer(Log_rank_df$TILE)
  
  if(Condition_DEBUG == 1)
  {
    cat("Log_rank_df_1\n")
    cat(str(Log_rank_df))
    cat("\n")
    cat(str(unique(Log_rank_df$VAR)))
    cat("\n")
    cat(str(unique(Log_rank_df$KEY_Plus_carried_variants)))
    cat("\n")
  }
  
  
  
  
  
  ################# Classification of TILES ----------------
  
  Log_rank_df$ACTIVE[which(Log_rank_df$TILE <= Log_rank_df$E_Plus_ASE_CLASS_TILES)]<-"1"
  Log_rank_df$ACTIVE[which(Log_rank_df$TILE > Log_rank_df$E_Plus_ASE_CLASS_TILES)]<-"0"
  
  Log_rank_df$ACTIVE<-as.numeric(Log_rank_df$ACTIVE)
  
  
  
  if(Condition_DEBUG == 1)
  {
    cat("Log_rank_df_2\n")
    cat(str(Log_rank_df))
    cat("\n")
    cat(str(unique(Log_rank_df$VAR)))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(Log_rank_df$ACTIVE))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(Log_rank_df$ACTIVE)))))
    cat("\n")
  }
  
  Log_rank_df<-Log_rank_df[order(Log_rank_df$KEY_Plus_carried_variants,Log_rank_df$Cell_Type, Log_rank_df$TILE),]
  
  
  setwd(out2)
  
  write.table(Log_rank_df, file="check2.tsv",sep="\t",quote=F,row.names = F)
  
  # ################################################################
  # quit(status = 1)
  
  
  ############## Freq Table ----------------------------
  
  ### TOTAL ---
  
  CUMMULATIVE_CLASSES_restricted.dt<-data.table(CUMMULATIVE_CLASSES_restricted, key=c("Cell_Type","CLASS_gran_mono"))
  
  Freq_TOTAL<-as.data.frame(CUMMULATIVE_CLASSES_restricted.dt[,.(TOTAL=.N),by=key(CUMMULATIVE_CLASSES_restricted.dt)], stringsAsFactors=F)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_TOTAL_0\n")
    cat(str(Freq_TOTAL))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_TOTAL.dt<-data.table(Freq_TOTAL, key=c("CLASS_gran_mono"))
  
  
  MAXED_df<-as.data.frame(Freq_TOTAL.dt[,.(MAX=max(TOTAL)),by=key(Freq_TOTAL.dt)], stringsAsFactors=F)
  
  
  
  if(Condition_DEBUG == 1)
  {
    cat("MAXED_df_0\n")
    cat(str(MAXED_df))
    cat("\n")
    #quit(status = 1)
  }
  
  ########## Active _tiles --------------------
  
  Log_rank_df.dt<-data.table(Log_rank_df, key=c("Cell_Type","CLASS_gran_mono","TILE"))
  
  Freq_table<-as.data.frame(Log_rank_df.dt[,.(Freq=sum(ACTIVE)),by=key(Log_rank_df.dt)], stringsAsFactors=F)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_table_0\n")
    cat(str(Freq_table))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_table<-merge(MAXED_df,
                    Freq_table,
                    by=c("CLASS_gran_mono"),
                    all=T)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_table_1\n")
    cat(str(Freq_table))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_table$Perc<-round(100*(Freq_table$Freq/Freq_table$MAX),2)
  
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_table_2\n")
    cat(str(Freq_table))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_table$CLASS_gran_mono<-factor(Freq_table$CLASS_gran_mono,
                                     levels=levels_CLASS_gran_mono,
                                     ordered=T)
  
  Freq_table$Cell_Type<-factor(Freq_table$Cell_Type,
                               levels=Cell_Type_levels,
                               ordered=T)
  
  Freq_table$TILE<-as.integer(Freq_table$TILE)
  
  Freq_table<-Freq_table[order(Freq_table$CLASS_gran_mono,Freq_table$Cell_Type, Freq_table$TILE),]
  
  
  setwd(out2)
  
  write.table(Freq_table, file="check3.tsv", row.names = F,quote = F,sep="\t")
  
  # ####################################################################
  # quit(status = 1)
  
  ############### Per CT LOOP ----
  
  
  CT_levels<-Cell_Type_levels
  
  list_ABSOLUTE_DEF<-list()
  list_ABSOLUTE_DEF_Freq<-list()
  
  Condition_DEBUG <- 1
  
  for(i in 1:length(CT_levels))
  {
    CT_sel<-CT_levels[i]
    
    cat("CT:---------------------------------------->\t")
    cat(sprintf(as.character(CT_sel)))
    cat("\n")
    
    Log_rank_df_sel<-Log_rank_df[which(Log_rank_df$Cell_Type == CT_sel),]
    
    if(Condition_DEBUG == 1)
    {
      cat("Log_rank_df_sel_0\n")
      cat(str(Log_rank_df_sel))
      cat("\n")
      cat(str(unique(Log_rank_df_sel$VAR)))
      cat("\n")
    }
    
    
    #### Log Rank test ----
    
    Condition_DEBUG <- 0
    
    list_DEF<-list()
    
    CLASS_gran_mono_labels<-levels(Log_rank_df_sel$CLASS_gran_mono)
    
    for(k in 1 :length(CLASS_gran_mono_labels))
    {
      BAIT<-CLASS_gran_mono_labels[k]
      
      cat("-----BAIT------>\t")
      cat(sprintf(as.character(BAIT)))
      cat("\t")
      
      list_PREY<-list()
      
      for(h in 1 :length(CLASS_gran_mono_labels))
      {
        
        PREY<-CLASS_gran_mono_labels[h]
        
        cat("------PREY----->\t")
        cat(sprintf(as.character(PREY)))
        cat("\n")
        
        Mini_comparison<-unique(c(BAIT,PREY))
        
        if(Condition_DEBUG == 1)
        {
          cat("Mini_comparison_0\n")
          cat(str(Mini_comparison))
          cat("\n")
        }
        
        Log_rank_df_sel_subset<-Log_rank_df_sel[which(Log_rank_df_sel$CLASS_gran_mono%in%Mini_comparison),]
        
        if(Condition_DEBUG == 1)
        {
          cat("Log_rank_df_sel_subset_0\n")
          cat(str(Log_rank_df_sel_subset))
          cat("\n")
        }
        
        FLAG_diversity<-as.integer(length(unique(Log_rank_df_sel_subset$CLASS_gran_mono)))
        
        if(Condition_DEBUG == 1)
        {
          cat("FLAG_diversity_0\n")
          cat(sprintf(as.character(FLAG_diversity)))
          cat("\n")
          
          cat(str(FLAG_diversity))
          cat("\n")
        }
        
        if(FLAG_diversity >1)
        {
          cat("Hello_world\n")
          
          Mini_comparison_c1<-as.character(Mini_comparison[1])
          Mini_comparison_c2<-as.character(Mini_comparison[2])
          
          FLAG_diversity_2<-length(unique(Log_rank_df_sel_subset$ACTIVE))
          
          if(Condition_DEBUG == 1)
          {
            cat("FLAG_diversity_2_0\n")
            cat(str(FLAG_diversity_2))
            cat("\n")
          }
          
          if(FLAG_diversity_2 >1)
          {
            if(Condition_DEBUG == 1)
            {
              cat("----------->\t")
              cat(sprintf(as.character(Mini_comparison_c1)))
              cat("\t")
              cat(sprintf(as.character(Mini_comparison_c2)))
              cat("\n")
              
              cat("Log_rank_df_sel_subset_0\n")
              cat(str(Log_rank_df_sel_subset))
              cat("\n")
              cat(str(unique(Log_rank_df_sel_subset$VAR)))
              cat("\n")
              cat(sprintf(as.character(names(summary(as.factor(Log_rank_df_sel_subset$ACTIVE))))))
              cat("\n")
              cat(sprintf(as.character(summary(as.factor(Log_rank_df_sel_subset$ACTIVE)))))
              cat("\n")
            }
            
            #, "less", "greater"
            
            STATS_test<-logrank.test(
              Log_rank_df_sel_subset$TILE,
              Log_rank_df_sel_subset$ACTIVE,
              Log_rank_df_sel_subset$CLASS_gran_mono,
              alternative = c("two.sided"),
              rho = 0,
              gamma = 0,
              event_time_weights = NULL
            )
            
            # if(Condition_DEBUG == 1)
            # {
            #   cat("STATS_test_0\n")
            #   cat(str(STATS_test))
            #   cat("\n")
            # }
            
            
            p_value<-STATS_test$test$p
            minus_log_val<-round(-1*log10(p_value),2)
            
            if(Condition_DEBUG == 1)
            {
              cat("minus_log_val_0\n")
              cat(str(minus_log_val))
              cat("\n")
            }
            
            
            a.df<-as.data.frame(cbind(Mini_comparison_c1,Mini_comparison_c2,minus_log_val), stringsAsFactors=F)
            
            
            colnames(a.df)<-c("CLASS_gran_mono_c1","CLASS_gran_mono_c2","minus_logpval")
            a.df$minus_logpval<-as.numeric(a.df$minus_logpval)
            
            if(Condition_DEBUG == 1)
            {
              cat("a.df\n")
              cat(str(a.df))
              cat("\n")
            }
            
            list_PREY[[h]]<-a.df
            
            # ###############################################
            # quit(status = 1)
            
          }#FLAG_diversity_2 >1
        }#FLAG_diversity >1
      }# h in 1 :length(unique(Log_rank_df_sel$CLASS_gran_mono)
      
      Condition_DEBUG <- 0
      
      if(length(list_PREY) >0)
      {
        
        PREY_df = unique(as.data.frame(data.table::rbindlist(list_PREY, fill=T), stringsAsFactors=F))
        
        if(Condition_DEBUG == 1)
        {
          cat("PREY_df_0\n")
          cat(str(PREY_df))
          cat("\n")
          #quit(status = 1)
        }
        
        list_DEF[[k]]<-PREY_df
        
      }# length(list_PREY) >0
      
      Condition_DEBUG <- 0
      
      
    }#k in 1 :length(unique(Log_rank_df_sel$CLASS_gran_mono))
    
    
    Condition_DEBUG <- 0
    
    if(length(list_DEF) >0)
    {
      FINAL_df_per_CT = unique(as.data.frame(data.table::rbindlist(list_DEF, fill=T), stringsAsFactors=F))
      
      FINAL_df_per_CT$Cell_Type<-CT_sel
      
      if(Condition_DEBUG == 1)
      {
        cat("FINAL_df_per_CT_0\n")
        cat(str(FINAL_df_per_CT))
        cat("\n")
        #quit(status = 1)
      }
      
      list_ABSOLUTE_DEF[[i]]<-FINAL_df_per_CT
      
    }#length(list_DEF) >0
    
    Condition_DEBUG <- 0
    
    #### GRAPH -----
    
    Condition_DEBUG <- 1
    
    vector_colors_CLASS_gran_mono<-c("gray","magenta")
    
    Freq_table_sel<-Freq_table[which(Freq_table$Cell_Type == CT_sel),]
    
    if(Condition_DEBUG == 1)
    {
      cat("Freq_table_sel_0\n")
      cat(str(Freq_table_sel))
      cat("\n")
      #quit(status = 1)
    }
    
    
    
    
    path5<-paste(out2,'LINEAGE_CUMULATIVE_CURVES','/', sep='')
    
    if (file.exists(path5)){
      
      
      
      
    } else {
      dir.create(file.path(path5))
      
    }
    
    setwd(path5)
    
    v_parameter<-"NA" 
    
    indx_CLASS_gran_mono<-which(colnames(Freq_table_sel) == "CLASS_gran_mono")
    indx.yaxis<-which(colnames(Freq_table_sel) == "Perc")
    mycols <- vector_colors_CLASS_gran_mono
    
    pdfname<-paste(CT_sel,"_CLASS_gran_mono","_Cummulative_frequency_","E_Plus_ASE_activity",".pdf",sep='')
    pdf(file=pdfname, width=5, height=4, pointsize=12)
    
    par(mai=c(0.9,0.9,0.3,0.2))
    lab <- as.character(unique(Freq_table_sel[,indx_CLASS_gran_mono]))
    
    cat("lab\n")
    cat(sprintf(as.character(lab)))
    cat("\n")
    
    
    plot(Freq_table_sel$TILE, Freq_table_sel[,indx.yaxis],
         ty="n", xlab="TILES with E_Plus_ASE activity",
         ylab="Cummulative % of variants",
         axes=F, cex.lab=1.2, cex.lab=1.3, ylim=c(0,100))
    
    
    
    # cat("Hello_world1\n")
    
    # points(Freq_table_sel$TILE, Freq_table_sel[,indx.yaxis], col="darkgrey", pch=19)
    
    # cat("Hello_world2\n")
    # 
    # 
    for (iteration_graph in 1:length(lab))  {
      
      cat("Hello_world3\t")
      cat(sprintf(as.character(lab[iteration_graph])))
      cat("\n")
      
      ind <- which(Freq_table_sel[,indx_CLASS_gran_mono]==lab[iteration_graph])
      
      cat("------------------->ind\n")
      cat(str(ind))
      cat("\n")
      color_sel<-mycols[iteration_graph]
      
      cat(sprintf(as.character(color_sel)))
      cat("\n")
      
      points(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind],
             pch=21,lwd=5, col=color_sel)
      lines(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind],
            lty=1,lwd=5, col=color_sel)
      
      cat("END_3\n")
      
      # points(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind], pch=19, col=color_sel)
    }#iteration_graph in 1:length(lab)
    
    cat("END_4\n")
    
    
    legend("topright", legend=lab, fill=mycols, border=mycols, bty="n")
    axis(1, at=seq(0,100))
    axis(2, las=1)
    
    dev.off()
    
    
    
    list_ABSOLUTE_DEF_Freq[[i]]<-Freq_table_sel
    
    
    Condition_DEBUG <- 0
    
    # ############################################################### HERE HERE
    # quit(status = 1)
    
  }################################ i 1:length(CT_levels)
  
  
  
  
  
  
  
  
  
  Condition_DEBUG <- 1
  
  if(length(list_ABSOLUTE_DEF) >0)
  {
    FINAL_df_ABSOLUTE = unique(as.data.frame(data.table::rbindlist(list_ABSOLUTE_DEF, fill=T), stringsAsFactors=F))
    
    
    if(Condition_DEBUG == 1)
    {
      cat("FINAL_df_ABSOLUTE_0\n")
      cat(str(FINAL_df_ABSOLUTE))
      cat("\n")
      #quit(status = 1)
    }
    
    path5<-paste(out2,'LINEAGE_CUMULATIVE_CURVES','/', sep='')
    
    
    if (file.exists(path5)){
      
      
      
      
    } else {
      dir.create(file.path(path5))
      
    }
    
    setwd(path5)
    
    
    
    write.table(FINAL_df_ABSOLUTE, file="E_Plus_ASE_CLASS_gran_mono_log_Rank_test_STAT.tsv", sep="\t", row.names = F, quote = F)
    
    
  }#length(list_ABSOLUTE_DEF) >0
  
  if(length(list_ABSOLUTE_DEF_Freq) >0)
  {
    FINAL_ABSOLUTE_DEF_Freq = unique(as.data.frame(data.table::rbindlist(list_ABSOLUTE_DEF_Freq, fill=T), stringsAsFactors=F))
    
    
    if(Condition_DEBUG == 1)
    {
      cat("FINAL_ABSOLUTE_DEF_Freq_0\n")
      cat(str(FINAL_ABSOLUTE_DEF_Freq))
      cat("\n")
      #quit(status = 1)
    }
    
    path5<-paste(out2,'LINEAGE_CUMULATIVE_CURVES','/', sep='')
    
    if (file.exists(path5)){
      
      
      
      
    } else {
      dir.create(file.path(path5))
      
    }
    
    setwd(path5)
    
    write.table(FINAL_ABSOLUTE_DEF_Freq, file="E_Plus_ASE_CLASS_gran_mono_Cummulative_Freq_table.tsv", sep="\t", quote=F,row.names = F)
    
    
  }#length(list_ABSOLUTE_DEF_Freq) >0
  
  
  
}

cummulative_CT_CLASS_lymph_and_Log_Rank_test_E_Plus_ASE = function(option_list)
{
  suppressMessages(library("nph", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ and transform out2 ----
  
  out2 = opt$out2
  
  cat("out2_\n")
  cat(sprintf(as.character(out2)))
  cat("\n")
  
  setwd(out)
  
  
  filename<-paste('df_Cell_colors','.rds',sep='')
  
  df_Cell_colors<-readRDS(file=filename)
  
  cat("df_Cell_colors_0\n")
  cat(str(df_Cell_colors))
  cat("\n")
  
  
  ### CUMMULATIVE_CLASSES #----
  
  path5<-paste(out2,'LINEAGE_CUMULATIVE_CURVES','/', sep='')
  
  cat("path5\n")
  cat(sprintf(as.character(path5)))
  cat("\n")
  
  setwd(path5)
  
  filename_1<-paste("cummulative_classification_ONLY_ASSAYED_with_Lineage_classification",".rds", sep='')
  
  
  CUMMULATIVE_CLASSES_restricted<-readRDS(file=filename_1)
  
  
  cat("CUMMULATIVE_CLASSES_restricted\n")
  cat(str(CUMMULATIVE_CLASSES_restricted))
  cat("\n")
  cat(str(unique(CUMMULATIVE_CLASSES_restricted$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(CUMMULATIVE_CLASSES_restricted$CLASS_lymph)))))
  cat("\n")
  cat(sprintf(as.character(summary(CUMMULATIVE_CLASSES_restricted$CLASS_lymph))))
  cat("\n")
  
  Condition_DEBUG <- 1
  
  ########### Build df ---------------
  
  Cell_Type_levels<-levels(CUMMULATIVE_CLASSES_restricted$Cell_Type)
  n_Cell_Type_levels<-length(levels(CUMMULATIVE_CLASSES_restricted$Cell_Type))
  n_MPRA_TILES<-5
  VAR_vector<-as.character(unique(CUMMULATIVE_CLASSES_restricted$VAR))
  KEY_Plus_carried_variants_vector<-as.character(unique(CUMMULATIVE_CLASSES_restricted$KEY_Plus_carried_variants))
  
  
  
  n_VAR<-length(VAR_vector)
  
  if(Condition_DEBUG == 1)
  {
    cat("Parameters\n")
    cat(str(Cell_Type_levels))
    cat("\n")
    cat(str(n_Cell_Type_levels))
    cat("\n")
    cat(str(n_MPRA_TILES))
    cat("\n")
    cat(str(VAR_vector))
    cat("\n")
    cat(str(KEY_Plus_carried_variants_vector))
    cat("\n")
    
    
  }
  
  
  #### Build LogRank matrix ----
  
  
  
  Log_rank_df<-as.data.frame(cbind(rep(VAR_vector,n_MPRA_TILES),
                                   rep(KEY_Plus_carried_variants_vector,n_MPRA_TILES),
                                   c(rep(1,n_VAR),rep(2,n_VAR),rep(3,n_VAR),rep(4,n_VAR),rep(5,n_VAR)),
                                   rep("NA",n_MPRA_TILES*n_VAR)),
                             stringsAsFactors=F)
  
  colnames(Log_rank_df)<-c("VAR","KEY_Plus_carried_variants","TILE","ACTIVE")
  
  
  if(Condition_DEBUG == 1)
  {
    cat("Log_rank_df_0\n")
    cat(str(Log_rank_df))
    cat("\n")
    cat(str(unique(Log_rank_df$VAR)))
    cat("\n")
    cat(str(unique(Log_rank_df$KEY_Plus_carried_variants)))
    cat("\n")
  }
  
  
  
  CUMMULATIVE_CLASSES_restricted_subset<-unique(CUMMULATIVE_CLASSES_restricted[,c(which(colnames(CUMMULATIVE_CLASSES_restricted) == "VAR"),
                                                                                  which(colnames(CUMMULATIVE_CLASSES_restricted) == "KEY_Plus_carried_variants"),
                                                                                  which(colnames(CUMMULATIVE_CLASSES_restricted) == "CLASS_lymph"),
                                                                                  which(colnames(CUMMULATIVE_CLASSES_restricted) == "Cell_Type"),
                                                                                  which(colnames(CUMMULATIVE_CLASSES_restricted) == "E_Plus_ASE_CLASS_TILES"))])
  levels_CLASS_lymph<-levels(CUMMULATIVE_CLASSES_restricted$CLASS_lymph)
  
  if(Condition_DEBUG == 1)
  {
    cat("CUMMULATIVE_CLASSES_restricted_subset_0\n")
    cat(str(CUMMULATIVE_CLASSES_restricted_subset))
    cat("\n")
    cat(str(unique(CUMMULATIVE_CLASSES_restricted_subset$VAR)))
    cat("\n")
    cat(str(unique(CUMMULATIVE_CLASSES_restricted_subset$KEY_Plus_carried_variants)))
    cat("\n")
  }
  
  Log_rank_df<-merge(Log_rank_df,
                     CUMMULATIVE_CLASSES_restricted_subset,
                     by=c("VAR","KEY_Plus_carried_variants"),
                     all=T)
  
  Log_rank_df$CLASS_lymph<-Log_rank_df$CLASS_lymph
  
  Log_rank_df$CLASS_lymph<-factor(Log_rank_df$CLASS_lymph,
                                      levels=levels_CLASS_lymph,
                                      ordered=T)
  
  Log_rank_df$Cell_Type<-factor(Log_rank_df$Cell_Type,
                                levels=Cell_Type_levels,
                                ordered=T)
  
  Log_rank_df$TILE<-as.integer(Log_rank_df$TILE)
  
  if(Condition_DEBUG == 1)
  {
    cat("Log_rank_df_1\n")
    cat(str(Log_rank_df))
    cat("\n")
    cat(str(unique(Log_rank_df$VAR)))
    cat("\n")
    cat(str(unique(Log_rank_df$KEY_Plus_carried_variants)))
    cat("\n")
  }
  
  
  
  
  
  ################# Classification of TILES ----------------
  
  Log_rank_df$ACTIVE[which(Log_rank_df$TILE <= Log_rank_df$E_Plus_ASE_CLASS_TILES)]<-"1"
  Log_rank_df$ACTIVE[which(Log_rank_df$TILE > Log_rank_df$E_Plus_ASE_CLASS_TILES)]<-"0"
  
  Log_rank_df$ACTIVE<-as.numeric(Log_rank_df$ACTIVE)
  
  
  
  if(Condition_DEBUG == 1)
  {
    cat("Log_rank_df_2\n")
    cat(str(Log_rank_df))
    cat("\n")
    cat(str(unique(Log_rank_df$VAR)))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(Log_rank_df$ACTIVE))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(Log_rank_df$ACTIVE)))))
    cat("\n")
  }
  
  Log_rank_df<-Log_rank_df[order(Log_rank_df$KEY_Plus_carried_variants,Log_rank_df$Cell_Type, Log_rank_df$TILE),]
  
  
  setwd(out2)
  
  write.table(Log_rank_df, file="check2.tsv",sep="\t",quote=F,row.names = F)
  
  # ################################################################
  # quit(status = 1)
  
  
  ############## Freq Table ----------------------------
  
  ### TOTAL ---
  
  CUMMULATIVE_CLASSES_restricted.dt<-data.table(CUMMULATIVE_CLASSES_restricted, key=c("Cell_Type","CLASS_lymph"))
  
  Freq_TOTAL<-as.data.frame(CUMMULATIVE_CLASSES_restricted.dt[,.(TOTAL=.N),by=key(CUMMULATIVE_CLASSES_restricted.dt)], stringsAsFactors=F)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_TOTAL_0\n")
    cat(str(Freq_TOTAL))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_TOTAL.dt<-data.table(Freq_TOTAL, key=c("CLASS_lymph"))
  
  
  MAXED_df<-as.data.frame(Freq_TOTAL.dt[,.(MAX=max(TOTAL)),by=key(Freq_TOTAL.dt)], stringsAsFactors=F)
  
  
  
  if(Condition_DEBUG == 1)
  {
    cat("MAXED_df_0\n")
    cat(str(MAXED_df))
    cat("\n")
    #quit(status = 1)
  }
  
  ########## Active _tiles --------------------
  
  Log_rank_df.dt<-data.table(Log_rank_df, key=c("Cell_Type","CLASS_lymph","TILE"))
  
  Freq_table<-as.data.frame(Log_rank_df.dt[,.(Freq=sum(ACTIVE)),by=key(Log_rank_df.dt)], stringsAsFactors=F)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_table_0\n")
    cat(str(Freq_table))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_table<-merge(MAXED_df,
                    Freq_table,
                    by=c("CLASS_lymph"),
                    all=T)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_table_1\n")
    cat(str(Freq_table))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_table$Perc<-round(100*(Freq_table$Freq/Freq_table$MAX),2)
  
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_table_2\n")
    cat(str(Freq_table))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_table$CLASS_lymph<-factor(Freq_table$CLASS_lymph,
                                     levels=levels_CLASS_lymph,
                                     ordered=T)
  
  Freq_table$Cell_Type<-factor(Freq_table$Cell_Type,
                               levels=Cell_Type_levels,
                               ordered=T)
  
  Freq_table$TILE<-as.integer(Freq_table$TILE)
  
  Freq_table<-Freq_table[order(Freq_table$CLASS_lymph,Freq_table$Cell_Type, Freq_table$TILE),]
  
  
  setwd(out2)
  
  write.table(Freq_table, file="check3.tsv", row.names = F,quote = F,sep="\t")
  
  # ####################################################################
  # quit(status = 1)
  
  ############### Per CT LOOP ----
  
  
  CT_levels<-Cell_Type_levels
  
  list_ABSOLUTE_DEF<-list()
  list_ABSOLUTE_DEF_Freq<-list()
  
  Condition_DEBUG <- 1
  
  for(i in 1:length(CT_levels))
  {
    CT_sel<-CT_levels[i]
    
    cat("CT:---------------------------------------->\t")
    cat(sprintf(as.character(CT_sel)))
    cat("\n")
    
    Log_rank_df_sel<-Log_rank_df[which(Log_rank_df$Cell_Type == CT_sel),]
    
    if(Condition_DEBUG == 1)
    {
      cat("Log_rank_df_sel_0\n")
      cat(str(Log_rank_df_sel))
      cat("\n")
      cat(str(unique(Log_rank_df_sel$VAR)))
      cat("\n")
    }
    
    
    #### Log Rank test ----
    
    Condition_DEBUG <- 0
    
    list_DEF<-list()
    
    CLASS_lymph_labels<-levels(Log_rank_df_sel$CLASS_lymph)
    
    for(k in 1 :length(CLASS_lymph_labels))
    {
      BAIT<-CLASS_lymph_labels[k]
      
      cat("-----BAIT------>\t")
      cat(sprintf(as.character(BAIT)))
      cat("\t")
      
      list_PREY<-list()
      
      for(h in 1 :length(CLASS_lymph_labels))
      {
        
        PREY<-CLASS_lymph_labels[h]
        
        cat("------PREY----->\t")
        cat(sprintf(as.character(PREY)))
        cat("\n")
        
        Mini_comparison<-unique(c(BAIT,PREY))
        
        if(Condition_DEBUG == 1)
        {
          cat("Mini_comparison_0\n")
          cat(str(Mini_comparison))
          cat("\n")
        }
        
        Log_rank_df_sel_subset<-Log_rank_df_sel[which(Log_rank_df_sel$CLASS_lymph%in%Mini_comparison),]
        
        if(Condition_DEBUG == 1)
        {
          cat("Log_rank_df_sel_subset_0\n")
          cat(str(Log_rank_df_sel_subset))
          cat("\n")
        }
        
        FLAG_diversity<-as.integer(length(unique(Log_rank_df_sel_subset$CLASS_lymph)))
        
        if(Condition_DEBUG == 1)
        {
          cat("FLAG_diversity_0\n")
          cat(sprintf(as.character(FLAG_diversity)))
          cat("\n")
          
          cat(str(FLAG_diversity))
          cat("\n")
        }
        
        if(FLAG_diversity >1)
        {
          cat("Hello_world\n")
          
          Mini_comparison_c1<-as.character(Mini_comparison[1])
          Mini_comparison_c2<-as.character(Mini_comparison[2])
          
          FLAG_diversity_2<-length(unique(Log_rank_df_sel_subset$ACTIVE))
          
          if(Condition_DEBUG == 1)
          {
            cat("FLAG_diversity_2_0\n")
            cat(str(FLAG_diversity_2))
            cat("\n")
          }
          
          if(FLAG_diversity_2 >1)
          {
            if(Condition_DEBUG == 1)
            {
              cat("----------->\t")
              cat(sprintf(as.character(Mini_comparison_c1)))
              cat("\t")
              cat(sprintf(as.character(Mini_comparison_c2)))
              cat("\n")
              
              cat("Log_rank_df_sel_subset_0\n")
              cat(str(Log_rank_df_sel_subset))
              cat("\n")
              cat(str(unique(Log_rank_df_sel_subset$VAR)))
              cat("\n")
              cat(sprintf(as.character(names(summary(as.factor(Log_rank_df_sel_subset$ACTIVE))))))
              cat("\n")
              cat(sprintf(as.character(summary(as.factor(Log_rank_df_sel_subset$ACTIVE)))))
              cat("\n")
            }
            
            #, "less", "greater"
            
            STATS_test<-logrank.test(
              Log_rank_df_sel_subset$TILE,
              Log_rank_df_sel_subset$ACTIVE,
              Log_rank_df_sel_subset$CLASS_lymph,
              alternative = c("two.sided"),
              rho = 0,
              gamma = 0,
              event_time_weights = NULL
            )
            
            # if(Condition_DEBUG == 1)
            # {
            #   cat("STATS_test_0\n")
            #   cat(str(STATS_test))
            #   cat("\n")
            # }
            
            
            p_value<-STATS_test$test$p
            minus_log_val<-round(-1*log10(p_value),2)
            
            if(Condition_DEBUG == 1)
            {
              cat("minus_log_val_0\n")
              cat(str(minus_log_val))
              cat("\n")
            }
            
            
            a.df<-as.data.frame(cbind(Mini_comparison_c1,Mini_comparison_c2,minus_log_val), stringsAsFactors=F)
            
            
            colnames(a.df)<-c("CLASS_lymph_c1","CLASS_lymph_c2","minus_logpval")
            a.df$minus_logpval<-as.numeric(a.df$minus_logpval)
            
            if(Condition_DEBUG == 1)
            {
              cat("a.df\n")
              cat(str(a.df))
              cat("\n")
            }
            
            list_PREY[[h]]<-a.df
            
            # ###############################################
            # quit(status = 1)
            
          }#FLAG_diversity_2 >1
        }#FLAG_diversity >1
      }# h in 1 :length(unique(Log_rank_df_sel$CLASS_lymph)
      
      Condition_DEBUG <- 0
      
      if(length(list_PREY) >0)
      {
        
        PREY_df = unique(as.data.frame(data.table::rbindlist(list_PREY, fill=T), stringsAsFactors=F))
        
        if(Condition_DEBUG == 1)
        {
          cat("PREY_df_0\n")
          cat(str(PREY_df))
          cat("\n")
          #quit(status = 1)
        }
        
        list_DEF[[k]]<-PREY_df
        
      }# length(list_PREY) >0
      
      Condition_DEBUG <- 0
      
      
    }#k in 1 :length(unique(Log_rank_df_sel$CLASS_lymph))
    
    
    Condition_DEBUG <- 0
    
    if(length(list_DEF) >0)
    {
      FINAL_df_per_CT = unique(as.data.frame(data.table::rbindlist(list_DEF, fill=T), stringsAsFactors=F))
      
      FINAL_df_per_CT$Cell_Type<-CT_sel
      
      if(Condition_DEBUG == 1)
      {
        cat("FINAL_df_per_CT_0\n")
        cat(str(FINAL_df_per_CT))
        cat("\n")
        #quit(status = 1)
      }
      
      list_ABSOLUTE_DEF[[i]]<-FINAL_df_per_CT
      
    }#length(list_DEF) >0
    
    Condition_DEBUG <- 0
    
    #### GRAPH -----
    
    Condition_DEBUG <- 1
    
    vector_colors_CLASS_lymph<-c("gray","dodgerblue")
    
    Freq_table_sel<-Freq_table[which(Freq_table$Cell_Type == CT_sel),]
    
    if(Condition_DEBUG == 1)
    {
      cat("Freq_table_sel_0\n")
      cat(str(Freq_table_sel))
      cat("\n")
      #quit(status = 1)
    }
    
    
    
    
    path5<-paste(out2,'LINEAGE_CUMULATIVE_CURVES','/', sep='')
    
    if (file.exists(path5)){
      
      
      
      
    } else {
      dir.create(file.path(path5))
      
    }
    
    setwd(path5)
    
    v_parameter<-"NA" 
    
    indx_CLASS_lymph<-which(colnames(Freq_table_sel) == "CLASS_lymph")
    indx.yaxis<-which(colnames(Freq_table_sel) == "Perc")
    mycols <- vector_colors_CLASS_lymph
    
    pdfname<-paste(CT_sel,"_CLASS_lymph","_Cummulative_frequency_","E_Plus_ASE_activity",".pdf",sep='')
    pdf(file=pdfname, width=5, height=4, pointsize=12)
    
    par(mai=c(0.9,0.9,0.3,0.2))
    lab <- as.character(unique(Freq_table_sel[,indx_CLASS_lymph]))
    
    cat("lab\n")
    cat(sprintf(as.character(lab)))
    cat("\n")
    
    
    plot(Freq_table_sel$TILE, Freq_table_sel[,indx.yaxis],
         ty="n", xlab="TILES with E_Plus_ASE activity",
         ylab="Cummulative % of variants",
         axes=F, cex.lab=1.2, cex.lab=1.3, ylim=c(0,100))
    
    
    
    # cat("Hello_world1\n")
    
    # points(Freq_table_sel$TILE, Freq_table_sel[,indx.yaxis], col="darkgrey", pch=19)
    
    # cat("Hello_world2\n")
    # 
    # 
    for (iteration_graph in 1:length(lab))  {
      
      cat("Hello_world3\t")
      cat(sprintf(as.character(lab[iteration_graph])))
      cat("\n")
      
      ind <- which(Freq_table_sel[,indx_CLASS_lymph]==lab[iteration_graph])
      
      cat("------------------->ind\n")
      cat(str(ind))
      cat("\n")
      color_sel<-mycols[iteration_graph]
      
      cat(sprintf(as.character(color_sel)))
      cat("\n")
      
      points(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind],
             pch=21,lwd=5, col=color_sel)
      lines(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind],
            lty=1,lwd=5, col=color_sel)
      
      cat("END_3\n")
      
      # points(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind], pch=19, col=color_sel)
    }#iteration_graph in 1:length(lab)
    
    cat("END_4\n")
    
    
    legend("topright", legend=lab, fill=mycols, border=mycols, bty="n")
    axis(1, at=seq(0,100))
    axis(2, las=1)
    
    dev.off()
    
    
    
    list_ABSOLUTE_DEF_Freq[[i]]<-Freq_table_sel
    
    
    Condition_DEBUG <- 0
    
    # ############################################################### HERE HERE
    # quit(status = 1)
    
  }################################ i 1:length(CT_levels)
  
  
  
  
  
  
  
  
  
  Condition_DEBUG <- 1
  
  if(length(list_ABSOLUTE_DEF) >0)
  {
    FINAL_df_ABSOLUTE = unique(as.data.frame(data.table::rbindlist(list_ABSOLUTE_DEF, fill=T), stringsAsFactors=F))
    
    
    if(Condition_DEBUG == 1)
    {
      cat("FINAL_df_ABSOLUTE_0\n")
      cat(str(FINAL_df_ABSOLUTE))
      cat("\n")
      #quit(status = 1)
    }
    
    path5<-paste(out2,'LINEAGE_CUMULATIVE_CURVES','/', sep='')
    
    
    if (file.exists(path5)){
      
      
      
      
    } else {
      dir.create(file.path(path5))
      
    }
    
    setwd(path5)
    
    
    
    write.table(FINAL_df_ABSOLUTE, file="E_Plus_ASE_CLASS_lymph_log_Rank_test_STAT.tsv", sep="\t", row.names = F, quote = F)
    
    
  }#length(list_ABSOLUTE_DEF) >0
  
  if(length(list_ABSOLUTE_DEF_Freq) >0)
  {
    FINAL_ABSOLUTE_DEF_Freq = unique(as.data.frame(data.table::rbindlist(list_ABSOLUTE_DEF_Freq, fill=T), stringsAsFactors=F))
    
    
    if(Condition_DEBUG == 1)
    {
      cat("FINAL_ABSOLUTE_DEF_Freq_0\n")
      cat(str(FINAL_ABSOLUTE_DEF_Freq))
      cat("\n")
      #quit(status = 1)
    }
    
    path5<-paste(out2,'LINEAGE_CUMULATIVE_CURVES','/', sep='')
    
    if (file.exists(path5)){
      
      
      
      
    } else {
      dir.create(file.path(path5))
      
    }
    
    setwd(path5)
    
    write.table(FINAL_ABSOLUTE_DEF_Freq, file="E_Plus_ASE_CLASS_lymph_Cummulative_Freq_table.tsv", sep="\t", quote=F,row.names = F)
    
    
  }#length(list_ABSOLUTE_DEF_Freq) >0
  
  
  
}

cummulative_CT_CLASS_ALL_wbc_and_Log_Rank_test_E_Plus_ASE = function(option_list)
{
  suppressMessages(library("nph", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ and transform out2 ----
  
  out2 = opt$out2
  
  cat("out2_\n")
  cat(sprintf(as.character(out2)))
  cat("\n")
  
  setwd(out)
  
  
  filename<-paste('df_Cell_colors','.rds',sep='')
  
  df_Cell_colors<-readRDS(file=filename)
  
  cat("df_Cell_colors_0\n")
  cat(str(df_Cell_colors))
  cat("\n")
  
  
  ### CUMMULATIVE_CLASSES #----
  
  path5<-paste(out2,'LINEAGE_CUMULATIVE_CURVES','/', sep='')
  
  cat("path5\n")
  cat(sprintf(as.character(path5)))
  cat("\n")
  
  setwd(path5)
  
  filename_1<-paste("cummulative_classification_ONLY_ASSAYED_with_Lineage_classification",".rds", sep='')
  
  
  CUMMULATIVE_CLASSES_restricted<-readRDS(file=filename_1)
  
  
  cat("CUMMULATIVE_CLASSES_restricted\n")
  cat(str(CUMMULATIVE_CLASSES_restricted))
  cat("\n")
  cat(str(unique(CUMMULATIVE_CLASSES_restricted$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(CUMMULATIVE_CLASSES_restricted$CLASS_ALL_wbc)))))
  cat("\n")
  cat(sprintf(as.character(summary(CUMMULATIVE_CLASSES_restricted$CLASS_ALL_wbc))))
  cat("\n")
  
  Condition_DEBUG <- 1
  
  ########### Build df ---------------
  
  Cell_Type_levels<-levels(CUMMULATIVE_CLASSES_restricted$Cell_Type)
  n_Cell_Type_levels<-length(levels(CUMMULATIVE_CLASSES_restricted$Cell_Type))
  n_MPRA_TILES<-5
  VAR_vector<-as.character(unique(CUMMULATIVE_CLASSES_restricted$VAR))
  KEY_Plus_carried_variants_vector<-as.character(unique(CUMMULATIVE_CLASSES_restricted$KEY_Plus_carried_variants))
  
  
  
  n_VAR<-length(VAR_vector)
  
  if(Condition_DEBUG == 1)
  {
    cat("Parameters\n")
    cat(str(Cell_Type_levels))
    cat("\n")
    cat(str(n_Cell_Type_levels))
    cat("\n")
    cat(str(n_MPRA_TILES))
    cat("\n")
    cat(str(VAR_vector))
    cat("\n")
    cat(str(KEY_Plus_carried_variants_vector))
    cat("\n")
    
    
  }
  
  
  #### Build LogRank matrix ----
  
  
  
  Log_rank_df<-as.data.frame(cbind(rep(VAR_vector,n_MPRA_TILES),
                                   rep(KEY_Plus_carried_variants_vector,n_MPRA_TILES),
                                   c(rep(1,n_VAR),rep(2,n_VAR),rep(3,n_VAR),rep(4,n_VAR),rep(5,n_VAR)),
                                   rep("NA",n_MPRA_TILES*n_VAR)),
                             stringsAsFactors=F)
  
  colnames(Log_rank_df)<-c("VAR","KEY_Plus_carried_variants","TILE","ACTIVE")
  
  
  if(Condition_DEBUG == 1)
  {
    cat("Log_rank_df_0\n")
    cat(str(Log_rank_df))
    cat("\n")
    cat(str(unique(Log_rank_df$VAR)))
    cat("\n")
    cat(str(unique(Log_rank_df$KEY_Plus_carried_variants)))
    cat("\n")
  }
  
  
  
  CUMMULATIVE_CLASSES_restricted_subset<-unique(CUMMULATIVE_CLASSES_restricted[,c(which(colnames(CUMMULATIVE_CLASSES_restricted) == "VAR"),
                                                                                  which(colnames(CUMMULATIVE_CLASSES_restricted) == "KEY_Plus_carried_variants"),
                                                                                  which(colnames(CUMMULATIVE_CLASSES_restricted) == "CLASS_ALL_wbc"),
                                                                                  which(colnames(CUMMULATIVE_CLASSES_restricted) == "Cell_Type"),
                                                                                  which(colnames(CUMMULATIVE_CLASSES_restricted) == "E_Plus_ASE_CLASS_TILES"))])
  levels_CLASS_ALL_wbc<-levels(CUMMULATIVE_CLASSES_restricted$CLASS_ALL_wbc)
  
  if(Condition_DEBUG == 1)
  {
    cat("CUMMULATIVE_CLASSES_restricted_subset_0\n")
    cat(str(CUMMULATIVE_CLASSES_restricted_subset))
    cat("\n")
    cat(str(unique(CUMMULATIVE_CLASSES_restricted_subset$VAR)))
    cat("\n")
    cat(str(unique(CUMMULATIVE_CLASSES_restricted_subset$KEY_Plus_carried_variants)))
    cat("\n")
  }
  
  Log_rank_df<-merge(Log_rank_df,
                     CUMMULATIVE_CLASSES_restricted_subset,
                     by=c("VAR","KEY_Plus_carried_variants"),
                     all=T)
  
  Log_rank_df$CLASS_ALL_wbc<-Log_rank_df$CLASS_ALL_wbc
  
  Log_rank_df$CLASS_ALL_wbc<-factor(Log_rank_df$CLASS_ALL_wbc,
                                  levels=levels_CLASS_ALL_wbc,
                                  ordered=T)
  
  Log_rank_df$Cell_Type<-factor(Log_rank_df$Cell_Type,
                                levels=Cell_Type_levels,
                                ordered=T)
  
  Log_rank_df$TILE<-as.integer(Log_rank_df$TILE)
  
  if(Condition_DEBUG == 1)
  {
    cat("Log_rank_df_1\n")
    cat(str(Log_rank_df))
    cat("\n")
    cat(str(unique(Log_rank_df$VAR)))
    cat("\n")
    cat(str(unique(Log_rank_df$KEY_Plus_carried_variants)))
    cat("\n")
  }
  
  
  
  
  
  ################# Classification of TILES ----------------
  
  Log_rank_df$ACTIVE[which(Log_rank_df$TILE <= Log_rank_df$E_Plus_ASE_CLASS_TILES)]<-"1"
  Log_rank_df$ACTIVE[which(Log_rank_df$TILE > Log_rank_df$E_Plus_ASE_CLASS_TILES)]<-"0"
  
  Log_rank_df$ACTIVE<-as.numeric(Log_rank_df$ACTIVE)
  
  
  
  if(Condition_DEBUG == 1)
  {
    cat("Log_rank_df_2\n")
    cat(str(Log_rank_df))
    cat("\n")
    cat(str(unique(Log_rank_df$VAR)))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(Log_rank_df$ACTIVE))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(Log_rank_df$ACTIVE)))))
    cat("\n")
  }
  
  Log_rank_df<-Log_rank_df[order(Log_rank_df$KEY_Plus_carried_variants,Log_rank_df$Cell_Type, Log_rank_df$TILE),]
  
  
  setwd(out2)
  
  write.table(Log_rank_df, file="check2.tsv",sep="\t",quote=F,row.names = F)
  
  # ################################################################
  # quit(status = 1)
  
  
  ############## Freq Table ----------------------------
  
  ### TOTAL ---
  
  CUMMULATIVE_CLASSES_restricted.dt<-data.table(CUMMULATIVE_CLASSES_restricted, key=c("Cell_Type","CLASS_ALL_wbc"))
  
  Freq_TOTAL<-as.data.frame(CUMMULATIVE_CLASSES_restricted.dt[,.(TOTAL=.N),by=key(CUMMULATIVE_CLASSES_restricted.dt)], stringsAsFactors=F)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_TOTAL_0\n")
    cat(str(Freq_TOTAL))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_TOTAL.dt<-data.table(Freq_TOTAL, key=c("CLASS_ALL_wbc"))
  
  
  MAXED_df<-as.data.frame(Freq_TOTAL.dt[,.(MAX=max(TOTAL)),by=key(Freq_TOTAL.dt)], stringsAsFactors=F)
  
  
  
  if(Condition_DEBUG == 1)
  {
    cat("MAXED_df_0\n")
    cat(str(MAXED_df))
    cat("\n")
    #quit(status = 1)
  }
  
  ########## Active _tiles --------------------
  
  Log_rank_df.dt<-data.table(Log_rank_df, key=c("Cell_Type","CLASS_ALL_wbc","TILE"))
  
  Freq_table<-as.data.frame(Log_rank_df.dt[,.(Freq=sum(ACTIVE)),by=key(Log_rank_df.dt)], stringsAsFactors=F)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_table_0\n")
    cat(str(Freq_table))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_table<-merge(MAXED_df,
                    Freq_table,
                    by=c("CLASS_ALL_wbc"),
                    all=T)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_table_1\n")
    cat(str(Freq_table))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_table$Perc<-round(100*(Freq_table$Freq/Freq_table$MAX),2)
  
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_table_2\n")
    cat(str(Freq_table))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_table$CLASS_ALL_wbc<-factor(Freq_table$CLASS_ALL_wbc,
                                 levels=levels_CLASS_ALL_wbc,
                                 ordered=T)
  
  Freq_table$Cell_Type<-factor(Freq_table$Cell_Type,
                               levels=Cell_Type_levels,
                               ordered=T)
  
  Freq_table$TILE<-as.integer(Freq_table$TILE)
  
  Freq_table<-Freq_table[order(Freq_table$CLASS_ALL_wbc,Freq_table$Cell_Type, Freq_table$TILE),]
  
  
  setwd(out2)
  
  write.table(Freq_table, file="check3.tsv", row.names = F,quote = F,sep="\t")
  
  # ####################################################################
  # quit(status = 1)
  
  ############### Per CT LOOP ----
  
  
  CT_levels<-Cell_Type_levels
  
  list_ABSOLUTE_DEF<-list()
  list_ABSOLUTE_DEF_Freq<-list()
  
  Condition_DEBUG <- 1
  
  for(i in 1:length(CT_levels))
  {
    CT_sel<-CT_levels[i]
    
    cat("CT:---------------------------------------->\t")
    cat(sprintf(as.character(CT_sel)))
    cat("\n")
    
    Log_rank_df_sel<-Log_rank_df[which(Log_rank_df$Cell_Type == CT_sel),]
    
    if(Condition_DEBUG == 1)
    {
      cat("Log_rank_df_sel_0\n")
      cat(str(Log_rank_df_sel))
      cat("\n")
      cat(str(unique(Log_rank_df_sel$VAR)))
      cat("\n")
    }
    
    
    #### Log Rank test ----
    
    Condition_DEBUG <- 0
    
    list_DEF<-list()
    
    CLASS_ALL_wbc_labels<-levels(Log_rank_df_sel$CLASS_ALL_wbc)
    
    for(k in 1 :length(CLASS_ALL_wbc_labels))
    {
      BAIT<-CLASS_ALL_wbc_labels[k]
      
      cat("-----BAIT------>\t")
      cat(sprintf(as.character(BAIT)))
      cat("\t")
      
      list_PREY<-list()
      
      for(h in 1 :length(CLASS_ALL_wbc_labels))
      {
        
        PREY<-CLASS_ALL_wbc_labels[h]
        
        cat("------PREY----->\t")
        cat(sprintf(as.character(PREY)))
        cat("\n")
        
        Mini_comparison<-unique(c(BAIT,PREY))
        
        if(Condition_DEBUG == 1)
        {
          cat("Mini_comparison_0\n")
          cat(str(Mini_comparison))
          cat("\n")
        }
        
        Log_rank_df_sel_subset<-Log_rank_df_sel[which(Log_rank_df_sel$CLASS_ALL_wbc%in%Mini_comparison),]
        
        if(Condition_DEBUG == 1)
        {
          cat("Log_rank_df_sel_subset_0\n")
          cat(str(Log_rank_df_sel_subset))
          cat("\n")
        }
        
        FLAG_diversity<-as.integer(length(unique(Log_rank_df_sel_subset$CLASS_ALL_wbc)))
        
        if(Condition_DEBUG == 1)
        {
          cat("FLAG_diversity_0\n")
          cat(sprintf(as.character(FLAG_diversity)))
          cat("\n")
          
          cat(str(FLAG_diversity))
          cat("\n")
        }
        
        if(FLAG_diversity >1)
        {
          cat("Hello_world\n")
          
          Mini_comparison_c1<-as.character(Mini_comparison[1])
          Mini_comparison_c2<-as.character(Mini_comparison[2])
          
          FLAG_diversity_2<-length(unique(Log_rank_df_sel_subset$ACTIVE))
          
          if(Condition_DEBUG == 1)
          {
            cat("FLAG_diversity_2_0\n")
            cat(str(FLAG_diversity_2))
            cat("\n")
          }
          
          if(FLAG_diversity_2 >1)
          {
            if(Condition_DEBUG == 1)
            {
              cat("----------->\t")
              cat(sprintf(as.character(Mini_comparison_c1)))
              cat("\t")
              cat(sprintf(as.character(Mini_comparison_c2)))
              cat("\n")
              
              cat("Log_rank_df_sel_subset_0\n")
              cat(str(Log_rank_df_sel_subset))
              cat("\n")
              cat(str(unique(Log_rank_df_sel_subset$VAR)))
              cat("\n")
              cat(sprintf(as.character(names(summary(as.factor(Log_rank_df_sel_subset$ACTIVE))))))
              cat("\n")
              cat(sprintf(as.character(summary(as.factor(Log_rank_df_sel_subset$ACTIVE)))))
              cat("\n")
            }
            
            #, "less", "greater"
            
            STATS_test<-logrank.test(
              Log_rank_df_sel_subset$TILE,
              Log_rank_df_sel_subset$ACTIVE,
              Log_rank_df_sel_subset$CLASS_ALL_wbc,
              alternative = c("two.sided"),
              rho = 0,
              gamma = 0,
              event_time_weights = NULL
            )
            
            # if(Condition_DEBUG == 1)
            # {
            #   cat("STATS_test_0\n")
            #   cat(str(STATS_test))
            #   cat("\n")
            # }
            
            
            p_value<-STATS_test$test$p
            minus_log_val<-round(-1*log10(p_value),2)
            
            if(Condition_DEBUG == 1)
            {
              cat("minus_log_val_0\n")
              cat(str(minus_log_val))
              cat("\n")
            }
            
            
            a.df<-as.data.frame(cbind(Mini_comparison_c1,Mini_comparison_c2,minus_log_val), stringsAsFactors=F)
            
            
            colnames(a.df)<-c("CLASS_ALL_wbc_c1","CLASS_ALL_wbc_c2","minus_logpval")
            a.df$minus_logpval<-as.numeric(a.df$minus_logpval)
            
            if(Condition_DEBUG == 1)
            {
              cat("a.df\n")
              cat(str(a.df))
              cat("\n")
            }
            
            list_PREY[[h]]<-a.df
            
            # ###############################################
            # quit(status = 1)
            
          }#FLAG_diversity_2 >1
        }#FLAG_diversity >1
      }# h in 1 :length(unique(Log_rank_df_sel$CLASS_ALL_wbc)
      
      Condition_DEBUG <- 0
      
      if(length(list_PREY) >0)
      {
        
        PREY_df = unique(as.data.frame(data.table::rbindlist(list_PREY, fill=T), stringsAsFactors=F))
        
        if(Condition_DEBUG == 1)
        {
          cat("PREY_df_0\n")
          cat(str(PREY_df))
          cat("\n")
          #quit(status = 1)
        }
        
        list_DEF[[k]]<-PREY_df
        
      }# length(list_PREY) >0
      
      Condition_DEBUG <- 0
      
      
    }#k in 1 :length(unique(Log_rank_df_sel$CLASS_ALL_wbc))
    
    
    Condition_DEBUG <- 0
    
    if(length(list_DEF) >0)
    {
      FINAL_df_per_CT = unique(as.data.frame(data.table::rbindlist(list_DEF, fill=T), stringsAsFactors=F))
      
      FINAL_df_per_CT$Cell_Type<-CT_sel
      
      if(Condition_DEBUG == 1)
      {
        cat("FINAL_df_per_CT_0\n")
        cat(str(FINAL_df_per_CT))
        cat("\n")
        #quit(status = 1)
      }
      
      list_ABSOLUTE_DEF[[i]]<-FINAL_df_per_CT
      
    }#length(list_DEF) >0
    
    Condition_DEBUG <- 0
    
    #### GRAPH -----
    
    Condition_DEBUG <- 1
    
    vector_colors_CLASS_ALL_wbc<-c("gray","yellow2")
    
    Freq_table_sel<-Freq_table[which(Freq_table$Cell_Type == CT_sel),]
    
    if(Condition_DEBUG == 1)
    {
      cat("Freq_table_sel_0\n")
      cat(str(Freq_table_sel))
      cat("\n")
      #quit(status = 1)
    }
    
    
    
    
    path5<-paste(out2,'LINEAGE_CUMULATIVE_CURVES','/', sep='')
    
    if (file.exists(path5)){
      
      
      
      
    } else {
      dir.create(file.path(path5))
      
    }
    
    setwd(path5)
    
    v_parameter<-"NA" 
    
    indx_CLASS_ALL_wbc<-which(colnames(Freq_table_sel) == "CLASS_ALL_wbc")
    indx.yaxis<-which(colnames(Freq_table_sel) == "Perc")
    mycols <- vector_colors_CLASS_ALL_wbc
    
    pdfname<-paste(CT_sel,"_CLASS_ALL_wbc","_Cummulative_frequency_","E_Plus_ASE_activity",".pdf",sep='')
    pdf(file=pdfname, width=5, height=4, pointsize=12)
    
    par(mai=c(0.9,0.9,0.3,0.2))
    lab <- as.character(unique(Freq_table_sel[,indx_CLASS_ALL_wbc]))
    
    cat("lab\n")
    cat(sprintf(as.character(lab)))
    cat("\n")
    
    
    plot(Freq_table_sel$TILE, Freq_table_sel[,indx.yaxis],
         ty="n", xlab="TILES with E_Plus_ASE activity",
         ylab="Cummulative % of variants",
         axes=F, cex.lab=1.2, cex.lab=1.3, ylim=c(0,100))
    
    
    
    # cat("Hello_world1\n")
    
    # points(Freq_table_sel$TILE, Freq_table_sel[,indx.yaxis], col="darkgrey", pch=19)
    
    # cat("Hello_world2\n")
    # 
    # 
    for (iteration_graph in 1:length(lab))  {
      
      cat("Hello_world3\t")
      cat(sprintf(as.character(lab[iteration_graph])))
      cat("\n")
      
      ind <- which(Freq_table_sel[,indx_CLASS_ALL_wbc]==lab[iteration_graph])
      
      cat("------------------->ind\n")
      cat(str(ind))
      cat("\n")
      color_sel<-mycols[iteration_graph]
      
      cat(sprintf(as.character(color_sel)))
      cat("\n")
      
      points(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind],
             pch=21,lwd=5, col=color_sel)
      lines(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind],
            lty=1,lwd=5, col=color_sel)
      
      cat("END_3\n")
      
      # points(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind], pch=19, col=color_sel)
    }#iteration_graph in 1:length(lab)
    
    cat("END_4\n")
    
    
    legend("topright", legend=lab, fill=mycols, border=mycols, bty="n")
    axis(1, at=seq(0,100))
    axis(2, las=1)
    
    dev.off()
    
    
    
    list_ABSOLUTE_DEF_Freq[[i]]<-Freq_table_sel
    
    
    Condition_DEBUG <- 0
    
    # ############################################################### HERE HERE
    # quit(status = 1)
    
  }################################ i 1:length(CT_levels)
  
  
  
  
  
  
  
  
  
  Condition_DEBUG <- 1
  
  if(length(list_ABSOLUTE_DEF) >0)
  {
    FINAL_df_ABSOLUTE = unique(as.data.frame(data.table::rbindlist(list_ABSOLUTE_DEF, fill=T), stringsAsFactors=F))
    
    
    if(Condition_DEBUG == 1)
    {
      cat("FINAL_df_ABSOLUTE_0\n")
      cat(str(FINAL_df_ABSOLUTE))
      cat("\n")
      #quit(status = 1)
    }
    
    path5<-paste(out2,'LINEAGE_CUMULATIVE_CURVES','/', sep='')
    
    
    if (file.exists(path5)){
      
      
      
      
    } else {
      dir.create(file.path(path5))
      
    }
    
    setwd(path5)
    
    
    
    write.table(FINAL_df_ABSOLUTE, file="E_Plus_ASE_CLASS_ALL_wbc_log_Rank_test_STAT.tsv", sep="\t", row.names = F, quote = F)
    
    
  }#length(list_ABSOLUTE_DEF) >0
  
  if(length(list_ABSOLUTE_DEF_Freq) >0)
  {
    FINAL_ABSOLUTE_DEF_Freq = unique(as.data.frame(data.table::rbindlist(list_ABSOLUTE_DEF_Freq, fill=T), stringsAsFactors=F))
    
    
    if(Condition_DEBUG == 1)
    {
      cat("FINAL_ABSOLUTE_DEF_Freq_0\n")
      cat(str(FINAL_ABSOLUTE_DEF_Freq))
      cat("\n")
      #quit(status = 1)
    }
    
    path5<-paste(out2,'LINEAGE_CUMULATIVE_CURVES','/', sep='')
    
    if (file.exists(path5)){
      
      
      
      
    } else {
      dir.create(file.path(path5))
      
    }
    
    setwd(path5)
    
    write.table(FINAL_ABSOLUTE_DEF_Freq, file="E_Plus_ASE_CLASS_ALL_wbc_Cummulative_Freq_table.tsv", sep="\t", quote=F,row.names = F)
    
    
  }#length(list_ABSOLUTE_DEF_Freq) >0
  
  
  
}


printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}


#### main script ----

main = function() {
  
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  option_list <- list(
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out2"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--CUMMULATIVE_CLASSES"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--dB"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--TOME_correspondence"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--finemap_prob_Threshold"), type="numeric", default=NULL, 
                metavar="filename", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "137_MPRA_normalization_and_filtering_Rscript_v2.R
                        --regular_table FILE.txt
                        --replicas charac
                        --type1 type1
                        --type2 type2
                        --pvalThreshold integer
                        --FDRThreshold integer
                        --EquivalenceTable FILE.txt
                        --sharpr2Threshold charac",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  
  Classification_of_variants_by_lineages(opt)
  cummulative_CT_CLASS_mega_and_Log_Rank_test_E_Plus_ASE(opt)
  cummulative_CT_CLASS_erythroid_and_Log_Rank_test_E_Plus_ASE(opt)
  cummulative_CT_CLASS_gran_mono_and_Log_Rank_test_E_Plus_ASE(opt)
  cummulative_CT_CLASS_lymph_and_Log_Rank_test_E_Plus_ASE(opt)
  cummulative_CT_CLASS_ALL_wbc_and_Log_Rank_test_E_Plus_ASE(opt)
  
}
  
  
  
 

###########################################################################

system.time( main() )
