#### WARNING HARDCODED INDEXES -----


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
library("ggeasy",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("reshape2",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("viridisLite",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")


opt = NULL

options(warn=1)

data_wrangling = function(option_list)
{
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
  
  #### READ and transform out ----
  
  finemap_prob_Threshold = opt$finemap_prob_Threshold
  
  cat("finemap_prob_Threshold\n")
  cat(sprintf(as.character(finemap_prob_Threshold)))
  cat("\n")
  
  #### Categories colors ----
  
  
  df_CSQ_colors<-readRDS(file=opt$CSQ_colors)
  
  df_CSQ_colors$color[df_CSQ_colors$VEP_DEF_LABELS == "TFBS"]<-"red"
  
  
  cat("df_CSQ_colors_0\n")
  cat(str(df_CSQ_colors))
  cat("\n")
  
  #### READ DATA----
  
  setwd(out)
  
  filename<-paste('Element_collapse','.rds',sep='')
  
  
  KEY_collapse<-readRDS(file=filename)
  
  cat("KEY_collapse_0\n")
  cat(str(KEY_collapse))
  cat("\n")
  
  
  filename<-paste('df_Cell_colors','.rds',sep='')
  
  df_Cell_colors<-readRDS(file=filename)
  
  cat("df_Cell_colors_0\n")
  cat(str(df_Cell_colors))
  cat("\n")
  
  #### Categories colors ----
  
  
  df_CSQ_colors<-readRDS(file=opt$CSQ_colors)
  
  df_CSQ_colors$color[df_CSQ_colors$VEP_DEF_LABELS == "TFBS"]<-"red"
  
  
  cat("df_CSQ_colors_0\n")
  cat(str(df_CSQ_colors))
  cat("\n")
  
  color_vector<-c("black",'#1877C9','gray')
  
  A.df<-as.data.frame(cbind(c("NCGR","ASE_CTRL","enhancer_CTRL"),color_vector), stringsAsFactors=F)
  
  colnames(A.df)<-c("VEP_DEF_LABELS","color")
  
  cat("A.df_1\n")
  cat(str(A.df))
  cat("\n")
  
  
  
  df_CSQ_colors<-rbind(df_CSQ_colors,A.df)
  
  
  df_CSQ_colors$VEP_DEF_LABELS<-factor(df_CSQ_colors$VEP_DEF_LABELS,
                                       levels=c("LOF","MISS","SYN","UTR5","UTR3",
                                                "INTRON","INTERGENIC","UPSTREAM","DOWNSTREAM","REGULATORY",
                                                "TFBS","SPLICE",
                                                "OTHER","NMD","NCT","PCHIC_Relevant_link",
                                                "NCGR","ASE_CTRL","enhancer_CTRL","Kousik_variant"),
                                       ordered=T)
  
  
  cat("df_CSQ_colors_2\n")
  cat(str(df_CSQ_colors))
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
  

  dB_subset<-dB[which(dB$VAR%in%KEY_collapse$VAR),]
  
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
  
  KEY_collapse<-merge(KEY_collapse,
                    VAR_lineage_CLASSIF,
                    by="VAR",
                    all.x=T)
  
  cat("KEY_collapse_\n")
  cat(str(KEY_collapse))
  cat("\n")
  cat(sprintf(as.character(names(summary(KEY_collapse$CLASS_erythroid)))))
  cat("\n")
  cat(sprintf(as.character(summary(KEY_collapse$CLASS_erythroid))))
  cat("\n")
  cat(sprintf(as.character(names(summary(KEY_collapse$CLASS_mega)))))
  cat("\n")
  cat(sprintf(as.character(summary(KEY_collapse$CLASS_mega))))
  cat("\n")
  cat(sprintf(as.character(names(summary(KEY_collapse$CLASS_gran_mono)))))
  cat("\n")
  cat(sprintf(as.character(summary(KEY_collapse$CLASS_gran_mono))))
  cat("\n")
  cat(sprintf(as.character(names(summary(KEY_collapse$CLASS_lymph)))))
  cat("\n")
  cat(sprintf(as.character(summary(KEY_collapse$CLASS_lymph))))
  cat("\n")
  cat(sprintf(as.character(names(summary(KEY_collapse$CLASS_ALL_wbc)))))
  cat("\n")
  cat(sprintf(as.character(summary(KEY_collapse$CLASS_ALL_wbc))))
  cat("\n")
  
  
  ##### LINEAGE_plots  ----
  
  
  path5<-paste(out,'LINEAGE_plots','/', sep='')
  
  cat("path5\n")
  cat(sprintf(as.character(path5)))
  cat("\n")
  
  
  if (file.exists(path5)){
    
    
    
    
  } else {
    dir.create(file.path(path5))
    
  }
  
  
  ####    SAVE ----
  
  setwd(path5)
  
  # quit(status = 1)
  

  filename_1<-paste("KEY_collpase_Plus_Variant_lineage_CLASSIFICATION",".rds", sep='')
  saveRDS(file=filename_1,KEY_collapse)
  

}


phenotype_Lineage_ER_stacked_barplots_E_Plus_ASE_CLASS = function(option_list)
{
  
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
  
  #### Categories colors ----
  
  
  df_CSQ_colors<-readRDS(file=opt$CSQ_colors)
  
  df_CSQ_colors$color[df_CSQ_colors$VEP_DEF_LABELS == "TFBS"]<-"red"
  
  
  cat("df_CSQ_colors_0\n")
  cat(str(df_CSQ_colors))
  cat("\n")
  
  #### READ DATA----
  
  
  path5<-paste(out,'LINEAGE_plots','/', sep='')
  
  cat("path5\n")
  cat(sprintf(as.character(path5)))
  cat("\n")
  
  
  if (file.exists(path5)){
    
    
    
    
  } else {
    dir.create(file.path(path5))
    
  }
  
  
  
  setwd(path5)
  
  # quit(status = 1)
  
  
  filename_1<-paste("KEY_collpase_Plus_Variant_lineage_CLASSIFICATION",".rds", sep='')
  
  
  KEY_collapse<-readRDS(file=filename_1)
  
  cat("KEY_collapse_0\n")
  cat(str(KEY_collapse))
  cat("\n")
  
  
  check<-KEY_collapse[is.na(KEY_collapse$CLASS_erythroid),]
  
  cat("check_0\n")
  cat(str(check))
  cat("\n")
  cat(sprintf(as.character(names(summary(check$VEP_DEF_LABELS)))))
  cat("\n")
  cat(sprintf(as.character(summary(check$VEP_DEF_LABELS))))
  cat("\n")
  
  
  ##### exclude CTRLS -----
  
  KEY_collapse<-KEY_collapse[which(KEY_collapse$Label_2 == "ASSAYED_VARIANT"),]
  
  cat("KEY_collapse_excluded_CTRLS\n")
  cat(str(KEY_collapse))
  cat("\n")
  
  
  check<-KEY_collapse[is.na(KEY_collapse$CLASS_erythroid),]
  
  cat("check_0\n")
  cat(str(check))
  cat("\n")
  
  
  
  
  setwd(out)
  
  
  filename<-paste('df_Cell_colors','.rds',sep='')
  
  df_Cell_colors<-readRDS(file=filename)
  
  cat("df_Cell_colors_0\n")
  cat(str(df_Cell_colors))
  cat("\n")
  
  #### Categories colors ----
  
  
  df_CSQ_colors<-readRDS(file=opt$CSQ_colors)
  
  df_CSQ_colors$color[df_CSQ_colors$VEP_DEF_LABELS == "TFBS"]<-"red"
  
  
  cat("df_CSQ_colors_0\n")
  cat(str(df_CSQ_colors))
  cat("\n")
  
  ############################# TILES   ---------------------------
  
  
  CLASSES_array<-c("CLASS_erythroid","CLASS_mega","CLASS_gran_mono","CLASS_lymph","CLASS_ALL_wbc")
  
  
  List_SUPER_DEF_STATS<-list()
  
  List_CLASSIF_FOR_INTER_ANALYSIS<-list()
  
  

  for(iteration_CLASSES_array in 1:length(CLASSES_array))
  {
    CLASSES_array_sel<-CLASSES_array[iteration_CLASSES_array]
    
    cat("------------------------------------------------------------------------------------------------------------------------->\t")
    cat(sprintf(as.character(CLASSES_array_sel)))
    cat("\n")
    
    
    
    KEY_collapse_E_Plus_ASE_CLASS.dt<-data.table(KEY_collapse, key=c("Cell_Type","E_Plus_ASE_CLASS",CLASSES_array_sel))
    
    
    KEY_collapse_E_Plus_ASE_Fq<-as.data.frame(KEY_collapse_E_Plus_ASE_CLASS.dt[,.N,by=key(KEY_collapse_E_Plus_ASE_CLASS.dt)], stringsAsFactors=F)
    
    colnames(KEY_collapse_E_Plus_ASE_Fq)[which(colnames(KEY_collapse_E_Plus_ASE_Fq) == "N")]<-"instances"
    
    cat("KEY_collapse_E_Plus_ASE_Fq_\n")
    cat(str(KEY_collapse_E_Plus_ASE_Fq))
    cat("\n")
    
    
    
    KEY_collapse_E_Plus_ASE_CLASS_TOTAL.dt<-data.table(KEY_collapse, key=c("Cell_Type","E_Plus_ASE_CLASS"))
    
    
    KEY_collapse_E_Plus_ASE_CLASS_TOTAL<-as.data.frame(KEY_collapse_E_Plus_ASE_CLASS_TOTAL.dt[,.N,by=key(KEY_collapse_E_Plus_ASE_CLASS_TOTAL.dt)], stringsAsFactors=F)
    
    colnames(KEY_collapse_E_Plus_ASE_CLASS_TOTAL)[which(colnames(KEY_collapse_E_Plus_ASE_CLASS_TOTAL) == "N")]<-"TOTAL"
    
    cat("KEY_collapse_E_Plus_ASE_CLASS_TOTAL_\n")
    cat(str(KEY_collapse_E_Plus_ASE_CLASS_TOTAL))
    cat("\n")
    
    KEY_collapse_E_Plus_ASE_Fq<-merge(KEY_collapse_E_Plus_ASE_Fq,
                                    KEY_collapse_E_Plus_ASE_CLASS_TOTAL,
                                    by=c("Cell_Type","E_Plus_ASE_CLASS"))
    
    # KEY_collapse_E_Plus_ASE_Fq$Perc<-round(100*(KEY_collapse_E_Plus_ASE_Fq$instances/KEY_collapse_E_Plus_ASE_Fq$TOTAL),1)
    
    
    cat("KEY_collapse_E_Plus_ASE_Fq_\n")
    cat(str(KEY_collapse_E_Plus_ASE_Fq))
    cat("\n")
    
    
    
    E_Plus_ASE_CLASS_array<-levels(KEY_collapse_E_Plus_ASE_Fq$E_Plus_ASE_CLASS)
    
    
    cat("E_Plus_ASE_CLASS_array_\n")
    cat(str(E_Plus_ASE_CLASS_array))
    cat("\n")
    
    
    
    accepted_values<-NULL
    AT_LEAST_category<-NULL
    
    list_results<-list()
    
    for(iteration_E_Plus_ASE_CLASS_array in 1:length(E_Plus_ASE_CLASS_array))
    {
      E_Plus_ASE_CLASS_array_sel<-E_Plus_ASE_CLASS_array[iteration_E_Plus_ASE_CLASS_array]
      
      cat("------------------------>\t")
      cat(sprintf(as.character(E_Plus_ASE_CLASS_array_sel)))
      cat("\n")
      
      if(E_Plus_ASE_CLASS_array_sel == "0")
      {
        
        accepted_values<-"0"
        AT_LEAST_category<-"0"
      }
      if(E_Plus_ASE_CLASS_array_sel == "1")
      {
        
        accepted_values<-c("1","2","3 or >")
        AT_LEAST_category<-"AT LEAST 1"
        
      }
      if(E_Plus_ASE_CLASS_array_sel == "2")
      {
        
        accepted_values<-c("2","3 or >")
        AT_LEAST_category<-"AT LEAST 2"
      }
      if(E_Plus_ASE_CLASS_array_sel == "3 or >")
      {
        
        accepted_values<-c("3 or >")
        AT_LEAST_category<-"AT LEAST 3"
      }
      
      cat("accepted_values_\n")
      cat(str(accepted_values))
      cat("\n")
      
      KEY_collapse_E_Plus_ASE_Fq_sel<-KEY_collapse_E_Plus_ASE_Fq[which(KEY_collapse_E_Plus_ASE_Fq$E_Plus_ASE_CLASS%in%accepted_values),]
      
      cat("KEY_collapse_E_Plus_ASE_Fq_sel_\n")
      cat(str(KEY_collapse_E_Plus_ASE_Fq_sel))
      cat("\n")
      
      
      #### collapse AT LEAST ----
      
      AT_LEAST_COLLAPSE.dt<-data.table(KEY_collapse_E_Plus_ASE_Fq_sel, key=c("Cell_Type",CLASSES_array_sel))
      
      
      AT_LEAST_COLLAPSE_Fq<-as.data.frame(AT_LEAST_COLLAPSE.dt[,.(AT_LEAST_instances=sum(instances)),by=key(AT_LEAST_COLLAPSE.dt)], stringsAsFactors=F)
      
      
      cat("AT_LEAST_COLLAPSE_Fq_\n")
      cat(str(AT_LEAST_COLLAPSE_Fq))
      cat("\n")
      
      
      
      AT_LEAST_COLLAPSE_TOTAL.dt<-data.table(KEY_collapse_E_Plus_ASE_Fq_sel, key=c("Cell_Type"))
      
      
      AT_LEAST_COLLAPSE_TOTAL<-as.data.frame(AT_LEAST_COLLAPSE_TOTAL.dt[,.(AT_LEAST_TOTAL=sum(instances)),by=key(AT_LEAST_COLLAPSE_TOTAL.dt)], stringsAsFactors=F)
      
      
      cat("AT_LEAST_COLLAPSE_TOTAL_\n")
      cat(str(AT_LEAST_COLLAPSE_TOTAL))
      cat("\n")
      
      AT_LEAST_COLLAPSE_Fq<-merge(AT_LEAST_COLLAPSE_Fq,
                                  AT_LEAST_COLLAPSE_TOTAL,
                                  by=c("Cell_Type"))
      
      AT_LEAST_COLLAPSE_Fq$AT_LEAST_category<-AT_LEAST_category
      
      cat("AT_LEAST_COLLAPSE_Fq_\n")
      cat(str(AT_LEAST_COLLAPSE_Fq))
      cat("\n")
      
      list_results[[iteration_E_Plus_ASE_CLASS_array]]<-AT_LEAST_COLLAPSE_Fq
      
      
      
      
    }# iteration_E_Plus_ASE_CLASS_array
    
    
    
    
    
    TABLE_AT_LEAST_E_Plus_ASE = unique(as.data.frame(data.table::rbindlist(list_results, fill = T)))
    
    TABLE_AT_LEAST_E_Plus_ASE$AT_LEAST_category<-factor(TABLE_AT_LEAST_E_Plus_ASE$AT_LEAST_category,
                                                      levels = c("0","AT LEAST 1","AT LEAST 2","AT LEAST 3"),
                                                      ordered=T)
    
    TABLE_AT_LEAST_E_Plus_ASE$Perc<-round(100*(TABLE_AT_LEAST_E_Plus_ASE$AT_LEAST_instances/TABLE_AT_LEAST_E_Plus_ASE$AT_LEAST_TOTAL),1)
    
    TABLE_AT_LEAST_E_Plus_ASE<-TABLE_AT_LEAST_E_Plus_ASE[order(TABLE_AT_LEAST_E_Plus_ASE$Cell_Type,TABLE_AT_LEAST_E_Plus_ASE$AT_LEAST_category),]
    
    cat("TABLE_AT_LEAST_E_Plus_ASE_0\n")
    cat(str(TABLE_AT_LEAST_E_Plus_ASE))
    cat("\n")
    
    List_CLASSIF_FOR_INTER_ANALYSIS[[iteration_CLASSES_array]]<-TABLE_AT_LEAST_E_Plus_ASE
    
    
    #### CT LOOP for graphs and STATS ----
    
    Cell_Type_array<-levels(TABLE_AT_LEAST_E_Plus_ASE$Cell_Type)
    
    
    cat("Cell_Type_array_0\n")
    cat(str(Cell_Type_array))
    cat("\n")
    
    
    list_DEF<-list()
    
    
    list_ChiSq_DEF<-list()
    
    for(i in 1:length(Cell_Type_array))
    {
      Cell_Type_array_sel<-Cell_Type_array[i]
      
      cat("------------------------------------------------------------------------------------------------>\t")
      cat(sprintf(as.character(Cell_Type_array_sel)))
      cat("\n")
      
      TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel<-TABLE_AT_LEAST_E_Plus_ASE[which(TABLE_AT_LEAST_E_Plus_ASE$Cell_Type == Cell_Type_array_sel),]
      TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel<-droplevels(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel)
      
      indx.sel<-which(colnames(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel) == CLASSES_array_sel)
      
      cat("indx.sel_0\n")
      cat(str(indx.sel))
      cat("\n")
      
      colnames.sel<-colnames(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel)[which(colnames(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel) == CLASSES_array_sel)]
      
      cat("colnames.sel\n")
      cat(str(colnames.sel))
      cat("\n")
      
      
      
      #### STATS ----
      
      
      indx.int_2<-which(colnames(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel) =="AT_LEAST_category")
      
      
      if(dim(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel)[1] > 0)
      {
        
        list_ChiSq_category<-list()
        for(k in 1:length(levels(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel[,indx.int_2])))
        {
          level_sel<-levels(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel[,indx.int_2])[k]
          
          cat("--->\t")
          cat(sprintf(as.character(level_sel)))
          cat("\n")
          
          if(level_sel != "0")
          {
            
            TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0<-TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel[which(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel[,indx.int_2] == "0"),]
            
            # cat("TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0\n")
            # cat(str(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0))
            # cat("\n")
            
            TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel<-TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel[which(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel[,indx.int_2] == level_sel),]
            
            # cat("TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel\n")
            # cat(str(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel))
            # cat("\n")
            
            levels_LINEAGE_CLASSES<-levels(TABLE_AT_LEAST_E_Plus_ASE[,indx.sel])
            
            list_LINEAGE_CLASSES<-list()
            
            for(iteration_levels_LINEAGE_CLASSES in 1:length(levels_LINEAGE_CLASSES))
            {
              
              
              levels_LINEAGE_CLASSES_sel<-levels_LINEAGE_CLASSES[iteration_levels_LINEAGE_CLASSES]
              
              cat("->\t")
              cat(sprintf(as.character(levels_LINEAGE_CLASSES_sel)))
              cat("\n")
              
              
              TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0_LINEAGE_CLASSES_sel<-TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0[which(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0[,indx.sel] == levels_LINEAGE_CLASSES_sel),]
              
              # cat("TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0_LINEAGE_CLASSES_sel\n")
              # cat(str(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0_LINEAGE_CLASSES_sel))
              # cat("\n")
              
              TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_LINEAGE_CLASSES_sel<-TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel[which(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel[,indx.sel] == levels_LINEAGE_CLASSES_sel),]
              
              # cat("TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_LINEAGE_CLASSES_sel\n")
              # cat(str(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_LINEAGE_CLASSES_sel))
              # cat("\n")
              
              
              TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0_NOT_LINEAGE_CLASSES_sel<-TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0[which(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0[,indx.sel] != levels_LINEAGE_CLASSES_sel),]
              
              # cat("TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0_NOT_LINEAGE_CLASSES_sel\n")
              # cat(str(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0_NOT_LINEAGE_CLASSES_sel))
              # cat("\n")
              
              TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_NOT_LINEAGE_CLASSES_sel<-TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel[which(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel[,indx.sel] != levels_LINEAGE_CLASSES_sel),]
              
              # cat("TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_NOT_LINEAGE_CLASSES_sel\n")
              # cat(str(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_NOT_LINEAGE_CLASSES_sel))
              # cat("\n")
              
              
              
              if(dim(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0_LINEAGE_CLASSES_sel)[1]>0)
              {
                if(dim(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_LINEAGE_CLASSES_sel)[1]>0)
                {
                  cat("->Hello_world_1\n")
                  
                  
                  A.df<-as.data.frame(rbind(cbind(sum(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0_NOT_LINEAGE_CLASSES_sel$AT_LEAST_instances),sum(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_NOT_LINEAGE_CLASSES_sel$AT_LEAST_instances)),
                                            cbind(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0_LINEAGE_CLASSES_sel$AT_LEAST_instances,TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_LINEAGE_CLASSES_sel$AT_LEAST_instances)),
                                      stringsAsFactors=F)
                  
                  colnames(A.df)<-c("0",level_sel)
                  row.names(A.df)<-c(paste("NOT",levels_LINEAGE_CLASSES_sel,sep='_'),
                                     levels_LINEAGE_CLASSES_sel)
                  
                  
                  # cat("A.df\n")
                  # cat(str(A.df))
                  # cat("\n")
                  
                  tab.chisq.test<-chisq.test(A.df,correct = TRUE)
                  
                  # cat("tab.chisq.test\n")
                  # cat(str(tab.chisq.test))
                  # cat("\n")
                  
                  log_pval<-as.numeric(round(-1*log10(tab.chisq.test$p.value),2))
                  
                  # cat("log_pval\n")
                  # cat(str(log_pval))
                  # cat("\n")
                  
                  a.df<-as.data.frame(cbind(level_sel,levels_LINEAGE_CLASSES_sel,log_pval), stringsAsFactors=F)
                  
                  colnames(a.df)<-c("Category","COMPARISON_category","log_pval")
                  
                  # cat("a.df\n")
                  # cat(str(a.df))
                  # cat("\n")
                  
                  list_LINEAGE_CLASSES[[iteration_levels_LINEAGE_CLASSES]]<-a.df
                  
                  
                }else{
                  cat("->Hello_world_2\n")
                  
                  
                  A.df<-as.data.frame(rbind(cbind(sum(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0_NOT_LINEAGE_CLASSES_sel$AT_LEAST_instances),sum(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_NOT_LINEAGE_CLASSES_sel$AT_LEAST_instances)),
                                            cbind(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0_LINEAGE_CLASSES_sel$AT_LEAST_instances,0)),
                                      stringsAsFactors=F)
                  
                  colnames(A.df)<-c("0",level_sel)
                  row.names(A.df)<-c(paste("NOT",levels_LINEAGE_CLASSES_sel,sep='_'),
                                     levels_LINEAGE_CLASSES_sel)
                  
                  
                  # cat("A.df\n")
                  # cat(str(A.df))
                  # cat("\n")
                  
                  tab.chisq.test<-chisq.test(A.df,correct = TRUE)
                  
                  # cat("tab.chisq.test\n")
                  # cat(str(tab.chisq.test))
                  # cat("\n")
                  
                  log_pval<-as.numeric(round(-1*log10(tab.chisq.test$p.value),2))
                  
                  # cat("log_pval\n")
                  # cat(str(log_pval))
                  # cat("\n")
                  
                  a.df<-as.data.frame(cbind(level_sel,levels_LINEAGE_CLASSES_sel,log_pval), stringsAsFactors=F)
                  
                  colnames(a.df)<-c("Category","COMPARISON_category","log_pval")
                  
                  # cat("a.df\n")
                  # cat(str(a.df))
                  # cat("\n")
                  
                  list_LINEAGE_CLASSES[[iteration_levels_LINEAGE_CLASSES]]<-a.df
                  
                  
                }# dim(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_LINEAGE_CLASSES_sel)[1]>0
                
              }else{
                
                if(dim(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_LINEAGE_CLASSES_sel)[1]>0)
                {
                  
                  cat("->Hello_world_3\n")
                  
                  A.df<-as.data.frame(rbind(cbind(sum(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0_NOT_LINEAGE_CLASSES_sel$AT_LEAST_instances),sum(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_NOT_LINEAGE_CLASSES_sel$AT_LEAST_instances)),
                                            cbind(0,TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_LINEAGE_CLASSES_sel$AT_LEAST_instances)),
                                      stringsAsFactors=F)
                  
                  colnames(A.df)<-c("0",level_sel)
                  row.names(A.df)<-c(paste("NOT",levels_LINEAGE_CLASSES_sel,sep='_'),
                                     levels_LINEAGE_CLASSES_sel)
                  
                  
                  # cat("A.df\n")
                  # cat(str(A.df))
                  # cat("\n")
                  
                  tab.chisq.test<-chisq.test(A.df,correct = TRUE)
                  
                  # cat("tab.chisq.test\n")
                  # cat(str(tab.chisq.test))
                  # cat("\n")
                  
                  log_pval<-as.numeric(round(-1*log10(tab.chisq.test$p.value),2))
                  
                  # cat("log_pval\n")
                  # cat(str(log_pval))
                  # cat("\n")
                  
                  a.df<-as.data.frame(cbind(level_sel,levels_LINEAGE_CLASSES_sel,log_pval), stringsAsFactors=F)
                  
                  colnames(a.df)<-c("Category","COMPARISON_category","log_pval")
                  
                  # cat("a.df\n")
                  # cat(str(a.df))
                  # cat("\n")
                  
                  list_LINEAGE_CLASSES[[iteration_levels_LINEAGE_CLASSES]]<-a.df
                }else{
                  
                  cat("->Hello_world_4\n")
                  
                  A.df<-as.data.frame(rbind(cbind(sum(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0_NOT_LINEAGE_CLASSES_sel$AT_LEAST_instances),sum(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_NOT_LINEAGE_CLASSES_sel$AT_LEAST_instances)),
                                            cbind(0,0)),
                                      stringsAsFactors=F)
                  
                  colnames(A.df)<-c("0",level_sel)
                  row.names(A.df)<-c(paste("NOT",levels_LINEAGE_CLASSES_sel,sep='_'),
                                     levels_LINEAGE_CLASSES_sel)
                  
                  
                  # cat("A.df\n")
                  # cat(str(A.df))
                  # cat("\n")
                  
                  tab.chisq.test<-chisq.test(A.df,correct = TRUE)
                  
                  # cat("tab.chisq.test\n")
                  # cat(str(tab.chisq.test))
                  # cat("\n")
                  
                  log_pval<-as.numeric(round(-1*log10(tab.chisq.test$p.value),2))
                  
                  # cat("log_pval\n")
                  # cat(str(log_pval))
                  # cat("\n")
                  
                  a.df<-as.data.frame(cbind(level_sel,levels_LINEAGE_CLASSES_sel,log_pval), stringsAsFactors=F)
                  
                  colnames(a.df)<-c("Category","COMPARISON_category","log_pval")
                  
                  # cat("a.df\n")
                  # cat(str(a.df))
                  # cat("\n")
                  
                  list_LINEAGE_CLASSES[[iteration_levels_LINEAGE_CLASSES]]<-a.df
                  
                }#dim(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_LINEAGE_CLASSES_sel)[1]>0
                
                
                
              }# dim(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0_LINEAGE_CLASSES_sel)[1]>0
              
              
            }# iteration_levels_LINEAGE_CLASSES
            
            CHI_Sq_LINEAGE_CLASSES = unique(as.data.frame(data.table::rbindlist(list_LINEAGE_CLASSES, fill = T)))
            
            # cat("CHI_Sq_LINEAGE_CLASSES\n")
            # cat(str(CHI_Sq_LINEAGE_CLASSES))
            # cat("\n")
            
            list_ChiSq_category[[k]]<-CHI_Sq_LINEAGE_CLASSES
            
            # ############################################################################################################
            # quit(status = 1) 
            
            
          }#level_sel != "0"
        }# k
      }# dim(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel)[1] > 0
      
      CHI_Sq_CSQ_CATEGORY = unique(as.data.frame(data.table::rbindlist(list_ChiSq_category, fill = T)))
      
      CHI_Sq_CSQ_CATEGORY$Cell_Type<-Cell_Type_array_sel
      
      # cat("CHI_Sq_CSQ_CATEGORY\n")
      # cat(str(CHI_Sq_CSQ_CATEGORY))
      # cat("\n")
      # 
      
      
      list_ChiSq_DEF[[i]]<-CHI_Sq_CSQ_CATEGORY
      
      # ############################################################################################################
      # quit(status = 1) 
      
    }# i Cell_Type
    

    CHI_Sq_DEF = unique(as.data.frame(data.table::rbindlist(list_ChiSq_DEF, fill = T)))
    
    
    cat("CHI_Sq_DEF_0\n")
    cat(str(CHI_Sq_DEF))
    cat("\n")
    
    
    CHI_Sq_DEF$GROUP<-CLASSES_array_sel
    
    
    cat("CHI_Sq_DEF_1\n")
    cat(str(CHI_Sq_DEF))
    cat("\n")
    
    
    
    List_SUPER_DEF_STATS[[iteration_CLASSES_array]]<-CHI_Sq_DEF
  
    
  }#iteration_CLASSES_array
  
  
  

  CLASSIF_FOR_INTER = unique(as.data.frame(data.table::rbindlist(List_CLASSIF_FOR_INTER_ANALYSIS, fill = T)))
  
  cat("CLASSIF_FOR_INTER_0\n")
  cat(str(CLASSIF_FOR_INTER))
  cat("\n")
  
  
  CLASSIF_FOR_INTER.m<-melt(CLASSIF_FOR_INTER, id.vars=c("Cell_Type","AT_LEAST_instances","AT_LEAST_TOTAL","AT_LEAST_category","Perc"))
  
  # measure.vars=c("CLASS_erythroid","CLASS_mega","CLASS_gran_mono","CLASS_lymph","CLASS_ALL_wbc"),
  # value.name="Perc"
  
  
  CLASSIF_FOR_INTER.m<-CLASSIF_FOR_INTER.m[!is.na(CLASSIF_FOR_INTER.m$value),]
  
  CLASSIF_FOR_INTER.m<-droplevels(CLASSIF_FOR_INTER.m[-which(CLASSIF_FOR_INTER.m$variable == "CLASS_ALL_wbc"),])
  
  CLASSIF_FOR_INTER.m$value<-factor(CLASSIF_FOR_INTER.m$value,
                                    levels=c("OTHER","erythroid_lineage","mega_lineage","gran_mono_lineage","lymph_lineage"),
                                    ordered=T)
  
  
  
  cat("CLASSIF_FOR_INTER.m_0\n")
  cat(str(CLASSIF_FOR_INTER.m))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(CLASSIF_FOR_INTER.m$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(CLASSIF_FOR_INTER.m$variable)))))
  cat("\n")
  
  cat(sprintf(as.character(names(summary(as.factor(CLASSIF_FOR_INTER.m$value))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(CLASSIF_FOR_INTER.m$value)))))
  cat("\n")
  
 
    

#### PANEL GRAPH ----    
    
  breaks.y<-seq(0,100,by=25)
  labels.y<-as.character(breaks.y)
 
    
  
  colors<-c("gray","firebrick1","brown","magenta","dodgerblue")
    


  cat("GRAPH_START\n")


  graph<-ggplot(data=CLASSIF_FOR_INTER.m,
                aes(x=AT_LEAST_category,
                    y=Perc,
                    fill=value)) +
    geom_bar(stat="identity",colour='black')+
    theme(axis.title.y=element_text(size=18, family="sans"),
          axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
          axis.text.x=element_text(angle=45,size=18, color="black", family="sans"))+
    scale_y_continuous(name=paste("Perentage"),breaks=breaks.y,labels=labels.y,
                       limits=c(breaks.y[1],breaks.y[length(breaks.y)]+1))+
    scale_x_discrete(name=NULL, drop=F)+
    theme_classic()+
    scale_fill_manual(values=colors, drop=F)+
    theme(legend.position="bottom")+
    facet_grid(Cell_Type ~ variable)+
    ggeasy::easy_center_title()
  
  cat("GRAPH_END\n")
  

  # setwd(out)
  # 
  # svgname<-paste(paste("TEST",".svg",sep=''))
  # makesvg = TRUE
  # 
  # if (makesvg == TRUE)
  # {
  #   ggsave(svgname, plot= graph,
  #          device="svg",
  #          height=10, width=12)
  # }
  # 
  # ########################################################################################
  # quit(status = 1)
#   

  setwd(path5)

  svgname<-paste(paste("E_Plus_ASE_single_graph",".svg",sep=''))
  makesvg = TRUE

  if (makesvg == TRUE)
  {
    ggsave(svgname, plot= graph,
           device="svg",
           height=10, width=12)
  }
  
  CHI_Sq_FINAL_DEF = unique(as.data.frame(data.table::rbindlist(List_SUPER_DEF_STATS, fill = T)))
  
  cat("CHI_Sq_FINAL_DEF_0\n")
  cat(str(CHI_Sq_FINAL_DEF))
  cat("\n")
  
  CHI_Sq_FINAL_DEF$COMPARISON_category<-factor(CHI_Sq_FINAL_DEF$COMPARISON_category,
                                               levels=c("OTHER","erythroid_lineage","mega_lineage","gran_mono_lineage","lymph_lineage",
                                                        "CLASS_ALL_wbc"),
                                               ordered=T)
  
  CHI_Sq_FINAL_DEF$Category<-factor(CHI_Sq_FINAL_DEF$Category,
                                    levels = c("0","AT LEAST 1","AT LEAST 2","AT LEAST 3"),
                                    ordered=T)
  
  CHI_Sq_FINAL_DEF$log_pval<-as.numeric(CHI_Sq_FINAL_DEF$log_pval)
  
  CHI_Sq_FINAL_DEF$Cell_Type<-factor(CHI_Sq_FINAL_DEF$Cell_Type,
                                     levels = c("K562","CHRF","HL60","THP1"),
                                     ordered=T)
  CHI_Sq_FINAL_DEF$GROUP<-factor(CHI_Sq_FINAL_DEF$GROUP,
                                 levels=c("CLASS_erythroid","CLASS_mega","CLASS_gran_mono","CLASS_lymph",
                                          "CLASS_ALL_wbc"),
                                 ordered=T)
  
  
  CHI_Sq_FINAL_DEF<-CHI_Sq_FINAL_DEF[order(CHI_Sq_FINAL_DEF$Cell_Type,CHI_Sq_FINAL_DEF$GROUP,CHI_Sq_FINAL_DEF$Category,CHI_Sq_FINAL_DEF$COMPARISON_category),]
  
  cat("CHI_Sq_FINAL_DEF_1\n")
  cat(str(CHI_Sq_FINAL_DEF))
  cat("\n")
  
  CHI_Sq_FINAL_DEF<-CHI_Sq_FINAL_DEF[-which(CHI_Sq_FINAL_DEF$COMPARISON_category == "OTHER"),]
  
  CHI_Sq_FINAL_DEF<-droplevels(CHI_Sq_FINAL_DEF)
  
  cat("CHI_Sq_FINAL_DEF_2\n")
  cat(str(CHI_Sq_FINAL_DEF))
  cat("\n")
  
  #### SAVE----
  
  
  setwd(path5)
  
  saveRDS(CLASSIF_FOR_INTER,file="CLASSIF_FOR_INTER_E_Plus_ASE_CLASS_stats.rds")
  
  
  
  saveRDS(CHI_Sq_FINAL_DEF,file="Lineage_E_Plus_ASE_CLASS_stats.rds")
  
  # saveRDS(TABLE_AT_LEAST_E_Plus_ASE,file="Lineage_E_Plus_ASE_df.rds")
  
}

phenotype_Lineage_ER_stacked_barplots_enhancer_CLASS = function(option_list)
{
  
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
  
  #### Categories colors ----
  
  
  df_CSQ_colors<-readRDS(file=opt$CSQ_colors)
  
  df_CSQ_colors$color[df_CSQ_colors$VEP_DEF_LABELS == "TFBS"]<-"red"
  
  
  cat("df_CSQ_colors_0\n")
  cat(str(df_CSQ_colors))
  cat("\n")
  
  #### READ DATA----
  
  
  path5<-paste(out,'LINEAGE_plots','/', sep='')
  
  cat("path5\n")
  cat(sprintf(as.character(path5)))
  cat("\n")
  
  
  if (file.exists(path5)){
    
    
    
    
  } else {
    dir.create(file.path(path5))
    
  }
  
  
  
  setwd(path5)
  
  # quit(status = 1)
  
  
  filename_1<-paste("KEY_collpase_Plus_Variant_lineage_CLASSIFICATION",".rds", sep='')
  
  
  KEY_collapse<-readRDS(file=filename_1)
  
  cat("KEY_collapse_0\n")
  cat(str(KEY_collapse))
  cat("\n")
  
  
  check<-KEY_collapse[is.na(KEY_collapse$CLASS_erythroid),]
  
  cat("check_0\n")
  cat(str(check))
  cat("\n")
  cat(sprintf(as.character(names(summary(check$VEP_DEF_LABELS)))))
  cat("\n")
  cat(sprintf(as.character(summary(check$VEP_DEF_LABELS))))
  cat("\n")
  
  
  ##### exclude CTRLS -----
  
  KEY_collapse<-KEY_collapse[which(KEY_collapse$Label_2 == "ASSAYED_VARIANT"),]
  
  cat("KEY_collapse_excluded_CTRLS\n")
  cat(str(KEY_collapse))
  cat("\n")
  
  
  check<-KEY_collapse[is.na(KEY_collapse$CLASS_erythroid),]
  
  cat("check_0\n")
  cat(str(check))
  cat("\n")
  
  
  
  
  setwd(out)
  
  
  filename<-paste('df_Cell_colors','.rds',sep='')
  
  df_Cell_colors<-readRDS(file=filename)
  
  cat("df_Cell_colors_0\n")
  cat(str(df_Cell_colors))
  cat("\n")
  
  #### Categories colors ----
  
  
  df_CSQ_colors<-readRDS(file=opt$CSQ_colors)
  
  df_CSQ_colors$color[df_CSQ_colors$VEP_DEF_LABELS == "TFBS"]<-"red"
  
  
  cat("df_CSQ_colors_0\n")
  cat(str(df_CSQ_colors))
  cat("\n")
  
  ############################# TILES   ---------------------------
  
  
  CLASSES_array<-c("CLASS_erythroid","CLASS_mega","CLASS_gran_mono","CLASS_lymph","CLASS_ALL_wbc")
  
  
  List_SUPER_DEF_STATS<-list()
  
  List_CLASSIF_FOR_INTER_ANALYSIS<-list()
  
  
  
  for(iteration_CLASSES_array in 1:length(CLASSES_array))
  {
    CLASSES_array_sel<-CLASSES_array[iteration_CLASSES_array]
    
    cat("------------------------------------------------------------------------------------------------------------------------->\t")
    cat(sprintf(as.character(CLASSES_array_sel)))
    cat("\n")
    
    
    
    KEY_collapse_enhancer_CLASS.dt<-data.table(KEY_collapse, key=c("Cell_Type","enhancer_CLASS",CLASSES_array_sel))
    
    
    KEY_collapse_enhancer_Fq<-as.data.frame(KEY_collapse_enhancer_CLASS.dt[,.N,by=key(KEY_collapse_enhancer_CLASS.dt)], stringsAsFactors=F)
    
    colnames(KEY_collapse_enhancer_Fq)[which(colnames(KEY_collapse_enhancer_Fq) == "N")]<-"instances"
    
    cat("KEY_collapse_enhancer_Fq_\n")
    cat(str(KEY_collapse_enhancer_Fq))
    cat("\n")
    
    
    
    KEY_collapse_enhancer_CLASS_TOTAL.dt<-data.table(KEY_collapse, key=c("Cell_Type","enhancer_CLASS"))
    
    
    KEY_collapse_enhancer_CLASS_TOTAL<-as.data.frame(KEY_collapse_enhancer_CLASS_TOTAL.dt[,.N,by=key(KEY_collapse_enhancer_CLASS_TOTAL.dt)], stringsAsFactors=F)
    
    colnames(KEY_collapse_enhancer_CLASS_TOTAL)[which(colnames(KEY_collapse_enhancer_CLASS_TOTAL) == "N")]<-"TOTAL"
    
    cat("KEY_collapse_enhancer_CLASS_TOTAL_\n")
    cat(str(KEY_collapse_enhancer_CLASS_TOTAL))
    cat("\n")
    
    KEY_collapse_enhancer_Fq<-merge(KEY_collapse_enhancer_Fq,
                                      KEY_collapse_enhancer_CLASS_TOTAL,
                                      by=c("Cell_Type","enhancer_CLASS"))
    
    # KEY_collapse_enhancer_Fq$Perc<-round(100*(KEY_collapse_enhancer_Fq$instances/KEY_collapse_enhancer_Fq$TOTAL),1)
    
    
    cat("KEY_collapse_enhancer_Fq_\n")
    cat(str(KEY_collapse_enhancer_Fq))
    cat("\n")
    
    
    
    enhancer_CLASS_array<-levels(KEY_collapse_enhancer_Fq$enhancer_CLASS)
    
    
    cat("enhancer_CLASS_array_\n")
    cat(str(enhancer_CLASS_array))
    cat("\n")
    
    
    
    accepted_values<-NULL
    AT_LEAST_category<-NULL
    
    list_results<-list()
    
    for(iteration_enhancer_CLASS_array in 1:length(enhancer_CLASS_array))
    {
      enhancer_CLASS_array_sel<-enhancer_CLASS_array[iteration_enhancer_CLASS_array]
      
      cat("------------------------>\t")
      cat(sprintf(as.character(enhancer_CLASS_array_sel)))
      cat("\n")
      
      if(enhancer_CLASS_array_sel == "0")
      {
        
        accepted_values<-"0"
        AT_LEAST_category<-"0"
      }
      if(enhancer_CLASS_array_sel == "1")
      {
        
        accepted_values<-c("1","2","3 or >")
        AT_LEAST_category<-"AT LEAST 1"
        
      }
      if(enhancer_CLASS_array_sel == "2")
      {
        
        accepted_values<-c("2","3 or >")
        AT_LEAST_category<-"AT LEAST 2"
      }
      if(enhancer_CLASS_array_sel == "3 or >")
      {
        
        accepted_values<-c("3 or >")
        AT_LEAST_category<-"AT LEAST 3"
      }
      
      cat("accepted_values_\n")
      cat(str(accepted_values))
      cat("\n")
      
      KEY_collapse_enhancer_Fq_sel<-KEY_collapse_enhancer_Fq[which(KEY_collapse_enhancer_Fq$enhancer_CLASS%in%accepted_values),]
      
      cat("KEY_collapse_enhancer_Fq_sel_\n")
      cat(str(KEY_collapse_enhancer_Fq_sel))
      cat("\n")
      
      
      #### collapse AT LEAST ----
      
      AT_LEAST_COLLAPSE.dt<-data.table(KEY_collapse_enhancer_Fq_sel, key=c("Cell_Type",CLASSES_array_sel))
      
      
      AT_LEAST_COLLAPSE_Fq<-as.data.frame(AT_LEAST_COLLAPSE.dt[,.(AT_LEAST_instances=sum(instances)),by=key(AT_LEAST_COLLAPSE.dt)], stringsAsFactors=F)
      
      
      cat("AT_LEAST_COLLAPSE_Fq_\n")
      cat(str(AT_LEAST_COLLAPSE_Fq))
      cat("\n")
      
      
      
      AT_LEAST_COLLAPSE_TOTAL.dt<-data.table(KEY_collapse_enhancer_Fq_sel, key=c("Cell_Type"))
      
      
      AT_LEAST_COLLAPSE_TOTAL<-as.data.frame(AT_LEAST_COLLAPSE_TOTAL.dt[,.(AT_LEAST_TOTAL=sum(instances)),by=key(AT_LEAST_COLLAPSE_TOTAL.dt)], stringsAsFactors=F)
      
      
      cat("AT_LEAST_COLLAPSE_TOTAL_\n")
      cat(str(AT_LEAST_COLLAPSE_TOTAL))
      cat("\n")
      
      AT_LEAST_COLLAPSE_Fq<-merge(AT_LEAST_COLLAPSE_Fq,
                                  AT_LEAST_COLLAPSE_TOTAL,
                                  by=c("Cell_Type"))
      
      AT_LEAST_COLLAPSE_Fq$AT_LEAST_category<-AT_LEAST_category
      
      cat("AT_LEAST_COLLAPSE_Fq_\n")
      cat(str(AT_LEAST_COLLAPSE_Fq))
      cat("\n")
      
      list_results[[iteration_enhancer_CLASS_array]]<-AT_LEAST_COLLAPSE_Fq
      
      
      
      
    }# iteration_enhancer_CLASS_array
    
    
    
    
    
    TABLE_AT_LEAST_enhancer = unique(as.data.frame(data.table::rbindlist(list_results, fill = T)))
    
    TABLE_AT_LEAST_enhancer$AT_LEAST_category<-factor(TABLE_AT_LEAST_enhancer$AT_LEAST_category,
                                                        levels = c("0","AT LEAST 1","AT LEAST 2","AT LEAST 3"),
                                                        ordered=T)
    
    TABLE_AT_LEAST_enhancer$Perc<-round(100*(TABLE_AT_LEAST_enhancer$AT_LEAST_instances/TABLE_AT_LEAST_enhancer$AT_LEAST_TOTAL),1)
    
    TABLE_AT_LEAST_enhancer<-TABLE_AT_LEAST_enhancer[order(TABLE_AT_LEAST_enhancer$Cell_Type,TABLE_AT_LEAST_enhancer$AT_LEAST_category),]
    
    cat("TABLE_AT_LEAST_enhancer_0\n")
    cat(str(TABLE_AT_LEAST_enhancer))
    cat("\n")
    
    List_CLASSIF_FOR_INTER_ANALYSIS[[iteration_CLASSES_array]]<-TABLE_AT_LEAST_enhancer
    
    
    #### CT LOOP for graphs and STATS ----
    
    Cell_Type_array<-levels(TABLE_AT_LEAST_enhancer$Cell_Type)
    
    
    cat("Cell_Type_array_0\n")
    cat(str(Cell_Type_array))
    cat("\n")
    
    
    list_DEF<-list()
    
    
    list_ChiSq_DEF<-list()
    
    for(i in 1:length(Cell_Type_array))
    {
      Cell_Type_array_sel<-Cell_Type_array[i]
      
      cat("------------------------------------------------------------------------------------------------>\t")
      cat(sprintf(as.character(Cell_Type_array_sel)))
      cat("\n")
      
      TABLE_AT_LEAST_enhancer_Cell_Type_sel<-TABLE_AT_LEAST_enhancer[which(TABLE_AT_LEAST_enhancer$Cell_Type == Cell_Type_array_sel),]
      TABLE_AT_LEAST_enhancer_Cell_Type_sel<-droplevels(TABLE_AT_LEAST_enhancer_Cell_Type_sel)
      
      indx.sel<-which(colnames(TABLE_AT_LEAST_enhancer_Cell_Type_sel) == CLASSES_array_sel)
      
      cat("indx.sel_0\n")
      cat(str(indx.sel))
      cat("\n")
      
      colnames.sel<-colnames(TABLE_AT_LEAST_enhancer_Cell_Type_sel)[which(colnames(TABLE_AT_LEAST_enhancer_Cell_Type_sel) == CLASSES_array_sel)]
      
      cat("colnames.sel\n")
      cat(str(colnames.sel))
      cat("\n")
      
      
      
      #### STATS ----
      
      
      indx.int_2<-which(colnames(TABLE_AT_LEAST_enhancer_Cell_Type_sel) =="AT_LEAST_category")
      
      
      if(dim(TABLE_AT_LEAST_enhancer_Cell_Type_sel)[1] > 0)
      {
        
        list_ChiSq_category<-list()
        for(k in 1:length(levels(TABLE_AT_LEAST_enhancer_Cell_Type_sel[,indx.int_2])))
        {
          level_sel<-levels(TABLE_AT_LEAST_enhancer_Cell_Type_sel[,indx.int_2])[k]
          
          cat("--->\t")
          cat(sprintf(as.character(level_sel)))
          cat("\n")
          
          if(level_sel != "0")
          {
            
            TABLE_AT_LEAST_enhancer_Cell_Type_sel_0<-TABLE_AT_LEAST_enhancer_Cell_Type_sel[which(TABLE_AT_LEAST_enhancer_Cell_Type_sel[,indx.int_2] == "0"),]
            
            # cat("TABLE_AT_LEAST_enhancer_Cell_Type_sel_0\n")
            # cat(str(TABLE_AT_LEAST_enhancer_Cell_Type_sel_0))
            # cat("\n")
            
            TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel<-TABLE_AT_LEAST_enhancer_Cell_Type_sel[which(TABLE_AT_LEAST_enhancer_Cell_Type_sel[,indx.int_2] == level_sel),]
            
            # cat("TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel\n")
            # cat(str(TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel))
            # cat("\n")
            
            levels_LINEAGE_CLASSES<-levels(TABLE_AT_LEAST_enhancer[,indx.sel])
            
            list_LINEAGE_CLASSES<-list()
            
            for(iteration_levels_LINEAGE_CLASSES in 1:length(levels_LINEAGE_CLASSES))
            {
              
              
              levels_LINEAGE_CLASSES_sel<-levels_LINEAGE_CLASSES[iteration_levels_LINEAGE_CLASSES]
              
              cat("->\t")
              cat(sprintf(as.character(levels_LINEAGE_CLASSES_sel)))
              cat("\n")
              
              
              TABLE_AT_LEAST_enhancer_Cell_Type_sel_0_LINEAGE_CLASSES_sel<-TABLE_AT_LEAST_enhancer_Cell_Type_sel_0[which(TABLE_AT_LEAST_enhancer_Cell_Type_sel_0[,indx.sel] == levels_LINEAGE_CLASSES_sel),]
              
              # cat("TABLE_AT_LEAST_enhancer_Cell_Type_sel_0_LINEAGE_CLASSES_sel\n")
              # cat(str(TABLE_AT_LEAST_enhancer_Cell_Type_sel_0_LINEAGE_CLASSES_sel))
              # cat("\n")
              
              TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_LINEAGE_CLASSES_sel<-TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel[which(TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel[,indx.sel] == levels_LINEAGE_CLASSES_sel),]
              
              # cat("TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_LINEAGE_CLASSES_sel\n")
              # cat(str(TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_LINEAGE_CLASSES_sel))
              # cat("\n")
              
              
              TABLE_AT_LEAST_enhancer_Cell_Type_sel_0_NOT_LINEAGE_CLASSES_sel<-TABLE_AT_LEAST_enhancer_Cell_Type_sel_0[which(TABLE_AT_LEAST_enhancer_Cell_Type_sel_0[,indx.sel] != levels_LINEAGE_CLASSES_sel),]
              
              # cat("TABLE_AT_LEAST_enhancer_Cell_Type_sel_0_NOT_LINEAGE_CLASSES_sel\n")
              # cat(str(TABLE_AT_LEAST_enhancer_Cell_Type_sel_0_NOT_LINEAGE_CLASSES_sel))
              # cat("\n")
              
              TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_NOT_LINEAGE_CLASSES_sel<-TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel[which(TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel[,indx.sel] != levels_LINEAGE_CLASSES_sel),]
              
              # cat("TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_NOT_LINEAGE_CLASSES_sel\n")
              # cat(str(TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_NOT_LINEAGE_CLASSES_sel))
              # cat("\n")
              
              
              
              if(dim(TABLE_AT_LEAST_enhancer_Cell_Type_sel_0_LINEAGE_CLASSES_sel)[1]>0)
              {
                if(dim(TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_LINEAGE_CLASSES_sel)[1]>0)
                {
                  cat("->Hello_world_1\n")
                  
                  
                  A.df<-as.data.frame(rbind(cbind(sum(TABLE_AT_LEAST_enhancer_Cell_Type_sel_0_NOT_LINEAGE_CLASSES_sel$AT_LEAST_instances),sum(TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_NOT_LINEAGE_CLASSES_sel$AT_LEAST_instances)),
                                            cbind(TABLE_AT_LEAST_enhancer_Cell_Type_sel_0_LINEAGE_CLASSES_sel$AT_LEAST_instances,TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_LINEAGE_CLASSES_sel$AT_LEAST_instances)),
                                      stringsAsFactors=F)
                  
                  colnames(A.df)<-c("0",level_sel)
                  row.names(A.df)<-c(paste("NOT",levels_LINEAGE_CLASSES_sel,sep='_'),
                                     levels_LINEAGE_CLASSES_sel)
                  
                  
                  # cat("A.df\n")
                  # cat(str(A.df))
                  # cat("\n")
                  
                  tab.chisq.test<-chisq.test(A.df,correct = TRUE)
                  
                  # cat("tab.chisq.test\n")
                  # cat(str(tab.chisq.test))
                  # cat("\n")
                  
                  log_pval<-as.numeric(round(-1*log10(tab.chisq.test$p.value),2))
                  
                  # cat("log_pval\n")
                  # cat(str(log_pval))
                  # cat("\n")
                  
                  a.df<-as.data.frame(cbind(level_sel,levels_LINEAGE_CLASSES_sel,log_pval), stringsAsFactors=F)
                  
                  colnames(a.df)<-c("Category","COMPARISON_category","log_pval")
                  
                  # cat("a.df\n")
                  # cat(str(a.df))
                  # cat("\n")
                  
                  list_LINEAGE_CLASSES[[iteration_levels_LINEAGE_CLASSES]]<-a.df
                  
                  
                }else{
                  cat("->Hello_world_2\n")
                  
                  
                  A.df<-as.data.frame(rbind(cbind(sum(TABLE_AT_LEAST_enhancer_Cell_Type_sel_0_NOT_LINEAGE_CLASSES_sel$AT_LEAST_instances),sum(TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_NOT_LINEAGE_CLASSES_sel$AT_LEAST_instances)),
                                            cbind(TABLE_AT_LEAST_enhancer_Cell_Type_sel_0_LINEAGE_CLASSES_sel$AT_LEAST_instances,0)),
                                      stringsAsFactors=F)
                  
                  colnames(A.df)<-c("0",level_sel)
                  row.names(A.df)<-c(paste("NOT",levels_LINEAGE_CLASSES_sel,sep='_'),
                                     levels_LINEAGE_CLASSES_sel)
                  
                  
                  # cat("A.df\n")
                  # cat(str(A.df))
                  # cat("\n")
                  
                  tab.chisq.test<-chisq.test(A.df,correct = TRUE)
                  
                  # cat("tab.chisq.test\n")
                  # cat(str(tab.chisq.test))
                  # cat("\n")
                  
                  log_pval<-as.numeric(round(-1*log10(tab.chisq.test$p.value),2))
                  
                  # cat("log_pval\n")
                  # cat(str(log_pval))
                  # cat("\n")
                  
                  a.df<-as.data.frame(cbind(level_sel,levels_LINEAGE_CLASSES_sel,log_pval), stringsAsFactors=F)
                  
                  colnames(a.df)<-c("Category","COMPARISON_category","log_pval")
                  
                  # cat("a.df\n")
                  # cat(str(a.df))
                  # cat("\n")
                  
                  list_LINEAGE_CLASSES[[iteration_levels_LINEAGE_CLASSES]]<-a.df
                  
                  
                }# dim(TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_LINEAGE_CLASSES_sel)[1]>0
                
              }else{
                
                if(dim(TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_LINEAGE_CLASSES_sel)[1]>0)
                {
                  
                  cat("->Hello_world_3\n")
                  
                  A.df<-as.data.frame(rbind(cbind(sum(TABLE_AT_LEAST_enhancer_Cell_Type_sel_0_NOT_LINEAGE_CLASSES_sel$AT_LEAST_instances),sum(TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_NOT_LINEAGE_CLASSES_sel$AT_LEAST_instances)),
                                            cbind(0,TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_LINEAGE_CLASSES_sel$AT_LEAST_instances)),
                                      stringsAsFactors=F)
                  
                  colnames(A.df)<-c("0",level_sel)
                  row.names(A.df)<-c(paste("NOT",levels_LINEAGE_CLASSES_sel,sep='_'),
                                     levels_LINEAGE_CLASSES_sel)
                  
                  
                  # cat("A.df\n")
                  # cat(str(A.df))
                  # cat("\n")
                  
                  tab.chisq.test<-chisq.test(A.df,correct = TRUE)
                  
                  # cat("tab.chisq.test\n")
                  # cat(str(tab.chisq.test))
                  # cat("\n")
                  
                  log_pval<-as.numeric(round(-1*log10(tab.chisq.test$p.value),2))
                  
                  # cat("log_pval\n")
                  # cat(str(log_pval))
                  # cat("\n")
                  
                  a.df<-as.data.frame(cbind(level_sel,levels_LINEAGE_CLASSES_sel,log_pval), stringsAsFactors=F)
                  
                  colnames(a.df)<-c("Category","COMPARISON_category","log_pval")
                  
                  # cat("a.df\n")
                  # cat(str(a.df))
                  # cat("\n")
                  
                  list_LINEAGE_CLASSES[[iteration_levels_LINEAGE_CLASSES]]<-a.df
                }else{
                  
                  cat("->Hello_world_4\n")
                  
                  A.df<-as.data.frame(rbind(cbind(sum(TABLE_AT_LEAST_enhancer_Cell_Type_sel_0_NOT_LINEAGE_CLASSES_sel$AT_LEAST_instances),sum(TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_NOT_LINEAGE_CLASSES_sel$AT_LEAST_instances)),
                                            cbind(0,0)),
                                      stringsAsFactors=F)
                  
                  colnames(A.df)<-c("0",level_sel)
                  row.names(A.df)<-c(paste("NOT",levels_LINEAGE_CLASSES_sel,sep='_'),
                                     levels_LINEAGE_CLASSES_sel)
                  
                  
                  # cat("A.df\n")
                  # cat(str(A.df))
                  # cat("\n")
                  
                  tab.chisq.test<-chisq.test(A.df,correct = TRUE)
                  
                  # cat("tab.chisq.test\n")
                  # cat(str(tab.chisq.test))
                  # cat("\n")
                  
                  log_pval<-as.numeric(round(-1*log10(tab.chisq.test$p.value),2))
                  
                  # cat("log_pval\n")
                  # cat(str(log_pval))
                  # cat("\n")
                  
                  a.df<-as.data.frame(cbind(level_sel,levels_LINEAGE_CLASSES_sel,log_pval), stringsAsFactors=F)
                  
                  colnames(a.df)<-c("Category","COMPARISON_category","log_pval")
                  
                  # cat("a.df\n")
                  # cat(str(a.df))
                  # cat("\n")
                  
                  list_LINEAGE_CLASSES[[iteration_levels_LINEAGE_CLASSES]]<-a.df
                  
                }#dim(TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_LINEAGE_CLASSES_sel)[1]>0
                
                
                
              }# dim(TABLE_AT_LEAST_enhancer_Cell_Type_sel_0_LINEAGE_CLASSES_sel)[1]>0
              
              
            }# iteration_levels_LINEAGE_CLASSES
            
            CHI_Sq_LINEAGE_CLASSES = unique(as.data.frame(data.table::rbindlist(list_LINEAGE_CLASSES, fill = T)))
            
            # cat("CHI_Sq_LINEAGE_CLASSES\n")
            # cat(str(CHI_Sq_LINEAGE_CLASSES))
            # cat("\n")
            
            list_ChiSq_category[[k]]<-CHI_Sq_LINEAGE_CLASSES
            
            # ############################################################################################################
            # quit(status = 1) 
            
            
          }#level_sel != "0"
        }# k
      }# dim(TABLE_AT_LEAST_enhancer_Cell_Type_sel)[1] > 0
      
      CHI_Sq_CSQ_CATEGORY = unique(as.data.frame(data.table::rbindlist(list_ChiSq_category, fill = T)))
      
      CHI_Sq_CSQ_CATEGORY$Cell_Type<-Cell_Type_array_sel
      
      # cat("CHI_Sq_CSQ_CATEGORY\n")
      # cat(str(CHI_Sq_CSQ_CATEGORY))
      # cat("\n")
      # 
      
      
      list_ChiSq_DEF[[i]]<-CHI_Sq_CSQ_CATEGORY
      
      # ############################################################################################################
      # quit(status = 1) 
      
    }# i Cell_Type
    
    
    CHI_Sq_DEF = unique(as.data.frame(data.table::rbindlist(list_ChiSq_DEF, fill = T)))
    
    
    cat("CHI_Sq_DEF_0\n")
    cat(str(CHI_Sq_DEF))
    cat("\n")
    
    
    CHI_Sq_DEF$GROUP<-CLASSES_array_sel
    
    
    cat("CHI_Sq_DEF_1\n")
    cat(str(CHI_Sq_DEF))
    cat("\n")
    
    
    
    List_SUPER_DEF_STATS[[iteration_CLASSES_array]]<-CHI_Sq_DEF
    
    
  }#iteration_CLASSES_array
  
  
  
  
  CLASSIF_FOR_INTER = unique(as.data.frame(data.table::rbindlist(List_CLASSIF_FOR_INTER_ANALYSIS, fill = T)))
  
  cat("CLASSIF_FOR_INTER_0\n")
  cat(str(CLASSIF_FOR_INTER))
  cat("\n")
  
  
  CLASSIF_FOR_INTER.m<-melt(CLASSIF_FOR_INTER, id.vars=c("Cell_Type","AT_LEAST_instances","AT_LEAST_TOTAL","AT_LEAST_category","Perc"))
  
  # measure.vars=c("CLASS_erythroid","CLASS_mega","CLASS_gran_mono","CLASS_lymph","CLASS_ALL_wbc"),
  # value.name="Perc"
  
  
  CLASSIF_FOR_INTER.m<-CLASSIF_FOR_INTER.m[!is.na(CLASSIF_FOR_INTER.m$value),]
  
  CLASSIF_FOR_INTER.m<-droplevels(CLASSIF_FOR_INTER.m[-which(CLASSIF_FOR_INTER.m$variable == "CLASS_ALL_wbc"),])
  
  CLASSIF_FOR_INTER.m$value<-factor(CLASSIF_FOR_INTER.m$value,
                                    levels=c("OTHER","erythroid_lineage","mega_lineage","gran_mono_lineage","lymph_lineage"),
                                    ordered=T)
  
  
  
  cat("CLASSIF_FOR_INTER.m_0\n")
  cat(str(CLASSIF_FOR_INTER.m))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(CLASSIF_FOR_INTER.m$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(CLASSIF_FOR_INTER.m$variable)))))
  cat("\n")
  
  cat(sprintf(as.character(names(summary(as.factor(CLASSIF_FOR_INTER.m$value))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(CLASSIF_FOR_INTER.m$value)))))
  cat("\n")
  
  
  
  
  #### PANEL GRAPH ----    
  
  breaks.y<-seq(0,100,by=25)
  labels.y<-as.character(breaks.y)
  
  
  
  colors<-c("gray","firebrick1","brown","magenta","dodgerblue")
  
  
  
  cat("GRAPH_START\n")
  
  
  graph<-ggplot(data=CLASSIF_FOR_INTER.m,
                aes(x=AT_LEAST_category,
                    y=Perc,
                    fill=value)) +
    geom_bar(stat="identity",colour='black')+
    theme(axis.title.y=element_text(size=18, family="sans"),
          axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
          axis.text.x=element_text(angle=45,size=18, color="black", family="sans"))+
    scale_y_continuous(name=paste("Perentage"),breaks=breaks.y,labels=labels.y,
                       limits=c(breaks.y[1],breaks.y[length(breaks.y)]+1))+
    scale_x_discrete(name=NULL, drop=F)+
    theme_classic()+
    scale_fill_manual(values=colors, drop=F)+
    theme(legend.position="bottom")+
    facet_grid(Cell_Type ~ variable)+
    ggeasy::easy_center_title()
  
  cat("GRAPH_END\n")
  
  
  # setwd(out)
  # 
  # svgname<-paste(paste("TEST",".svg",sep=''))
  # makesvg = TRUE
  # 
  # if (makesvg == TRUE)
  # {
  #   ggsave(svgname, plot= graph,
  #          device="svg",
  #          height=10, width=12)
  # }
  # 
  # ########################################################################################
  # quit(status = 1)
  #   
  
  setwd(path5)
  
  svgname<-paste(paste("enhancer_single_graph",".svg",sep=''))
  makesvg = TRUE
  
  if (makesvg == TRUE)
  {
    ggsave(svgname, plot= graph,
           device="svg",
           height=10, width=12)
  }
  
  CHI_Sq_FINAL_DEF = unique(as.data.frame(data.table::rbindlist(List_SUPER_DEF_STATS, fill = T)))
  
  cat("CHI_Sq_FINAL_DEF_0\n")
  cat(str(CHI_Sq_FINAL_DEF))
  cat("\n")
  
  CHI_Sq_FINAL_DEF$COMPARISON_category<-factor(CHI_Sq_FINAL_DEF$COMPARISON_category,
                                               levels=c("OTHER","erythroid_lineage","mega_lineage","gran_mono_lineage","lymph_lineage",
                                                        "CLASS_ALL_wbc"),
                                               ordered=T)
  
  CHI_Sq_FINAL_DEF$Category<-factor(CHI_Sq_FINAL_DEF$Category,
                                    levels = c("0","AT LEAST 1","AT LEAST 2","AT LEAST 3"),
                                    ordered=T)
  
  CHI_Sq_FINAL_DEF$log_pval<-as.numeric(CHI_Sq_FINAL_DEF$log_pval)
  
  CHI_Sq_FINAL_DEF$Cell_Type<-factor(CHI_Sq_FINAL_DEF$Cell_Type,
                                     levels = c("K562","CHRF","HL60","THP1"),
                                     ordered=T)
  CHI_Sq_FINAL_DEF$GROUP<-factor(CHI_Sq_FINAL_DEF$GROUP,
                                 levels=c("CLASS_erythroid","CLASS_mega","CLASS_gran_mono","CLASS_lymph",
                                          "CLASS_ALL_wbc"),
                                 ordered=T)
  
  
  CHI_Sq_FINAL_DEF<-CHI_Sq_FINAL_DEF[order(CHI_Sq_FINAL_DEF$Cell_Type,CHI_Sq_FINAL_DEF$GROUP,CHI_Sq_FINAL_DEF$Category,CHI_Sq_FINAL_DEF$COMPARISON_category),]
  
  cat("CHI_Sq_FINAL_DEF_1\n")
  cat(str(CHI_Sq_FINAL_DEF))
  cat("\n")
  
  CHI_Sq_FINAL_DEF<-CHI_Sq_FINAL_DEF[-which(CHI_Sq_FINAL_DEF$COMPARISON_category == "OTHER"),]
  
  CHI_Sq_FINAL_DEF<-droplevels(CHI_Sq_FINAL_DEF)
  
  cat("CHI_Sq_FINAL_DEF_2\n")
  cat(str(CHI_Sq_FINAL_DEF))
  cat("\n")
  
  #### SAVE----
  
  
  setwd(path5)
  
  saveRDS(CLASSIF_FOR_INTER,file="CLASSIF_FOR_INTER_enhancer_CLASS_stats.rds")
  
  
  
  saveRDS(CHI_Sq_FINAL_DEF,file="Lineage_enhancer_CLASS_stats.rds")
  
  # saveRDS(TABLE_AT_LEAST_enhancer,file="Lineage_enhancer_df.rds")
  
}


print_join = function(option_list)
{
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
  
  #### READ DATA----
  
  
  path5<-paste(out,'LINEAGE_plots','/', sep='')
  
  cat("path5\n")
  cat(sprintf(as.character(path5)))
  cat("\n")
  
  
  if (file.exists(path5)){
    
    
    
    
  } else {
    dir.create(file.path(path5))
    
  }
  
  
  
  #### READ DATA ----
  
  setwd(path5)
  
 
  
  # cat("enhancer_CLASS_graphs\n")
  # cat(str(enhancer_CLASS_graphs))
  # cat("\n")
  
  enhancer_CLASS_stats<-readRDS(file="Lineage_enhancer_CLASS_stats.rds")
  enhancer_CLASS_stats$ER_Type<-"enhancer_CLASS"
  
  cat("enhancer_CLASS_stats\n")
  cat(str(enhancer_CLASS_stats))
  cat("\n")
  
  
  
  E_Plus_ASE_CLASS_stats<-readRDS(file="Lineage_E_Plus_ASE_CLASS_stats.rds")
  E_Plus_ASE_CLASS_stats$ER_Type<-"E_Plus_ASE_CLASS"
  
  cat("E_Plus_ASE_CLASS_stats\n")
  cat(str(E_Plus_ASE_CLASS_stats))
  cat("\n")
  
  Cell_Type_array<-c("K562","CHRF","HL60","THP1")
  
  cat("Cell_Type_array\n")
  cat(str(Cell_Type_array))
  cat("\n")
  
  
  DEF_stats<-rbind(enhancer_CLASS_stats,
                   E_Plus_ASE_CLASS_stats)
  
  cat("DEF_stats\n")
  cat(str(DEF_stats))
  cat("\n")
  
  setwd(path5)
  
  write.table(DEF_stats,file="stats_object_ER_Chi.tsv",sep="\t",quote=F,row.names = F)
  
  
  
  # 
  # graph_enhancer_CLASS_erythroid<-readRDS(file=paste("Lineage_enhancer_CLASS_graphs_","CLASS_erythroid",".rds",sep=''))
  # graph_E_Plus_ASE_CLASS_erythroid<-readRDS(file=paste("Lineage_E_Plus_ASE_CLASS_graphs_","CLASS_erythroid",".rds",sep=''))
  # 
  # graph_enhancer_CLASS_mega<-readRDS(file=paste("Lineage_enhancer_CLASS_graphs_","CLASS_mega",".rds",sep=''))
  # graph_E_Plus_ASE_CLASS_mega<-readRDS(file=paste("Lineage_E_Plus_ASE_CLASS_graphs_","CLASS_mega",".rds",sep=''))
  # 
  # graph_enhancer_CLASS_gran_mono<-readRDS(file=paste("Lineage_enhancer_CLASS_graphs_","CLASS_gran_mono",".rds",sep=''))
  # graph_E_Plus_ASE_CLASS_gran_mono<-readRDS(file=paste("Lineage_E_Plus_ASE_CLASS_graphs_","CLASS_gran_mono",".rds",sep=''))
  # 
  # graph_enhancer_CLASS_lymph<-readRDS(file=paste("Lineage_enhancer_CLASS_graphs_","CLASS_lymph",".rds",sep=''))
  # graph_E_Plus_ASE_CLASS_lymph<-readRDS(file=paste("Lineage_E_Plus_ASE_CLASS_graphs_","CLASS_lymph",".rds",sep=''))
  # 
  # graph_enhancer_CLASS_ALL_wbc<-readRDS(file=paste("Lineage_enhancer_CLASS_graphs_","CLASS_ALL_wbc",".rds",sep=''))
  # graph_E_Plus_ASE_CLASS_ALL_wbc<-readRDS(file=paste("Lineage_E_Plus_ASE_CLASS_graphs_","CLASS_ALL_wbc",".rds",sep=''))
  # 
  # 
  # 
  # for(i in 1:length(Cell_Type_array))
  # {
  # 
  #   Cell_Type_array_sel<-Cell_Type_array[i]
  # 
  #   cat("---------------------------------------------->\t")
  #   cat(sprintf(as.character(Cell_Type_array_sel)))
  #   cat("\n")
  # 
  #   graph_enhancer_CLASS_erythroid_sel<-graph_enhancer_CLASS_erythroid[[Cell_Type_array_sel]]
  # 
  #   graph_enhancer_CLASS_erythroid_sel<-graph_enhancer_CLASS_erythroid_sel+
  #     theme(axis.title.y=element_text(size=18, family="sans"),
  #           axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
  #           axis.text.x=element_text(angle=45,vjust=1,hjust=1, size=16, color="black", family="sans"))+
  #     theme(legend.position="hidden")
  # 
  #   graph_E_Plus_ASE_CLASS_erythroid_sel<-graph_E_Plus_ASE_CLASS_erythroid[[Cell_Type_array_sel]]
  # 
  #   graph_E_Plus_ASE_CLASS_erythroid_sel<-graph_E_Plus_ASE_CLASS_erythroid_sel+
  #     theme(axis.title.y=element_text(size=18, family="sans"),
  #           axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
  #           axis.text.x=element_text(angle=45,vjust=1,hjust=1, size=16, color="black", family="sans"))+
  #     theme(legend.position="hidden")
  # 
  #   # cat("graph_E_Plus_ASE_CLASS_erythroid_sel\n")
  #   # cat(str(graph_E_Plus_ASE_CLASS_erythroid_sel))
  #   # cat("\n")
  # 
  # 
  #   graph_DEF2<-plot_grid(graph_enhancer_CLASS_erythroid_sel,graph_E_Plus_ASE_CLASS_erythroid_sel,
  #                         nrow = 1,
  #                         ncol=2,
  #                         labels = c("ACTIVE enhancer TILES)","ACTIVE E+ASE TILES)"),
  #                         label_size = 12,
  #                         align = "h",
  #                         rel_widths=c(1,1))
  # 
  #   cat("graph_DEF2\n")
  # 
  # 
  # 
  #   title <- ggdraw() + draw_label(paste("ER","CLASS_erythroid",Cell_Type_array_sel,sep="  "))
  # 
  # 
  #   graph_DEF3<- plot_grid(title,graph_DEF2,
  #                          ncol=1,
  #                          rel_heights = c(0.05,1))
  # 
  #   cat("graph_DEF3\n")
  # 
  # 
  #   setwd(path5)
  # 
  #   svgname<-paste("ER","_","CLASS_erythroid","_",Cell_Type_array_sel,".svg", sep='')
  #   makesvg = TRUE
  # 
  #   if (makesvg == TRUE)
  #   {
  # 
  #     ggsave(svgname, plot= graph_DEF3,
  #            device="svg",
  #            height=10, width=12)
  #   }
  # 
  #   cat("svg_WITHOUT_LEGEND\n")
  # 
  # 
  # 
  #   graph_enhancer_CLASS_erythroid_sel<-graph_enhancer_CLASS_erythroid_sel+
  #     theme(axis.title.y=element_text(size=18, family="sans"),
  #           axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
  #           axis.text.x=element_text(angle=45,vjust=1,hjust=1, size=16, color="black", family="sans"))+
  #     theme(legend.position="right")
  # 
  # 
  #   graph_E_Plus_ASE_CLASS_erythroid_sel<-graph_E_Plus_ASE_CLASS_erythroid_sel+
  #     theme(axis.title.y=element_text(size=18, family="sans"),
  #           axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
  #           axis.text.x=element_text(angle=45,vjust=1,hjust=1, size=16, color="black", family="sans"))+
  #     theme(legend.position="right")
  # 
  #   # cat("graph_E_Plus_ASE_CLASS_erythroid_sel\n")
  #   # cat(str(graph_E_Plus_ASE_CLASS_erythroid_sel))
  #   # cat("\n")
  # 
  # 
  #   graph_DEF2<-plot_grid(graph_enhancer_CLASS_erythroid_sel,graph_E_Plus_ASE_CLASS_erythroid_sel,
  #                         nrow = 1,
  #                         ncol=2,
  #                         labels = c("ACTIVE enhancer TILES)","ACTIVE E+ASE TILES)"),
  #                         label_size = 12,
  #                         align = "h",
  #                         rel_widths=c(1,1))
  # 
  #   cat("graph_DEF2\n")
  # 
  # 
  # 
  #   title <- ggdraw() + draw_label(paste("ER","CLASS_erythroid",Cell_Type_array_sel,sep="  "))
  # 
  # 
  #   graph_DEF3<- plot_grid(title,graph_DEF2,
  #                          ncol=1,
  #                          rel_heights = c(0.05,1))
  # 
  #   cat("graph_DEF3\n")
  # 
  # 
  #   setwd(path5)
  # 
  #   svgname<-paste("ER","_","CLASS_erythroid","_",Cell_Type_array_sel,".svg", sep='')
  #   makesvg = TRUE
  # 
  #   if (makesvg == TRUE)
  #   {
  # 
  #     ggsave(svgname, plot= graph_DEF3,
  #            device="svg",
  #            height=10, width=12)
  #   }
  # 
  #   cat("svg_WITH_LEGEND\n")
  # 
  #   # #############################################################################################################################################
  #   # quit(status=1)
  # 
  # }#i
  
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
    make_option(c("--CSQ_colors"), type="character", default=NULL, 
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
  
  
  data_wrangling(opt)
  phenotype_Lineage_ER_stacked_barplots_E_Plus_ASE_CLASS(opt)
  phenotype_Lineage_ER_stacked_barplots_enhancer_CLASS(opt)
  print_join(opt)
  
}
  
  
  
 

###########################################################################

system.time( main() )
