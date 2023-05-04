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
suppressMessages(library("ggrepel", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
library("gtools", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")


library("reshape2", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")


opt = NULL



Parameter_recalculation = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  

  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
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
  
  #### READ and transform K562_replicates ----
  
  K562_replicates<-unlist(strsplit(opt$K562_replicates_QC_PASS, ","))
  
  cat("K562_replicates\n")
  cat(sprintf(as.character(K562_replicates)))
  cat("\n")
  
  K562_df<-as.data.frame(cbind(rep("K562",length(K562_replicates)),K562_replicates), stringsAsFactors=F)
  colnames(K562_df)<-c("Cell_Type","master_sample")
  
  
  cat("K562_df\n")
  str(K562_df)
  cat("\n")
  
  
  #### READ and transform CHRF_replicates ----
  
  CHRF_replicates<-unlist(strsplit(opt$CHRF_replicates_QC_PASS, ","))
  
  cat("CHRF_replicates\n")
  cat(sprintf(as.character(CHRF_replicates)))
  cat("\n")
  
  CHRF_df<-as.data.frame(cbind(rep("CHRF",length(CHRF_replicates)),CHRF_replicates), stringsAsFactors=F)
  colnames(CHRF_df)<-c("Cell_Type","master_sample")
  
  
  cat("CHRF_df\n")
  str(CHRF_df)
  cat("\n")
  
  #### READ and transform HL60_replicates ----
  
  HL60_replicates<-unlist(strsplit(opt$HL60_replicates_QC_PASS, ","))
  
  cat("HL60_replicates\n")
  cat(sprintf(as.character(HL60_replicates)))
  cat("\n")
  
  HL60_df<-as.data.frame(cbind(rep("HL60",length(HL60_replicates)),HL60_replicates), stringsAsFactors=F)
  
  colnames(HL60_df)<-c("Cell_Type","master_sample")
  
  cat("HL60_df\n")
  str(HL60_df)
  cat("\n")
  
  #### READ and transform THP1_replicates ----
  
  THP1_replicates<-unlist(strsplit(opt$THP1_replicates_QC_PASS, ","))
  
  cat("THP1_replicates\n")
  cat(sprintf(as.character(THP1_replicates)))
  cat("\n")
  
  #quit(status = 1)
  
  THP1_df<-as.data.frame(cbind(rep("THP1",length(THP1_replicates)),THP1_replicates), stringsAsFactors=F)
  
  colnames(THP1_df)<-c("Cell_Type","master_sample")
  cat("THP1_df\n")
  str(THP1_df)
  cat("\n")
  
  #### Replicates merge ----
  
  Replicates_merge<-rbind(K562_df,CHRF_df,HL60_df,THP1_df)
  
  Replicates_merge$Cell_Type<-factor(Replicates_merge$Cell_Type,
                                     levels=c("K562","CHRF","HL60","THP1"),
                                     ordered=T)
  
  cat("Replicates_merge\n")
  str(Replicates_merge)
  cat("\n")
  
  #### READ data ----
  
  
  
  
  
  setwd(out)
  
  filename=paste('Rosetta_extended','.rds',sep='')
  Rosetta_extended <-readRDS(file=filename)
  
  
  cat("Rosetta_extended_0\n")
  cat(str(Rosetta_extended))
  cat("\n")
  cat(sprintf(as.character(names(summary(Rosetta_extended$batch)))))
  cat("\n")
  cat(sprintf(as.character(summary(Rosetta_extended$batch))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Rosetta_extended$type)))))
  cat("\n")
  cat(sprintf(as.character(summary(Rosetta_extended$type))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Rosetta_extended$condition)))))
  cat("\n")
  cat(sprintf(as.character(summary(Rosetta_extended$condition))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Rosetta_extended$master_sample)))))
  cat("\n")
  cat(sprintf(as.character(summary(Rosetta_extended$master_sample))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Rosetta_extended$bc)))))
  cat("\n")
  cat(sprintf(as.character(summary(Rosetta_extended$bc))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Rosetta_extended$Cell_Type)))))
  cat("\n")
  cat(sprintf(as.character(summary(Rosetta_extended$Cell_Type))))
  cat("\n")
  
  Deconvolve_table<-unique(Rosetta_extended[,c(which(colnames(Rosetta_extended) == "Cell_Type"),
                                               which(colnames(Rosetta_extended) == "master_sample"),
                                               which(colnames(Rosetta_extended) == "batch"))])
  
  cat("Deconvolve_table_0\n")
  cat(str(Deconvolve_table))
  cat("\n")
  
  
  filename=paste('Rosetta_df','.rds',sep='')
  Rosetta <-readRDS(file=filename)
  
  
  cat("Rosetta_0\n")
  cat(str(Rosetta))
  cat("\n")
  
  rm(Rosetta_extended)
  
  
  filename=paste('matrix_gDNA_NOR','.rds',sep='')
  gDNA_NOR <-readRDS(file=filename)
  
  cat("gDNA_NOR\n")
  cat(str(gDNA_NOR))
  cat("\n")
  
  # matrix_gDNA_NOR<-as.matrix(gDNA_NOR[,-which(colnames(gDNA_NOR) == "REAL_TILE_Plus_carried_variants")])
  # 
  # row.names(matrix_gDNA_NOR)<-gDNA_NOR$REAL_TILE_Plus_carried_variants
  # 
  # cat("matrix_gDNA_NOR\n")
  # cat(str(matrix_gDNA_NOR))
  # cat("\n")
  # 
  # rm(gDNA_NOR)
  
  filename=paste('matrix_cDNA_NOR','.rds',sep='')
  #cDNA_NOR <-as.data.frame(fread(file=filename, sep="\t",header=T), stringsAsFactors = F)
  cDNA_NOR <-readRDS(file=filename)
  
  
  cat("cDNA_NOR\n")
  cat(str(cDNA_NOR))
  cat("\n")
  
  # matrix_cDNA_NOR<-as.matrix(cDNA_NOR[,-which(colnames(cDNA_NOR) == "REAL_TILE_Plus_carried_variants")])
  # 
  # row.names(matrix_cDNA_NOR)<-cDNA_NOR$REAL_TILE_Plus_carried_variants
  # 
  # cat("matrix_cDNA_NOR\n")
  # cat(str(matrix_cDNA_NOR))
  # cat("\n")
  # 
  # rm(cDNA_NOR)
  
  
  filename_1<-paste("Annot_gDNA_df",".rds", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  Annot_gDNA_df<-readRDS(file=filename_1)#,sep="\t")#,header=T), stringsAsFactors=F)
  
  cat("Annot_gDNA_df\n")
  cat(str(Annot_gDNA_df))
  cat("\n")
  
  
  #### Important all the Annot columns should be factors
  
  REF_bc<-paste("bc",seq(1,15,by=1),sep='')
  ALT_bc<-paste("bc",seq(16,30,by=1),sep='')
  
  Annot_gDNA_df$batch<-factor(Annot_gDNA_df$batch,
                              levels=levels(Annot_gDNA_df$batch))
  
  Annot_gDNA_df$bc<-factor(Annot_gDNA_df$bc,
                           levels=c(REF_bc,ALT_bc))
  
  Annot_gDNA_df$condition<-factor(Annot_gDNA_df$condition,
                                  levels=c("REF","ALT"))
  
  
  
  filename_1<-paste("Annot_cDNA_df",".rds", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  Annot_cDNA_df<-readRDS(file=filename_1)#,sep="\t")#,header=T), stringsAsFactors=F)
  
  cat("Annot_cDNA_df\n")
  cat(str(Annot_cDNA_df))
  cat("\n")
  
  
  #### Important all the Annot columns should be factors
  
  REF_bc<-paste("bc",seq(1,15,by=1),sep='')
  ALT_bc<-paste("bc",seq(16,30,by=1),sep='')
  
  Annot_cDNA_df$batch<-factor(Annot_cDNA_df$batch,
                              levels=levels(Annot_cDNA_df$batch))
  
  Annot_cDNA_df$bc<-factor(Annot_cDNA_df$bc,
                           levels=c(REF_bc,ALT_bc))
  
  Annot_cDNA_df$condition<-factor(Annot_cDNA_df$condition,
                                  levels=c("REF","ALT"))
  
  
  
  
  #### keep only accepted replicates ----
  
  Annot_gDNA_df_QC_PASS<-Annot_gDNA_df[which(Annot_gDNA_df$master_sample%in%Replicates_merge$master_sample),]
  Annot_gDNA_df_QC_PASS$accepted_colnames<-paste(Annot_gDNA_df_QC_PASS$batch,Annot_gDNA_df_QC_PASS$bc,sep="_")
  
  Annot_gDNA_df_QC_PASS<-Annot_gDNA_df_QC_PASS[order(Annot_gDNA_df_QC_PASS$master_sample),]
  
  cat("Annot_gDNA_df_QC_PASS_0\n")
  cat(str(Annot_gDNA_df_QC_PASS))
  cat("\n")
  
  Annot_gDNA_df_QC_PASS<-merge(Annot_gDNA_df_QC_PASS,
                               Deconvolve_table,
                               by=c("batch","master_sample"),
                               all.x=T)
  cat("Annot_gDNA_df_QC_PASS_1\n")
  cat(str(Annot_gDNA_df_QC_PASS))
  cat("\n")
  
  
  
  
  
  
  Annot_cDNA_df_QC_PASS<-Annot_cDNA_df[which(Annot_cDNA_df$master_sample%in%Replicates_merge$master_sample),]
  Annot_cDNA_df_QC_PASS$accepted_colnames<-paste(Annot_cDNA_df_QC_PASS$batch,Annot_cDNA_df_QC_PASS$bc,sep="_")
  
  Annot_cDNA_df_QC_PASS<-Annot_cDNA_df_QC_PASS[order(Annot_cDNA_df_QC_PASS$master_sample),]
  
  
  cat("Annot_cDNA_df_QC_PASS_0\n")
  cat(str(Annot_cDNA_df_QC_PASS))
  cat("\n")
  
  Annot_cDNA_df_QC_PASS<-merge(Annot_cDNA_df_QC_PASS,
                               Deconvolve_table,
                               by=c("batch","master_sample"),
                               all.x=T)
  cat("Annot_cDNA_df_QC_PASS_1\n")
  cat(str(Annot_cDNA_df_QC_PASS))
  cat("\n")
  
  
  
  
  
  #### QC_PASS df ----
  
  indx.gDNA<-c(which(colnames(gDNA_NOR) == "REAL_TILE_Plus_carried_variants"),
               which(colnames(gDNA_NOR)%in%Annot_gDNA_df_QC_PASS$accepted_colnames))
  
  cat("indx.gDNA\n")
  cat(str(indx.gDNA))
  cat("\n")
  
  gDNA_NOR_QC_PASS<-gDNA_NOR[,indx.gDNA]
  
  cat("gDNA_NOR_QC_PASS\n")
  cat(str(gDNA_NOR_QC_PASS))
  cat("\n")
  
    
  
  indx.cDNA<-c(which(colnames(cDNA_NOR) == "REAL_TILE_Plus_carried_variants"),
               which(colnames(cDNA_NOR)%in%Annot_cDNA_df_QC_PASS$accepted_colnames))
  
  cat("indx.cDNA\n")
  cat(str(indx.cDNA))
  cat("\n")
  
  cDNA_NOR_QC_PASS<-cDNA_NOR[,indx.cDNA]
  
  cat("cDNA_NOR_QC_PASS\n")
  cat(str(cDNA_NOR_QC_PASS))
  cat("\n")
  
 
  
  ##### Annot_gDNA_df_QC_PASS and Annot_cDNA_df_QC_PASS MERGE and ordered levels ----
  
  
  
  
  Annot_gDNA_df_QC_PASS$Cell_Type<-gsub("_.+$","",Annot_gDNA_df_QC_PASS$master_sample)
  Annot_gDNA_df_QC_PASS$Replicate<-gsub("^[^_]+_","",Annot_gDNA_df_QC_PASS$master_sample)
  Annot_gDNA_df_QC_PASS$type<-"gDNA"
  
  
  Annot_cDNA_df_QC_PASS$Cell_Type<-gsub("_.+$","",Annot_cDNA_df_QC_PASS$master_sample)
  Annot_cDNA_df_QC_PASS$Replicate<-gsub("^[^_]+_","",Annot_cDNA_df_QC_PASS$master_sample)
  Annot_cDNA_df_QC_PASS$type<-"cDNA"
  
  Annot_DEF<-rbind(Annot_gDNA_df_QC_PASS,Annot_cDNA_df_QC_PASS)
  
  
  Annot_DEF$type<-factor(Annot_DEF$type,
                         c("gDNA","cDNA"),
                         ordered=T)
  Annot_DEF$Cell_Type<-factor(Annot_DEF$Cell_Type,
                              c("K562","CHRF","HL60","THP1"),
                              ordered=T)
  # Annot_DEF$Replicate<-factor(Annot_DEF$Replicate,
  #                                 c("Rep1","Rep2","Rep3","Rep4","Rep5","Rep6","Rep7","Rep8","Rep9"),
  #                                 ordered=T)
  
  
  Annot_DEF<-Annot_DEF[order(Annot_DEF$Cell_Type,Annot_DEF$master_sample,Annot_DEF$type),]
  
  
  
  levels_master_samples<-levels(Annot_DEF$master_sample)
  # 
  # Annot_DEF$master_sample<-factor(Annot_DEF$master_sample,
  #                          levels=levels_master_samples,
  #                          ordered=T)
  # 
  # levels_batchs<-unique(as.character(Annot_DEF$batch))
  # 
  # Annot_DEF$batch<-factor(Annot_DEF$batch,
  #                                 levels=levels_batchs,
  #                                 ordered=T)
  
  cat("Annot_DEF_2\n")
  cat(str(Annot_DEF))
  cat("\n")
  
  Annot_DEF_subset<-unique(Annot_DEF[,c(which(colnames(Annot_DEF) == "Cell_Type"),
                                        which(colnames(Annot_DEF) == "master_sample"),
                                        which(colnames(Annot_DEF) == "batch"),
                                        which(colnames(Annot_DEF) == "type"),
                                        which(colnames(Annot_DEF) == "Replicate"))])
  
  cat("Annot_DEF_subset\n")
  cat(str(Annot_DEF_subset))
  cat("\n")
  
  levels_batchs<-levels(Annot_DEF_subset$batch)
  
  cat("levels_batchs\n")
  cat(str(levels_batchs))
  cat("\n")
  
  #### melt of df ----
  
  gDNA_NOR_QC_PASS.m<-melt(gDNA_NOR_QC_PASS, id.variables="REAL_TILE_Plus_carried_variants")
  
  cat("gDNA_NOR_QC_PASS.m_0\n")
  cat(str(gDNA_NOR_QC_PASS.m))
  cat("\n")
  
  #gDNA_NOR_QC_PASS.m$variable<-gsub("^aggregated_counts_","",gDNA_NOR_QC_PASS.m$variable)
  
  
  
  
  
  gDNA_NOR_QC_PASS.m$batch<-gsub("_.+$","",gDNA_NOR_QC_PASS.m$variable)
  gDNA_NOR_QC_PASS.m$bc<-gsub("^[^_]+_","",gDNA_NOR_QC_PASS.m$variable)
  
  
  gDNA_NOR_QC_PASS.m$type<-"gDNA"
  
  
  
  gDNA_NOR_QC_PASS.m<-gDNA_NOR_QC_PASS.m[,-which(colnames(gDNA_NOR_QC_PASS.m) == "variable")]
  
  cat("gDNA_NOR_QC_PASS.m_1\n")
  cat(str(gDNA_NOR_QC_PASS.m))
  cat("\n")
  
  cDNA_NOR_QC_PASS.m<-melt(cDNA_NOR_QC_PASS, id.variables="REAL_TILE_Plus_carried_variants")
  
  cat("cDNA_NOR_QC_PASS.m_0\n")
  cat(str(cDNA_NOR_QC_PASS.m))
  cat("\n")
  
  #cDNA_NOR_QC_PASS.m$variable<-gsub("^aggregated_counts_","",cDNA_NOR_QC_PASS.m$variable)
  
  
  cDNA_NOR_QC_PASS.m$batch<-gsub("_.+$","",cDNA_NOR_QC_PASS.m$variable)
  cDNA_NOR_QC_PASS.m$bc<-gsub("^[^_]+_","",cDNA_NOR_QC_PASS.m$variable)
  
  
  cDNA_NOR_QC_PASS.m$type<-"cDNA"
  
  
  
  cDNA_NOR_QC_PASS.m<-cDNA_NOR_QC_PASS.m[,-which(colnames(cDNA_NOR_QC_PASS.m) == "variable")]
  
  cat("cDNA_NOR_QC_PASS.m_1\n")
  cat(str(cDNA_NOR_QC_PASS.m))
  cat("\n")
  
  
  
  
  #### rbind of matrixes ----
  
  
  REP_df<-rbind(gDNA_NOR_QC_PASS.m,cDNA_NOR_QC_PASS.m)
  
  
  REP_df$type<-factor(REP_df$type,
                      levels=c("gDNA","cDNA"),
                      ordered=T)
  

  
  
  REP_df$batch<-factor(REP_df$batch,
                       levels=levels_batchs,
                       ordered=T)
  
  cat("REP_df_1\n")
  cat(str(REP_df))
  cat("\n")
  
  #### merge REP with Annot_DEF ----
  
  
  REP_df<-merge(REP_df,
                Annot_DEF,
                by=c("batch","bc","type"),
                all.x=T)
  
  cat("REP_df_2\n")
  cat(str(REP_df))
  cat("\n")
  
  REP_df<-REP_df[order(REP_df$master_sample,REP_df$type),]
  
  
  REP_df$sample<-paste(REP_df$master_sample,REP_df$type,sep="_")
  
  
  levels_sample<-unique(as.character(REP_df$sample))
  
  cat("levels_sample\n")
  cat(str(levels_sample))
  cat("\n")
  
  REP_df$sample<-factor(REP_df$sample,
                        levels=levels_sample,
                        ordered=T)
  
  cat("REP_df_3\n")
  cat(str(REP_df))
  cat("\n")
  
  
 
  
  #### LISTS DEFINITION ----
  
  List_values<-list()
  List_medians<-list()
  
  #### POOL for cDNA and gDNA calculation -----
  
  
  REP_df.dt<-data.table(REP_df,
                        key=c("REAL_TILE_Plus_carried_variants","sample"))
  
  
  cat("----------------------------------------------------------------------------------------------------------------->REP_df.dt\n")
  cat(str(REP_df.dt))
  cat("\n")
  
  
  dna_rna_pool_df<-unique(as.data.frame(REP_df.dt[,.(Median_value=median(value),
                                                     type=type,
                                                     Cell_Type=Cell_Type,
                                                     batch=batch,
                                                     master_sample=master_sample,
                                                     Replicate=Replicate),by=key(REP_df.dt)],stringsAsFactors=F))
  
  
  cat("dna_rna_pool_df_1\n")
  cat(str(dna_rna_pool_df))
  cat("\n")
  
  dna_rna_pool_df_NOR<-dna_rna_pool_df[,-which(colnames(dna_rna_pool_df) == "sample")]
  
  
  # cat("dna_rna_pool_df_NOR_1\n")
  # cat(str(dna_rna_pool_df_NOR))
  # cat("\n")
  
  dna_rna_pool_df_NOR<-droplevels(dna_rna_pool_df_NOR)
  
  # cat("dna_rna_pool_df_NOR_2\n")
  # cat(str(dna_rna_pool_df_NOR))
  # cat("\n")
  
  
  dna_pool_df_NOR<-dna_rna_pool_df_NOR[which(dna_rna_pool_df_NOR$type == "gDNA"),]
  
  # cat("dna_pool_df_NOR_1\n")
  # cat(str(dna_pool_df_NOR))
  # cat("\n")
  
  dna_pool_df_NOR<-droplevels(dna_pool_df_NOR)
  
  cat("dna_pool_df_NOR_2\n")
  cat(str(dna_pool_df_NOR))
  cat("\n")
  
  List_values[['dna_pool']]<-dna_pool_df_NOR
  
  
  dna_pool_df_NOR.dt<-data.table(dna_pool_df_NOR,
                                 key=c("REAL_TILE_Plus_carried_variants","Cell_Type"))
  
  dna_pool_median<-unique(as.data.frame(dna_pool_df_NOR.dt[,.(median_Median_value=median(Median_value)),
                                                           by=key(dna_pool_df_NOR.dt)],stringsAsFactors=F))
  
  
  cat("dna_pool_median_0\n")
  cat(str(dna_pool_median))
  cat("\n")
  
  List_medians[['dna_pool']]<-dna_pool_median
  
  
  dna_pool_df_NOR_wide<-pivot_wider(dna_pool_df_NOR,
                                    id_cols=REAL_TILE_Plus_carried_variants,
                                    names_from=master_sample,
                                    values_from=Median_value)
  
  cat("dna_pool_df_NOR_wide\n")
  cat(str(dna_pool_df_NOR_wide))
  cat("\n")
  
  dna_pool_df_NOR_wide.matrix<-as.matrix(dna_pool_df_NOR_wide[,-1])
  
  row.names(dna_pool_df_NOR_wide.matrix)<-dna_pool_df_NOR_wide$REAL_TILE_Plus_carried_variants
  
  cat("dna_pool_df_NOR_wide.matrix\n")
  cat(str(dna_pool_df_NOR_wide.matrix))
  cat("\n")
  
  
  rna_pool_df_NOR<-dna_rna_pool_df_NOR[which(dna_rna_pool_df_NOR$type == "cDNA"),]
  
  # cat("rna_pool_df_NOR_1\n")
  # cat(str(rna_pool_df_NOR))
  # cat("\n")
  
  rna_pool_df_NOR<-droplevels(rna_pool_df_NOR)
  
  cat("rna_pool_df_NOR_2\n")
  cat(str(rna_pool_df_NOR))
  cat("\n")
  
  List_values[['rna_pool']]<-rna_pool_df_NOR
  
  rna_pool_df_NOR.dt<-data.table(rna_pool_df_NOR,
                                 key=c("REAL_TILE_Plus_carried_variants","Cell_Type"))
  
  rna_pool_median<-unique(as.data.frame(rna_pool_df_NOR.dt[,.(median_Median_value=median(Median_value)),
                                                           by=key(rna_pool_df_NOR.dt)],stringsAsFactors=F))
  
  
  cat("rna_pool_median_0\n")
  cat(str(rna_pool_median))
  cat("\n")
  
  List_medians[['rna_pool']]<-rna_pool_median
  
  
  rna_pool_df_NOR_wide<-pivot_wider(rna_pool_df_NOR,
                                    id_cols=REAL_TILE_Plus_carried_variants,
                                    names_from=master_sample,
                                    values_from=Median_value)
  
  cat("rna_pool_df_NOR_wide\n")
  cat(str(rna_pool_df_NOR_wide))
  cat("\n")
  
  rna_pool_df_NOR_wide.matrix<-as.matrix(rna_pool_df_NOR_wide[,-1])
  
  row.names(rna_pool_df_NOR_wide.matrix)<-rna_pool_df_NOR_wide$REAL_TILE_Plus_carried_variants
  
  cat("rna_pool_df_NOR_wide.matrix\n")
  cat(str(rna_pool_df_NOR_wide.matrix))
  cat("\n")
  
  #### LogFC ----
  
  
  LogFC_matrix<- foldchange2logratio(foldchange(rna_pool_df_NOR_wide.matrix+0.0000001,dna_pool_df_NOR_wide.matrix+0.0000001))
  
  cat("LogFC_matrix\n")
  cat(str(LogFC_matrix))
  cat("\n")
  
  LogFC_df<-as.data.frame(LogFC_matrix, stringsAsFactors=F)
  
  
  LogFC_df$REAL_TILE_Plus_carried_variants<-row.names(LogFC_df)
  
  cat("LogFC_df\n")
  cat(str(LogFC_df))
  cat("\n")
  
  LogFC_df.m<-melt(LogFC_df, id.cols="REAL_TILE_Plus_carried_variants")
  
  cat("LogFC_df.m_1\n")
  cat(str(LogFC_df.m))
  cat("\n")
  
  colnames(LogFC_df.m)[which(colnames(LogFC_df.m) == "variable")]<-"master_sample"
  colnames(LogFC_df.m)[which(colnames(LogFC_df.m) == "value")]<-"LogFC"
  
  
  LogFC_df.m$master_sample<-factor(LogFC_df.m$master_sample,
                                   levels=levels_master_samples,ordered=T)
  
  
  cat("LogFC_df.m_2\n")
  cat(str(LogFC_df.m))
  cat("\n")
  
  LogFC_df.m<-merge(LogFC_df.m,
                    Annot_DEF_subset,
                    by="master_sample")
  
  
  cat("LogFC_df.m_3\n")
  cat(str(LogFC_df.m))
  cat("\n")
  
  
  ### NO NA
  
  LogFC_df.m_NO_NA<-LogFC_df.m[!is.na(LogFC_df.m$LogFC),]
  
  
  cat("LogFC_df.m_NO_NA_1\n")
  cat(str(LogFC_df.m_NO_NA))
  cat("\n")
  
  ### Infinite reverted
  
  LogFC_df.m_NO_NA$LogFC[(is.infinite(LogFC_df.m_NO_NA$LogFC) &  LogFC_df.m_NO_NA$LogFC > 0)]<-60
  # LogFC_df.m_NO_NA<-LogFC_df.m_NO_NA[-(is.infinite(LogFC_df.m_NO_NA$LogFC) &  LogFC_df.m_NO_NA$LogFC < 0),]
  #LogFC_df.m_NO_NA$LogFC[(is.infinite(LogFC_df.m_NO_NA$LogFC) &  LogFC_df.m_NO_NA$LogFC < 0)]<--1*10000
  
 
  
  cat("LogFC_df.m_NO_NA_2_infinite_reverted\n")
  cat(str(LogFC_df.m_NO_NA))
  cat("\n")
  
  
  
  
  
  A<-summary(LogFC_df.m_NO_NA$LogFC)
  
  cat("summary_LogFC\n")
  cat(sprintf(as.character(names(A))))
  cat("\n")
  cat(sprintf(as.character(A)))
  cat("\n")
  
  
  check<-LogFC_df.m_NO_NA[which(LogFC_df.m_NO_NA$LogFC >0),]
  
  cat("check_LogFC\n")
  cat(str(check))
  cat("\n")
  cat(sprintf(as.character(names(summary(check$master_sample)))))
  cat("\n")
  cat(sprintf(as.character(summary(check$master_sample))))
  cat("\n")
  
  A<-summary(check$LogFC)
  
  cat("summary_LogFC\n")
  cat(sprintf(as.character(names(A))))
  cat("\n")
  cat(sprintf(as.character(A)))
  cat("\n")
  
  
  check<-LogFC_df.m_NO_NA[is.na(LogFC_df.m_NO_NA$LogFC),]
  
  cat("check_2\n")
  cat(str(check))
  cat("\n")
  
  
  List_values[['LogFC']]<-LogFC_df.m_NO_NA
  
  LogFC_df.m_NO_NA.dt<-data.table(LogFC_df.m_NO_NA,
                                  key=c("REAL_TILE_Plus_carried_variants","Cell_Type"))
  
  LogFC_median<-unique(as.data.frame(LogFC_df.m_NO_NA.dt[,.(median_LogFC=median(LogFC)),
                                                         by=key(LogFC_df.m_NO_NA.dt)],stringsAsFactors=F))
  
  
  cat("LogFC_median_0\n")
  cat(str(LogFC_median))
  cat("\n")
  
  List_medians[['LogFC']]<-LogFC_median
  
  
  check_LogFC<-LogFC_median[which(LogFC_median$median_LogFC >0),]
  
  cat("check_LogFC_LogFC\n")
  cat(str(check_LogFC))
  cat("\n")
  cat(sprintf(as.character(names(summary(check_LogFC$Cell_Type)))))
  cat("\n")
  cat(sprintf(as.character(summary(check_LogFC$Cell_Type))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(check_LogFC$REAL_TILE_Plus_carried_variants))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(check_LogFC$REAL_TILE_Plus_carried_variants)))))
  cat("\n")
  
  A<-summary(check_LogFC$median_LogFC)
  
  cat("summary_median_LogFC\n")
  cat(sprintf(as.character(names(A))))
  cat("\n")
  cat(sprintf(as.character(A)))
  cat("\n")
  
  
  #### REF and ALT dna and rna ----
  
  REP_df.dt<-data.table(REP_df,
                        key=c("REAL_TILE_Plus_carried_variants","sample","condition"))
  
  
  cat("-------------------------------------------------------------------------------------------------------------------------->condition\n")
  cat(str(REP_df.dt))
  cat("\n")
  
  
  dna_rna_condition_df<-unique(as.data.frame(REP_df.dt[,.(Median_value=median(value),
                                                          type=type,
                                                          Cell_Type=Cell_Type,
                                                          batch=batch,
                                                          master_sample=master_sample,
                                                          Replicate=Replicate),by=key(REP_df.dt)],stringsAsFactors=F))
  
  
  cat("dna_rna_condition_df_1\n")
  cat(str(dna_rna_condition_df))
  cat("\n")
  
  ### REF
  
  
  dna_rna_condition_df_NOR_REF<-unique(dna_rna_condition_df[which(dna_rna_condition_df$condition == "REF"),-c(which(colnames(dna_rna_condition_df) == "sample"),
                                                                                                                which(colnames(dna_rna_condition_df) == "condition"))])
  
  
  # cat("dna_rna_condition_df_NOR_REF_1\n")
  # cat(str(dna_rna_condition_df_NOR_REF))
  # cat("\n")
  
  dna_rna_condition_df_NOR_REF<-droplevels(dna_rna_condition_df_NOR_REF)
  
  cat("dna_rna_condition_df_NOR_REF_2\n")
  cat(str(dna_rna_condition_df_NOR_REF))
  cat("\n")
  
  
  
  dna_condition_df_NOR_REF<-dna_rna_condition_df_NOR_REF[which(dna_rna_condition_df_NOR_REF$type == "gDNA"),]
  
  # cat("dna_condition_df_NOR_REF_1\n")
  # cat(str(dna_condition_df_NOR_REF))
  # cat("\n")
  
  dna_condition_df_NOR_REF<-droplevels(dna_condition_df_NOR_REF)
  
  cat("dna_condition_df_NOR_REF_2\n")
  cat(str(dna_condition_df_NOR_REF))
  cat("\n")
  
  List_values[['dna_REF']]<-dna_condition_df_NOR_REF
  
  dna_condition_df_NOR_REF.dt<-data.table(dna_condition_df_NOR_REF,
                                          key=c("REAL_TILE_Plus_carried_variants","Cell_Type"))
  
  dna_REF_median<-unique(as.data.frame(dna_condition_df_NOR_REF.dt[,.(median_Median_value=median(Median_value)),
                                                                   by=key(dna_condition_df_NOR_REF.dt)],stringsAsFactors=F))
  
  
  cat("dna_REF_median_0\n")
  cat(str(dna_REF_median))
  cat("\n")
  
  List_medians[['dna_REF']]<-dna_REF_median
  
  
  dna_condition_df_NOR_REF_wide<-pivot_wider(dna_condition_df_NOR_REF,
                                             id_cols=REAL_TILE_Plus_carried_variants,
                                             names_from=master_sample,
                                             values_from=Median_value)
  
  cat("dna_condition_df_NOR_REF_wide\n")
  cat(str(dna_condition_df_NOR_REF_wide))
  cat("\n")
  
  dna_condition_df_NOR_REF_wide.matrix<-as.matrix(dna_condition_df_NOR_REF_wide[,-1])
  
  row.names(dna_condition_df_NOR_REF_wide.matrix)<-dna_condition_df_NOR_REF_wide$REAL_TILE_Plus_carried_variants
  
  cat("dna_condition_df_NOR_REF_wide.matrix\n")
  cat(str(dna_condition_df_NOR_REF_wide.matrix))
  cat("\n")
  
  
  
  
  rna_condition_df_NOR_REF<-dna_rna_condition_df_NOR_REF[which(dna_rna_condition_df_NOR_REF$type == "cDNA"),]
  
  # cat("rna_condition_df_NOR_REF_1\n")
  # cat(str(rna_condition_df_NOR_REF))
  # cat("\n")
  
  rna_condition_df_NOR_REF<-droplevels(rna_condition_df_NOR_REF)
  
  cat("rna_condition_df_NOR_REF_2\n")
  cat(str(rna_condition_df_NOR_REF))
  cat("\n")
  
  List_values[['rna_REF']]<-rna_condition_df_NOR_REF
  
  rna_condition_df_NOR_REF.dt<-data.table(rna_condition_df_NOR_REF,
                                          key=c("REAL_TILE_Plus_carried_variants","Cell_Type"))
  
  rna_REF_median<-unique(as.data.frame(rna_condition_df_NOR_REF.dt[,.(median_Median_value=median(Median_value)),
                                                                   by=key(rna_condition_df_NOR_REF.dt)],stringsAsFactors=F))
  
  
  cat("rna_REF_median_0\n")
  cat(str(rna_REF_median))
  cat("\n")
  
  List_medians[['rna_REF']]<-rna_REF_median
  
  
  rna_condition_df_NOR_REF_wide<-pivot_wider(rna_condition_df_NOR_REF,
                                             id_cols=REAL_TILE_Plus_carried_variants,
                                             names_from=master_sample,
                                             values_from=Median_value)
  
  cat("rna_condition_df_NOR_REF_wide\n")
  cat(str(rna_condition_df_NOR_REF_wide))
  cat("\n")
  
  rna_condition_df_NOR_REF_wide.matrix<-as.matrix(rna_condition_df_NOR_REF_wide[,-1])
  
  row.names(rna_condition_df_NOR_REF_wide.matrix)<-rna_condition_df_NOR_REF_wide$REAL_TILE_Plus_carried_variants
  
  cat("rna_condition_df_NOR_REF_wide.matrix\n")
  cat(str(rna_condition_df_NOR_REF_wide.matrix))
  cat("\n")
  
  
  
  ### ALT
  
  
  dna_rna_condition_df_NOR_ALT<-unique(dna_rna_condition_df[which(dna_rna_condition_df$condition == "ALT"),-c(which(colnames(dna_rna_condition_df) == "sample"),
                                                                                                                which(colnames(dna_rna_condition_df) == "condition"))])
  
  
  # cat("dna_rna_condition_df_NOR_ALT_1\n")
  # cat(str(dna_rna_condition_df_NOR_ALT))
  # cat("\n")
  
  dna_rna_condition_df_NOR_ALT<-droplevels(dna_rna_condition_df_NOR_ALT)
  
  cat("dna_rna_condition_df_NOR_ALT_2\n")
  cat(str(dna_rna_condition_df_NOR_ALT))
  cat("\n")
  
  
  
  dna_condition_df_NOR_ALT<-dna_rna_condition_df_NOR_ALT[which(dna_rna_condition_df_NOR_ALT$type == "gDNA"),]
  
  # cat("dna_condition_df_NOR_ALT_1\n")
  # cat(str(dna_condition_df_NOR_ALT))
  # cat("\n")
  
  dna_condition_df_NOR_ALT<-droplevels(dna_condition_df_NOR_ALT)
  
  cat("dna_condition_df_NOR_ALT_2\n")
  cat(str(dna_condition_df_NOR_ALT))
  cat("\n")
  
  List_values[['dna_ALT']]<-dna_condition_df_NOR_ALT
  
  dna_condition_df_NOR_ALT.dt<-data.table(dna_condition_df_NOR_ALT,
                                          key=c("REAL_TILE_Plus_carried_variants","Cell_Type"))
  
  dna_ALT_median<-unique(as.data.frame(dna_condition_df_NOR_ALT.dt[,.(median_Median_value=median(Median_value)),
                                                                   by=key(dna_condition_df_NOR_ALT.dt)],stringsAsFactors=F))
  
  
  cat("dna_ALT_median_0\n")
  cat(str(dna_ALT_median))
  cat("\n")
  
  List_medians[['dna_ALT']]<-dna_ALT_median
  
  
  dna_condition_df_NOR_ALT_wide<-pivot_wider(dna_condition_df_NOR_ALT,
                                             id_cols=REAL_TILE_Plus_carried_variants,
                                             names_from=master_sample,
                                             values_from=Median_value)
  
  cat("dna_condition_df_NOR_ALT_wide\n")
  cat(str(dna_condition_df_NOR_ALT_wide))
  cat("\n")
  
  dna_condition_df_NOR_ALT_wide.matrix<-as.matrix(dna_condition_df_NOR_ALT_wide[,-1])
  
  row.names(dna_condition_df_NOR_ALT_wide.matrix)<-dna_condition_df_NOR_ALT_wide$REAL_TILE_Plus_carried_variants
  
  cat("dna_condition_df_NOR_ALT_wide.matrix\n")
  cat(str(dna_condition_df_NOR_ALT_wide.matrix))
  cat("\n")
  
  
  
  
  rna_condition_df_NOR_ALT<-dna_rna_condition_df_NOR_ALT[which(dna_rna_condition_df_NOR_ALT$type == "cDNA"),]
  
  # cat("rna_condition_df_NOR_ALT_1\n")
  # cat(str(rna_condition_df_NOR_ALT))
  # cat("\n")
  
  rna_condition_df_NOR_ALT<-droplevels(rna_condition_df_NOR_ALT)
  
  cat("rna_condition_df_NOR_ALT_2\n")
  cat(str(rna_condition_df_NOR_ALT))
  cat("\n")
  
  List_values[['rna_ALT']]<-rna_condition_df_NOR_ALT
  
  rna_condition_df_NOR_ALT.dt<-data.table(rna_condition_df_NOR_ALT,
                                          key=c("REAL_TILE_Plus_carried_variants","Cell_Type"))
  
  rna_ALT_median<-unique(as.data.frame(rna_condition_df_NOR_ALT.dt[,.(median_Median_value=median(Median_value)),
                                                                   by=key(rna_condition_df_NOR_ALT.dt)],stringsAsFactors=F))
  
  
  cat("rna_ALT_median_0\n")
  cat(str(rna_ALT_median))
  cat("\n")
  
  List_medians[['rna_ALT']]<-rna_ALT_median
  
  
  rna_condition_df_NOR_ALT_wide<-pivot_wider(rna_condition_df_NOR_ALT,
                                             id_cols=REAL_TILE_Plus_carried_variants,
                                             names_from=master_sample,
                                             values_from=Median_value)
  
  cat("rna_condition_df_NOR_ALT_wide\n")
  cat(str(rna_condition_df_NOR_ALT_wide))
  cat("\n")
  
  rna_condition_df_NOR_ALT_wide.matrix<-as.matrix(rna_condition_df_NOR_ALT_wide[,-1])
  
  row.names(rna_condition_df_NOR_ALT_wide.matrix)<-rna_condition_df_NOR_ALT_wide$REAL_TILE_Plus_carried_variants
  
  cat("rna_condition_df_NOR_ALT_wide.matrix\n")
  cat(str(rna_condition_df_NOR_ALT_wide.matrix))
  cat("\n")
  
  
  #### ASE ----
  
  
  
  ASE_matrix<-(rna_condition_df_NOR_ALT_wide.matrix+0.0000001/dna_condition_df_NOR_ALT_wide.matrix+0.0000001)/(rna_condition_df_NOR_REF_wide.matrix+0.0000001/dna_condition_df_NOR_REF_wide.matrix+0.0000001)
  
  cat("--------------------------------------------------------------------------------------------->ASE_matrix\n")
  cat(str(ASE_matrix))
  cat("\n")
  
  ASE_df<-as.data.frame(ASE_matrix, stringsAsFactors=F)
  
  
  ASE_df$REAL_TILE_Plus_carried_variants<-row.names(ASE_df)
  
  cat("ASE_df\n")
  cat(str(ASE_df))
  cat("\n")
  
  ASE_df.m<-melt(ASE_df, id.cols="REAL_TILE_Plus_carried_variants")
  
  cat("ASE_df.m_1\n")
  cat(str(ASE_df.m))
  cat("\n")
  
  colnames(ASE_df.m)[which(colnames(ASE_df.m) == "variable")]<-"master_sample"
  colnames(ASE_df.m)[which(colnames(ASE_df.m) == "value")]<-"ASE"
  
  
  ASE_df.m$master_sample<-factor(ASE_df.m$master_sample,
                                         levels=levels_master_samples,ordered=T)
  
  
  cat("ASE_df.m_2\n")
  cat(str(ASE_df.m))
  cat("\n")
  
  ASE_df.m<-merge(ASE_df.m,
                          Annot_DEF_subset,
                          by="master_sample")
  
  
  cat("ASE_df.m_3\n")
  cat(str(ASE_df.m))
  cat("\n")
  
  
  ### NO NA
  
  ASE_df.m_NO_NA<-ASE_df.m[!is.na(ASE_df.m$ASE),]
  
  
  cat("ASE_df.m_NO_NA_1\n")
  cat(str(ASE_df.m_NO_NA))
  cat("\n")
  
  ### Infinite reverted
  
  ASE_df.m_NO_NA$ASE[(is.infinite(ASE_df.m_NO_NA$ASE) &  ASE_df.m_NO_NA$ASE > 0)]<-60
  
  # ASE_df.m_NO_NA<-ASE_df.m_NO_NA[-(is.infinite(ASE_df.m_NO_NA$ASE) &  ASE_df.m_NO_NA$ASE < 0),]
  
  
  
  # ASE_df.m_NO_NA$ASE[(is.infinite(ASE_df.m_NO_NA$ASE) &  ASE_df.m_NO_NA$ASE < 0)]<--1*10000
  

  cat("ASE_df.m_NO_NA_2_infinite_reverted\n")
  cat(str(ASE_df.m_NO_NA))
  cat("\n")
  
  
  
  
  
  A<-summary(ASE_df.m_NO_NA$ASE)
  
  cat("summary_ASE\n")
  cat(sprintf(as.character(names(A))))
  cat("\n")
  cat(sprintf(as.character(A)))
  cat("\n")
  
  
  check<-ASE_df.m_NO_NA[which(ASE_df.m_NO_NA$ASE >1.2 | ASE_df.m_NO_NA$ASE <0.8),]
  
  cat("check_ASE\n")
  cat(str(check))
  cat("\n")
  cat(sprintf(as.character(names(summary(check$master_sample)))))
  cat("\n")
  cat(sprintf(as.character(summary(check$master_sample))))
  cat("\n")
  
  A<-summary(check$ASE)
  
  cat("summary_ASE\n")
  cat(sprintf(as.character(names(A))))
  cat("\n")
  cat(sprintf(as.character(A)))
  cat("\n")
  
  
  check<-ASE_df.m_NO_NA[is.na(ASE_df.m_NO_NA$ASE),]
  
  cat("check_2\n")
  cat(str(check))
  cat("\n")
  
  
  List_values[['ASE']]<-ASE_df.m_NO_NA
  
  ASE_df.m_NO_NA.dt<-data.table(ASE_df.m_NO_NA,
                                        key=c("REAL_TILE_Plus_carried_variants","Cell_Type"))
  
  ASE_median<-unique(as.data.frame(ASE_df.m_NO_NA.dt[,.(median_ASE=median(ASE)),
                                                                     by=key(ASE_df.m_NO_NA.dt)],stringsAsFactors=F))
  
  
  cat("ASE_median_0\n")
  cat(str(ASE_median))
  cat("\n")
  
  List_medians[['ASE']]<-ASE_median
  
  
  check_ASE<-ASE_median[which(ASE_median$median_ASE >1.2 | ASE_median$median_ASE <0.8),]
  
  cat("check_ASE_ASE\n")
  cat(str(check_ASE))
  cat("\n")
  cat(sprintf(as.character(names(summary(check_ASE$Cell_Type)))))
  cat("\n")
  cat(sprintf(as.character(summary(check_ASE$Cell_Type))))
  cat("\n")
  # cat(sprintf(as.character(names(summary(as.factor(check_ASE$REAL_TILE_Plus_carried_variants))))))
  # cat("\n")
  # cat(sprintf(as.character(summary(as.factor(check_ASE$REAL_TILE_Plus_carried_variants)))))
  # cat("\n")
  
  A<-summary(check_ASE$median_ASE)
  
  cat("summary_median_ASE\n")
  cat(sprintf(as.character(names(A))))
  cat("\n")
  cat(sprintf(as.character(A)))
  cat("\n")
  
  
  
  
  check_DEF<-check_LogFC[which(check_LogFC$REAL_TILE_Plus_carried_variants%in%check_ASE$REAL_TILE_Plus_carried_variants),]
  
  
  
  cat("check_DEF_ASE\n")
  cat(str(check_DEF))
  cat("\n")
  
  
  
  #### dna.REF.vs.dna.ALT ----
  
  
  
  dna.REF.vs.dna.ALT_matrix<-(dna_condition_df_NOR_REF_wide.matrix+0.0000001)/(dna_condition_df_NOR_ALT_wide.matrix+0.0000001)
  
  cat("--------------------------------------------------------------------------------------------->dna.REF.vs.dna.ALT_matrix\n")
  cat(str(dna.REF.vs.dna.ALT_matrix))
  cat("\n")
  
  dna.REF.vs.dna.ALT_df<-as.data.frame(dna.REF.vs.dna.ALT_matrix, stringsAsFactors=F)
  
  
  dna.REF.vs.dna.ALT_df$REAL_TILE_Plus_carried_variants<-row.names(dna.REF.vs.dna.ALT_df)
  
  cat("dna.REF.vs.dna.ALT_df\n")
  cat(str(dna.REF.vs.dna.ALT_df))
  cat("\n")
  
  dna.REF.vs.dna.ALT_df.m<-melt(dna.REF.vs.dna.ALT_df, id.cols="REAL_TILE_Plus_carried_variants")
  
  cat("dna.REF.vs.dna.ALT_df.m_1\n")
  cat(str(dna.REF.vs.dna.ALT_df.m))
  cat("\n")
  
  colnames(dna.REF.vs.dna.ALT_df.m)[which(colnames(dna.REF.vs.dna.ALT_df.m) == "variable")]<-"master_sample"
  colnames(dna.REF.vs.dna.ALT_df.m)[which(colnames(dna.REF.vs.dna.ALT_df.m) == "value")]<-"dna.REF.vs.dna.ALT"
  
  
  dna.REF.vs.dna.ALT_df.m$master_sample<-factor(dna.REF.vs.dna.ALT_df.m$master_sample,
                                                levels=levels_master_samples,ordered=T)
  
  
  cat("dna.REF.vs.dna.ALT_df.m_2\n")
  cat(str(dna.REF.vs.dna.ALT_df.m))
  cat("\n")
  
  dna.REF.vs.dna.ALT_df.m<-merge(dna.REF.vs.dna.ALT_df.m,
                                 Annot_DEF_subset,
                                 by="master_sample")
  
  
  cat("dna.REF.vs.dna.ALT_df.m_3\n")
  cat(str(dna.REF.vs.dna.ALT_df.m))
  cat("\n")
  
  
  ### NO NA
  
  dna.REF.vs.dna.ALT_df.m_NO_NA<-dna.REF.vs.dna.ALT_df.m[!is.na(dna.REF.vs.dna.ALT_df.m$dna.REF.vs.dna.ALT),]
  
  
  cat("dna.REF.vs.dna.ALT_df.m_NO_NA_1\n")
  cat(str(dna.REF.vs.dna.ALT_df.m_NO_NA))
  cat("\n")
  
  ### Infinite reverted
  
  dna.REF.vs.dna.ALT_df.m_NO_NA$dna.REF.vs.dna.ALT[(is.infinite(dna.REF.vs.dna.ALT_df.m_NO_NA$dna.REF.vs.dna.ALT) &  dna.REF.vs.dna.ALT_df.m_NO_NA$dna.REF.vs.dna.ALT > 0)]<-60
  #dna.REF.vs.dna.ALT_df.m_NO_NA$dna.REF.vs.dna.ALT[(is.infinite(dna.REF.vs.dna.ALT_df.m_NO_NA$dna.REF.vs.dna.ALT) &  dna.REF.vs.dna.ALT_df.m_NO_NA$dna.REF.vs.dna.ALT < 0)]<--1*10000
  
  # dna.REF.vs.dna.ALT_df.m_NO_NA<-dna.REF.vs.dna.ALT_df.m_NO_NA[-(is.infinite(dna.REF.vs.dna.ALT_df.m_NO_NA$dna.REF.vs.dna.ALT) &  dna.REF.vs.dna.ALT_df.m_NO_NA$dna.REF.vs.dna.ALT < 0),]
  
  
  
  

  cat("dna.REF.vs.dna.ALT_df.m_NO_NA_2_infinite_reverted\n")
  cat(str(dna.REF.vs.dna.ALT_df.m_NO_NA))
  cat("\n")
  
  
  
  
  
  A<-summary(dna.REF.vs.dna.ALT_df.m_NO_NA$dna.REF.vs.dna.ALT)
  
  cat("summary_dna.REF.vs.dna.ALT\n")
  cat(sprintf(as.character(names(A))))
  cat("\n")
  cat(sprintf(as.character(A)))
  cat("\n")
  
  
  # check<-dna.REF.vs.dna.ALT_df.m_NO_NA[which(dna.REF.vs.dna.ALT_df.m_NO_NA$dna.REF.vs.dna.ALT >1.2 | dna.REF.vs.dna.ALT_df.m_NO_NA$dna.REF.vs.dna.ALT <0.8),]
  # 
  # cat("check_dna.REF.vs.dna.ALT\n")
  # cat(str(check))
  # cat("\n")
  # cat(sprintf(as.character(names(summary(check$master_sample)))))
  # cat("\n")
  # cat(sprintf(as.character(summary(check$master_sample))))
  # cat("\n")
  # 
  # A<-summary(check$dna.REF.vs.dna.ALT)
  # 
  # cat("summary_dna.REF.vs.dna.ALT\n")
  # cat(sprintf(as.character(names(A))))
  # cat("\n")
  # cat(sprintf(as.character(A)))
  # cat("\n")
  
  
  check<-dna.REF.vs.dna.ALT_df.m_NO_NA[is.na(dna.REF.vs.dna.ALT_df.m_NO_NA$dna.REF.vs.dna.ALT),]
  
  cat("check_2\n")
  cat(str(check))
  cat("\n")
  
  
  List_values[['dna.REF.vs.dna.ALT']]<-dna.REF.vs.dna.ALT_df.m_NO_NA
  
  dna.REF.vs.dna.ALT_df.m_NO_NA.dt<-data.table(dna.REF.vs.dna.ALT_df.m_NO_NA,
                                               key=c("REAL_TILE_Plus_carried_variants","Cell_Type"))
  
  dna.REF.vs.dna.ALT_median<-unique(as.data.frame(dna.REF.vs.dna.ALT_df.m_NO_NA.dt[,.(median_dna.REF.vs.dna.ALT=median(dna.REF.vs.dna.ALT)),
                                                                                   by=key(dna.REF.vs.dna.ALT_df.m_NO_NA.dt)],stringsAsFactors=F))
  
  
  cat("dna.REF.vs.dna.ALT_median_0\n")
  cat(str(dna.REF.vs.dna.ALT_median))
  cat("\n")
  
  List_medians[['dna.REF.vs.dna.ALT']]<-dna.REF.vs.dna.ALT_median
  
  
  # check_dna.REF.vs.dna.ALT<-dna.REF.vs.dna.ALT_median[which(dna.REF.vs.dna.ALT_median$median_dna.REF.vs.dna.ALT >1.2 | dna.REF.vs.dna.ALT_median$median_dna.REF.vs.dna.ALT <0.8),]
  # 
  # cat("check_dna.REF.vs.dna.ALT_dna.REF.vs.dna.ALT\n")
  # cat(str(check_dna.REF.vs.dna.ALT))
  # cat("\n")
  # cat(sprintf(as.character(names(summary(check_dna.REF.vs.dna.ALT$Cell_Type)))))
  # cat("\n")
  # cat(sprintf(as.character(summary(check_dna.REF.vs.dna.ALT$Cell_Type))))
  # cat("\n")
  # # cat(sprintf(as.character(names(summary(as.factor(check_dna.REF.vs.dna.ALT$REAL_TILE_Plus_carried_variants))))))
  # # cat("\n")
  # # cat(sprintf(as.character(summary(as.factor(check_dna.REF.vs.dna.ALT$REAL_TILE_Plus_carried_variants)))))
  # # cat("\n")
  # 
  # A<-summary(check_dna.REF.vs.dna.ALT$median_dna.REF.vs.dna.ALT)
  # 
  # cat("summary_median_dna.REF.vs.dna.ALT\n")
  # cat(sprintf(as.character(names(A))))
  # cat("\n")
  # cat(sprintf(as.character(A)))
  # cat("\n")
  
  
  #### rna.REF.vs.rna.ALT ----
  
  
  
  rna.REF.vs.rna.ALT_matrix<-(rna_condition_df_NOR_REF_wide.matrix+0.0000001)/(rna_condition_df_NOR_ALT_wide.matrix+0.0000001)
  
  cat("--------------------------------------------------------------------------------------------->rna.REF.vs.rna.ALT_matrix\n")
  cat(str(rna.REF.vs.rna.ALT_matrix))
  cat("\n")
  
  rna.REF.vs.rna.ALT_df<-as.data.frame(rna.REF.vs.rna.ALT_matrix, stringsAsFactors=F)
  
  
  rna.REF.vs.rna.ALT_df$REAL_TILE_Plus_carried_variants<-row.names(rna.REF.vs.rna.ALT_df)
  
  cat("rna.REF.vs.rna.ALT_df\n")
  cat(str(rna.REF.vs.rna.ALT_df))
  cat("\n")
  
  rna.REF.vs.rna.ALT_df.m<-melt(rna.REF.vs.rna.ALT_df, id.cols="REAL_TILE_Plus_carried_variants")
  
  cat("rna.REF.vs.rna.ALT_df.m_1\n")
  cat(str(rna.REF.vs.rna.ALT_df.m))
  cat("\n")
  
  colnames(rna.REF.vs.rna.ALT_df.m)[which(colnames(rna.REF.vs.rna.ALT_df.m) == "variable")]<-"master_sample"
  colnames(rna.REF.vs.rna.ALT_df.m)[which(colnames(rna.REF.vs.rna.ALT_df.m) == "value")]<-"rna.REF.vs.rna.ALT"
  
  
  rna.REF.vs.rna.ALT_df.m$master_sample<-factor(rna.REF.vs.rna.ALT_df.m$master_sample,
                                                levels=levels_master_samples,ordered=T)
  
  
  cat("rna.REF.vs.rna.ALT_df.m_2\n")
  cat(str(rna.REF.vs.rna.ALT_df.m))
  cat("\n")
  
  rna.REF.vs.rna.ALT_df.m<-merge(rna.REF.vs.rna.ALT_df.m,
                                 Annot_DEF_subset,
                                 by="master_sample")
  
  
  cat("rna.REF.vs.rna.ALT_df.m_3\n")
  cat(str(rna.REF.vs.rna.ALT_df.m))
  cat("\n")
  
  
  ### NO NA
  
  rna.REF.vs.rna.ALT_df.m_NO_NA<-rna.REF.vs.rna.ALT_df.m[!is.na(rna.REF.vs.rna.ALT_df.m$rna.REF.vs.rna.ALT),]
  
  
  cat("rna.REF.vs.rna.ALT_df.m_NO_NA_1\n")
  cat(str(rna.REF.vs.rna.ALT_df.m_NO_NA))
  cat("\n")
  
  ### Infinite reverted
  
  rna.REF.vs.rna.ALT_df.m_NO_NA$rna.REF.vs.rna.ALT[(is.infinite(rna.REF.vs.rna.ALT_df.m_NO_NA$rna.REF.vs.rna.ALT) &  rna.REF.vs.rna.ALT_df.m_NO_NA$rna.REF.vs.rna.ALT > 0)]<-60
  # rna.REF.vs.rna.ALT_df.m_NO_NA<-rna.REF.vs.rna.ALT_df.m_NO_NA[-(is.infinite(rna.REF.vs.rna.ALT_df.m_NO_NA$rna.REF.vs.rna.ALT) &  rna.REF.vs.rna.ALT_df.m_NO_NA$rna.REF.vs.rna.ALT < 0),]
  
  #  rna.REF.vs.rna.ALT_df.m_NO_NA$rna.REF.vs.rna.ALT[(is.infinite(rna.REF.vs.rna.ALT_df.m_NO_NA$rna.REF.vs.rna.ALT) &  rna.REF.vs.rna.ALT_df.m_NO_NA$rna.REF.vs.rna.ALT < 0)]<--1*10000
  

  cat("rna.REF.vs.rna.ALT_df.m_NO_NA_2_infinite_reverted\n")
  cat(str(rna.REF.vs.rna.ALT_df.m_NO_NA))
  cat("\n")
  
  
  
  
  
  A<-summary(rna.REF.vs.rna.ALT_df.m_NO_NA$rna.REF.vs.rna.ALT)
  
  cat("summary_rna.REF.vs.rna.ALT\n")
  cat(sprintf(as.character(names(A))))
  cat("\n")
  cat(sprintf(as.character(A)))
  cat("\n")
  
  
  # check<-rna.REF.vs.rna.ALT_df.m_NO_NA[which(rna.REF.vs.rna.ALT_df.m_NO_NA$rna.REF.vs.rna.ALT >1.2 | rna.REF.vs.rna.ALT_df.m_NO_NA$rna.REF.vs.rna.ALT <0.8),]
  # 
  # cat("check_rna.REF.vs.rna.ALT\n")
  # cat(str(check))
  # cat("\n")
  # cat(sprintf(as.character(names(summary(check$master_sample)))))
  # cat("\n")
  # cat(sprintf(as.character(summary(check$master_sample))))
  # cat("\n")
  # 
  # A<-summary(check$rna.REF.vs.rna.ALT)
  # 
  # cat("summary_rna.REF.vs.rna.ALT\n")
  # cat(sprintf(as.character(names(A))))
  # cat("\n")
  # cat(sprintf(as.character(A)))
  # cat("\n")
  
  
  check<-rna.REF.vs.rna.ALT_df.m_NO_NA[is.na(rna.REF.vs.rna.ALT_df.m_NO_NA$rna.REF.vs.rna.ALT),]
  
  cat("check_2\n")
  cat(str(check))
  cat("\n")
  
  
  List_values[['rna.REF.vs.rna.ALT']]<-rna.REF.vs.rna.ALT_df.m_NO_NA
  
  rna.REF.vs.rna.ALT_df.m_NO_NA.dt<-data.table(rna.REF.vs.rna.ALT_df.m_NO_NA,
                                               key=c("REAL_TILE_Plus_carried_variants","Cell_Type"))
  
  rna.REF.vs.rna.ALT_median<-unique(as.data.frame(rna.REF.vs.rna.ALT_df.m_NO_NA.dt[,.(median_rna.REF.vs.rna.ALT=median(rna.REF.vs.rna.ALT)),
                                                                                   by=key(rna.REF.vs.rna.ALT_df.m_NO_NA.dt)],stringsAsFactors=F))
  
  
  cat("rna.REF.vs.rna.ALT_median_0\n")
  cat(str(rna.REF.vs.rna.ALT_median))
  cat("\n")
  
  List_medians[['rna.REF.vs.rna.ALT']]<-rna.REF.vs.rna.ALT_median
  
  
  # check_rna.REF.vs.rna.ALT<-rna.REF.vs.rna.ALT_median[which(rna.REF.vs.rna.ALT_median$median_rna.REF.vs.rna.ALT >1.2 | rna.REF.vs.rna.ALT_median$median_rna.REF.vs.rna.ALT <0.8),]
  # 
  # cat("check_rna.REF.vs.rna.ALT_rna.REF.vs.rna.ALT\n")
  # cat(str(check_rna.REF.vs.rna.ALT))
  # cat("\n")
  # cat(sprintf(as.character(names(summary(check_rna.REF.vs.rna.ALT$Cell_Type)))))
  # cat("\n")
  # cat(sprintf(as.character(summary(check_rna.REF.vs.rna.ALT$Cell_Type))))
  # cat("\n")
  # # cat(sprintf(as.character(names(summary(as.factor(check_rna.REF.vs.rna.ALT$REAL_TILE_Plus_carried_variants))))))
  # # cat("\n")
  # # cat(sprintf(as.character(summary(as.factor(check_rna.REF.vs.rna.ALT$REAL_TILE_Plus_carried_variants)))))
  # # cat("\n")
  # 
  # A<-summary(check_rna.REF.vs.rna.ALT$median_rna.REF.vs.rna.ALT)
  # 
  # cat("summary_median_rna.REF.vs.rna.ALT\n")
  # cat(sprintf(as.character(names(A))))
  # cat("\n")
  # cat(sprintf(as.character(A)))
  # cat("\n")
  
  
  
  
  ##### SAVE LISTS with parameters ----
  
  setwd(out2)
  
  filename_1<-paste("list_values_post_QC",".rds", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  saveRDS(List_values, file=filename_1)
  
  
  filename_1<-paste("list_medians_post_QC",".rds", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  saveRDS(List_medians, file=filename_1)
}

# #####################################################################################################################################################################################################################
# quit(status=1)


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
    make_option(c("--K562_replicates_QC_PASS"), type="character", default=NULL, 
                metavar="NEG.CTRLS", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--CHRF_replicates_QC_PASS"), type="character", default=NULL, 
                metavar="NEG.CTRLS", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--HL60_replicates_QC_PASS"), type="character", default=NULL, 
                metavar="NEG.CTRLS", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--THP1_replicates_QC_PASS"), type="character", default=NULL, 
                metavar="NEG.CTRLS",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="out", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out2"), type="character", default=NULL, 
                metavar="out", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "137_MPRA_normalization_and_filtering_Rscript_v2.R
                        --CHARAC_TABLE FILE.txt
                        --replicas charac
                        --type1 type1
                        --type2 type2
                        --barcodes_per_tile integer
                        --pvalThreshold integer
                        --FDRThreshold integer
                        --EquivalenceTable FILE.txt
                        --sharpr2Threshold charac",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  Parameter_recalculation(opt)
 
  
  
  
}


###########################################################################

system.time( main() )
