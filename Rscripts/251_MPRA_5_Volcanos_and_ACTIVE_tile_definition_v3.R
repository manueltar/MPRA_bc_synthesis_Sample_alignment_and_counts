

suppressMessages(library("zoo", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("biomaRt", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("Sushi", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
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
library("gtools", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
suppressMessages(library("labeling", lib.loc="/nfs/team151/software/manuel_R_libs_4_1//"))
suppressMessages(library("cowplot", lib.loc="/nfs/team151/software/manuel_R_libs_4_1//"))
suppressMessages(library("farver", lib.loc="/nfs/team151/software/manuel_R_libs_4_1//"))
suppressMessages(library("ggeasy", lib.loc="/nfs/team151/software/manuel_R_libs_4_1//"))


opt = NULL

options(warn = 1)

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
  
  #### READ and transform out2 ----
  
  out2 = opt$out2
  
  cat("out2_\n")
  cat(sprintf(as.character(out2)))
  cat("\n")
  
  #### enhancer_result_K562 ----
  
  enhancer_result_K562<-as.data.frame(fread(file=opt$enhancer_result_K562, sep="\t", header = T, fill=TRUE), stringsAsFactors = F)
  
  
  colnames(enhancer_result_K562)<-c("REAL_TILE_Plus_carried_variants","statistic","control","zscore","mad.score","pval.empirical","pval.mad","pval.zscore")
  
  
  enhancer_result_K562$logpval_enhancer_pval.empirical<--1*log10(enhancer_result_K562$pval.empirical)
  enhancer_result_K562$enhancer_logpval<--1*log10(enhancer_result_K562$pval.mad)
  
  
  cat("enhancer_result_K562\n")
  cat(str(enhancer_result_K562))
  cat("\n")
  
  enhancer_result_K562_subset<-unique(enhancer_result_K562[,c(which(colnames(enhancer_result_K562) == "REAL_TILE_Plus_carried_variants"),
                                                       which(colnames(enhancer_result_K562) == "logpval_enhancer_pval.empirical"),
                                                       which(colnames(enhancer_result_K562) == "enhancer_logpval"))])
  
  
  cat("enhancer_result_K562_subset\n")
  cat(str(enhancer_result_K562_subset))
  cat("\n")
 
  #### ASE_result_K562 ----
  
  ASE_result_K562<-as.data.frame(fread(file=opt$ASE_result_K562, sep="\t", header = T, fill=TRUE), stringsAsFactors = F)
  
  colnames(ASE_result_K562)<-c("REAL_TILE_Plus_carried_variants","statistic","pval","fdr","df.test","df.dna","df.rna.full","df.rna.red","logFC")
  ASE_result_K562$ASE_logpval<--1*log10(ASE_result_K562$pval)
  
  
  cat("ASE_result_K562\n")
  cat(str(ASE_result_K562))
  cat("\n")
  
  ASE_result_K562_subset<-unique(ASE_result_K562[,c(which(colnames(ASE_result_K562) == "REAL_TILE_Plus_carried_variants"),
                                                              which(colnames(ASE_result_K562) == "fdr"),
                                                              which(colnames(ASE_result_K562) == "ASE_logpval"))])
  
  
  cat("ASE_result_K562_subset\n")
  cat(str(ASE_result_K562_subset))
  cat("\n")
  
  K562_result<-merge(enhancer_result_K562_subset,
                     ASE_result_K562_subset,
                     by="REAL_TILE_Plus_carried_variants")
  
  
  K562_result$Cell_Type<-"K562"
  
  
  cat("K562_result_0\n")
  cat(str(K562_result))
  cat("\n")
  cat(sprintf(as.character(names(summary(K562_result$logpval_enhancer_pval.empirical)))))
  cat("\n")
  cat(sprintf(as.character(summary(K562_result$logpval_enhancer_pval.empirical))))
  cat("\n")
  cat(sprintf(as.character(names(summary(K562_result$enhancer_logpval)))))
  cat("\n")
  cat(sprintf(as.character(summary(K562_result$enhancer_logpval))))
  cat("\n")
  cat(sprintf(as.character(names(summary(K562_result$fdr)))))
  cat("\n")
  cat(sprintf(as.character(summary(K562_result$fdr))))
  cat("\n")
  cat(sprintf(as.character(names(summary(K562_result$ASE_logpval)))))
  cat("\n")
  cat(sprintf(as.character(summary(K562_result$ASE_logpval))))
  cat("\n")
  
  
  K562_result$logpval_enhancer_pval.empirical[!is.finite(K562_result$logpval_enhancer_pval.empirical)]<-90
  K562_result$enhancer_logpval[!is.finite(K562_result$enhancer_logpval)]<-300
  K562_result$ASE_logpval[!is.finite(K562_result$ASE_logpval)]<-90
  
  
  cat("K562_result_1\n")
  cat(str(K562_result))
  cat("\n")
  cat(sprintf(as.character(names(summary(K562_result$logpval_enhancer_pval.empirical)))))
  cat("\n")
  cat(sprintf(as.character(summary(K562_result$logpval_enhancer_pval.empirical))))
  cat("\n")
  cat(sprintf(as.character(names(summary(K562_result$enhancer_logpval)))))
  cat("\n")
  cat(sprintf(as.character(summary(K562_result$enhancer_logpval))))
  cat("\n")
  cat(sprintf(as.character(names(summary(K562_result$fdr)))))
  cat("\n")
  cat(sprintf(as.character(summary(K562_result$fdr))))
  cat("\n")
  cat(sprintf(as.character(names(summary(K562_result$ASE_logpval)))))
  cat("\n")
  cat(sprintf(as.character(summary(K562_result$ASE_logpval))))
  cat("\n")
 
  
  
  
  #### enhancer_result_CHRF ----
  
  enhancer_result_CHRF<-as.data.frame(fread(file=opt$enhancer_result_CHRF, sep="\t", header = T, fill=TRUE), stringsAsFactors = F)
  
  
  colnames(enhancer_result_CHRF)<-c("REAL_TILE_Plus_carried_variants","statistic","control","zscore","mad.score","pval.empirical","pval.mad","pval.zscore")
  
  
  enhancer_result_CHRF$logpval_enhancer_pval.empirical<--1*log10(enhancer_result_CHRF$pval.empirical)
  enhancer_result_CHRF$enhancer_logpval<--1*log10(enhancer_result_CHRF$pval.mad)
  
  
  cat("enhancer_result_CHRF\n")
  cat(str(enhancer_result_CHRF))
  cat("\n")
  
  enhancer_result_CHRF_subset<-unique(enhancer_result_CHRF[,c(which(colnames(enhancer_result_CHRF) == "REAL_TILE_Plus_carried_variants"),
                                                              which(colnames(enhancer_result_CHRF) == "logpval_enhancer_pval.empirical"),
                                                              which(colnames(enhancer_result_CHRF) == "enhancer_logpval"))])
  
  
  cat("enhancer_result_CHRF_subset\n")
  cat(str(enhancer_result_CHRF_subset))
  cat("\n")
  
  #### ASE_result_CHRF ----
  
  ASE_result_CHRF<-as.data.frame(fread(file=opt$ASE_result_CHRF, sep="\t", header = T, fill=TRUE), stringsAsFactors = F)
  
  colnames(ASE_result_CHRF)<-c("REAL_TILE_Plus_carried_variants","statistic","pval","fdr","df.test","df.dna","df.rna.full","df.rna.red","logFC")
  ASE_result_CHRF$ASE_logpval<--1*log10(ASE_result_CHRF$pval)
  
  
  cat("ASE_result_CHRF\n")
  cat(str(ASE_result_CHRF))
  cat("\n")
  
  ASE_result_CHRF_subset<-unique(ASE_result_CHRF[,c(which(colnames(ASE_result_CHRF) == "REAL_TILE_Plus_carried_variants"),
                                                    which(colnames(ASE_result_CHRF) == "fdr"),
                                                    which(colnames(ASE_result_CHRF) == "ASE_logpval"))])
  
  
  cat("ASE_result_CHRF_subset\n")
  cat(str(ASE_result_CHRF_subset))
  cat("\n")
  
  CHRF_result<-merge(enhancer_result_CHRF_subset,
                     ASE_result_CHRF_subset,
                     by="REAL_TILE_Plus_carried_variants")
  
  
  CHRF_result$Cell_Type<-"CHRF"
  
  
  cat("CHRF_result_0\n")
  cat(str(CHRF_result))
  cat("\n")
  cat(sprintf(as.character(names(summary(CHRF_result$logpval_enhancer_pval.empirical)))))
  cat("\n")
  cat(sprintf(as.character(summary(CHRF_result$logpval_enhancer_pval.empirical))))
  cat("\n")
  cat(sprintf(as.character(names(summary(CHRF_result$enhancer_logpval)))))
  cat("\n")
  cat(sprintf(as.character(summary(CHRF_result$enhancer_logpval))))
  cat("\n")
  cat(sprintf(as.character(names(summary(CHRF_result$fdr)))))
  cat("\n")
  cat(sprintf(as.character(summary(CHRF_result$fdr))))
  cat("\n")
  cat(sprintf(as.character(names(summary(CHRF_result$ASE_logpval)))))
  cat("\n")
  cat(sprintf(as.character(summary(CHRF_result$ASE_logpval))))
  cat("\n")
  
  
  CHRF_result$logpval_enhancer_pval.empirical[!is.finite(CHRF_result$logpval_enhancer_pval.empirical)]<-90
  CHRF_result$enhancer_logpval[!is.finite(CHRF_result$enhancer_logpval)]<-300
  CHRF_result$ASE_logpval[!is.finite(CHRF_result$ASE_logpval)]<-90
  
  
  cat("CHRF_result_1\n")
  cat(str(CHRF_result))
  cat("\n")
  cat(sprintf(as.character(names(summary(CHRF_result$logpval_enhancer_pval.empirical)))))
  cat("\n")
  cat(sprintf(as.character(summary(CHRF_result$logpval_enhancer_pval.empirical))))
  cat("\n")
  cat(sprintf(as.character(names(summary(CHRF_result$enhancer_logpval)))))
  cat("\n")
  cat(sprintf(as.character(summary(CHRF_result$enhancer_logpval))))
  cat("\n")
  cat(sprintf(as.character(names(summary(CHRF_result$fdr)))))
  cat("\n")
  cat(sprintf(as.character(summary(CHRF_result$fdr))))
  cat("\n")
  cat(sprintf(as.character(names(summary(CHRF_result$ASE_logpval)))))
  cat("\n")
  cat(sprintf(as.character(summary(CHRF_result$ASE_logpval))))
  cat("\n")
  
  
  
  
  #### enhancer_result_HL60 ----

  enhancer_result_HL60<-as.data.frame(fread(file=opt$enhancer_result_HL60, sep="\t", header = T, fill=TRUE), stringsAsFactors = F)


  colnames(enhancer_result_HL60)<-c("REAL_TILE_Plus_carried_variants","statistic","control","zscore","mad.score","pval.empirical","pval.mad","pval.zscore")


  enhancer_result_HL60$logpval_enhancer_pval.empirical<--1*log10(enhancer_result_HL60$pval.empirical)
  enhancer_result_HL60$enhancer_logpval<--1*log10(enhancer_result_HL60$pval.mad)


  cat("enhancer_result_HL60\n")
  cat(str(enhancer_result_HL60))
  cat("\n")

  enhancer_result_HL60_subset<-unique(enhancer_result_HL60[,c(which(colnames(enhancer_result_HL60) == "REAL_TILE_Plus_carried_variants"),
                                                              which(colnames(enhancer_result_HL60) == "logpval_enhancer_pval.empirical"),
                                                              which(colnames(enhancer_result_HL60) == "enhancer_logpval"))])


  cat("enhancer_result_HL60_subset\n")
  cat(str(enhancer_result_HL60_subset))
  cat("\n")

  #### ASE_result_HL60 ----

  ASE_result_HL60<-as.data.frame(fread(file=opt$ASE_result_HL60, sep="\t", header = T, fill=TRUE), stringsAsFactors = F)

  colnames(ASE_result_HL60)<-c("REAL_TILE_Plus_carried_variants","statistic","pval","fdr","df.test","df.dna","df.rna.full","df.rna.red","logFC")
  ASE_result_HL60$ASE_logpval<--1*log10(ASE_result_HL60$pval)


  cat("ASE_result_HL60\n")
  cat(str(ASE_result_HL60))
  cat("\n")

  ASE_result_HL60_subset<-unique(ASE_result_HL60[,c(which(colnames(ASE_result_HL60) == "REAL_TILE_Plus_carried_variants"),
                                                    which(colnames(ASE_result_HL60) == "fdr"),
                                                    which(colnames(ASE_result_HL60) == "ASE_logpval"))])


  cat("ASE_result_HL60_subset\n")
  cat(str(ASE_result_HL60_subset))
  cat("\n")

  HL60_result<-merge(enhancer_result_HL60_subset,
                     ASE_result_HL60_subset,
                     by="REAL_TILE_Plus_carried_variants")


  HL60_result$Cell_Type<-"HL60"


  cat("HL60_result_0\n")
  cat(str(HL60_result))
  cat("\n")
  cat(sprintf(as.character(names(summary(HL60_result$logpval_enhancer_pval.empirical)))))
  cat("\n")
  cat(sprintf(as.character(summary(HL60_result$logpval_enhancer_pval.empirical))))
  cat("\n")
  cat(sprintf(as.character(names(summary(HL60_result$enhancer_logpval)))))
  cat("\n")
  cat(sprintf(as.character(summary(HL60_result$enhancer_logpval))))
  cat("\n")
  cat(sprintf(as.character(names(summary(HL60_result$fdr)))))
  cat("\n")
  cat(sprintf(as.character(summary(HL60_result$fdr))))
  cat("\n")
  cat(sprintf(as.character(names(summary(HL60_result$ASE_logpval)))))
  cat("\n")
  cat(sprintf(as.character(summary(HL60_result$ASE_logpval))))
  cat("\n")


  HL60_result$logpval_enhancer_pval.empirical[!is.finite(HL60_result$logpval_enhancer_pval.empirical)]<-90
  HL60_result$enhancer_logpval[!is.finite(HL60_result$enhancer_logpval)]<-300
  HL60_result$ASE_logpval[!is.finite(HL60_result$ASE_logpval)]<-90


  cat("HL60_result_1\n")
  cat(str(HL60_result))
  cat("\n")
  cat(sprintf(as.character(names(summary(HL60_result$logpval_enhancer_pval.empirical)))))
  cat("\n")
  cat(sprintf(as.character(summary(HL60_result$logpval_enhancer_pval.empirical))))
  cat("\n")
  cat(sprintf(as.character(names(summary(HL60_result$enhancer_logpval)))))
  cat("\n")
  cat(sprintf(as.character(summary(HL60_result$enhancer_logpval))))
  cat("\n")
  cat(sprintf(as.character(names(summary(HL60_result$fdr)))))
  cat("\n")
  cat(sprintf(as.character(summary(HL60_result$fdr))))
  cat("\n")
  cat(sprintf(as.character(names(summary(HL60_result$ASE_logpval)))))
  cat("\n")
  cat(sprintf(as.character(summary(HL60_result$ASE_logpval))))
  cat("\n")




  #### enhancer_result_THP1 ----

  enhancer_result_THP1<-as.data.frame(fread(file=opt$enhancer_result_THP1, sep="\t", header = T, fill=TRUE), stringsAsFactors = F)


  colnames(enhancer_result_THP1)<-c("REAL_TILE_Plus_carried_variants","statistic","control","zscore","mad.score","pval.empirical","pval.mad","pval.zscore")


  enhancer_result_THP1$logpval_enhancer_pval.empirical<--1*log10(enhancer_result_THP1$pval.empirical)
  enhancer_result_THP1$enhancer_logpval<--1*log10(enhancer_result_THP1$pval.mad)


  cat("enhancer_result_THP1\n")
  cat(str(enhancer_result_THP1))
  cat("\n")

  enhancer_result_THP1_subset<-unique(enhancer_result_THP1[,c(which(colnames(enhancer_result_THP1) == "REAL_TILE_Plus_carried_variants"),
                                                              which(colnames(enhancer_result_THP1) == "logpval_enhancer_pval.empirical"),
                                                              which(colnames(enhancer_result_THP1) == "enhancer_logpval"))])


  cat("enhancer_result_THP1_subset\n")
  cat(str(enhancer_result_THP1_subset))
  cat("\n")

  #### ASE_result_THP1 ----

  ASE_result_THP1<-as.data.frame(fread(file=opt$ASE_result_THP1, sep="\t", header = T, fill=TRUE), stringsAsFactors = F)

  colnames(ASE_result_THP1)<-c("REAL_TILE_Plus_carried_variants","statistic","pval","fdr","df.test","df.dna","df.rna.full","df.rna.red","logFC")
  ASE_result_THP1$ASE_logpval<--1*log10(ASE_result_THP1$pval)


  cat("ASE_result_THP1\n")
  cat(str(ASE_result_THP1))
  cat("\n")

  ASE_result_THP1_subset<-unique(ASE_result_THP1[,c(which(colnames(ASE_result_THP1) == "REAL_TILE_Plus_carried_variants"),
                                                    which(colnames(ASE_result_THP1) == "fdr"),
                                                    which(colnames(ASE_result_THP1) == "ASE_logpval"))])


  cat("ASE_result_THP1_subset\n")
  cat(str(ASE_result_THP1_subset))
  cat("\n")

  THP1_result<-merge(enhancer_result_THP1_subset,
                     ASE_result_THP1_subset,
                     by="REAL_TILE_Plus_carried_variants")


  THP1_result$Cell_Type<-"THP1"


  cat("THP1_result_0\n")
  cat(str(THP1_result))
  cat("\n")
  cat(sprintf(as.character(names(summary(THP1_result$logpval_enhancer_pval.empirical)))))
  cat("\n")
  cat(sprintf(as.character(summary(THP1_result$logpval_enhancer_pval.empirical))))
  cat("\n")
  cat(sprintf(as.character(names(summary(THP1_result$enhancer_logpval)))))
  cat("\n")
  cat(sprintf(as.character(summary(THP1_result$enhancer_logpval))))
  cat("\n")
  cat(sprintf(as.character(names(summary(THP1_result$fdr)))))
  cat("\n")
  cat(sprintf(as.character(summary(THP1_result$fdr))))
  cat("\n")
  cat(sprintf(as.character(names(summary(THP1_result$ASE_logpval)))))
  cat("\n")
  cat(sprintf(as.character(summary(THP1_result$ASE_logpval))))
  cat("\n")


  THP1_result$logpval_enhancer_pval.empirical[!is.finite(THP1_result$logpval_enhancer_pval.empirical)]<-90
  THP1_result$enhancer_logpval[!is.finite(THP1_result$enhancer_logpval)]<-300
  THP1_result$ASE_logpval[!is.finite(THP1_result$ASE_logpval)]<-90


  cat("THP1_result_1\n")
  cat(str(THP1_result))
  cat("\n")
  cat(sprintf(as.character(names(summary(THP1_result$logpval_enhancer_pval.empirical)))))
  cat("\n")
  cat(sprintf(as.character(summary(THP1_result$logpval_enhancer_pval.empirical))))
  cat("\n")
  cat(sprintf(as.character(names(summary(THP1_result$enhancer_logpval)))))
  cat("\n")
  cat(sprintf(as.character(summary(THP1_result$enhancer_logpval))))
  cat("\n")
  cat(sprintf(as.character(names(summary(THP1_result$fdr)))))
  cat("\n")
  cat(sprintf(as.character(summary(THP1_result$fdr))))
  cat("\n")
  cat(sprintf(as.character(names(summary(THP1_result$ASE_logpval)))))
  cat("\n")
  cat(sprintf(as.character(summary(THP1_result$ASE_logpval))))
  cat("\n")




  #### READ data ----
  
  
  setwd(out2)
  
  filename_1<-paste("list_medians_post_QC",".rds", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  medians_object<-readRDS(file=filename_1)
  
  cat("medians_object_\n")
  str(medians_object)
  cat("\n")
  
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
  Rosetta_df <-readRDS(file=filename)
  
  
  cat("Rosetta_df_0\n")
  cat(str(Rosetta_df))
  cat("\n")
  cat(str(unique(Rosetta_df$VAR)))
  cat("\n")
  
  
 
    
  Rosetta_df$Tile<-"NA"
  Rosetta_df$Tile[which(Rosetta_df$TILE == "TILE_1")]<-"NINE"
  Rosetta_df$Tile[which(Rosetta_df$TILE == "TILE_2")]<-"TWO_THIRDS"
  Rosetta_df$Tile[which(Rosetta_df$TILE == "TILE_3")]<-"HALF"
  Rosetta_df$Tile[which(Rosetta_df$TILE == "TILE_4")]<-"ONE_THIRD"
  Rosetta_df$Tile[which(Rosetta_df$TILE == "TILE_5")]<-"A_TENTH"
  
  Rosetta_df$Tile<-factor(Rosetta_df$Tile,
                             levels=c("NINE",
                                      "TWO_THIRDS","HALF","ONE_THIRD","A_TENTH"),
                             ordered=T)
  

  
  Rosetta_df$Label_2<-"NA"
  Rosetta_df$Label_2[which(Rosetta_df$Label == "downstream_gene_variant")]<-"ASSAYED_VARIANT"
  Rosetta_df$Label_2[which(Rosetta_df$Label == "intergenic_variant")]<-"ASSAYED_VARIANT"
  Rosetta_df$Label_2[which(Rosetta_df$Label == "intron_variant")]<-"ASSAYED_VARIANT"
  Rosetta_df$Label_2[which(Rosetta_df$Label == "Kousik_variant")]<-"Kousik_variant"
  Rosetta_df$Label_2[which(Rosetta_df$Label == "Negative_Control_Genomic_Regions")]<-"NCGR"
  Rosetta_df$Label_2[which(Rosetta_df$Label == "POSITIVE_CTRL")]<-"ASE_CTRL"
  Rosetta_df$Label_2[which(Rosetta_df$Label == "regulatory_region_variant")]<-"ASSAYED_VARIANT"
  Rosetta_df$Label_2[which(Rosetta_df$Label == "TF_binding_site_variant")]<-"ASSAYED_VARIANT"
  Rosetta_df$Label_2[which(Rosetta_df$Label == "UNDETERMINED_CTRL")]<-"Enhancer_CTRL"
  Rosetta_df$Label_2[which(Rosetta_df$Label == "upstream_gene_variant")]<-"ASSAYED_VARIANT"
 
  
  Rosetta_df$Label_2<-factor(Rosetta_df$Label_2,
                          levels=c("ASSAYED_VARIANT",
                                   "NCGR","ASE_CTRL","Enhancer_CTRL","Kousik_variant"),
                          ordered=T)
  
  Rosetta_df$Label<-factor(Rosetta_df$Label,
                             levels=c("intron_variant","intergenic_variant","upstream_gene_variant","downstream_gene_variant",
                                      "regulatory_region_variant","TF_binding_site_variant",
                                      "Negative_Control_Genomic_Regions","POSITIVE_CTRL","UNDETERMINED_CTRL","Kousik_variant"),
                             ordered=T)
  
  
  cat("Rosetta_df_1\n")
  cat(str(Rosetta_df))
  cat("\n")
  cat(sprintf(as.character(names(summary(Rosetta_df$Label)))))
  cat("\n")
  cat(sprintf(as.character(summary(Rosetta_df$Label))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Rosetta_df$Label_2)))))
  cat("\n")
  cat(sprintf(as.character(summary(Rosetta_df$Label_2))))
  cat("\n")
  
  # intron_variant intergenic_variant upstream_gene_variant downstream_gene_variant regulatory_region_variant TF_binding_site_variant Negative_Control_Genomic_Regions POSITIVE_CTRL UNDETERMINED_CTRL Kousik_variant
  # 409 74 15 25 20 18 10 37 60 55
  # ASSAYED_VARIANT Negative_Control_Genomic_Regions ASE_CTRL UNDETERMINED_CTRL Kousik_variant
  # 561 10 37 60 55
  # NP LIBRARY ------------------> ASE_CTRL ASSAYED_VARIANT Negative_Control_Genomic_Regions UNDETERMINED_CTRL
  # NP LIBRARY ------------------> 52 11499 11 59
  
  
  REAL_TILE_DF_NEG_CTRLS<-Rosetta_df[which(Rosetta_df$Label == 'Negative_Control_Genomic_Regions'),]
  
  cat("REAL_TILE_DF_NEG_CTRLS\n")
  cat(str(REAL_TILE_DF_NEG_CTRLS))
  cat("\n")
  
  
  cat(sprintf(as.character(names(summary(as.factor(REAL_TILE_DF_NEG_CTRLS$KEY))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(REAL_TILE_DF_NEG_CTRLS$KEY)))))
  cat("\n")
  
  

  
  # ##############################################################
  # quit(status = 1)
  
  #### undo the list median object ----
  
  
  names_medians_object<-names(medians_object)
  
  # cat("names_medians_object\n")
  # cat(str(names_medians_object))
  # cat("\n")
  
  temp<-data.frame(matrix(vector(), 0, 4,
                          dimnames=list(c(),
                                        c("REAL_TILE_Plus_carried_variants","Cell_Type","value","variable"))
                                        ))
  
  
  for(i in 1:length(names_medians_object))
  {
    names_medians_object_sel<-names_medians_object[i]
    
    cat("names_medians_object_sel\n")
    cat(str(names_medians_object_sel))
    cat("\n")
    
    
    medians_object_sel<-medians_object[[names_medians_object_sel]]
    
    # cat("medians_object_sel\n")
    # cat(str(medians_object_sel))
    # cat("\n")
    
    colnames(medians_object_sel)[which(colnames(medians_object_sel) == "median_Median_value")]<-"value"
    colnames(medians_object_sel)[which(colnames(medians_object_sel) == "median_LogFC")]<-"value"
    colnames(medians_object_sel)[which(colnames(medians_object_sel) == "median_ASE")]<-"value"
    colnames(medians_object_sel)[which(colnames(medians_object_sel) == "median_dna.REF.vs.dna.ALT")]<-"value"
    colnames(medians_object_sel)[which(colnames(medians_object_sel) == "median_rna.REF.vs.rna.ALT")]<-"value"
    
    
    
    
    
    
    
    medians_object_sel$variable<-names_medians_object_sel
    
    temp<-rbind(temp,medians_object_sel)
    
  }# i 
  
  
  cat("temp_0\n")
  cat(str(temp))
  cat("\n")
 
  
  #### merge MPRA_Real_Tile & subset ----
  
  MPRA_Real_Tile<-merge(Rosetta_df,
                        temp,
                        by=c("REAL_TILE_Plus_carried_variants"), stringsAsFactors=F)
  
  cat("MPRA_Real_Tile_0\n")
  cat(str(MPRA_Real_Tile))
  cat("\n")
  
  MPRA_Real_Tile$chr<-gsub("_.+$","",MPRA_Real_Tile$VAR)
  
  MPRA_Real_Tile$chr<-factor(MPRA_Real_Tile$chr,
                             levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                                      "chr17","chr18","chr19","chr20","chr21","chr22","chr23","chrX","chrY"),
                             ordered=T)
  MPRA_Real_Tile$pos<-gsub("^[^_]+_","",MPRA_Real_Tile$VAR)
  MPRA_Real_Tile$pos<-as.integer(gsub("_.+$","",MPRA_Real_Tile$pos))
  
  
  
  cat("MPRA_Real_Tile_1\n")
  cat(str(MPRA_Real_Tile))
  cat("\n")
  
  #### RMV Kousik ----
  
  MPRA_Real_Tile<-MPRA_Real_Tile[which(MPRA_Real_Tile$Label != "Kousik_variant"),]
  
  cat("Rosetta_RMV_Kousik\n")
  cat(str(MPRA_Real_Tile))
  cat("\n")
 
  
  # quit(status=1)
  
 
 
  
  
  MPRA_Real_Tile_subset<-unique(MPRA_Real_Tile[,c(which(colnames(MPRA_Real_Tile) == "REAL_TILE_Plus_carried_variants"),
                                           which(colnames(MPRA_Real_Tile) == "Cell_Type"),
                                           which(colnames(MPRA_Real_Tile) == "factor4"),
                                           which(colnames(MPRA_Real_Tile) == "Label"),
                                           which(colnames(MPRA_Real_Tile) == "Label_2"),
                                           which(colnames(MPRA_Real_Tile) == "KEY"),
                                           which(colnames(MPRA_Real_Tile) == "TILE"),
                                           which(colnames(MPRA_Real_Tile) == "Tile"),
                                           which(colnames(MPRA_Real_Tile) == "VAR"),
                                           which(colnames(MPRA_Real_Tile) == "carried_variants"),
                                           which(colnames(MPRA_Real_Tile) == "chr"),
                                           which(colnames(MPRA_Real_Tile) == "pos"),
                                           which(colnames(MPRA_Real_Tile) == "variable"),
                                           which(colnames(MPRA_Real_Tile) == "value"))])
  
  cat("MPRA_Real_Tile_subset\n")
  cat(str(MPRA_Real_Tile_subset))
  cat("\n")
  
  
  
  
  ##### merge with enhancer and ASE ----
  
  
  RESULTS<-rbind(K562_result,CHRF_result,HL60_result,THP1_result)
  
  cat("RESULTS\n")
  cat(str(RESULTS))
  cat("\n")
  
  
  ####  SAVE ----
  
  
  setwd(out2)
  
  filename_1<-paste("MPRA_Real_Tile",".rds", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  saveRDS(MPRA_Real_Tile_subset,file=filename_1)
  
  filename_1<-paste("RESULTS_MPRAnalyze",".rds", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  saveRDS(RESULTS,file=filename_1)
  
  
 
  
}

QC2_barplots_plots = function(option_list)
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
  
  #### READ and transform out2 ----
  
  out2 = opt$out2
  
  cat("out2_\n")
  cat(sprintf(as.character(out2)))
  cat("\n")
  
  #### READ and transform fdr_Threshold ----
  
  fdr_Threshold = opt$fdr_Threshold
  
  cat("fdr_Threshold_\n")
  cat(sprintf(as.character(fdr_Threshold)))
  cat("\n")
  
  
  #### read RDS----
  
  setwd(out2)
  
  filename_1<-paste("MPRA_Real_Tile",".rds", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  MPRA_Real_Tile<-readRDS(file=filename_1)
  
  cat("MPRA_Real_Tile_0\n")
  cat(str(MPRA_Real_Tile))
  cat("\n")
  
  
  filename_1<-paste("RESULTS_MPRAnalyze",".rds", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  RESULTS<-readRDS(file=filename_1)
  
  cat("RESULTS\n")
  cat(str(RESULTS))
  cat("\n")
  
  #### path2 ----
  
  path2<-paste(out2,'QC2_Barplots','/', sep='')
  
  cat("path2\n")
  cat(sprintf(as.character(path2)))
  cat("\n")
  
  if (file.exists(path2)){
    
    unlink(path2)
    dir.create(file.path(path2))
    
    
  } else {
    dir.create(file.path(path2))
    
  }
  
  
  
  ##### LOOP graphs ----
  
    
  Cell_Types_vector<-levels(MPRA_Real_Tile$Cell_Type)
  
  cat("Cell_Types_vector_3\n")
  cat(str(Cell_Types_vector))
  cat("\n")
  
  for(i in 1:length(Cell_Types_vector))
  {
    Cell_Types_vector_sel<-Cell_Types_vector[i]
    
    cat("----------------------------------------------------------------------------------------------------------------------------------------------------->\t")
    cat(sprintf(as.character(Cell_Types_vector_sel)))
    cat("\n")
    
    
    MPRA_Real_Tile_sel<-MPRA_Real_Tile[which(MPRA_Real_Tile$Cell_Type == Cell_Types_vector_sel),]
    
    cat("MPRA_Real_Tile_sel\n")
    cat(str(MPRA_Real_Tile_sel))
    cat("\n")
    
    
    
    # #### path3 ----
    # 
    # path3<-paste(out,'QC2_Barplots','/',Cell_Types_vector_sel,'/', sep='')
    # 
    # cat("path3\n")
    # cat(sprintf(as.character(path3)))
    # cat("\n")
    # 
    # if (file.exists(path3)){
    #   
    #   
    #   
    #   
    # } else {
    #   dir.create(file.path(path3))
    #   
    # }
    
    
    #### select features that are part of QC2 ----
    
    
    features_sel<-c("dna_pool","rna_pool","dna.REF.vs.dna.ALT")
    
    
    REP<-MPRA_Real_Tile_sel[which(MPRA_Real_Tile_sel$variable%in%features_sel),]
    
    
    cat("REP_0\n")
    cat(str(REP))
    cat("\n")
    
    #### Discretize select features that are part of QC2 ----
    
    
    REP$variable_CLASS<-"NA"
    
    REP$variable_CLASS[which(REP$value >0)]<-"> 0"
    REP$variable_CLASS[which(REP$value == 0)]<-"0"
    
    REP$variable_CLASS<-factor(REP$variable_CLASS,
                                                           levels=c("0","> 0"),
                                                           ordered = T)
    
    cat("REP_1\n")
    cat(str(REP))
    cat("\n")
    cat(sprintf(as.character(names(summary(REP$variable_CLASS)))))
    cat("\n")
    cat(sprintf(as.character(summary(REP$variable_CLASS))))
    cat("\n")
    
    breaks.Rank<-seq(0,110,by=10)
    labels.Rank<-as.character(breaks.Rank)
    
    cat(sprintf(as.character(labels.Rank)))
    cat("\n")
    
    
    #### Freq.table_TOTAL
    
    REP.dt<-data.table(REP, key=c("Cell_Type","variable","variable_CLASS"))
    
    
    
    Freq.table_TOTAL<-as.data.frame(REP.dt[,.N,by=key(REP.dt)], stringsAsFactors=F)
    
    colnames(Freq.table_TOTAL)[which(colnames(Freq.table_TOTAL) == "N")]<-"TOTAL"
    
    cat("Freq.table_TOTAL_\n")
    cat(str(Freq.table_TOTAL))
    cat("\n")
    
    setwd(out2)
    
    write.table(Freq.table_TOTAL, file="TOTAL.txt",sep="\t",quote=F, row.names = F)
    
    #### Freq.table_pool
    
    REP.dt<-data.table(REP, key=c("Cell_Type","variable","variable_CLASS","factor4"))
    
    
    
    Freq.table_pool<-as.data.frame(REP.dt[,.N,by=key(REP.dt)], stringsAsFactors=F)
    
    colnames(Freq.table_pool)[which(colnames(Freq.table_pool) == "N")]<-"instances"
    
    cat("Freq.table_pool_0\n")
    cat(str(Freq.table_pool))
    cat("\n")
    
    setwd(out2)
    
    write.table(Freq.table_pool, file="Pool.txt",sep="\t",quote=F, row.names = F)
    
    #### Freq.table_Label
    
    REP.dt<-data.table(REP, key=c("Cell_Type","variable","variable_CLASS","Label"))
    
    
    
    Freq.table_Label<-as.data.frame(REP.dt[,.N,by=key(REP.dt)], stringsAsFactors=F)
    
    colnames(Freq.table_Label)[which(colnames(Freq.table_Label) == "N")]<-"instances"
    
    cat("Freq.table_Label_0\n")
    cat(str(Freq.table_Label))
    cat("\n")
    
    #### Barplot # pool
    
    Freq.table_pool<-merge(Freq.table_pool,
                           Freq.table_TOTAL,
                           by=c("Cell_Type","variable","variable_CLASS"))
    
    Freq.table_pool$Perc<-round(100*(Freq.table_pool$instances/Freq.table_pool$TOTAL), 1)
    
    
    cat("Freq.table_pool_1\n")
    cat(str(Freq.table_pool))
    cat("\n")
    
    setwd(out2)
    
    write.table(Freq.table_pool, file="Pool_DEF.txt",sep="\t",quote=F, row.names = F)
    
   
    
    QC_2_plot<-Freq.table_pool %>%
      mutate(myaxis = paste0(Freq.table_pool$variable_CLASS, "\n", "n=", TOTAL), drop=F) %>%
      mutate(myaxis=fct_reorder(myaxis,as.numeric(Freq.table_pool$variable_CLASS)), drop=F) %>%
      ggplot(aes(x=myaxis, y=Perc, fill=factor4)) +
      geom_bar(stat="identity",colour='black')+
      theme_bw()+
      theme(axis.title.y=element_text(size=24, family="sans"),
            axis.title.x=element_text(size=24, family="sans"),
            axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
            axis.text.x=element_text(angle=45,size=12, color="black", family="sans"),
            legend.title=element_text(size=16,color="black", family="sans"),
            legend.text=element_text(size=12,color="black", family="sans"))+
      scale_x_discrete(name=NULL, drop=F)+
      scale_y_continuous(name="Percentage",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
      scale_fill_manual(values=c('red','#6DB2EE','#62D07F','#1877C9','#C9244B','#D45E85','#87447B','#553B68','#D6B8E6',"#89C5DA",
                                  "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C","#D7C1B1", "#AD6F3B", "#689030", "red"), drop=F)+
      theme(legend.position="bottom")+
      ggeasy::easy_center_title()
    
      
    QC_2_plot <- QC_2_plot + facet_grid(. ~ variable, drop=F)
      
    setwd(path2)
    
    svgname<-paste("QC2_Barplot_pool_",Cell_Types_vector_sel,".svg",sep='')
   
    makesvg = TRUE
    
    if (makesvg == TRUE)
    {
      ggsave(svgname, plot= QC_2_plot,
             device="svg",
             height=10, width=12)
    }
    
    #### Barplot # Label
    
    Freq.table_Label<-merge(Freq.table_Label,
                           Freq.table_TOTAL,
                           by=c("Cell_Type","variable","variable_CLASS"))
    
    Freq.table_Label$Perc<-round(100*(Freq.table_Label$instances/Freq.table_Label$TOTAL), 1)
    
    
    cat("Freq.table_Label_1\n")
    cat(str(Freq.table_Label))
    cat("\n")
    
    setwd(out2)
    
    write.table(Freq.table_Label, file="Label_DEF.txt",sep="\t",quote=F, row.names = F)
    
    
    
    QC_2_plot<-Freq.table_Label %>%
      mutate(myaxis = paste0(Freq.table_Label$variable_CLASS, "\n", "n=", TOTAL), drop=F) %>%
      mutate(myaxis=fct_reorder(myaxis,as.numeric(Freq.table_Label$variable_CLASS)), drop=F) %>%
      ggplot(aes(x=myaxis, y=Perc, fill=Label)) +
      geom_bar(stat="identity",colour='black')+
      theme_bw()+
      theme(axis.title.y=element_text(size=24, family="sans"),
            axis.title.x=element_text(size=24, family="sans"),
            axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
            axis.text.x=element_text(angle=45,size=12, color="black", family="sans"),
            legend.title=element_text(size=16,color="black", family="sans"),
            legend.text=element_text(size=12,color="black", family="sans"))+
      scale_x_discrete(name=NULL, drop=F)+
      scale_y_continuous(name="Percentage",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
      scale_fill_manual(values=c('red','#6DB2EE','#62D07F','#1877C9','#C9244B','#D45E85','#87447B','#553B68','#D6B8E6',"#89C5DA",
                                 "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C","#D7C1B1", "#AD6F3B", "#689030", "red"), drop=F)+
      theme(legend.position="bottom")+
      ggeasy::easy_center_title()
    
    
    QC_2_plot <- QC_2_plot + facet_grid(. ~ variable, drop=F)
    
    setwd(path2)
    
    svgname<-paste("QC2_Barplot_Label_",Cell_Types_vector_sel,".svg",sep='')
    
    makesvg = TRUE
    
    if (makesvg == TRUE)
    {
      ggsave(svgname, plot= QC_2_plot,
             device="svg",
             height=10, width=12)
    }
   
  }#i
  
  # ############################################################################################################################################################################################################################################
  # quit(status=1)
}

data_wrangling_after_QC2_plots = function(option_list)
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
  
  #### READ and transform out2 ----
  
  out2 = opt$out2
  
  cat("out2_\n")
  cat(sprintf(as.character(out2)))
  cat("\n")
  
  #### READ and transform fdr_Threshold ----
  
  fdr_Threshold = opt$fdr_Threshold
  
  cat("fdr_Threshold_\n")
  cat(sprintf(as.character(fdr_Threshold)))
  cat("\n")
  
  #### READ and transform enhancer_logval_Threshold ----
  
  enhancer_logval_Threshold = opt$enhancer_logval_Threshold
  
  cat("enhancer_logval_Threshold_\n")
  cat(sprintf(as.character(enhancer_logval_Threshold)))
  cat("\n")
  
  #### READ and transform FC_Threshold ----
  
  FC_Threshold = opt$FC_Threshold
  
  cat("FC_Threshold_\n")
  cat(sprintf(as.character(FC_Threshold)))
  cat("\n")
  
 
  
  #### READ and transform ASE_log_pval_Threshold ----
  
  ASE_log_pval_Threshold = opt$ASE_log_pval_Threshold
  
  cat("ASE_log_pval_Threshold_\n")
  cat(sprintf(as.character(ASE_log_pval_Threshold)))
  cat("\n")
  
  #### READ and transform ASE_Threshold ----
  
  ASE_Threshold = as.numeric(unlist(strsplit(opt$ASE_Threshold, split=",")))
  
  cat("ASE_Threshold_\n")
  cat(sprintf(as.character(ASE_Threshold)))
  cat("\n")
  
  ASE_LOW<-ASE_Threshold[1]
  
  cat("ASE_LOW_\n")
  cat(sprintf(as.character(ASE_LOW)))
  cat("\n")
  
  ASE_HIGH<-ASE_Threshold[2]
  
  cat("ASE_HIGH_\n")
  cat(sprintf(as.character(ASE_HIGH)))
  cat("\n")
  
  
  
  #### read RDS----
  
  setwd(out2)
  
  filename_1<-paste("MPRA_Real_Tile",".rds", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  MPRA_Real_Tile<-readRDS(file=filename_1)
  
  cat("MPRA_Real_Tile_0\n")
  cat(str(MPRA_Real_Tile))
  cat("\n")
  
  
  filename_1<-paste("RESULTS_MPRAnalyze",".rds", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  RESULTS<-readRDS(file=filename_1)
  
  cat("RESULTS\n")
  cat(str(RESULTS))
  cat("\n")
  
  #### merge MPRA_Real_Tile and RESULTS ----
  
  MPRA_Real_Tile<-merge(MPRA_Real_Tile,
                        RESULTS,
                        by=c("REAL_TILE_Plus_carried_variants","Cell_Type"),
                        all=T)
  
  cat("MPRA_Real_Tile_1\n")
  cat(str(MPRA_Real_Tile))
  cat("\n")
  cat("variable\n")
  cat(sprintf(as.character(names(summary(as.factor(MPRA_Real_Tile$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(MPRA_Real_Tile$variable)))))
  cat("\n")
  
  
  ##### QC2 PASS ----
  
  # QC2_PASS_CLASS
  
 
    
    
  indx.int<-c(which(colnames(MPRA_Real_Tile) == "REAL_TILE_Plus_carried_variants"),
              which(colnames(MPRA_Real_Tile) == "Cell_Type"),
              which(colnames(MPRA_Real_Tile) == "carried_variants"),
              which(colnames(MPRA_Real_Tile) == "factor4"),
              which(colnames(MPRA_Real_Tile) == "Label"),
              which(colnames(MPRA_Real_Tile) == "Label_2"),
              which(colnames(MPRA_Real_Tile) == "KEY"),
              which(colnames(MPRA_Real_Tile) == "TILE"),
              which(colnames(MPRA_Real_Tile) == "Tile"),
              which(colnames(MPRA_Real_Tile) == "VAR"),
              which(colnames(MPRA_Real_Tile) == "chr"),
              which(colnames(MPRA_Real_Tile) == "pos"),
              which(colnames(MPRA_Real_Tile) == "logpval_enhancer_pval.empirical"),
              which(colnames(MPRA_Real_Tile) == "enhancer_logpval"),
              which(colnames(MPRA_Real_Tile) == "fdr"),
              which(colnames(MPRA_Real_Tile) == "ASE_logpval"))
  
  
  
  MPRA_Real_Tile_QC2_PASS<-MPRA_Real_Tile[which(MPRA_Real_Tile$variable == "dna_pool" |
                                             MPRA_Real_Tile$variable == "rna_pool" |
                                             MPRA_Real_Tile$variable == "dna.REF.vs.dna.ALT"),]
  
  cat("MPRA_Real_Tile_QC2_PASS\n")
  cat(str(MPRA_Real_Tile_QC2_PASS))
  cat("\n")
  
  MPRA_Real_Tile_QC2_PASS_wide<-as.data.frame(pivot_wider(MPRA_Real_Tile_QC2_PASS,
                                                     id_cols=colnames(MPRA_Real_Tile)[indx.int],
                                                     names_from=variable,
                                                     values_from=value),stringsAsFactors=F)
  
  cat("MPRA_Real_Tile_QC2_PASS_wide_0\n")
  cat(str(MPRA_Real_Tile_QC2_PASS_wide))
  cat("\n")
  
  
  
  
  MPRA_Real_Tile_QC2_PASS_wide$QC2_CLASS<-"NA"
  
  
  
  indx.int<-c(which(colnames(MPRA_Real_Tile_QC2_PASS_wide) == "rna_pool"),
              which(colnames(MPRA_Real_Tile_QC2_PASS_wide) == 'dna.REF.vs.dna.ALT'),
              which(colnames(MPRA_Real_Tile_QC2_PASS_wide) == 'dna_pool'))
  
  
  MPRA_Real_Tile_QC2_PASS_wide$dna_pool_CLASS<-"NA"
  
  MPRA_Real_Tile_QC2_PASS_wide$dna_pool_CLASS[MPRA_Real_Tile_QC2_PASS_wide[,indx.int[3]]>0]<-"PASS"
  
  MPRA_Real_Tile_QC2_PASS_wide$dna_ratio_CLASS<-"NA"
  
  MPRA_Real_Tile_QC2_PASS_wide$dna_ratio_CLASS[MPRA_Real_Tile_QC2_PASS_wide[,indx.int[2]]<= 5]<-"PASS"
  MPRA_Real_Tile_QC2_PASS_wide$dna_ratio_CLASS[MPRA_Real_Tile_QC2_PASS_wide[,indx.int[2]]<= 0.2]<-"NA"
  
  
  MPRA_Real_Tile_QC2_PASS_wide$QC2_CLASS[which(MPRA_Real_Tile_QC2_PASS_wide$dna_pool_CLASS == "PASS" &
                                                 MPRA_Real_Tile_QC2_PASS_wide$dna_ratio_CLASS == "PASS")]<-"PASS"
  
  
  cat("MPRA_Real_Tile_QC2_PASS_wide_1\n")
  cat(str(MPRA_Real_Tile_QC2_PASS_wide))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(MPRA_Real_Tile_QC2_PASS_wide$dna_pool_CLASS))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(MPRA_Real_Tile_QC2_PASS_wide$dna_pool_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(MPRA_Real_Tile_QC2_PASS_wide$dna_ratio_CLASS))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(MPRA_Real_Tile_QC2_PASS_wide$dna_ratio_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(MPRA_Real_Tile_QC2_PASS_wide$QC2_CLASS))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(MPRA_Real_Tile_QC2_PASS_wide$QC2_CLASS)))))
  cat("\n")
  
  
  MPRA_Real_Tile_QC2_PASS_wide_subset<-MPRA_Real_Tile_QC2_PASS_wide[,c(which(colnames(MPRA_Real_Tile_QC2_PASS_wide) == "REAL_TILE_Plus_carried_variants"),
                                                                       which(colnames(MPRA_Real_Tile_QC2_PASS_wide) == "Cell_Type"),
                                                                      which(colnames(MPRA_Real_Tile_QC2_PASS_wide) == "QC2_CLASS"))]
  
  cat("MPRA_Real_Tile_QC2_PASS_wide_subset\n")
  cat(str(MPRA_Real_Tile_QC2_PASS_wide_subset))
  cat("\n")
  
  ########################################## NO FILTER 
  
 
  MPRA_Real_Tile_NO_FILTER<-merge(MPRA_Real_Tile,
                                  MPRA_Real_Tile_QC2_PASS_wide_subset,
                                  by=c("REAL_TILE_Plus_carried_variants","Cell_Type"),
                                  all=T)
  
  cat("MPRA_Real_Tile_NO_FILTER_0\n")
  cat(str(MPRA_Real_Tile_NO_FILTER))
  cat("\n")
  
  
  
  ########################################## FILTER 
  
  
  
  
  
  MPRA_Real_Tile_QC2_PASS_wide_subset_PASS<-MPRA_Real_Tile_QC2_PASS_wide_subset[which(MPRA_Real_Tile_QC2_PASS_wide_subset$QC2_CLASS == "PASS"),]
  
  cat("MPRA_Real_Tile_QC2_PASS_wide_subset_PASS\n")
  cat(str(MPRA_Real_Tile_QC2_PASS_wide_subset_PASS))
  cat("\n")
  
  cat("MPRA_Real_Tile_PRE\n")
  cat(str(MPRA_Real_Tile))
  cat("\n")
  
  MPRA_Real_Tile<-merge(MPRA_Real_Tile,
                        MPRA_Real_Tile_QC2_PASS_wide_subset_PASS,
                        by=c("REAL_TILE_Plus_carried_variants","Cell_Type"),
                        all.y=T)
  
  cat("MPRA_Real_Tile_2\n")
  cat(str(MPRA_Real_Tile))
  cat("\n")
  

  
  #### RESULTS_SIG ----

  RESULTS_SIG_enhancer<-MPRA_Real_Tile[which(MPRA_Real_Tile$variable == "LogFC"),]
  
  cat("RESULTS_SIG_enhancer_1\n")
  cat(str(RESULTS_SIG_enhancer))
  cat("\n")
  
  RESULTS_SIG_enhancer<-RESULTS_SIG_enhancer[which(RESULTS_SIG_enhancer$enhancer_logpval >= enhancer_logval_Threshold),]
  
  
  
  cat("RESULTS_SIG_enhancer_2\n")
  cat(str(RESULTS_SIG_enhancer))
  cat("\n")
  
  RESULTS_SIG_enhancer<-RESULTS_SIG_enhancer[which(RESULTS_SIG_enhancer$value > FC_Threshold),]
  
  
  RESULTS_SIG_enhancer$CLASS_enhancer<-"enhancer"
  
  
  cat("RESULTS_SIG_enhancer_3\n")
  cat(str(RESULTS_SIG_enhancer))
  cat("\n")
  
  cat("distribution_SIG_enhancer_LogFC\n")
  cat(sprintf(as.character(names(summary(RESULTS_SIG_enhancer$value)))))
  cat("\n")
  cat(sprintf(as.character(summary(RESULTS_SIG_enhancer$value))))
  cat("\n")
  cat("distribution_SIG_enhancer_enhancer_logpval\n")
  cat(sprintf(as.character(names(summary(RESULTS_SIG_enhancer$enhancer_logpval)))))
  cat("\n")
  cat(sprintf(as.character(summary(RESULTS_SIG_enhancer$enhancer_logpval))))
  cat("\n")
  cat("distribution_SIG_enhancer_pool\n")
  cat(sprintf(as.character(names(summary(as.factor(RESULTS_SIG_enhancer$factor4))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(RESULTS_SIG_enhancer$factor4)))))
  cat("\n")
  cat("distribution_SIG_enhancer_Label\n")
  cat(sprintf(as.character(names(summary(as.factor(RESULTS_SIG_enhancer$Label))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(RESULTS_SIG_enhancer$Label)))))
  cat("\n")
  cat("distribution_SIG_enhancer_Cell_Type\n")
  cat(sprintf(as.character(names(summary(as.factor(RESULTS_SIG_enhancer$Cell_Type))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(RESULTS_SIG_enhancer$Cell_Type)))))
  cat("\n")
  
  # distribution_SIG_ASE_Cell_Type
  # K562 CHRF HL60 THP1
  # 148 153 138 147
  
  
  RESULTS_SIG_ASE<-MPRA_Real_Tile[which(MPRA_Real_Tile$variable == "ASE"),]
  
  cat("RESULTS_SIG_ASE_1\n")
  cat(str(RESULTS_SIG_ASE))
  cat("\n")
  
  RESULTS_SIG_ASE<-RESULTS_SIG_ASE[which(RESULTS_SIG_ASE$ASE_logpval >= ASE_log_pval_Threshold),]
  
  
  
  cat("RESULTS_SIG_ASE_2\n")
  cat(str(RESULTS_SIG_ASE))
  cat("\n")
  
  RESULTS_SIG_ASE<-RESULTS_SIG_ASE[which(RESULTS_SIG_ASE$fdr <= fdr_Threshold),]
  
  
  cat("RESULTS_SIG_ASE_3\n")
  cat(str(RESULTS_SIG_ASE))
  cat("\n")
  
  RESULTS_SIG_ASE_1<-RESULTS_SIG_ASE[which(RESULTS_SIG_ASE$value <= ASE_LOW),]
  
  
  RESULTS_SIG_ASE_2<-RESULTS_SIG_ASE[which(RESULTS_SIG_ASE$value >= ASE_HIGH),]
 
  RESULTS_SIG_ASE<-rbind(RESULTS_SIG_ASE_1,RESULTS_SIG_ASE_2)
  
  RESULTS_SIG_ASE$CLASS_ASE<-"ASE"
  
  
  cat("RESULTS_SIG_ASE_4\n")
  cat(str(RESULTS_SIG_ASE))
  cat("\n")
  
 
  
  
  cat("distribution_SIG_ASE_LogFC\n")
  cat(sprintf(as.character(names(summary(RESULTS_SIG_ASE$value)))))
  cat("\n")
  cat(sprintf(as.character(summary(RESULTS_SIG_ASE$value))))
  cat("\n")
  cat("distribution_SIG_ASE_ASE_logpval\n")
  cat(sprintf(as.character(names(summary(RESULTS_SIG_ASE$ASE_logpval)))))
  cat("\n")
  cat(sprintf(as.character(summary(RESULTS_SIG_ASE$ASE_logpval))))
  cat("\n")
  cat("distribution_SIG_ASE_pool\n")
  cat(sprintf(as.character(names(summary(as.factor(RESULTS_SIG_ASE$factor4))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(RESULTS_SIG_ASE$factor4)))))
  cat("\n")
  cat("distribution_SIG_ASE_Label\n")
  cat(sprintf(as.character(names(summary(as.factor(RESULTS_SIG_ASE$Label))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(RESULTS_SIG_ASE$Label)))))
  cat("\n")
  cat("distribution_SIG_ASE_Cell_Type\n")
  cat(sprintf(as.character(names(summary(as.factor(RESULTS_SIG_ASE$Cell_Type))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(RESULTS_SIG_ASE$Cell_Type)))))
  cat("\n")
  
  # distribution_SIG_ASE_Cell_Type
  # K562 CHRF HL60 THP1
  # 138 259 175 2
  
  
 
  
  #### INTERSECTION SIG ----
  
  
  
  indx.dep<-c(which(colnames(RESULTS_SIG_enhancer) == "carried_variants"),
              which(colnames(RESULTS_SIG_enhancer) == "factor4"),
              which(colnames(RESULTS_SIG_enhancer) == "Label"),
              which(colnames(RESULTS_SIG_enhancer) == "Label_2"),
              which(colnames(RESULTS_SIG_enhancer) == "KEY"),
              which(colnames(RESULTS_SIG_enhancer) == "TILE"),
              which(colnames(RESULTS_SIG_enhancer) == "Tile"),
              which(colnames(RESULTS_SIG_enhancer) == "VAR"),
              which(colnames(RESULTS_SIG_enhancer) == "chr"),
              which(colnames(RESULTS_SIG_enhancer) == "pos"))
  



  RESULTS_E_Plus_ASE_LogFC<-MPRA_Real_Tile[which(MPRA_Real_Tile$variable == "LogFC"),-c(indx.dep, which(colnames(MPRA_Real_Tile) == "variable"))]

  colnames(RESULTS_E_Plus_ASE_LogFC)[which(colnames(RESULTS_E_Plus_ASE_LogFC) == "value")]<-"LogFC"
  
  cat("RESULTS_E_Plus_ASE_LogFC\n")
  cat(str(RESULTS_E_Plus_ASE_LogFC))
  cat("\n")
  
  RESULTS_E_Plus_ASE_ASE<-MPRA_Real_Tile[which(MPRA_Real_Tile$variable == "ASE"),-c(indx.dep, which(colnames(MPRA_Real_Tile) == "variable"))]
  
  colnames(RESULTS_E_Plus_ASE_ASE)[which(colnames(RESULTS_E_Plus_ASE_ASE) == "value")]<-"ASE"
  
  cat("RESULTS_E_Plus_ASE_ASE\n")
  cat(str(RESULTS_E_Plus_ASE_ASE))
  cat("\n")
  
  
  indx.int<-c(which(colnames(RESULTS_E_Plus_ASE_ASE) == "REAL_TILE_Plus_carried_variants"),
              which(colnames(RESULTS_E_Plus_ASE_ASE) == "Cell_Type"),
              which(colnames(RESULTS_E_Plus_ASE_ASE) == "logpval_enhancer_pval.empirical"),
              which(colnames(RESULTS_E_Plus_ASE_ASE) == "enhancer_logpval"),
              which(colnames(RESULTS_E_Plus_ASE_ASE) == "fdr"),
              which(colnames(RESULTS_E_Plus_ASE_ASE) == "ASE_logpval"),
              which(colnames(RESULTS_E_Plus_ASE_ASE) == "QC2_CLASS"))
  
  
  # RESULTS_E_Plus_ASE_wide<-as.data.frame(pivot_wider(RESULTS_E_Plus_ASE,
  #                                              id_cols=colnames(RESULTS_E_Plus_ASE)[indx.int],
  #                                                 names_from=variable,
  #                                                 values_from=value),stringsAsFactors=F)
  
  RESULTS_E_Plus_ASE_wide<-merge(RESULTS_E_Plus_ASE_LogFC,
                                 RESULTS_E_Plus_ASE_ASE,
                                 by=colnames(RESULTS_E_Plus_ASE_ASE)[indx.int])

  cat("RESULTS_E_Plus_ASE_wide_0\n")
  cat(str(RESULTS_E_Plus_ASE_wide))
  cat("\n")
  
  

  RESULTS_E_Plus_ASE_wide<-RESULTS_E_Plus_ASE_wide[which(RESULTS_E_Plus_ASE_wide$enhancer_logpval >= enhancer_logval_Threshold),]
  
  cat("RESULTS_E_Plus_ASE_wide_1\n")
  cat(str(RESULTS_E_Plus_ASE_wide))
  cat("\n")
  
  RESULTS_E_Plus_ASE_wide<-RESULTS_E_Plus_ASE_wide[which(RESULTS_E_Plus_ASE_wide$LogFC > FC_Threshold),]
  
  cat("RESULTS_E_Plus_ASE_wide_2\n")
  cat(str(RESULTS_E_Plus_ASE_wide))
  cat("\n")
  
  RESULTS_E_Plus_ASE_wide<-RESULTS_E_Plus_ASE_wide[which(RESULTS_E_Plus_ASE_wide$ASE_logpval >= ASE_log_pval_Threshold),]

  cat("RESULTS_E_Plus_ASE_wide_3\n")
  cat(str(RESULTS_E_Plus_ASE_wide))
  cat("\n")
  
  
  if(dim(RESULTS_E_Plus_ASE_wide)[1] ==0)
  {
    
    cat("RESULTS_E_Plus_ASE_wide_4\n")
    cat(str(RESULTS_E_Plus_ASE_wide))
    cat("\n")
    cat("distribution_SIG_ASE_pool\n")
    cat(sprintf(as.character(names(summary(as.factor(RESULTS_E_Plus_ASE_wide$factor4))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(RESULTS_E_Plus_ASE_wide$factor4)))))
    cat("\n")
    cat("distribution_SIG_ASE_Label\n")
    cat(sprintf(as.character(names(summary(as.factor(RESULTS_E_Plus_ASE_wide$Label))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(RESULTS_E_Plus_ASE_wide$Label)))))
    cat("\n")
    cat("distribution_SIG_ASE_Cell_Type\n")
    cat(sprintf(as.character(names(summary(as.factor(RESULTS_E_Plus_ASE_wide$Cell_Type))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(RESULTS_E_Plus_ASE_wide$Cell_Type)))))
    cat("\n")
    
  }else{
    # quit(status = 1)
    
    RESULTS_E_Plus_ASE_wide_1<-RESULTS_E_Plus_ASE_wide[which(RESULTS_E_Plus_ASE_wide$ASE >= ASE_HIGH),]
    
    RESULTS_E_Plus_ASE_wide_2<-RESULTS_E_Plus_ASE_wide[which(RESULTS_E_Plus_ASE_wide$ASE <= ASE_LOW),]
    
    
    RESULTS_E_Plus_ASE_wide<-rbind(RESULTS_E_Plus_ASE_wide_1,RESULTS_E_Plus_ASE_wide_2)
    
    RESULTS_E_Plus_ASE_wide$CLASS_E_Plus_ASE<-"E_Plus_ASE"
    
    
    cat("RESULTS_E_Plus_ASE_wide_4_Hello_world\n")
    cat(str(RESULTS_E_Plus_ASE_wide))
    cat("\n")
    cat("distribution_SIG_ASE_pool\n")
    cat(sprintf(as.character(names(summary(as.factor(RESULTS_E_Plus_ASE_wide$factor4))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(RESULTS_E_Plus_ASE_wide$factor4)))))
    cat("\n")
    cat("distribution_SIG_ASE_Label\n")
    cat(sprintf(as.character(names(summary(as.factor(RESULTS_E_Plus_ASE_wide$Label))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(RESULTS_E_Plus_ASE_wide$Label)))))
    cat("\n")
    cat("distribution_SIG_ASE_Cell_Type\n")
    cat(sprintf(as.character(names(summary(as.factor(RESULTS_E_Plus_ASE_wide$Cell_Type))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(RESULTS_E_Plus_ASE_wide$Cell_Type)))))
    cat("\n")
    
  }
  
  
  
  # distribution_SIG_ASE_Cell_Type
  # K562 CHRF HL60 THP1
  # 48 74 48 13
  
 
  
  #### merge CLASSES----
  
  indx.int<-c(which(colnames(RESULTS_E_Plus_ASE_wide) == "REAL_TILE_Plus_carried_variants"),
              which(colnames(RESULTS_E_Plus_ASE_wide) == "Cell_Type"),
              which(colnames(RESULTS_E_Plus_ASE_wide) == "CLASS_E_Plus_ASE"))
  
  
  RESULTS_E_Plus_ASE_wide_subset<-unique(RESULTS_E_Plus_ASE_wide[,indx.int])
  
  cat("RESULTS_E_Plus_ASE_wide_subset\n")
  cat(str(RESULTS_E_Plus_ASE_wide_subset))
  cat("\n")
  
  indx.int<-c(which(colnames(RESULTS_SIG_ASE) == "REAL_TILE_Plus_carried_variants"),
              which(colnames(RESULTS_SIG_ASE) == "Cell_Type"),
              which(colnames(RESULTS_SIG_ASE) == "CLASS_ASE"))
  
  
  RESULTS_SIG_ASE_subset<-unique(RESULTS_SIG_ASE[,indx.int])
  
  cat("RESULTS_SIG_ASE_subset\n")
  cat(str(RESULTS_SIG_ASE_subset))
  cat("\n")
  
  indx.int<-c(which(colnames(RESULTS_SIG_enhancer) == "REAL_TILE_Plus_carried_variants"),
              which(colnames(RESULTS_SIG_enhancer) == "Cell_Type"),
              which(colnames(RESULTS_SIG_enhancer) == "CLASS_enhancer"))
  
  
  RESULTS_SIG_enhancer_subset<-unique(RESULTS_SIG_enhancer[,indx.int])
  
  cat("RESULTS_SIG_enhancer_subset\n")
  cat(str(RESULTS_SIG_enhancer_subset))
  cat("\n")
  
  
  MERGE_CLASS<-merge(RESULTS_E_Plus_ASE_wide_subset,
                     RESULTS_SIG_ASE_subset,
                     by=c("REAL_TILE_Plus_carried_variants","Cell_Type"),
                     all=T)
  
  cat("MERGE_CLASS_1\n")
  cat(str(MERGE_CLASS))
  cat("\n")
  
  MERGE_CLASS<-merge(MERGE_CLASS,
                     RESULTS_SIG_enhancer_subset,
                     by=c("REAL_TILE_Plus_carried_variants","Cell_Type"),
                     all=T)
  
  cat("MERGE_CLASS_2\n")
  cat(str(MERGE_CLASS))
  cat("\n")
  
  
  MERGE_CLASS$DEF_CLASS<-"NA"
  
  MERGE_CLASS$DEF_CLASS[!is.na(MERGE_CLASS$CLASS_ASE)]<-MERGE_CLASS$CLASS_ASE[!is.na(MERGE_CLASS$CLASS_ASE)]
  
  MERGE_CLASS$DEF_CLASS[!is.na(MERGE_CLASS$CLASS_enhancer)]<-MERGE_CLASS$CLASS_enhancer[!is.na(MERGE_CLASS$CLASS_enhancer)]
  
  
  MERGE_CLASS$DEF_CLASS[!is.na(MERGE_CLASS$CLASS_E_Plus_ASE)]<-MERGE_CLASS$CLASS_E_Plus_ASE[!is.na(MERGE_CLASS$CLASS_E_Plus_ASE)]
  
  cat("MERGE_CLASS_3\n")
  cat(str(MERGE_CLASS))
  cat("\n")
  cat("distribution_SIG_ASE_Cell_Type\n")
  cat(sprintf(as.character(names(summary(as.factor(MERGE_CLASS$DEF_CLASS))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(MERGE_CLASS$DEF_CLASS)))))
  cat("\n")
  
  
  MERGE_CLASS$DEF_CLASS<-factor(MERGE_CLASS$DEF_CLASS,
                                levels=c("E_Plus_ASE","enhancer","ASE"),
                                ordered=T)
  
  indx.int<-c(which(colnames(MERGE_CLASS) == "REAL_TILE_Plus_carried_variants"),
              which(colnames(MERGE_CLASS) == "Cell_Type"),
              which(colnames(MERGE_CLASS) == "DEF_CLASS"))
  
  
  MERGE_CLASS_subset<-unique(MERGE_CLASS[,indx.int])
  
  cat("MERGE_CLASS_subset\n")
  cat(str(MERGE_CLASS_subset))
  cat("\n")
  
  ##### ALL TILES VERSION -----
  
  
 
  MPRA_Real_Tile_NO_FILTER<-merge(MPRA_Real_Tile_NO_FILTER,
                        MERGE_CLASS_subset,
                        by=c("REAL_TILE_Plus_carried_variants","Cell_Type"),
                        all=T)
  
  cat("MPRA_Real_Tile_NO_FILTER_1\n")
  cat(str(MPRA_Real_Tile_NO_FILTER))
  cat("\n")
  
  setwd(out2)
  
  saveRDS(MPRA_Real_Tile_NO_FILTER,file="MPRA_Real_Tile_QC2_NO_FILTERED.rds")
  
  
  
  #### SAVE ACTIVE TILES AND MPRA_Real_Tile_QC2_PASS----
  
  setwd(out2)
  
  write.table(MERGE_CLASS_subset,file="ACTIVE_TILES.txt", sep="\t", quote=F,row.names = F)
  
  write.table(MERGE_CLASS,file="ACTIVE_TILES_for_cummulative.txt", sep="\t", quote=F,row.names = F)
  
  
  MPRA_Real_Tile<-droplevels(MPRA_Real_Tile)
  
  MPRA_Real_Tile<-merge(MPRA_Real_Tile,
                        MERGE_CLASS_subset,
                     by=c("REAL_TILE_Plus_carried_variants","Cell_Type"),
                     all=T)
  
  cat("MPRA_Real_Tile_DEF\n")
  cat(str(MPRA_Real_Tile))
  cat("\n")
  
  
  saveRDS(MPRA_Real_Tile,file="MPRA_Real_Tile_QC2_PASS.rds")
  
}

########################################################################################  ACTIVATE EVERYTHING UPSTREAM


volcano_graphs= function(option_list)
{
  
  library("BiocManager", lib="/nfs/team151/software/manuel_R_libs_4_1//")
  library("memoise", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
  library("BiocGenerics", lib="/nfs/team151/software/manuel_R_libs_4_1/")
  library("Biobase", lib="/nfs/team151/software/manuel_R_libs_4_1/")
  library("GWASTools", lib="/nfs/team151/software/manuel_R_libs_4_1/")
  library("GENESIS", lib="/nfs/team151/software/manuel_R_libs_4_1/")
  
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
 
  #### Read ACTIVE TILES AND MPRA_Real_Tile_QC2_PASS----
  
  setwd(out2)
  
  MERGE_CLASS_subset<-read.table(file="ACTIVE_TILES.txt", sep="\t")
  
  cat("MERGE_CLASS_subset_DEF\n")
  cat(str(MERGE_CLASS_subset))
  cat("\n")
  
  MPRA_Real_Tile<-readRDS(file="MPRA_Real_Tile_QC2_PASS.rds")
 
  
  cat("MPRA_Real_Tile_DEF\n")
  cat(str(MPRA_Real_Tile))
  cat("\n")
  
  
 
  
  
  #### path2 ----
  
  path2<-paste(out2,'Volcano','/', sep='')
  
  cat("path2\n")
  cat(sprintf(as.character(path2)))
  cat("\n")
  
  if (file.exists(path2)){
    
    unlink(path2)
    dir.create(file.path(path2))
    
    
  } else {
    dir.create(file.path(path2))
    
  }
  
  
  
  
  ##### LOOP graphs ----
  
  list_graphs<-list()
  
  
  Cell_Types_vector<-levels(MPRA_Real_Tile$Cell_Type)
  
  cat("Cell_Types_vector_3\n")
  cat(str(Cell_Types_vector))
  cat("\n")
  
  represented_features<-c("dna_pool","rna_pool","dna.REF.vs.dna.ALT","rna.REF.vs.rna.ALT","LogFC","ASE")
  
  features<-unique(MPRA_Real_Tile$variable)
  
  features<-features[which(features%in%represented_features)]
  
  
  cat("features\n")
  cat(str(features))
  cat("\n")
  cat(sprintf(as.character(features)))
  cat("\n")
  
  indx.int<-c(which(features == "dna_pool"),which(features == "rna_pool"),which(features == "LogFC"))
  
  df_representation_enhancer_logpval<-as.data.frame(cbind(c(features[indx.int]),rep("enhancer_logpval",length(indx.int))), stringsAsFactors=F)
  
  colnames(df_representation_enhancer_logpval)<-c("xaxis","yaxis")
  
  cat("df_representation_enhancer_logpval\n")
  cat(str(df_representation_enhancer_logpval))
  cat("\n")
  
  indx.int<-c(which(features == "dna.REF.vs.dna.ALT"),which(features == "rna.REF.vs.rna.ALT"),which(features == "ASE"))
  
  
  df_representation_ASE_logpval<-as.data.frame(cbind(c(features[c(indx.int)]),rep("ASE_logpval",length(indx.int))), stringsAsFactors=F)
  
  colnames(df_representation_ASE_logpval)<-c("xaxis","yaxis")
  
  cat("df_representation_ASE_logpval\n")
  cat(str(df_representation_ASE_logpval))
  cat("\n")
  
  
  df_representation<-rbind(df_representation_enhancer_logpval,
                           df_representation_ASE_logpval)
  
  
  cat("df_representation\n")
  cat(str(df_representation))
  cat("\n")
  
  
  features<-unique(df_representation$xaxis)
  cat(sprintf(as.character(features)))
  cat("\n")
  
  for(i in 1:length(Cell_Types_vector))
  {
    Cell_Types_vector_sel<-Cell_Types_vector[i]
    
    cat("----------------------------------------------------------------------------------------------------------------------------------------------------->\t")
    cat(sprintf(as.character(Cell_Types_vector_sel)))
    cat("\n")
    
    
    MPRA_Real_Tile_sel<-MPRA_Real_Tile[which(MPRA_Real_Tile$Cell_Type == Cell_Types_vector_sel),]
    
    cat("MPRA_Real_Tile_sel\n")
    cat(str(MPRA_Real_Tile_sel))
    cat("\n")
    
    #### path3 ----
    
    path3<-paste(out2,'Volcano','/',Cell_Types_vector_sel,'/', sep='')
    
    cat("path3\n")
    cat(sprintf(as.character(path3)))
    cat("\n")
    
    if (file.exists(path3)){
      
      
      
      
    } else {
      dir.create(file.path(path3))
      
    }
    
    
    for(l in 1:length(features))
    {
      features_sel<-features[l]
      
      cat("------------------------------------------------------------------------------------------------->\t")
      cat(sprintf(as.character(features_sel)))
      cat("\n")
      
      df_representation_sel<-df_representation[which(df_representation$xaxis == features_sel),]
      
      cat("df_representation_sel\n")
      cat(str(df_representation_sel))
      cat("\n")
      
      
      yaxis_sel<-unique(df_representation_sel$yaxis)
      
      cat("yaxis_sel\n")
      cat(str(yaxis_sel))
      cat("\n")
      
      REP<-MPRA_Real_Tile_sel[which(MPRA_Real_Tile_sel$variable == features_sel),]
      
      cat("REP\n")
      cat(str(REP))
      cat("\n")
      
      indx.yaxis<-which(colnames(REP) == yaxis_sel)
      
      cat("indx.yaxis\n")
      cat(str(indx.yaxis))
      cat("\n")
     
      
      A<-summary(REP[,indx.yaxis])
      
      
      
      cat("A\n")
      cat(sprintf(as.character(names(A))))
      cat("\n")
      cat(sprintf(as.character(A)))
      cat("\n")
      
      max_value<-A[6]
      min_value<-A[1]
      
      
      step<-round((max_value-min_value)/5,0)
      
      cat("step_logpval\n")
      cat(sprintf(as.character(step)))
      cat("\n")
      # if(step)
      
      breaks.logpval<-unique(seq(min_value,max_value+step, by=step))
      labels.logpval<-as.character(round(breaks.logpval),1)
      
      
      cat("labels.logpval\n")
      cat(sprintf(as.character(labels.logpval)))
      cat("\n")
      
      
      if(features_sel == "LogFC")
      {
        cat("HELLO_WORLD_1\n")
        
        A<-summary(REP$value)
        
        
        
        cat("A\n")
        cat(sprintf(as.character(names(A))))
        cat("\n")
        cat(sprintf(as.character(A)))
        cat("\n")
        
        
        
        max_log_value<-A[6]
        min_log_value<--1.5
        
        #min_log_value
        breaks.log_value<-sort(c(0,seq(min_log_value,max_log_value+0.5, by=0.5)))
        
        
        #labels_value<-as.character(round(logratio2foldchange(breaks.log_value, base=2),2))
        labels_value<-as.character(round(breaks.log_value),2)
        
        
      }
      if(features_sel == "ASE")
      {
        cat("HELLO_WORLD_2\n")
        
        
        A<-summary(REP$value)
        
        
        
        cat("A\n")
        cat(sprintf(as.character(names(A))))
        cat("\n")
        cat(sprintf(as.character(A)))
        cat("\n")
        
        
        
        # max_value<-A[6]
        #min_value<-A[1]
        min_value<-0.4
        max_value<-3
        
        step<-(max_value-min_value)/5
        
        breaks.value<-sort(c(1,seq(min_value,max_value,by=step)))
        labels_value<-as.character(round(breaks.value,1))
        
      }
      
      OTHER_FEATURES<-c("dna_pool","rna_pool","dna.REF.vs.dna.ALT","rna.REF.vs.rna.ALT")
      
      if(features_sel%in%OTHER_FEATURES)
      {
        cat("HELLO_WORLD_3\n")
        
        
        A<-summary(REP$value)
        
        
        
        cat("A\n")
        cat(sprintf(as.character(names(A))))
        cat("\n")
        cat(sprintf(as.character(A)))
        cat("\n")
        
        
        
        max_value<-A[6]
        min_value<-A[1]
        
        step<-(max_value-min_value)/5
        
        breaks.log_value<-log10(seq(min_value,max_value+step,by=step)+0.000000001)
        #breaks.log_value<-c(seq(-5,25, by=0.5))
        
        labels_value<-as.character(round(10^(breaks.log_value),1))
        
      }
      
      cat("labels_value\n")
      cat(sprintf(as.character(labels_value)))
      cat("\n")
        
     
      
      #### Print Klaudia file ----
      #### path4
      
      path4<-paste(out2,'Volcano','/',Cell_Types_vector_sel,'/',yaxis_sel,'/', sep='')
      
      cat("path4\n")
      cat(sprintf(as.character(path4)))
      cat("\n")
      
      
      if (file.exists(path4)){
        
        
        
        
      } else {
        dir.create(file.path(path4))
        
      }
      
      setwd(path4)

      write.table(REP, file=paste("Table_for_volcano_",yaxis_sel,"_",features_sel,".tsv",sep=''),
                  sep="\t",
                  quote=F,
                  row.names = F)
       
        
        
        
      X_parameter_array<-c("DEF_CLASS","factor4","Label","Label_2")
      
      cat("X_parameter_array:\t")
      cat(str(X_parameter_array))
      cat("\n")
      
      
      for(iteration_X_parameter_array in 1:length(X_parameter_array))
      {
        
        
        ####
        
        X_parameter_array_sel<-X_parameter_array[iteration_X_parameter_array]
        
        cat("X_parameter_array_sel:\t")
        cat(sprintf(as.character(X_parameter_array_sel)))
        cat("\n")
        
        indx_x_parameter_sel<-which(colnames(REP) == X_parameter_array_sel)
        
        
        
        cat("indx_x_parameter_sel:\t")
        cat(sprintf(as.character(indx_x_parameter_sel)))
        cat("\n")
        
        ####
        
        selection_labels<-levels(REP[,indx_x_parameter_sel])
        
               
        if(X_parameter_array_sel == "DEF_CLASS")
        {
          mycols <- c("blue","dark cyan")
          
          if(yaxis_sel == "enhancer_logpval")
          {
            selection_labels<-c("enhancer","E_Plus_ASE")
          }
          if(yaxis_sel == "ASE_logpval")
          {
            selection_labels<-c("ASE","E_Plus_ASE")
          }
          
        }
        if(X_parameter_array_sel == "factor4")
        {
          mycols <- c("red","black","blue","dark cyan")
        }
        if(X_parameter_array_sel == "Label_2")
        {
          mycols <- c("red","black","blue","dark cyan")
        }
        if(X_parameter_array_sel == "Label")
        {
          mycols <- c("palevioletred2","black",'#1877C9','#32A852','#C9244B','#D45E85','#87447B','#553B68','#D6B8E6',
                             "#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C","#D7C1B1", "#AD6F3B", "#689030", "red")
        }
        
        
        #### Freq table 
        
        REP.dt<-data.table(REP, key=colnames(REP)[indx_x_parameter_sel])
        
        
        
        Freq.table<-as.data.frame(REP.dt[,.N,by=key(REP.dt)], stringsAsFactors=F)
        
        colnames(Freq.table)[which(colnames(Freq.table) == "N")]<-"Total"
        
        cat("Freq.table_\n")
        cat(str(Freq.table))
        cat("\n")
        
        indx_x_parameter_sel_Freq_table<-which(colnames(Freq.table) == X_parameter_array_sel)
        
        indx.yaxis<-which(colnames(REP) == yaxis_sel)
        
        cat("indx.yaxis\n")
        cat(str(indx.yaxis))
        cat("\n")
        
        #### volcanos
        
        ylabel="NA"
          
        if(yaxis_sel == "enhancer_logpval")
        {
          ylabel<-"-log10pval_Enhancer"
        }
        if(yaxis_sel == "ASE_logpval")
        {
          ylabel<-"-log10pval_ASE"
        }
        
        OTHER_FEATURES<-c("dna_pool","rna_pool","dna.REF.vs.dna.ALT","rna.REF.vs.rna.ALT")
        
        cat("selection_labels\n")
        cat(str(selection_labels))
        cat("\n")
        
        space_x_parameter<-REP[,indx_x_parameter_sel]
        
        cat("space_x_parameter\n")
        cat(str(space_x_parameter))
        cat("\n")
        
        
        
        indx.selection<-!is.na(space_x_parameter)
        
        cat("indx.selection\n")
        cat(str(indx.selection))
        cat("\n")
        
        
        
        if(features_sel%in%OTHER_FEATURES)
        {
         
          cat("HELLO_WORLD_3_3\n")
          
          volcano<-ggplot(data=REP,
                          aes(x=log10(value +0.000000001), 
                              y=REP[,indx.yaxis])) +
            geom_point(size=4, color="gray")+
            geom_point(data=REP[indx.selection,],
                       aes(x=log10(value +0.000000001),
                           y=REP[indx.selection,indx.yaxis],
                           color=REP[indx.selection,indx_x_parameter_sel]),
                       size=4)+
            theme_bw()+
            theme(axis.title.y=element_text(size=24, family="sans"),
                  axis.title.x=element_text(size=24, family="sans"),
                  axis.text.y=element_text(angle=0,size=24, color="black", family="sans"),
                  axis.text.x=element_text(angle=0, size=24, color="black", family="sans"),
                  legend.title=element_text(size=16,color="black", family="sans"),
                  legend.text=element_text(size=12,color="black", family="sans"))+
            scale_x_continuous(name=features_sel, breaks=breaks.log_value,labels=labels_value, limits=c(breaks.log_value[1],breaks.log_value[length(breaks.log_value)]))+
            scale_y_continuous(name=ylabel,
                               breaks=breaks.logpval,labels=labels.logpval,
                               limits=c(breaks.logpval[1],breaks.logpval[length(breaks.logpval)]))+
            scale_color_manual(values=c("palevioletred2","black",'#1877C9','#32A852','#C9244B','#D45E85','#87447B','#553B68','#D6B8E6',
                                        "#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C","#D7C1B1", "#AD6F3B", "#689030", "red"), drop=F,
                               name="Totals",breaks=Freq.table[,indx_x_parameter_sel_Freq_table],
                               labels=paste(Freq.table[,indx_x_parameter_sel_Freq_table],
                                            Freq.table$Total, sep =' n= '))+
            theme(legend.position="bottom",legend.title=element_blank(), legend.text = element_text(size=10))+
            ggeasy::easy_center_title()
          
          #### PDF Klaudia
          
          v_parameter<-"NA" 
          
          
          if(yaxis_sel == "enhancer_logpval")
          {
           v_parameter<-0 
          }
          if(yaxis_sel == "ASE_logpval")
          {
            v_parameter<-1 
            
          }
          
          # pdfname<-paste("volcano_",yaxis_sel,"_",features_sel,"_",X_parameter_array_sel,".pdf",sep='')
          # pdf(file=pdfname, width=5, height=4, pointsize=12)
          # 
          # par(mai=c(0.9,0.9,0.3,0.2))
          # lab <- as.character(unique(REP[,indx_x_parameter_sel]))
          # 
          # lab<-lab[which(lab%in%selection_labels)]
          # 
          # # cat("lab\n")
          # # cat(sprintf(as.character(lab)))
          # # cat("\n")
          # 
          # ind <- which(REP[,indx.yaxis] <= -log10(0.05)); length(ind)
          # 
          # plot(log10(REP$value +0.000000001), REP[,indx.yaxis], ty="n", xlab=features_sel, ylab=ylabel, axes=F, cex.lab=1.2, cex.lab=1.3, xlim=c(breaks.log_value[1],breaks.log_value[length(breaks.log_value)]))
          # abline(v=v_parameter, col="black",lty=2,lwd=2)
          # 
          # points(REP$value[ind], REP[,indx.yaxis][ind], col="darkgrey", pch=19,lwd=4)
          # 
          # sREP <- REP[-ind,]
          # 
          # for (i in 1:length(lab))  {
          #   
          #   cat(sprintf(as.character(lab[i])))
          #   cat("\n")
          #   
          #   ind <- which(sREP[,indx_x_parameter_sel]==lab[i])
          #   
          #   cat("------------------->ind\n")
          #   cat(str(ind))
          #   cat("\n")
          #   
          #   
          #   points(sREP$value[ind], sREP[,indx.yaxis][ind], pch=19,lwd=4, col=mycols[i])
          # }
          # legend("topright", legend=lab, fill=mycols, border=mycols, bty="n")
          # axis(1, at=seq(breaks.log_value[1],breaks.log_value[length(breaks.log_value)]))
          # axis(2, las=1)
          # 
          # dev.off()
          
        }# length(features_sel[which(features_sel%in%OTHER_FEATURES)])>1
        
        if(features_sel == "ASE")
        {
          cat("HELLO_WORLD_2_2\n")
          
          
          volcano<-ggplot(data=REP,
                          aes(x=value, 
                              y=REP[,indx.yaxis])) +
            geom_point(size=4, color="gray")+
            geom_point(data=REP[indx.selection,],
                       aes(x=value, 
                           y=REP[indx.selection,indx.yaxis],
                           color=REP[indx.selection,indx_x_parameter_sel]),
                       size=4)+
            theme_bw()+
            theme(axis.title.y=element_text(size=24, family="sans"),
                  axis.title.x=element_text(size=24, family="sans"),
                  axis.text.y=element_text(angle=0,size=24, color="black", family="sans"),
                  axis.text.x=element_text(angle=0, size=24, color="black", family="sans"),
                  legend.title=element_text(size=16,color="black", family="sans"),
                  legend.text=element_text(size=12,color="black", family="sans"))+
            scale_x_continuous(name=features_sel, breaks=breaks.value,labels=labels_value, limits=c(breaks.value[1],breaks.value[length(breaks.value)]))+
            scale_y_continuous(name=ylabel,
                               breaks=breaks.logpval,labels=labels.logpval, 
                               limits=c(breaks.logpval[1],breaks.logpval[length(breaks.logpval)]))+
            scale_color_manual(values=c("palevioletred2","black",'#1877C9','#32A852','#C9244B','#D45E85','#87447B','#553B68','#D6B8E6',
                                        "#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C","#D7C1B1", "#AD6F3B", "#689030", "red"), drop=F,
                               name="Totals",breaks=Freq.table[,indx_x_parameter_sel_Freq_table],
                               labels=paste(Freq.table[,indx_x_parameter_sel_Freq_table],
                                            Freq.table$Total, sep =' n= '))+
            theme(legend.position="bottom",legend.title=element_blank(), legend.text = element_text(size=10))+
            geom_vline(xintercept=1,linetype="solid")+
            ggeasy::easy_center_title()
          
          
        
          
          
          
          
          #### PDF Klaudia
          
          pdfname<-paste("volcano_",yaxis_sel,"_",features_sel,"_",X_parameter_array_sel,".pdf",sep='')
          pdf(file=pdfname, width=5, height=4, pointsize=12)
          
          par(mai=c(0.9,0.9,0.3,0.2))
          lab <- as.character(unique(REP[,indx_x_parameter_sel]))
          
          lab<-lab[which(lab%in%selection_labels)]
          
          # cat("lab\n")
          # cat(sprintf(as.character(lab)))
          # cat("\n")
        
          ASE_accepted<-c("ASE","E_Plus_ASE")
          
          ind <- -which(REP$DEF_CLASS%in%ASE_accepted); length(ind)
          
          plot(REP$value, REP[,indx.yaxis], ty="n", xlab=features_sel, ylab=ylabel, axes=F, cex.lab=1.2, cex.lab=1.3, xlim=c(breaks.value[1],breaks.value[length(breaks.value)]))
          abline(v=1, col="black",lty=2,lwd=2)
          
          points(REP$value[ind], REP[,indx.yaxis][ind], col="darkgrey", pch=19,lwd=4)
          
          sREP <- REP[-ind,]
          
          for (i in 1:length(lab))  {
            
            cat(sprintf(as.character(lab[i])))
            cat("\n")
            
            ind <- which(sREP[,indx_x_parameter_sel]==lab[i])
            
            cat("------------------->ind\n")
            cat(str(ind))
            cat("\n")
            
            
            points(sREP$value[ind], sREP[,indx.yaxis][ind], pch=19,lwd=4, col=mycols[i])
          }
          legend("topright", legend=lab, fill=mycols, border=mycols, bty="n")
          axis(1, at=seq(breaks.value[1],breaks.value[length(breaks.value)]))
          axis(2, las=1)
          
          dev.off()
          
          
          
        }#features_sel == "ASE"
        
        if(features_sel == "LogFC")
        {
          cat("HELLO_WORLD_1_1\n")
          
          
          volcano<-ggplot(data=REP,
                          aes(x=value,
                              y=REP[,indx.yaxis])) +
            geom_point(size=4, color="gray")+
            geom_point(data=REP[indx.selection,],
                       aes(x=value,
                           y=REP[indx.selection,indx.yaxis],
                           color=REP[indx.selection,indx_x_parameter_sel]),
                       size=4)+
            theme_bw()+
            theme(axis.title.y=element_text(size=24, family="sans"),
                  axis.title.x=element_text(size=24, family="sans"),
                  axis.text.y=element_text(angle=0,size=24, color="black", family="sans"),
                  axis.text.x=element_text(angle=0, size=24, color="black", family="sans"),
                  legend.title=element_text(size=16,color="black", family="sans"),
                  legend.text=element_text(size=12,color="black", family="sans"))+
            scale_x_continuous(name=features_sel, breaks=breaks.log_value,labels=labels_value, limits=c(breaks.log_value[1],breaks.log_value[length(breaks.log_value)]))+
            scale_y_continuous(name=ylabel,
                               breaks=breaks.logpval,labels=labels.logpval,
                               limits=c(breaks.logpval[1],breaks.logpval[length(breaks.logpval)]))+
            scale_color_manual(values=c("palevioletred2","black",'#1877C9','#32A852','#C9244B','#D45E85','#87447B','#553B68','#D6B8E6',
                                        "#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C","#D7C1B1", "#AD6F3B", "#689030", "red"), drop=F,
                               name="Totals",breaks=Freq.table[,indx_x_parameter_sel_Freq_table],
                               labels=paste(Freq.table[,indx_x_parameter_sel_Freq_table],
                                            Freq.table$Total, sep =' n= '))+
            theme(legend.position="bottom",legend.title=element_blank(), legend.text = element_text(size=10))+
            ggeasy::easy_center_title()
          
          #### PDF Klaudia
          
          v_parameter<-"NA" 
          
          
          if(yaxis_sel == "enhancer_logpval")
          {
            v_parameter<-0 
          }
          if(yaxis_sel == "ASE_logpval")
          {
            v_parameter<-1 
            
          }
          
          pdfname<-paste("volcano_",yaxis_sel,"_",features_sel,"_",X_parameter_array_sel,".pdf",sep='')
          pdf(file=pdfname, width=5, height=4, pointsize=12)
          
          par(mai=c(0.9,0.9,0.3,0.2))
          lab <- as.character(unique(REP[,indx_x_parameter_sel]))
          
          lab<-lab[which(lab%in%selection_labels)]
          
          # cat("lab\n")
          # cat(sprintf(as.character(lab)))
          # cat("\n")
          
          
          enhancer_accepted<-c("enhancer","E_Plus_ASE")
          
          ind <- -which(REP$DEF_CLASS%in%enhancer_accepted); length(ind)
          
          
          plot(REP$value, REP[,indx.yaxis], ty="n", xlab=features_sel, ylab=ylabel, axes=F, cex.lab=1.2, cex.lab=1.3, xlim=c(breaks.log_value[1],breaks.log_value[length(breaks.log_value)]))
          abline(v=v_parameter, col="black",lty=2,lwd=2)
          
          points(REP$value[ind], REP[,indx.yaxis][ind], col="darkgrey", pch=19,lwd=4)
          
          sREP <- REP[-ind,]
          
          for (i in 1:length(lab))  {
            
            cat(sprintf(as.character(lab[i])))
            cat("\n")
            
            ind <- which(sREP[,indx_x_parameter_sel]==lab[i])
            
            cat("------------------->ind\n")
            cat(str(ind))
            cat("\n")
            
            
            points(sREP$value[ind], sREP[,indx.yaxis][ind], pch=19,lwd=4, col=mycols[i])
          }
          legend("topleft", legend=lab, fill=mycols, border=mycols, bty="n")
          axis(1, at=seq(breaks.log_value[1],breaks.log_value[length(breaks.log_value)]))
          axis(2, las=1)
          
          dev.off()
          
         
        }
          
        
        
        
        svgname<-paste("volcano_",yaxis_sel,"_",features_sel,"_",X_parameter_array_sel,".svg",sep='')
        makesvg = TRUE
        
        if (makesvg == TRUE)
        {
          ggsave(svgname, plot= volcano,
                 device="svg",
                 height=10, width=12)
        }
        
        # if(X_parameter_array_sel == "factor4")
        # {
        #   setwd(out2)
        #   
        #   svgname<-paste("volcano_",yaxis_sel,"_",features_sel,"_",X_parameter_array_sel,".svg",sep='')
        #   makesvg = TRUE
        #   
        #   if (makesvg == TRUE)
        #   {
        #     ggsave(svgname, plot= volcano,
        #            device="svg",
        #            height=10, width=12)
        #   }
        #   
        #   quit(status = 1)
        #   
        # }#
        
        
        
      }# iteration_X_parameter_array
    }#l features
  }#i
  
  
  # if(features_sel == "LogFC")
  # {
  #   ######################################################################
  #   quit(status = 1)
  # }
  # 
  # 
  # # If Nicole wants a box instead of the two axes, you can just remove the flag 'axes=F' in the plotting command, and then you also don't need the separate functions 'axis'.
 
  
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
    make_option(c("--enhancer_result_K562"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ASE_result_K562"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--enhancer_result_CHRF"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ASE_result_CHRF"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--enhancer_result_HL60"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ASE_result_HL60"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--enhancer_result_THP1"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ASE_result_THP1"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--fdr_Threshold"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--FC_Threshold"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ASE_Threshold"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--enhancer_logval_Threshold"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ASE_log_pval_Threshold"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Rosetta_df"), type="character", default=NULL,
                metavar="FILE.txt",
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
  
 data_wrangling(opt)
 QC2_barplots_plots(opt)
  data_wrangling_after_QC2_plots(opt)
  volcano_graphs(opt)
  
}


###########################################################################

system.time( main() )
