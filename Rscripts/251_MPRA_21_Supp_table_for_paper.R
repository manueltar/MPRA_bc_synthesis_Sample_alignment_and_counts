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
suppressMessages(library("reshape2", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("tzdb", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("cli", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("tidyverse", lib.loc="/nfs/team151/software/manuel_R_libs_4_1//"))


opt = NULL

options(warn = 1)



data_wrangling = function(option_list)
{
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
  
  #### MPRA Real Tile ----
  
  MPRA_Real_tile = readRDS(opt$MPRA_Real_Tile_QC2_NO_FILTERED)#, sep="\t", header = T), stringsAsFactors = F)
  
  cat("MPRA_Real_tile_\n")
  cat(str(MPRA_Real_tile))
  cat("\n")
  cat(str(unique(MPRA_Real_tile$VAR)))
  cat("\n")
  cat(str(unique(MPRA_Real_tile$carried_variants)))
  cat("\n")
  
  # cat(sprintf(as.character(names(summary(as.factor(MPRA_Real_tile$Label))))))
  # cat("\n")
  # cat(sprintf(as.character(summary(as.factor(MPRA_Real_tile$Label)))))
  # cat("\n")
  # cat(sprintf(as.character(names(summary(as.factor(MPRA_Real_tile$DEF_CLASS))))))
  # cat("\n")
  # cat(sprintf(as.character(summary(as.factor(MPRA_Real_tile$DEF_CLASS)))))
  # cat("\n")
  
  indx.int<-c(which(colnames(MPRA_Real_tile) == "Cell_Type"),
              which(colnames(MPRA_Real_tile) == "VAR"),which(colnames(MPRA_Real_tile) == "carried_variants"),which(colnames(MPRA_Real_tile) == "KEY"),which(colnames(MPRA_Real_tile) == "TILE"),which(colnames(MPRA_Real_tile) == "Tile"),
              which(colnames(MPRA_Real_tile) == "Label_2"),
              which(colnames(MPRA_Real_tile) == "REAL_TILE_Plus_carried_variants"),which(colnames(MPRA_Real_tile) == "QC2_CLASS"),
              which(colnames(MPRA_Real_tile) == "variable"),which(colnames(MPRA_Real_tile) == "value"),
              which(colnames(MPRA_Real_tile) == "enhancer_logpval"),which(colnames(MPRA_Real_tile) == "logpval_enhancer_pval.empirical"),
              which(colnames(MPRA_Real_tile) == "fdr"),which(colnames(MPRA_Real_tile) == "ASE_logpval"),
              which(colnames(MPRA_Real_tile) == "DEF_CLASS"))
  
  
  MPRA_Real_tile_subset<-unique(MPRA_Real_tile[,indx.int])
  
  
  cat("MPRA_Real_tile_subset_\n")
  cat(str(MPRA_Real_tile_subset))
  cat("\n")
  
  
  #### Reported parameters ----
 
  variable_accepted<-c("LogFC","ASE")
  
  
  MPRA_Real_tile_subset_variable_accepted<-unique(MPRA_Real_tile_subset[which(MPRA_Real_tile_subset$variable%in%variable_accepted),])
  
  cat("MPRA_Real_tile_subset_variable_accepted_\n")
  cat(str(MPRA_Real_tile_subset_variable_accepted))
  cat("\n")
  
  indx.keep<-c(which(colnames(MPRA_Real_tile_subset_variable_accepted) == "Cell_Type"),
               which(colnames(MPRA_Real_tile_subset_variable_accepted) == "VAR"),which(colnames(MPRA_Real_tile_subset_variable_accepted) == "carried_variants"),which(colnames(MPRA_Real_tile_subset_variable_accepted) == "KEY"),which(colnames(MPRA_Real_tile_subset_variable_accepted) == "TILE"),which(colnames(MPRA_Real_tile_subset_variable_accepted) == "Tile"),
               which(colnames(MPRA_Real_tile_subset_variable_accepted) == "Label_2"),
               which(colnames(MPRA_Real_tile_subset_variable_accepted) == "REAL_TILE_Plus_carried_variants"),which(colnames(MPRA_Real_tile_subset_variable_accepted) == "QC2_CLASS"),
               which(colnames(MPRA_Real_tile_subset_variable_accepted) == "enhancer_logpval"),which(colnames(MPRA_Real_tile_subset_variable_accepted) == "logpval_enhancer_pval.empirical"),
               which(colnames(MPRA_Real_tile_subset_variable_accepted) == "fdr"),which(colnames(MPRA_Real_tile_subset_variable_accepted) == "ASE_logpval"),
               which(colnames(MPRA_Real_tile_subset_variable_accepted) == "DEF_CLASS"))
  
  MPRA_Real_tile_subset_variable_accepted_wide<-as.data.frame(pivot_wider(MPRA_Real_tile_subset_variable_accepted,
                                                          id_cols=colnames(MPRA_Real_tile_subset_variable_accepted)[indx.keep],
                                                          names_from=variable,
                                                          values_from=value),stringsAsFactors=F)
  
 
  indx.reorder<-c(which(colnames(MPRA_Real_tile_subset_variable_accepted_wide) == "Cell_Type"),which(colnames(MPRA_Real_tile_subset_variable_accepted_wide) == "VAR"),which(colnames(MPRA_Real_tile_subset_variable_accepted_wide) == "carried_variants"),which(colnames(MPRA_Real_tile_subset_variable_accepted_wide) == "KEY"),which(colnames(MPRA_Real_tile_subset_variable_accepted_wide) == "TILE"),which(colnames(MPRA_Real_tile_subset_variable_accepted_wide) == "Tile"),
               which(colnames(MPRA_Real_tile_subset_variable_accepted_wide) == "Label_2"),
               which(colnames(MPRA_Real_tile_subset_variable_accepted_wide) == "REAL_TILE_Plus_carried_variants"),which(colnames(MPRA_Real_tile_subset_variable_accepted_wide) == "QC2_CLASS"),
               which(colnames(MPRA_Real_tile_subset_variable_accepted_wide) == "LogFC"),
               which(colnames(MPRA_Real_tile_subset_variable_accepted_wide) == "enhancer_logpval"),which(colnames(MPRA_Real_tile_subset_variable_accepted_wide) == "logpval_enhancer_pval.empirical"),
               which(colnames(MPRA_Real_tile_subset_variable_accepted_wide) == "ASE"),
               which(colnames(MPRA_Real_tile_subset_variable_accepted_wide) == "fdr"),which(colnames(MPRA_Real_tile_subset_variable_accepted_wide) == "ASE_logpval"),
               which(colnames(MPRA_Real_tile_subset_variable_accepted_wide) == "DEF_CLASS"))
  
  MPRA_Real_tile_subset_variable_accepted_wide<-MPRA_Real_tile_subset_variable_accepted_wide[,indx.reorder]
  
  cat("MPRA_Real_tile_subset_variable_accepted_wide_0\n")
  cat(str(MPRA_Real_tile_subset_variable_accepted_wide))
  cat("\n")
  cat(str(unique(MPRA_Real_tile_subset_variable_accepted_wide$VAR)))
  cat("\n")
  cat(str(unique(MPRA_Real_tile_subset_variable_accepted_wide$carried_variants)))
  cat("\n")
  cat(sprintf(as.character(names(summary(MPRA_Real_tile_subset_variable_accepted_wide$Label_2)))))
  cat("\n")
  cat(sprintf(as.character(summary(MPRA_Real_tile_subset_variable_accepted_wide$Label_2))))
  cat("\n")
  
  
  #### READ and transform LONG_MATRIX ----
  
  LONG_MATRIX = read.table(opt$LONG_MATRIX, sep="\t", stringsAsFactors = F, header = T)
  
  
  colnames(LONG_MATRIX)[which(colnames(LONG_MATRIX) == "Real_Tile")]<-"REAL_TILE"
  colnames(LONG_MATRIX)[which(colnames(LONG_MATRIX) == "Carried_variants")]<-"carried_variants"
  
  
  cat("----->LONG_MATRIX_0\n")
  cat(str(LONG_MATRIX))
  cat("\n")
  cat(str(unique(LONG_MATRIX$carried_variants)))
  cat("\n")
  
  
  indx.int<-c(which(colnames(LONG_MATRIX) == "carried_variants"),
              which(colnames(LONG_MATRIX) == "REAL_TILE"),
              which(colnames(LONG_MATRIX) == "chr"),which(colnames(LONG_MATRIX) == "start"),which(colnames(LONG_MATRIX) == "stop"),
              which(colnames(LONG_MATRIX) == "seq_name"),
              which(colnames(LONG_MATRIX) == "REST_no_bc"),which(colnames(LONG_MATRIX) == "Enzymes"))
  
  
  LONG_MATRIX_subset<-unique(LONG_MATRIX[,indx.int])
  
  
  cat("LONG_MATRIX_subset_0\n")
  cat(str(LONG_MATRIX_subset))
  cat("\n")
  
  LONG_MATRIX_subset$KEY<-gsub(";.+$","",LONG_MATRIX_subset$seq_name)
  
  LONG_MATRIX_subset$Tile<-gsub("^[^;]+;","",LONG_MATRIX_subset$seq_name)
  LONG_MATRIX_subset$Tile<-gsub(";.+$","",LONG_MATRIX_subset$Tile)
  
  LONG_MATRIX_subset$Allele<-gsub("^[^;]+;[^;]+;","",LONG_MATRIX_subset$seq_name)
  LONG_MATRIX_subset$Allele<-gsub(";.+$","",LONG_MATRIX_subset$Allele)
  
  LONG_MATRIX_subset$TILE<-"NA"
  LONG_MATRIX_subset$TILE[which(LONG_MATRIX_subset$Tile == "NINE")]<-"TILE_1"
  LONG_MATRIX_subset$TILE[which(LONG_MATRIX_subset$Tile == "TWO_THIRDS")]<-"TILE_2"
  LONG_MATRIX_subset$TILE[which(LONG_MATRIX_subset$Tile == "HALF")]<-"TILE_3"
  LONG_MATRIX_subset$TILE[which(LONG_MATRIX_subset$Tile == "ONE_THIRD")]<-"TILE_4"
  LONG_MATRIX_subset$TILE[which(LONG_MATRIX_subset$Tile == "A_TENTH")]<-"TILE_5"
  
  cat("LONG_MATRIX_subset_1\n")
  cat(str(LONG_MATRIX_subset))
  cat("\n")
  
  LONG_MATRIX_subset$REAL_TILE_Plus_carried_variants<-paste(paste(LONG_MATRIX_subset$KEY,LONG_MATRIX_subset$TILE, sep="__"),LONG_MATRIX_subset$carried_variants,sep=";")
  
  cat("----->LONG_MATRIX_subset_0\n")
  cat(str(LONG_MATRIX_subset))
  cat("\n")
  cat(str(unique(LONG_MATRIX_subset$carried_variants)))
  cat("\n")
  cat(str(unique(LONG_MATRIX_subset$REAL_TILE_Plus_carried_variants)))
  cat("\n")
  
  check<-LONG_MATRIX_subset$REAL_TILE_Plus_carried_variants[which(LONG_MATRIX_subset$REAL_TILE_Plus_carried_variants%in%MPRA_Real_tile_subset_variable_accepted_wide$REAL_TILE_Plus_carried_variants)]
  
  cat("check\n")
  cat(str(check))
  cat("\n")
  
  indx.int_2<-c(which(colnames(LONG_MATRIX_subset) == "REAL_TILE_Plus_carried_variants"),
                which(colnames(LONG_MATRIX_subset) == "KEY"),which(colnames(LONG_MATRIX_subset) == "TILE"), which(colnames(LONG_MATRIX_subset) == "Allele"),
                which(colnames(LONG_MATRIX_subset) == "chr"),which(colnames(LONG_MATRIX_subset) == "start"),which(colnames(LONG_MATRIX_subset) == "stop"),
              which(colnames(LONG_MATRIX_subset) == "REST_no_bc"),which(colnames(LONG_MATRIX_subset) == "Enzymes"))
  
  LONG_MATRIX_subset_2<-unique(LONG_MATRIX_subset[,indx.int_2])
  
  LONG_MATRIX_subset_2$nchar_initial<-nchar(LONG_MATRIX_subset_2$REST_no_bc)
  
  cat("LONG_MATRIX_subset_2\n")
  cat(str(LONG_MATRIX_subset_2))
  cat("\n")
  cat(sprintf(as.character(names(summary(LONG_MATRIX_subset_2$nchar_initial)))))
  cat("\n")
  cat(sprintf(as.character(summary(LONG_MATRIX_subset_2$nchar_initial))))
  cat("\n")
  
  
  LONG_MATRIX_subset_2$check_enzymes<-substr(LONG_MATRIX_subset_2$REST_no_bc,1,12)
  LONG_MATRIX_subset_2$remaining<-substr(LONG_MATRIX_subset_2$REST_no_bc,13,nchar(LONG_MATRIX_subset_2$REST_no_bc))
  
  LONG_MATRIX_subset_2$nchar_two<-nchar(LONG_MATRIX_subset_2$remaining)
  
  
  cat("LONG_MATRIX_subset_3\n")
  cat(str(LONG_MATRIX_subset_2))
  cat("\n")
  cat(sprintf(as.character(names(summary(LONG_MATRIX_subset_2$nchar_two)))))
  cat("\n")
  cat(sprintf(as.character(summary(LONG_MATRIX_subset_2$nchar_two))))
  cat("\n")
  
  
  check_2<-sum(LONG_MATRIX_subset_2$Enzymes != LONG_MATRIX_subset_2$check_enzymes)
  
  cat("check_2\n")
  cat(str(check_2))
  cat("\n")
  
  ### REV arme has 14 nucleotides ----
  
  
  LONG_MATRIX_subset_2$check_rev_arm<-substr(LONG_MATRIX_subset_2$remaining,nchar(LONG_MATRIX_subset_2$remaining)-14+1,nchar(LONG_MATRIX_subset_2$remaining))
  LONG_MATRIX_subset_2$remaining_2<-substr(LONG_MATRIX_subset_2$remaining,1,nchar(LONG_MATRIX_subset_2$remaining)-14)
  
  LONG_MATRIX_subset_2$nchar_three<-nchar(LONG_MATRIX_subset_2$remaining_2)
  
  cat("LONG_MATRIX_subset_4\n")
  cat(str(LONG_MATRIX_subset_2))
  cat("\n")
  cat(sprintf(as.character(names(summary(LONG_MATRIX_subset_2$nchar_three)))))
  cat("\n")
  cat(sprintf(as.character(summary(LONG_MATRIX_subset_2$nchar_three))))
  cat("\n")
  
  
  check_3<-as.factor(LONG_MATRIX_subset_2$check_rev_arm)
  
  cat("check_3\n")
  cat(sprintf(as.character(names(summary(check_3)))))
  cat("\n")
  cat(sprintf(as.character(summary(check_3))))
  cat("\n")
  
  
  indx.int_3<-c(which(colnames(LONG_MATRIX_subset_2) == "REAL_TILE_Plus_carried_variants"),
                which(colnames(LONG_MATRIX_subset_2) == "KEY"),which(colnames(LONG_MATRIX_subset_2) == "TILE"), which(colnames(LONG_MATRIX_subset_2) == "Allele"),
                which(colnames(LONG_MATRIX_subset_2) == "chr"),which(colnames(LONG_MATRIX_subset_2) == "start"),which(colnames(LONG_MATRIX_subset_2) == "stop"),
                which(colnames(LONG_MATRIX_subset_2) == "remaining_2"))
  
  LONG_MATRIX_subset_3<-unique(LONG_MATRIX_subset_2[,indx.int_3])
  
  colnames(LONG_MATRIX_subset_3)[which(colnames(LONG_MATRIX_subset_3) == "remaining_2")]<-"genomic_sequence"
  
  cat("LONG_MATRIX_subset_3\n")
  cat(str(LONG_MATRIX_subset_3))
  cat("\n")
  
  ###### MERGE with MPRA results ----------------
  
  cat("MPRA_Real_tile_subset_variable_accepted_wide_REMEMBER\n")
  cat(str(MPRA_Real_tile_subset_variable_accepted_wide))
  cat("\n")
  cat(str(unique(MPRA_Real_tile_subset_variable_accepted_wide$VAR)))
  cat("\n")
  cat(str(unique(MPRA_Real_tile_subset_variable_accepted_wide$carried_variants)))
  cat("\n")
  cat(sprintf(as.character(names(summary(MPRA_Real_tile_subset_variable_accepted_wide$Label_2)))))
  cat("\n")
  cat(sprintf(as.character(summary(MPRA_Real_tile_subset_variable_accepted_wide$Label_2))))
  cat("\n")
  
  Table_Supp<-merge(LONG_MATRIX_subset_3,
                    MPRA_Real_tile_subset_variable_accepted_wide,
                    by=c("REAL_TILE_Plus_carried_variants","KEY","TILE"),
                    all.y=T)
  
  cat("Table_Supp_0\n")
  cat(str(Table_Supp))
  cat("\n")
  cat(str(unique(Table_Supp$VAR)))
  cat("\n")
  cat(str(unique(Table_Supp$carried_variants)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_Supp$Label_2)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_Supp$Label_2))))
  cat("\n")
  
  indx.dep<-c(which(colnames(Table_Supp) == "Tile"))
  
  
  
  Table_Supp_subset<-unique(Table_Supp[,-indx.dep])
  
  cat("Table_Supp_subset_1\n")
  cat(str(Table_Supp_subset))
  cat("\n")
  
  indx.reorder<-c(which(colnames(Table_Supp_subset) == "VAR"),
                  which(colnames(Table_Supp_subset) == "carried_variants"),
                  which(colnames(Table_Supp_subset) == "REAL_TILE_Plus_carried_variants"),
                  which(colnames(Table_Supp_subset) == "KEY"),
                  which(colnames(Table_Supp_subset) == "TILE"),
                  which(colnames(Table_Supp_subset) == "chr"),
                  which(colnames(Table_Supp_subset) == "start"),
                  which(colnames(Table_Supp_subset) == "stop"),
                  which(colnames(Table_Supp_subset) == "Allele"),
                  which(colnames(Table_Supp_subset) == "genomic_sequence"),
                  which(colnames(Table_Supp_subset) == "Cell_Type"),
                  which(colnames(Table_Supp_subset) == "Label_2"),
                  which(colnames(Table_Supp_subset) == "QC2_CLASS"),
                  which(colnames(Table_Supp_subset) == "LogFC"),
                  which(colnames(Table_Supp_subset) == "enhancer_logpval"),which(colnames(Table_Supp_subset) == "logpval_enhancer_pval.empirical"),
                  which(colnames(Table_Supp_subset) == "ASE"),
                  which(colnames(Table_Supp_subset) == "fdr"),which(colnames(Table_Supp_subset) == "ASE_logpval"),
                  which(colnames(Table_Supp_subset) == "DEF_CLASS"))
  
  Table_Supp_subset<-droplevels(Table_Supp_subset[,indx.reorder])
  
  cat("Table_Supp_subset_2\n")
  cat(str(Table_Supp_subset))
  cat("\n")
  
  
  colnames(Table_Supp_subset)[which(colnames(Table_Supp_subset) == "Label_2")]<-"ASSAY_CLASS"
  colnames(Table_Supp_subset)[which(colnames(Table_Supp_subset) == "QC2_CLASS")]<-"QC_FILTER"
  colnames(Table_Supp_subset)[which(colnames(Table_Supp_subset) == "DEF_CLASS")]<-"Per_tile_experimental_class"
  
  Table_Supp_subset$Per_tile_experimental_class<-as.character(Table_Supp_subset$Per_tile_experimental_class)
  
  
  Table_Supp_subset$QC_FILTER[which(Table_Supp_subset$QC_FILTER == "NA")]<-"FAIL"
  Table_Supp_subset$Per_tile_experimental_class[is.na(Table_Supp_subset$Per_tile_experimental_class)]<-"NOT_ACTIVE"
  Table_Supp_subset$Per_tile_experimental_class[which(Table_Supp_subset$QC_FILTER == "FAIL")]<-"FAIL"
  
  
 
  check_FAIL<-Table_Supp_subset[which(Table_Supp_subset$QC_FILTER == "FAIL"),]
  
  cat("check_FAIL_3\n")
  cat(str(check_FAIL))
  cat("\n")
  cat(str(unique(check_FAIL$VAR)))
  cat("\n")
  cat(str(unique(check_FAIL$carried_variants)))
  cat("\n")
  cat(str(unique(check_FAIL$REAL_TILE_Plus_carried_variants)))
  cat("\n")
  cat(sprintf(as.character(names(summary(check_FAIL$ASSAY_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(check_FAIL$ASSAY_CLASS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(check_FAIL$QC_FILTER))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(check_FAIL$QC_FILTER)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(check_FAIL$Per_tile_experimental_class))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(check_FAIL$Per_tile_experimental_class)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(check_FAIL$Allele))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(check_FAIL$Allele)))))
  cat("\n")
  
  ############################### SAVE -------------------------------------------------------
  
  Table_Supp_subset$chr<-factor(Table_Supp_subset$chr,
                             levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                                      "chr17","chr18","chr19","chr20","chr21","chr22","chr23","chrX","chrY"),
                             ordered=T)
  
  Table_Supp_subset$Allele<-factor(Table_Supp_subset$Allele,
                                levels=c("REF","ALT"),
                                ordered=T)
  

  cat("Table_Supp_subset_3\n")
  cat(str(Table_Supp_subset))
  cat("\n")
  cat(str(unique(Table_Supp_subset$VAR)))
  cat("\n")
  cat(str(unique(Table_Supp_subset$carried_variants)))
  cat("\n")
  cat(str(unique(Table_Supp_subset$REAL_TILE_Plus_carried_variants)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_Supp_subset$ASSAY_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_Supp_subset$ASSAY_CLASS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Table_Supp_subset$QC_FILTER))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Table_Supp_subset$QC_FILTER)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Table_Supp_subset$Per_tile_experimental_class))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Table_Supp_subset$Per_tile_experimental_class)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Table_Supp_subset$Allele))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Table_Supp_subset$Allele)))))
  cat("\n")
  
  Table_Supp_subset<-Table_Supp_subset[order(Table_Supp_subset$chr,Table_Supp_subset$start,Table_Supp_subset$carried_variants,Table_Supp_subset$Allele),]
  
  cat("Table_Supp_subset_FINAL\n")
  cat(str(Table_Supp_subset))
  cat("\n")
  cat(str(unique(Table_Supp_subset$VAR)))
  cat("\n")
  cat(str(unique(Table_Supp_subset$carried_variants)))
  cat("\n")
  cat(str(unique(Table_Supp_subset$REAL_TILE_Plus_carried_variants)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_Supp_subset$ASSAY_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_Supp_subset$ASSAY_CLASS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Table_Supp_subset$QC_FILTER))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Table_Supp_subset$QC_FILTER)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Table_Supp_subset$Per_tile_experimental_class))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Table_Supp_subset$Per_tile_experimental_class)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Table_Supp_subset$Allele))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Table_Supp_subset$Allele)))))
  cat("\n")
  
  setwd(out2)
  
  saveRDS(file="Table_S2_MPRA_Results.rds",
              Table_Supp_subset)
 
}

Add_Rsid_and_cumulative_classification = function(option_list)
{
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
  
  #### MPRA Real Tile ----
  
  setwd(out2)
  
  Table_S2 = readRDS(file="Table_S2_MPRA_Results.rds") #, sep="\t", header = T), stringsAsFactors = F)
  
  cat("Table_S2_0\n")
  cat(str(Table_S2))
  cat("\n")
  cat(str(unique(Table_S2$VAR)))
  cat("\n")
  cat(str(unique(Table_S2$carried_variants)))
  cat("\n")
  
  #### Read ALL_dB file ----
  
  ALL_dB<-as.data.frame(fread(file=opt$ALL_dB,sep="\t") , stringsAsFactors=F)
  
  cat("ALL_dB_0\n")
  cat(str(ALL_dB))
  cat("\n")
  
  
  indx.dep<-c(which(colnames(ALL_dB) == "maf_origin"))
  
  ALL_dB_subset<-unique(ALL_dB[,-indx.dep])
  
  # cat("ALL_dB_subset\n")
  # cat(str(ALL_dB_subset))
  # cat("\n")
  
  rm(ALL_dB)
  
  ALL_dB_subset$Allelic_Series_ID<-paste(ALL_dB_subset$phenotype,ALL_dB_subset$block_no,sep='__')
  
  cat("ALL_dB_subset_1\n")
  cat(str(ALL_dB_subset))
  cat("\n")
  
  
  indx.int<-c(which(colnames(ALL_dB_subset) == "VAR"),which(colnames(ALL_dB_subset) == "rs"))
  
  ALL_FINAL<-unique(ALL_dB_subset[,indx.int])
  
  cat("ALL_FINAL\n")
  cat(str(ALL_FINAL))
  cat("\n")
  cat(str(unique(ALL_FINAL$VAR)))
  cat("\n")
  
  rm(ALL_dB_subset)
  
  ### Add Rsid ----
  
  Table_S2<-merge(ALL_FINAL,
                  Table_S2,
                  by="VAR",
                  all.y=T)
  
  colnames(Table_S2)[which(colnames(Table_S2) == "rs")]<-"Rsid"
  
  cat("Table_S2_1\n")
  cat(str(Table_S2))
  cat("\n")
  cat(str(unique(Table_S2$VAR)))
  cat("\n")
  cat(str(unique(Table_S2$carried_variants)))
  cat("\n")
  
  #### Common KEY with CUMMULATIVE_CLASSES ----
  
  Table_S2$KEY_Plus_carried_variants<-paste(Table_S2$KEY,Table_S2$carried_variants, sep=";")
  
  cat("Table_S2_2\n")
  cat(str(Table_S2))
  cat("\n")
  cat(str(unique(Table_S2$VAR)))
  cat("\n")
  cat(str(unique(Table_S2$carried_variants)))
  cat("\n")
  cat(str(unique(Table_S2$KEY_Plus_carried_variants)))
  cat("\n")
  
  ### CUMMULATIVE_CLASSES #----
  
  CUMMULATIVE_CLASSES<-readRDS(file=opt$CUMMULATIVE_CLASSES)
  
  
  cat("CUMMULATIVE_CLASSES\n")
  cat(str(CUMMULATIVE_CLASSES))
  cat("\n")
  cat(str(unique(CUMMULATIVE_CLASSES$KEY_Plus_carried_variants)))
  cat("\n")
 
  indx.int<-c(which(colnames(CUMMULATIVE_CLASSES) == "KEY_Plus_carried_variants"),
              which(colnames(CUMMULATIVE_CLASSES) == "Cell_Type"),which(colnames(CUMMULATIVE_CLASSES) == "enhancer_CLASS_TILES"),which(colnames(CUMMULATIVE_CLASSES) == "E_Plus_ASE_CLASS_TILES"))
  
  
  
  CUMMULATIVE_CLASSES_subset<-unique(CUMMULATIVE_CLASSES[which(CUMMULATIVE_CLASSES$Cell_Type != "ALL_CT"),indx.int])
 
  cat("CUMMULATIVE_CLASSES_subset\n")
  cat(str(CUMMULATIVE_CLASSES_subset))
  cat("\n")
  cat(str(unique(CUMMULATIVE_CLASSES_subset$KEY_Plus_carried_variants)))
  cat("\n")
  
  
  # CUMMULATIVE_CLASSES_subset_wide<-as.data.frame(pivot_wider(CUMMULATIVE_CLASSES_subset,
  #                                                                         id_cols=KEY_Plus_carried_variants,
  #                                                                         names_from=Cell_Type,
  #                                                                         values_from=c("enhancer_CLASS_TILES","E_Plus_ASE_CLASS_TILES")),stringsAsFactors=F)
  # 
  # cat("CUMMULATIVE_CLASSES_subset_wide\n")
  # cat(str(CUMMULATIVE_CLASSES_subset_wide))
  # cat("\n")
  # cat(str(unique(CUMMULATIVE_CLASSES_subset_wide$KEY_Plus_carried_variants)))
  # cat("\n")
  
  
  # Add cummulative classes with controls; keep KEY and add by KEY the TILES with enhancer and the TILES with E Plus ASE
  
  
  # Table_S2<-merge(Table_S2,
  #                 CUMMULATIVE_CLASSES_subset_wide,
  #                 by=c("KEY_Plus_carried_variants","Cell_Type"),
  #                 all=T)
  
  Table_S2<-merge(Table_S2,
                  CUMMULATIVE_CLASSES_subset,
                  by=c("KEY_Plus_carried_variants","Cell_Type"),
                  all=T)
  
  
  cat("Table_S2_3\n")
  cat(str(Table_S2))
  cat("\n")
  cat(str(unique(Table_S2$VAR)))
  cat("\n")
  cat(str(unique(Table_S2$carried_variants)))
  cat("\n")
  cat(str(unique(Table_S2$KEY_Plus_carried_variants)))
  cat("\n")
  
  
  
  
  
  ########## SAVE ----
  
  setwd(out2)

  Table_S2<-Table_S2[order(Table_S2$chr,Table_S2$start,Table_S2$carried_variants,Table_S2$Cell_Type,Table_S2$Allele),]

  cat("Table_S2_FINAL\n")
  cat(str(Table_S2))
  cat("\n")
  cat(str(unique(Table_S2$VAR)))
  cat("\n")
  cat(str(unique(Table_S2$carried_variants)))
  cat("\n")
  cat(str(unique(Table_S2$REAL_TILE_Plus_carried_variants)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S2$ASSAY_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S2$ASSAY_CLASS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Table_S2$QC_FILTER))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Table_S2$QC_FILTER)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Table_S2$Per_tile_experimental_class))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Table_S2$Per_tile_experimental_class)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Table_S2$Allele))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Table_S2$Allele)))))
  cat("\n")

  #
  write.table(file="Table_S2_MPRA_Results.tsv",
              Table_S2, sep="\t", quote=F, row.names = F)
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
    make_option(c("--LONG_MATRIX"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--MPRA_Real_Tile_QC2_NO_FILTERED"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--CUMMULATIVE_CLASSES"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ALL_dB"), type="character", default=NULL,
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
                        --REAL_TILE_DF FILE.txt
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
  
  
  # data_wrangling(opt)
  Add_Rsid_and_cumulative_classification(opt)
 
  
  
}


###########################################################################

system.time( main() )
