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


library("MPRAnalyze", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("ttutils", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/")








# library("ps", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
# library("usethis", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
# library("withr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
# library("desc", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
# library("devtools", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")

opt = NULL

options(warn = 1)

QC_PASS_CT = function(option_list)
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
  
  ### Read initial selection----
  
  
  Initial_Selection<-as.data.frame(readRDS(file=opt$Initial_Selection) , stringsAsFactors=F)
  
  cat("Initial_Selection\n")
  cat(str(Initial_Selection))
  cat("\n")
  
  Initial_Selection<-unique(Initial_Selection[,-which(colnames(Initial_Selection) == "maf_origin")])
  
  cat("Initial_Selection_2\n")
  cat(str(Initial_Selection))
  cat("\n")
  
  Screened_variants<-Initial_Selection[!is.na(Initial_Selection$Fig2_Annot_Category),]
  
  cat("Screened_variants\n")
  cat(str(Screened_variants))
  cat("\n")
  cat(str(unique(Screened_variants$VAR)))
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
  cat(str(unique(Rosetta$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Rosetta$Label))))))
  cat("\n")
  
  
  
  
  
  # ###################################################
  # quit(status = 1)
  
  rm(Rosetta_extended)
  
  
  filename=paste('matrix_gDNA_PRE','.rds',sep='')
  gDNA_PRE <-readRDS(file=filename)
  
  # cat("gDNA_PRE\n")
  # cat(str(gDNA_PRE))
  # cat("\n")
  
  matrix_gDNA_PRE<-as.matrix(gDNA_PRE[,-which(colnames(gDNA_PRE) == "REAL_TILE_Plus_carried_variants")])

  row.names(matrix_gDNA_PRE)<-gDNA_PRE$REAL_TILE_Plus_carried_variants

  # matrix_gDNA_PRE<-matrix_gDNA_PRE+0.000000001
  
  
  cat("matrix_gDNA_PRE\n")
  cat(str(matrix_gDNA_PRE))
  cat("\n")

  rm(gDNA_PRE)
  
  filename=paste('matrix_cDNA_PRE','.rds',sep='')
  #cDNA_PRE <-as.data.frame(fread(file=filename, sep="\t",header=T), stringsAsFactors = F)
  cDNA_PRE <-readRDS(file=filename)
  
  
  # cat("cDNA_PRE\n")
  # cat(str(cDNA_PRE))
  # cat("\n")
  
  matrix_cDNA_PRE<-as.matrix(cDNA_PRE[,-which(colnames(cDNA_PRE) == "REAL_TILE_Plus_carried_variants")])

  row.names(matrix_cDNA_PRE)<-cDNA_PRE$REAL_TILE_Plus_carried_variants
  
  # matrix_cDNA_PRE<-matrix_cDNA_PRE+0.000000001

  cat("matrix_cDNA_PRE\n")
  cat(str(matrix_cDNA_PRE))
  cat("\n")

  rm(cDNA_PRE)
  
  
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
  
  # ##################################################################################################
  # quit(status = 1)
  
  #### Important all the Annot columns should be factors
  
  REF_bc<-paste("bc",seq(1,15,by=1),sep='')
  ALT_bc<-paste("bc",seq(16,30,by=1),sep='')
  
  Annot_cDNA_df$batch<-factor(Annot_cDNA_df$batch,
                              levels=levels(Annot_cDNA_df$batch))
  
  Annot_cDNA_df$bc<-factor(Annot_cDNA_df$bc,
                           levels=c(REF_bc,ALT_bc))
  
  Annot_cDNA_df$condition<-factor(Annot_cDNA_df$condition,
                                  levels=c("REF","ALT"))
  
  
  filename_1<-paste("dnaDepth_df",".tsv", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  dnaDepth_df<-as.data.frame(fread(file=filename_1,sep="\t",header=T), stringsAsFactors=F)
  
  cat("dnaDepth_df\n")
  cat(str(dnaDepth_df))
  cat("\n")
  
  filename_1<-paste("rnaDepth_df",".tsv", sep='')
  
  # cat("filename_1\n")
  # cat(sprintf(as.character(filename_1)))
  # cat("\n")
  
  rnaDepth_df<-as.data.frame(fread(file=filename_1,sep="\t",header=T), stringsAsFactors=F)
  
  cat("rnaDepth_df\n")
  cat(str(rnaDepth_df))
  cat("\n")
  
  
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
  
  
  
 
  
  #### QC_PASS matrixes ----
  
  indx.gDNA<-which(colnames(matrix_gDNA_PRE)%in%Annot_gDNA_df_QC_PASS$accepted_colnames)
  
  cat("indx.gDNA\n")
  cat(str(indx.gDNA))
  cat("\n")
  
  matrix_gDNA_PRE_QC_PASS<-as.matrix(matrix_gDNA_PRE[,indx.gDNA])
  row.names(matrix_gDNA_PRE_QC_PASS)<-row.names(matrix_gDNA_PRE)#$REAL_TILE_Plus_carried_variants
  
  
  cat("matrix_gDNA_PRE_QC_PASS\n")
  cat(str(matrix_gDNA_PRE_QC_PASS))
  cat("\n")
  
  
  cat(sprintf(as.character(head(row.names(matrix_gDNA_PRE_QC_PASS)))))
  cat("\n")
  
  
 
  indx.cDNA<-which(colnames(matrix_cDNA_PRE)%in%Annot_cDNA_df_QC_PASS$accepted_colnames)
  
  cat("indx.cDNA\n")
  cat(str(indx.cDNA))
  cat("\n")
  
  matrix_cDNA_PRE_QC_PASS<-as.matrix(matrix_cDNA_PRE[,indx.cDNA])
  row.names(matrix_cDNA_PRE_QC_PASS)<-row.names(matrix_cDNA_PRE)#$REAL_TILE_Plus_carried_variants
  
  cat("matrix_cDNA_PRE_QC_PASS\n")
  cat(str(matrix_cDNA_PRE_QC_PASS))
  cat("\n")
  
  #### QC_PASS dna/rna_Depth -----
  
  dnaDepth_df_QC_PASS<-dnaDepth_df[which(dnaDepth_df$matrix_colnames%in%Annot_gDNA_df_QC_PASS$accepted_colnames),]
  
  cat("dnaDepth_df_QC_PASS\n")
  cat(str(dnaDepth_df_QC_PASS))
  cat("\n")
  
  rnaDepth_df_QC_PASS<-rnaDepth_df[which(rnaDepth_df$matrix_colnames%in%Annot_cDNA_df_QC_PASS$accepted_colnames),]
  
  cat("rnaDepth_df_QC_PASS\n")
  cat(str(rnaDepth_df_QC_PASS))
  cat("\n")
  
  
  
  
  
  #### indexes ctrls ----
  
  
 
  
  Rosetta_NCGR<-Rosetta[which(Rosetta$Label == 'Negative_Control_Genomic_Regions'),]
  
  cat("Rosetta_NCGR\n")
  cat(str(Rosetta_NCGR))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Rosetta_NCGR$KEY))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Rosetta_NCGR$KEY)))))
  cat("\n")
  
  indx.ctrls<-which(row.names(matrix_gDNA_PRE_QC_PASS)%in%Rosetta_NCGR$REAL_TILE_Plus_carried_variants)
  
  
  
  check.matrix_gDNA_PRE_QC_PASS<-matrix_gDNA_PRE_QC_PASS[indx.ctrls,]
  
  cat("check.matrix_gDNA_PRE_0\n")
  cat(str(check.matrix_gDNA_PRE_QC_PASS))
  cat("\n")
  
  
  check.matrix_cDNA_PRE_QC_PASS<-matrix_cDNA_PRE_QC_PASS[indx.ctrls,]
  
  cat("check.matrix_cDNA_PRE_0\n")
  cat(str(check.matrix_cDNA_PRE_QC_PASS))
  cat("\n")
  
  
  FLAG_check_ctrls<-sum(row.names(check.matrix_gDNA_PRE_QC_PASS) == row.names(check.matrix_cDNA_PRE_QC_PASS))
  
  cat("FLAG_check_ctrls_0\n")
  cat(str(FLAG_check_ctrls))
  cat("\n")
  
  
  if(FLAG_check_ctrls == length(row.names(check.matrix_gDNA_PRE_QC_PASS)))
  {
    
    
  }else{
    
    cat(sprintf("Inconcordant ctrls\n"))
    
    quit(status=1)
  }
  
  
  #### LOOP cell type and analysis ----
  
  
  for(i in 1:length(levels(Replicates_merge$Cell_Type)))
  {
    Cell_Type_sel<-levels(Replicates_merge$Cell_Type)[i]
    
    cat("--->\t")
    cat(sprintf(as.character(Cell_Type_sel)))
    cat("\n")
    
    
    ### subset annot by cell type
    
    Annot_gDNA_df_QC_PASS_sel<-Annot_gDNA_df_QC_PASS[which(Annot_gDNA_df_QC_PASS$Cell_Type == Cell_Type_sel),]
    
    
    cat("Annot_gDNA_df_QC_PASS_sel\n")
    cat(str(Annot_gDNA_df_QC_PASS_sel))
    cat("\n")
    
    Annot_cDNA_df_QC_PASS_sel<-Annot_cDNA_df_QC_PASS[which(Annot_cDNA_df_QC_PASS$Cell_Type == Cell_Type_sel),]
    
    
    cat("Annot_cDNA_df_QC_PASS_sel\n")
    cat(str(Annot_cDNA_df_QC_PASS_sel))
    cat("\n")
    
    
    
    
    ### subset annot for obj
    
    
    
    Annot_gDNA_df_QC_PASS_sel_subset<-unique(Annot_gDNA_df_QC_PASS_sel[,c(which(colnames(Annot_gDNA_df_QC_PASS_sel) == "batch"),
                                                                  which(colnames(Annot_gDNA_df_QC_PASS_sel) == "bc"),
                                                                  which(colnames(Annot_gDNA_df_QC_PASS_sel) == "condition"))])
    
    cat("Annot_gDNA_df_QC_PASS_sel_subset\n")
    cat(str(Annot_gDNA_df_QC_PASS_sel_subset))
    cat("\n")
    
    
    
    Annot_cDNA_df_QC_PASS_sel_subset<-unique(Annot_cDNA_df_QC_PASS_sel[,c(which(colnames(Annot_cDNA_df_QC_PASS_sel) == "batch"),
                                                                          which(colnames(Annot_cDNA_df_QC_PASS_sel) == "bc"),
                                                                          which(colnames(Annot_cDNA_df_QC_PASS_sel) == "condition"))])
    
    cat("Annot_cDNA_df_QC_PASS_sel_subset\n")
    cat(str(Annot_cDNA_df_QC_PASS_sel_subset))
    cat("\n")
    
    
    
    
    ### Cell_type_sel matrixes and Depth
    
    indx.gDNA<-which(colnames(matrix_gDNA_PRE_QC_PASS)%in%Annot_gDNA_df_QC_PASS_sel$accepted_colnames)
    
    cat("indx.gDNA\n")
    cat(str(indx.gDNA))
    cat("\n")
    
    matrix_gDNA_PRE_QC_PASS_sel<-as.data.frame(matrix_gDNA_PRE_QC_PASS[,indx.gDNA])
    
    
    cat("matrix_gDNA_PRE_QC_PASS_sel_0\n")
    cat(str(matrix_gDNA_PRE_QC_PASS_sel))
    cat("\n")
    
    
    
    
    matrix_gDNA_PRE_QC_PASS_sel$REAL_TILE_Plus_carried_variants<-row.names(matrix_gDNA_PRE_QC_PASS_sel)
    row.names(matrix_gDNA_PRE_QC_PASS_sel)<-NULL
    
    
    cat("matrix_gDNA_PRE_QC_PASS_sel_1\n")
    cat(str(matrix_gDNA_PRE_QC_PASS_sel))
    cat("\n")
    
    
    
    
    indx.cDNA<-which(colnames(matrix_cDNA_PRE_QC_PASS)%in%Annot_cDNA_df_QC_PASS_sel$accepted_colnames)
    
    cat("indx.cDNA\n")
    cat(str(indx.cDNA))
    cat("\n")
    
    matrix_cDNA_PRE_QC_PASS_sel<-as.data.frame(matrix_cDNA_PRE_QC_PASS[,indx.cDNA])
    
    
    cat("matrix_cDNA_PRE_QC_PASS_sel_0\n")
    cat(str(matrix_cDNA_PRE_QC_PASS_sel))
    cat("\n")
    
    
    
    
    matrix_cDNA_PRE_QC_PASS_sel$REAL_TILE_Plus_carried_variants<-row.names(matrix_cDNA_PRE_QC_PASS_sel)
    row.names(matrix_cDNA_PRE_QC_PASS_sel)<-NULL
    
    
    cat("matrix_cDNA_PRE_QC_PASS_sel_1\n")
    cat(str(matrix_cDNA_PRE_QC_PASS_sel))
    cat("\n")
    
    
    
    
    #### QC_PASS dna/rna_Depth -----
    
    dnaDepth_df_QC_PASS_sel<-dnaDepth_df_QC_PASS[which(dnaDepth_df_QC_PASS$matrix_colnames%in%Annot_gDNA_df_QC_PASS_sel$accepted_colnames),]

    cat("dnaDepth_df_QC_PASS_sel\n")
    cat(str(dnaDepth_df_QC_PASS_sel))
    cat("\n")

    rnaDepth_df_QC_PASS_sel<-rnaDepth_df_QC_PASS[which(rnaDepth_df_QC_PASS$matrix_colnames%in%Annot_cDNA_df_QC_PASS_sel$accepted_colnames),]

    cat("rnaDepth_df_QC_PASS_sel\n")
    cat(str(rnaDepth_df_QC_PASS_sel))
    cat("\n")
    
    #### Annot_DEF----
    
    Annot_DEF<-droplevels(Annot_gDNA_df_QC_PASS_sel_subset)
    
    levels_batch<-unique(as.character(Annot_DEF$batch))
    
    Annot_DEF$batch<-factor(Annot_DEF$batch,
                            levels=levels_batch)
    
    cat("Annot_DEF\n")
    cat(str(Annot_DEF))
    cat("\n")
    cat(sprintf(as.character(unique(Annot_DEF$batch))))
    cat("\n")
  
    
    
    #### SAVE ----

    setwd(out)

    filename10<-paste("Annot_DEF_",type,"_",Cell_Type_sel,".tsv",sep='')
    write.table(Annot_DEF,file=filename10,sep="\t",row.names = F,quote = F)

    filename10<-paste("matrix_gDNA_PRE_",type,"_",Cell_Type_sel,".tsv",sep='')
    write.table(matrix_gDNA_PRE_QC_PASS_sel,file=filename10,sep="\t",row.names = F,quote = F)
    
    filename10<-paste("matrix_cDNA_PRE_",type,"_",Cell_Type_sel,".tsv",sep='')
    write.table(matrix_cDNA_PRE_QC_PASS_sel,file=filename10,sep="\t",row.names = F,quote = F)
    
    filename10<-paste("indx.ctrls_",type,"_",Cell_Type_sel,".tsv",sep='')
    write.table(indx.ctrls,file=filename10,sep="\t",row.names = F,quote = F)
    
    filename10<-paste("dnaDepth_df_",type,"_",Cell_Type_sel,".tsv",sep='')
    write.table(dnaDepth_df_QC_PASS_sel,file=filename10,sep="\t",row.names = F,quote = F)
    
    filename10<-paste("rnaDepth_df_",type,"_",Cell_Type_sel,".tsv",sep='')
    write.table(rnaDepth_df_QC_PASS_sel,file=filename10,sep="\t",row.names = F,quote = F)
    
    

  }#i Cell_Type_sel
}


# ###############################################################################################################################################################################################################################################
# quit(status = 1)



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
    make_option(c("--REAL_TILE_DF"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--LONG_MATRIX"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Initial_Selection"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--barcodes_per_tile"), type="character", default=NULL, 
                metavar="type1", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Threshold_EXP"), type="numeric", default=NULL, 
                metavar="Threshold_EXP", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--indir"), type="character", default=NULL, 
                metavar="indir", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--swapped_NEG.CTRLS"), type="character", default=NULL, 
                metavar="NEG.CTRLS", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
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
    make_option(c("--out"), type="character", default=NULL, 
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
  
  QC_PASS_CT(opt)
   
}


###########################################################################

system.time( main() )
