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

opt = NULL

options(warn = 1)




MPRAnalyze_model = function(option_list)
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
  
  #### READ and transform Cell_Type ----
  
  Cell_Type_sel = opt$Cell_Type
  

  
  
  ##### LOOP TO READ PER CT FILES ----
  
  setwd(out)
  
 
    
    cat("-------------->\t")
    cat(sprintf(as.character(Cell_Type_sel)))
    cat("\n")
    
    filename10<-paste("Annot_DEF_QC_pass_",Cell_Type_sel,".tsv",sep='')
    Annot_DEF<-as.data.frame(fread(file=filename10,sep="\t",header=T), stringsAsFactors=F)
    
    cat("Annot_DEF_0\n")
    cat(str(Annot_DEF))
    cat("\n")
    
    levels_batch<-unique(as.character(Annot_DEF$batch))
    
    cat("levels_batch\n")
    cat(str(levels_batch))
    cat("\n")
    
    levels_bc<-unique(as.character(Annot_DEF$bc))
    
    cat("levels_bc\n")
    cat(str(levels_bc))
    cat("\n")
    
    
    levels_condition<-unique(as.character(Annot_DEF$condition))
    
    cat("levels_condition\n")
    cat(str(levels_condition))
    cat("\n")
    
    Annot_DEF$batch<-factor(Annot_DEF$batch,
                            levels=levels_batch)
    
    Annot_DEF$bc<-factor(Annot_DEF$bc,
                            levels=levels_bc)
    
    Annot_DEF$condition<-factor(Annot_DEF$condition,
                            levels=levels_condition)
    
    cat("Annot_DEF_1\n")
    cat(str(Annot_DEF))
    cat("\n")
    
    
    filename10<-paste("matrix_gDNA_PRE_QC_pass_",Cell_Type_sel,".tsv",sep='')
    
    gDNA<-as.data.frame(fread(file=filename10,sep="\t",header=T), stringsAsFactors=F)
    
    # cat("gDNA\n")
    # cat(str(gDNA))
    # cat("\n")
    
    matrix_gDNA<-as.matrix(gDNA[,-which(colnames(gDNA) == "REAL_TILE_Plus_carried_variants")])
    
    row.names(matrix_gDNA)<-gDNA$REAL_TILE_Plus_carried_variants
    
    cat("matrix_gDNA\n")
    cat(str(matrix_gDNA))
    cat("\n")
    
    
    filename10<-paste("matrix_cDNA_PRE_QC_pass_",Cell_Type_sel,".tsv",sep='')
    
    cDNA<-as.data.frame(fread(file=filename10,sep="\t",header=T), stringsAsFactors=F)
    
    # cat("cDNA\n")
    # cat(str(cDNA))
    # cat("\n")
    
    matrix_cDNA<-as.matrix(cDNA[,-which(colnames(cDNA) == "REAL_TILE_Plus_carried_variants")])
    
    row.names(matrix_cDNA)<-cDNA$REAL_TILE_Plus_carried_variants
    
    cat("matrix_cDNA\n")
    cat(str(matrix_cDNA))
    cat("\n")
    
    
    filename10<-paste("indx.ctrls_QC_pass_",Cell_Type_sel,".tsv",sep='')
    
    indx_ctrls<-as.data.frame(fread(file=filename10,sep="\t",header=T), stringsAsFactors=F)
    
    indx_ctrls<-as.integer(indx_ctrls$x)
    
    cat("indx_ctrls\n")
    cat(str(indx_ctrls))
    cat("\n")
    
    
    
    
    
    filename10<-paste("dnaDepth_df_QC_pass_",Cell_Type_sel,".tsv",sep='')
    
    dnaDepth<-as.data.frame(fread(file=filename10,sep="\t",header=T), stringsAsFactors=F)
    
    dnaDepth<-as.numeric(dnaDepth$dnaDepth)
    
    cat("dnaDepth\n")
    cat(sprintf(as.character(dnaDepth)))
    cat("\n")
    
    
    
    
    filename10<-paste("rnaDepth_df_QC_pass_",Cell_Type_sel,".tsv",sep='')
    
    rnaDepth<-as.data.frame(fread(file=filename10,sep="\t",header=T), stringsAsFactors=F)
    
    rnaDepth<-as.numeric(rnaDepth$rnaDepth)
    
    cat("rnaDepth\n")
    cat(sprintf(as.character(rnaDepth)))
    cat("\n")
    
    
    #### MPRA normalization by UQ ----
    
    obj <- MpraObject(dnaCounts = matrix_gDNA,
                      rnaCounts = matrix_cDNA,
                      colAnnot = Annot_DEF,
                      controls = indx_ctrls)
    
    cat("obj_1\n")
    cat(str(obj))
    cat("\n")
    
    # obj <- estimateDepthFactors(obj,
    #                             lib.factor = "batch",
    #                             which.lib = "both",
    #                             depth.estimator = "uq")
    # cat("obj_2A\n")
    # cat(str(obj))
    # cat("\n")
    
    obj <- setDepthFactors(obj, dnaDepth = dnaDepth,
                           rnaDepth = rnaDepth)

    cat("obj_2B\n")
    cat(str(obj))
    cat("\n")
    
    ######## ASE ANALYSIS ----
    
    obj <- analyzeComparative(obj, dnaDesign = ~ batch + bc + condition,
                              rnaDesign = ~ condition,
                              reducedDesign = ~ 1)
    
    cat("obj_3\n")
    cat(str(obj))
    cat("\n")
    
    results_ASE <- testLrt(obj)
    
    cat("results_ASE\n")
    cat(str(results_ASE))
    cat("\n")
    
    
   
    
    ######## enhancer ANALYSIS ----
    
    
    obj <- analyzeQuantification(obj, dnaDesign = ~ batch + bc,
                                 rnaDesign = ~1)
    
    cat("obj_3\n")
    cat(str(obj))
    cat("\n")
    
    #### RESULTS ----
    
    
    alpha <- getAlpha(obj)
    
    cat("alpha\n")
    cat(str(alpha))
    cat("\n")
    
    results_enhancer <- testEmpirical(obj,
                             useControls = TRUE)
    
    cat("results_enhancer\n")
    cat(str(results_enhancer))
    cat("\n")
    
    #### SAVE ----
    
    setwd(out)
    
    filename10<-paste("MPRAnalyze_enhancer_results_",Cell_Type_sel,".txt",sep='')
    write.table(results_enhancer,file=filename10,sep="\t",row.names = T,quote = F)
    
    filename11<-paste("MPRAnalyze_enhancer_alphas_",Cell_Type_sel,".txt",sep='')
    write.table(alpha,file=filename11,sep="\t",row.names = T,quote = F)
   
    
    filename10<-paste("MPRAnalyze_ASE_results_ASE_",Cell_Type_sel,".txt",sep='')
    write.table(results_ASE,file=filename10,sep="\t",row.names = T,quote = F)
    
    
    
    
    # quit(status=1)
    
    
    
    
  
  
 
  
  
  
  
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
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Cell_Type"), type="character", default=NULL, 
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
  
  
  MPRAnalyze_model(opt)
  
}


###########################################################################

system.time( main() )
