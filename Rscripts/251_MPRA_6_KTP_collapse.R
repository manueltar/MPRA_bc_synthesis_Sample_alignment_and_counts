
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



opt = NULL


opt = NULL

CLASSIFICATION_AT = function(option_list)
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
  
  
  #### ACTIVE_TILES ----
  
  ACTIVE_TILES<-as.data.frame(fread(file=opt$ACTIVE_TILES, sep="\t", header = T), stringsAsFactors = F)
 
  
  cat("ACTIVE_TILES\n")
  cat(str(ACTIVE_TILES))
  cat("\n")
  
  
  ACTIVE_TILES$KEY<-gsub(";.+$","",ACTIVE_TILES$REAL_TILE_Plus_carried_variants)
  ACTIVE_TILES$KEY<-gsub("__.+$","",ACTIVE_TILES$KEY)
  
 
  cat("ACTIVE_TILES_sandwich\n")
  cat(str(ACTIVE_TILES))
  cat("\n")
  
  
  ACTIVE_TILES$carried_variants<-gsub("[^;]+;","",ACTIVE_TILES$REAL_TILE_Plus_carried_variants)
  
  cat("ACTIVE_TILES_2\n")
  cat(str(ACTIVE_TILES))
  cat("\n")
  
  ACTIVE_TILES$KEY_Plus_carried_variants<-paste(ACTIVE_TILES$KEY,ACTIVE_TILES$carried_variants,sep=";")
  
  cat("ACTIVE_TILES_3\n")
  cat(str(ACTIVE_TILES))
  cat("\n")
  
  
  ACTIVE_TILES.dt<-data.table(ACTIVE_TILES, key=c("KEY_Plus_carried_variants","Cell_Type","DEF_CLASS"))
  
  
  ACTIVE_TILES_Fq<-as.data.frame(ACTIVE_TILES.dt[,.N,by=key(ACTIVE_TILES.dt)], stringsAsFactors=F)
  
  colnames(ACTIVE_TILES_Fq)[which(colnames(ACTIVE_TILES_Fq) == "N")]<-"instances"
  
  cat("ACTIVE_TILES_Fq_\n")
  cat(str(ACTIVE_TILES_Fq))
  cat("\n")
  
  
  
  ### enhancer classification
  
  enhancer_labels<-c("enhancer","E_Plus_ASE")
  
  
  ACTIVE_TILES_Fq_enhancer<-ACTIVE_TILES_Fq[which(ACTIVE_TILES_Fq$DEF_CLASS%in%enhancer_labels),]
  
  cat("ACTIVE_TILES_Fq_enhancer\n")
  cat(str(ACTIVE_TILES_Fq_enhancer))
  cat("\n")
  
  ACTIVE_TILES_Fq_enhancer.dt<-data.table(ACTIVE_TILES_Fq_enhancer, key=c("KEY_Plus_carried_variants","Cell_Type"))
  
  
  ACTIVE_TILES_Fq_enhancer_Fq<-as.data.frame(ACTIVE_TILES_Fq_enhancer.dt[,.(enhancer_CLASS_TILES=sum(instances)),
                                                                         by=key(ACTIVE_TILES_Fq_enhancer.dt)], stringsAsFactors=F)
  
  
  
  ACTIVE_TILES_Fq_enhancer_Fq$enhancer_CLASS<-"NA"
  
  ACTIVE_TILES_Fq_enhancer_Fq$enhancer_CLASS[ACTIVE_TILES_Fq_enhancer_Fq$enhancer_CLASS_TILES == 1]<-"1"
  ACTIVE_TILES_Fq_enhancer_Fq$enhancer_CLASS[ACTIVE_TILES_Fq_enhancer_Fq$enhancer_CLASS_TILES == 2]<-"2"
  ACTIVE_TILES_Fq_enhancer_Fq$enhancer_CLASS[ACTIVE_TILES_Fq_enhancer_Fq$enhancer_CLASS_TILES >= 3]<-"3 or >"
  
  cat("ACTIVE_TILES_Fq_enhancer_Fq_\n")
  cat(str(ACTIVE_TILES_Fq_enhancer_Fq))
  cat("\n")
  
  ### ASE classification
  
  ASE_labels<-c("ASE","E_Plus_ASE")
  
  
  ACTIVE_TILES_Fq_ASE<-ACTIVE_TILES_Fq[which(ACTIVE_TILES_Fq$DEF_CLASS%in%ASE_labels),]
  
  cat("ACTIVE_TILES_Fq_ASE\n")
  cat(str(ACTIVE_TILES_Fq_ASE))
  cat("\n")
  
  ACTIVE_TILES_Fq_ASE.dt<-data.table(ACTIVE_TILES_Fq_ASE, key=c("KEY_Plus_carried_variants","Cell_Type"))
  
  
  ACTIVE_TILES_Fq_ASE_Fq<-as.data.frame(ACTIVE_TILES_Fq_ASE.dt[,.(ASE_CLASS_TILES=sum(instances)),
                                                                         by=key(ACTIVE_TILES_Fq_ASE.dt)], stringsAsFactors=F)
  
  
  
  ACTIVE_TILES_Fq_ASE_Fq$ASE_CLASS<-"NA"
  
  ACTIVE_TILES_Fq_ASE_Fq$ASE_CLASS[ACTIVE_TILES_Fq_ASE_Fq$ASE_CLASS_TILES == 1]<-"1"
  ACTIVE_TILES_Fq_ASE_Fq$ASE_CLASS[ACTIVE_TILES_Fq_ASE_Fq$ASE_CLASS_TILES == 2]<-"2"
  ACTIVE_TILES_Fq_ASE_Fq$ASE_CLASS[ACTIVE_TILES_Fq_ASE_Fq$ASE_CLASS_TILES >= 3]<-"3 or >"
  
  cat("ACTIVE_TILES_Fq_ASE_Fq_\n")
  cat(str(ACTIVE_TILES_Fq_ASE_Fq))
  cat("\n")
  
  ### E_Plus_ASE classification
  
  E_Plus_ASE_labels<-c("E_Plus_ASE")
  
  
  ACTIVE_TILES_Fq_E_Plus_ASE<-ACTIVE_TILES_Fq[which(ACTIVE_TILES_Fq$DEF_CLASS%in%E_Plus_ASE_labels),]
  
  cat("ACTIVE_TILES_Fq_E_Plus_ASE\n")
  cat(str(ACTIVE_TILES_Fq_E_Plus_ASE))
  cat("\n")
  
  ACTIVE_TILES_Fq_E_Plus_ASE.dt<-data.table(ACTIVE_TILES_Fq_E_Plus_ASE, key=c("KEY_Plus_carried_variants","Cell_Type"))
  
  
  ACTIVE_TILES_Fq_E_Plus_ASE_Fq<-as.data.frame(ACTIVE_TILES_Fq_E_Plus_ASE.dt[,.(E_Plus_ASE_CLASS_TILES=sum(instances)),
                                                                             by=key(ACTIVE_TILES_Fq_E_Plus_ASE.dt)], stringsAsFactors=F)
  
  E_Plus_ASE_CLASS<-"NA"
  
  ACTIVE_TILES_Fq_E_Plus_ASE_Fq$E_Plus_ASE_CLASS[ACTIVE_TILES_Fq_E_Plus_ASE_Fq$E_Plus_ASE_CLASS_TILES == 1]<-"1"
  ACTIVE_TILES_Fq_E_Plus_ASE_Fq$E_Plus_ASE_CLASS[ACTIVE_TILES_Fq_E_Plus_ASE_Fq$E_Plus_ASE_CLASS_TILES == 2]<-"2"
  ACTIVE_TILES_Fq_E_Plus_ASE_Fq$E_Plus_ASE_CLASS[ACTIVE_TILES_Fq_E_Plus_ASE_Fq$E_Plus_ASE_CLASS_TILES >= 3]<-"3 or >"
  
  cat("ACTIVE_TILES_Fq_E_Plus_ASE_Fq_\n")
  cat(str(ACTIVE_TILES_Fq_E_Plus_ASE_Fq))
  cat("\n")
  
  
  
  setwd(out)
  
  write.table(ACTIVE_TILES_Fq_E_Plus_ASE_Fq,file="test.tsv",sep="\t",quote=F,row.names=F)
  
  
  #### MPRA Real Tile ----
  
  MPRA_Real_tile = readRDS(opt$MPRA_Real_tile_QC2_PASS)#, sep="\t", header = T), stringsAsFactors = F)
  
  cat("MPRA_Real_tile_\n")
  cat(str(MPRA_Real_tile))
  cat("\n")
  
  cat(sprintf(as.character(names(summary(as.factor(MPRA_Real_tile$Label))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(MPRA_Real_tile$Label)))))
  cat("\n")
  
  
  
  MPRA_Real_tile$VEP_DEF_LABELS<-revalue(MPRA_Real_tile$Label, c("intron_variant"="INTRON",
                                                         "intergenic_variant" = "INTERGENIC",
                                                         "TF_binding_site_variant"= "TFBS",
                                                         "regulatory_region_variant"= "REGULATORY",
                                                         "Kousik_variant" = "Kousik_variant",
                                                         "upstream_gene_variant" = "UPSTREAM",
                                                         "downstream_gene_variant" = "DOWNSTREAM",
                                                         "Negative_Control_Genomic_Regions" = "NCGR",
                                                         "POSITIVE_CTRL" = "ASE_CTRL",
                                                         "UNDETERMINED_CTRL"= "enhancer_CTRL"
  ))
  
  MPRA_Real_tile$VEP_DEF_LABELS<-factor(MPRA_Real_tile$VEP_DEF_LABELS,
                                levels=c("LOF","MISS","SYN","UTR5","UTR3",
                                         "INTRON","INTERGENIC","UPSTREAM","DOWNSTREAM","REGULATORY",
                                         "TFBS","SPLICE",
                                         "OTHER","NMD","NCT",
                                         "NCGR","ASE_CTRL","enhancer_CTRL","Kousik_variant"),
                                ordered=T)
  
  MPRA_Real_tile<-droplevels(MPRA_Real_tile)
  
  cat("MPRA_Real_tile\n")
  cat(str(MPRA_Real_tile))
  cat("\n")
  
  
  KEY_collapse<-unique(MPRA_Real_tile[,c(which(colnames(MPRA_Real_tile) == "KEY"),
                                         which(colnames(MPRA_Real_tile) == "carried_variants"),
                                         which(colnames(MPRA_Real_tile) == "VAR"),
                                         which(colnames(MPRA_Real_tile) == "chr"),
                                         which(colnames(MPRA_Real_tile) == "Cell_Type"),
                                         which(colnames(MPRA_Real_tile) == "factor4"),
                                         which(colnames(MPRA_Real_tile) == "Label"),
                                         which(colnames(MPRA_Real_tile) == "Label_2"),
                                         which(colnames(MPRA_Real_tile) == "VEP_DEF_LABELS"))])
  
  KEY_collapse$KEY_Plus_carried_variants<-paste(KEY_collapse$KEY,KEY_collapse$carried_variants,sep=";")
  
  
  cat("KEY_collapse_0\n")
  cat(str(KEY_collapse))
  cat("\n")
  
  ######### Merge ACTIVE TILES KEYED CLASS with KEY collapse----------
  
  KEY_collapse<-merge(KEY_collapse,
                      ACTIVE_TILES_Fq_enhancer_Fq,
                      by=c("KEY_Plus_carried_variants","Cell_Type"),
                      all=T)
  
  KEY_collapse$enhancer_CLASS_TILES[is.na(KEY_collapse$enhancer_CLASS_TILES)]<-0
  KEY_collapse$enhancer_CLASS[is.na(KEY_collapse$enhancer_CLASS)]<-"0"
  
  KEY_collapse$enhancer_CLASS<-factor(KEY_collapse$enhancer_CLASS,
                                        levels = c("0","1","2","3 or >"),
                                        ordered=T)
  
  
  cat("KEY_collapse_1\n")
  cat(str(KEY_collapse))
  cat("\n")
  
  KEY_collapse<-merge(KEY_collapse,
                      ACTIVE_TILES_Fq_ASE_Fq,
                      by=c("KEY_Plus_carried_variants","Cell_Type"),
                      all=T)
  
  KEY_collapse$ASE_CLASS_TILES[is.na(KEY_collapse$ASE_CLASS_TILES)]<-0
  KEY_collapse$ASE_CLASS[is.na(KEY_collapse$ASE_CLASS)]<-"0"
  
  KEY_collapse$ASE_CLASS<-factor(KEY_collapse$ASE_CLASS,
                                      levels = c("0","1","2","3 or >"),
                                      ordered=T)
  
  
  cat("KEY_collapse_1_5\n")
  cat(str(KEY_collapse))
  cat("\n")
  
  KEY_collapse<-merge(KEY_collapse,
                      ACTIVE_TILES_Fq_E_Plus_ASE_Fq,
                      by=c("KEY_Plus_carried_variants","Cell_Type"),
                      all=T)
  
  KEY_collapse$E_Plus_ASE_CLASS_TILES[is.na(KEY_collapse$E_Plus_ASE_CLASS_TILES)]<-0
  KEY_collapse$E_Plus_ASE_CLASS[is.na(KEY_collapse$E_Plus_ASE_CLASS)]<-"0"
  
  
  KEY_collapse$E_Plus_ASE_CLASS<-factor(KEY_collapse$E_Plus_ASE_CLASS,
                                        levels = c("0","1","2","3 or >"),
                                        ordered=T)
  
  
  cat("KEY_collapse_2\n")
  cat(str(KEY_collapse))
  cat("\n")
  
 
  
  KEY_collapse<-unique(KEY_collapse[,-which(colnames(KEY_collapse) == "factor4")])
  
  
  KEY_collapse_NCGR<-KEY_collapse[which(KEY_collapse$VEP_DEF_LABELS == "NCGR"),]
  
  cat("KEY_collapse_NCGR_0\n")
  cat(str(KEY_collapse_NCGR))
  cat("\n")
  cat(str(unique(KEY_collapse_NCGR$VAR)))
  cat("\n")
  cat(str(unique(KEY_collapse_NCGR$KEY_Plus_carried_variants[KEY_collapse_NCGR$Label_2 == "ASSAYED_VARIANT"])))
  cat("\n")
  cat(str(unique(KEY_collapse_NCGR$VAR[KEY_collapse_NCGR$Label_2 == "ASSAYED_VARIANT"])))
  cat("\n")
  cat(sprintf(as.character(names(summary(KEY_collapse_NCGR$VEP_DEF_LABELS)))))
  cat("\n")
  cat(sprintf(as.character(summary(KEY_collapse_NCGR$VEP_DEF_LABELS))))
  cat("\n")
  
  
  
  Element_collapse<-KEY_collapse
  
  Element_collapse$carried_variant_derived<-gsub("^chr","",Element_collapse$VAR)

  cat("Element_collapse_0\n")
  cat(str(Element_collapse))
  cat("\n")
  
  Element_collapse<-unique(Element_collapse[which(Element_collapse$carried_variant_derived == Element_collapse$carried_variants),])
  
  
  Element_collapse<-droplevels(Element_collapse)
  
  Element_collapse<-unique(Element_collapse[,-which(colnames(Element_collapse) == "carried_variant_derived")])
  
  
  cat("Element_collapse_1\n")
  cat(str(Element_collapse))
  cat("\n")
  cat(str(unique(Element_collapse$VAR)))
  cat("\n")
  cat(str(unique(Element_collapse$KEY_Plus_carried_variants[Element_collapse$Label_2 == "ASSAYED_VARIANT"])))
  cat("\n")
  cat(str(unique(Element_collapse$VAR[Element_collapse$Label_2 == "ASSAYED_VARIANT"])))
  cat("\n")
  cat(sprintf(as.character(names(summary(Element_collapse$VEP_DEF_LABELS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Element_collapse$VEP_DEF_LABELS))))
  cat("\n")
  
  Element_collapse<-rbind(Element_collapse,
                          KEY_collapse_NCGR)
  
  
  
  cat("Element_collapse_2\n")
  cat(str(Element_collapse))
  cat("\n")
  cat(str(unique(Element_collapse$VAR)))
  cat("\n")
  cat(str(unique(Element_collapse$KEY_Plus_carried_variants[Element_collapse$Label_2 == "ASSAYED_VARIANT"])))
  cat("\n")
  cat(str(unique(Element_collapse$VAR[Element_collapse$Label_2 == "ASSAYED_VARIANT"])))
  cat("\n")
  cat(sprintf(as.character(names(summary(Element_collapse$VEP_DEF_LABELS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Element_collapse$VEP_DEF_LABELS))))
  cat("\n")
  
  setwd(out)
  
  write.table(KEY_collapse,file="test_KEY_collapse.tsv",sep="\t",quote=F,row.names=F)
  
  
  ### The condition to element collapse eliminates NCGR ----> recover NCGR
  
  
  # ########################################################################################################################################################
  # quit(status = 1)
#   
  #### Cell_Type colors ----
  
  
  Cell_Type_levels<-levels(KEY_collapse$Cell_Type)
  colors_Cell_Type_levels<-c('#32A852','#1877C9','#553B68','#D45E85','#6DB2EE','#62D07F','#C9244B','#87447B','#D6B8E6')
  
  df.color_Cell_Type<-as.data.frame(cbind(Cell_Type_levels,colors_Cell_Type_levels[1:length(Cell_Type_levels)]), stringsAsFactors=F)
  colnames(df.color_Cell_Type)<-c("Cell_Type","colors")
  
  df.color_Cell_Type$Cell_Type<-factor(df.color_Cell_Type$Cell_Type,
                                       levels=levels(KEY_collapse$Cell_Type),
                                       ordered=T)
  cat("df.color_Cell_Type_\n")
  cat(str(df.color_Cell_Type))
  cat("\n")
  
  
  ######### SAVE THE COLLAPSED KEY -----
  
  setwd(out)
  
  filename<-paste('Element_collapse','.rds',sep='')
  
  
  saveRDS(Element_collapse,file=filename)
  
  filename<-paste('KEY_collapse','.rds',sep='')
  
  
  saveRDS(KEY_collapse,file=filename)
  
  
  filename<-paste('df_Cell_colors','.rds',sep='')
  
  
  saveRDS(df.color_Cell_Type,file=filename)
 
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
    make_option(c("--MPRA_Real_tile_QC2_PASS"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--CSQ_colors"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ACTIVE_TILES"), type="character", default=NULL, 
                metavar="type", 
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
  
  
  CLASSIFICATION_AT(opt)
 
  
}
  
  
  
 

###########################################################################

system.time( main() )
