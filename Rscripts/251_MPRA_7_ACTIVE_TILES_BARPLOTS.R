
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

options(warn=1)


Tier_stacked_barplots = function(option_list)
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
  
  ##### enhancer stacked barplot ----
  
 
  path5<-paste(out,'Tier_plots','/', sep='')
  
  cat("path5\n")
  cat(sprintf(as.character(path5)))
  cat("\n")
  
  
  if (file.exists(path5)){
    
    
    
    
  } else {
    dir.create(file.path(path5))
    
  }
  

  
  
  
  KEY_collapse_enhancer_CLASS.dt<-data.table(KEY_collapse, key=c("Cell_Type","enhancer_CLASS"))
  
  
  KEY_collapse_enhancer_Fq<-as.data.frame(KEY_collapse_enhancer_CLASS.dt[,.N,by=key(KEY_collapse_enhancer_CLASS.dt)], stringsAsFactors=F)
  
  colnames(KEY_collapse_enhancer_Fq)[which(colnames(KEY_collapse_enhancer_Fq) == "N")]<-"instances"
  
  cat("KEY_collapse_enhancer_Fq_\n")
  cat(str(KEY_collapse_enhancer_Fq))
  cat("\n")
  
  
  
  KEY_collapse_enhancer_CLASS_TOTAL.dt<-data.table(KEY_collapse, key=c("Cell_Type"))
  
  
  KEY_collapse_enhancer_CLASS_TOTAL<-as.data.frame(KEY_collapse_enhancer_CLASS_TOTAL.dt[,.N,by=key(KEY_collapse_enhancer_CLASS_TOTAL.dt)], stringsAsFactors=F)
  
  colnames(KEY_collapse_enhancer_CLASS_TOTAL)[which(colnames(KEY_collapse_enhancer_CLASS_TOTAL) == "N")]<-"TOTAL"
  
  cat("KEY_collapse_enhancer_CLASS_TOTAL_\n")
  cat(str(KEY_collapse_enhancer_CLASS_TOTAL))
  cat("\n")
  
  KEY_collapse_enhancer_Fq<-merge(KEY_collapse_enhancer_Fq,
                                  KEY_collapse_enhancer_CLASS_TOTAL,
                                   by="Cell_Type")
  
  KEY_collapse_enhancer_Fq$Perc<-round(100*(KEY_collapse_enhancer_Fq$instances/KEY_collapse_enhancer_Fq$TOTAL),1)
  
  
  cat("KEY_collapse_enhancer_Fq_\n")
  cat(str(KEY_collapse_enhancer_Fq))
  cat("\n")
  
  setwd(out)
  
  write.table(KEY_collapse_enhancer_Fq,file="test.tsv",sep="\t",quote=F,row.names=F)
  
  
  ##### stacked barplot
  
  levels_CT<-levels(KEY_collapse_enhancer_Fq$Cell_Type)
  
  KEY_collapse_enhancer_Fq$Cell_Type<-factor(KEY_collapse_enhancer_Fq$Cell_Type,
                                             levels=rev(levels_CT),
                                             ordered=T)
  
  breaks.Rank<-seq(0,110,by=10)
  labels.Rank<-as.character(breaks.Rank)
  #   
  
  
  graph<-KEY_collapse_enhancer_Fq %>%
    mutate(myaxis = paste0(Cell_Type, "\n", "n=", TOTAL)) %>%
    mutate(myaxis=fct_reorder(myaxis,as.numeric(Cell_Type))) %>%
    ggplot(aes(x=Perc, y=myaxis, fill=enhancer_CLASS)) +
    geom_bar(stat="identity",colour='black')+
    theme_bw()+
    theme(axis.title.y=element_text(size=24, family="sans"),
          axis.title.x=element_text(size=24, family="sans"),
          axis.text.y=element_text(angle=0,size=18, color=df_Cell_colors$color, family="sans"),
          axis.text.x=element_text(angle=0,size=18, color="black", family="sans"),
          legend.title=element_text(size=16,color="black", family="sans"),
          legend.text=element_text(size=12,color="black", family="sans"))+
    scale_y_discrete(name=NULL, drop=F)+
    scale_x_continuous(name="Percentage",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
    scale_fill_manual(values=c("white",'#32A852','#6DB2EE','#62D07F','#1877C9','#C9244B','#D45E85'),drop=F)+
    theme(legend.position="bottom")+
    guides(fill=guide_legend(title=paste("ACTIVE","enhancer TILES","per variant assayed", sep="\n")))+
    ggeasy::easy_center_title()
  
  # graph<-graph +coord_flip()
  
  path6<-paste(out,'Tier_plots','/','WITH_CTRLS','/', sep='')
  
  cat("path6\n")
  cat(sprintf(as.character(path6)))
  cat("\n")
  
  
  if (file.exists(path6)){
    
    
    
    
  } else {
    dir.create(file.path(path6))
    
  }
  
  
  
  setwd(path6)
  
  svglite(paste('enhancer_ACTIVE','.svg',sep=''), width = 8, height = 8)
  print(graph)
  dev.off()
  
  #### without CTRLS ----
  

  KEY_collapse_EXCLUDED_other_labels<-KEY_collapse[which(KEY_collapse$Label_2 == "ASSAYED_VARIANT"),]

  KEY_collapse_EXCLUDED_other_labels<-droplevels(KEY_collapse_EXCLUDED_other_labels)


  cat("KEY_collapse_EXCLUDED_other_labels\n")
  cat(str(KEY_collapse_EXCLUDED_other_labels))
  cat("\n")
  
  
  
  KEY_collapse_EXCLUDED_other_labels_enhancer_CLASS.dt<-data.table(KEY_collapse_EXCLUDED_other_labels, key=c("Cell_Type","enhancer_CLASS"))
  
  
  KEY_collapse_EXCLUDED_other_labels_enhancer_Fq<-as.data.frame(KEY_collapse_EXCLUDED_other_labels_enhancer_CLASS.dt[,.N,by=key(KEY_collapse_EXCLUDED_other_labels_enhancer_CLASS.dt)], stringsAsFactors=F)
  
  colnames(KEY_collapse_EXCLUDED_other_labels_enhancer_Fq)[which(colnames(KEY_collapse_EXCLUDED_other_labels_enhancer_Fq) == "N")]<-"instances"
  
  cat("KEY_collapse_EXCLUDED_other_labels_enhancer_Fq_\n")
  cat(str(KEY_collapse_EXCLUDED_other_labels_enhancer_Fq))
  cat("\n")
  
  
  
  KEY_collapse_EXCLUDED_other_labels_enhancer_CLASS_TOTAL.dt<-data.table(KEY_collapse_EXCLUDED_other_labels, key=c("Cell_Type"))
  
  
  KEY_collapse_EXCLUDED_other_labels_enhancer_CLASS_TOTAL<-as.data.frame(KEY_collapse_EXCLUDED_other_labels_enhancer_CLASS_TOTAL.dt[,.N,by=key(KEY_collapse_EXCLUDED_other_labels_enhancer_CLASS_TOTAL.dt)], stringsAsFactors=F)
  
  colnames(KEY_collapse_EXCLUDED_other_labels_enhancer_CLASS_TOTAL)[which(colnames(KEY_collapse_EXCLUDED_other_labels_enhancer_CLASS_TOTAL) == "N")]<-"TOTAL"
  
  cat("KEY_collapse_EXCLUDED_other_labels_enhancer_CLASS_TOTAL_\n")
  cat(str(KEY_collapse_EXCLUDED_other_labels_enhancer_CLASS_TOTAL))
  cat("\n")
  
  KEY_collapse_EXCLUDED_other_labels_enhancer_Fq<-merge(KEY_collapse_EXCLUDED_other_labels_enhancer_Fq,
                                  KEY_collapse_EXCLUDED_other_labels_enhancer_CLASS_TOTAL,
                                  by="Cell_Type")
  
  KEY_collapse_EXCLUDED_other_labels_enhancer_Fq$Perc<-round(100*(KEY_collapse_EXCLUDED_other_labels_enhancer_Fq$instances/KEY_collapse_EXCLUDED_other_labels_enhancer_Fq$TOTAL),1)
  
  
  cat("KEY_collapse_EXCLUDED_other_labels_enhancer_Fq_\n")
  cat(str(KEY_collapse_EXCLUDED_other_labels_enhancer_Fq))
  cat("\n")
  
  setwd(out)
  
  write.table(KEY_collapse_EXCLUDED_other_labels_enhancer_Fq,file="test.tsv",sep="\t",quote=F,row.names=F)
  
  
  ##### stacked barplot
  
  levels_CT<-levels(KEY_collapse_EXCLUDED_other_labels_enhancer_Fq$Cell_Type)
  
  KEY_collapse_EXCLUDED_other_labels_enhancer_Fq$Cell_Type<-factor(KEY_collapse_EXCLUDED_other_labels_enhancer_Fq$Cell_Type,
                                             levels=rev(levels_CT),
                                             ordered=T)
  
  breaks.Rank<-seq(0,110,by=10)
  labels.Rank<-as.character(breaks.Rank)
  #   
  
  
  graph<-KEY_collapse_EXCLUDED_other_labels_enhancer_Fq %>%
    mutate(myaxis = paste0(Cell_Type, "\n", "n=", TOTAL)) %>%
    mutate(myaxis=fct_reorder(myaxis,as.numeric(Cell_Type))) %>%
    ggplot(aes(x=Perc, y=myaxis, fill=enhancer_CLASS)) +
    geom_bar(stat="identity",colour='black')+
    theme_bw()+
    theme(axis.title.y=element_text(size=24, family="sans"),
          axis.title.x=element_text(size=24, family="sans"),
          axis.text.y=element_text(angle=0,size=18, color=df_Cell_colors$color, family="sans"),
          axis.text.x=element_text(angle=0,size=18, color="black", family="sans"),
          legend.title=element_text(size=16,color="black", family="sans"),
          legend.text=element_text(size=12,color="black", family="sans"))+
    scale_y_discrete(name=NULL, drop=F)+
    scale_x_continuous(name="Percentage",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
    scale_fill_manual(values=c("white",'#32A852','#6DB2EE','#62D07F','#1877C9','#C9244B','#D45E85'),drop=F)+
    theme(legend.position="bottom")+
    guides(fill=guide_legend(title=paste("ACTIVE","enhancer TILES","per variant assayed", sep="\n")))+
    ggeasy::easy_center_title()
  
  # graph<-graph +coord_flip()
  
  path6<-paste(out,'Tier_plots','/','WITHOUT_CTRLS','/', sep='')
  
  cat("path6\n")
  cat(sprintf(as.character(path6)))
  cat("\n")
  
  
  if (file.exists(path6)){
    
    
    
    
  } else {
    dir.create(file.path(path6))
    
  }
  
  
  
  setwd(path6)
  
  svglite(paste('enhancer_ACTIVE','.svg',sep=''), width = 8, height = 8)
  print(graph)
  dev.off()
  
  
  

  ##### E_Plus_ASE stacked barplot ----
  
  
  path5<-paste(out,'Tier_plots','/', sep='')
  
  cat("path5\n")
  cat(sprintf(as.character(path5)))
  cat("\n")
  
  
  if (file.exists(path5)){
    
    
    
    
  } else {
    dir.create(file.path(path5))
    
  }
  
  
  
  
  
  KEY_collapse_E_Plus_ASE_CLASS.dt<-data.table(KEY_collapse, key=c("Cell_Type","E_Plus_ASE_CLASS"))
  
  
  KEY_collapse_E_Plus_ASE_Fq<-as.data.frame(KEY_collapse_E_Plus_ASE_CLASS.dt[,.N,by=key(KEY_collapse_E_Plus_ASE_CLASS.dt)], stringsAsFactors=F)
  
  colnames(KEY_collapse_E_Plus_ASE_Fq)[which(colnames(KEY_collapse_E_Plus_ASE_Fq) == "N")]<-"instances"
  
  cat("KEY_collapse_E_Plus_ASE_Fq_\n")
  cat(str(KEY_collapse_E_Plus_ASE_Fq))
  cat("\n")
  
  
  
  KEY_collapse_E_Plus_ASE_CLASS_TOTAL.dt<-data.table(KEY_collapse, key=c("Cell_Type"))
  
  
  KEY_collapse_E_Plus_ASE_CLASS_TOTAL<-as.data.frame(KEY_collapse_E_Plus_ASE_CLASS_TOTAL.dt[,.N,by=key(KEY_collapse_E_Plus_ASE_CLASS_TOTAL.dt)], stringsAsFactors=F)
  
  colnames(KEY_collapse_E_Plus_ASE_CLASS_TOTAL)[which(colnames(KEY_collapse_E_Plus_ASE_CLASS_TOTAL) == "N")]<-"TOTAL"
  
  cat("KEY_collapse_E_Plus_ASE_CLASS_TOTAL_\n")
  cat(str(KEY_collapse_E_Plus_ASE_CLASS_TOTAL))
  cat("\n")
  
  KEY_collapse_E_Plus_ASE_Fq<-merge(KEY_collapse_E_Plus_ASE_Fq,
                                  KEY_collapse_E_Plus_ASE_CLASS_TOTAL,
                                  by="Cell_Type")
  
  KEY_collapse_E_Plus_ASE_Fq$Perc<-round(100*(KEY_collapse_E_Plus_ASE_Fq$instances/KEY_collapse_E_Plus_ASE_Fq$TOTAL),1)
  
  
  cat("KEY_collapse_E_Plus_ASE_Fq_\n")
  cat(str(KEY_collapse_E_Plus_ASE_Fq))
  cat("\n")
  
  setwd(out)
  
  write.table(KEY_collapse_E_Plus_ASE_Fq,file="test.tsv",sep="\t",quote=F,row.names=F)
  
  
  ##### stacked barplot
  
  levels_CT<-levels(KEY_collapse_E_Plus_ASE_Fq$Cell_Type)
  
  KEY_collapse_E_Plus_ASE_Fq$Cell_Type<-factor(KEY_collapse_E_Plus_ASE_Fq$Cell_Type,
                                             levels=rev(levels_CT),
                                             ordered=T)
  
  breaks.Rank<-seq(0,110,by=10)
  labels.Rank<-as.character(breaks.Rank)
  #   
  
  
  graph<-KEY_collapse_E_Plus_ASE_Fq %>%
    mutate(myaxis = paste0(Cell_Type, "\n", "n=", TOTAL)) %>%
    mutate(myaxis=fct_reorder(myaxis,as.numeric(Cell_Type))) %>%
    ggplot(aes(x=Perc, y=myaxis, fill=E_Plus_ASE_CLASS)) +
    geom_bar(stat="identity",colour='black')+
    theme_bw()+
    theme(axis.title.y=element_text(size=24, family="sans"),
          axis.title.x=element_text(size=24, family="sans"),
          axis.text.y=element_text(angle=0,size=18, color=df_Cell_colors$color, family="sans"),
          axis.text.x=element_text(angle=0,size=18, color="black", family="sans"),
          legend.title=element_text(size=16,color="black", family="sans"),
          legend.text=element_text(size=12,color="black", family="sans"))+
    scale_y_discrete(name=NULL, drop=F)+
    scale_x_continuous(name="Percentage",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
    scale_fill_manual(values=c("white",'#32A852','#6DB2EE','#62D07F','#1877C9','#C9244B','#D45E85'),drop=F)+
    theme(legend.position="bottom")+
    guides(fill=guide_legend(title=paste("ACTIVE","E_Plus_ASE TILES","per variant assayed", sep="\n")))+
    ggeasy::easy_center_title()
  
  # graph<-graph +coord_flip()
  
  path6<-paste(out,'Tier_plots','/','WITH_CTRLS','/', sep='')
  
  cat("path6\n")
  cat(sprintf(as.character(path6)))
  cat("\n")
  
  
  if (file.exists(path6)){
    
    
    
    
  } else {
    dir.create(file.path(path6))
    
  }
  
  
  
  setwd(path6)
  
  svglite(paste('E_Plus_ASE_ACTIVE','.svg',sep=''), width = 8, height = 8)
  print(graph)
  dev.off()
  
  #### without CTRLS ----
  
  
  KEY_collapse_EXCLUDED_other_labels<-KEY_collapse[which(KEY_collapse$Label_2 == "ASSAYED_VARIANT"),]
  
  KEY_collapse_EXCLUDED_other_labels<-droplevels(KEY_collapse_EXCLUDED_other_labels)
  
  
  cat("KEY_collapse_EXCLUDED_other_labels\n")
  cat(str(KEY_collapse_EXCLUDED_other_labels))
  cat("\n")
  
  
  
  KEY_collapse_EXCLUDED_other_labels_E_Plus_ASE_CLASS.dt<-data.table(KEY_collapse_EXCLUDED_other_labels, key=c("Cell_Type","E_Plus_ASE_CLASS"))
  
  
  KEY_collapse_EXCLUDED_other_labels_E_Plus_ASE_Fq<-as.data.frame(KEY_collapse_EXCLUDED_other_labels_E_Plus_ASE_CLASS.dt[,.N,by=key(KEY_collapse_EXCLUDED_other_labels_E_Plus_ASE_CLASS.dt)], stringsAsFactors=F)
  
  colnames(KEY_collapse_EXCLUDED_other_labels_E_Plus_ASE_Fq)[which(colnames(KEY_collapse_EXCLUDED_other_labels_E_Plus_ASE_Fq) == "N")]<-"instances"
  
  cat("KEY_collapse_EXCLUDED_other_labels_E_Plus_ASE_Fq_\n")
  cat(str(KEY_collapse_EXCLUDED_other_labels_E_Plus_ASE_Fq))
  cat("\n")
  
  
  
  KEY_collapse_EXCLUDED_other_labels_E_Plus_ASE_CLASS_TOTAL.dt<-data.table(KEY_collapse_EXCLUDED_other_labels, key=c("Cell_Type"))
  
  
  KEY_collapse_EXCLUDED_other_labels_E_Plus_ASE_CLASS_TOTAL<-as.data.frame(KEY_collapse_EXCLUDED_other_labels_E_Plus_ASE_CLASS_TOTAL.dt[,.N,by=key(KEY_collapse_EXCLUDED_other_labels_E_Plus_ASE_CLASS_TOTAL.dt)], stringsAsFactors=F)
  
  colnames(KEY_collapse_EXCLUDED_other_labels_E_Plus_ASE_CLASS_TOTAL)[which(colnames(KEY_collapse_EXCLUDED_other_labels_E_Plus_ASE_CLASS_TOTAL) == "N")]<-"TOTAL"
  
  cat("KEY_collapse_EXCLUDED_other_labels_E_Plus_ASE_CLASS_TOTAL_\n")
  cat(str(KEY_collapse_EXCLUDED_other_labels_E_Plus_ASE_CLASS_TOTAL))
  cat("\n")
  
  KEY_collapse_EXCLUDED_other_labels_E_Plus_ASE_Fq<-merge(KEY_collapse_EXCLUDED_other_labels_E_Plus_ASE_Fq,
                                                        KEY_collapse_EXCLUDED_other_labels_E_Plus_ASE_CLASS_TOTAL,
                                                        by="Cell_Type")
  
  KEY_collapse_EXCLUDED_other_labels_E_Plus_ASE_Fq$Perc<-round(100*(KEY_collapse_EXCLUDED_other_labels_E_Plus_ASE_Fq$instances/KEY_collapse_EXCLUDED_other_labels_E_Plus_ASE_Fq$TOTAL),1)
  
  
  cat("KEY_collapse_EXCLUDED_other_labels_E_Plus_ASE_Fq_\n")
  cat(str(KEY_collapse_EXCLUDED_other_labels_E_Plus_ASE_Fq))
  cat("\n")
  
  setwd(out)
  
  write.table(KEY_collapse_EXCLUDED_other_labels_E_Plus_ASE_Fq,file="test.tsv",sep="\t",quote=F,row.names=F)
  
  
  ##### stacked barplot
  
  levels_CT<-levels(KEY_collapse_EXCLUDED_other_labels_E_Plus_ASE_Fq$Cell_Type)
  
  KEY_collapse_EXCLUDED_other_labels_E_Plus_ASE_Fq$Cell_Type<-factor(KEY_collapse_EXCLUDED_other_labels_E_Plus_ASE_Fq$Cell_Type,
                                                                   levels=rev(levels_CT),
                                                                   ordered=T)
  
  breaks.Rank<-seq(0,110,by=10)
  labels.Rank<-as.character(breaks.Rank)
  #   
  
  
  graph<-KEY_collapse_EXCLUDED_other_labels_E_Plus_ASE_Fq %>%
    mutate(myaxis = paste0(Cell_Type, "\n", "n=", TOTAL)) %>%
    mutate(myaxis=fct_reorder(myaxis,as.numeric(Cell_Type))) %>%
    ggplot(aes(x=Perc, y=myaxis, fill=E_Plus_ASE_CLASS)) +
    geom_bar(stat="identity",colour='black')+
    theme_bw()+
    theme(axis.title.y=element_text(size=24, family="sans"),
          axis.title.x=element_text(size=24, family="sans"),
          axis.text.y=element_text(angle=0,size=18, color=df_Cell_colors$color, family="sans"),
          axis.text.x=element_text(angle=0,size=18, color="black", family="sans"),
          legend.title=element_text(size=16,color="black", family="sans"),
          legend.text=element_text(size=12,color="black", family="sans"))+
    scale_y_discrete(name=NULL, drop=F)+
    scale_x_continuous(name="Percentage",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
    scale_fill_manual(values=c("white",'#32A852','#6DB2EE','#62D07F','#1877C9','#C9244B','#D45E85'),drop=F)+
    theme(legend.position="bottom")+
    guides(fill=guide_legend(title=paste("ACTIVE","E_Plus_ASE TILES","per variant assayed", sep="\n")))+
    ggeasy::easy_center_title()
  
  # graph<-graph +coord_flip()
  
  path6<-paste(out,'Tier_plots','/','WITHOUT_CTRLS','/', sep='')
  
  cat("path6\n")
  cat(sprintf(as.character(path6)))
  cat("\n")
  
  
  if (file.exists(path6)){
    
    
    
    
  } else {
    dir.create(file.path(path6))
    
  }
  
  
  
  setwd(path6)
  
  svglite(paste('E_Plus_ASE_ACTIVE','.svg',sep=''), width = 8, height = 8)
  print(graph)
  dev.off()
  
  
  
  
}

VEP_LABEL_ER_stacked_barplots_enhancer_CLASS = function(option_list)
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
  cat(sprintf(as.character(df_CSQ_colors$VEP_DEF_LABELS)))
  cat("\n")
  cat(sprintf(as.character(df_CSQ_colors$color)))
  cat("\n")
  
  
  
  
  
  ##### enhancer stacked barplot ----
  
  
  path5<-paste(out,'Tier_plots_ER','/', sep='')
  
  cat("path5\n")
  cat(sprintf(as.character(path5)))
  cat("\n")
  
  
  if (file.exists(path5)){
    
    
    
    
  } else {
    dir.create(file.path(path5))
    
  }
  
  
  
  
  
  KEY_collapse_enhancer_CLASS.dt<-data.table(KEY_collapse, key=c("Cell_Type","enhancer_CLASS","VEP_DEF_LABELS"))
  
  
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
    
    AT_LEAST_COLLAPSE.dt<-data.table(KEY_collapse_enhancer_Fq_sel, key=c("Cell_Type","VEP_DEF_LABELS"))
    
    
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
  
  
  # setwd(out)
  # 
  # write.table(file="test.tsv",TABLE_AT_LEAST_enhancer,quote=F, sep="\t", row.names = F)
  # 
  # ####################################################################################################################################
  # quit(status = 1)
  
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
    
    #### graph ----
    
    cat("TABLE_AT_LEAST_enhancer_Cell_Type_sel_0\n")
    cat(str(TABLE_AT_LEAST_enhancer_Cell_Type_sel))
    cat("\n")
    cat(sprintf(as.character(names(summary(TABLE_AT_LEAST_enhancer_Cell_Type_sel$VEP_DEF_LABELS)))))
    cat("\n")
    cat(sprintf(as.character(summary(TABLE_AT_LEAST_enhancer_Cell_Type_sel$VEP_DEF_LABELS))))
    cat("\n")
    
    
    df_CSQ_colors_sel<-df_CSQ_colors[which(df_CSQ_colors$VEP_DEF_LABELS%in%TABLE_AT_LEAST_enhancer_Cell_Type_sel$VEP_DEF_LABELS),]
    
    df_CSQ_colors_sel<-droplevels(df_CSQ_colors_sel)
    
    df_CSQ_colors_sel<-df_CSQ_colors_sel[order(df_CSQ_colors_sel$VEP_DEF_LABELS),]
    
    cat("df_CSQ_colors_sel_0\n")
    cat(str(df_CSQ_colors_sel))
    cat("\n")
    
    
    # df_CSQ_colors_sel$color<-factor(df_CSQ_colors_sel$color,
    #                                 levels=unique(as.character(df_CSQ_colors_sel$color)),
    #                                 ordered=T)
    
    
    cat("df_CSQ_colors_sel_1\n")
    cat(str(df_CSQ_colors_sel))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_CSQ_colors_sel$VEP_DEF_LABELS)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_CSQ_colors_sel$VEP_DEF_LABELS))))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(df_CSQ_colors_sel$color))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(df_CSQ_colors_sel$color)))))
    cat("\n")
    
    
    # write.table(df_CSQ_colors_sel,file="test_color.tsv",sep="\t",quote=F,row.names = F)
    
    
    
    breaks.y<-seq(0,110,by=10)
    labels.y<-as.character(breaks.y)
    
    
    
    
    
    
    graph<-TABLE_AT_LEAST_enhancer_Cell_Type_sel%>%
      mutate(myaxis = paste0(AT_LEAST_category, "\n", "n=", AT_LEAST_TOTAL)) %>%
      mutate(myaxis=fct_reorder(myaxis,as.numeric(AT_LEAST_category))) %>%
      ggplot(aes(x=myaxis,
                 y=Perc,
                 fill=VEP_DEF_LABELS)) +
      geom_bar(stat="identity",colour='black')+
      theme_bw()+
      theme(axis.title.y=element_text(size=18, family="sans"),
            axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
            axis.text.x=element_text(angle=0,size=18, color="black", family="sans"))+
      scale_y_continuous(name=paste("Percentage"),breaks=breaks.y,labels=labels.y,
                         limits=c(breaks.y[1],breaks.y[length(breaks.y)]))+
      scale_x_discrete(name=NULL, drop=F)+
      scale_fill_manual(breaks=df_CSQ_colors_sel$VEP_DEF_LABELS,
                        values=df_CSQ_colors_sel$color, drop=F)+
      theme(legend.position="right")+
      ggeasy::easy_center_title()
    
    
    setwd(out)
    
    svgname<-paste(paste("test","_",Cell_Type_array_sel,".svg",sep=''))
    makesvg = TRUE
    
    if (makesvg == TRUE)
    {
      ggsave(svgname, plot= graph,
             device="svg",
             height=10, width=12)
    }
    
    list_DEF[[Cell_Type_array_sel]]<-graph
    
    
    
    
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
          
          levels_CSQ<-levels(TABLE_AT_LEAST_enhancer$VEP_DEF_LABELS)
          
          list_CSQ<-list()
          
          for(iteration_levels_CSQ in 1:length(levels_CSQ))
          {
            
            
            levels_CSQ_sel<-levels_CSQ[iteration_levels_CSQ]
            
            cat("->\t")
            cat(sprintf(as.character(levels_CSQ_sel)))
            cat("\n")
            
            
            TABLE_AT_LEAST_enhancer_Cell_Type_sel_0_CSQ_sel<-TABLE_AT_LEAST_enhancer_Cell_Type_sel_0[which(TABLE_AT_LEAST_enhancer_Cell_Type_sel_0$VEP_DEF_LABELS == levels_CSQ_sel),]
            
            # cat("TABLE_AT_LEAST_enhancer_Cell_Type_sel_0_CSQ_sel\n")
            # cat(str(TABLE_AT_LEAST_enhancer_Cell_Type_sel_0_CSQ_sel))
            # cat("\n")
            
            TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_CSQ_sel<-TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel[which(TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel$VEP_DEF_LABELS == levels_CSQ_sel),]
            
            # cat("TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_CSQ_sel\n")
            # cat(str(TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_CSQ_sel))
            # cat("\n")
            
            
            TABLE_AT_LEAST_enhancer_Cell_Type_sel_0_NOT_CSQ_sel<-TABLE_AT_LEAST_enhancer_Cell_Type_sel_0[which(TABLE_AT_LEAST_enhancer_Cell_Type_sel_0$VEP_DEF_LABELS != levels_CSQ_sel),]
            
            # cat("TABLE_AT_LEAST_enhancer_Cell_Type_sel_0_NOT_CSQ_sel\n")
            # cat(str(TABLE_AT_LEAST_enhancer_Cell_Type_sel_0_NOT_CSQ_sel))
            # cat("\n")
            
            TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_NOT_CSQ_sel<-TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel[which(TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel$VEP_DEF_LABELS != levels_CSQ_sel),]
            
            # cat("TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_NOT_CSQ_sel\n")
            # cat(str(TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_NOT_CSQ_sel))
            # cat("\n")
            
            
            if(dim(TABLE_AT_LEAST_enhancer_Cell_Type_sel_0_CSQ_sel)[1]>0)
            {
              if(dim(TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_CSQ_sel)[1]>0)
              {
                cat("->Hello_world_1\n")
                
                
                A.df<-as.data.frame(rbind(cbind(sum(TABLE_AT_LEAST_enhancer_Cell_Type_sel_0_NOT_CSQ_sel$AT_LEAST_instances),sum(TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_NOT_CSQ_sel$AT_LEAST_instances)),
                                          cbind(TABLE_AT_LEAST_enhancer_Cell_Type_sel_0_CSQ_sel$AT_LEAST_instances,TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_CSQ_sel$AT_LEAST_instances)),
                                    stringsAsFactors=F)
                
                colnames(A.df)<-c("0",level_sel)
                row.names(A.df)<-c(paste("NOT",levels_CSQ_sel,sep='_'),
                                   levels_CSQ_sel)
                
                
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
                
                a.df<-as.data.frame(cbind(level_sel,levels_CSQ_sel,log_pval), stringsAsFactors=F)
                
                colnames(a.df)<-c("Category","COMPARISON_category","log_pval")
                
                # cat("a.df\n")
                # cat(str(a.df))
                # cat("\n")
                
                list_CSQ[[iteration_levels_CSQ]]<-a.df
                
                
              }else{
                cat("->Hello_world_2\n")
                
                
                A.df<-as.data.frame(rbind(cbind(sum(TABLE_AT_LEAST_enhancer_Cell_Type_sel_0_NOT_CSQ_sel$AT_LEAST_instances),sum(TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_NOT_CSQ_sel$AT_LEAST_instances)),
                                          cbind(TABLE_AT_LEAST_enhancer_Cell_Type_sel_0_CSQ_sel$AT_LEAST_instances,0)),
                                    stringsAsFactors=F)
                
                colnames(A.df)<-c("0",level_sel)
                row.names(A.df)<-c(paste("NOT",levels_CSQ_sel,sep='_'),
                                   levels_CSQ_sel)
                
                
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
                
                a.df<-as.data.frame(cbind(level_sel,levels_CSQ_sel,log_pval), stringsAsFactors=F)
                
                colnames(a.df)<-c("Category","COMPARISON_category","log_pval")
                
                # cat("a.df\n")
                # cat(str(a.df))
                # cat("\n")
                
                list_CSQ[[iteration_levels_CSQ]]<-a.df
                
                
              }# dim(TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_CSQ_sel)[1]>0
              
            }else{
              
              if(dim(TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_CSQ_sel)[1]>0)
              {
                
                cat("->Hello_world_3\n")
                
                A.df<-as.data.frame(rbind(cbind(sum(TABLE_AT_LEAST_enhancer_Cell_Type_sel_0_NOT_CSQ_sel$AT_LEAST_instances),sum(TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_NOT_CSQ_sel$AT_LEAST_instances)),
                                          cbind(0,TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_CSQ_sel$AT_LEAST_instances)),
                                    stringsAsFactors=F)
                
                colnames(A.df)<-c("0",level_sel)
                row.names(A.df)<-c(paste("NOT",levels_CSQ_sel,sep='_'),
                                   levels_CSQ_sel)
                
                
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
                
                a.df<-as.data.frame(cbind(level_sel,levels_CSQ_sel,log_pval), stringsAsFactors=F)
                
                colnames(a.df)<-c("Category","COMPARISON_category","log_pval")
                
                # cat("a.df\n")
                # cat(str(a.df))
                # cat("\n")
                
                list_CSQ[[iteration_levels_CSQ]]<-a.df
              }else{
                
                cat("->Hello_world_4\n")
                
                A.df<-as.data.frame(rbind(cbind(sum(TABLE_AT_LEAST_enhancer_Cell_Type_sel_0_NOT_CSQ_sel$AT_LEAST_instances),sum(TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_NOT_CSQ_sel$AT_LEAST_instances)),
                                          cbind(0,0)),
                                    stringsAsFactors=F)
                
                colnames(A.df)<-c("0",level_sel)
                row.names(A.df)<-c(paste("NOT",levels_CSQ_sel,sep='_'),
                                   levels_CSQ_sel)
                
                
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
                
                a.df<-as.data.frame(cbind(level_sel,levels_CSQ_sel,log_pval), stringsAsFactors=F)
                
                colnames(a.df)<-c("Category","COMPARISON_category","log_pval")
                
                # cat("a.df\n")
                # cat(str(a.df))
                # cat("\n")
                
                list_CSQ[[iteration_levels_CSQ]]<-a.df
                
              }#dim(TABLE_AT_LEAST_enhancer_Cell_Type_sel_sel_CSQ_sel)[1]>0
              
              
              
            }# dim(TABLE_AT_LEAST_enhancer_Cell_Type_sel_0_CSQ_sel)[1]>0
            
            
          }# iteration_levels_CSQ
          
          CHI_Sq_CSQ = unique(as.data.frame(data.table::rbindlist(list_CSQ, fill = T)))
          
          # cat("CHI_Sq_CSQ\n")
          # cat(str(CHI_Sq_CSQ))
          # cat("\n")
          
          list_ChiSq_category[[k]]<-CHI_Sq_CSQ
          
          # ############################################################################################################
          # quit(status = 1) 
          
          
        }#level_sel != "0"
      }# k
    }# dim(TABLE_AT_LEAST_enhancer_Cell_Type_sel)[1] > 0
    
    CHI_Sq_CSQ_CATEGORY = unique(as.data.frame(data.table::rbindlist(list_ChiSq_category, fill = T)))
    
    CHI_Sq_CSQ_CATEGORY$Cell_Type<-Cell_Type_array_sel
    
    cat("CHI_Sq_CSQ_CATEGORY\n")
    cat(str(CHI_Sq_CSQ_CATEGORY))
    cat("\n")
    
    list_ChiSq_DEF[[i]]<-CHI_Sq_CSQ_CATEGORY
    
    # ############################################################################################################
    # quit(status = 1) 
    
  }# i Cell_Type
  
  CHI_Sq_DEF = unique(as.data.frame(data.table::rbindlist(list_ChiSq_DEF, fill = T)))
  
  
  cat("CHI_Sq_DEF_0\n")
  cat(str(CHI_Sq_DEF))
  cat("\n")
  
  
  
  CHI_Sq_DEF$COMPARISON_category<-factor(CHI_Sq_DEF$COMPARISON_category,
                                         levels=c("LOF","MISS","SYN","UTR5","UTR3",
                                                  "INTRON","INTERGENIC","UPSTREAM","DOWNSTREAM","REGULATORY",
                                                  "TFBS","SPLICE",
                                                  "OTHER","NMD","NCT",
                                                  "NCGR","ASE_CTRL","enhancer_CTRL","Kousik_variant"),
                                         ordered=T)
  
  CHI_Sq_DEF$Category<-factor(CHI_Sq_DEF$Category,
                              levels = c("0","AT LEAST 1","AT LEAST 2","AT LEAST 3"),
                              ordered=T)
  
  CHI_Sq_DEF$log_pval<-as.numeric(CHI_Sq_DEF$log_pval)
  
  CHI_Sq_DEF$Cell_Type<-factor(CHI_Sq_DEF$Cell_Type,
                               levels = c("K562","CHRF","HL60","THP1"),
                               ordered=T)
  
  
  cat("CHI_Sq_DEF_1\n")
  cat(str(CHI_Sq_DEF))
  cat("\n")
  
  #### SAVE----
  
  
  setwd(path5)
  
  saveRDS(list_DEF,file="worstCSQ_enhancer_CLASS_graphs.rds")
  
  saveRDS(CHI_Sq_DEF,file="worstCSQ_enhancer_CLASS_stats.rds")
  
}

VEP_LABEL_ER_stacked_barplots_E_Plus_ASE_CLASS = function(option_list)
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
  cat(sprintf(as.character(df_CSQ_colors$VEP_DEF_LABELS)))
  cat("\n")
  cat(sprintf(as.character(df_CSQ_colors$color)))
  cat("\n")
  
  
  
  
  
  ##### E_Plus_ASE stacked barplot ----
  
  
  path5<-paste(out,'Tier_plots_ER','/', sep='')
  
  cat("path5\n")
  cat(sprintf(as.character(path5)))
  cat("\n")
  
  
  if (file.exists(path5)){
    
    
    
    
  } else {
    dir.create(file.path(path5))
    
  }
  
  
  
  
  
  KEY_collapse_E_Plus_ASE_CLASS.dt<-data.table(KEY_collapse, key=c("Cell_Type","E_Plus_ASE_CLASS","VEP_DEF_LABELS"))
  
  
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
    
    AT_LEAST_COLLAPSE.dt<-data.table(KEY_collapse_E_Plus_ASE_Fq_sel, key=c("Cell_Type","VEP_DEF_LABELS"))
    
    
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
  
  
  # setwd(out)
  # 
  # write.table(file="test.tsv",TABLE_AT_LEAST_E_Plus_ASE,quote=F, sep="\t", row.names = F)
  # 
  # ####################################################################################################################################
  # quit(status = 1)
  
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
    
    #### graph ----
    
    cat("TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0\n")
    cat(str(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel))
    cat("\n")
    cat(sprintf(as.character(names(summary(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel$VEP_DEF_LABELS)))))
    cat("\n")
    cat(sprintf(as.character(summary(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel$VEP_DEF_LABELS))))
    cat("\n")
    
    
    df_CSQ_colors_sel<-df_CSQ_colors[which(df_CSQ_colors$VEP_DEF_LABELS%in%TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel$VEP_DEF_LABELS),]
    
    df_CSQ_colors_sel<-droplevels(df_CSQ_colors_sel)
    
    df_CSQ_colors_sel<-df_CSQ_colors_sel[order(df_CSQ_colors_sel$VEP_DEF_LABELS),]
    
    cat("df_CSQ_colors_sel_0\n")
    cat(str(df_CSQ_colors_sel))
    cat("\n")
    
    
    # df_CSQ_colors_sel$color<-factor(df_CSQ_colors_sel$color,
    #                                 levels=unique(as.character(df_CSQ_colors_sel$color)),
    #                                 ordered=T)
    
    
    cat("df_CSQ_colors_sel_1\n")
    cat(str(df_CSQ_colors_sel))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_CSQ_colors_sel$VEP_DEF_LABELS)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_CSQ_colors_sel$VEP_DEF_LABELS))))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(df_CSQ_colors_sel$color))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(df_CSQ_colors_sel$color)))))
    cat("\n")
    
    
    # write.table(df_CSQ_colors_sel,file="test_color.tsv",sep="\t",quote=F,row.names = F)
    
    
    
    breaks.y<-seq(0,110,by=10)
    labels.y<-as.character(breaks.y)
    
  
    
    
    
    
    graph<-TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel%>%
      mutate(myaxis = paste0(AT_LEAST_category, "\n", "n=", AT_LEAST_TOTAL)) %>%
      mutate(myaxis=fct_reorder(myaxis,as.numeric(AT_LEAST_category))) %>%
      ggplot(aes(x=myaxis,
                 y=Perc,
                 fill=VEP_DEF_LABELS)) +
      geom_bar(stat="identity",colour='black')+
      theme_bw()+
      theme(axis.title.y=element_text(size=18, family="sans"),
            axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
            axis.text.x=element_text(angle=0,size=18, color="black", family="sans"))+
      scale_y_continuous(name=paste("Percentage"),breaks=breaks.y,labels=labels.y,
                         limits=c(breaks.y[1],breaks.y[length(breaks.y)]))+
      scale_x_discrete(name=NULL, drop=F)+
      scale_fill_manual(breaks=df_CSQ_colors_sel$VEP_DEF_LABELS,
                        values=df_CSQ_colors_sel$color, drop=F)+
      theme(legend.position="right")+
      ggeasy::easy_center_title()
    
    
    setwd(out)
    
    svgname<-paste(paste("test","_",Cell_Type_array_sel,".svg",sep=''))
    makesvg = TRUE
    
    if (makesvg == TRUE)
    {
      ggsave(svgname, plot= graph,
             device="svg",
             height=10, width=12)
    }
    
    list_DEF[[Cell_Type_array_sel]]<-graph
    
    
   
    
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
          
          levels_CSQ<-levels(TABLE_AT_LEAST_E_Plus_ASE$VEP_DEF_LABELS)
          
          list_CSQ<-list()
          
          for(iteration_levels_CSQ in 1:length(levels_CSQ))
          {
            
            
            levels_CSQ_sel<-levels_CSQ[iteration_levels_CSQ]
            
            cat("->\t")
            cat(sprintf(as.character(levels_CSQ_sel)))
            cat("\n")
            
            
            TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0_CSQ_sel<-TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0[which(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0$VEP_DEF_LABELS == levels_CSQ_sel),]
            
            # cat("TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0_CSQ_sel\n")
            # cat(str(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0_CSQ_sel))
            # cat("\n")
            
            TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_CSQ_sel<-TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel[which(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel$VEP_DEF_LABELS == levels_CSQ_sel),]
            
            # cat("TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_CSQ_sel\n")
            # cat(str(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_CSQ_sel))
            # cat("\n")
            
            
            TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0_NOT_CSQ_sel<-TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0[which(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0$VEP_DEF_LABELS != levels_CSQ_sel),]
            
            # cat("TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0_NOT_CSQ_sel\n")
            # cat(str(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0_NOT_CSQ_sel))
            # cat("\n")
            
            TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_NOT_CSQ_sel<-TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel[which(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel$VEP_DEF_LABELS != levels_CSQ_sel),]
            
            # cat("TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_NOT_CSQ_sel\n")
            # cat(str(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_NOT_CSQ_sel))
            # cat("\n")
            
            
            if(dim(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0_CSQ_sel)[1]>0)
            {
              if(dim(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_CSQ_sel)[1]>0)
              {
                cat("->Hello_world_1\n")
                
                
                A.df<-as.data.frame(rbind(cbind(sum(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0_NOT_CSQ_sel$AT_LEAST_instances),sum(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_NOT_CSQ_sel$AT_LEAST_instances)),
                                          cbind(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0_CSQ_sel$AT_LEAST_instances,TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_CSQ_sel$AT_LEAST_instances)),
                                    stringsAsFactors=F)
                
                colnames(A.df)<-c("0",level_sel)
                row.names(A.df)<-c(paste("NOT",levels_CSQ_sel,sep='_'),
                                   levels_CSQ_sel)
                
                
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
                
                a.df<-as.data.frame(cbind(level_sel,levels_CSQ_sel,log_pval), stringsAsFactors=F)
                
                colnames(a.df)<-c("Category","COMPARISON_category","log_pval")
                
                # cat("a.df\n")
                # cat(str(a.df))
                # cat("\n")
                
                list_CSQ[[iteration_levels_CSQ]]<-a.df
                
                
              }else{
                cat("->Hello_world_2\n")
                
                
                A.df<-as.data.frame(rbind(cbind(sum(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0_NOT_CSQ_sel$AT_LEAST_instances),sum(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_NOT_CSQ_sel$AT_LEAST_instances)),
                                          cbind(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0_CSQ_sel$AT_LEAST_instances,0)),
                                    stringsAsFactors=F)
                
                colnames(A.df)<-c("0",level_sel)
                row.names(A.df)<-c(paste("NOT",levels_CSQ_sel,sep='_'),
                                   levels_CSQ_sel)
                
                
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
                
                a.df<-as.data.frame(cbind(level_sel,levels_CSQ_sel,log_pval), stringsAsFactors=F)
                
                colnames(a.df)<-c("Category","COMPARISON_category","log_pval")
                
                # cat("a.df\n")
                # cat(str(a.df))
                # cat("\n")
                
                list_CSQ[[iteration_levels_CSQ]]<-a.df
                
                
              }# dim(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_CSQ_sel)[1]>0
              
            }else{
              
              if(dim(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_CSQ_sel)[1]>0)
              {
                
                cat("->Hello_world_3\n")
                
                A.df<-as.data.frame(rbind(cbind(sum(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0_NOT_CSQ_sel$AT_LEAST_instances),sum(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_NOT_CSQ_sel$AT_LEAST_instances)),
                                          cbind(0,TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_CSQ_sel$AT_LEAST_instances)),
                                    stringsAsFactors=F)
                
                colnames(A.df)<-c("0",level_sel)
                row.names(A.df)<-c(paste("NOT",levels_CSQ_sel,sep='_'),
                                   levels_CSQ_sel)
                
                
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
                
                a.df<-as.data.frame(cbind(level_sel,levels_CSQ_sel,log_pval), stringsAsFactors=F)
                
                colnames(a.df)<-c("Category","COMPARISON_category","log_pval")
                
                # cat("a.df\n")
                # cat(str(a.df))
                # cat("\n")
                
                list_CSQ[[iteration_levels_CSQ]]<-a.df
              }else{
                
                cat("->Hello_world_4\n")
                
                A.df<-as.data.frame(rbind(cbind(sum(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0_NOT_CSQ_sel$AT_LEAST_instances),sum(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_NOT_CSQ_sel$AT_LEAST_instances)),
                                          cbind(0,0)),
                                    stringsAsFactors=F)
                
                colnames(A.df)<-c("0",level_sel)
                row.names(A.df)<-c(paste("NOT",levels_CSQ_sel,sep='_'),
                                   levels_CSQ_sel)
                
                
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
                
                a.df<-as.data.frame(cbind(level_sel,levels_CSQ_sel,log_pval), stringsAsFactors=F)
                
                colnames(a.df)<-c("Category","COMPARISON_category","log_pval")
                
                # cat("a.df\n")
                # cat(str(a.df))
                # cat("\n")
                
                list_CSQ[[iteration_levels_CSQ]]<-a.df
                  
              }#dim(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_sel_CSQ_sel)[1]>0
              
              
              
            }# dim(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel_0_CSQ_sel)[1]>0
            
            
          }# iteration_levels_CSQ
          
          CHI_Sq_CSQ = unique(as.data.frame(data.table::rbindlist(list_CSQ, fill = T)))
          
          # cat("CHI_Sq_CSQ\n")
          # cat(str(CHI_Sq_CSQ))
          # cat("\n")
          
          list_ChiSq_category[[k]]<-CHI_Sq_CSQ
          
          # ############################################################################################################
          # quit(status = 1) 
          
          
        }#level_sel != "0"
      }# k
    }# dim(TABLE_AT_LEAST_E_Plus_ASE_Cell_Type_sel)[1] > 0
    
    CHI_Sq_CSQ_CATEGORY = unique(as.data.frame(data.table::rbindlist(list_ChiSq_category, fill = T)))
    
    CHI_Sq_CSQ_CATEGORY$Cell_Type<-Cell_Type_array_sel
    
    cat("CHI_Sq_CSQ_CATEGORY\n")
    cat(str(CHI_Sq_CSQ_CATEGORY))
    cat("\n")
    
    list_ChiSq_DEF[[i]]<-CHI_Sq_CSQ_CATEGORY
    
    # ############################################################################################################
    # quit(status = 1) 
    
  }# i Cell_Type
  
  CHI_Sq_DEF = unique(as.data.frame(data.table::rbindlist(list_ChiSq_DEF, fill = T)))
  
  
  cat("CHI_Sq_DEF_0\n")
  cat(str(CHI_Sq_DEF))
  cat("\n")
  
  
  
  CHI_Sq_DEF$COMPARISON_category<-factor(CHI_Sq_DEF$COMPARISON_category,
                                         levels=c("LOF","MISS","SYN","UTR5","UTR3",
                                                  "INTRON","INTERGENIC","UPSTREAM","DOWNSTREAM","REGULATORY",
                                                  "TFBS","SPLICE",
                                                  "OTHER","NMD","NCT",
                                                  "NCGR","ASE_CTRL","enhancer_CTRL","Kousik_variant"),
                                         ordered=T)
  
  CHI_Sq_DEF$Category<-factor(CHI_Sq_DEF$Category,
                              levels = c("0","AT LEAST 1","AT LEAST 2","AT LEAST 3"),
                              ordered=T)
  
  CHI_Sq_DEF$log_pval<-as.numeric(CHI_Sq_DEF$log_pval)
  
  CHI_Sq_DEF$Cell_Type<-factor(CHI_Sq_DEF$Cell_Type,
                               levels = c("K562","CHRF","HL60","THP1"),
                               ordered=T)
  
  
  cat("CHI_Sq_DEF_1\n")
  cat(str(CHI_Sq_DEF))
  cat("\n")
  
  #### SAVE----
  
  
  setwd(path5)
  
  saveRDS(list_DEF,file="worstCSQ_E_Plus_ASE_CLASS_graphs.rds")
  
  saveRDS(CHI_Sq_DEF,file="worstCSQ_E_Plus_ASE_CLASS_stats.rds")
  
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
  
  ##### E_Plus_ASE stacked barplot ----
  
  
  path5<-paste(out,'Tier_plots_ER','/', sep='')
  
  cat("path5\n")
  cat(sprintf(as.character(path5)))
  cat("\n")
  
  
  if (file.exists(path5)){
    
    
    
    
  } else {
    dir.create(file.path(path5))
    
  }
  
  
  setwd(path5)
  
  
 #### READ DATA ----
 
 enhancer_CLASS_graphs<-readRDS(file="worstCSQ_enhancer_CLASS_graphs.rds")  

  # cat("enhancer_CLASS_graphs\n")
  # cat(str(enhancer_CLASS_graphs))
  # cat("\n")
  
 enhancer_CLASS_stats<-readRDS(file="worstCSQ_enhancer_CLASS_stats.rds")
 enhancer_CLASS_stats$ER_Type<-"enhancer_CLASS"
 
 cat("enhancer_CLASS_stats\n")
 cat(str(enhancer_CLASS_stats))
 cat("\n")
 
 
 E_Plus_ASE_CLASS_graphs<-readRDS(file="worstCSQ_E_Plus_ASE_CLASS_graphs.rds")  
 
 E_Plus_ASE_CLASS_stats<-readRDS(file="worstCSQ_E_Plus_ASE_CLASS_stats.rds")
 E_Plus_ASE_CLASS_stats$ER_Type<-"E_Plus_ASE_CLASS"
 
 cat("E_Plus_ASE_CLASS_stats\n")
 cat(str(E_Plus_ASE_CLASS_stats))
 cat("\n")
 
 Cell_Type_array<-names(enhancer_CLASS_graphs)
 
 cat("Cell_Type_array\n")
 cat(str(Cell_Type_array))
 cat("\n")

 
 DEF_stats<-rbind(enhancer_CLASS_stats,
                  E_Plus_ASE_CLASS_stats)
 
 cat("DEF_stats\n")
 cat(str(DEF_stats))
 cat("\n")
 
 for(i in 1:length(Cell_Type_array))
 {
   
   Cell_Type_array_sel<-Cell_Type_array[i]
   
   cat("---------------------------------------------->\t")
   cat(sprintf(as.character(Cell_Type_array_sel)))
   cat("\n")
   
   graph_enhancer<-enhancer_CLASS_graphs[[Cell_Type_array_sel]]
   
   graph_enhancer<-graph_enhancer+
                    theme(axis.title.y=element_text(size=18, family="sans"),
                    axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
                    axis.text.x=element_text(angle=45,vjust=1,hjust=1, size=16, color="black", family="sans"))+
                    theme(legend.position="hidden")
   
   # cat("graph_enhancer\n")
   # cat(str(graph_enhancer))
   # cat("\n")
   
   graph_E_Plus_ASE<-E_Plus_ASE_CLASS_graphs[[Cell_Type_array_sel]]
   
   
   graph_E_Plus_ASE<-graph_E_Plus_ASE+
     theme(axis.title.y=element_text(size=18, family="sans"),
           axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
           axis.text.x=element_text(angle=45,vjust=1,hjust=1, size=16, color="black", family="sans"))+
     theme(legend.position="hidden")
    
    
   
   graph_DEF2<-plot_grid(graph_enhancer,graph_E_Plus_ASE,
                          nrow = 1,
                          ncol=2,
                          labels = c("ACTIVE enhancer TILES)","ACTIVE E+ASE TILES)"),
                          label_size = 12,
                          align = "h",
                          rel_widths=c(1,1))
   
   title <- ggdraw() + draw_label(paste("ER","VEP_wCSQ+CTRLs",Cell_Type_array_sel,sep="  "))


   graph_DEF3<- plot_grid(title,graph_DEF2,
             ncol=1,
             rel_heights = c(0.05,1))
   
   
   svgname<-paste("ER","_VEP_wCSQ_",Cell_Type_array_sel,".svg", sep='')
   makesvg = TRUE
   
   if (makesvg == TRUE)
   {
     
     ggsave(svgname, plot= graph_DEF3,
            device="svg",
            height=10, width=12)
   }
   
   
   graph_enhancer<-graph_enhancer+
                  theme(legend.position="right")
   
      
   graph_E_Plus_ASE<-graph_E_Plus_ASE+theme(legend.position="right")
   
   
   
   graph_DEF2<-plot_grid(graph_enhancer,graph_E_Plus_ASE,
                         nrow = 1,
                         ncol=2,
                         labels = c("ACTIVE enhancer TILES)","ACTIVE E+ASE TILES)"),
                         label_size = 12,
                         align = "h",
                         rel_widths=c(1,1))
   
   title <- ggdraw() + draw_label(paste("ER","VEP_wCSQ+CTRLs",Cell_Type_array_sel,sep="  "))
   
   
   graph_DEF3<- plot_grid(title,graph_DEF2,
                          ncol=1,
                          rel_heights = c(0.05,1))
   
   
   svgname<-paste("ER","_VEP_wCSQ_",Cell_Type_array_sel,"_WITH_LEGENDS.svg", sep='')
   makesvg = TRUE
   
   if (makesvg == TRUE)
   {
     
     ggsave(svgname, plot= graph_DEF3,
            device="svg",
            height=10, width=12)
   }
   
   
   
 }#i
 
 setwd(path5)
 
 write.table(DEF_stats,file="stats_object_ER_Chi.tsv",sep="\t",quote=F,row.names = F)
  
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
  
  
  
  Tier_stacked_barplots(opt)
  VEP_LABEL_ER_stacked_barplots_E_Plus_ASE_CLASS(opt)
  VEP_LABEL_ER_stacked_barplots_enhancer_CLASS(opt)
  print_join(opt)
  
}
  
  
  
 

###########################################################################

system.time( main() )
