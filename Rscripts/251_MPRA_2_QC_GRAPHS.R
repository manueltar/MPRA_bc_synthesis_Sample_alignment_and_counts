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
library("reshape2", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("gtools", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")

library("extrafont", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("extrafontdb", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("sysfonts", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("showtextdb", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("showtext", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("jsonlite", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("curl", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("ggridges", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")



# quit(status = 1)

opt = NULL

options(warn = 1)


ggridge.graphs = function(option_list)
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
  
  
  #### READ INPUT FILES ----
  
  setwd(out)
  
  
  
  filename_1<-paste("list_values",".rds", sep='')
  
  list_values<-readRDS(file=filename_1)
  
  features<-names(list_values)
  
  cat("features\n")
  cat(sprintf(as.character(features)))
  cat("\n")
  
  
  
  filename_1<-paste("list_medians",".rds", sep='')
  
  list_medians<-readRDS(file=filename_1)
  
  # cat("list_medians\n")
  # cat(str(list_medians))
  # cat("\n")
  
  filename_1<-paste("Rosetta_df",".rds", sep='')
  
  
  Rosetta_df<-readRDS(file=filename_1)
  
  cat("Rosetta_df\n")
  cat(str(Rosetta_df))
  cat("\n")
  
  
  #### path2 ----
  
  path2<-paste(out,'QC_graphs','/', sep='')
  
  cat("path2\n")
  cat(sprintf(as.character(path2)))
  cat("\n")
  
  if (file.exists(path2)){
    
    
    
    
  } else {
    dir.create(file.path(path2))
    
  }
  
  
  
  ##### LOOP ggridge ----
  
  
  for(i in 1:length(features))
  {
    features_sel<-features[i]
    
    cat("----------------------------------------------------------------------------------------->\t")
    cat(sprintf(as.character(features_sel)))
    cat("\n")
    
    list_values_sel<-unique(list_values[[features_sel]])
    
    

    cat("list_values_sel_0\n")
    cat(str(list_values_sel))
    cat("\n")
    
    list_values_sel<-droplevels(unique(list_values_sel[,-c(which(colnames(list_values_sel) == "type"))]))
    
    
    cat("list_values_sel_1\n")
    cat(str(list_values_sel))
    cat("\n")
    
    ##### Elements per master sample
    
    list_values_sel.dt<-data.table(list_values_sel, key="master_sample")
    
    
    
    Freq.table<-as.data.frame(list_values_sel.dt[,.N,by=key(list_values_sel.dt)], stringsAsFactors=F)
    
    colnames(Freq.table)[which(colnames(Freq.table) == "N")]<-"instances"
    
    cat("Freq.table_\n")
    cat(str(Freq.table))
    cat("\n")
    
    
    list_values_sel<-merge(list_values_sel,Freq.table,
                           by="master_sample",
                           all.x=T)
    
    cat("list_values_sel_2\n")
    cat(str(list_values_sel))
    cat("\n")
    
    #### Cell_Type colors
    
    Cell_Type_levels<-levels(list_values_sel$Cell_Type)
    colors_Cell_Type_levels<-c('#32A852','#1877C9','#553B68','#D45E85','#6DB2EE','#62D07F','#C9244B','#87447B','#D6B8E6')
    
    df.color_Cell_Type<-as.data.frame(cbind(Cell_Type_levels,colors_Cell_Type_levels[1:length(Cell_Type_levels)]), stringsAsFactors=F)
    colnames(df.color_Cell_Type)<-c("Cell_Type","colors")
    
    df.color_Cell_Type$Cell_Type<-factor(df.color_Cell_Type$Cell_Type,
                                         levels=levels(list_values_sel$Cell_Type),
                                         ordered=T)
    
    cat("df.color_Cell_Type_\n")
    cat(str(df.color_Cell_Type))
    cat("\n")
    
    
   
    # 
    # quit(status = 1)
    # 
    list_values_sel<-merge(list_values_sel,
                           df.color_Cell_Type,
                           by="Cell_Type")
    
    
    cat("list_values_sel_3\n")
    cat(str(list_values_sel))
    cat("\n")
    
    
    special_features<-c("LogFC","Vockley_REF","dna.REF.vs.dna.ALT","rna.REF.vs.rna.ALT")
    
    indx.feature<-"NA"
    
    
    FLAG_feature<-length(features_sel[features_sel%in%special_features])
    
    cat("FLAG_feature\n")
    cat(sprintf(as.character(FLAG_feature)))
    cat("\n")
    
      
      
    if(FLAG_feature >0)
    {
      
      indx.feature<-which(colnames(list_values_sel) == features_sel)
      
    }
    if(FLAG_feature  == 0)
    {
      
      indx.feature<-which(colnames(list_values_sel) == "Median_value")
      
    }
    
    cat("indx.feature\n")
    cat(sprintf(as.character(indx.feature)))
    cat("\n")
    
    
    list_values_sel_colors<-unique(list_values_sel[,c(which(colnames(list_values_sel) == "master_sample"),
                                                      which(colnames(list_values_sel) == "colors"))])
    
    cat("list_values_sel_colors\n")
    cat(str(list_values_sel_colors))
    cat("\n")
    
   
    
    if(features_sel == "LogFC")
    {
      
      
      A<-summary(list_values_sel[,indx.feature])
      
      
      
      cat("A\n")
      cat(sprintf(as.character(names(A))))
      cat("\n")
      cat(sprintf(as.character(A)))
      cat("\n")
      
      max_log_value<-A[6]
      min_log_value<-A[1]
      
      breaks.log_value<-c(min_log_value,seq(-5,5, by=0.5),max_log_value)
      labels_value<-as.character(round(logratio2foldchange(breaks.log_value, base=2),2))
      
      
      cat("labels_value\n")
      cat(sprintf(as.character(labels_value)))
      cat("\n")
      
      
     
    
      
      ggridge_plot<-list_values_sel[order(list_values_sel$master_sample),] %>%
        mutate(myaxis = paste0(list_values_sel$master_sample, " ", "n=", list_values_sel$instances)) %>%
        mutate(myaxis=fct_reorder(myaxis,as.numeric(list_values_sel$master_sample))) %>%
        ggplot(aes(y=myaxis, x=list_values_sel[,indx.feature], fill=Cell_Type,alpha=Cell_Type)) +
        stat_density_ridges(quantile_lines = TRUE, quantiles = 2, geom="density_ridges",scale=5)+
        theme_bw()+
        theme(axis.title.y=element_text(size=24, family="sans"),
              axis.title.x=element_text(size=24, family="sans"),
              axis.text.y=element_text(angle=0,size=10, color="black", family="sans"),
              axis.text.x=element_text(angle=45, vjust=1,hjust=1,size=18, color="black", family="sans"),
              legend.title=element_text(size=16,color="black", family="sans"),
              legend.text=element_text(size=12,color="black", family="sans"))+
        scale_y_discrete(name=NULL, drop=F)+
        scale_x_continuous(name="LogFC",breaks=breaks.log_value,labels=labels_value, limits=c(breaks.log_value[1],breaks.log_value[length(breaks.log_value)]))+
        scale_fill_manual(values = df.color_Cell_Type$colors,
                          drop=F)+
        scale_alpha_manual(values=c(2, 0.5, 2,0.5)) +
        theme(legend.position = "none")+
        ggeasy::easy_center_title()
      
      setwd(path2)
      
      svgname<-paste("ggridge_",features_sel,".svg",sep='')
      makesvg = TRUE
      
      if (makesvg == TRUE)
      {
        ggsave(svgname, plot= ggridge_plot,
               device="svg",
               height=10, width=12)
      }
      
      
      
      list_values_sel_subset<-list_values_sel[which(list_values_sel[,indx.feature] > 0),]
      
      
      
      cat("----------------------------->list_values_sel_subset_0\n")
      cat(str(list_values_sel_subset))
      cat("\n")
      
      list_values_sel_subset<-droplevels(list_values_sel_subset)
      
      cat("list_values_sel_subset_1\n")
      cat(str(list_values_sel_subset))
      cat("\n")
      
      
      list_values_sel_subset.dt<-data.table(list_values_sel_subset, key="master_sample")
      
      
      
      Freq.table<-as.data.frame(list_values_sel_subset.dt[,.N,by=key(list_values_sel_subset.dt)], stringsAsFactors=F)
      
      colnames(Freq.table)[which(colnames(Freq.table) == "N")]<-"instances_2"
      
      cat("Freq.table_\n")
      cat(str(Freq.table))
      cat("\n")
      
      
      list_values_sel_subset<-merge(list_values_sel_subset,Freq.table,
                             by="master_sample",
                             all.x=T)
      
      cat("list_values_sel_subset_2\n")
      cat(str(list_values_sel_subset))
      cat("\n")
      
      
      A<-summary(list_values_sel_subset[,indx.feature])
      
      
      
      cat("A\n")
      cat(sprintf(as.character(names(A))))
      cat("\n")
      cat(sprintf(as.character(A)))
      cat("\n")
      
      max_log_value<-A[6]
      min_log_value<-A[1]
      
      breaks.log_value<-c(min_log_value,seq(0,2, by=0.25),5)
      labels_value<-as.character(round(logratio2foldchange(breaks.log_value, base=2),2))
      
      
      cat("labels_value\n")
      cat(sprintf(as.character(labels_value)))
      cat("\n")
      
      
      
    
      
      ggridge_plot<-list_values_sel_subset[order(list_values_sel_subset$master_sample),] %>%
        mutate(myaxis = paste0(list_values_sel_subset$master_sample, " ", "n=", list_values_sel_subset$instances_2)) %>%
        mutate(myaxis=fct_reorder(myaxis,as.numeric(list_values_sel_subset$master_sample))) %>%
        ggplot(aes(y=myaxis, x=list_values_sel_subset[,indx.feature], fill=Cell_Type,alpha=Cell_Type)) +
        stat_density_ridges(quantile_lines = TRUE, quantiles = 2, geom="density_ridges",scale=5)+
        theme_bw()+
        theme(axis.title.y=element_text(size=24, family="sans"),
              axis.title.x=element_text(size=24, family="sans"),
              axis.text.y=element_text(angle=0,size=10, color="black", family="sans"),
              axis.text.x=element_text(angle=45, vjust=1,hjust=1,size=18, color="black", family="sans"),
              legend.title=element_text(size=16,color="black", family="sans"),
              legend.text=element_text(size=12,color="black", family="sans"))+
        scale_y_discrete(name=NULL, drop=F)+
        scale_x_continuous(name="LogFC",breaks=breaks.log_value,labels=labels_value, limits=c(breaks.log_value[1],breaks.log_value[length(breaks.log_value)]))+
        scale_fill_manual(values = df.color_Cell_Type$colors,
                          drop=F)+
        scale_alpha_manual(values=c(2, 1, 2,1)) +
        theme(legend.position = "none")+
        ggeasy::easy_center_title()
      
      setwd(path2)
      
      svgname<-paste("ggridge_SUBSET_",features_sel,".svg",sep='')
      makesvg = TRUE
      
      if (makesvg == TRUE)
      {
        ggsave(svgname, plot= ggridge_plot,
               device="svg",
               height=10, width=12)
      }
      
      
    
      
   
   
    }else{
      
      A<-summary(log10(list_values_sel[,indx.feature]+ 0.00001))
      
      
      
      cat("A\n")
      cat(sprintf(as.character(names(A))))
      cat("\n")
      cat(sprintf(as.character(A)))
      cat("\n")
      
      max_log_value<-A[6]
      min_log_value<-A[1]
      
      breaks.log_value<-seq(min_log_value,max_log_value+0.5, by=0.5)
      labels_value<-10^breaks.log_value
      
      labels_value[which(labels_value < 1)]<-round(labels_value[which(labels_value < 1)],5)
      labels_value[which(labels_value >= 1)]<-round(labels_value[which(labels_value >= 1)],0)
      
      labels_value<-as.character(labels_value)
      
      cat("labels_value\n")
      cat(sprintf(as.character(labels_value)))
      cat("\n")
      
      
      ggridge_plot<-list_values_sel[order(list_values_sel$master_sample),] %>%
        mutate(myaxis = paste0(list_values_sel$master_sample, " ", "n=", list_values_sel$instances)) %>%
        mutate(myaxis=fct_reorder(myaxis,as.numeric(list_values_sel$master_sample))) %>%
        ggplot(aes(y=myaxis, x=log10(list_values_sel[,indx.feature]+ 0.00001), fill=Cell_Type,alpha=Cell_Type)) +
        stat_density_ridges(quantile_lines = TRUE, quantiles = 2, geom="density_ridges",scale=5)+
        theme_bw()+
        theme(axis.title.y=element_text(size=24, family="sans"),
              axis.title.x=element_text(size=24, family="sans"),
              axis.text.y=element_text(angle=0,size=10, color="black", family="sans"),
              axis.text.x=element_text(angle=45, vjust=1,hjust=1,size=18, color="black", family="sans"),
              legend.title=element_text(size=16,color="black", family="sans"),
              legend.text=element_text(size=12,color="black", family="sans"))+
        scale_y_discrete(name=NULL, drop=F)+
        scale_x_continuous(name=paste(features_sel,"log10(value +0.00001)",sep = ' '),breaks=breaks.log_value,labels=labels_value, limits=c(breaks.log_value[1],breaks.log_value[length(breaks.log_value)]))+
        scale_fill_manual(values = df.color_Cell_Type$colors,
                          drop=F)+
        scale_alpha_manual(values=c(2, 0.5, 2,0.5)) +
        theme(legend.position = "none")+
        ggeasy::easy_center_title()
      
      setwd(path2)
      
      svgname<-paste("ggridge_",features_sel,".svg",sep='')
      makesvg = TRUE
      
      if (makesvg == TRUE)
      {
        ggsave(svgname, plot= ggridge_plot,
               device="svg",
               height=10, width=12)
      }
      
      if(features_sel == "Vockley_REF")
      {
        
        list_values_sel_subset_1<-list_values_sel[which(list_values_sel[,indx.feature] >= 1.1 & list_values_sel[,indx.feature] < 10),]
        list_values_sel_subset_2<-list_values_sel[which(list_values_sel[,indx.feature] <= 0.8 & list_values_sel[,indx.feature] > 0.1),]
        
        list_values_sel_subset<-rbind(list_values_sel_subset_1,list_values_sel_subset_2)
        
        
        
        cat("----------------------------->list_values_sel_subset_0\n")
        cat(str(list_values_sel_subset))
        cat("\n")
        
        list_values_sel_subset<-droplevels(list_values_sel_subset)
        
        cat("list_values_sel_subset_1\n")
        cat(str(list_values_sel_subset))
        cat("\n")
        
        
        list_values_sel_subset.dt<-data.table(list_values_sel_subset, key="master_sample")
        
        
        
        Freq.table<-as.data.frame(list_values_sel_subset.dt[,.N,by=key(list_values_sel_subset.dt)], stringsAsFactors=F)
        
        colnames(Freq.table)[which(colnames(Freq.table) == "N")]<-"instances_2"
        
        cat("Freq.table_\n")
        cat(str(Freq.table))
        cat("\n")
        
        
        list_values_sel_subset<-merge(list_values_sel_subset,Freq.table,
                                      by="master_sample",
                                      all.x=T)
        
        cat("list_values_sel_subset_2\n")
        cat(str(list_values_sel_subset))
        cat("\n")
        
        
        A<-summary(log10(list_values_sel_subset[,indx.feature]+ 0.00001))
        
        
        
        cat("A\n")
        cat(sprintf(as.character(names(A))))
        cat("\n")
        cat(sprintf(as.character(A)))
        cat("\n")
        
        max_log_value<-A[6]
        min_log_value<-A[1]
        
        breaks.log_value<-seq(min_log_value,max_log_value+0.25, by=0.25)
        labels_value<-10^breaks.log_value
        
        labels_value[which(labels_value < 1)]<-round(labels_value[which(labels_value < 1)],5)
        labels_value[which(labels_value >= 1)]<-round(labels_value[which(labels_value >= 1)],0)
        
        labels_value<-as.character(labels_value)
        
        cat("labels_value\n")
        cat(sprintf(as.character(labels_value)))
        cat("\n")
        
        
        ggridge_plot<-list_values_sel_subset[order(list_values_sel_subset$master_sample),] %>%
          mutate(myaxis = paste0(list_values_sel_subset$master_sample, " ", "n=", list_values_sel_subset$instances_2)) %>%
          mutate(myaxis=fct_reorder(myaxis,as.numeric(list_values_sel_subset$master_sample))) %>%
          ggplot(aes(y=myaxis, x=log10(list_values_sel_subset[,indx.feature]+ 0.00001), fill=Cell_Type,alpha=Cell_Type)) +
          stat_density_ridges(quantile_lines = TRUE, quantiles = 2, geom="density_ridges",scale=5)+
          theme_bw()+
          theme(axis.title.y=element_text(size=24, family="sans"),
                axis.title.x=element_text(size=24, family="sans"),
                axis.text.y=element_text(angle=0,size=10, color="black", family="sans"),
                axis.text.x=element_text(angle=45, vjust=1,hjust=1,size=18, color="black", family="sans"),
                legend.title=element_text(size=16,color="black", family="sans"),
                legend.text=element_text(size=12,color="black", family="sans"))+
          scale_y_discrete(name=NULL, drop=F)+
          scale_x_continuous(name=paste(features_sel,"log10(value +0.00001)",sep = ' '),breaks=breaks.log_value,labels=labels_value, limits=c(breaks.log_value[1],breaks.log_value[length(breaks.log_value)]))+
          scale_fill_manual(values = df.color_Cell_Type$colors,
                            drop=F)+
          scale_alpha_manual(values=c(2, 0.5, 2,0.5)) +
          theme(legend.position = "none")+
          ggeasy::easy_center_title()
        
        setwd(path2)
        
        svgname<-paste("ggridge_SUBSET_",features_sel,".svg",sep='')
        makesvg = TRUE
        
        if (makesvg == TRUE)
        {
          ggsave(svgname, plot= ggridge_plot,
                 device="svg",
                 height=10, width=12)
        }
      }
      
    }
    
  }#i features
}


corr.graphs = function(option_list)
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
  
  
  #### READ INPUT FILES ----
  
  setwd(out)
  
  
  
  filename_1<-paste("list_values",".rds", sep='')
  
  list_values<-readRDS(file=filename_1)
  
  features<-names(list_values)
  
  cat("features\n")
  cat(sprintf(as.character(features)))
  cat("\n")
  
  
  
  filename_1<-paste("list_medians",".rds", sep='')
  
  list_medians<-readRDS(file=filename_1)
  
  # cat("list_medians\n")
  # cat(str(list_medians))
  # cat("\n")
  
  filename_1<-paste("Rosetta_df",".rds", sep='')
  
  
  Rosetta_df<-readRDS(file=filename_1)
  
  cat("Rosetta_df\n")
  cat(str(Rosetta_df))
  cat("\n")
  
  
  #### path2 ----
  
  path2<-paste(out,'QC_graphs','/', sep='')
  
  cat("path2\n")
  cat(sprintf(as.character(path2)))
  cat("\n")
  
  if (file.exists(path2)){
    
    
    
    
  } else {
    dir.create(file.path(path2))
    
  }
  

  
  ##### LOOP ggridge ----
  
  
  for(i in 1:length(features))
  {
    features_sel<-features[i]
    
    cat("----------------------------------------------------------------------------------------->\t")
    cat(sprintf(as.character(features_sel)))
    cat("\n")
    
    list_values_sel<-unique(list_values[[features_sel]])
    
    
    
    
    cat("list_values_sel_0\n")
    cat(str(list_values_sel))
    cat("\n")
    
    
    list_values_sel<-droplevels(unique(list_values_sel[,-c(which(colnames(list_values_sel) == "type"))]))
    
    cat("list_values_sel_1\n")
    cat(str(list_values_sel))
    cat("\n")
    
    special_features<-c("LogFC","Vockley_REF","dna.REF.vs.dna.ALT","rna.REF.vs.rna.ALT")
    
    indx.feature<-"NA"
    
    
    FLAG_feature<-length(features_sel[features_sel%in%special_features])
    
    cat("FLAG_feature\n")
    cat(sprintf(as.character(FLAG_feature)))
    cat("\n")
    
    
    
    if(FLAG_feature >0)
    {
      
      indx.feature<-which(colnames(list_values_sel) == features_sel)
      
    }
    if(FLAG_feature  == 0)
    {
      
      indx.feature<-which(colnames(list_values_sel) == "Median_value")
      
    }
    
    cat("indx.feature\n")
    cat(sprintf(as.character(indx.feature)))
    cat("\n")
    
    list_values_sel_NO_NA<-list_values_sel[!is.na(list_values_sel[,indx.feature]),]
    
    
    cat("list_values_sel_NO_NA\n")
    cat(str(list_values_sel_NO_NA))
    cat("\n")
    
    
    
    list_values_sel_wide<-as.data.frame(pivot_wider(list_values_sel_NO_NA,
                                       id_cols=REAL_TILE_Plus_carried_variants,
                                       names_from=master_sample,
                                       values_from=colnames(list_values_sel)[indx.feature]),stringsAsFactors=F)
    
    cat("list_values_sel_wide_0\n")
    cat(str(list_values_sel_wide))
    cat("\n")
    
    
    list_values_sel_wide_matrix<-as.matrix(list_values_sel_wide[,-1])
    row.names(list_values_sel_wide_matrix)<-list_values_sel_wide$REAL_TILE_Plus_carried_variants
    
    cat("list_values_sel_wide_matrix_2\n")
    cat(str(list_values_sel_wide_matrix))
    cat("\n")
    
    
    # list_values_sel_wide_matrix[is.na(list_values_sel_wide_matrix)]<-0
    # 
    # cat("list_values_sel_wide_matrix_3\n")
    # cat(str(list_values_sel_wide_matrix))
    # cat("\n")
    
    
    setwd(path2)
    write.table(list_values_sel_wide_matrix,file="test.tsv", sep="\t")
    # quit(status=1)
    
    Corr.Matrix<-cor(list_values_sel_wide_matrix, method = c("pearson"))
    
    cat("Corr.Matrix\n")
    str(Corr.Matrix)
    cat("\n")
    
    Corr.Matrix[is.na(Corr.Matrix)]<-0
    
    cat("Corr.Matrix_POST\n")
    str(Corr.Matrix)
    cat("\n")
    
    setwd(path2)
    write.table(Corr.Matrix,file="test_2.tsv", sep="\t")
    
    # Reorder the correlation matrix according to clusters
    
    cormat <- reorder_cormat(Corr.Matrix)
    cormat<-Corr.Matrix
    
    
    check<-colnames(cormat)
    
    cat("check_reorder\n")
    cat(sprintf(as.character(check)))
    cat("\n")
    
    
    
    # quit(status=1)
    
    # cormat<-Corr.Matrix
    lower_tri <- get_lower_tri(cormat)
    
    cat("lower_tri\n")
    str(lower_tri)
    cat("\n")
    
    new_order<-row.names(lower_tri)
    
    cat("new_order\n")
    cat(sprintf(as.character(new_order)))
    cat("\n")
    
    ##### h2 and h3
    
    h2_h3_df<-unique(list_values_sel[,c(which(colnames(list_values_sel) == "master_sample"),
                                 which(colnames(list_values_sel) == "Cell_Type"),
                                 which(colnames(list_values_sel) == "Replicate"))])
    
    cat("h2_h3_df_0\n")
    str(h2_h3_df)
    cat("\n")
    
    
    h2_h3_df$master_sample<-factor(h2_h3_df$master_sample,
                                             levels=new_order,
                                             ordered = T)
    
    h2_h3_df<-h2_h3_df[order(h2_h3_df$Cell_Type,h2_h3_df$master_sample),]
    
   
    h2_h3_df$MOCK_COORD<-"MOCK"
    
    cat("h2_h3_df\n")
    str(h2_h3_df)
    cat("\n")
    cat(sprintf(as.character(h2_h3_df$Cell_Type)))
    cat("\n")
    cat(sprintf(as.character(levels(h2_h3_df$Cell_Type))))
    cat("\n")
    cat(sprintf(as.character(h2_h3_df$master_sample)))
    cat("\n")
    cat(sprintf(as.character(levels(h2_h3_df$master_sample))))
    cat("\n")
    
    
    #### Cell_Type colors
    
    Cell_Type_levels<-levels(h2_h3_df$Cell_Type)
    colors_Cell_Type_levels<-c('#32A852','#1877C9','#553B68','#D45E85','#6DB2EE','#62D07F','#C9244B','#87447B','#D6B8E6')
    
    df.color_Cell_Type<-as.data.frame(cbind(Cell_Type_levels,colors_Cell_Type_levels[1:length(Cell_Type_levels)]), stringsAsFactors=F)
    colnames(df.color_Cell_Type)<-c("Cell_Type","colors")
    
    df.color_Cell_Type$Cell_Type<-factor(df.color_Cell_Type$Cell_Type,
                                         levels=levels(h2_h3_df$Cell_Type),
                                         ordered=T)
    
    cat("df.color_Cell_Type_\n")
    cat(str(df.color_Cell_Type))
    cat("\n")
    
    
   
    
    h2 <- ggplot(data=h2_h3_df, aes(x=MOCK_COORD,
                                                     y=master_sample,
                                                     fill = Cell_Type))+
      geom_tile(color = "white")+
      scale_fill_manual(values = df.color_Cell_Type$colors,
                        drop=F)+
      ggeasy::easy_center_title()+
      theme_minimal()+ # minimal theme
      scale_x_discrete(name=element_blank(), drop =F) +
      scale_y_discrete(name=NULL, drop=F)+
      theme(legend.position = "hidden", legend.title = element_blank())+
      theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
      theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    
    
    h3 <- ggplot(data=h2_h3_df, aes(x=master_sample,
                                                     y=MOCK_COORD,
                                                     fill = Cell_Type))+
      geom_tile(color = "white")+
      scale_fill_manual(values = df.color_Cell_Type$colors,
                        drop=F)+
      ggeasy::easy_center_title()+
      theme_minimal()+ # minimal theme
      scale_x_discrete(name=element_blank(), drop =F) +
      scale_y_discrete(name=NULL, drop=F)+
      theme(legend.position = "hidden", legend.title = element_blank())+
      theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
      theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    
    
    
  
    
    #################### Melt the correlation matrix
    
    cat("lower_tri\n")
    str(lower_tri)
    cat("\n")
    
    melted_cormat <- melt(lower_tri, na.rm = TRUE)
    
    cat("melted_cormat\n")
    str(melted_cormat)
    cat("\n")
    
    colnames(melted_cormat)[which(colnames(melted_cormat) == "Var2")]<-"master_sample"
    
    
    melted_cormat<-merge(melted_cormat,
                         h2_h3_df,
                         by="master_sample")
    
    
    
    cat("melted_cormat\n")
    str(melted_cormat)
    cat("\n")
    
    
    levels_Cell_Type<-levels(melted_cormat$Cell_Type)
    
    vector_COORDS<-0
    
    FINAL_COORDS<-NULL
    
    for(i in 1:length(levels_Cell_Type))
    {
      
      Cell_Type_sel<-levels_Cell_Type[i]
      
      cat("--->\t")
      cat(sprintf(as.character(Cell_Type_sel)))
      cat("\n")
      
      melted_cormat_sel<-melted_cormat[which(melted_cormat$Cell_Type == Cell_Type_sel),]
      
      cat("melted_cormat_sel\n")
      str(melted_cormat_sel)
      cat("\n")
      
      melted_cormat_sel_master_samples<-unique(as.character(melted_cormat_sel$master_sample))
      
      cat("melted_cormat_sel_master_samples\n")
      str(melted_cormat_sel_master_samples)
      cat("\n")
      
      COORD_X<-length(melted_cormat_sel_master_samples)
      
      cat("COORD_X--->\t")
      cat(sprintf(as.character(COORD_X)))
      cat("\n")
      
      COORD_DEF<-sum(vector_COORDS,COORD_X)
      
      cat("COORD_DEF--->\t")
      cat(sprintf(as.character(COORD_DEF)))
      cat("\n")
      
      vector_COORDS[i]<-COORD_X
      
      
      # COORD_Y<-length(melted_cormat_sel_master_samples)
      # 
      # cat("COORD_Y--->\t")
      # cat(sprintf(as.character(COORD_Y)))
      # cat("\n")
      
      FINAL_COORDS[i]<-as.numeric(COORD_DEF)
      
    }
    
    FINAL_COORDS<-FINAL_COORDS+0.5
    
    FINAL_COORDS<-FINAL_COORDS[-length(FINAL_COORDS)]
    
    cat("FINAL_COORDS--->\t")
    cat(sprintf(as.character(FINAL_COORDS)))
    cat("\n")
    
    #### heatmap graph
    
    breaks.Rank<-seq(-1,1,by=0.5)
    labels.Rank<-as.character(breaks.Rank)
    
    
    ggheatmap <- ggplot(melted_cormat, aes(master_sample, Var1, fill = value))+
      geom_tile(color = "white")+
      scale_fill_gradient2(name = paste("Pearson","Correlation","Coefficient", sep="\n"),
                           low = "red", high = "green",mid="white",midpoint=0,
                           na.value = NA,
                           breaks=breaks.Rank,
                           labels=labels.Rank,
                           limits=c(breaks.Rank[1],
                                    breaks.Rank[length(breaks.Rank)]))+
      theme_minimal()+ # minimal theme
      scale_x_discrete(name=element_blank(), drop =F) +
      scale_y_discrete(name=NULL, drop=F)+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                       size = 12, hjust = 1))+
      geom_hline(yintercept = FINAL_COORDS, linetype=3, color="black")+
      geom_vline(xintercept = FINAL_COORDS, linetype=3, color="black")+
      coord_fixed()
    
    
    
    ggheatmap2<-ggheatmap+
      theme(axis.text.y  = element_blank())+
      theme(axis.text.x  = element_blank())+
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.8, 0.2),
        legend.direction = "horizontal")+
      guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                   title.position = "top", title.hjust = 0.5))+
      theme(text=element_text(size=24,  family="sans"))
    
    
    h3<-plot_grid(NULL,h3,NULL,
                  nrow = 1,
                  ncol = 3,
                  align = "hv",
                  rel_widths = c(0.18, 1,0.1),
                  rel_heights = 1)
    
    graph_DEF<-plot_grid(h2,ggheatmap2,
                         nrow = 1,
                         ncol = 2,
                         align = "hv",
                         rel_widths = c(0.065, 1),
                         rel_heights = 1)
    
    graph_DEF<-plot_grid(h3,graph_DEF,
                         nrow = 2,
                         ncol = 1,
                         align = "hv",
                         rel_widths = 1,
                         rel_heights = c(0.065,1))
    
    
    
    title <- ggdraw() + draw_label(paste("Corr plot",features_sel,sep="  "))
    
    
  
    
    cat("svg_graph\n")
    
    setwd(path2)
    
    
    
    svgname<-paste("Graph_corr_",features_sel,".svg", sep='')
    makesvg = TRUE
    
    if (makesvg == TRUE)
    {
      
      ggsave(svgname, plot= graph_DEF,
             device="svg",
             height=10, width=12)
    }
    
    
   
    graph_DEF2<-plot_grid(h2,ggheatmap,
                          nrow = 2,
                          ncol = 2,
                          align = "hv",
                          rel_widths = c(0.065, 1),
                          rel_heights = c(1,0.05))
    
    graph_DEF3<-plot_grid(title,graph_DEF2,
                          ncol=1,
                          rel_heights = c(0.1,1))
    
    
    svgname<-paste("Graph_corr_NAMED_",features_sel,".svg", sep='')
    makesvg = TRUE
    
    if (makesvg == TRUE)
    {
      
      ggsave(svgname, plot= graph_DEF3,
             device="svg",
             height=10, width=12)
    }
    
    
    
  
    
    if(features_sel == "LogFC" | features_sel == "Vockley_REF")
    {
      if(features_sel == "LogFC")
      {
        list_values_sel_subset<-list_values_sel[which(list_values_sel[,indx.feature] > 0),]
        
      }
      
      if(features_sel == "Vockley_REF")
      {
        list_values_sel_subset_1<-list_values_sel[which(list_values_sel[,indx.feature] >= 1.1 & list_values_sel[,indx.feature] < 10),]
        list_values_sel_subset_2<-list_values_sel[which(list_values_sel[,indx.feature] <= 0.8 & list_values_sel[,indx.feature] > 0.1),]
        
        list_values_sel_subset<-rbind(list_values_sel_subset_1,list_values_sel_subset_2)
      }
      
      
      
      
      cat("----------------------------->list_values_sel_subset_0\n")
      cat(str(list_values_sel_subset))
      cat("\n")
      
      list_values_sel_subset<-droplevels(list_values_sel_subset)
      
      cat("list_values_sel_subset_1\n")
      cat(str(list_values_sel_subset))
      cat("\n")
      
      
        list_values_sel_subset_wide<-as.data.frame(pivot_wider(list_values_sel_subset,
                                          id_cols=REAL_TILE_Plus_carried_variants,
                                          names_from=master_sample,
                                          values_from=colnames(list_values_sel_subset)[indx.feature]),stringsAsFactors=F)
        
        cat("list_values_sel_subset_wide_0\n")
        cat(str(list_values_sel_subset_wide))
        cat("\n")
        
        
        list_values_sel_subset_wide_matrix<-as.matrix(list_values_sel_subset_wide[,-1])
        row.names(list_values_sel_subset_wide_matrix)<-list_values_sel_subset_wide$REAL_TILE_Plus_carried_variants
        
        cat("list_values_sel_subset_wide_matrix_2\n")
        cat(str(list_values_sel_subset_wide_matrix))
        cat("\n")
        
        
        
        list_values_sel_subset_wide_matrix[is.na(list_values_sel_subset_wide_matrix)]<-0
        
        cat("list_values_sel_subset_wide_matrix\n")
        cat(str(list_values_sel_subset_wide_matrix))
        cat("\n")
        
        Corr.Matrix<-cor(list_values_sel_subset_wide_matrix, method = c("pearson"))
        
        cat("Corr.Matrix\n")
        str(Corr.Matrix)
        cat("\n")
        
     
        
        # Reorder the correlation matrix according to clusters
        
        cormat <- reorder_cormat(Corr.Matrix)
        cormat<-Corr.Matrix
        
        
        check<-colnames(cormat)
        
        cat("check_reorder\n")
        cat(sprintf(as.character(check)))
        cat("\n")
        
        
        
        # quit(status=1)
        
        # cormat<-Corr.Matrix
        lower_tri <- get_lower_tri(cormat)
        
        cat("lower_tri\n")
        str(lower_tri)
        cat("\n")
        
        new_order<-row.names(lower_tri)
        
        cat("new_order\n")
        cat(sprintf(as.character(new_order)))
        cat("\n")
        
        ##### h2 and h3
        
        h2_h3_df<-unique(list_values_sel_subset[,c(which(colnames(list_values_sel_subset) == "master_sample"),
                                            which(colnames(list_values_sel_subset) == "Cell_Type"),
                                            which(colnames(list_values_sel_subset) == "Replicate"))])
        
        cat("h2_h3_df_0\n")
        str(h2_h3_df)
        cat("\n")
        
        
        h2_h3_df$master_sample<-factor(h2_h3_df$master_sample,
                                       levels=new_order,
                                       ordered = T)
        
        h2_h3_df<-h2_h3_df[order(h2_h3_df$Cell_Type,h2_h3_df$master_sample),]
        
        
        h2_h3_df$MOCK_COORD<-"MOCK"
        
        cat("h2_h3_df\n")
        str(h2_h3_df)
        cat("\n")
        cat(sprintf(as.character(h2_h3_df$Cell_Type)))
        cat("\n")
        cat(sprintf(as.character(levels(h2_h3_df$Cell_Type))))
        cat("\n")
        cat(sprintf(as.character(h2_h3_df$master_sample)))
        cat("\n")
        cat(sprintf(as.character(levels(h2_h3_df$master_sample))))
        cat("\n")
        
        #### Cell_Type colors
        
        Cell_Type_levels<-levels(h2_h3_df$Cell_Type)
        colors_Cell_Type_levels<-c('#32A852','#1877C9','#553B68','#D45E85','#6DB2EE','#62D07F','#C9244B','#87447B','#D6B8E6')
        
        df.color_Cell_Type<-as.data.frame(cbind(Cell_Type_levels,colors_Cell_Type_levels[1:length(Cell_Type_levels)]), stringsAsFactors=F)
        colnames(df.color_Cell_Type)<-c("Cell_Type","colors")
        
        df.color_Cell_Type$Cell_Type<-factor(df.color_Cell_Type$Cell_Type,
                                             levels=levels(h2_h3_df$Cell_Type),
                                             ordered=T)
        
        cat("df.color_Cell_Type_\n")
        cat(str(df.color_Cell_Type))
        cat("\n")
        
       
        
        
        
        h2 <- ggplot(data=h2_h3_df, aes(x=MOCK_COORD,
                                        y=master_sample,
                                        fill = Cell_Type))+
          geom_tile(color = "white")+
          scale_fill_manual(values = df.color_Cell_Type$colors,
                            drop=F)+
          ggeasy::easy_center_title()+
          theme_minimal()+ # minimal theme
          scale_x_discrete(name=element_blank(), drop =F) +
          scale_y_discrete(name=NULL, drop=F)+
          theme(legend.position = "hidden", legend.title = element_blank())+
          theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
          theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
        
        
        h3 <- ggplot(data=h2_h3_df, aes(x=master_sample,
                                        y=MOCK_COORD,
                                        fill = Cell_Type))+
          geom_tile(color = "white")+
          scale_fill_manual(values = df.color_Cell_Type$colors,
                            drop=F)+
          ggeasy::easy_center_title()+
          theme_minimal()+ # minimal theme
          scale_x_discrete(name=element_blank(), drop =F) +
          scale_y_discrete(name=NULL, drop=F)+
          theme(legend.position = "hidden", legend.title = element_blank())+
          theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
          theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
        
        
        #################### Melt the correlation matrix
        
        cat("lower_tri\n")
        str(lower_tri)
        cat("\n")
        
        melted_cormat <- melt(lower_tri, na.rm = TRUE)
        
        cat("melted_cormat\n")
        str(melted_cormat)
        cat("\n")
        
        colnames(melted_cormat)[which(colnames(melted_cormat) == "Var2")]<-"master_sample"
        
        
        melted_cormat<-merge(melted_cormat,
                             h2_h3_df,
                             by="master_sample")
        
        
        
        cat("melted_cormat\n")
        str(melted_cormat)
        cat("\n")
        
        
        levels_Cell_Type<-levels(melted_cormat$Cell_Type)
        
        vector_COORDS<-0
        
        FINAL_COORDS<-NULL
        
        for(i in 1:length(levels_Cell_Type))
        {
          
          Cell_Type_sel<-levels_Cell_Type[i]
          
          cat("--->\t")
          cat(sprintf(as.character(Cell_Type_sel)))
          cat("\n")
          
          melted_cormat_sel<-melted_cormat[which(melted_cormat$Cell_Type == Cell_Type_sel),]
          
          cat("melted_cormat_sel\n")
          str(melted_cormat_sel)
          cat("\n")
          
          melted_cormat_sel_master_samples<-unique(as.character(melted_cormat_sel$master_sample))
          
          cat("melted_cormat_sel_master_samples\n")
          str(melted_cormat_sel_master_samples)
          cat("\n")
          
          COORD_X<-length(melted_cormat_sel_master_samples)
          
          cat("COORD_X--->\t")
          cat(sprintf(as.character(COORD_X)))
          cat("\n")
          
          COORD_DEF<-sum(vector_COORDS,COORD_X)
          
          cat("COORD_DEF--->\t")
          cat(sprintf(as.character(COORD_DEF)))
          cat("\n")
          
          vector_COORDS[i]<-COORD_X
          
          
          # COORD_Y<-length(melted_cormat_sel_master_samples)
          # 
          # cat("COORD_Y--->\t")
          # cat(sprintf(as.character(COORD_Y)))
          # cat("\n")
          
          FINAL_COORDS[i]<-as.numeric(COORD_DEF)
          
        }
        
        FINAL_COORDS<-FINAL_COORDS+0.5
        
        FINAL_COORDS<-FINAL_COORDS[-length(FINAL_COORDS)]
        
        cat("FINAL_COORDS--->\t")
        cat(sprintf(as.character(FINAL_COORDS)))
        cat("\n")
        
        #### heatmap graph
        
        breaks.Rank<-seq(-1,1,by=0.5)
        labels.Rank<-as.character(breaks.Rank)
        
        
        ggheatmap <- ggplot(melted_cormat, aes(master_sample, Var1, fill = value))+
          geom_tile(color = "white")+
          scale_fill_gradient2(name = paste("Pearson","Correlation","Coefficient", sep="\n"),
                               low = "red", high = "green",mid="white",midpoint=0,
                               na.value = NA,
                               breaks=breaks.Rank,
                               labels=labels.Rank,
                               limits=c(breaks.Rank[1],
                                        breaks.Rank[length(breaks.Rank)]))+
          theme_minimal()+ # minimal theme
          scale_x_discrete(name=element_blank(), drop =F) +
          scale_y_discrete(name=NULL, drop=F)+
          theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                           size = 12, hjust = 1))+
          geom_hline(yintercept = FINAL_COORDS, linetype=3, color="black")+
          geom_vline(xintercept = FINAL_COORDS, linetype=3, color="black")+
          coord_fixed()
        
        
        
        ggheatmap2<-ggheatmap+
          theme(axis.text.y  = element_blank())+
          theme(axis.text.x  = element_blank())+
          theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.ticks = element_blank(),
            legend.justification = c(1, 0),
            legend.position = c(0.8, 0.2),
            legend.direction = "horizontal")+
          guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                       title.position = "top", title.hjust = 0.5))+
          theme(text=element_text(size=24,  family="sans"))
        
        
        h3<-plot_grid(NULL,h3,NULL,
                      nrow = 1,
                      ncol = 3,
                      align = "hv",
                      rel_widths = c(0.18, 1,0.1),
                      rel_heights = 1)
        
        graph_DEF<-plot_grid(h2,ggheatmap2,
                             nrow = 1,
                             ncol = 2,
                             align = "hv",
                             rel_widths = c(0.065, 1),
                             rel_heights = 1)
        
        graph_DEF<-plot_grid(h3,graph_DEF,
                             nrow = 2,
                             ncol = 1,
                             align = "hv",
                             rel_widths = 1,
                             rel_heights = c(0.065,1))
        
        
        
        title <- ggdraw() + draw_label(paste("Corr plot",features_sel,sep="  "))
        
        
        
        
        cat("svg_graph\n")
        
        setwd(path2)
        
        
        
        svgname<-paste("Graph_corr_SUBSET_",features_sel,".svg", sep='')
        makesvg = TRUE
        
        if (makesvg == TRUE)
        {
          
          ggsave(svgname, plot= graph_DEF,
                 device="svg",
                 height=10, width=12)
        }
        
        
        graph_DEF2<-plot_grid(h2,ggheatmap,
                              nrow = 2,
                              ncol = 2,
                              align = "hv",
                              rel_widths = c(0.065, 1),
                              rel_heights = c(1,0.05))
        
        graph_DEF3<-plot_grid(title,graph_DEF2,
                              ncol=1,
                              rel_heights = c(0.1,1))
        
        
        svgname<-paste("Graph_corr_SUBSET_NAMED_",features_sel,".svg", sep='')
        makesvg = TRUE
        
        if (makesvg == TRUE)
        {
          
          ggsave(svgname, plot= graph_DEF3,
                 device="svg",
                 height=10, width=12)
        }
        
      
      
      
    
    }# features_sel == "LogFC"
  }#i features
  
  # #################################################################################################################################################################################################################################################################################
  # quit(status=1)
}



# 



get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
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
  
  ggridge.graphs(opt)
  corr.graphs(opt)
  
  
}


###########################################################################

system.time( main() )
