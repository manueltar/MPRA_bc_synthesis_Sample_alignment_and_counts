
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

cummulative_CLASSIF_with_ctrls = function(option_list)
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
  
  # ##############################################################
  # quit(status = 1)
  
  
  #### READ DATA----
  
  setwd(out2)
  
  filename<-paste('Element_collapse','.rds',sep='')
  
  KEY_collapse<-readRDS(file=filename)
  
  cat("KEY_collapse_0\n")
  cat(str(KEY_collapse))
  cat("\n")
  cat(str(unique(KEY_collapse$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(KEY_collapse$Label_2)))))
  cat("\n")
  cat(sprintf(as.character(summary(KEY_collapse$Label_2))))
  cat("\n")
  
  
  check_NCGR<-KEY_collapse[which(KEY_collapse$Label_2 == "NCGR"),]
  
  cat("check_NCGR_0\n")
  cat(str(check_NCGR))
  cat("\n")
  cat(str(unique(check_NCGR$VAR)))
  cat("\n")
  cat(str(unique(check_NCGR$KEY_Plus_carried_variants)))
  cat("\n")
  cat(sprintf(as.character(names(summary(check_NCGR$Label_2)))))
  cat("\n")
  cat(sprintf(as.character(summary(check_NCGR$Label_2))))
  cat("\n")
  
  
  
  setwd(out)
  
  filename<-paste('df_Cell_colors','.rds',sep='')
  
  df_Cell_colors<-readRDS(file=filename)
  
  cat("df_Cell_colors_0\n")
  cat(str(df_Cell_colors))
  cat("\n")
  
  # #########################################
  # quit(status = 1)
  
  ############ ALL CT CLASSIFICATION ----
  
  
  array_Elements<-KEY_collapse$KEY_Plus_carried_variants
  
  cat("array_Elements\n")
  cat(str(array_Elements))
  cat("\n")
  
  
  check_NCGR_2<-array_Elements[which(array_Elements%in%check_NCGR$KEY_Plus_carried_variants)]
  
  cat("check_NCGR\n")
  cat(str(check_NCGR))
  cat("\n")
  
  
  
  
  Condition_DEBUG <- 0
  
  list_ALL_CT<-list()
  
 
  
  for(i in 1:length(array_Elements))
  {
    array_Elements_sel<-array_Elements[i]
    
    cat("---------------------------------------->\t")
    cat(sprintf(as.character(array_Elements_sel)))
    cat("\n")
    
    if(array_Elements_sel == "Element_87;NCGR")
    {
      
      Condition_DEBUG <- 1
    }
    
  
    
    KEY_collapse_sel<-KEY_collapse[which(KEY_collapse$KEY_Plus_carried_variants == array_Elements_sel),]
    
    if(Condition_DEBUG == 1)
    {
      cat("KEY_collapse_sel\n")
      cat(str(KEY_collapse_sel))
      cat("\n")
      cat(sprintf(as.character(colnames(KEY_collapse_sel))))
      cat("\n")
    }
    
    
    indx.int <- which(colnames(KEY_collapse_sel)%in%c("KEY_Plus_carried_variants","KEY","carried_variants","VAR","chr","Label","Label_2","VEP_DEF_LABELS","factor4_CLASS"))
    
    KEY_subset<-unique(KEY_collapse_sel[,indx.int])
    
    if(Condition_DEBUG == 1)
    {
      cat("KEY_subset\n")
      cat(str(KEY_subset))
      cat("\n")
    }
    
    ACTIVE_TILES_sel<-ACTIVE_TILES[which(ACTIVE_TILES$KEY_Plus_carried_variants%in%KEY_collapse_sel$KEY_Plus_carried_variants),]
    
    if(Condition_DEBUG == 1)
    {
      cat("ACTIVE_TILES_sel\n")
      cat(str(ACTIVE_TILES_sel))
      cat("\n")
    }
    
   
    
    if(dim(ACTIVE_TILES_sel)[1] >0)
    {
      
      ### ALL_CT_enhancer
      
      enhancer_labels<-c("enhancer")
      
      
      ACTIVE_TILES_sel_enhancer<-ACTIVE_TILES_sel[which(ACTIVE_TILES_sel$CLASS_enhancer%in%enhancer_labels),]
      
      
      # cat("ACTIVE_TILES_sel_enhancer\n")
      # cat(str(ACTIVE_TILES_sel_enhancer))
      # cat("\n")
      
      if(dim(ACTIVE_TILES_sel_enhancer)[1] >0)
      {
        
        
        EVERY_CT_enhancer<-ACTIVE_TILES_sel_enhancer
        
        # cat("EVERY_CT_enhancer\n")
        # cat(str(EVERY_CT_enhancer))
        # cat("\n")
        
        EVERY_CT_enhancer.dt<-data.table(EVERY_CT_enhancer, key=c("KEY_Plus_carried_variants","Cell_Type"))
        
        
        EVERY_CT_enhancer_Fq<-as.data.frame(EVERY_CT_enhancer.dt[,.(enhancer_CLASS_TILES=.N),
                                                                 by=key(EVERY_CT_enhancer.dt)], stringsAsFactors=F)
        
        
        
        
        
        # cat("EVERY_CT_enhancer_Fq_\n")
        # cat(str(EVERY_CT_enhancer_Fq))
        # cat("\n")
        
        
        
        
        ALL_CT_enhancer.dt<-data.table(EVERY_CT_enhancer, key=c("REAL_TILE_Plus_carried_variants"))
        
        
        ALL_CT_enhancer_Fq<-as.data.frame(ALL_CT_enhancer.dt[,.(enhancer_CLASS_TILES=.N),
                                                             by=key(ALL_CT_enhancer.dt)], stringsAsFactors=F)
        
        
        
        ALL_enhancer<-as.data.frame(cbind(dim(ALL_CT_enhancer_Fq)[1],"ALL_CT"), stringAsFactors=F)
        
        colnames(ALL_enhancer)<-c("enhancer_CLASS_TILES","Cell_Type")
        
        ALL_enhancer$KEY_Plus_carried_variants<-array_Elements_sel
        
        # cat("ALL_enhancer_\n")
        # cat(str(ALL_enhancer))
        # cat("\n")
        
        FINAL_enhancer<-rbind(EVERY_CT_enhancer_Fq,ALL_enhancer)
        
        # cat("FINAL_enhancer_\n")
        # cat(str(FINAL_enhancer))
        # cat("\n")
        
        
        # ############################################################
        # quit(status = 1)   
        
        
        
        
      }#dim(ACTIVE_TILES_sel_enhancer)[1] >0
      else{
        
        FINAL_enhancer<-as.data.frame(cbind(rep(array_Elements_sel,5),rep(0,5),c("K562","CHRF","HL60","THP1","ALL_CT")))
        colnames(FINAL_enhancer)<-c("KEY_Plus_carried_variants","enhancer_CLASS_TILES","Cell_Type")
        
        
      }
      
      # cat("FINAL_enhancer_\n")
      # cat(str(FINAL_enhancer))
      # cat("\n")
      
      ### ALL_CT_ASE
      
      ASE_labels<-c("ASE")
      
      
      ACTIVE_TILES_sel_ASE<-ACTIVE_TILES_sel[which(ACTIVE_TILES_sel$CLASS_ASE%in%ASE_labels),]
      
      
      # cat("ACTIVE_TILES_sel_ASE\n")
      # cat(str(ACTIVE_TILES_sel_ASE))
      # cat("\n")
      
      if(dim(ACTIVE_TILES_sel_ASE)[1] >0)
      {
        
        
        EVERY_CT_ASE<-ACTIVE_TILES_sel_ASE
        
        # cat("EVERY_CT_ASE\n")
        # cat(str(EVERY_CT_ASE))
        # cat("\n")
        
        EVERY_CT_ASE.dt<-data.table(EVERY_CT_ASE, key=c("KEY_Plus_carried_variants","Cell_Type"))
        
        
        EVERY_CT_ASE_Fq<-as.data.frame(EVERY_CT_ASE.dt[,.(ASE_CLASS_TILES=.N),
                                                       by=key(EVERY_CT_ASE.dt)], stringsAsFactors=F)
        
        
        
        
        
        # cat("EVERY_CT_ASE_Fq_\n")
        # cat(str(EVERY_CT_ASE_Fq))
        # cat("\n")
        
        
        
        
        ALL_CT_ASE.dt<-data.table(EVERY_CT_ASE, key=c("REAL_TILE_Plus_carried_variants"))
        
        
        ALL_CT_ASE_Fq<-as.data.frame(ALL_CT_ASE.dt[,.(ASE_CLASS_TILES=.N),
                                                   by=key(ALL_CT_ASE.dt)], stringsAsFactors=F)
        
        
        
        ALL_ASE<-as.data.frame(cbind(dim(ALL_CT_ASE_Fq)[1],"ALL_CT"), stringAsFactors=F)
        
        colnames(ALL_ASE)<-c("ASE_CLASS_TILES","Cell_Type")
        
        ALL_ASE$KEY_Plus_carried_variants<-array_Elements_sel
        
        # cat("ALL_ASE_\n")
        # cat(str(ALL_ASE))
        # cat("\n")
        
        FINAL_ASE<-rbind(EVERY_CT_ASE_Fq,ALL_ASE)
        
        # cat("FINAL_ASE_\n")
        # cat(str(FINAL_ASE))
        # cat("\n")
        
        
        # ############################################################
        # quit(status = 1)   
        
        
        
        
      }#dim(ACTIVE_TILES_sel_ASE)[1] >0
      else{
        
        FINAL_ASE<-as.data.frame(cbind(rep(array_Elements_sel,5),rep(0,5),c("K562","CHRF","HL60","THP1","ALL_CT")))
        colnames(FINAL_ASE)<-c("KEY_Plus_carried_variants","ASE_CLASS_TILES","Cell_Type")
        
        
      }
      
      # cat("FINAL_ASE_\n")
      # cat(str(FINAL_ASE))
      # cat("\n")
      
      ### ALL_CT_E_Plus_ASE
      
      E_Plus_ASE_labels<-c("E_Plus_ASE")
      
      
      ACTIVE_TILES_sel_E_Plus_ASE<-ACTIVE_TILES_sel[which(ACTIVE_TILES_sel$CLASS_E_Plus_ASE%in%E_Plus_ASE_labels),]
      
      
      # cat("ACTIVE_TILES_sel_E_Plus_ASE\n")
      # cat(str(ACTIVE_TILES_sel_E_Plus_ASE))
      # cat("\n")
      
      if(dim(ACTIVE_TILES_sel_E_Plus_ASE)[1] >0)
      {
        
        
        EVERY_CT_E_Plus_ASE<-ACTIVE_TILES_sel_E_Plus_ASE
        
        # cat("EVERY_CT_E_Plus_ASE\n")
        # cat(str(EVERY_CT_E_Plus_ASE))
        # cat("\n")
        
        EVERY_CT_E_Plus_ASE.dt<-data.table(EVERY_CT_E_Plus_ASE, key=c("KEY_Plus_carried_variants","Cell_Type"))
        
        
        EVERY_CT_E_Plus_ASE_Fq<-as.data.frame(EVERY_CT_E_Plus_ASE.dt[,.(E_Plus_ASE_CLASS_TILES=.N),
                                                                     by=key(EVERY_CT_E_Plus_ASE.dt)], stringsAsFactors=F)
        
        
        
        
        
        # cat("EVERY_CT_E_Plus_ASE_Fq_\n")
        # cat(str(EVERY_CT_E_Plus_ASE_Fq))
        # cat("\n")
        
        
        
        
        ALL_CT_E_Plus_ASE.dt<-data.table(EVERY_CT_E_Plus_ASE, key=c("REAL_TILE_Plus_carried_variants"))
        
        
        ALL_CT_E_Plus_ASE_Fq<-as.data.frame(ALL_CT_E_Plus_ASE.dt[,.(E_Plus_ASE_CLASS_TILES=.N),
                                                                 by=key(ALL_CT_E_Plus_ASE.dt)], stringsAsFactors=F)
        
        
        
        ALL_E_Plus_ASE<-as.data.frame(cbind(dim(ALL_CT_E_Plus_ASE_Fq)[1],"ALL_CT"), stringAsFactors=F)
        
        colnames(ALL_E_Plus_ASE)<-c("E_Plus_ASE_CLASS_TILES","Cell_Type")
        
        ALL_E_Plus_ASE$KEY_Plus_carried_variants<-array_Elements_sel
        
        # cat("ALL_E_Plus_ASE_\n")
        # cat(str(ALL_E_Plus_ASE))
        # cat("\n")
        
        FINAL_E_Plus_ASE<-rbind(EVERY_CT_E_Plus_ASE_Fq,ALL_E_Plus_ASE)
        
        # cat("FINAL_E_Plus_ASE_\n")
        # cat(str(FINAL_E_Plus_ASE))
        # cat("\n")
        
        
        # ############################################################
        # quit(status = 1)   
        
        
        
        
      }#dim(ACTIVE_TILES_sel_E_Plus_ASE)[1] >0
      else{
        
        FINAL_E_Plus_ASE<-as.data.frame(cbind(rep(array_Elements_sel,5),rep(0,5),c("K562","CHRF","HL60","THP1","ALL_CT")))
        colnames(FINAL_E_Plus_ASE)<-c("KEY_Plus_carried_variants","E_Plus_ASE_CLASS_TILES","Cell_Type")
        
        
      }
      
      # cat("FINAL_E_Plus_ASE_\n")
      # cat(str(FINAL_E_Plus_ASE))
      # cat("\n")
      
      #### ALL 
      
      ALL_DEF<-merge(FINAL_enhancer,FINAL_ASE,
                     by=c("KEY_Plus_carried_variants","Cell_Type"),
                     all=T)
      
      ALL_DEF<-merge(ALL_DEF,FINAL_E_Plus_ASE,
                     by=c("KEY_Plus_carried_variants","Cell_Type"),
                     all=T)
      
      
      ALL_DEF$enhancer_CLASS_TILES<-as.numeric(ALL_DEF$enhancer_CLASS_TILES)
      ALL_DEF$ASE_CLASS_TILES<-as.numeric(ALL_DEF$ASE_CLASS_TILES)
      ALL_DEF$E_Plus_ASE_CLASS_TILES<-as.numeric(ALL_DEF$E_Plus_ASE_CLASS_TILES)
      
      ALL_DEF[is.na(ALL_DEF)]<-0
      
      # cat("ALL_DEF_0\n")
      # cat(str(ALL_DEF))
      # cat("\n")
      
      
      ALL_DEF<-merge(KEY_subset,
                     ALL_DEF,
                     by="KEY_Plus_carried_variants",
                     all.y=T)
      
      # cat("ALL_DEF_1\n")
      # cat(str(ALL_DEF))
      # cat("\n")
      
      
      # ######################################################################
      # quit(status = 1)
      
      check_1<-max(ALL_DEF$enhancer_CLASS_TILES)
      
      if(check_1 > 5)
      {
        
        ######################################################################
        quit(status = 1)
      }
      
      
      check_2<-max(ALL_DEF$ASE_CLASS_TILES)
      
      if(check_2 > 5)
      {
        setwd(out)
        write.table(ACTIVE_TILES_sel_ASE, file="test.tsv",sep="\t",quote=F,row.names = F)
        
        ######################################################################
        quit(status = 1)
      }
      
      check_3<-max(ALL_DEF$E_Plus_ASE_CLASS_TILES)
      
      if(check_3 > 5)
      {
        
        ######################################################################
        quit(status = 1)
      }
      
      
      
      
      list_ALL_CT[[i]]<-ALL_DEF
      
    }else{
      
      FINAL_enhancer<-as.data.frame(cbind(rep(array_Elements_sel,5),rep(0,5),c("K562","CHRF","HL60","THP1","ALL_CT")))
      colnames(FINAL_enhancer)<-c("KEY_Plus_carried_variants","enhancer_CLASS_TILES","Cell_Type")
      
      FINAL_ASE<-as.data.frame(cbind(rep(array_Elements_sel,5),rep(0,5),c("K562","CHRF","HL60","THP1","ALL_CT")))
      colnames(FINAL_ASE)<-c("KEY_Plus_carried_variants","ASE_CLASS_TILES","Cell_Type")
      
      FINAL_E_Plus_ASE<-as.data.frame(cbind(rep(array_Elements_sel,5),rep(0,5),c("K562","CHRF","HL60","THP1","ALL_CT")))
      colnames(FINAL_E_Plus_ASE)<-c("KEY_Plus_carried_variants","E_Plus_ASE_CLASS_TILES","Cell_Type")
      
      #### ALL 
      
      ALL_DEF<-merge(FINAL_enhancer,FINAL_ASE,
                     by=c("KEY_Plus_carried_variants","Cell_Type"),
                     all=T)
      
      ALL_DEF<-merge(ALL_DEF,FINAL_E_Plus_ASE,
                     by=c("KEY_Plus_carried_variants","Cell_Type"),
                     all=T)
      
      
      ALL_DEF$enhancer_CLASS_TILES<-as.numeric(ALL_DEF$enhancer_CLASS_TILES)
      ALL_DEF$ASE_CLASS_TILES<-as.numeric(ALL_DEF$ASE_CLASS_TILES)
      ALL_DEF$E_Plus_ASE_CLASS_TILES<-as.numeric(ALL_DEF$E_Plus_ASE_CLASS_TILES)
      
      ALL_DEF[is.na(ALL_DEF)]<-0
      
      if(Condition_DEBUG == 1)
      {
        cat("ALL_DEF_0\n")
        cat(str(ALL_DEF))
        cat("\n")
      }
      
      
      ALL_DEF<-merge(KEY_subset,
                     ALL_DEF,
                     by="KEY_Plus_carried_variants",
                     all.y=T)
      
      
      list_ALL_CT[[i]]<-ALL_DEF
      
      
      }# dim(ACTIVE_TILES_sel)[1] >0
    
    if(array_Elements_sel == "Element_87;NCGR")
    {
      # quit(status = 1)
    }
  }#i
  
  
  
  ALL_CT_df = unique(as.data.frame(data.table::rbindlist(list_ALL_CT, fill=T), stringsAsFactors=F))
  
  
  
  cat("ALL_CT_df_0\n")
  cat(str(ALL_CT_df))
  cat("\n")
  
  ALL_CT_df$Cell_Type<-factor(ALL_CT_df$Cell_Type,
                              levels=c("ALL_CT","K562","CHRF","HL60","THP1"),
                              ordered=T)
  
  ##################### CORRECTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ALL_CT_df$ASE_CLASS_TILES[which(ALL_CT_df$ASE_CLASS_TILES < ALL_CT_df$E_Plus_ASE_CLASS_TILES)]<-ALL_CT_df$E_Plus_ASE_CLASS_TILES[which(ALL_CT_df$ASE_CLASS_TILES < ALL_CT_df$E_Plus_ASE_CLASS_TILES)]
  
  cat("ALL_CT_df_1\n")
  cat(str(ALL_CT_df))
  cat("\n")
  cat(sprintf(as.character(names(summary(ALL_CT_df$Label_2)))))
  cat("\n")
  cat(sprintf(as.character(summary(ALL_CT_df$Label_2))))
  cat("\n")
  
  path5<-paste(out2,'cummulative_plots','/', sep='')
  
  cat("path5\n")
  cat(sprintf(as.character(path5)))
  cat("\n")
  
  
  if (file.exists(path5)){
    
    
    
    
  } else {
    dir.create(file.path(path5))
    
  }
  
  setwd(path5)
  
  saveRDS(ALL_CT_df,file="cummulative_classification_with_CTRLS.rds")
  
  
  # quit(status = 1)
}

cummulative_CLASSIF = function(option_list)
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
  
  # ##############################################################
  # quit(status = 1)
  
 
  #### READ DATA----
 
  setwd(out2)
  
  filename<-paste('Element_collapse','.rds',sep='')
  
  KEY_collapse<-readRDS(file=filename)
  
  cat("KEY_collapse_0\n")
  cat(str(KEY_collapse))
  cat("\n")
  
  setwd(out)
  
  filename<-paste('df_Cell_colors','.rds',sep='')
  
  df_Cell_colors<-readRDS(file=filename)
  
  cat("df_Cell_colors_0\n")
  cat(str(df_Cell_colors))
  cat("\n")
  
  
  #### without CTRLS ----
  

  KEY_collapse_EXCLUDED_other_labels<-KEY_collapse[which(KEY_collapse$Label_2 == "ASSAYED_VARIANT"),]

  KEY_collapse_EXCLUDED_other_labels<-droplevels(KEY_collapse_EXCLUDED_other_labels)


  cat("KEY_collapse_EXCLUDED_other_labels\n")
  cat(str(KEY_collapse_EXCLUDED_other_labels))
  cat("\n")
  
  ############ ALL CT CLASSIFICATION ----
  
  
  array_Elements<-KEY_collapse_EXCLUDED_other_labels$KEY_Plus_carried_variants
  
  cat("array_Elements\n")
  cat(str(array_Elements))
  cat("\n")
  
  list_ALL_CT<-list()
  
  
  for(i in 1:length(array_Elements))
  {
    array_Elements_sel<-array_Elements[i]
    
    cat("---------------------------------------->\t")
    cat(sprintf(as.character(array_Elements_sel)))
    cat("\n")
    
    
    
    
    
    KEY_collapse_EXCLUDED_other_labels_sel<-KEY_collapse_EXCLUDED_other_labels[which(KEY_collapse_EXCLUDED_other_labels$KEY_Plus_carried_variants == array_Elements_sel),]
    
    # cat("KEY_collapse_EXCLUDED_other_labels_sel\n")
    # cat(str(KEY_collapse_EXCLUDED_other_labels_sel))
    # cat("\n")
    # cat(sprintf(as.character(colnames(KEY_collapse_EXCLUDED_other_labels_sel))))
    # cat("\n")
    
    
    indx.int <- which(colnames(KEY_collapse_EXCLUDED_other_labels_sel)%in%c("KEY_Plus_carried_variants","KEY","carried_variants","VAR","chr","Label","Label_2","VEP_DEF_LABELS","factor4_CLASS"))

    KEY_subset<-unique(KEY_collapse_EXCLUDED_other_labels_sel[,indx.int])

    # cat("KEY_subset\n")
    # cat(str(KEY_subset))
    # cat("\n")
    
    ACTIVE_TILES_sel<-ACTIVE_TILES[which(ACTIVE_TILES$KEY_Plus_carried_variants%in%KEY_collapse_EXCLUDED_other_labels_sel$KEY_Plus_carried_variants),]
    
    # cat("ACTIVE_TILES_sel\n")
    # cat(str(ACTIVE_TILES_sel))
    # cat("\n")
    
    
    
    if(dim(ACTIVE_TILES_sel)[1] >0)
    {
      
      ### ALL_CT_enhancer
      
      enhancer_labels<-c("enhancer")
      
      
      ACTIVE_TILES_sel_enhancer<-ACTIVE_TILES_sel[which(ACTIVE_TILES_sel$CLASS_enhancer%in%enhancer_labels),]
      
      
      # cat("ACTIVE_TILES_sel_enhancer\n")
      # cat(str(ACTIVE_TILES_sel_enhancer))
      # cat("\n")
      
      if(dim(ACTIVE_TILES_sel_enhancer)[1] >0)
      {
        
        
        EVERY_CT_enhancer<-ACTIVE_TILES_sel_enhancer
        
        # cat("EVERY_CT_enhancer\n")
        # cat(str(EVERY_CT_enhancer))
        # cat("\n")
        
        EVERY_CT_enhancer.dt<-data.table(EVERY_CT_enhancer, key=c("KEY_Plus_carried_variants","Cell_Type"))
        
        
        EVERY_CT_enhancer_Fq<-as.data.frame(EVERY_CT_enhancer.dt[,.(enhancer_CLASS_TILES=.N),
                                                       by=key(EVERY_CT_enhancer.dt)], stringsAsFactors=F)
        
        
        
        
        
        # cat("EVERY_CT_enhancer_Fq_\n")
        # cat(str(EVERY_CT_enhancer_Fq))
        # cat("\n")
        
        
        
        
        ALL_CT_enhancer.dt<-data.table(EVERY_CT_enhancer, key=c("REAL_TILE_Plus_carried_variants"))
        
        
        ALL_CT_enhancer_Fq<-as.data.frame(ALL_CT_enhancer.dt[,.(enhancer_CLASS_TILES=.N),
                                                   by=key(ALL_CT_enhancer.dt)], stringsAsFactors=F)
        
        
        
        ALL_enhancer<-as.data.frame(cbind(dim(ALL_CT_enhancer_Fq)[1],"ALL_CT"), stringAsFactors=F)
        
        colnames(ALL_enhancer)<-c("enhancer_CLASS_TILES","Cell_Type")
        
        ALL_enhancer$KEY_Plus_carried_variants<-array_Elements_sel
        
        # cat("ALL_enhancer_\n")
        # cat(str(ALL_enhancer))
        # cat("\n")
        
        FINAL_enhancer<-rbind(EVERY_CT_enhancer_Fq,ALL_enhancer)
        
        # cat("FINAL_enhancer_\n")
        # cat(str(FINAL_enhancer))
        # cat("\n")
        
        
        # ############################################################
        # quit(status = 1)   
        
        
        
        
      }#dim(ACTIVE_TILES_sel_enhancer)[1] >0
      else{
        
        FINAL_enhancer<-as.data.frame(cbind(rep(array_Elements_sel,5),rep(0,5),c("K562","CHRF","HL60","THP1","ALL_CT")))
        colnames(FINAL_enhancer)<-c("KEY_Plus_carried_variants","enhancer_CLASS_TILES","Cell_Type")
        
       
      }
      
      # cat("FINAL_enhancer_\n")
      # cat(str(FINAL_enhancer))
      # cat("\n")
      
      ### ALL_CT_ASE
      
      ASE_labels<-c("ASE")
      
      
      ACTIVE_TILES_sel_ASE<-ACTIVE_TILES_sel[which(ACTIVE_TILES_sel$CLASS_ASE%in%ASE_labels),]
      
      
      # cat("ACTIVE_TILES_sel_ASE\n")
      # cat(str(ACTIVE_TILES_sel_ASE))
      # cat("\n")
      
      if(dim(ACTIVE_TILES_sel_ASE)[1] >0)
      {
        
      
        EVERY_CT_ASE<-ACTIVE_TILES_sel_ASE
        
        # cat("EVERY_CT_ASE\n")
        # cat(str(EVERY_CT_ASE))
        # cat("\n")
        
        EVERY_CT_ASE.dt<-data.table(EVERY_CT_ASE, key=c("KEY_Plus_carried_variants","Cell_Type"))
        
        
        EVERY_CT_ASE_Fq<-as.data.frame(EVERY_CT_ASE.dt[,.(ASE_CLASS_TILES=.N),
                                                   by=key(EVERY_CT_ASE.dt)], stringsAsFactors=F)
        
        
        
        
        
        # cat("EVERY_CT_ASE_Fq_\n")
        # cat(str(EVERY_CT_ASE_Fq))
        # cat("\n")
        
        
     
        
        ALL_CT_ASE.dt<-data.table(EVERY_CT_ASE, key=c("REAL_TILE_Plus_carried_variants"))
        
        
        ALL_CT_ASE_Fq<-as.data.frame(ALL_CT_ASE.dt[,.(ASE_CLASS_TILES=.N),
                                                   by=key(ALL_CT_ASE.dt)], stringsAsFactors=F)
        
        
        
        ALL_ASE<-as.data.frame(cbind(dim(ALL_CT_ASE_Fq)[1],"ALL_CT"), stringAsFactors=F)
        
        colnames(ALL_ASE)<-c("ASE_CLASS_TILES","Cell_Type")
        
        ALL_ASE$KEY_Plus_carried_variants<-array_Elements_sel
        
        # cat("ALL_ASE_\n")
        # cat(str(ALL_ASE))
        # cat("\n")
        
        FINAL_ASE<-rbind(EVERY_CT_ASE_Fq,ALL_ASE)
        
        # cat("FINAL_ASE_\n")
        # cat(str(FINAL_ASE))
        # cat("\n")
        
        
        # ############################################################
        # quit(status = 1)   
        
        
       
        
      }#dim(ACTIVE_TILES_sel_ASE)[1] >0
      else{
        
        FINAL_ASE<-as.data.frame(cbind(rep(array_Elements_sel,5),rep(0,5),c("K562","CHRF","HL60","THP1","ALL_CT")))
        colnames(FINAL_ASE)<-c("KEY_Plus_carried_variants","ASE_CLASS_TILES","Cell_Type")
        
       
      }
      
      # cat("FINAL_ASE_\n")
      # cat(str(FINAL_ASE))
      # cat("\n")
      
      ### ALL_CT_E_Plus_ASE
      
      E_Plus_ASE_labels<-c("E_Plus_ASE")
      
      
      ACTIVE_TILES_sel_E_Plus_ASE<-ACTIVE_TILES_sel[which(ACTIVE_TILES_sel$CLASS_E_Plus_ASE%in%E_Plus_ASE_labels),]
      
      
      # cat("ACTIVE_TILES_sel_E_Plus_ASE\n")
      # cat(str(ACTIVE_TILES_sel_E_Plus_ASE))
      # cat("\n")
      
      if(dim(ACTIVE_TILES_sel_E_Plus_ASE)[1] >0)
      {
        
        
        EVERY_CT_E_Plus_ASE<-ACTIVE_TILES_sel_E_Plus_ASE
        
        # cat("EVERY_CT_E_Plus_ASE\n")
        # cat(str(EVERY_CT_E_Plus_ASE))
        # cat("\n")
        
        EVERY_CT_E_Plus_ASE.dt<-data.table(EVERY_CT_E_Plus_ASE, key=c("KEY_Plus_carried_variants","Cell_Type"))
        
        
        EVERY_CT_E_Plus_ASE_Fq<-as.data.frame(EVERY_CT_E_Plus_ASE.dt[,.(E_Plus_ASE_CLASS_TILES=.N),
                                                       by=key(EVERY_CT_E_Plus_ASE.dt)], stringsAsFactors=F)
        
        
        
        
        
        # cat("EVERY_CT_E_Plus_ASE_Fq_\n")
        # cat(str(EVERY_CT_E_Plus_ASE_Fq))
        # cat("\n")
        
        
        
        
        ALL_CT_E_Plus_ASE.dt<-data.table(EVERY_CT_E_Plus_ASE, key=c("REAL_TILE_Plus_carried_variants"))
        
        
        ALL_CT_E_Plus_ASE_Fq<-as.data.frame(ALL_CT_E_Plus_ASE.dt[,.(E_Plus_ASE_CLASS_TILES=.N),
                                                   by=key(ALL_CT_E_Plus_ASE.dt)], stringsAsFactors=F)
        
        
        
        ALL_E_Plus_ASE<-as.data.frame(cbind(dim(ALL_CT_E_Plus_ASE_Fq)[1],"ALL_CT"), stringAsFactors=F)
        
        colnames(ALL_E_Plus_ASE)<-c("E_Plus_ASE_CLASS_TILES","Cell_Type")
        
        ALL_E_Plus_ASE$KEY_Plus_carried_variants<-array_Elements_sel
        
        # cat("ALL_E_Plus_ASE_\n")
        # cat(str(ALL_E_Plus_ASE))
        # cat("\n")
        
        FINAL_E_Plus_ASE<-rbind(EVERY_CT_E_Plus_ASE_Fq,ALL_E_Plus_ASE)
        
        # cat("FINAL_E_Plus_ASE_\n")
        # cat(str(FINAL_E_Plus_ASE))
        # cat("\n")
        
        
        # ############################################################
        # quit(status = 1)   
        
        
        
        
      }#dim(ACTIVE_TILES_sel_E_Plus_ASE)[1] >0
      else{
        
        FINAL_E_Plus_ASE<-as.data.frame(cbind(rep(array_Elements_sel,5),rep(0,5),c("K562","CHRF","HL60","THP1","ALL_CT")))
        colnames(FINAL_E_Plus_ASE)<-c("KEY_Plus_carried_variants","E_Plus_ASE_CLASS_TILES","Cell_Type")
        
        
      }
      
      # cat("FINAL_E_Plus_ASE_\n")
      # cat(str(FINAL_E_Plus_ASE))
      # cat("\n")
      
      #### ALL 
      
      ALL_DEF<-merge(FINAL_enhancer,FINAL_ASE,
                     by=c("KEY_Plus_carried_variants","Cell_Type"),
                     all=T)
      
      ALL_DEF<-merge(ALL_DEF,FINAL_E_Plus_ASE,
                     by=c("KEY_Plus_carried_variants","Cell_Type"),
                     all=T)
      
      
      ALL_DEF$enhancer_CLASS_TILES<-as.numeric(ALL_DEF$enhancer_CLASS_TILES)
      ALL_DEF$ASE_CLASS_TILES<-as.numeric(ALL_DEF$ASE_CLASS_TILES)
      ALL_DEF$E_Plus_ASE_CLASS_TILES<-as.numeric(ALL_DEF$E_Plus_ASE_CLASS_TILES)
      
      ALL_DEF[is.na(ALL_DEF)]<-0
      
      # cat("ALL_DEF_0\n")
      # cat(str(ALL_DEF))
      # cat("\n")
      
      
      ALL_DEF<-merge(KEY_subset,
            ALL_DEF,
            by="KEY_Plus_carried_variants",
            all.y=T)
      
      # cat("ALL_DEF_1\n")
      # cat(str(ALL_DEF))
      # cat("\n")
      
      
      # ######################################################################
      # quit(status = 1)
      
      check_1<-max(ALL_DEF$enhancer_CLASS_TILES)
      
      if(check_1 > 5)
      {
        
        ######################################################################
        quit(status = 1)
      }
      
      
      check_2<-max(ALL_DEF$ASE_CLASS_TILES)
      
      if(check_2 > 5)
      {
        setwd(out)
        write.table(ACTIVE_TILES_sel_ASE, file="test.tsv",sep="\t",quote=F,row.names = F)
        
        ######################################################################
        quit(status = 1)
      }
      
      check_3<-max(ALL_DEF$E_Plus_ASE_CLASS_TILES)
      
      if(check_3 > 5)
      {
        
        ######################################################################
        quit(status = 1)
      }
      
      
      
      
      list_ALL_CT[[i]]<-ALL_DEF
      
    }# dim(ACTIVE_TILES_sel)[1] >0
  }#i
  
  
 
    ALL_CT_df = unique(as.data.frame(data.table::rbindlist(list_ALL_CT, fill=T), stringsAsFactors=F))
    
    
    
    cat("ALL_CT_df_0\n")
    cat(str(ALL_CT_df))
    cat("\n")
    
   ALL_CT_df$Cell_Type<-factor(ALL_CT_df$Cell_Type,
                               levels=c("ALL_CT","K562","CHRF","HL60","THP1"),
                               ordered=T)
   
   ##################### CORRECTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   ALL_CT_df$ASE_CLASS_TILES[which(ALL_CT_df$ASE_CLASS_TILES < ALL_CT_df$E_Plus_ASE_CLASS_TILES)]<-ALL_CT_df$E_Plus_ASE_CLASS_TILES[which(ALL_CT_df$ASE_CLASS_TILES < ALL_CT_df$E_Plus_ASE_CLASS_TILES)]
   
   cat("ALL_CT_df_1\n")
   cat(str(ALL_CT_df))
   cat("\n")
   
   path5<-paste(out2,'cummulative_plots','/', sep='')
   
   cat("path5\n")
   cat(sprintf(as.character(path5)))
   cat("\n")
   
   
   if (file.exists(path5)){
     
     
     
     
   } else {
     dir.create(file.path(path5))
     
   }
    
   setwd(path5)
   
   saveRDS(ALL_CT_df,file="cummulative_classification.rds")
  
}

cummulative_PLOTS_enhancer = function(option_list)
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
  
  setwd(out)
  
  
  filename<-paste('df_Cell_colors','.rds',sep='')
  
  df_Cell_colors<-readRDS(file=filename)
  
  cat("df_Cell_colors_0\n")
  cat(str(df_Cell_colors))
  cat("\n")
  
  
  #### rEAD DATA ----
  
  path5<-paste(out2,'cummulative_plots','/', sep='')
  
  cat("path5\n")
  cat(sprintf(as.character(path5)))
  cat("\n")
  
  
  if (file.exists(path5)){
    
    
    
    
  } else {
    dir.create(file.path(path5))
    
  }
  
  setwd(path5)
  
  ALL_CT_df<-readRDS(file="cummulative_classification.rds")
  
  cat("ALL_CT_df_0\n")
  cat(str(ALL_CT_df))
  cat("\n")
  cat(str(unique(ALL_CT_df$VAR)))
  cat("\n")
  
  subset<-ALL_CT_df[which(ALL_CT_df$enhancer_CLASS_TILES >0),]
  
  subset$enhancer_CLASS_TILES<-factor(subset$enhancer_CLASS_TILES,
                                 levels=c("5","4","3","2","1"),
                                 ordered=T)
  
  subset<-subset[order(subset$enhancer_CLASS_TILES),]
  
  cat("subset_0\n")
  cat(str(subset))
  cat("\n")
  
  # enhancer.dt<-data.table(subset, key=c("Cell_Type","enhancer_CLASS_TILES"))
  # 
  # cat("enhancer.dt_0\n")
  # cat(str(enhancer.dt))
  # cat("\n")
  
  
  TABLE_enhancer<-as.data.frame(table(subset$Cell_Type,subset$enhancer_CLASS_TILES), stringsAsFactors=F)
  
  colnames(TABLE_enhancer)[which(colnames(TABLE_enhancer) == "Var1")]<-"Cell_Type"
  colnames(TABLE_enhancer)[which(colnames(TABLE_enhancer) == "Var2")]<-"enhancer_CLASS_TILES"
  colnames(TABLE_enhancer)[which(colnames(TABLE_enhancer) == "Freq")]<-"n"
  
  
  TABLE_enhancer$enhancer_CLASS_TILES<-factor(TABLE_enhancer$enhancer_CLASS_TILES,
                                    levels=c("5","4","3","2","1"),
                                    ordered=T)
  
  TABLE_enhancer$enhancer_CLASS_TILES<-factor(TABLE_enhancer$enhancer_CLASS_TILES,
                                    levels=c("5","4","3","2","1"),
                                    ordered=T)
  
  TABLE_enhancer$Cell_Type<-factor(TABLE_enhancer$Cell_Type,
                              levels=c("ALL_CT","K562","CHRF","HL60","THP1"),
                              ordered=T)
  
  cat("TABLE_enhancer_1\n")
  cat(str(TABLE_enhancer))
  cat("\n")
  
  
  
  TABLE_enhancer.dt<-setDT(TABLE_enhancer)[,Cum.Sum := cumsum(n),by=.(Cell_Type)]
  
  TABLE_enhancer<-as.data.frame(setDT(TABLE_enhancer)[,Total := sum(n),by=.(Cell_Type)], stringsAsFactors=F)
  
  cat("TABLE_enhancer_1\n")
  cat(str(TABLE_enhancer))
  cat("\n")
  
  check<-TABLE_enhancer[(TABLE_enhancer$Cell_Type == "K562"),]
  
  
  cat("check_1\n")
  cat(str(check))
  cat("\n")
  
  
  TABLE_enhancer$Total<-97
  TABLE_enhancer$Perc<-round(100*(TABLE_enhancer$Cum.Sum/TABLE_enhancer$Total),1)
  
  cat("TABLE_enhancer_0\n")
  cat(str(TABLE_enhancer))
  cat("\n")
  
  # setwd(out)
  # 
  # write.table(TABLE_enhancer,file="test.tsv",sep="\t",quote=F, row.names = F)
  # 
  # ##############################################################
  # quit(status = 1)
  
  setwd(path5)
  
  #### PDF Klaudia
  
  v_parameter<-"NA" 
  
  indx_Cell_Type<-which(colnames(TABLE_enhancer) == "Cell_Type")
  indx.yaxis<-which(colnames(TABLE_enhancer) == "Perc")
  mycols <- c("black",df_Cell_colors$colors)
  
  
  
  pdfname<-paste("Cummulative_frequency_","enhancer_activity",".pdf",sep='')
  pdf(file=pdfname, width=5, height=4, pointsize=12)
  
  par(mai=c(0.9,0.9,0.3,0.2))
  lab <- as.character(unique(TABLE_enhancer[,indx_Cell_Type]))
  
  cat("lab\n")
  cat(sprintf(as.character(lab)))
  cat("\n")
  
  
  plot(rev(TABLE_enhancer$enhancer_CLASS_TILES), TABLE_enhancer[,indx.yaxis],
       ty="n", xlab="TILES with enhancer activity",
       ylab="Cummulative % of variants",
       axes=F, cex.lab=1.2, cex.lab=1.3, ylim=c(0,100))
  
  
  
  # cat("Hello_world1\n")
  
  # points(TABLE_enhancer$enhancer_CLASS_TILES, TABLE_enhancer[,indx.yaxis], col="darkgrey", pch=19)
  
  # cat("Hello_world2\n")
  # 
  # 
  for (i in 1:length(lab))  {
    
    cat("Hello_world3\t")
    cat(sprintf(as.character(lab[i])))
    cat("\n")
    
    ind <- which(TABLE_enhancer[,indx_Cell_Type]==lab[i])
    
    cat("------------------->ind\n")
    cat(str(ind))
    cat("\n")
    color_sel<-mycols[i]
    
    cat(sprintf(as.character(color_sel)))
    cat("\n")
    
    points(rev(TABLE_enhancer$enhancer_CLASS_TILES[ind]), TABLE_enhancer[,indx.yaxis][ind],
           pch=19, col=color_sel)
    lines(rev(TABLE_enhancer$enhancer_CLASS_TILES[ind]), TABLE_enhancer[,indx.yaxis][ind],
          lty=1,lwd=3, col=color_sel)
    
    cat("END_3\n")
    
    # points(TABLE_enhancer$enhancer_CLASS_TILES[ind], TABLE_enhancer[,indx.yaxis][ind], pch=19, col=color_sel)
  }
  
  cat("END_4\n")
  
  
  legend("topright", legend=lab, fill=mycols, border=mycols, bty="n")
  axis(1, at=seq(0,100))
  axis(2, las=1)
  
  dev.off()
  
  
  setwd(path5)
  
  write.table(TABLE_enhancer, file="Cummulative_plot_enhancer_activity_table.tsv",sep="\t",quote=F,row.names = F)
  
  
  
}

cummulative_PLOTS_ASE = function(option_list)
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
  
  setwd(out)
  
  
  filename<-paste('df_Cell_colors','.rds',sep='')
  
  df_Cell_colors<-readRDS(file=filename)
  
  cat("df_Cell_colors_0\n")
  cat(str(df_Cell_colors))
  cat("\n")
  
  
  #### rEAD DATA ----
  
  path5<-paste(out2,'cummulative_plots','/', sep='')
  
  cat("path5\n")
  cat(sprintf(as.character(path5)))
  cat("\n")
  
  
  if (file.exists(path5)){
    
    
    
    
  } else {
    dir.create(file.path(path5))
    
  }
  
  setwd(path5)
  
  ALL_CT_df<-readRDS(file="cummulative_classification.rds")
  
  cat("ALL_CT_df_0\n")
  cat(str(ALL_CT_df))
  cat("\n")
  cat(str(unique(ALL_CT_df$VAR)))
  cat("\n")
  
  subset<-ALL_CT_df[which(ALL_CT_df$ASE_CLASS_TILES >0),]
  
  subset$ASE_CLASS_TILES<-factor(subset$ASE_CLASS_TILES,
                                      levels=c("5","4","3","2","1"),
                                      ordered=T)
  
  subset<-subset[order(subset$ASE_CLASS_TILES),]
  
  cat("subset_0\n")
  cat(str(subset))
  cat("\n")
  
  # ASE.dt<-data.table(subset, key=c("Cell_Type","ASE_CLASS_TILES"))
  # 
  # cat("ASE.dt_0\n")
  # cat(str(ASE.dt))
  # cat("\n")
  
  
  TABLE_ASE<-as.data.frame(table(subset$Cell_Type,subset$ASE_CLASS_TILES), stringsAsFactors=F)
  
  colnames(TABLE_ASE)[which(colnames(TABLE_ASE) == "Var1")]<-"Cell_Type"
  colnames(TABLE_ASE)[which(colnames(TABLE_ASE) == "Var2")]<-"ASE_CLASS_TILES"
  colnames(TABLE_ASE)[which(colnames(TABLE_ASE) == "Freq")]<-"n"
  
  
  TABLE_ASE$ASE_CLASS_TILES<-factor(TABLE_ASE$ASE_CLASS_TILES,
                                    levels=c("5","4","3","2","1"),
                                    ordered=T)
  
  TABLE_ASE$ASE_CLASS_TILES<-factor(TABLE_ASE$ASE_CLASS_TILES,
                                 levels=c("5","4","3","2","1"),
                                 ordered=T)
  
  TABLE_ASE$Cell_Type<-factor(TABLE_ASE$Cell_Type,
                              levels=c("ALL_CT","K562","CHRF","HL60","THP1"),
                              ordered=T)
  
  cat("TABLE_ASE_1\n")
  cat(str(TABLE_ASE))
  cat("\n")
  
  
  
  TABLE_ASE.dt<-setDT(TABLE_ASE)[,Cum.Sum := cumsum(n),by=.(Cell_Type)]
  
  TABLE_ASE<-as.data.frame(setDT(TABLE_ASE)[,Total := sum(n),by=.(Cell_Type)], stringsAsFactors=F)
  
  cat("TABLE_ASE_1\n")
  cat(str(TABLE_ASE))
  cat("\n")
  
  check<-TABLE_ASE[(TABLE_ASE$Cell_Type == "K562"),]
  
  
  cat("check_1\n")
  cat(str(check))
  cat("\n")
 
  
  TABLE_ASE$Total<-97
  TABLE_ASE$Perc<-round(100*(TABLE_ASE$Cum.Sum/TABLE_ASE$Total),1)
  
  cat("TABLE_ASE_0\n")
  cat(str(TABLE_ASE))
  cat("\n")
  
  # setwd(out)
  # 
  # write.table(TABLE_ASE,file="test.tsv",sep="\t",quote=F, row.names = F)
  # 
  # ##############################################################
  # quit(status = 1)
  
  setwd(path5)
  
  #### PDF Klaudia
  
  v_parameter<-"NA" 
  
  indx_Cell_Type<-which(colnames(TABLE_ASE) == "Cell_Type")
  indx.yaxis<-which(colnames(TABLE_ASE) == "Perc")
  mycols <- c("black",df_Cell_colors$colors)
  
  
  
  pdfname<-paste("Cummulative_frequency_","ASE_activity",".pdf",sep='')
  pdf(file=pdfname, width=5, height=4, pointsize=12)
  
  par(mai=c(0.9,0.9,0.3,0.2))
  lab <- as.character(unique(TABLE_ASE[,indx_Cell_Type]))
  
  cat("lab\n")
  cat(sprintf(as.character(lab)))
  cat("\n")
  
  
  plot(rev(TABLE_ASE$ASE_CLASS_TILES), TABLE_ASE[,indx.yaxis],
       ty="n", xlab="TILES with ASE activity",
       ylab="Cummulative % of variants",
       axes=F, cex.lab=1.2, cex.lab=1.3, ylim=c(0,100))
  
  
  
  # cat("Hello_world1\n")
  
  # points(TABLE_ASE$ASE_CLASS_TILES, TABLE_ASE[,indx.yaxis], col="darkgrey", pch=19)
  
  # cat("Hello_world2\n")
  # 
  # 
  for (i in 1:length(lab))  {
    
    cat("Hello_world3\t")
    cat(sprintf(as.character(lab[i])))
    cat("\n")
    
    ind <- which(TABLE_ASE[,indx_Cell_Type]==lab[i])
    
    cat("------------------->ind\n")
    cat(str(ind))
    cat("\n")
    color_sel<-mycols[i]
    
    cat(sprintf(as.character(color_sel)))
    cat("\n")
    
    points(rev(TABLE_ASE$ASE_CLASS_TILES[ind]), TABLE_ASE[,indx.yaxis][ind],
           pch=19, col=color_sel)
    lines(rev(TABLE_ASE$ASE_CLASS_TILES[ind]), TABLE_ASE[,indx.yaxis][ind],
          lty=1,lwd=3, col=color_sel)
    
    cat("END_3\n")
    
    # points(TABLE_ASE$ASE_CLASS_TILES[ind], TABLE_ASE[,indx.yaxis][ind], pch=19, col=color_sel)
  }
  
  cat("END_4\n")
  
  
  legend("topright", legend=lab, fill=mycols, border=mycols, bty="n")
  axis(1, at=seq(0,100))
  axis(2, las=1)
  
  dev.off()
  
  
  setwd(path5)
  
  write.table(TABLE_ASE, file="Cummulative_plot_ASE_activity_table.tsv",sep="\t",quote=F,row.names = F)
  
  
  
}

cummulative_PLOTS_E_Plus_ASE = function(option_list)
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
  setwd(out)
  
  
  filename<-paste('df_Cell_colors','.rds',sep='')
  
  df_Cell_colors<-readRDS(file=filename)
  
  cat("df_Cell_colors_0\n")
  cat(str(df_Cell_colors))
  cat("\n")
  
  
  #### rEAD DATA ----
  
  path5<-paste(out2,'cummulative_plots','/', sep='')
  
  cat("path5\n")
  cat(sprintf(as.character(path5)))
  cat("\n")
  
  
  if (file.exists(path5)){
    
    
    
    
  } else {
    dir.create(file.path(path5))
    
  }
  
  setwd(path5)
  
  ALL_CT_df<-readRDS(file="cummulative_classification.rds")
  
  cat("ALL_CT_df_0\n")
  cat(str(ALL_CT_df))
  cat("\n")
  cat(str(unique(ALL_CT_df$VAR)))
  cat("\n")
  
  subset<-ALL_CT_df[which(ALL_CT_df$E_Plus_ASE_CLASS_TILES >0),]
  
  subset$E_Plus_ASE_CLASS_TILES<-factor(subset$E_Plus_ASE_CLASS_TILES,
                                 levels=c("5","4","3","2","1"),
                                 ordered=T)
  
  subset<-subset[order(subset$E_Plus_ASE_CLASS_TILES),]
  
  cat("subset_0\n")
  cat(str(subset))
  cat("\n")
  
  # E_Plus_ASE.dt<-data.table(subset, key=c("Cell_Type","E_Plus_ASE_CLASS_TILES"))
  # 
  # cat("E_Plus_ASE.dt_0\n")
  # cat(str(E_Plus_ASE.dt))
  # cat("\n")
  
  
  TABLE_E_Plus_ASE<-as.data.frame(table(subset$Cell_Type,subset$E_Plus_ASE_CLASS_TILES), stringsAsFactors=F)
  
  colnames(TABLE_E_Plus_ASE)[which(colnames(TABLE_E_Plus_ASE) == "Var1")]<-"Cell_Type"
  colnames(TABLE_E_Plus_ASE)[which(colnames(TABLE_E_Plus_ASE) == "Var2")]<-"E_Plus_ASE_CLASS_TILES"
  colnames(TABLE_E_Plus_ASE)[which(colnames(TABLE_E_Plus_ASE) == "Freq")]<-"n"
  
  
  TABLE_E_Plus_ASE$E_Plus_ASE_CLASS_TILES<-factor(TABLE_E_Plus_ASE$E_Plus_ASE_CLASS_TILES,
                                    levels=c("5","4","3","2","1"),
                                    ordered=T)
  
  TABLE_E_Plus_ASE$E_Plus_ASE_CLASS_TILES<-factor(TABLE_E_Plus_ASE$E_Plus_ASE_CLASS_TILES,
                                    levels=c("5","4","3","2","1"),
                                    ordered=T)
  
  TABLE_E_Plus_ASE$Cell_Type<-factor(TABLE_E_Plus_ASE$Cell_Type,
                              levels=c("ALL_CT","K562","CHRF","HL60","THP1"),
                              ordered=T)
  
  cat("TABLE_E_Plus_ASE_1\n")
  cat(str(TABLE_E_Plus_ASE))
  cat("\n")
  
  
  
  TABLE_E_Plus_ASE.dt<-setDT(TABLE_E_Plus_ASE)[,Cum.Sum := cumsum(n),by=.(Cell_Type)]
  
  TABLE_E_Plus_ASE<-as.data.frame(setDT(TABLE_E_Plus_ASE)[,Total := sum(n),by=.(Cell_Type)], stringsAsFactors=F)
  
  cat("TABLE_E_Plus_ASE_1\n")
  cat(str(TABLE_E_Plus_ASE))
  cat("\n")
  
  check<-TABLE_E_Plus_ASE[(TABLE_E_Plus_ASE$Cell_Type == "K562"),]
  
  
  cat("check_1\n")
  cat(str(check))
  cat("\n")
  
  
  TABLE_E_Plus_ASE$Total<-97
  TABLE_E_Plus_ASE$Perc<-round(100*(TABLE_E_Plus_ASE$Cum.Sum/TABLE_E_Plus_ASE$Total),1)
  
  cat("TABLE_E_Plus_ASE_0\n")
  cat(str(TABLE_E_Plus_ASE))
  cat("\n")
  
  # setwd(out)
  # 
  # write.table(TABLE_E_Plus_ASE,file="test.tsv",sep="\t",quote=F, row.names = F)
  # 
  # ##############################################################
  # quit(status = 1)
  
  setwd(path5)
  
  #### PDF Klaudia
  
  v_parameter<-"NA" 
  
  indx_Cell_Type<-which(colnames(TABLE_E_Plus_ASE) == "Cell_Type")
  indx.yaxis<-which(colnames(TABLE_E_Plus_ASE) == "Perc")
  mycols <- c("black",df_Cell_colors$colors)
  
  
  
  pdfname<-paste("Cummulative_frequency_","E_Plus_ASE_activity",".pdf",sep='')
  pdf(file=pdfname, width=5, height=4, pointsize=12)
  
  par(mai=c(0.9,0.9,0.3,0.2))
  lab <- as.character(unique(TABLE_E_Plus_ASE[,indx_Cell_Type]))
  
  cat("lab\n")
  cat(sprintf(as.character(lab)))
  cat("\n")
  
  
  plot(rev(TABLE_E_Plus_ASE$E_Plus_ASE_CLASS_TILES), TABLE_E_Plus_ASE[,indx.yaxis],
       ty="n", xlab="TILES with E_Plus_ASE activity",
       ylab="Cummulative % of variants",
       axes=F, cex.lab=1.2, cex.lab=1.3, ylim=c(0,100))
  
  
  
  # cat("Hello_world1\n")
  
  # points(TABLE_E_Plus_ASE$E_Plus_ASE_CLASS_TILES, TABLE_E_Plus_ASE[,indx.yaxis], col="darkgrey", pch=19)
  
  # cat("Hello_world2\n")
  # 
  # 
  for (i in 1:length(lab))  {
    
    cat("Hello_world3\t")
    cat(sprintf(as.character(lab[i])))
    cat("\n")
    
    ind <- which(TABLE_E_Plus_ASE[,indx_Cell_Type]==lab[i])
    
    cat("------------------->ind\n")
    cat(str(ind))
    cat("\n")
    color_sel<-mycols[i]
    
    cat(sprintf(as.character(color_sel)))
    cat("\n")
    
    points(rev(TABLE_E_Plus_ASE$E_Plus_ASE_CLASS_TILES[ind]), TABLE_E_Plus_ASE[,indx.yaxis][ind],
           pch=19, col=color_sel)
    lines(rev(TABLE_E_Plus_ASE$E_Plus_ASE_CLASS_TILES[ind]), TABLE_E_Plus_ASE[,indx.yaxis][ind],
          lty=1,lwd=3, col=color_sel)
    
    cat("END_3\n")
    
    # points(TABLE_E_Plus_ASE$E_Plus_ASE_CLASS_TILES[ind], TABLE_E_Plus_ASE[,indx.yaxis][ind], pch=19, col=color_sel)
  }
  
  cat("END_4\n")
  
  
  legend("topright", legend=lab, fill=mycols, border=mycols, bty="n")
  axis(1, at=seq(0,100))
  axis(2, las=1)
  
  dev.off()
  
  
  setwd(path5)
  
  write.table(TABLE_E_Plus_ASE, file="Cummulative_plot_E_Plus_ASE_activity_table.tsv",sep="\t",quote=F,row.names = F)
  
  
  
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
  
  
  cummulative_CLASSIF_with_ctrls(opt)
  cummulative_CLASSIF(opt)
  cummulative_PLOTS_enhancer(opt)
  cummulative_PLOTS_ASE(opt)
  cummulative_PLOTS_E_Plus_ASE(opt)
  
  
  
}
  
  
  
 

###########################################################################

system.time( main() )
