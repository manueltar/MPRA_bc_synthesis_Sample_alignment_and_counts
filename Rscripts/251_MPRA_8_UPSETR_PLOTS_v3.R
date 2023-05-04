
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


UpSetR_like_E_Plus_ASE = function(option_list)
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
  
  #### Categories colors ----
  
  
  df_CSQ_colors<-readRDS(file=opt$CSQ_colors)
  
  df_CSQ_colors$color[df_CSQ_colors$VEP_DEF_LABELS == "TFBS"]<-"red"
  
  
  cat("df_CSQ_colors_0\n")
  cat(str(df_CSQ_colors))
  cat("\n")
  
  #### READ DATA----
  
  setwd(out2)
  
  filename<-paste('Element_collapse','.rds',sep='')

  KEY_collapse<-readRDS(file=filename)
  
  
  KEY_collapse<-KEY_collapse[order(KEY_collapse$KEY_Plus_carried_variants,KEY_collapse$Cell_Type),]
  
  cat("KEY_collapse_0\n")
  cat(str(KEY_collapse))
  cat("\n")
  
  
  filename<-paste('df_Cell_colors','.rds',sep='')
  
  df_Cell_colors<-readRDS(file=filename)
  
  cat("df_Cell_colors_0\n")
  cat(str(df_Cell_colors))
  cat("\n")
  
  ##### enhancer UpSetR_like_plots ----
  
  
  path5<-paste(out2,'UpSetR_like_plots','/', sep='')
  
  cat("path5\n")
  cat(sprintf(as.character(path5)))
  cat("\n")
  
  
  if (file.exists(path5)){
    
    
    
    
  } else {
    dir.create(file.path(path5))
    
  }
  
  
  #### E_Plus_ASE_CLASS ----
  
  
  path6<-paste(out2,'UpSetR_like_plots','/','E_Plus_ASE','/', sep='')
  
  cat("path6\n")
  cat(sprintf(as.character(path6)))
  cat("\n")
  
  
  if (file.exists(path6)){
    
    
    
    
  } else {
    dir.create(file.path(path6))
    
  }
  
  KEY_collapse_melt<-melt(KEY_collapse,id.vars=c("KEY_Plus_carried_variants","KEY","carried_variants","VAR","chr","Label","Label_2","VEP_DEF_LABELS","Cell_Type"))
  
  cat("KEY_collapse_melt\n")
  cat(str(KEY_collapse_melt))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(KEY_collapse_melt$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(KEY_collapse_melt$variable)))))
  cat("\n")
  
  
  
  KEY_collapse_melt_sel_E_Plus_ASE_CLASS<-KEY_collapse_melt[which(KEY_collapse_melt$variable == "E_Plus_ASE_CLASS"),]
  
  
  cat("KEY_collapse_melt_sel_E_Plus_ASE_CLASS\n")
  cat(str(KEY_collapse_melt_sel_E_Plus_ASE_CLASS))
  cat("\n")
  
  
  KEY_collapse_melt_sel_E_Plus_ASE_CLASS$value<-factor(as.character(KEY_collapse_melt_sel_E_Plus_ASE_CLASS$value),
                                      levels = c("0","1","2","3 or >"),
                                      ordered=T)
  
  
  cat("KEY_collapse_melt_sel_E_Plus_ASE_CLASS_2\n")
  cat(str(KEY_collapse_melt_sel_E_Plus_ASE_CLASS))
  cat("\n")
  
  
  value_arrays<-levels(KEY_collapse_melt_sel_E_Plus_ASE_CLASS$value)
  
  cat("value_arrays\n")
  cat(str(value_arrays))
  cat("\n")
  
  
  for(iteration_value_arrays in 1:length(value_arrays))
  {
    
    value_arrays_sel<-value_arrays[iteration_value_arrays]
    
    cat("------------------------------------------------------------------>\t")
    cat(sprintf(as.character(value_arrays_sel)))
    cat("\n")
    
    
    if(value_arrays_sel == "0")
    {
      
      accepted_values<-"0"
    }
    if(value_arrays_sel == "1")
    {
      
      accepted_values<-c("1","2","3 or >")
    }
    if(value_arrays_sel == "2")
    {
      
      accepted_values<-c("2","3 or >")
    }
    if(value_arrays_sel == "3 or >")
    {
      
      accepted_values<-c("3 or >")
    }
    
    
    
    
    df_sel<-KEY_collapse_melt_sel_E_Plus_ASE_CLASS[which(KEY_collapse_melt_sel_E_Plus_ASE_CLASS$value%in%accepted_values),]
    
    cat("df_sel\n")
    cat(str(df_sel))
    cat("\n")
      
    df_sel.dt<-data.table(df_sel,key=c("KEY_Plus_carried_variants","KEY","carried_variants","VAR","chr","Label","Label_2","VEP_DEF_LABELS"))
    
    cat("df_sel.dt\n")
    cat(str(df_sel.dt))
    cat("\n")
    
    
    
    df_sel_CT_summarised<-as.data.frame(df_sel.dt[,.(Cell_Type_string=paste(Cell_Type,collapse="|")), by=key(df_sel.dt)],stringsAsFactors=F)
    
    cat("df_sel_CT_summarised\n")
    cat(str(df_sel_CT_summarised))
    cat("\n")
    
    
    df_sel_CT_summarised$Classif_DEF<-"NA"
    
    df_sel_CT_summarised$Classif_DEF<-factor(df_sel_CT_summarised$Cell_Type_string,
                                                                 levels = c("K562|CHRF|HL60|THP1",
                                                                            "K562|CHRF|HL60","K562|CHRF|THP1","K562|HL60|THP1","CHRF|HL60|THP1",
                                                                            "K562|CHRF","K562|HL60","K562|THP1",
                                                                            "CHRF|THP1","CHRF|HL60",
                                                                            "HL60|THP1",
                                                                            "K562","CHRF","HL60","THP1"),
                                                                 ordered=T)
    
    cat(sprintf(as.character(names(summary(df_sel_CT_summarised$Classif_DEF)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_sel_CT_summarised$Classif_DEF))))
    cat("\n")
    
    df_sel_CT_summarised<-droplevels(df_sel_CT_summarised)


    cat("df_sel_CT_summarised_droplevels\n")
    cat(str(df_sel_CT_summarised))
    cat("\n")

    cat(sprintf(as.character(names(summary(df_sel_CT_summarised$Classif_DEF)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_sel_CT_summarised$Classif_DEF))))
    cat("\n")
    
    
    
    #### Freq table ----
    
    
    
    df_sel_CT_summarised_Classif_DEF.dt<-data.table(df_sel_CT_summarised, key=c("VEP_DEF_LABELS","Classif_DEF"))
    
    
    df_sel_CT_summarised_Classif_DEF_Fq<-as.data.frame(df_sel_CT_summarised_Classif_DEF.dt[,.N,by=key(df_sel_CT_summarised_Classif_DEF.dt)], stringsAsFactors=F)
    
    colnames(df_sel_CT_summarised_Classif_DEF_Fq)[which(colnames(df_sel_CT_summarised_Classif_DEF_Fq) == "N")]<-"instances"
    
    cat("df_sel_CT_summarised_Classif_DEF_Fq_0\n")
    cat(str(df_sel_CT_summarised_Classif_DEF_Fq))
    cat("\n")
    
    
    
    df_sel_CT_summarised_Classif_DEF_TOTAL.dt<-data.table(df_sel_CT_summarised, key=c("Classif_DEF"))
    
    
    df_sel_CT_summarised_Classif_DEF_TOTAL<-as.data.frame(df_sel_CT_summarised_Classif_DEF_TOTAL.dt[,.N,by=key(df_sel_CT_summarised_Classif_DEF_TOTAL.dt)], stringsAsFactors=F)
    
    colnames(df_sel_CT_summarised_Classif_DEF_TOTAL)[which(colnames(df_sel_CT_summarised_Classif_DEF_TOTAL) == "N")]<-"TOTAL"
    
    cat("df_sel_CT_summarised_Classif_DEF_TOTAL_\n")
    cat(str(df_sel_CT_summarised_Classif_DEF_TOTAL))
    cat("\n")
    
    df_sel_CT_summarised_Classif_DEF_Fq<-merge(df_sel_CT_summarised_Classif_DEF_Fq,
                                    df_sel_CT_summarised_Classif_DEF_TOTAL,
                                    by="Classif_DEF")
    
    df_sel_CT_summarised_Classif_DEF_Fq$Perc<-round(100*(df_sel_CT_summarised_Classif_DEF_Fq$instances/df_sel_CT_summarised_Classif_DEF_Fq$TOTAL),1)
    
    
    cat("df_sel_CT_summarised_Classif_DEF_Fq_1\n")
    cat(str(df_sel_CT_summarised_Classif_DEF_Fq))
    cat("\n")
    
  
    
    
    VEP_LABELS_df.dt<-data.table(df_sel_CT_summarised, key=c("VEP_DEF_LABELS"))
    
    
    VEP_LABELS_df<-as.data.frame(VEP_LABELS_df.dt[,.N,by=key(VEP_LABELS_df.dt)], stringsAsFactors=F)
    
    colnames(VEP_LABELS_df)[which(colnames(VEP_LABELS_df) == "N")]<-"TOTAL"
    
    cat("VEP_LABELS_df_\n")
    cat(str(VEP_LABELS_df))
    cat("\n")
    
    
    
    
    step<-round(max(df_sel_CT_summarised_Classif_DEF_Fq$TOTAL)/10,0)

    cat("--------->\t")
    cat(sprintf(as.character(step)))
    cat("\n")

    if(step == 0){

      step=1
    }


    breaks.y<-seq(0,max(df_sel_CT_summarised_Classif_DEF_Fq$TOTAL)+step, by=step)
    labels.y<-as.character(breaks.y)


    cat(sprintf(as.character(breaks.y)))
    cat("\n")
    
    ### df_CSQ_colors

    
    df_CSQ_colors_sel<-df_CSQ_colors[which(df_CSQ_colors$VEP_DEF_LABELS%in%df_sel_CT_summarised_Classif_DEF_Fq$VEP_DEF_LABELS),]
    
    cat("df_CSQ_colors_sel_1\n")
    cat(str(df_CSQ_colors_sel))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_CSQ_colors_sel$VEP_DEF_LABELS)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_CSQ_colors_sel$VEP_DEF_LABELS))))
    cat("\n")
    
    missing_CTRL_levels<-levels(df_sel_CT_summarised_Classif_DEF_Fq$VEP_DEF_LABELS)[-which(levels(df_sel_CT_summarised_Classif_DEF_Fq$VEP_DEF_LABELS)%in%levels(df_CSQ_colors$VEP_DEF_LABELS))]
    
    cat("missing_CTRL_levels_1\n")
    cat(str(missing_CTRL_levels))
    cat("\n")
    
    if(length(missing_CTRL_levels) >0)
    {
      color_vector<-c("black",'#1877C9','gray')
      
      A.df<-as.data.frame(cbind(c("NCGR","ASE_CTRL","enhancer_CTRL"),color_vector), stringsAsFactors=F)
      
      colnames(A.df)<-c("VEP_DEF_LABELS","color")
      
      cat("A.df_1\n")
      cat(str(A.df))
      cat("\n")
      
      A.df_sel<-A.df[which(A.df$VEP_DEF_LABELS%in%missing_CTRL_levels),]
      
      cat("A.df_sel_1\n")
      cat(str(A.df_sel))
      cat("\n")
      
      
      df_CSQ_colors_sel<-rbind(df_CSQ_colors_sel,A.df_sel)
      
      
      df_CSQ_colors_sel$VEP_DEF_LABELS<-factor(df_CSQ_colors_sel$VEP_DEF_LABELS,
                                            levels=c("LOF","MISS","SYN","UTR5","UTR3",
                                                     "INTRON","INTERGENIC","UPSTREAM","DOWNSTREAM","REGULATORY",
                                                     "TFBS","SPLICE",
                                                     "OTHER","NMD","NCT","PCHIC_Relevant_link",
                                                     "NCGR","ASE_CTRL","enhancer_CTRL","Kousik_variant"),
                                            ordered=T)
      
      
      df_CSQ_colors_sel<-droplevels(df_CSQ_colors_sel)
      
      
      
    }
    
    cat("df_CSQ_colors_sel_2\n")
    cat(str(df_CSQ_colors_sel))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_CSQ_colors_sel$VEP_DEF_LABELS)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_CSQ_colors_sel$VEP_DEF_LABELS))))
    cat("\n")
    
    df_CSQ_colors_sel<-merge(VEP_LABELS_df,
                             df_CSQ_colors_sel,
                             by="VEP_DEF_LABELS")
    
    df_CSQ_colors_sel<-df_CSQ_colors_sel[order(df_CSQ_colors_sel$VEP_DEF_LABELS),]
    
    cat("df_CSQ_colors_sel_3\n")
    cat(str(df_CSQ_colors_sel))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_CSQ_colors_sel$VEP_DEF_LABELS)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_CSQ_colors_sel$VEP_DEF_LABELS))))
    cat("\n")
    
    ### graph
   
    # ,
    # fill=VEP_DEF_LABELS
    # scale_fill_manual(values=df_CSQ_colors_sel$color, drop=F,
    #                   name="Totals", breaks=df_CSQ_colors_sel$VEP_DEF_LABELS,
    #                   labels=paste(df_CSQ_colors_sel$VEP_DEF_LABELS,
    #                                df_CSQ_colors_sel$TOTAL, sep =' n= '))+
    
    
    cat("---------------------------->df_sel_CT_summarised_Classif_DEF_TOTAL\n")
    cat(str(df_sel_CT_summarised_Classif_DEF_TOTAL))
    cat("\n")

    graph<-ggplot(data=df_sel_CT_summarised_Classif_DEF_TOTAL,
                  aes(x=Classif_DEF,
                      y=TOTAL)) +
      geom_bar(stat="identity",colour='black')+
      theme_bw()+
      theme(axis.title.y=element_text(size=18, family="sans"),
            axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
            axis.text.x=element_text(angle=90,vjust=1,hjust=1,size=8, color="black", family="sans"))+
      scale_y_continuous(name=paste("Variants with at least",value_arrays_sel,"ACTIVE TILES",sep=" "),breaks=breaks.y,labels=labels.y,
                         limits=c(breaks.y[1],breaks.y[length(breaks.y)]))+
      scale_x_discrete(name=NULL, drop=F)+
      theme(legend.position="right")+
      theme_classic()+
      ggeasy::easy_center_title()
    
    ### subgraph
    
    Cell_Type.dt<-data.table(df_sel,key=c("Cell_Type"))
    
    cat("Cell_Type.dt\n")
    cat(str(Cell_Type.dt))
    cat("\n")
    
    
    
    Cell_Type_df<-as.data.frame(Cell_Type.dt[,.N, by=key(Cell_Type.dt)],stringsAsFactors=F)
    colnames(Cell_Type_df)[which(colnames(Cell_Type_df) == "N")]<-"Total"
   
    
    Cell_Type_df<-merge(Cell_Type_df,
                        df_Cell_colors,
                        by="Cell_Type",
                        all=T)
    

    cat("Cell_Type_df_1\n")
    cat(str(Cell_Type_df))
    cat("\n")
    
    Cell_Type_df$Total[is.na(Cell_Type_df$Total)]<-0
    
    
    cat("Cell_Type_df_WITHOUT_CTRLS0\n")
    cat(str(Cell_Type_df))
    cat("\n")
   


    step<-round(max(Cell_Type_df$Total)/5,0)

    cat("--------->\t")
    cat(sprintf(as.character(step)))
    cat("\n")

    if(step == 0){

      step=1
    }
    breaks.x<-rev(c(seq(0,max(Cell_Type_df$Total)+step, by=step)))
    labels.x<-as.character(breaks.x)


    cat(sprintf(as.character(breaks.x)))
    cat("\n")
    
    levels_CT<-rev(levels(Cell_Type_df$Cell_Type))
    
    cat("levels_CT_2\n")
    cat(str(levels_CT))
    cat("\n")
    
    Cell_Type_df$Cell_type<-factor(Cell_Type_df$Cell_Type,
                                   levels=levels_CT,
                                   ordered=T)
    
    # Cell_Type_df<-Cell_Type_df[order(Cell_Type_df$Cell_Type),]
    
    cat("Cell_Type_df_2\n")
    cat(str(Cell_Type_df))
    cat("\n")
    
    
    subgraph<-ggplot(data=Cell_Type_df,
                     aes(y=Cell_Type,
                         x=Total,
                         fill=Cell_Type)) +
      geom_bar(stat="identity",colour='black')+
      theme_bw()+
      theme(axis.title.y=element_text(size=16, family="sans"),
            axis.text.y=element_text(angle=0,size=16, color="black", family="sans"),
            axis.text.x=element_text(angle=0,size=10, color="black", family="sans"))+
      scale_x_reverse(name=paste("Variants with at least",value_arrays_sel,"ACTIVE TILES",sep=" "),breaks=breaks.x,labels=labels.x,
                      limits=c(breaks.x[1],breaks.x[length(breaks.x)]))+
      scale_y_discrete(position="right",name=NULL, drop=F)+
      scale_fill_manual(values=Cell_Type_df$colors, drop=F)+
      theme(legend.position="hidden")+
      theme_classic()+
      ggeasy::easy_center_title()
    
    cat("subgraph DONE\n")
    
    
  


    graph_FINAL<-plot_grid(NULL,graph,subgraph,NULL,
                           nrow = 2,
                           ncol=2,
                           rel_heights = c(1, 0.25),
                           rel_widths=c(0.5,1))
   
    path7<-paste(out2,'UpSetR_like_plots','/','E_Plus_ASE','/','WITH_CTRLS','/', sep='')
    
    cat("path7\n")
    cat(sprintf(as.character(path7)))
    cat("\n")
    
    
    if (file.exists(path7)){
      
      
      
      
    } else {
      dir.create(file.path(path7))
      
    }
    
    
    setwd(path7)
    
    svglite(paste('ACTIVE_TILES_AT_LEAST_',value_arrays_sel,'.svg',sep=''), width = 8, height = 8)
    print(graph_FINAL)
    dev.off()
    
    
    
    # ###################################################################################################
    # quit(status = 1)
    
    
    #### WITHOUT CTRLS ----
    
    
    df_sel<-KEY_collapse_melt_sel_E_Plus_ASE_CLASS[which(KEY_collapse_melt_sel_E_Plus_ASE_CLASS$value%in%accepted_values),]
    
    cat("df_sel_0\n")
    cat(str(df_sel))
    cat("\n")
    
    df_sel<-df_sel[which(df_sel$Label_2 == "ASSAYED_VARIANT"),]
    
    cat("df_sel_1\n")
    cat(str(df_sel))
    cat("\n")
    
    df_sel.dt<-data.table(df_sel,key=c("KEY_Plus_carried_variants","KEY","carried_variants","VAR","chr","Label","Label_2","VEP_DEF_LABELS"))
    
    cat("df_sel.dt\n")
    cat(str(df_sel.dt))
    cat("\n")
    
    
    
    df_sel_CT_summarised<-as.data.frame(df_sel.dt[,.(Cell_Type_string=paste(Cell_Type,collapse="|")), by=key(df_sel.dt)],stringsAsFactors=F)
    
    cat("df_sel_CT_summarised\n")
    cat(str(df_sel_CT_summarised))
    cat("\n")
    
    
    df_sel_CT_summarised$Classif_DEF<-"NA"
    
    df_sel_CT_summarised$Classif_DEF<-factor(df_sel_CT_summarised$Cell_Type_string,
                                             levels = c("K562|CHRF|HL60|THP1",
                                                        "K562|CHRF|HL60","K562|CHRF|THP1","K562|HL60|THP1","CHRF|HL60|THP1",
                                                        "K562|CHRF","K562|HL60","K562|THP1",
                                                        "CHRF|THP1","CHRF|HL60",
                                                        "HL60|THP1",
                                                        "K562","CHRF","HL60","THP1"),
                                             ordered=T)
    
    cat(sprintf(as.character(names(summary(df_sel_CT_summarised$Classif_DEF)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_sel_CT_summarised$Classif_DEF))))
    cat("\n")
    
    df_sel_CT_summarised<-droplevels(df_sel_CT_summarised)
    
    
    cat("df_sel_CT_summarised_droplevels\n")
    cat(str(df_sel_CT_summarised))
    cat("\n")
    
    cat(sprintf(as.character(names(summary(df_sel_CT_summarised$Classif_DEF)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_sel_CT_summarised$Classif_DEF))))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_sel_CT_summarised$VEP_DEF_LABELS)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_sel_CT_summarised$VEP_DEF_LABELS))))
    cat("\n")
    
    ####
    
    df_CSQ_colors_sel<-df_CSQ_colors_sel[which(df_CSQ_colors_sel$VEP_DEF_LABELS%in%df_sel_CT_summarised$VEP_DEF_LABELS),]
    df_CSQ_colors_sel<-droplevels(df_CSQ_colors_sel)
    
    
    cat("df_CSQ_colors_sel\n")
    cat(str(df_CSQ_colors_sel))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_CSQ_colors_sel$VEP_DEF_LABELS)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_CSQ_colors_sel$VEP_DEF_LABELS))))
    cat("\n")
      
      
    df_sel<-droplevels(df_sel)
    
    
    cat("df_sel\n")
    cat(str(df_sel))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_sel$VEP_DEF_LABELS)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_sel$VEP_DEF_LABELS))))
    cat("\n")
    
    
    
    #### Freq table ----
    
    
    
    df_sel_CT_summarised_Classif_DEF.dt<-data.table(df_sel_CT_summarised, key=c("VEP_DEF_LABELS","Classif_DEF"))
    
    
    df_sel_CT_summarised_Classif_DEF_Fq<-as.data.frame(df_sel_CT_summarised_Classif_DEF.dt[,.N,by=key(df_sel_CT_summarised_Classif_DEF.dt)], stringsAsFactors=F)
    
    colnames(df_sel_CT_summarised_Classif_DEF_Fq)[which(colnames(df_sel_CT_summarised_Classif_DEF_Fq) == "N")]<-"instances"
    
    cat("df_sel_CT_summarised_Classif_DEF_Fq_0\n")
    cat(str(df_sel_CT_summarised_Classif_DEF_Fq))
    cat("\n")
    
    
    
    df_sel_CT_summarised_Classif_DEF_TOTAL.dt<-data.table(df_sel_CT_summarised, key=c("Classif_DEF"))
    
    
    df_sel_CT_summarised_Classif_DEF_TOTAL<-as.data.frame(df_sel_CT_summarised_Classif_DEF_TOTAL.dt[,.N,by=key(df_sel_CT_summarised_Classif_DEF_TOTAL.dt)], stringsAsFactors=F)
    
    colnames(df_sel_CT_summarised_Classif_DEF_TOTAL)[which(colnames(df_sel_CT_summarised_Classif_DEF_TOTAL) == "N")]<-"TOTAL"
    
    cat("df_sel_CT_summarised_Classif_DEF_TOTAL_\n")
    cat(str(df_sel_CT_summarised_Classif_DEF_TOTAL))
    cat("\n")
    
    df_sel_CT_summarised_Classif_DEF_Fq<-merge(df_sel_CT_summarised_Classif_DEF_Fq,
                                               df_sel_CT_summarised_Classif_DEF_TOTAL,
                                               by="Classif_DEF")
    
    df_sel_CT_summarised_Classif_DEF_Fq$Perc<-round(100*(df_sel_CT_summarised_Classif_DEF_Fq$instances/df_sel_CT_summarised_Classif_DEF_Fq$TOTAL),1)
    
    
    cat("df_sel_CT_summarised_Classif_DEF_Fq_1\n")
    cat(str(df_sel_CT_summarised_Classif_DEF_Fq))
    cat("\n")
    
    
    
    
    VEP_LABELS_df.dt<-data.table(df_sel_CT_summarised, key=c("VEP_DEF_LABELS"))
    
    
    VEP_LABELS_df<-as.data.frame(VEP_LABELS_df.dt[,.N,by=key(VEP_LABELS_df.dt)], stringsAsFactors=F)
    
    colnames(VEP_LABELS_df)[which(colnames(VEP_LABELS_df) == "N")]<-"TOTAL"
    
    cat("VEP_LABELS_df_\n")
    cat(str(VEP_LABELS_df))
    cat("\n")
    
    
    
    
    step<-round(max(df_sel_CT_summarised_Classif_DEF_Fq$TOTAL)/10,0)
    
    cat("--------->\t")
    cat(sprintf(as.character(step)))
    cat("\n")
    
    if(step == 0){
      
      step=1
    }
    
    
    breaks.y<-seq(0,max(df_sel_CT_summarised_Classif_DEF_Fq$TOTAL)+step, by=step)
    labels.y<-as.character(breaks.y)
    
    
    cat(sprintf(as.character(breaks.y)))
    cat("\n")
    
    ### graph
    
    # fill=VEP_DEF_LABELS
    # scale_fill_manual(values=df_CSQ_colors_sel$color, drop=F,
    #                   name="Totals", breaks=df_CSQ_colors_sel$VEP_DEF_LABELS,
    #                   labels=paste(df_CSQ_colors_sel$VEP_DEF_LABELS,
    #                                df_CSQ_colors_sel$TOTAL, sep =' n= '))+
    
    graph<-ggplot(data=df_sel_CT_summarised_Classif_DEF_TOTAL,
                  aes(x=Classif_DEF,
                      y=TOTAL)) +
      geom_bar(stat="identity",colour='black')+
      theme_bw()+
      theme_classic()+
      theme(axis.title.y=element_text(size=18, family="sans"),
            axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
            axis.text.x=element_text(angle=90,vjust=1,hjust=1,size=8, color="black", family="sans"))+
      scale_y_continuous(name=paste("Variants with at least",value_arrays_sel,"ACTIVE TILES",sep=" "),breaks=breaks.y,labels=labels.y,
                         limits=c(breaks.y[1],breaks.y[length(breaks.y)]))+
      scale_x_discrete(name=NULL, drop=F)+
      theme(legend.position="hidden")+
      ggeasy::easy_center_title()
    
    ### subgraph
    
    Cell_Type.dt<-data.table(df_sel,key=c("Cell_Type"))
    
    cat("Cell_Type.dt\n")
    cat(str(Cell_Type.dt))
    cat("\n")
    
    
    
    Cell_Type_df<-as.data.frame(Cell_Type.dt[,.N, by=key(Cell_Type.dt)],stringsAsFactors=F)
    colnames(Cell_Type_df)[which(colnames(Cell_Type_df) == "N")]<-"Total"
   
    
    
    cat("Cell_Type_df_0\n")
    cat(str(Cell_Type_df))
    cat("\n")
    
    Cell_Type_df<-merge(Cell_Type_df,
                        df_Cell_colors,
                        by="Cell_Type",
                        all=T)
    
    Cell_Type_df$Total[is.na(Cell_Type_df$Total)]<-0
    
    
    cat("Cell_Type_df_1\n")
    cat(str(Cell_Type_df))
    cat("\n")
    
   
    step<-round(max(Cell_Type_df$Total)/5,0)
    
    cat("--------->\t")
    cat(sprintf(as.character(step)))
    cat("\n")
    
    if(step == 0){
      
      step=1
    }
    breaks.x<-rev(c(seq(0,max(Cell_Type_df$Total)+step, by=step)))
    labels.x<-as.character(breaks.x)
    
    
    cat(sprintf(as.character(breaks.x)))
    cat("\n")
    
    levels_CT<-rev(levels(Cell_Type_df$Cell_Type))
    
    cat("levels_CT_2\n")
    cat(str(levels_CT))
    cat("\n")
    
    Cell_Type_df$Cell_type<-factor(Cell_Type_df$Cell_Type,
                                   levels=levels_CT,
                                   ordered=T)
    
    # Cell_Type_df<-Cell_Type_df[order(Cell_Type_df$Cell_Type),]
    
    
    
    cat("Cell_Type_df_2\n")
    cat(str(Cell_Type_df))
    cat("\n")
    
    
    subgraph<-ggplot(data=Cell_Type_df,
                     aes(y=Cell_Type,
                         x=Total,
                         fill=Cell_Type)) +
      geom_bar(stat="identity",colour='black')+
      theme_bw()+
      theme_classic()+
      theme(axis.title.y=element_text(size=16, family="sans"),
            axis.text.y=element_text(angle=0,size=16, color="black", family="sans"),
            axis.text.x=element_text(angle=0,size=10, color="black", family="sans"))+
      scale_x_reverse(name=paste("Variants with at least",value_arrays_sel,"ACTIVE TILES",sep=" "),breaks=breaks.x,labels=labels.x,
                      limits=c(breaks.x[1],breaks.x[length(breaks.x)]))+
      scale_y_discrete(position="right",name=NULL, drop=F)+
      scale_fill_manual(values=Cell_Type_df$colors, drop=F)+
      theme(legend.position="hidden")+
      ggeasy::easy_center_title()
    
    cat("subgraph DONE WITHOUT_CTRLS\n")
    
    
    
    
    
    graph_FINAL<-plot_grid(NULL,graph,subgraph,NULL,
                           nrow = 2,
                           ncol=2,
                           rel_heights = c(1, 0.25),
                           rel_widths=c(0.5,1))
    
    path8<-paste(out2,'UpSetR_like_plots','/','E_Plus_ASE','/','WITHOUT_CTRLS','/', sep='')
    
    cat("path8\n")
    cat(sprintf(as.character(path8)))
    cat("\n")
    
    
    if (file.exists(path8)){
      
      
      
      
    } else {
      dir.create(file.path(path8))
      
    }
    
    
    setwd(path8)
    
    svglite(paste('ACTIVE_TILES_AT_LEAST_',value_arrays_sel,'.svg',sep=''), width = 8, height = 8)
    print(graph_FINAL)
    dev.off()
    
    
    
    
    
    

  }
  
 
  
  
  
  
 
}

UpSetR_like_enhancer = function(option_list)
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
  
  #### Categories colors ----
  
  
  df_CSQ_colors<-readRDS(file=opt$CSQ_colors)
  
  df_CSQ_colors$color[df_CSQ_colors$VEP_DEF_LABELS == "TFBS"]<-"red"
  
  
  cat("df_CSQ_colors_0\n")
  cat(str(df_CSQ_colors))
  cat("\n")
  
  #### READ DATA----
  
  setwd(out2)
  
  
  
  
  filename<-paste('Element_collapse','.rds',sep='')
  
  
  KEY_collapse<-readRDS(file=filename)
  
  
  KEY_collapse<-KEY_collapse[order(KEY_collapse$KEY_Plus_carried_variants,KEY_collapse$Cell_Type),]
  
  cat("KEY_collapse_0\n")
  cat(str(KEY_collapse))
  cat("\n")
  
  
  filename<-paste('df_Cell_colors','.rds',sep='')
  
  df_Cell_colors<-readRDS(file=filename)
  
  cat("df_Cell_colors_0\n")
  cat(str(df_Cell_colors))
  cat("\n")
  
  ##### enhancer UpSetR_like_plots ----
  
  
  path5<-paste(out2,'UpSetR_like_plots','/', sep='')
  
  cat("path5\n")
  cat(sprintf(as.character(path5)))
  cat("\n")
  
  
  if (file.exists(path5)){
    
    
    
    
  } else {
    dir.create(file.path(path5))
    
  }
  
  
  #### enhancer_CLASS ----
  
  
  path6<-paste(out2,'UpSetR_like_plots','/','enhancer','/', sep='')
  
  cat("path6\n")
  cat(sprintf(as.character(path6)))
  cat("\n")
  
  
  if (file.exists(path6)){
    
    
    
    
  } else {
    dir.create(file.path(path6))
    
  }
  
  KEY_collapse_melt<-melt(KEY_collapse,id.vars=c("KEY_Plus_carried_variants","KEY","carried_variants","VAR","chr","Label","Label_2","VEP_DEF_LABELS","Cell_Type"))
  
  cat("KEY_collapse_melt\n")
  cat(str(KEY_collapse_melt))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(KEY_collapse_melt$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(KEY_collapse_melt$variable)))))
  cat("\n")
  
  
  
  KEY_collapse_melt_sel_enhancer_CLASS<-KEY_collapse_melt[which(KEY_collapse_melt$variable == "enhancer_CLASS"),]
  
  
  cat("KEY_collapse_melt_sel_enhancer_CLASS\n")
  cat(str(KEY_collapse_melt_sel_enhancer_CLASS))
  cat("\n")
  
  
  KEY_collapse_melt_sel_enhancer_CLASS$value<-factor(as.character(KEY_collapse_melt_sel_enhancer_CLASS$value),
                                                       levels = c("0","1","2","3 or >"),
                                                       ordered=T)
  
  
  cat("KEY_collapse_melt_sel_enhancer_CLASS_2\n")
  cat(str(KEY_collapse_melt_sel_enhancer_CLASS))
  cat("\n")
  
  
  value_arrays<-levels(KEY_collapse_melt_sel_enhancer_CLASS$value)
  
  cat("value_arrays\n")
  cat(str(value_arrays))
  cat("\n")
  
  
  for(iteration_value_arrays in 1:length(value_arrays))
  {
    
    value_arrays_sel<-value_arrays[iteration_value_arrays]
    
    cat("------------------------------------------------------------------>\t")
    cat(sprintf(as.character(value_arrays_sel)))
    cat("\n")
    
    
    if(value_arrays_sel == "0")
    {
      
      accepted_values<-"0"
    }
    if(value_arrays_sel == "1")
    {
      
      accepted_values<-c("1","2","3 or >")
    }
    if(value_arrays_sel == "2")
    {
      
      accepted_values<-c("2","3 or >")
    }
    if(value_arrays_sel == "3 or >")
    {
      
      accepted_values<-c("3 or >")
    }
    
    
    
    
    df_sel<-KEY_collapse_melt_sel_enhancer_CLASS[which(KEY_collapse_melt_sel_enhancer_CLASS$value%in%accepted_values),]
    
    cat("df_sel\n")
    cat(str(df_sel))
    cat("\n")
    
    df_sel.dt<-data.table(df_sel,key=c("KEY_Plus_carried_variants","KEY","carried_variants","VAR","chr","Label","Label_2","VEP_DEF_LABELS"))
    
    cat("df_sel.dt\n")
    cat(str(df_sel.dt))
    cat("\n")
    
    
    
    df_sel_CT_summarised<-as.data.frame(df_sel.dt[,.(Cell_Type_string=paste(Cell_Type,collapse="|")), by=key(df_sel.dt)],stringsAsFactors=F)
    
    cat("df_sel_CT_summarised\n")
    cat(str(df_sel_CT_summarised))
    cat("\n")
    
    
    df_sel_CT_summarised$Classif_DEF<-"NA"
    
    df_sel_CT_summarised$Classif_DEF<-factor(df_sel_CT_summarised$Cell_Type_string,
                                             levels = c("K562|CHRF|HL60|THP1",
                                                        "K562|CHRF|HL60","K562|CHRF|THP1","K562|HL60|THP1","CHRF|HL60|THP1",
                                                        "K562|CHRF","K562|HL60","K562|THP1",
                                                        "CHRF|THP1","CHRF|HL60",
                                                        "HL60|THP1",
                                                        "K562","CHRF","HL60","THP1"),
                                             ordered=T)
    
    cat(sprintf(as.character(names(summary(df_sel_CT_summarised$Classif_DEF)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_sel_CT_summarised$Classif_DEF))))
    cat("\n")
    
    df_sel_CT_summarised<-droplevels(df_sel_CT_summarised)
    
    
    cat("df_sel_CT_summarised_droplevels\n")
    cat(str(df_sel_CT_summarised))
    cat("\n")
    
    cat(sprintf(as.character(names(summary(df_sel_CT_summarised$Classif_DEF)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_sel_CT_summarised$Classif_DEF))))
    cat("\n")
    
    
    
    #### Freq table ----
    
    
    
    df_sel_CT_summarised_Classif_DEF.dt<-data.table(df_sel_CT_summarised, key=c("VEP_DEF_LABELS","Classif_DEF"))
    
    
    df_sel_CT_summarised_Classif_DEF_Fq<-as.data.frame(df_sel_CT_summarised_Classif_DEF.dt[,.N,by=key(df_sel_CT_summarised_Classif_DEF.dt)], stringsAsFactors=F)
    
    colnames(df_sel_CT_summarised_Classif_DEF_Fq)[which(colnames(df_sel_CT_summarised_Classif_DEF_Fq) == "N")]<-"instances"
    
    cat("df_sel_CT_summarised_Classif_DEF_Fq_0\n")
    cat(str(df_sel_CT_summarised_Classif_DEF_Fq))
    cat("\n")
    
    
    
    df_sel_CT_summarised_Classif_DEF_TOTAL.dt<-data.table(df_sel_CT_summarised, key=c("Classif_DEF"))
    
    
    df_sel_CT_summarised_Classif_DEF_TOTAL<-as.data.frame(df_sel_CT_summarised_Classif_DEF_TOTAL.dt[,.N,by=key(df_sel_CT_summarised_Classif_DEF_TOTAL.dt)], stringsAsFactors=F)
    
    colnames(df_sel_CT_summarised_Classif_DEF_TOTAL)[which(colnames(df_sel_CT_summarised_Classif_DEF_TOTAL) == "N")]<-"TOTAL"
    
    cat("df_sel_CT_summarised_Classif_DEF_TOTAL_\n")
    cat(str(df_sel_CT_summarised_Classif_DEF_TOTAL))
    cat("\n")
    
    df_sel_CT_summarised_Classif_DEF_Fq<-merge(df_sel_CT_summarised_Classif_DEF_Fq,
                                               df_sel_CT_summarised_Classif_DEF_TOTAL,
                                               by="Classif_DEF")
    
    df_sel_CT_summarised_Classif_DEF_Fq$Perc<-round(100*(df_sel_CT_summarised_Classif_DEF_Fq$instances/df_sel_CT_summarised_Classif_DEF_Fq$TOTAL),1)
    
    
    cat("df_sel_CT_summarised_Classif_DEF_Fq_1\n")
    cat(str(df_sel_CT_summarised_Classif_DEF_Fq))
    cat("\n")
    
    
    
    
    VEP_LABELS_df.dt<-data.table(df_sel_CT_summarised, key=c("VEP_DEF_LABELS"))
    
    
    VEP_LABELS_df<-as.data.frame(VEP_LABELS_df.dt[,.N,by=key(VEP_LABELS_df.dt)], stringsAsFactors=F)
    
    colnames(VEP_LABELS_df)[which(colnames(VEP_LABELS_df) == "N")]<-"TOTAL"
    
    cat("VEP_LABELS_df_\n")
    cat(str(VEP_LABELS_df))
    cat("\n")
    
    
    
    
    step<-round(max(df_sel_CT_summarised_Classif_DEF_Fq$TOTAL)/10,0)
    
    cat("--------->\t")
    cat(sprintf(as.character(step)))
    cat("\n")
    
    if(step == 0){
      
      step=1
    }
    
    
    breaks.y<-seq(0,max(df_sel_CT_summarised_Classif_DEF_Fq$TOTAL)+step, by=step)
    labels.y<-as.character(breaks.y)
    
    
    cat(sprintf(as.character(breaks.y)))
    cat("\n")
    
    ### df_CSQ_colors
    
    
    df_CSQ_colors_sel<-df_CSQ_colors[which(df_CSQ_colors$VEP_DEF_LABELS%in%df_sel_CT_summarised_Classif_DEF_Fq$VEP_DEF_LABELS),]
    
    cat("df_CSQ_colors_sel_1\n")
    cat(str(df_CSQ_colors_sel))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_CSQ_colors_sel$VEP_DEF_LABELS)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_CSQ_colors_sel$VEP_DEF_LABELS))))
    cat("\n")
    
    missing_CTRL_levels<-levels(df_sel_CT_summarised_Classif_DEF_Fq$VEP_DEF_LABELS)[-which(levels(df_sel_CT_summarised_Classif_DEF_Fq$VEP_DEF_LABELS)%in%levels(df_CSQ_colors$VEP_DEF_LABELS))]
    
    cat("missing_CTRL_levels_1\n")
    cat(str(missing_CTRL_levels))
    cat("\n")
    
    if(length(missing_CTRL_levels) >0)
    {
      color_vector<-c("black",'#1877C9','gray')
      
      A.df<-as.data.frame(cbind(c("NCGR","ASE_CTRL","enhancer_CTRL"),color_vector), stringsAsFactors=F)
      
      colnames(A.df)<-c("VEP_DEF_LABELS","color")
      
      cat("A.df_1\n")
      cat(str(A.df))
      cat("\n")
      
      A.df_sel<-A.df[which(A.df$VEP_DEF_LABELS%in%missing_CTRL_levels),]
      
      cat("A.df_sel_1\n")
      cat(str(A.df_sel))
      cat("\n")
      
      
      df_CSQ_colors_sel<-rbind(df_CSQ_colors_sel,A.df_sel)
      
      
      df_CSQ_colors_sel$VEP_DEF_LABELS<-factor(df_CSQ_colors_sel$VEP_DEF_LABELS,
                                               levels=c("LOF","MISS","SYN","UTR5","UTR3",
                                                        "INTRON","INTERGENIC","UPSTREAM","DOWNSTREAM","REGULATORY",
                                                        "TFBS","SPLICE",
                                                        "OTHER","NMD","NCT","PCHIC_Relevant_link",
                                                        "NCGR","ASE_CTRL","enhancer_CTRL","Kousik_variant"),
                                               ordered=T)
      
      
      df_CSQ_colors_sel<-droplevels(df_CSQ_colors_sel)
      
      
      
    }
    
    cat("df_CSQ_colors_sel_2\n")
    cat(str(df_CSQ_colors_sel))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_CSQ_colors_sel$VEP_DEF_LABELS)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_CSQ_colors_sel$VEP_DEF_LABELS))))
    cat("\n")
    
    df_CSQ_colors_sel<-merge(VEP_LABELS_df,
                             df_CSQ_colors_sel,
                             by="VEP_DEF_LABELS")
    
    df_CSQ_colors_sel<-df_CSQ_colors_sel[order(df_CSQ_colors_sel$VEP_DEF_LABELS),]
    
    cat("df_CSQ_colors_sel_3\n")
    cat(str(df_CSQ_colors_sel))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_CSQ_colors_sel$VEP_DEF_LABELS)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_CSQ_colors_sel$VEP_DEF_LABELS))))
    cat("\n")
    
    ### graph
    
    # ,
    # fill=VEP_DEF_LABELS
    # scale_fill_manual(values=df_CSQ_colors_sel$color, drop=F,
    #                   name="Totals", breaks=df_CSQ_colors_sel$VEP_DEF_LABELS,
    #                   labels=paste(df_CSQ_colors_sel$VEP_DEF_LABELS,
    #                                df_CSQ_colors_sel$TOTAL, sep =' n= '))+
    
    
    graph<-ggplot(data=df_sel_CT_summarised_Classif_DEF_TOTAL,
                  aes(x=Classif_DEF,
                      y=TOTAL)) +
      geom_bar(stat="identity",colour='black')+
      theme_bw()+
      theme_classic()+
      theme(axis.title.y=element_text(size=18, family="sans"),
            axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
            axis.text.x=element_text(angle=90,vjust=1,hjust=1,size=8, color="black", family="sans"))+
      scale_y_continuous(name=paste("Variants with at least",value_arrays_sel,"ACTIVE TILES",sep=" "),breaks=breaks.y,labels=labels.y,
                         limits=c(breaks.y[1],breaks.y[length(breaks.y)]))+
      scale_x_discrete(name=NULL, drop=F)+
      theme(legend.position="hidden")+
      ggeasy::easy_center_title()
    
    ### subgraph
    
    Cell_Type.dt<-data.table(df_sel,key=c("Cell_Type"))
    
    cat("Cell_Type.dt\n")
    cat(str(Cell_Type.dt))
    cat("\n")
    
    
    
    Cell_Type_df<-as.data.frame(Cell_Type.dt[,.N, by=key(Cell_Type.dt)],stringsAsFactors=F)
    colnames(Cell_Type_df)[which(colnames(Cell_Type_df) == "N")]<-"Total"
    
    
    Cell_Type_df<-merge(Cell_Type_df,
                        df_Cell_colors,
                        by="Cell_Type",
                        all=T)
    
    
    cat("Cell_Type_df_1\n")
    cat(str(Cell_Type_df))
    cat("\n")
    
    Cell_Type_df$Total[is.na(Cell_Type_df$Total)]<-0
    
    
    cat("Cell_Type_df_WITHOUT_CTRLS0\n")
    cat(str(Cell_Type_df))
    cat("\n")
    
    
    
    step<-round(max(Cell_Type_df$Total)/5,0)
    
    cat("--------->\t")
    cat(sprintf(as.character(step)))
    cat("\n")
    
    if(step == 0){
      
      step=1
    }
    breaks.x<-rev(c(seq(0,max(Cell_Type_df$Total)+step, by=step)))
    labels.x<-as.character(breaks.x)
    
    
    cat(sprintf(as.character(breaks.x)))
    cat("\n")
    
    levels_CT<-rev(levels(Cell_Type_df$Cell_Type))
    
    cat("levels_CT_2\n")
    cat(str(levels_CT))
    cat("\n")
    
    Cell_Type_df$Cell_type<-factor(Cell_Type_df$Cell_Type,
                                   levels=levels_CT,
                                   ordered=T)
    
    # Cell_Type_df<-Cell_Type_df[order(Cell_Type_df$Cell_Type),]
    
    cat("Cell_Type_df_2\n")
    cat(str(Cell_Type_df))
    cat("\n")
    
    # ,
    # fill=Cell_Type
    # scale_fill_manual(values=Cell_Type_df$colors, drop=F)+
      
    subgraph<-ggplot(data=Cell_Type_df,
                     aes(y=Cell_Type,
                         x=Total)) +
      geom_bar(stat="identity",colour='black')+
      theme_bw()+
      theme_classic()+
      theme(axis.title.y=element_text(size=16, family="sans"),
            axis.text.y=element_text(angle=0,size=16, color="black", family="sans"),
            axis.text.x=element_text(angle=0,size=10, color="black", family="sans"))+
      scale_x_reverse(name=paste("Variants with at least",value_arrays_sel,"ACTIVE TILES",sep=" "),breaks=breaks.x,labels=labels.x,
                      limits=c(breaks.x[1],breaks.x[length(breaks.x)]))+
      scale_y_discrete(position="right",name=NULL, drop=F)+
      theme(legend.position="hidden")+
      ggeasy::easy_center_title()
    
    cat("subgraph DONE\n")
    
    
    
    
    
    graph_FINAL<-plot_grid(NULL,graph,subgraph,NULL,
                           nrow = 2,
                           ncol=2,
                           rel_heights = c(1, 0.25),
                           rel_widths=c(0.5,1))
    
    path7<-paste(out2,'UpSetR_like_plots','/','enhancer','/','WITH_CTRLS','/', sep='')
    
    cat("path7\n")
    cat(sprintf(as.character(path7)))
    cat("\n")
    
    
    if (file.exists(path7)){
      
      
      
      
    } else {
      dir.create(file.path(path7))
      
    }
    
    
    setwd(path7)
    
    svglite(paste('ACTIVE_TILES_AT_LEAST_',value_arrays_sel,'.svg',sep=''), width = 8, height = 8)
    print(graph_FINAL)
    dev.off()
    
    
    #### WITHOUT CTRLS ----
    
    
    df_sel<-KEY_collapse_melt_sel_enhancer_CLASS[which(KEY_collapse_melt_sel_enhancer_CLASS$value%in%accepted_values),]
    
    cat("df_sel_0\n")
    cat(str(df_sel))
    cat("\n")
    
    df_sel<-df_sel[which(df_sel$Label_2 == "ASSAYED_VARIANT"),]
    
    cat("df_sel_1\n")
    cat(str(df_sel))
    cat("\n")
    
    df_sel.dt<-data.table(df_sel,key=c("KEY_Plus_carried_variants","KEY","carried_variants","VAR","chr","Label","Label_2","VEP_DEF_LABELS"))
    
    cat("df_sel.dt\n")
    cat(str(df_sel.dt))
    cat("\n")
    
    
    
    df_sel_CT_summarised<-as.data.frame(df_sel.dt[,.(Cell_Type_string=paste(Cell_Type,collapse="|")), by=key(df_sel.dt)],stringsAsFactors=F)
    
    cat("df_sel_CT_summarised\n")
    cat(str(df_sel_CT_summarised))
    cat("\n")
    
    
    df_sel_CT_summarised$Classif_DEF<-"NA"
    
    df_sel_CT_summarised$Classif_DEF<-factor(df_sel_CT_summarised$Cell_Type_string,
                                             levels = c("K562|CHRF|HL60|THP1",
                                                        "K562|CHRF|HL60","K562|CHRF|THP1","K562|HL60|THP1","CHRF|HL60|THP1",
                                                        "K562|CHRF","K562|HL60","K562|THP1",
                                                        "CHRF|THP1","CHRF|HL60",
                                                        "HL60|THP1",
                                                        "K562","CHRF","HL60","THP1"),
                                             ordered=T)
    
    cat(sprintf(as.character(names(summary(df_sel_CT_summarised$Classif_DEF)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_sel_CT_summarised$Classif_DEF))))
    cat("\n")
    
    df_sel_CT_summarised<-droplevels(df_sel_CT_summarised)
    
    
    cat("df_sel_CT_summarised_droplevels\n")
    cat(str(df_sel_CT_summarised))
    cat("\n")
    
    cat(sprintf(as.character(names(summary(df_sel_CT_summarised$Classif_DEF)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_sel_CT_summarised$Classif_DEF))))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_sel_CT_summarised$VEP_DEF_LABELS)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_sel_CT_summarised$VEP_DEF_LABELS))))
    cat("\n")
    
    ####
    
    df_CSQ_colors_sel<-df_CSQ_colors_sel[which(df_CSQ_colors_sel$VEP_DEF_LABELS%in%df_sel_CT_summarised$VEP_DEF_LABELS),]
    df_CSQ_colors_sel<-droplevels(df_CSQ_colors_sel)
    
    
    cat("df_CSQ_colors_sel\n")
    cat(str(df_CSQ_colors_sel))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_CSQ_colors_sel$VEP_DEF_LABELS)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_CSQ_colors_sel$VEP_DEF_LABELS))))
    cat("\n")
    
    
    df_sel<-droplevels(df_sel)
    
    
    cat("df_sel\n")
    cat(str(df_sel))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_sel$VEP_DEF_LABELS)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_sel$VEP_DEF_LABELS))))
    cat("\n")
    
    
    
    #### Freq table ----
    
    
    
    df_sel_CT_summarised_Classif_DEF.dt<-data.table(df_sel_CT_summarised, key=c("VEP_DEF_LABELS","Classif_DEF"))
    
    
    df_sel_CT_summarised_Classif_DEF_Fq<-as.data.frame(df_sel_CT_summarised_Classif_DEF.dt[,.N,by=key(df_sel_CT_summarised_Classif_DEF.dt)], stringsAsFactors=F)
    
    colnames(df_sel_CT_summarised_Classif_DEF_Fq)[which(colnames(df_sel_CT_summarised_Classif_DEF_Fq) == "N")]<-"instances"
    
    cat("df_sel_CT_summarised_Classif_DEF_Fq_0\n")
    cat(str(df_sel_CT_summarised_Classif_DEF_Fq))
    cat("\n")
    
    
    
    df_sel_CT_summarised_Classif_DEF_TOTAL.dt<-data.table(df_sel_CT_summarised, key=c("Classif_DEF"))
    
    
    df_sel_CT_summarised_Classif_DEF_TOTAL<-as.data.frame(df_sel_CT_summarised_Classif_DEF_TOTAL.dt[,.N,by=key(df_sel_CT_summarised_Classif_DEF_TOTAL.dt)], stringsAsFactors=F)
    
    colnames(df_sel_CT_summarised_Classif_DEF_TOTAL)[which(colnames(df_sel_CT_summarised_Classif_DEF_TOTAL) == "N")]<-"TOTAL"
    
    cat("df_sel_CT_summarised_Classif_DEF_TOTAL_\n")
    cat(str(df_sel_CT_summarised_Classif_DEF_TOTAL))
    cat("\n")
    
    df_sel_CT_summarised_Classif_DEF_Fq<-merge(df_sel_CT_summarised_Classif_DEF_Fq,
                                               df_sel_CT_summarised_Classif_DEF_TOTAL,
                                               by="Classif_DEF")
    
    df_sel_CT_summarised_Classif_DEF_Fq$Perc<-round(100*(df_sel_CT_summarised_Classif_DEF_Fq$instances/df_sel_CT_summarised_Classif_DEF_Fq$TOTAL),1)
    
    
    cat("df_sel_CT_summarised_Classif_DEF_Fq_1\n")
    cat(str(df_sel_CT_summarised_Classif_DEF_Fq))
    cat("\n")
    
    
    
    
    VEP_LABELS_df.dt<-data.table(df_sel_CT_summarised, key=c("VEP_DEF_LABELS"))
    
    
    VEP_LABELS_df<-as.data.frame(VEP_LABELS_df.dt[,.N,by=key(VEP_LABELS_df.dt)], stringsAsFactors=F)
    
    colnames(VEP_LABELS_df)[which(colnames(VEP_LABELS_df) == "N")]<-"TOTAL"
    
    cat("VEP_LABELS_df_\n")
    cat(str(VEP_LABELS_df))
    cat("\n")
    
    
    
    
    step<-round(max(df_sel_CT_summarised_Classif_DEF_Fq$TOTAL)/10,0)
    
    cat("--------->\t")
    cat(sprintf(as.character(step)))
    cat("\n")
    
    if(step == 0){
      
      step=1
    }
    
    
    breaks.y<-seq(0,max(df_sel_CT_summarised_Classif_DEF_Fq$TOTAL)+step, by=step)
    labels.y<-as.character(breaks.y)
    
    
    cat(sprintf(as.character(breaks.y)))
    cat("\n")
    
    ### graph
    
    # ,
    # fill=VEP_DEF_LABELS
    # scale_fill_manual(values=df_CSQ_colors_sel$color, drop=F,
    #                   name="Totals", breaks=df_CSQ_colors_sel$VEP_DEF_LABELS,
    #                   labels=paste(df_CSQ_colors_sel$VEP_DEF_LABELS,
    #                                df_CSQ_colors_sel$TOTAL, sep =' n= '))+
    
    graph<-ggplot(data=df_sel_CT_summarised_Classif_DEF_TOTAL,
                  aes(x=Classif_DEF,
                      y=TOTAL)) +
      geom_bar(stat="identity",colour='black')+
      theme_bw()+
      theme_classic()+
      theme(axis.title.y=element_text(size=18, family="sans"),
            axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
            axis.text.x=element_text(angle=90,vjust=1,hjust=1,size=8, color="black", family="sans"))+
      scale_y_continuous(name=paste("Variants with at least",value_arrays_sel,"ACTIVE TILES",sep=" "),breaks=breaks.y,labels=labels.y,
                         limits=c(breaks.y[1],breaks.y[length(breaks.y)]))+
      scale_x_discrete(name=NULL, drop=F)+
      theme(legend.position="hidden")+
      ggeasy::easy_center_title()
    
    ### subgraph
    
    Cell_Type.dt<-data.table(df_sel,key=c("Cell_Type"))
    
    cat("Cell_Type.dt\n")
    cat(str(Cell_Type.dt))
    cat("\n")
    
    
    
    Cell_Type_df<-as.data.frame(Cell_Type.dt[,.N, by=key(Cell_Type.dt)],stringsAsFactors=F)
    colnames(Cell_Type_df)[which(colnames(Cell_Type_df) == "N")]<-"Total"
    
    
    
    cat("Cell_Type_df_0\n")
    cat(str(Cell_Type_df))
    cat("\n")
    
    Cell_Type_df<-merge(Cell_Type_df,
                        df_Cell_colors,
                        by="Cell_Type",
                        all=T)
    
    Cell_Type_df$Total[is.na(Cell_Type_df$Total)]<-0
    
    
    cat("Cell_Type_df_1\n")
    cat(str(Cell_Type_df))
    cat("\n")
    
    
    step<-round(max(Cell_Type_df$Total)/5,0)
    
    cat("--------->\t")
    cat(sprintf(as.character(step)))
    cat("\n")
    
    if(step == 0){
      
      step=1
    }
    breaks.x<-rev(c(seq(0,max(Cell_Type_df$Total)+step, by=step)))
    labels.x<-as.character(breaks.x)
    
    
    cat(sprintf(as.character(breaks.x)))
    cat("\n")
    
    levels_CT<-rev(levels(Cell_Type_df$Cell_Type))
    
    cat("levels_CT_2\n")
    cat(str(levels_CT))
    cat("\n")
    
    Cell_Type_df$Cell_type<-factor(Cell_Type_df$Cell_Type,
                                   levels=levels_CT,
                                   ordered=T)
    
    # Cell_Type_df<-Cell_Type_df[order(Cell_Type_df$Cell_Type),]
    
    
    
    cat("Cell_Type_df_2\n")
    cat(str(Cell_Type_df))
    cat("\n")
    
    
    subgraph<-ggplot(data=Cell_Type_df,
                     aes(y=Cell_Type,
                         x=Total,
                         fill=Cell_Type)) +
      geom_bar(stat="identity",colour='black')+
      theme_bw()+
      theme_classic()+
      theme(axis.title.y=element_text(size=16, family="sans"),
            axis.text.y=element_text(angle=0,size=16, color="black", family="sans"),
            axis.text.x=element_text(angle=0,size=10, color="black", family="sans"))+
      scale_x_reverse(name=paste("Variants with at least",value_arrays_sel,"ACTIVE TILES",sep=" "),breaks=breaks.x,labels=labels.x,
                      limits=c(breaks.x[1],breaks.x[length(breaks.x)]))+
      scale_y_discrete(position="right",name=NULL, drop=F)+
      scale_fill_manual(values=Cell_Type_df$colors, drop=F)+
      theme(legend.position="hidden")+
      ggeasy::easy_center_title()
    
    cat("subgraph DONE WITHOUT_CTRLS\n")
    
    
    
    
    
    graph_FINAL<-plot_grid(NULL,graph,subgraph,NULL,
                           nrow = 2,
                           ncol=2,
                           rel_heights = c(1, 0.25),
                           rel_widths=c(0.5,1))
    
    path8<-paste(out2,'UpSetR_like_plots','/','enhancer','/','WITHOUT_CTRLS','/', sep='')
    
    cat("path8\n")
    cat(sprintf(as.character(path8)))
    cat("\n")
    
    
    if (file.exists(path8)){
      
      
      
      
    } else {
      dir.create(file.path(path8))
      
    }
    
    
    setwd(path8)
    
    svglite(paste('ACTIVE_TILES_AT_LEAST_',value_arrays_sel,'.svg',sep=''), width = 8, height = 8)
    print(graph_FINAL)
    dev.off()
    
    
    
    
    
    
    
  }
  
  
  
  
  
  
  
}

Tier_Printer_E_Plus_ASE = function(option_list)
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
  
  #### Categories colors ----
  
  
  df_CSQ_colors<-readRDS(file=opt$CSQ_colors)
  
  df_CSQ_colors$color[df_CSQ_colors$VEP_DEF_LABELS == "TFBS"]<-"red"
  
  
  cat("df_CSQ_colors_0\n")
  cat(str(df_CSQ_colors))
  cat("\n")
  
  #### READ DATA----
  
  setwd(out2)
  
  filename<-paste('Element_collapse','.rds',sep='')
  
  
  
  KEY_collapse<-readRDS(file=filename)
  
  
  KEY_collapse<-KEY_collapse[order(KEY_collapse$KEY_Plus_carried_variants,KEY_collapse$Cell_Type),]
  
  cat("KEY_collapse_0\n")
  cat(str(KEY_collapse))
  cat("\n")
  
  
  filename<-paste('df_Cell_colors','.rds',sep='')
  
  df_Cell_colors<-readRDS(file=filename)
  
  cat("df_Cell_colors_0\n")
  cat(str(df_Cell_colors))
  cat("\n")
  

  KEY_collapse_melt<-melt(KEY_collapse,id.vars=c("KEY_Plus_carried_variants","KEY","carried_variants","VAR","chr","Label","Label_2","VEP_DEF_LABELS","Cell_Type"))
  
  cat("KEY_collapse_melt\n")
  cat(str(KEY_collapse_melt))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(KEY_collapse_melt$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(KEY_collapse_melt$variable)))))
  cat("\n")
  
  
  
  KEY_collapse_melt_sel_E_Plus_ASE_CLASS<-KEY_collapse_melt[which(KEY_collapse_melt$variable == "E_Plus_ASE_CLASS"),]
  
  
  cat("KEY_collapse_melt_sel_E_Plus_ASE_CLASS\n")
  cat(str(KEY_collapse_melt_sel_E_Plus_ASE_CLASS))
  cat("\n")
  
  
  KEY_collapse_melt_sel_E_Plus_ASE_CLASS$value<-factor(as.character(KEY_collapse_melt_sel_E_Plus_ASE_CLASS$value),
                                                       levels = c("0","1","2","3 or >"),
                                                       ordered=T)
  
  
  cat("KEY_collapse_melt_sel_E_Plus_ASE_CLASS_2\n")
  cat(str(KEY_collapse_melt_sel_E_Plus_ASE_CLASS))
  cat("\n")
  
  
  value_arrays<-levels(KEY_collapse_melt_sel_E_Plus_ASE_CLASS$value)
  
  cat("value_arrays\n")
  cat(str(value_arrays))
  cat("\n")
  
  AT_LEAST_category<-NULL
  
  
  for(iteration_value_arrays in 1:length(value_arrays))
  {
    
    value_arrays_sel<-value_arrays[iteration_value_arrays]
    
    cat("------------------------------------------------------------------>\t")
    cat(sprintf(as.character(value_arrays_sel)))
    cat("\n")
    
    
    if(value_arrays_sel == "0")
    {
      
      accepted_values<-"0"
      AT_LEAST_category<-"0"
      
      
    }
    if(value_arrays_sel == "1")
    {
      
      accepted_values<-c("1","2","3 or >")
      AT_LEAST_category<-"AT_LEAST_1"
      
    }
    if(value_arrays_sel == "2")
    {
      
      accepted_values<-c("2","3 or >")
      AT_LEAST_category<-"AT_LEAST_2"
    }
    if(value_arrays_sel == "3 or >")
    {
      
      accepted_values<-c("3 or >")
      AT_LEAST_category<-"AT_LEAST_3"
    }
    
    
    
    
    df_sel<-KEY_collapse_melt_sel_E_Plus_ASE_CLASS[which(KEY_collapse_melt_sel_E_Plus_ASE_CLASS$value%in%accepted_values),]
    
    df_sel<-droplevels(df_sel)
    
    cat("df_sel_0\n")
    cat(str(df_sel))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_sel$VEP_DEF_LABELS)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_sel$VEP_DEF_LABELS))))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_sel$Label_2)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_sel$Label_2))))
    cat("\n")
    
    #### RMV CTRLS ----
    
    df_sel<-df_sel[which(df_sel$Label_2 == "ASSAYED_VARIANT"),]
    
    df_sel<-droplevels(df_sel)
    
    cat("df_sel_1\n")
    cat(str(df_sel))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_sel$VEP_DEF_LABELS)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_sel$VEP_DEF_LABELS))))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_sel$Label_2)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_sel$Label_2))))
    cat("\n")
    
    #### COLLAPSE TO ELEMENT ----
    
    df_sel$VAR_derived_carried_variants<-gsub("^chr","",df_sel$VAR)
    
    df_sel<-df_sel[which(df_sel$VAR_derived_carried_variants == df_sel$carried_variants),]
    
    df_sel<-droplevels(df_sel)
    
    cat("df_sel_2\n")
    cat(str(df_sel))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_sel$VEP_DEF_LABELS)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_sel$VEP_DEF_LABELS))))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_sel$Label_2)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_sel$Label_2))))
    cat("\n")
    
    indx.deplete<-c(which(colnames(df_sel) == "KEY_Plus_carried_variants"),which(colnames(df_sel) == "carried_variants"),which(colnames(df_sel) == "chr"),
                    which(colnames(df_sel) == "Label_2"),which(colnames(df_sel) == "Label"),which(colnames(df_sel) == "variable"),which(colnames(df_sel) == "value"))
    
    df_sel<-unique(df_sel[,-indx.deplete])
    
    indx.int<-c(which(colnames(df_sel) == "KEY"),which(colnames(df_sel) == "VAR"),which(colnames(df_sel) == "VEP_DEF_LABELS"))
    
    
    df_sel<-droplevels(df_sel)
    
    df_sel$value<-1
    
    cat("df_sel_3\n")
    cat(str(df_sel))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_sel$VEP_DEF_LABELS)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_sel$VEP_DEF_LABELS))))
    cat("\n")
    
    
    df_sel_wide<-as.data.frame(pivot_wider(df_sel,
                                           id_cols=colnames(df_sel)[indx.int],
                                           names_from=Cell_Type,
                                           values_from=value),stringsAsFactors=F)
    
    cat("df_sel_wide_0\n")
    cat(str(df_sel_wide))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_sel_wide$VEP_DEF_LABELS)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_sel_wide$VEP_DEF_LABELS))))
    cat("\n")
    
    
    
    df_sel_wide[is.na(df_sel_wide)]<-0
    
    
    cat("df_sel_wide_1\n")
    cat(str(df_sel_wide))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_sel_wide$VEP_DEF_LABELS)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_sel_wide$VEP_DEF_LABELS))))
    cat("\n")
    
    if(value_arrays_sel == "0")
    {
     
      
      
      df_sel_wide_ALL_0<-df_sel_wide[which(df_sel_wide$K562 != 0 &
                                             df_sel_wide$CHRF != 0 &
                                             df_sel_wide$HL60 != 0 &
                                             df_sel_wide$THP1 != 0),]
      
      
      cat("df_sel_wide_ALL_0\n")
      cat(str(df_sel_wide_ALL_0))
      cat("\n")
      cat(sprintf(as.character(names(summary(df_sel_wide_ALL_0$VEP_DEF_LABELS)))))
      cat("\n")
      cat(sprintf(as.character(summary(df_sel_wide_ALL_0$VEP_DEF_LABELS))))
      cat("\n")
      cat(str(unique(df_sel_wide_ALL_0$VAR)))
      cat("\n")
      
      df_sel_wide<-df_sel_wide_ALL_0
      
      
      cat("df_sel_wide_2\n")
      cat(str(df_sel_wide))
      cat("\n")
      cat(sprintf(as.character(names(summary(df_sel_wide$VEP_DEF_LABELS)))))
      cat("\n")
      cat(sprintf(as.character(summary(df_sel_wide$VEP_DEF_LABELS))))
      cat("\n")
      cat(str(unique(df_sel_wide$VAR)))
      cat("\n")
      
      # quit(status = 1)
      
      
    }
    
    #### SAVE----
    
    
    setwd(out2)
    
    
    saveRDS(df_sel_wide,file=paste("E_Plus_ASE_ACTIVE_",AT_LEAST_category,".rds",sep=''))
   
  }#iteration_value_arrays
}

Tier_Printer_enhancer = function(option_list)
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
  
  #### Categories colors ----
  
  
  df_CSQ_colors<-readRDS(file=opt$CSQ_colors)
  
  df_CSQ_colors$color[df_CSQ_colors$VEP_DEF_LABELS == "TFBS"]<-"red"
  
  
  cat("df_CSQ_colors_0\n")
  cat(str(df_CSQ_colors))
  cat("\n")
  
  #### READ DATA----
  
  setwd(out2)
  
  filename<-paste('Element_collapse','.rds',sep='')
  
  
  
  KEY_collapse<-readRDS(file=filename)
  
  
  KEY_collapse<-KEY_collapse[order(KEY_collapse$KEY_Plus_carried_variants,KEY_collapse$Cell_Type),]
  
  cat("KEY_collapse_0\n")
  cat(str(KEY_collapse))
  cat("\n")
  
  
  filename<-paste('df_Cell_colors','.rds',sep='')
  
  df_Cell_colors<-readRDS(file=filename)
  
  cat("df_Cell_colors_0\n")
  cat(str(df_Cell_colors))
  cat("\n")
  
  
  KEY_collapse_melt<-melt(KEY_collapse,id.vars=c("KEY_Plus_carried_variants","KEY","carried_variants","VAR","chr","Label","Label_2","VEP_DEF_LABELS","Cell_Type"))
  
  cat("KEY_collapse_melt\n")
  cat(str(KEY_collapse_melt))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(KEY_collapse_melt$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(KEY_collapse_melt$variable)))))
  cat("\n")
  
  
  
  KEY_collapse_melt_sel_enhancer_CLASS<-KEY_collapse_melt[which(KEY_collapse_melt$variable == "enhancer_CLASS"),]
  
  
  cat("KEY_collapse_melt_sel_enhancer_CLASS\n")
  cat(str(KEY_collapse_melt_sel_enhancer_CLASS))
  cat("\n")
  
  
  KEY_collapse_melt_sel_enhancer_CLASS$value<-factor(as.character(KEY_collapse_melt_sel_enhancer_CLASS$value),
                                                       levels = c("0","1","2","3 or >"),
                                                       ordered=T)
  
  
  cat("KEY_collapse_melt_sel_enhancer_CLASS_2\n")
  cat(str(KEY_collapse_melt_sel_enhancer_CLASS))
  cat("\n")
  
  
  value_arrays<-levels(KEY_collapse_melt_sel_enhancer_CLASS$value)
  
  cat("value_arrays\n")
  cat(str(value_arrays))
  cat("\n")
  
  AT_LEAST_category<-NULL
  
  
  for(iteration_value_arrays in 1:length(value_arrays))
  {
    
    value_arrays_sel<-value_arrays[iteration_value_arrays]
    
    cat("------------------------------------------------------------------>\t")
    cat(sprintf(as.character(value_arrays_sel)))
    cat("\n")
    
    
    if(value_arrays_sel == "0")
    {
      
      accepted_values<-"0"
      AT_LEAST_category<-"0"
      
      
    }
    if(value_arrays_sel == "1")
    {
      
      accepted_values<-c("1","2","3 or >")
      AT_LEAST_category<-"AT_LEAST_1"
      
    }
    if(value_arrays_sel == "2")
    {
      
      accepted_values<-c("2","3 or >")
      AT_LEAST_category<-"AT_LEAST_2"
    }
    if(value_arrays_sel == "3 or >")
    {
      
      accepted_values<-c("3 or >")
      AT_LEAST_category<-"AT_LEAST_3"
    }
    
    
    
    
    df_sel<-KEY_collapse_melt_sel_enhancer_CLASS[which(KEY_collapse_melt_sel_enhancer_CLASS$value%in%accepted_values),]
    
    df_sel<-droplevels(df_sel)
    
    cat("df_sel_0\n")
    cat(str(df_sel))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_sel$VEP_DEF_LABELS)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_sel$VEP_DEF_LABELS))))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_sel$Label_2)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_sel$Label_2))))
    cat("\n")
    
    #### RMV CTRLS ----
    
    df_sel<-df_sel[which(df_sel$Label_2 == "ASSAYED_VARIANT"),]
    
    df_sel<-droplevels(df_sel)
    
    cat("df_sel_1\n")
    cat(str(df_sel))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_sel$VEP_DEF_LABELS)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_sel$VEP_DEF_LABELS))))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_sel$Label_2)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_sel$Label_2))))
    cat("\n")
    
    #### COLLAPSE TO ELEMENT ----
    
    df_sel$VAR_derived_carried_variants<-gsub("^chr","",df_sel$VAR)
    
    df_sel<-df_sel[which(df_sel$VAR_derived_carried_variants == df_sel$carried_variants),]
    
    df_sel<-droplevels(df_sel)
    
    cat("df_sel_2\n")
    cat(str(df_sel))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_sel$VEP_DEF_LABELS)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_sel$VEP_DEF_LABELS))))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_sel$Label_2)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_sel$Label_2))))
    cat("\n")
    
    indx.deplete<-c(which(colnames(df_sel) == "KEY_Plus_carried_variants"),which(colnames(df_sel) == "carried_variants"),which(colnames(df_sel) == "chr"),
                    which(colnames(df_sel) == "Label_2"),which(colnames(df_sel) == "Label"),which(colnames(df_sel) == "variable"),which(colnames(df_sel) == "value"))
    
    df_sel<-unique(df_sel[,-indx.deplete])
    
    indx.int<-c(which(colnames(df_sel) == "KEY"),which(colnames(df_sel) == "VAR"),which(colnames(df_sel) == "VEP_DEF_LABELS"))
    
    
    df_sel<-droplevels(df_sel)
    
    df_sel$value<-1
    
    cat("df_sel_3\n")
    cat(str(df_sel))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_sel$VEP_DEF_LABELS)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_sel$VEP_DEF_LABELS))))
    cat("\n")
    
    
    df_sel_wide<-as.data.frame(pivot_wider(df_sel,
                                           id_cols=colnames(df_sel)[indx.int],
                                           names_from=Cell_Type,
                                           values_from=value),stringsAsFactors=F)
    
    cat("df_sel_wide_0\n")
    cat(str(df_sel_wide))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_sel_wide$VEP_DEF_LABELS)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_sel_wide$VEP_DEF_LABELS))))
    cat("\n")
    
    
    
    df_sel_wide[is.na(df_sel_wide)]<-0
    
    
    cat("df_sel_wide_1\n")
    cat(str(df_sel_wide))
    cat("\n")
    cat(sprintf(as.character(names(summary(df_sel_wide$VEP_DEF_LABELS)))))
    cat("\n")
    cat(sprintf(as.character(summary(df_sel_wide$VEP_DEF_LABELS))))
    cat("\n")
    
    if(value_arrays_sel == "0")
    {
      
      
      
      df_sel_wide_ALL_0<-df_sel_wide[which(df_sel_wide$K562 != 0 &
                                             df_sel_wide$CHRF != 0 &
                                             df_sel_wide$HL60 != 0 &
                                             df_sel_wide$THP1 != 0),]
      
      
      cat("df_sel_wide_ALL_0\n")
      cat(str(df_sel_wide_ALL_0))
      cat("\n")
      cat(sprintf(as.character(names(summary(df_sel_wide_ALL_0$VEP_DEF_LABELS)))))
      cat("\n")
      cat(sprintf(as.character(summary(df_sel_wide_ALL_0$VEP_DEF_LABELS))))
      cat("\n")
      cat(str(unique(df_sel_wide_ALL_0$VAR)))
      cat("\n")
      
      df_sel_wide<-df_sel_wide_ALL_0
      
      
      cat("df_sel_wide_2\n")
      cat(str(df_sel_wide))
      cat("\n")
      cat(sprintf(as.character(names(summary(df_sel_wide$VEP_DEF_LABELS)))))
      cat("\n")
      cat(sprintf(as.character(summary(df_sel_wide$VEP_DEF_LABELS))))
      cat("\n")
      cat(str(unique(df_sel_wide$VAR)))
      cat("\n")
      
      # quit(status = 1)
      
      
    }
    
    #### SAVE----
    
    
    setwd(out2)
    
    
    saveRDS(df_sel_wide,file=paste("enhancer_ACTIVE_",AT_LEAST_category,".rds",sep=''))
    
  }#iteration_value_arrays
}
# #########################################################################################################################################################
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
  
  
 
  UpSetR_like_E_Plus_ASE(opt)
  UpSetR_like_enhancer(opt)
  Tier_Printer_E_Plus_ASE(opt)
  Tier_Printer_enhancer(opt)
  
}
  
  
  
 

###########################################################################

system.time( main() )
