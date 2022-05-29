

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

library("ggforce", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("ggrepel", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")


opt = NULL

options(warn = 1)



volcano_after_QC2_plots_per_variant = function(option_list)
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

  #### MPRA Real Tile ----
  
  MPRA_Real_tile = readRDS(opt$MPRA_Real_tile_QC2_PASS)#, sep="\t", header = T), stringsAsFactors = F)
  
  cat("MPRA_Real_tile_\n")
  cat(str(MPRA_Real_tile))
  cat("\n")
  
  cat(sprintf(as.character(names(summary(as.factor(MPRA_Real_tile$Label))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(MPRA_Real_tile$Label)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(MPRA_Real_tile$DEF_CLASS))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(MPRA_Real_tile$DEF_CLASS)))))
  cat("\n")
  
  

  #### Cell_Type colors-----
  
  Cell_Type_levels<-c("K562","CHRF","HL60","THP1")
  colors_Cell_Type_levels<-c('#32A852','#1877C9','#553B68','#D45E85','#6DB2EE','#62D07F','#C9244B','#87447B','#D6B8E6')
  
  df.color_Cell_Type<-as.data.frame(cbind(Cell_Type_levels,colors_Cell_Type_levels[1:length(Cell_Type_levels)]), stringsAsFactors=F)
  colnames(df.color_Cell_Type)<-c("Cell_Type","colors")
  
  df.color_Cell_Type$Cell_Type<-factor(df.color_Cell_Type$Cell_Type,
                                       levels=Cell_Type_levels,
                                       ordered=T)
  
  cat("df.color_Cell_Type_\n")
  cat(str(df.color_Cell_Type))
  cat("\n")
  
  
  
  #### LOOP VAR and carried variant ----
  
  
  VARS<-unique(MPRA_Real_tile$VAR)
  
  
  cat("VARS\n")
  cat(str(VARS))
  cat("\n")
  
  #### path2 ----

  path2<-paste(out,'Per_variant_graphs','/', sep='')

  cat("path2\n")
  cat(sprintf(as.character(path2)))
  cat("\n")

  if (file.exists(path2)){

 

  } else {
    dir.create(file.path(path2))

  }
  
  for(i in 1:length(VARS))
  {
    VAR_sel<-as.character(VARS[i])
    
    
    cat("---------------------------------->\t")
    cat(sprintf(as.character(VAR_sel)))
    cat("\n")
    
    MPRA_Real_tile_sel<-MPRA_Real_tile[which(MPRA_Real_tile$VAR == VAR_sel),]
    
    cat("MPRA_Real_tile_sel\n")
    cat(str(MPRA_Real_tile_sel))
    cat("\n")
    
    carried_variants_array<-unique(as.character(MPRA_Real_tile_sel$carried_variants))
    
    cat("carried_variants_array\n")
    cat(str(carried_variants_array))
    cat("\n")
    
    Element_sel<-unique(as.character(MPRA_Real_tile_sel$KEY))
    Label_sel<-unique(as.character(MPRA_Real_tile_sel$Label))
    

    path3<-paste(out,'Per_variant_graphs','/',paste(Element_sel,Label_sel,VAR_sel, sep='_'),'/', sep='')
    
    cat("path3\n")
    cat(sprintf(as.character(path3)))
    cat("\n")
    
    if (file.exists(path3)){
      
    
      
      
    } else {
      dir.create(file.path(path3))
      
    }
    
    for(k in 1:length(carried_variants_array))
    {
     
      carried_variants_sel<-as.character(carried_variants_array[k])
      
      
      cat("----------->\t")
      cat(sprintf(as.character(carried_variants_sel)))
      cat("\n")
      
      path4<-paste(out,'Per_variant_graphs','/',paste(Element_sel,Label_sel,VAR_sel, sep='_'),'/',carried_variants_sel,'/', sep='')
      
      cat("path4\n")
      cat(sprintf(as.character(path4)))
      cat("\n")
      
      if (file.exists(path4)){
        
        
        
        
      } else {
        dir.create(file.path(path4))
        
      }
      
      path_variant_interpretation<-paste("/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/MPRA/",VAR_sel,'/',sep='')
      
      
      if (file.exists(path_variant_interpretation)){
        
        
        
        
      } else {
        dir.create(file.path(path_variant_interpretation))
        
      }
      
      path_variant_interpretation<-paste("/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/MPRA/",VAR_sel,'/',carried_variants_sel,'/',sep='')
      
      
      if (file.exists(path_variant_interpretation)){
        
        
        
        
      } else {
        dir.create(file.path(path_variant_interpretation))
        
      }
      
      
      MPRA_Real_tile_carried_variants_sel<-MPRA_Real_tile_sel[which(MPRA_Real_tile_sel$carried_variants == carried_variants_sel),]
      
      MPRA_Real_tile_carried_variants_sel<-droplevels(MPRA_Real_tile_carried_variants_sel)
      
      cat("MPRA_Real_tile_carried_variants_sel\n")
      cat(str(MPRA_Real_tile_carried_variants_sel))
      cat("\n")
    
      
      MPRA_Real_tile_carried_variants_sel_wide<-as.data.frame(pivot_wider(MPRA_Real_tile_carried_variants_sel,
                                                                                 names_from=variable,
                                                                                 values_from=value), stringsAsFactors=F)
      
      
      cat("MPRA_Real_tile_carried_variants_sel_wide\n")
      cat(str(MPRA_Real_tile_carried_variants_sel_wide))
      cat("\n")
      cat(sprintf(as.character(names(summary(MPRA_Real_tile_carried_variants_sel_wide$DEF_CLASS)))))
      cat("\n")
      cat(sprintf(as.character(summary(MPRA_Real_tile_carried_variants_sel_wide$DEF_CLASS))))
      cat("\n")
      
      A<-summary(MPRA_Real_tile_carried_variants_sel_wide$enhancer_logpval)
      
      cat("A\n")
      cat(sprintf(as.character(names(A))))
      cat("\n")
      cat(sprintf(as.character(A)))
      cat("\n")
      
      max_value<-A[6]
      min_value<-A[1]
      
      step<-round((max_value-min_value)/5,0)
      
      if(step ==0)
      {
        
        step<-1
      }
      
      cat("step_enhancer_logpval\n")
      cat(sprintf(as.character(step)))
      cat("\n")
      
      
      breaks.enhancer_logpval<-unique(seq(min_value,max_value+step, by=step))
      labels.enhancer_logpval<-as.character(round(breaks.enhancer_logpval),2)
      
      cat("labels.enhancer_logpval\n")
      cat(sprintf(as.character(labels.enhancer_logpval)))
      cat("\n")
      
      A<-summary(MPRA_Real_tile_carried_variants_sel_wide$ASE_logpval)
      
      cat("A\n")
      cat(sprintf(as.character(names(A))))
      cat("\n")
      cat(sprintf(as.character(A)))
      cat("\n")
      
      max_value<-A[6]
      min_value<-A[1]
      
      step<-round((max_value-min_value)/5,0)
      
      if(step ==0)
      {
        
        step<-1
      }
      
      cat("step_ASE_logpval\n")
      cat(sprintf(as.character(step)))
      cat("\n")
      
      
      breaks.ASE_logpval<-unique(seq(min_value,max_value+step, by=step))
      labels.ASE_logpval<-as.character(round(breaks.ASE_logpval),2)
      
      cat("labels.ASE_logpval\n")
      cat(sprintf(as.character(labels.ASE_logpval)))
      cat("\n")
      
      A<-summary(MPRA_Real_tile_carried_variants_sel_wide$LogFC)
      
      
      
      cat("A\n")
      cat(sprintf(as.character(names(A))))
      cat("\n")
      cat(sprintf(as.character(A)))
      cat("\n")
      
      
      
      max_log_LogFC<-A[6]
      min_log_LogFC<--1.5
      
      #min_log_LogFC
      breaks.log_LogFC<-sort(c(0,seq(min_log_LogFC,max_log_LogFC+0.5, by=0.5)))
      
      
      labels_LogFC<-as.character(round(logratio2foldchange(breaks.log_LogFC, base=2),2))
      
      
      cat("labels_LogFC\n")
      cat(sprintf(as.character(labels_LogFC)))
      cat("\n")
      
      # MPRA_Real_tile_carried_variants_sel_wide<-merge(MPRA_Real_tile_carried_variants_sel_wide,
      #                                                 df.color_Cell_Type,
      #                                                 by="Cell_Type")
      
      #MPRA_Real_tile_carried_variants_sel_wide<-droplevels(MPRA_Real_tile_carried_variants_sel_wide)
      
      
      shape_values<-c(15,16,17,18)
      
      
      accepted_CLASS_enhancer<-c("enhancer","E_Plus_ASE")
      
      volcano<-ggplot(data=MPRA_Real_tile_carried_variants_sel_wide,
                      aes(x=LogFC, 
                          y=enhancer_logpval,
                          shape=Cell_Type)) +
        geom_point(size=5, color="gray")+
        theme(legend.position="bottom",legend.title=element_blank(), legend.text = element_text(size=10))+
        geom_point(data=MPRA_Real_tile_carried_variants_sel_wide[which(MPRA_Real_tile_carried_variants_sel_wide$DEF_CLASS%in%accepted_CLASS_enhancer),],
                   aes(x=LogFC, 
                       y=enhancer_logpval,
                       color=Cell_Type,
                       shape=Cell_Type),
                   size=12)+
        theme_bw()+
        theme(axis.title.y=element_text(size=24, family="sans"),
              axis.title.x=element_text(size=24, family="sans"),
              axis.text.y=element_text(angle=0,size=24, color="black", family="sans"),
              axis.text.x=element_text(angle=0, size=24, color="black", family="sans"),
              legend.title=element_text(size=24,color="black", family="sans"),
              legend.text=element_text(size=24,color="black", family="sans"))+
        scale_x_continuous(name="FC",breaks=breaks.log_LogFC,labels=labels_LogFC, 
                           limits=c(breaks.log_LogFC[1],breaks.log_LogFC[length(breaks.log_LogFC)])) +
        scale_y_continuous(name="enhancer_logpval",breaks=breaks.enhancer_logpval,labels=labels.enhancer_logpval, 
                            limits=c(breaks.enhancer_logpval[1],breaks.enhancer_logpval[length(breaks.enhancer_logpval)]))+
        geom_vline(xintercept=0,linetype="dotted")+
        geom_hline(yintercept=1.3,linetype="dotted")+
        scale_shape_manual(values = shape_values)+
        scale_color_manual(values=df.color_Cell_Type$colors, drop=F)+
        ggeasy::easy_center_title()
      
   
      volcano<-volcano+
              geom_text_repel(data=MPRA_Real_tile_carried_variants_sel_wide[which(MPRA_Real_tile_carried_variants_sel_wide$DEF_CLASS%in%accepted_CLASS_enhancer),],
                              aes(x=LogFC, 
                                  y=enhancer_logpval,
                                  shape=Cell_Type,
                                  label=paste(TILE,DEF_CLASS,sep=";")),
                      nudge_x = .15,
                      box.padding = 0.5,
                      nudge_y = 1,
                      segment.curvature = -0.1,
                      segment.ncp = 3,
                      segment.angle = 20,
                      max.overlaps = Inf)
      
      setwd(path4)
      
      svgname<-paste("volcano_LogFC",".svg",sep='')
      makesvg = TRUE
      
      if (makesvg == TRUE)
      {
        ggsave(svgname, plot= volcano,
               device="svg",
               height=10, width=12)
      }
      
      setwd(path_variant_interpretation)
      
      
      svgname<-paste("volcano_LogFC",".svg",sep='')
      makesvg = TRUE
      
      if (makesvg == TRUE)
      {
        ggsave(svgname, plot= volcano,
               device="svg",
               height=10, width=12)
      }
      
      ### Vockley REF
      
      
      A<-summary(MPRA_Real_tile_carried_variants_sel_wide$Vockley_REF)
      
      
      
      cat("A\n")
      cat(sprintf(as.character(names(A))))
      cat("\n")
      cat(sprintf(as.character(A)))
      cat("\n")
      
      
      
      max_Vockley_REF<-A[6]
      min_Vockley_REF<-A[1]
      # min_Vockley_REF<-0.4
      # max_Vockley_REF<-3
      
      step<-(max_Vockley_REF-min_Vockley_REF)/5
      
      if(step ==0)
      {
        
        step<-1
      }
      
      breaks.Vockley_REF<-sort(c(1,seq(min_Vockley_REF,max_Vockley_REF,by=step)))
      labels_Vockley_REF<-as.character(round(breaks.Vockley_REF,2))
      
      
      cat("labels_Vockley_REF\n")
      cat(sprintf(as.character(labels_Vockley_REF)))
      cat("\n")
      
      accepted_CLASS_ASE<-c("ASE","E_Plus_ASE")
      
      
      volcano<-ggplot(data=MPRA_Real_tile_carried_variants_sel_wide,
                      aes(x=Vockley_REF, 
                          y=ASE_logpval,
                          shape=Cell_Type)) +
        geom_point(size=5, color="gray")+
        theme(legend.position="bottom",legend.title=element_blank(), legend.text = element_text(size=10))+
        geom_point(data=MPRA_Real_tile_carried_variants_sel_wide[which(MPRA_Real_tile_carried_variants_sel_wide$DEF_CLASS%in%accepted_CLASS_ASE),],
                   aes(x=Vockley_REF, 
                       y=ASE_logpval,
                       color=Cell_Type,
                       shape=Cell_Type),
                   size=12)+
        theme_bw()+
        theme(axis.title.y=element_text(size=24, family="sans"),
              axis.title.x=element_text(size=24, family="sans"),
              axis.text.y=element_text(angle=0,size=24, color="black", family="sans"),
              axis.text.x=element_text(angle=0, size=24, color="black", family="sans"),
              legend.title=element_text(size=24,color="black", family="sans"),
              legend.text=element_text(size=24,color="black", family="sans"))+
        scale_x_continuous(name="Vockley_REF",breaks=breaks.Vockley_REF,labels=labels_Vockley_REF, 
                           limits=c(breaks.Vockley_REF[1],breaks.Vockley_REF[length(breaks.Vockley_REF)])) +
        scale_y_continuous(name="ASE_logpval",breaks=breaks.ASE_logpval,labels=labels.ASE_logpval, 
                           limits=c(breaks.ASE_logpval[1],breaks.ASE_logpval[length(breaks.ASE_logpval)]))+
        geom_vline(xintercept=1,linetype="dotted")+
        geom_hline(yintercept=1.3,linetype="dotted")+
        scale_shape_manual(values = shape_values)+
        scale_color_manual(values=df.color_Cell_Type$colors, drop=F)+
        ggeasy::easy_center_title()
      
      
      volcano<-volcano+
        geom_text_repel(data=MPRA_Real_tile_carried_variants_sel_wide[which(MPRA_Real_tile_carried_variants_sel_wide$DEF_CLASS%in%accepted_CLASS_ASE),],
                        aes(x=Vockley_REF, 
                            y=ASE_logpval,
                            shape=Cell_Type,
                            label=paste(TILE,DEF_CLASS,sep=";")),
                        nudge_x = .15,
                        box.padding = 0.5,
                        nudge_y = 1,
                        segment.curvature = -0.1,
                        segment.ncp = 3,
                        segment.angle = 20,
                        max.overlaps = Inf)
      
      setwd(path4)
      
      svgname<-paste("volcano_Vockley_REF",".svg",sep='')
      makesvg = TRUE
      
      if (makesvg == TRUE)
      {
        ggsave(svgname, plot= volcano,
               device="svg",
               height=10, width=12)
      }
      
      setwd(path_variant_interpretation)
      
      svgname<-paste("volcano_Vockley_REF",".svg",sep='')
      makesvg = TRUE
      
      if (makesvg == TRUE)
      {
        ggsave(svgname, plot= volcano,
               device="svg",
               height=10, width=12)
      }
      
    }#k carried_variants_array[k]
  }# i VARS VARS[i]
  
  
  
  
  # 
  # 
  # 
  # ####################################################################################################################################################################################
  # quit(status = 1)
  # 
  # 
  # 
 
  
 
  
  
  
  # #####################################################################################################################################################################################################################################################################
  # quit(status = 1)
  
  
  
  
  
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
    make_option(c("--MPRA_Real_tile_QC2_PASS"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
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

  volcano_after_QC2_plots_per_variant(opt)
  
}


###########################################################################

system.time( main() )
