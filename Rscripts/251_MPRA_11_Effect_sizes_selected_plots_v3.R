#### WARNING HARDCODED INDEXES -----


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
library("reshape2",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("viridisLite",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("gtools", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/")


lm_eqn = function(m) {
  
  l <- list(a = format(coef(m)[1], digits = 2),
            b = format(abs(coef(m)[2]), digits = 2),
            R1_15 = format(summary(m)$r.squared, digits = 3));
  
  if (coef(m)[2] >= 0)  {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~R1_15,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~R1_15,l)    
  }
  
  as.character(as.expression(eq));                 
}

opt = NULL

options(warn=1)


graph_function_E_Plus_ASE = function(option_list)
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
  
  ### Read RAW data ----
  
  
  
  path5<-paste(out2,'GWAS_plots','/', sep='')
  
  cat("path5\n")
  cat(sprintf(as.character(path5)))
  cat("\n")
  
  
  if (file.exists(path5)){
    
    
    
    
  } else {
    dir.create(file.path(path5))
    
  }
  
  
  ####    READ ----
  
  setwd(path5)
  
  # quit(status = 1)
  
  
  filename_1<-paste("Element_collapse_Plus_Variant_Phenotypes_Plus_LogFC_ASE",".rds", sep='')
  
  
  Element_collapse<-readRDS(file=filename_1)
  
  # Reclassify WBC
  
  Element_collapse$Lineage[Element_collapse$phenotype == "wbc"]<-"gran_mono_lineage"
  
  cat("Element_collapse\n")
  cat(str(Element_collapse))
  cat("\n")
  
  
  
  
  #### AT_LEAST CLASSIFICATION ----
  
  
  
  Element_collapse$AT_LEAST_CLASS_1<-"NA"
  
  Element_collapse$AT_LEAST_CLASS_1[Element_collapse$E_Plus_ASE_CLASS == "0"]<-"0"
  Element_collapse$AT_LEAST_CLASS_1[Element_collapse$E_Plus_ASE_CLASS == "1"]<-"1"
  Element_collapse$AT_LEAST_CLASS_1[Element_collapse$E_Plus_ASE_CLASS == "2"]<-"1"
  Element_collapse$AT_LEAST_CLASS_1[Element_collapse$E_Plus_ASE_CLASS == "3 or >"]<-"1"
  
  Element_collapse$AT_LEAST_CLASS_2<-"NA"
  
  Element_collapse$AT_LEAST_CLASS_2[Element_collapse$E_Plus_ASE_CLASS == "0"]<-"0"
  Element_collapse$AT_LEAST_CLASS_2[Element_collapse$E_Plus_ASE_CLASS == "2"]<-"2"
  Element_collapse$AT_LEAST_CLASS_2[Element_collapse$E_Plus_ASE_CLASS == "3 or >"]<-"2"
  
  Element_collapse$AT_LEAST_CLASS_3<-"NA"
  
  Element_collapse$AT_LEAST_CLASS_3[Element_collapse$E_Plus_ASE_CLASS == "0"]<-"0"
  Element_collapse$AT_LEAST_CLASS_3[Element_collapse$E_Plus_ASE_CLASS == "3 or >"]<-"3 or >"
  
  
  
  Element_collapse$AT_LEAST_CLASS_1<-factor(Element_collapse$AT_LEAST_CLASS_1,
                                            levels=c("0","1"),
                                            ordered=T)
  
  
  Element_collapse$AT_LEAST_CLASS_2<-factor(Element_collapse$AT_LEAST_CLASS_2,
                                            levels=c("0","2"),
                                            ordered=T)
  
  Element_collapse$AT_LEAST_CLASS_3<-factor(Element_collapse$AT_LEAST_CLASS_3,
                                            levels=c("0","3 or >"),
                                            ordered=T)
  
  Element_collapse$FC<-logratio2foldchange(Element_collapse$LogFC, base=2)
  
  Element_collapse<-Element_collapse[order(Element_collapse$FC, decreasing = F),]
  
  
  cat("Element_collapse_2\n")
  cat(str(Element_collapse))
  cat("\n")
  cat(sprintf(as.character(names(summary(Element_collapse$AT_LEAST_CLASS_1)))))
  cat("\n")
  cat(sprintf(as.character(summary(Element_collapse$AT_LEAST_CLASS_1))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Element_collapse$AT_LEAST_CLASS_2)))))
  cat("\n")
  cat(sprintf(as.character(summary(Element_collapse$AT_LEAST_CLASS_2))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Element_collapse$AT_LEAST_CLASS_3)))))
  cat("\n")
  cat(sprintf(as.character(summary(Element_collapse$AT_LEAST_CLASS_3))))
  cat("\n")
  
  
  Element_collapse_sel<-droplevels(Element_collapse[which(Element_collapse$AT_LEAST_CLASS_1 == "1"),])
  
  cat("Element_collapse_sel_2\n")
  cat(str(Element_collapse_sel))
  cat("\n")
  cat(sprintf(as.character(names(summary(Element_collapse_sel$AT_LEAST_CLASS_1)))))
  cat("\n")
  cat(sprintf(as.character(summary(Element_collapse_sel$AT_LEAST_CLASS_1))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Element_collapse_sel$AT_LEAST_CLASS_2)))))
  cat("\n")
  cat(sprintf(as.character(summary(Element_collapse_sel$AT_LEAST_CLASS_2))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Element_collapse_sel$AT_LEAST_CLASS_3)))))
  cat("\n")
  cat(sprintf(as.character(summary(Element_collapse_sel$AT_LEAST_CLASS_3))))
  cat("\n")
  
  
  lineage_array<-levels(Element_collapse_sel$Lineage)
  
  cat("lineage_array_2\n")
  cat(str(lineage_array))
  cat("\n")
  cat(sprintf(as.character(lineage_array)))
  cat("\n")
  
  colors<-c("firebrick1","brown","magenta","dodgerblue")
  
  
  
  for(i in 1:length(lineage_array))
  {
    color_sel<-colors[i]

    lineage_array_sel<-lineage_array[i] 
    
    Element_collapse_sel_lineage_sel<-droplevels(Element_collapse_sel[which(Element_collapse_sel$Lineage == lineage_array_sel),])
    
    cat("Element_collapse_sel_lineage_sel_0\n")
    cat(str(Element_collapse_sel_lineage_sel))
    cat("\n")
    
    cat("GRAPH_START\n")
    
    
    graph<-ggplot(Element_collapse_sel_lineage_sel, 
                  aes(x=FC, 
                      y=abs_finemap_z)) +
      geom_point(size=4) +
      theme(plot.title = element_text(size=11)) +
      theme(plot.title=element_text(size=11))+
      scale_x_continuous(name="Fold Change")+
      scale_y_continuous(name="abs_finemap_z")+
      theme_bw()+
      theme(legend.position="bottom")+
      facet_grid(Cell_Type ~ phenotype_DEF)+
      ggeasy::easy_center_title()
    
    graph<-graph+
      theme(strip.background = element_rect(color="black", fill=color_sel, size=1.5, linetype="solid"))+
      theme(strip.text.x = element_text(size = 12, color = "white", face="bold"))+
      theme(strip.text.y = element_text(size = 12, color = "white", face="bold"))
      
    
    cat("GRAPH_END\n")
    
    
    setwd(path5)
    
    svgname<-paste("GWAS_effect_size_vs_Fold_change_",lineage_array_sel,"_E_Plus_ASE_AT_LEAST_ONE_CLASS",".svg",sep='')
    makesvg = TRUE
    
    if (makesvg == TRUE)
    {
      ggsave(svgname, plot= graph,
             device="svg",
             height=10, width=12)
    }
    
    graph<-ggplot(Element_collapse_sel_lineage_sel, 
                  aes(x=abs_ASE, 
                      y=abs_finemap_z)) +
      geom_point(size=4) +
      theme(plot.title = element_text(size=11)) +
      theme(plot.title=element_text(size=11))+
      scale_x_continuous(name="Fold Change")+
      scale_y_continuous(name="abs_finemap_z")+
      theme_bw()+
      theme(legend.position="bottom")+
      facet_grid(Cell_Type ~ phenotype_DEF)+
      ggeasy::easy_center_title()
    
    graph<-graph+
      theme(strip.background = element_rect(color="black", fill=color_sel, size=1.5, linetype="solid"))+
      theme(strip.text.x = element_text(size = 12, color = "white", face="bold"))+
      theme(strip.text.y = element_text(size = 12, color = "white", face="bold"))
    
    
    cat("GRAPH_END\n")
    
    
    setwd(path5)
    
    svgname<-paste("GWAS_effect_size_vs_abs_ASE_",lineage_array_sel,"_E_Plus_ASE_AT_LEAST_ONE_CLASS",".svg",sep='')
    makesvg = TRUE
    
    if (makesvg == TRUE)
    {
      ggsave(svgname, plot= graph,
             device="svg",
             height=10, width=12)
    }
    
    
  }#i lineage_array
  
}


graph_function_enhancer = function(option_list)
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
  
  ### Read RAW data ----
  
  #### READ and transform out2 ----
  
  out2 = opt$out2
  
  cat("out2_\n")
  cat(sprintf(as.character(out2)))
  cat("\n")
  
  path5<-paste(out2,'GWAS_plots','/', sep='')
  
  cat("path5\n")
  cat(sprintf(as.character(path5)))
  cat("\n")
  
  
  if (file.exists(path5)){
    
    
    
    
  } else {
    dir.create(file.path(path5))
    
  }
  
  
  ####    READ ----
  
  setwd(path5)
  
  # quit(status = 1)
  
  
  filename_1<-paste("Element_collapse_Plus_Variant_Phenotypes_Plus_LogFC_ASE",".rds", sep='')
  
  
  Element_collapse<-readRDS(file=filename_1)
  
  # Reclassify WBC
  
  Element_collapse$Lineage[Element_collapse$phenotype == "wbc"]<-"gran_mono_lineage"
  
  cat("Element_collapse\n")
  cat(str(Element_collapse))
  cat("\n")
  
  
  
  
  #### AT_LEAST CLASSIFICATION ----
  
  
  
  Element_collapse$AT_LEAST_CLASS_1<-"NA"
  
  Element_collapse$AT_LEAST_CLASS_1[Element_collapse$enhancer_CLASS == "0"]<-"0"
  Element_collapse$AT_LEAST_CLASS_1[Element_collapse$enhancer_CLASS == "1"]<-"1"
  Element_collapse$AT_LEAST_CLASS_1[Element_collapse$enhancer_CLASS == "2"]<-"1"
  Element_collapse$AT_LEAST_CLASS_1[Element_collapse$enhancer_CLASS == "3 or >"]<-"1"
  
  Element_collapse$AT_LEAST_CLASS_2<-"NA"
  
  Element_collapse$AT_LEAST_CLASS_2[Element_collapse$enhancer_CLASS == "0"]<-"0"
  Element_collapse$AT_LEAST_CLASS_2[Element_collapse$enhancer_CLASS == "2"]<-"2"
  Element_collapse$AT_LEAST_CLASS_2[Element_collapse$enhancer_CLASS == "3 or >"]<-"2"
  
  Element_collapse$AT_LEAST_CLASS_3<-"NA"
  
  Element_collapse$AT_LEAST_CLASS_3[Element_collapse$enhancer_CLASS == "0"]<-"0"
  Element_collapse$AT_LEAST_CLASS_3[Element_collapse$enhancer_CLASS == "3 or >"]<-"3 or >"
  
  
  
  Element_collapse$AT_LEAST_CLASS_1<-factor(Element_collapse$AT_LEAST_CLASS_1,
                                            levels=c("0","1"),
                                            ordered=T)
  
  
  Element_collapse$AT_LEAST_CLASS_2<-factor(Element_collapse$AT_LEAST_CLASS_2,
                                            levels=c("0","2"),
                                            ordered=T)
  
  Element_collapse$AT_LEAST_CLASS_3<-factor(Element_collapse$AT_LEAST_CLASS_3,
                                            levels=c("0","3 or >"),
                                            ordered=T)
  
  Element_collapse$FC<-logratio2foldchange(Element_collapse$LogFC, base=2)
  
  Element_collapse<-Element_collapse[order(Element_collapse$FC, decreasing = F),]
  
  
  cat("Element_collapse_2\n")
  cat(str(Element_collapse))
  cat("\n")
  cat(sprintf(as.character(names(summary(Element_collapse$AT_LEAST_CLASS_1)))))
  cat("\n")
  cat(sprintf(as.character(summary(Element_collapse$AT_LEAST_CLASS_1))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Element_collapse$AT_LEAST_CLASS_2)))))
  cat("\n")
  cat(sprintf(as.character(summary(Element_collapse$AT_LEAST_CLASS_2))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Element_collapse$AT_LEAST_CLASS_3)))))
  cat("\n")
  cat(sprintf(as.character(summary(Element_collapse$AT_LEAST_CLASS_3))))
  cat("\n")
  
  
  Element_collapse_sel<-droplevels(Element_collapse[which(Element_collapse$AT_LEAST_CLASS_1 == "1"),])
  
  cat("Element_collapse_sel_2\n")
  cat(str(Element_collapse_sel))
  cat("\n")
  cat(sprintf(as.character(names(summary(Element_collapse_sel$AT_LEAST_CLASS_1)))))
  cat("\n")
  cat(sprintf(as.character(summary(Element_collapse_sel$AT_LEAST_CLASS_1))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Element_collapse_sel$AT_LEAST_CLASS_2)))))
  cat("\n")
  cat(sprintf(as.character(summary(Element_collapse_sel$AT_LEAST_CLASS_2))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Element_collapse_sel$AT_LEAST_CLASS_3)))))
  cat("\n")
  cat(sprintf(as.character(summary(Element_collapse_sel$AT_LEAST_CLASS_3))))
  cat("\n")
  
  
  lineage_array<-levels(Element_collapse_sel$Lineage)
  
  cat("lineage_array_2\n")
  cat(str(lineage_array))
  cat("\n")
  cat(sprintf(as.character(lineage_array)))
  cat("\n")
  
  colors<-c("firebrick1","brown","magenta","dodgerblue")
  
  
  
  for(i in 1:length(lineage_array))
  {
    color_sel<-colors[i]
    
    lineage_array_sel<-lineage_array[i] 
    
    Element_collapse_sel_lineage_sel<-droplevels(Element_collapse_sel[which(Element_collapse_sel$Lineage == lineage_array_sel),])
    
    cat("Element_collapse_sel_lineage_sel_0\n")
    cat(str(Element_collapse_sel_lineage_sel))
    cat("\n")
    
    cat("GRAPH_START\n")
    
    
    graph<-ggplot(Element_collapse_sel_lineage_sel, 
                  aes(x=FC, 
                      y=abs_finemap_z)) +
      geom_point(size=4) +
      theme(plot.title = element_text(size=11)) +
      theme(plot.title=element_text(size=11))+
      scale_x_continuous(name="Fold Change")+
      scale_y_continuous(name="abs_finemap_z")+
      theme_bw()+
      theme(legend.position="bottom")+
      facet_grid(Cell_Type ~ phenotype_DEF)+
      ggeasy::easy_center_title()
    
    graph<-graph+
      theme(strip.background = element_rect(color="black", fill=color_sel, size=1.5, linetype="solid"))+
      theme(strip.text.x = element_text(size = 12, color = "white", face="bold"))+
      theme(strip.text.y = element_text(size = 12, color = "white", face="bold"))
    
    
    cat("GRAPH_END\n")
    
    
    setwd(path5)
    
    svgname<-paste("GWAS_effect_size_vs_Fold_change_",lineage_array_sel,"_enhancer_AT_LEAST_ONE_CLASS",".svg",sep='')
    makesvg = TRUE
    
    if (makesvg == TRUE)
    {
      ggsave(svgname, plot= graph,
             device="svg",
             height=10, width=12)
    }
    
    graph<-ggplot(Element_collapse_sel_lineage_sel, 
                  aes(x=abs_ASE, 
                      y=abs_finemap_z)) +
      geom_point(size=4) +
      theme(plot.title = element_text(size=11)) +
      theme(plot.title=element_text(size=11))+
      scale_x_continuous(name="Fold Change")+
      scale_y_continuous(name="abs_finemap_z")+
      theme_bw()+
      theme(legend.position="bottom")+
      facet_grid(Cell_Type ~ phenotype_DEF)+
      ggeasy::easy_center_title()
    
    graph<-graph+
      theme(strip.background = element_rect(color="black", fill=color_sel, size=1.5, linetype="solid"))+
      theme(strip.text.x = element_text(size = 12, color = "white", face="bold"))+
      theme(strip.text.y = element_text(size = 12, color = "white", face="bold"))
    
    
    cat("GRAPH_END\n")
    
    
    setwd(path5)
    
    svgname<-paste("GWAS_effect_size_vs_abs_ASE_",lineage_array_sel,"_enhancer_AT_LEAST_ONE_CLASS",".svg",sep='')
    makesvg = TRUE
    
    if (makesvg == TRUE)
    {
      ggsave(svgname, plot= graph,
             device="svg",
             height=10, width=12)
    }
    
    
  }#i lineage_array
  
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
    make_option(c("--TOME_correspondence"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--dB"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--MPRA_Real_tile_QC2_PASS"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--KEY_collpase_Plus_Variant_lineage_CLASSIFICATION"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--finemap_prob_Threshold"), type="numeric", default=NULL, 
                metavar="filename", 
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
  
  
 
  
  graph_function_E_Plus_ASE(opt)
  graph_function_enhancer(opt)
 
 
}
  
  
  
 

###########################################################################

system.time( main() )
