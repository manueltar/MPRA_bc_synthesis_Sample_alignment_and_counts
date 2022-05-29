
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
suppressMessages(library("ggrepel", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
library("gtools", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")


library("reshape2", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
suppressMessages(library("lemon", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("gridBase", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("gridExtra", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("grid", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))

opt = NULL

options(warn = 1)

Generating_count_matrix_NSCS = function(option_list)
{
  #### READ and transform type ----
  
  type = opt$type
  
  cat("HELLOW_WORLD_1_TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  
  #### Sankaran MPRA variants ----
  
  Sankaran_MPRA<-as.data.frame(fread(file=opt$Sankaran_MPRA, sep=",",
                                     header=T), stringsAsFactors=F)
  
  Sankaran_MPRA$Carried_variants<-paste(Sankaran_MPRA$chr,Sankaran_MPRA$pos,Sankaran_MPRA$ref,Sankaran_MPRA$alt, sep="_")
  
  
  
  cat("Sankaran_MPRA_\n")
  str(Sankaran_MPRA)
  cat("\n")
  
  Sankaran_MPRA_Sig<-Sankaran_MPRA#[which(Sankaran_MPRA$CTRL.padj <= 0.05),]
  
  
  Sankaran_MPRA_Sig$Tile<-"NA"
  
  Sankaran_MPRA_Sig$Tile[grep('1/2',Sankaran_MPRA_Sig$construct)]<-"HALF"
  Sankaran_MPRA_Sig$Tile[grep('1/3',Sankaran_MPRA_Sig$construct)]<-"ONE_THIRD"
  Sankaran_MPRA_Sig$Tile[grep('2/3',Sankaran_MPRA_Sig$construct)]<-"TWO_THIRDS"
  
  Sankaran_MPRA_Sig$ALLELE<-"NA"
  
  Sankaran_MPRA_Sig$ALLELE[grep('Mut',Sankaran_MPRA_Sig$type)]<-"ALT"
  Sankaran_MPRA_Sig$ALLELE[grep('Ref',Sankaran_MPRA_Sig$type)]<-"REF"
  
  
  cat("Sankaran_MPRA_Sig_\n")
  str(Sankaran_MPRA_Sig)
  cat("\n")
  
  indx.deplete<-grep("GATA1",colnames(Sankaran_MPRA_Sig))
  indx.keep<-grep("GATA1.fc",colnames(Sankaran_MPRA_Sig))
  
  indx.int<-c(which(colnames(Sankaran_MPRA_Sig) == "Carried_variants"),
              which(colnames(Sankaran_MPRA_Sig) == "dbSNP"),
              which(colnames(Sankaran_MPRA_Sig) == "Tile"),
              which(colnames(Sankaran_MPRA_Sig) == "CTRL.median"),
              which(colnames(Sankaran_MPRA_Sig) == "CTRL.padj"),
              which(colnames(Sankaran_MPRA_Sig) == "CTRL.mut.padj"),
              which(colnames(Sankaran_MPRA_Sig) == "CTRL.fc"),
              which(colnames(Sankaran_MPRA_Sig) == "GATA1.fc"),
              which(colnames(Sankaran_MPRA_Sig) == "ALLELE"))
  
  Sankaran_MPRA_Sig_subset<-unique(Sankaran_MPRA_Sig[,indx.int])
  
  cat("Sankaran_MPRA_Sig_subset_\n")
  str(Sankaran_MPRA_Sig_subset)
  cat("\n")
  
  
  Sankaran_MPRA_Sig_subset_wide<-as.data.frame(pivot_wider(Sankaran_MPRA_Sig_subset,names_from=ALLELE,
                                                           values_from=c("CTRL.median","CTRL.padj","CTRL.mut.padj","CTRL.fc","GATA1.fc"),
                                                           id_cols=c("Carried_variants","Tile","dbSNP")), stringsAsFactors=F)
  
  
  # Sankaran_MPRA_Sig_subset_wide$Delta_allelic_skew<-Sankaran_MPRA_Sig_subset_wide$CTRL.median_REF-Sankaran_MPRA_Sig_subset_wide$CTRL.median_ALT
  
  cat("Sankaran_MPRA_Sig_subset_wide_\n")
  str(Sankaran_MPRA_Sig_subset_wide)
  cat("\n")
  # 
  Sankaran_MPRA_Sig_subset_wide$Allelic_skew<-round(log2(Sankaran_MPRA_Sig_subset_wide$CTRL.median_REF/Sankaran_MPRA_Sig_subset_wide$CTRL.median_ALT),2)
  # 
  # 
  cat("Sankaran_MPRA_Sig_subset_wide_\n")
  str(Sankaran_MPRA_Sig_subset_wide)
  cat("\n")
  
  Sankaran_MPRA_Sig_subset_wide$CLASSIF<-"NA"
  
  Sankaran_MPRA_Sig_subset_wide$CLASSIF[which(Sankaran_MPRA_Sig_subset_wide$Allelic_skew < 0)]<-"ALT_biased"
  Sankaran_MPRA_Sig_subset_wide$CLASSIF[which(Sankaran_MPRA_Sig_subset_wide$Allelic_skew > 0)]<-"REF_biased"
  
  Sankaran_MPRA_Sig_subset_wide$CLASSIF<-factor(Sankaran_MPRA_Sig_subset_wide$CLASSIF,
                                                levels=c("REF_biased","ALT_biased"),
                                                ordered=T)
  
  mean_Allelic_skew_GLOBAL<-mean(Sankaran_MPRA_Sig_subset_wide$Allelic_skew,na.rm = T)
  sd_Allelic_skew_GLOBAL<-sd(Sankaran_MPRA_Sig_subset_wide$Allelic_skew,na.rm = T)
  
  cat("mean_Allelic_skew_GLOBAL\n")
  cat(str(mean_Allelic_skew_GLOBAL))
  cat("\n")
  
  cat("sd_Allelic_skew_GLOBAL\n")
  cat(str(sd_Allelic_skew_GLOBAL))
  cat("\n")
  
  Sankaran_MPRA_Sig_subset_wide.dt<-data.table(Sankaran_MPRA_Sig_subset_wide, key=c("Carried_variants","Tile"))
  
  cat("Sankaran_MPRA_Sig_subset_wide.dt\n")
  cat(str(Sankaran_MPRA_Sig_subset_wide.dt))
  cat("\n")
  
  Sankaran_MPRA_Sig_subset_wide.dt.Zscore<-as.data.frame(Sankaran_MPRA_Sig_subset_wide.dt[,.('Allelic_skew_Zscore'=((Allelic_skew-mean_Allelic_skew_GLOBAL)/sd_Allelic_skew_GLOBAL)),
                                                                                          by=.(Carried_variants,Tile)], stringsAsFactors=F)
  
  cat("Sankaran_MPRA_Sig_subset_wide.dt.Zscore\n")
  cat(str(Sankaran_MPRA_Sig_subset_wide.dt.Zscore))
  cat("\n")
  
  #quit(status=1)
  
  Sankaran_MPRA_Sig_subset_wide<-merge(Sankaran_MPRA_Sig_subset_wide,
                                       Sankaran_MPRA_Sig_subset_wide.dt.Zscore,
                                       by=c("Carried_variants","Tile"),
                                       all.x=T)
  
  
  Sankaran_MPRA_Sig_subset_wide$VAR<-paste('chr',Sankaran_MPRA_Sig_subset_wide$Carried_variants,sep='')
  
  cat("Sankaran_MPRA_Sig_subset_wide_\n")
  cat(str(Sankaran_MPRA_Sig_subset_wide))
  cat("\n")
  
  
  check<-Sankaran_MPRA_Sig_subset_wide[which(Sankaran_MPRA_Sig_subset_wide$dbSNP == "rs737092" | 
                                               Sankaran_MPRA_Sig_subset_wide$dbSNP == "rs4490057" | 
                                               Sankaran_MPRA_Sig_subset_wide$dbSNP == "rs1546723"),]
  
  
  cat("check_\n")
  str(check)
  cat("\n")
  
  #### KEY_collpase_Plus_Variant_lineage_CLASSIFICATION ----
  
  
  KEY_collapse<-readRDS(file=opt$KEY_collpase_Plus_Variant_lineage_CLASSIFICATION)
  
  cat("KEY_collapse_0\n")
  cat(str(KEY_collapse))
  cat("\n")
  cat(sprintf(as.character(names(summary(KEY_collapse$Label_2)))))
  cat("\n")
  cat(sprintf(as.character(summary(KEY_collapse$Label_2))))
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
  
  #### HALF selection ----
  
  MPRA_TILES<-MPRA_Real_tile#[which(MPRA_Real_tile$TILE == "TILE_3"),]

  cat("MPRA_TILES_\n")
  cat(str(MPRA_TILES))
  cat("\n")

  cat(sprintf(as.character(names(summary(as.factor(MPRA_TILES$Tile))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(MPRA_TILES$Tile)))))
  cat("\n")
  
  MPRA_TILES_in_Sankaran<-MPRA_TILES[which(MPRA_TILES$carried_variants%in%Sankaran_MPRA_Sig_subset_wide$Carried_variants),]
  
  cat("MPRA_TILES_in_Sankaran_\n")
  cat(str(MPRA_TILES_in_Sankaran))
  cat("\n")
  
  cat(sprintf(as.character(names(summary(as.factor(MPRA_TILES_in_Sankaran$Label_2))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(MPRA_TILES_in_Sankaran$Label_2)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(MPRA_TILES_in_Sankaran$VAR))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(MPRA_TILES_in_Sankaran$VAR)))))
  cat("\n")
  
  
  MPRA_TILES_in_Sankaran_Vockley_REF<-MPRA_TILES_in_Sankaran[which(MPRA_TILES_in_Sankaran$variable == "Vockley_REF"),]
  
  cat("MPRA_TILES_in_Sankaran_Vockley_REF_\n")
  cat(str(MPRA_TILES_in_Sankaran_Vockley_REF))
  cat("\n")
  
  Sankaran_DEF<-Sankaran_MPRA_Sig_subset_wide[which(Sankaran_MPRA_Sig_subset_wide$Carried_variants%in%MPRA_TILES_in_Sankaran_Vockley_REF$carried_variants),]
  
  cat("Sankaran_DEF_\n")
  cat(str(Sankaran_DEF))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Sankaran_DEF$Carried_variants))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Sankaran_DEF$Carried_variants)))))
  cat("\n")
  
  
  Sankaran_DEF_subset<-unique(Sankaran_DEF[,c(which(colnames(Sankaran_DEF) == "VAR"),
                                              which(colnames(Sankaran_DEF) == "Tile"),
                                       which(colnames(Sankaran_DEF) == "Allelic_skew"))])
  
  
  cat("Sankaran_DEF_subset_\n")
  cat(str(Sankaran_DEF_subset))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Sankaran_DEF_subset$VAR))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Sankaran_DEF_subset$VAR)))))
  cat("\n")
  
  ##### merge Sankaran and MPRA ----
  
  
  DEF<-merge(MPRA_TILES_in_Sankaran_Vockley_REF,
             Sankaran_DEF_subset,
             by=c("VAR","Tile"))
  
  
  cat("DEF_\n")
  cat(str(DEF))
  cat("\n")
  
  cat(sprintf(as.character(names(summary(as.factor(DEF$Label_2))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(DEF$Label_2)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(DEF$VAR))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(DEF$VAR)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(DEF$DEF_CLASS))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(DEF$DEF_CLASS)))))
  cat("\n")
  
  DEF<-DEF[!is.na(DEF$DEF_CLASS),]

  DEF<-DEF[which(DEF$DEF_CLASS != "enhancer"),]

  DEF<-droplevels(DEF)

  cat("DEF_2\n")
  cat(str(DEF))
  cat("\n")

  cat(sprintf(as.character(names(summary(as.factor(DEF$Label_2))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(DEF$Label_2)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(DEF$VAR))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(DEF$VAR)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(DEF$DEF_CLASS))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(DEF$DEF_CLASS)))))
  cat("\n")
  
  # quit(status = 1)
  #### LOOP TO PRINT ----
  
  Cell_Type_array<-levels(DEF$Cell_Type)
  
  
  path5<-paste(out,'VJ_directionality','/', sep='')
  
  cat("path5\n")
  cat(sprintf(as.character(path5)))
  cat("\n")
  
  
  if (file.exists(path5)){
    
    
    
    
  } else {
    dir.create(file.path(path5))
    
  }
  
  
  
  colors<-c('#32A852','#1877C9','#553B68','#D45E85','#6DB2EE','#62D07F','#C9244B','#87447B','#D6B8E6')
  
  
  list_graphs<-list()
  list_LM<-list()
  
  for(iteration_Cell_Type_array in 1:length(Cell_Type_array))
  {
    Cell_Type_array_sel<-Cell_Type_array[iteration_Cell_Type_array]
    
    cat("--->:\n")
    cat(sprintf(as.character(Cell_Type_array_sel)))
    cat("\n")
    
    color_sel=colors[iteration_Cell_Type_array]
    
    
    
    
    DEF_2<-DEF[which(DEF$Cell_Type == Cell_Type_array_sel),]
    
    DEF_2<-droplevels(DEF_2)
    
    cat("DEF_2_\n")
    cat(str(DEF_2))
    cat("\n")
    
    linearModel_VJ <- summary(lm(DEF_2$Allelic_skew ~ DEF_2$value,
                                               data=DEF_2))
    
    cat("linearModel_VJ\n")
    cat(str(linearModel_VJ))
    cat("\n")

    r.squared_linearModel_VJ<-round(linearModel_VJ$r.squared,2)
    
    linearModel_VJ_coeffcient_df.m<-melt(linearModel_VJ$coefficients)
    
    colnames(linearModel_VJ_coeffcient_df.m)[which(colnames(linearModel_VJ_coeffcient_df.m)=="Var1")]<-"Terms"
    colnames(linearModel_VJ_coeffcient_df.m)[which(colnames(linearModel_VJ_coeffcient_df.m)=="Var2")]<-"Parameters"
    
    # cat("linearModel_VJ_coeffcient_df.m\n")
    # cat(str(linearModel_VJ_coeffcient_df.m))
    # cat("\n")
    
    slope<-round(linearModel_VJ_coeffcient_df.m$value[which(linearModel_VJ_coeffcient_df.m$Parameters == "Estimate" &
                                                                            linearModel_VJ_coeffcient_df.m$Terms == "DEF_2$value")],2)
    
    cat("slope\n")
    cat(str(slope))
    cat("\n")
    
    
    slope_pval<-round(-1*log10(linearModel_VJ_coeffcient_df.m$value[which(linearModel_VJ_coeffcient_df.m$Parameters == "Pr(>|t|)" &
                                                                                          linearModel_VJ_coeffcient_df.m$Terms == "DEF_2$value")]),2)
    
    cat("slope_pval\n")
    cat(str(slope_pval))
    cat("\n")
    
    intercept<-round(linearModel_VJ_coeffcient_df.m$value[which(linearModel_VJ_coeffcient_df.m$Parameters == "Estimate" &
                                                                                linearModel_VJ_coeffcient_df.m$Terms == "(Intercept)")],2)
    
    cat("intercept\n")
    cat(str(intercept))
    cat("\n")
    
    
    intercept_pval<-linearModel_VJ_coeffcient_df.m$value[which(linearModel_VJ_coeffcient_df.m$Parameters == "Pr(>|t|)" &
                                                                               linearModel_VJ_coeffcient_df.m$Terms == "(Intercept)")]
    
    cat("intercept_pval\n")
    cat(str(intercept_pval))
    cat("\n")
    
    
    
    
    ALT_ALT<-DEF_2[which(DEF_2$value < 1 &
                           DEF_2$Allelic_skew < 0),]
    
    cat("ALT_ALT\n")
    cat(str(ALT_ALT))
    cat("\n")
    
    ALT_REF<-DEF_2[which(DEF_2$value < 1 &
                           DEF_2$Allelic_skew > 0),]
    
    cat("ALT_REF\n")
    cat(str(ALT_REF))
    cat("\n")
    
    REF_REF<-DEF_2[which(DEF_2$value > 1 &
                           DEF_2$Allelic_skew > 0),]
    
    cat("REF_REF\n")
    cat(str(REF_REF))
    cat("\n")
    
    REF_ALT<-DEF_2[which(DEF_2$value > 1 &
                           DEF_2$Allelic_skew < 0),]
    
    cat("REF_ALT\n")
    cat(str(REF_ALT))
    cat("\n")
    
    TOTAL_CLASSIFIED<-sum(dim(ALT_ALT)[1],dim(ALT_REF)[1],dim(REF_REF)[1],dim(REF_ALT)[1])
    
    cat(sprintf(as.character(TOTAL_CLASSIFIED)))
    cat("\n")
    
    CLASSIFIED_RIGHT<-sum(dim(ALT_ALT)[1],dim(REF_REF)[1])
    
    cat(sprintf(as.character(CLASSIFIED_RIGHT)))
    cat("\n")
    
    CLASSIFIED_WRONG<-sum(dim(ALT_REF)[1],dim(REF_ALT)[1])
    
    cat(sprintf(as.character(CLASSIFIED_WRONG)))
    cat("\n")
    
    
    A.df<-as.data.frame(cbind(Cell_Type_array_sel,slope,intercept,dim(DEF_2)[1],r.squared_linearModel_VJ,slope_pval,CLASSIFIED_RIGHT,CLASSIFIED_WRONG,TOTAL_CLASSIFIED))
    
    # cat("A.df\n")
    # cat(str(A.df))
    # cat("\n")
    
    colnames(A.df)<-c("Cell_Type","LM_slope","LM_intercept","n","r.squared","logpval","CLASSIFIED_RIGHT","CLASSIFIED_WRONG","TOTAL_CLASSIFIED")
    
    cat("A.df\n")
    cat(str(A.df))
    cat("\n")
    
    
    
    
    breaks.x<-summary(DEF_2$value)
    
    cat(sprintf(as.character(names(breaks.x))))
    cat("\n")
    cat(sprintf(as.character(breaks.x)))
    cat("\n")
    
    
    breaks.x<-unique(sort(c(as.numeric(breaks.x),as.numeric(breaks.x[length(breaks.x)])+0.2)))
    labels.x<-as.character(round(breaks.x,2))
    cat(sprintf(as.character(labels.x)))
    cat("\n")
    
    
    breaks.y<-summary(DEF_2$Allelic_skew)
    
    cat(sprintf(as.character(names(breaks.y))))
    cat("\n")
    cat(sprintf(as.character(breaks.y)))
    cat("\n")
    
    breaks.y<-unique(sort(c(as.numeric(breaks.y),as.numeric(breaks.y[length(breaks.y)])+0.2)))
    labels.y<-as.character(round(breaks.y,2))
    cat(sprintf(as.character(labels.y)))
    cat("\n")
    
    # ,
    # shape=DEF_CLASS
    #geom_abline(slope= slope, intercept= intercept, size=2)+
    
    # geom_text(x=1, y=0.25, label=paste("y=",slope,"x"," + ",intercept,sep=''),size=6, family="sans", color=color_sel)+
    #   geom_text(x=1, y=0, label=paste("R^2=",r.squared_linearModel_VJ,sep=''),size=6, family="sans", color=color_sel)+
    #   geom_text(x=1, y=0.5, label=paste("logpval_LM=",slope_pval,sep=''),size=6, family="sans", color=color_sel)+
    
       
    graph<-ggplot(DEF_2, 
                  aes(x=value, 
                      y=Allelic_skew)) +
      geom_point(size=4, color=color_sel)+
      theme_bw()+
      scale_x_continuous(name="MPRA Vockley REF", breaks=breaks.x,labels=labels.x, limits=c(breaks.x[1],breaks.x[length(breaks.x)]))+
      scale_y_continuous(name="MPRA Sankaran Allelic Skew",breaks=breaks.y,labels=labels.y, limits=c(breaks.y[1],breaks.y[length(breaks.y)]))+
      theme(axis.title.y=element_text(size=24, family="sans"),
            axis.title.x=element_text(size=24, family="sans"),
            axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
            axis.text.x=element_text(angle=45,vjust=1,hjust=1,size=18, color="black", family="sans"),
            legend.title=element_text(size=16,color="black", family="sans"),
            legend.text=element_text(size=16,color="black", family="sans"))+
      geom_vline(color="red",linetype="dotted", xintercept=1)+
      geom_hline(color="red",linetype="dotted", yintercept=0)+
      ggeasy::easy_center_title()
    
    list_graphs[[iteration_Cell_Type_array]]<-graph
    list_LM[[iteration_Cell_Type_array]]<-A.df
    
    setwd(path5)
    
    
    
    svgname<-paste("Scatter_Vockley_REF_Allelic_Skew_",Cell_Type_array_sel,".svg",sep='')
    makesvg = TRUE
    
    if (makesvg == TRUE)
    {
      ggsave(svgname, plot= graph,
             device="svg",
             height=10, width=12)
    }
    
   
    
    
    
  }# iteration_Cell_Type_array
 
 
  #### PRINT ALL ----
  
  
  setwd(path5)
  
  cat("svg_graph\n")
  
 
  
  K562_graph<-list_graphs[[1]]
  CHRF_graph<-list_graphs[[2]]
  HL60_graph<-list_graphs[[3]]
  # THP1_graph<-list_graphs[[4]]
  
  
  graph_DEF<-plot_grid(K562_graph,CHRF_graph,HL60_graph,NULL,
                       nrow = 2,
                       ncol = 2,
                       labels=c("A) K-562","B) CHRF","C) HL-60","D) THP-1"),
                       align = "hv")
  
  svgname<-paste("Graph_VJ_directionality",".svg", sep='')
  makesvg = TRUE
  
  if (makesvg == TRUE)
  {
    
    ggsave(svgname, plot= graph_DEF,
           device="svg",
           height=10, width=12)
  }
  
  if(length(list_LM) >0)
  {
    LM_DEF = unique(as.data.frame(data.table::rbindlist(list_LM, fill = T)))
    
    
    cat("LM_DEF\n")
    cat(str(LM_DEF))
    cat("\n")
    
    write.table(LM_DEF,file="LM_object.tsv",sep="\t",quote=F,row.names = F)
    
    
  }
  
  
}

# ###########################################################################################################################################################################################
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
    make_option(c("--Sankaran_MPRA"), type="character", default=NULL, 
                metavar="Sankaran_MPRA", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--MPRA_Real_tile_QC2_PASS"), type="character", default=NULL, 
                metavar="MPRA_Real_tile_QC2_PASS", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--KEY_collpase_Plus_Variant_lineage_CLASSIFICATION"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="out", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "140__Rscript_v106.R
                        --subset type
                        --TranscriptEXP FILE.txt
                        --cadd FILE.txt
                        --ncboost FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  Generating_count_matrix_NSCS(opt)
  
}


###########################################################################

system.time( main() )
