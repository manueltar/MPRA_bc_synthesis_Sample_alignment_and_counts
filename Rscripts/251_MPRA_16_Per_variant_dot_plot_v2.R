

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

suppressMessages(library("cowplot", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))


opt = NULL

options(warn = 1)



dot_plot_graphs = function(option_list)
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
  
  #### READ and transform LONG_MATRIX ----
  
  LONG_MATRIX = read.table(opt$LONG_MATRIX, sep="\t", stringsAsFactors = F, header = T)
  
  
  colnames(LONG_MATRIX)[which(colnames(LONG_MATRIX) == "Real_Tile")]<-"REAL_TILE"
  
  # cat("----->LONG_MATRIX_0\n")
  # cat(str(LONG_MATRIX))
  # cat("\n")
  
  
  
  
  
  
  LONG_MATRIX$KEY<-gsub(";.+$","",LONG_MATRIX$REAL_TILE)
  LONG_MATRIX$TILE<-gsub("^[^;]+;","",LONG_MATRIX$REAL_TILE)
  LONG_MATRIX$TILE<-gsub(";.+$","",LONG_MATRIX$TILE)
  LONG_MATRIX$ALLELE_VERSION<-gsub("[^;]+;[^;]+;","",LONG_MATRIX$REAL_TILE)
  
  LONG_MATRIX$Tile<-"NA"
  
  LONG_MATRIX$Tile[which(LONG_MATRIX$TILE == "NINE")]<-"TILE_1"
  LONG_MATRIX$Tile[which(LONG_MATRIX$TILE == "TWO_THIRDS")]<-"TILE_2"
  LONG_MATRIX$Tile[which(LONG_MATRIX$TILE == "HALF")]<-"TILE_3"
  LONG_MATRIX$Tile[which(LONG_MATRIX$TILE == "ONE_THIRD")]<-"TILE_4"
  LONG_MATRIX$Tile[which(LONG_MATRIX$TILE == "A_TENTH")]<-"TILE_5"
  
  
  
  
  
  LONG_MATRIX$REAL_TILE_Plus_carried_variants<-paste(paste(LONG_MATRIX$KEY,LONG_MATRIX$Tile,sep="__"),LONG_MATRIX$Carried_variants,sep=";")
  
  
  
  
  # cat("----->LONG_MATRIX_PRE\n")
  # cat(str(LONG_MATRIX))
  # cat("\n")
  
  indx.int<-c(which(colnames(LONG_MATRIX) == "REAL_TILE_Plus_carried_variants"),which(colnames(LONG_MATRIX) == "VAR"),which(colnames(LONG_MATRIX) == "Carried_variants"),
              which(colnames(LONG_MATRIX) == "Tile"),which(colnames(LONG_MATRIX) == "chr"),which(colnames(LONG_MATRIX) == "start"),
              which(colnames(LONG_MATRIX) == "stop"))
  
  LONG_MATRIX_subset<-unique(LONG_MATRIX[,indx.int])
  
  # cat("LONG_MATRIX_subset_0\n")
  # cat(str(LONG_MATRIX_subset))
  # cat("\n")
  
  colnames(LONG_MATRIX_subset)<-c("REAL_TILE_Plus_carried_variants","VAR","carried_variants","TILE","chr","start","stop")
  
  LONG_MATRIX_subset$TILE<-factor(LONG_MATRIX_subset$TILE,
                                  levels=c("TILE_1","TILE_2","TILE_3","TILE_4","TILE_5"),
                                  ordered=T)
  
  cat("LONG_MATRIX_subset_1\n")
  cat(str(LONG_MATRIX_subset))
  cat("\n")

  #### MPRA Real Tile ----
  
  MPRA_Real_tile = readRDS(opt$MPRA_Real_tile_QC2_PASS)#, sep="\t", header = T), stringsAsFactors = F)
  
  cat("MPRA_Real_tile_0\n")
  cat(str(MPRA_Real_tile))
  cat("\n")
  cat(sprintf(as.character(names(summary(MPRA_Real_tile$DEF_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(MPRA_Real_tile$DEF_CLASS))))
  cat("\n")
  
  MPRA_Real_tile$DEF_CLASS<-as.character(MPRA_Real_tile$DEF_CLASS)
  
  # cat("MPRA_Real_tile_0.5\n")
  # cat(str(MPRA_Real_tile))
  # cat("\n")
  
  indx<-is.na(MPRA_Real_tile$DEF_CLASS)
  
  # cat("indx\n")
  # cat(str(indx))
  # cat("\n")
  
  
  check<-sum(indx)
  
  # cat("check\n")
  # cat(str(check))
  # cat("\n")
  
 
  
  MPRA_Real_tile$DEF_CLASS[indx]<-"NOT_ACTIVE"

  MPRA_Real_tile$DEF_CLASS<-factor(MPRA_Real_tile$DEF_CLASS,
                                  levels=c("NOT_ACTIVE","enhancer","ASE","E_Plus_ASE"),
                                  ordered=T)
  
 
  

  cat("MPRA_Real_tile_1\n")
  cat(str(MPRA_Real_tile))
  cat("\n")
  
  cat(sprintf(as.character(names(summary(as.factor(MPRA_Real_tile$Label))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(MPRA_Real_tile$Label)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(MPRA_Real_tile$DEF_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(MPRA_Real_tile$DEF_CLASS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(MPRA_Real_tile$QC2_CLASS))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(MPRA_Real_tile$QC2_CLASS)))))
  cat("\n")
  
# quit(status = 1)
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
  
  
  #### Read ALL_dB file ----
  
  ALL_dB<-as.data.frame(fread(file=opt$ALL_dB,sep="\t") , stringsAsFactors=F)
  
  cat("ALL_dB_0\n")
  cat(str(ALL_dB))
  cat("\n")
  
  
  indx.dep<-c(which(colnames(ALL_dB) == "maf_origin"))
  
  ALL_dB_subset<-unique(ALL_dB[,-indx.dep])
  
  # cat("ALL_dB_subset\n")
  # cat(str(ALL_dB_subset))
  # cat("\n")
  
  
  ALL_dB_subset$Allelic_Series_ID<-paste(ALL_dB_subset$phenotype,ALL_dB_subset$block_no,sep='__')
  
  cat("ALL_dB_subset_1\n")
  cat(str(ALL_dB_subset))
  cat("\n")
  
  
  indx.int<-c(which(colnames(ALL_dB_subset) == "VAR"),which(colnames(ALL_dB_subset) == "rs"))
  
  ALL_FINAL<-unique(ALL_dB_subset[,indx.int])
  
  cat("ALL_FINAL\n")
  cat(str(ALL_FINAL))
  cat("\n")
  cat(str(unique(ALL_FINAL$VAR)))
  cat("\n")
  
  
  #### LOOP VAR and carried variant ----
  
  CONDITION_BUG<-"NA"
  
  
  VARS<-unique(MPRA_Real_tile$VAR)
  
  # VARS<-"chr9_5079248_T_C"
  
  
  cat("VARS\n")
  cat(str(VARS))
  cat("\n")
  
  #### Get limits for logpval enhancer and ASE ----
  
 
  A<-summary(unique(MPRA_Real_tile$enhancer_logpval))
  
  
  cat("A\n")
  cat(sprintf(as.character(names(A))))
  cat("\n")
  cat(sprintf(as.character(A)))
  cat("\n")
  
  
  max_value<-A[6]
  min_value<-0
  
  step<-round((max_value-min_value)/3,0)
  
  
  cat("max_value:min_value:step\n")
  cat(sprintf(as.character(max_value)))
  cat("\t")
  cat(sprintf(as.character(min_value)))
  cat("\t")
  cat(sprintf(as.character(step)))
  cat("\n")
  
  
  
  
  
  if(step == 0)
  {
    step<-0.1
    breaks.enhancer_logpval<-sort(unique(c(1,2,5,seq(min_value,max_value+step, by=step))))
    labels.enhancer_logpval<-as.character(round(breaks.enhancer_logpval,2))
  }else{
    breaks.enhancer_logpval<-sort(unique(c(1,2,5,seq(min_value,max_value+step, by=step))))
    labels.enhancer_logpval<-as.character(round(breaks.enhancer_logpval,2))
    
  }
  
  breaks.size_enhancer_logpval<-breaks.enhancer_logpval[c(1:4,length(breaks.enhancer_logpval)-2,length(breaks.enhancer_logpval))]
  labels.size_enhancer_logpval<-as.character(round(breaks.size_enhancer_logpval,2))
  
  
  
  cat("breaks.enhancer_logpval\n")
  cat(sprintf(as.character(breaks.enhancer_logpval)))
  cat("\n")
  cat("labels.enhancer_logpval\n")
  cat(sprintf(as.character(labels.enhancer_logpval)))
  cat("\n")
  
  cat("labels.size_enhancer_logpval\n")
  cat(str(labels.size_enhancer_logpval))
  cat("\n")
  
  
  A<-summary(unique(MPRA_Real_tile$ASE_logpval))
  
  
  cat("A\n")
  cat(sprintf(as.character(names(A))))
  cat("\n")
  cat(sprintf(as.character(A)))
  cat("\n")
  
  
  max_value<-A[6]
  min_value<-0
  
  step<-round((max_value-min_value)/3,0)
  
  
  cat("max_value:min_value:step\n")
  cat(sprintf(as.character(max_value)))
  cat("\t")
  cat(sprintf(as.character(min_value)))
  cat("\t")
  cat(sprintf(as.character(step)))
  cat("\n")
  
  
  
  
  
  if(step == 0)
  {
    step<-0.1
    breaks.ASE_logpval<-sort(unique(c(1,2,5,seq(min_value,max_value+step, by=step))))
    labels.ASE_logpval<-as.character(round(breaks.ASE_logpval,2))
  }else{
    breaks.ASE_logpval<-sort(unique(c(1,2,5,seq(min_value,max_value+step, by=step))))
    labels.ASE_logpval<-as.character(round(breaks.ASE_logpval,2))
    
  }
  
  
  breaks.size_ASE_logpval<-breaks.ASE_logpval[c(1:4,length(breaks.ASE_logpval)-2,length(breaks.ASE_logpval))]

  labels.size_ASE_logpval<-as.character(round(breaks.size_ASE_logpval,2))
  
  
  
  cat("breaks.ASE_logpval\n")
  cat(sprintf(as.character(breaks.ASE_logpval)))
  cat("\n")
  cat("labels.ASE_logpval\n")
  cat(sprintf(as.character(labels.ASE_logpval)))
  cat("\n")
  
  cat("labels.size_ASE_logpval\n")
  cat(str(labels.size_ASE_logpval))
  cat("\n")
  
  # quit(status = 1)
  
  
  #### path2 ----

  path2<-paste(out2,'Per_variant_graphs','/', sep='')

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
    
    rsid_sel<-ALL_FINAL$rs[ALL_FINAL$VAR == VAR_sel]
    
    cat("---------------------------------->\t")
    cat(sprintf(as.character(VAR_sel)))
    cat("\t")
    cat(sprintf(as.character(rsid_sel)))
    cat("\n")
    
    if(VAR_sel == "chr9_5079248_T_C")
    {
       #chr9_5079248_T_C #chr2_171570079_NCGR chr5_1093511_G_A
      CONDITION_BUG<-1
    }else{
      
      CONDITION_BUG<-0
      
    }
    
    MPRA_Real_tile_sel<-MPRA_Real_tile[which(MPRA_Real_tile$VAR == VAR_sel),]
    
    if(CONDITION_BUG == 1)
    {
      cat("MPRA_Real_tile_sel\n")
      cat(str(MPRA_Real_tile_sel))
      cat("\n")
    }
    
    carried_variants_array<-unique(as.character(MPRA_Real_tile_sel$carried_variants))
    
    if(CONDITION_BUG == 1)
    {
      cat("carried_variants_array\n")
      cat(str(carried_variants_array))
      cat("\n")
    }
    
    Element_sel<-unique(as.character(MPRA_Real_tile_sel$KEY))
    Label_sel<-unique(as.character(MPRA_Real_tile_sel$Label))
    
   

    path3<-paste(out2,'Per_variant_graphs','/',paste(Element_sel,Label_sel,VAR_sel, sep='_'),'/', sep='')
    
    if(CONDITION_BUG == 1)
    {
      cat("path3\n")
      cat(sprintf(as.character(path3)))
      cat("\n")
    }
    
    if (file.exists(path3)){
      
    
      
      
    } else {
      dir.create(file.path(path3))
      
    }
    
    for(k in 1:length(carried_variants_array))
    {
     
      carried_variants_sel<-as.character(carried_variants_array[k])
      
      if(carried_variants_sel != "NCGR")
      {
        #### segment plot of tiles ----
        carried_variants_sel<-carried_variants_sel
        MPRA_Real_tile_carried_variants_sel<-MPRA_Real_tile_sel[which(MPRA_Real_tile_sel$carried_variants == carried_variants_sel),]
        LONG_MATRIX_subset_sel<-LONG_MATRIX_subset[which(LONG_MATRIX_subset$carried_variants%in%carried_variants_sel),]
        LONG_MATRIX_subset_sel$TILE<-factor(LONG_MATRIX_subset_sel$TILE,
                                            levels=rev(c("TILE_1","TILE_2","TILE_3","TILE_4","TILE_5")),
                                            ordered=T)
        
        LONG_MATRIX_subset_sel<-LONG_MATRIX_subset_sel[order(LONG_MATRIX_subset_sel$TILE),]
        
        max_POS<-max(LONG_MATRIX_subset_sel$stop)
        min_POS<-min(LONG_MATRIX_subset_sel$start)
        
        
        breaks.POS<-seq(min_POS-50,max_POS+50,by=50)
        labels.POS<-as.character(breaks.POS)
        color.POS<-rep("transparent", length(labels.POS))
        color.POS[1]<-"black"
        color.POS[length(color.POS)]<-"black"
        
        max_TILES<-max(as.numeric(LONG_MATRIX_subset_sel$TILE))
        min_TILES<-min(as.numeric(LONG_MATRIX_subset_sel$TILE))
        
        
        breaks.TILES<-seq(1,5,by=1)
        labels.TILES<-gsub("_"," ",as.character(rev(c("TILE_1","TILE_2","TILE_3","TILE_4","TILE_5"))))
        
        
        name_carried_variant<-sub("_",":",unique(LONG_MATRIX_subset_sel$carried_variants))
        name_carried_variant<-sub("_"," ",name_carried_variant)
        name_carried_variant<-sub("_",">",name_carried_variant)
        
        
        
        name_MPRA<-paste(VAR_sel,rsid_sel,name_carried_variant,sep="  ")
        name_chr<-unique(LONG_MATRIX_subset_sel$chr)
        name_interval<-paste("Genomic coordinates (GRCh37)",paste(name_chr,paste(labels.POS[1],labels.POS[length(labels.POS)],sep="-"), sep=":"),sep=' ')
        
        if(CONDITION_BUG == 1)
        {
          
          cat("-----carried_variants_sel------>\t")
          cat(sprintf(as.character(carried_variants_sel)))
          cat("\n")
          cat("MPRA_Real_tile_carried_variants_sel\n")
          cat(str(MPRA_Real_tile_carried_variants_sel))
          cat("\n")
          cat("LONG_MATRIX_subset_sel\n")
          cat(str(LONG_MATRIX_subset_sel))
          cat("\n")
          
          cat(sprintf(as.character(name_MPRA)))
          cat("\t")
          cat(sprintf(as.character(name_chr)))
          cat("\t")
         
          cat("labels.POS\n")
          
          cat(sprintf(as.character(labels.POS)))
          cat("\t")
          cat(sprintf(as.character(color.POS)))
          cat("\n")
          
          
          
          cat("breaks.TILES\n")
          cat(sprintf(as.character(breaks.TILES)))
          cat("\n")
          cat(sprintf(as.character(labels.TILES)))
          cat("\n")
        }
        
        if(sum(grep("\\|",carried_variants_sel))>0)
        {
          carried_SNPS<-unlist(strsplit(carried_variants_sel, split="\\|"))
          
          BK_SNPS<-carried_SNPS
          
          carried_SNPS<-gsub("^[^_]+_","",carried_SNPS)
          carried_SNPS<-as.numeric(gsub("_.+","",carried_SNPS))
          df_lines<-as.data.frame(cbind(BK_SNPS,
                                        carried_SNPS), stringsAsFactors=F)
          
          
        }else{
          
          carried_SNPS<-gsub("^[^_]+_","",carried_variants_sel)
          carried_SNPS<-as.numeric(gsub("_.+","",carried_SNPS))
          df_lines<-as.data.frame(cbind(rep(carried_variants_sel, length(carried_SNPS)),
                                        carried_SNPS), stringsAsFactors=F)
          
        }
        
        
        
        if(CONDITION_BUG == 1)
        {
          
          cat("carried_SNPS_UNKNOWN:\t")
          cat(sprintf(as.character(carried_SNPS)))
          cat("\n")
        }
        
       
        
        TILE_selected<-rev(breaks.TILES)
        
        TILE_selected<-TILE_selected[c(1:dim(df_lines)[1])]
        
        ypos<-max(TILE_selected)
        
        colnames(df_lines)<-c("carried_variants","pos")
        df_lines$ypos<-ypos
        df_lines$ypos_offset<-ypos+0.5
        
        df_lines$xpos<-as.numeric(df_lines$pos)
        
        row_odd <- seq_len(nrow(df_lines)) %% 2
        
        df_lines$hjust_parameter<-0.1
        
        df_lines$hjust_parameter[which(row_odd == 1)]<-df_lines$hjust_parameter[which(row_odd == 1)]*-1
        
        
        df_lines$carried_variants<-sub("_",":",df_lines$carried_variants)
        df_lines$carried_variants<-sub("_"," ",df_lines$carried_variants)
        df_lines$carried_variants<-sub("_",">",df_lines$carried_variants)
        
        if(CONDITION_BUG == 1)
        {
          cat("df_lines_UNKNOWN\n")
          cat(str(df_lines))
          cat("\n")
          
          cat("row_odd\n")
          cat(str(row_odd))
          cat("\n")
          
          
          # quit(status = 1)
        }
        
      }else{
        
        carried_variants_sel<-VAR_sel
        MPRA_Real_tile_carried_variants_sel<-MPRA_Real_tile_sel[which(MPRA_Real_tile_sel$VAR == carried_variants_sel),]
        LONG_MATRIX_subset_sel<-LONG_MATRIX_subset[which(LONG_MATRIX_subset$VAR%in%carried_variants_sel),]
        LONG_MATRIX_subset_sel$TILE<-factor(LONG_MATRIX_subset_sel$TILE,
                                            levels=rev(c("TILE_1","TILE_2","TILE_3","TILE_4","TILE_5")),
                                            ordered=T)
        
        LONG_MATRIX_subset_sel<-LONG_MATRIX_subset_sel[order(LONG_MATRIX_subset_sel$TILE),]
        
        max_POS<-max(LONG_MATRIX_subset_sel$stop)
        min_POS<-min(LONG_MATRIX_subset_sel$start)
        
        
        breaks.POS<-seq(min_POS-50,max_POS+50,by=50)
        labels.POS<-as.character(breaks.POS)
        color.POS<-rep("transparent", length(labels.POS))
        color.POS[1]<-"black"
        color.POS[length(color.POS)]<-"black"
        
        max_TILES<-max(as.numeric(LONG_MATRIX_subset_sel$TILE))
        min_TILES<-min(as.numeric(LONG_MATRIX_subset_sel$TILE))
        
        
        breaks.TILES<-seq(1,5,by=1)
        labels.TILES<-gsub("_"," ", as.character(rev(c("TILE_1","TILE_2","TILE_3","TILE_4","TILE_5"))))
        
        name_carried_variant<-sub("_",":",unique(LONG_MATRIX_subset_sel$carried_variants))
        name_carried_variant<-sub("_"," ",name_carried_variant)
        name_carried_variant<-sub("_",">",name_carried_variant)
        
        
        
        name_MPRA<-paste(VAR_sel,rsid_sel,name_carried_variant,sep="  ")
        name_chr<-unique(LONG_MATRIX_subset_sel$chr)
        name_interval<-paste("Genomic coordinates (GRCh37)",paste(name_chr,paste(labels.POS[1],labels.POS[length(labels.POS)],sep="-"), sep=":"),sep=' ')
        
        if(CONDITION_BUG == 1)
        {
          
          cat("-----carried_variants_sel------>\t")
          cat(sprintf(as.character(carried_variants_sel)))
          cat("\n")
          cat("MPRA_Real_tile_carried_variants_sel\n")
          cat(str(MPRA_Real_tile_carried_variants_sel))
          cat("\n")
          cat("LONG_MATRIX_subset_sel\n")
          cat(str(LONG_MATRIX_subset_sel))
          cat("\n")
          
          cat(sprintf(as.character(name_MPRA)))
          cat("\t")
          cat(sprintf(as.character(name_chr)))
          cat("\t")
          
          
          cat(sprintf(as.character(labels.POS)))
          cat("\n")
          
          cat("breaks.TILES\n")
          cat(sprintf(as.character(breaks.TILES)))
          cat("\n")
          cat(sprintf(as.character(labels.TILES)))
          cat("\n")
        }
        
        
        
        carried_SNPS<-gsub("^[^_]+_","",carried_variants_sel)
        carried_SNPS<-as.numeric(gsub("_.+","",carried_SNPS))
        
        if(CONDITION_BUG == 1)
        {
          
          cat("carried_SNPS_NCGR:\t")
          cat(sprintf(as.character(carried_SNPS)))
          cat("\n")
        }
        
        df_lines<-as.data.frame(cbind(rep(carried_variants_sel, length(carried_SNPS)),
                                      carried_SNPS), stringsAsFactors=F)
        
        TILE_selected<-rev(breaks.TILES)
        
        TILE_selected<-TILE_selected[c(1:dim(df_lines)[1])]
        
        ypos<-max(TILE_selected)
        
        colnames(df_lines)<-c("carried_variants","pos")
        df_lines$ypos<-ypos
        df_lines$ypos_offset<-ypos+0.5
        
        df_lines$xpos<-as.numeric(df_lines$pos)
        
        row_odd <- seq_len(nrow(df_lines)) %% 2
        
        df_lines$hjust_parameter<-0.1
        
        df_lines$hjust_parameter[which(row_odd == 1)]<-df_lines$hjust_parameter[which(row_odd == 1)]*-1
        
        
        df_lines$carried_variants<-sub("_",":",df_lines$carried_variants)
        df_lines$carried_variants<-sub("_"," ",df_lines$carried_variants)
        df_lines$carried_variants<-sub("_",">",df_lines$carried_variants)
        
        if(CONDITION_BUG == 1)
        {
          cat("df_lines_NCGR\n")
          cat(str(df_lines))
          cat("\n")
        }
        
      }
     
      # geom_segment(color="steelblue",size=8)+
      #yend=as.numeric(TILE),
      # ,      xend=stop
      
      # element_blank(),
      
      
      track_plot<-ggplot(data=LONG_MATRIX_subset_sel,
                         aes(y=as.numeric(TILE),
                             x=start))+
        geom_rect(data=LONG_MATRIX_subset_sel,
                  aes(NULL,NULL,
                      xmin=start,
                      xmax=stop,
                      ymin=as.numeric(TILE),
                      ymax=as.numeric(TILE)+0.1), fill="white", color="black")+
        theme_classic()+
        scale_y_continuous(name=NULL,
                           breaks=breaks.TILES,labels=labels.TILES,
                           limits=c(0.65,breaks.TILES[length(breaks.TILES)]+0.5))+
        scale_x_continuous(name=name_interval,
                           breaks=breaks.POS,labels=labels.POS,
                           limits=c(breaks.POS[1],breaks.POS[length(breaks.POS)]))+
        theme(axis.title.x=element_blank(),
              axis.text.y=element_text(angle=0,size=10, color="black", family="sans"),
              axis.text.x=element_text(angle=45,size=6,hjust=1,vjust=1, color=color.POS, family="sans"))+
        ggeasy::easy_center_title()
        
        
        
        # 
      
      if(CONDITION_BUG == 1)
      {
        cat("Graph track_plot\n") 
        
        cat("df_lines_PREPLOT\n")
        cat(str(df_lines))
        cat("\n")
        
        # 
        # breaks.TILES[1]
      }
      
      
      track_plot<-track_plot+
        geom_segment(data=df_lines,
                     aes(x=xpos,
                         xend=xpos,
                         y=0.65,
                         yend=ypos_offset), linetype=2)
      
      if(CONDITION_BUG == 1)
      {
        cat("Graph with segments addded\n") 
      }
      
      
      
      # track_plot<-track_plot+
      #   geom_text(data=df_lines,
      #             aes(y=ypos_offset,
      #                 x=xpos,
      #                 label=carried_variants,
      #                 family="sans",
      #                 hjust=hjust_parameter),vjust=1,size=5)
      
      if(CONDITION_BUG == 1)
      {
        cat("Graph with carried variant text added\n") 
      }
      
      
    
      
      # setwd(out2)
      # 
      # #svglite(paste('TEST_track_plot_',carried_variants_sel,'.svg',sep=''), width = 8, height = 8)
      # svglite(paste('TEST_track_plot','.svg',sep=''), width = 8, height = 8)
      # print(track_plot)
      # dev.off()
      # #
      # quit(status = 1)
     
      ########################## LogFC/logpval PLOTS --------------
      
     
      
      
      selected_variables<-c("LogFC")
      
      LogFC_df<-MPRA_Real_tile_carried_variants_sel[which(MPRA_Real_tile_carried_variants_sel$variable%in%selected_variables),]
      
      levels_TILE<-levels(LogFC_df$TILE)
      
      LogFC_df$TILE<-factor(LogFC_df$TILE,
                            levels = rev(levels_TILE),
                            ordered=T)
      
      LogFC_df$REP_CLASS<-"NA"
      
      
      LogFC_df$REP_CLASS[which(LogFC_df$DEF_CLASS == "enhancer")]<-"ACTIVE"
      LogFC_df$REP_CLASS[which(LogFC_df$DEF_CLASS == "E_Plus_ASE")]<-"ACTIVE"
      
      LogFC_df$REP_CLASS[which(LogFC_df$DEF_CLASS == "NOT_ACTIVE")]<-"NOT_ACTIVE"
      LogFC_df$REP_CLASS[which(LogFC_df$DEF_CLASS == "ASE")]<-"NOT_ACTIVE"
      
      LogFC_df$REP_CLASS<-factor(LogFC_df$REP_CLASS,
                            levels =c("NOT_ACTIVE","ACTIVE"),
                            ordered=T)
      
      
      # cat("LogFC_df_0")
      # str(LogFC_df)
      # cat("\n")
      
      if(CONDITION_BUG == 1)
      {
        cat("LogFC_df_0")
        str(LogFC_df)
        cat("\n")
      }
      
      A<-summary(LogFC_df$value)
      
      
      if(CONDITION_BUG == 1)
      {
        cat("A\n")
        cat(sprintf(as.character(names(A))))
        cat("\n")
        cat(sprintf(as.character(A)))
        cat("\n")
      }
      
      max_log_value<-A[6]
      min_log_value<-A[1]
      
      step<-round((max_log_value-min_log_value)/2,2)
      
      if(CONDITION_BUG == 1)
      {
      
        cat("max_log_value:min_log_value:step\n")
        cat(sprintf(as.character(max_log_value)))
        cat("\t")
        cat(sprintf(as.character(min_log_value)))
        cat("\t")
        cat(sprintf(as.character(step)))
        cat("\n")
      }
      
      
      if(step == 0)
      {
        step<-0.1
       
      }else{
       
      }
      
      
      #min_log_value
      breaks.LogFC<-sort(unique(c(0,seq(min_log_value,max_log_value+step, by=step))))
     # breaks.LogFC<-c(which())
      
      labels.LogFC<-as.character(round(breaks.LogFC,2))
      
      if(CONDITION_BUG == 1)
      {
        cat("labels.LogFC\n")
        cat(sprintf(as.character(labels.LogFC)))
        cat("\n")
        
        # quit(status = 1)
      }
      
      
      dotplot_Log_FC<-ggplot(data=LogFC_df,
                      aes(y=TILE,
                          x=Cell_Type)) +
        geom_point(aes(fill=value,color=REP_CLASS,size=enhancer_logpval), stroke=1, shape=21)+
        scale_fill_gradient2(
          low = "gray", 
          mid = "white", 
          high = "green", 
          midpoint = .00,
          breaks=breaks.LogFC,labels=labels.LogFC,
          limits=c(breaks.LogFC[1]-0.1,breaks.LogFC[length(breaks.LogFC)] +0.1),name='LogFC',na.value = "gray")+
        scale_color_manual(values=c("gray","black"),name=paste('Significant','Activity',sep="\n"), drop=F)+
        scale_size(name=paste('-log10pval','Enhancer activity',sep="\n"), range = c(0,8), breaks=breaks.size_enhancer_logpval, 
                   labels=labels.size_enhancer_logpval, limits=c(breaks.size_enhancer_logpval[1],breaks.size_enhancer_logpval[length(breaks.size_enhancer_logpval)])) +
        scale_y_discrete(name=NULL, drop=F)+
        scale_x_discrete(name=NULL, drop=F)+
        theme_classic()+
        theme(axis.text.y=element_blank(),
              axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=8, color="black", family="sans"))+
        theme(legend.position="hidden",legend.title=element_text(size=6, color="black", family="sans", face="bold"), legend.text = element_text(size=5, color="black", family="sans", face="plain"))+
        ggeasy::easy_center_title()
      
      dotplot_Log_FC<-dotplot_Log_FC+ 
                      guides(col = guide_legend(order = 2), size = guide_legend(order = 3))
      
      # element_text(angle=0,size=8, color="black", family="sans")
      
      # setwd(out2)
      # 
      # svglite(paste('TEST_LogFc_dot_plot','.svg',sep=''), width = 8, height = 8)
      # print(dotplot_Log_FC)
      # dev.off()
      
      
      if(CONDITION_BUG == 1)
      {
        cat("dotplot_Log_FC DONE\n")
      }
      
      # if(CONDITION_BUG == 1)
      # {
      #   quit(status = 1)
      # }
     
      ########################## ASE/logpval PLOTS --------------
      
      
      
      
      selected_variables<-c("ASE")
      
      ASE_df<-MPRA_Real_tile_carried_variants_sel[which(MPRA_Real_tile_carried_variants_sel$variable%in%selected_variables),]
      
      levels_TILE<-levels(ASE_df$TILE)
      
      ASE_df$TILE<-factor(ASE_df$TILE,
                            levels = rev(levels_TILE),
                            ordered=T)
      
      ASE_df$REP_CLASS<-"NA"
      
      
      ASE_df$REP_CLASS[which(ASE_df$DEF_CLASS == "enhancer")]<-"NOT_ACTIVE"
      ASE_df$REP_CLASS[which(ASE_df$DEF_CLASS == "E_Plus_ASE")]<-"ACTIVE"
      
      ASE_df$REP_CLASS[which(ASE_df$DEF_CLASS == "NOT_ACTIVE")]<-"NOT_ACTIVE"
      ASE_df$REP_CLASS[which(ASE_df$DEF_CLASS == "ASE")]<-"ACTIVE"
      
      ASE_df$REP_CLASS<-factor(ASE_df$REP_CLASS,
                                 levels =c("NOT_ACTIVE","ACTIVE"),
                                 ordered=T)
      
      
      # cat("ASE_df_0")
      # str(ASE_df)
      # cat("\n")
      
      if(CONDITION_BUG == 1)
      {
        cat("ASE_df_0")
        str(ASE_df)
        cat("\n")
        
        # quit(status = 1)
      }
      
      A<-summary(ASE_df$value)
      
      
      if(CONDITION_BUG == 1)
      {
        cat("A\n")
        cat(sprintf(as.character(names(A))))
        cat("\n")
        cat(sprintf(as.character(A)))
        cat("\n")
      }
      
      max_value<-A[6]
      min_value<-A[1]
      
      step<-round((max_value-min_value)/2,2)
      
      if(CONDITION_BUG == 1)
      {
        
        cat("max_value:min_value:step\n")
        cat(sprintf(as.character(max_value)))
        cat("\t")
        cat(sprintf(as.character(min_value)))
        cat("\t")
        cat(sprintf(as.character(step)))
        cat("\n")
      }
      
      
      if(step == 0)
      {
        step<-0.1
        
      }else{
        
      }
      
      
      #min_log_value
      breaks.ASE<-sort(unique(c(1,seq(min_value,max_value+step, by=step))))
      
      
      labels.ASE<-as.character(round(breaks.ASE,2))
      
      if(CONDITION_BUG == 1)
      {
        cat("labels.ASE\n")
        cat(sprintf(as.character(labels.ASE)))
        cat("\n")
        
        # quit(status = 1)
      }
      
      
      dotplot_ASE<-ggplot(data=ASE_df,
                             aes(y=TILE,
                                 x=Cell_Type)) +
        geom_point(aes(color=REP_CLASS,fill=value,size=ASE_logpval), stroke=1, shape=21)+
        theme_classic()+
        scale_fill_gradient2(
          low = "blue", 
          mid = "white", 
          high = "red", 
          midpoint = 1,
          breaks=breaks.ASE,labels=labels.ASE,
          limits=c(breaks.ASE[1]-0.1,breaks.ASE[length(breaks.ASE)] +0.1),name='ASE',na.value = "gray")+
        scale_color_manual(values=c("gray","black"),name=paste('Significant','Activity', sep="\n"), drop=F)+
        scale_size(range = c(0,8), breaks=breaks.size_ASE_logpval, labels=labels.size_ASE_logpval, limits=c(breaks.size_ASE_logpval[1],breaks.size_ASE_logpval[length(breaks.size_ASE_logpval)]), 
                   name=paste('-log10pval','ASE activity', sep="\n")) +
        scale_y_discrete(name=NULL, drop=F)+
        scale_x_discrete(name=NULL, drop=F)+
        theme(axis.text.y=element_blank(),
              axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=8, color="black", family="sans"))+
        theme(legend.position="hidden",legend.title=element_text(size=6, color="black", family="sans", face="bold"), legend.text = element_text(size=5, color="black", family="sans", face="plain"))+
        ggeasy::easy_center_title()
      
      dotplot_ASE<-dotplot_ASE+ 
        guides(col = guide_legend(order = 2), size = guide_legend(order = 3))
      
     
        
      # setwd(out2)
      # 
      # svglite(paste('TEST_ASE_dot_plot','.svg',sep=''), width = 8, height = 8)
      # print(dotplot_ASE)
      # dev.off()
      
      if(CONDITION_BUG == 1)
      {
        cat("DONE dotplot_ASE\n")
       # quit(status = 1) 
      }
      
      ##### Graph DEF ----
      
      
      
      dotplot_ASE<-dotplot_ASE+ 
        theme(legend.position="right",legend.title=element_text(size=6, color="black", family="sans", face="bold"), legend.text = element_text(size=5, color="black", family="sans", face="plain"))
      
      dotplot_Log_FC<-dotplot_Log_FC+ 
        theme(legend.position="right",legend.title=element_text(size=6, color="black", family="sans", face="bold"), legend.text = element_text(size=5, color="black", family="sans", face="plain"))
      
      graph_DEF_1<-plot_grid(track_plot,dotplot_Log_FC,dotplot_ASE,
                             nrow = 1,
                             ncol = 3,
                             rel_widths=c(1,1,1),
                             align = "h")
      
      path4<-paste(out2,'Per_variant_graphs','/',paste(Element_sel,Label_sel,VAR_sel, sep='_'),'/',carried_variants_sel,'/', sep='')
      
      cat("path4\n")
      cat(sprintf(as.character(path4)))
      cat("\n")
      
      if (file.exists(path4)){
        
        
        
        
      } else {
        dir.create(file.path(path4))
        
      }
      
      setwd(path4)
      
      svglite(paste('MPRA_Dot_plot_combined_ALL','.svg',sep=''), width = 8, height = 8)
      print(graph_DEF_1)
      dev.off()
      
      saveRDS(graph_DEF_1,paste('MPRA_Dot_plot_combined_ALL_',carried_variants_sel,'.rds',sep=''))
      
      path_variant_interpretation<-paste("/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/MPRA/",VAR_sel,'/',sep='')
      
      
      if (file.exists(path_variant_interpretation)){
        
        
        
        
      } else {
        dir.create(file.path(path_variant_interpretation))
        
      }
      
      setwd(path_variant_interpretation)
      
      
      
      svglite(paste('MPRA_Dot_plot_combined_ALL',carried_variants_sel,'.svg',sep=''), width = 8, height = 8)
      print(graph_DEF_1)
      dev.off()
      
      saveRDS(graph_DEF_1,paste('MPRA_Dot_plot_combined_ALL_',carried_variants_sel,'.rds',sep=''))
      
      
      
    }#k carried_variants_array[k]
  }# i VARS VARS[i]
  
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
    make_option(c("--LONG_MATRIX"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ALL_dB"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="out", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out2"), type="character", default=NULL, 
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

  dot_plot_graphs(opt)
  
}


###########################################################################

system.time( main() )
