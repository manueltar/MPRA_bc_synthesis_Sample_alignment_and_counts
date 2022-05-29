

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
suppressMessages(library("BiocGenerics", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("S4Vectors", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("IRanges", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("GenomeInfoDb", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("GenomicRanges", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("Biobase", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("AnnotationDbi", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("GenomicFeatures", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("OrganismDbi", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("GO.db", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("org.Hs.eg.db", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("Homo.sapiens", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("gwascat", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("rtracklayer", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("liftOver",lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))

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
  
  #### READ and transform LONG_MATRIX ----
  
  LONG_MATRIX = read.table(opt$LONG_MATRIX, sep="\t", stringsAsFactors = F, header = T)
  
  
  colnames(LONG_MATRIX)[which(colnames(LONG_MATRIX) == "Real_Tile")]<-"REAL_TILE"
  
  cat("----->LONG_MATRIX_0\n")
  cat(str(LONG_MATRIX))
  cat("\n")
  
  
  
  
  
  
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
  
  
  
  
  cat("----->LONG_MATRIX_PRE\n")
  cat(str(LONG_MATRIX))
  cat("\n")
  
  indx.int<-c(which(colnames(LONG_MATRIX) == "REAL_TILE_Plus_carried_variants"),which(colnames(LONG_MATRIX) == "VAR"),which(colnames(LONG_MATRIX) == "Carried_variants"),
              which(colnames(LONG_MATRIX) == "Tile"),which(colnames(LONG_MATRIX) == "chr"),which(colnames(LONG_MATRIX) == "start"),
              which(colnames(LONG_MATRIX) == "stop"))
  
  LONG_MATRIX_subset<-unique(LONG_MATRIX[,indx.int])
  
  cat("LONG_MATRIX_subset_0\n")
  cat(str(LONG_MATRIX_subset))
  cat("\n")
  
  colnames(LONG_MATRIX_subset)<-c("REAL_TILE_Plus_carried_variants","VAR","carried_variants","TILE","chr","start","stop")
  
  cat("LONG_MATRIX_subset_1\n")
  cat(str(LONG_MATRIX_subset))
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
  
  #### X and Y GLOBAL points ----
  
  MPRA_Real_tile_wide<-as.data.frame(pivot_wider(MPRA_Real_tile,
                            names_from=variable,
                            values_from=value), stringsAsFactors=F)
  
  cat("MPRA_Real_tile_wide\n")
  cat(str(MPRA_Real_tile_wide))
  cat("\n")
  
  
  Vockley_REF<-summary(MPRA_Real_tile_wide$Vockley_REF, na.rm=T)
  
  cat("Vockley_REF\n")
  cat(sprintf(as.character(names(Vockley_REF))))
  cat("\n")
  cat(sprintf(as.character(Vockley_REF)))
  cat("\n")
  
  # max_Vockley_REF<-Vockley_REF[6]
  # min_Vockley_REF<-Vockley_REF[1]
  min_Vockley_REF<-0.4
  max_Vockley_REF<-3
  
  step<-(max_Vockley_REF-min_Vockley_REF)/5
  
  if(step ==0)
  {
    
    step<-1
  }
  
  breaks.Vockley_REF<-sort(c(1,seq(min_Vockley_REF,max_Vockley_REF,by=step)))
  labels_Vockley_REF<-as.numeric(round(breaks.Vockley_REF,2))
  
  
  cat("labels_Vockley_REF\n")
  cat(sprintf(as.character(labels_Vockley_REF)))
  cat("\n")
  
 
  LogFC<-summary(MPRA_Real_tile_wide$LogFC, na.rm=T)
  
  cat("LogFC\n")
  cat(sprintf(as.character(names(LogFC))))
  cat("\n")
  cat(sprintf(as.character(LogFC)))
  cat("\n")
  
  max_log_LogFC<-LogFC[6]
  min_log_LogFC<--1.5
  
  #min_log_LogFC
  breaks.log_LogFC<-sort(c(0,seq(min_log_LogFC,max_log_LogFC+0.5, by=0.5)))
  
  
  labels_LogFC<-as.numeric(round(logratio2foldchange(breaks.log_LogFC, base=2),2))
  
  
  cat("labels_LogFC\n")
  cat(sprintf(as.character(labels_LogFC)))
  cat("\n")
  
  
 
  list_BED_LogFC<-list()
  list_BED_Vockley_REF<-list()
  
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
      
      #### LONG MATRIX merge
      
      MPRA_Real_tile_carried_variants_sel_wide<-merge(MPRA_Real_tile_carried_variants_sel_wide,
                                                      LONG_MATRIX_subset,
                                                      by=colnames(LONG_MATRIX_subset)[which(colnames(LONG_MATRIX_subset)%in%colnames(MPRA_Real_tile_carried_variants_sel_wide))])
      
      cat("MPRA_Real_tile_carried_variants_sel_wide_2\n")
      cat(str(MPRA_Real_tile_carried_variants_sel_wide))
      cat("\n")
      
      MPRA_Real_tile_carried_variants_sel_wide<-MPRA_Real_tile_carried_variants_sel_wide[order(MPRA_Real_tile_carried_variants_sel_wide$Cell_Type,
                                                                                               MPRA_Real_tile_carried_variants_sel_wide$TILE),]
      
      
      cat("MPRA_Real_tile_carried_variants_sel_wide_3\n")
      cat(str(MPRA_Real_tile_carried_variants_sel_wide))
      cat("\n")
      
      POS_array<-c(MPRA_Real_tile_carried_variants_sel_wide$start,MPRA_Real_tile_carried_variants_sel_wide$stop)
      
      min_pos<-min(POS_array)-100
      max_pos<-max(POS_array)+100
      
      chrom2 = unique(as.character(MPRA_Real_tile_carried_variants_sel_wide$chr))
      chromstart2 = as.numeric(min_pos)
      chromend2 = as.numeric(max_pos)
      
      
      
      
     
      
      BED<- data.frame(matrix(vector(), dim(MPRA_Real_tile_carried_variants_sel_wide)[1], 7,
                              dimnames=list(c(),
                                            c("chrom","start","end","name","LogFC","strand","type"))),
                       stringsAsFactors=F)
      
      
      BED$chrom<-as.character(MPRA_Real_tile_carried_variants_sel_wide$chr)
      BED$start<-MPRA_Real_tile_carried_variants_sel_wide$start
      BED$end<-MPRA_Real_tile_carried_variants_sel_wide$stop
      BED$name<-paste(MPRA_Real_tile_carried_variants_sel_wide$REAL_TILE_Plus_carried_variants,MPRA_Real_tile_carried_variants_sel_wide$Cell_Type,MPRA_Real_tile_carried_variants_sel_wide$DEF_CLASS,sep="|")
      BED$LogFC<-MPRA_Real_tile_carried_variants_sel_wide$LogFC
      BED$strand<-'1'
      BED$type<-'utr'
      
      cat("BED\n")
      cat(str(BED))
      cat("\n")
      
      
     
      
      
      setwd(path4)
      
      svglite(paste('SUSHI','_FC_','.svg',sep=''), width = 8, height = 8)
      
      pg = plotGenes(BED[order(BED$start, decreasing=F),],chrom2,chromstart2,chromend2,
                     types = BED$type,
                     colorby=round(logratio2foldchange(BED$LogFC, base=2),1),
                     colorbycol= SushiColors(7),colorbyrange=c(labels_LogFC[1],labels_LogFC[length(labels_LogFC)]),
                     labeltext=TRUE,maxrows=50,height=0.4,plotgenetype="box")
      
      
      labelgenome( chrom2, chromstart2,chromend2,n=10,scale="bp")
      addlegend(pg[[1]],palette=pg[[2]],title=paste("FC"),side="right",
                bottominset=0.4,topinset=0,xoffset=-.035,labelside="left",
                width=0.025,title.offset=0.055)

      dev.off()
      
      list_BED_LogFC[[i]]<-BED
      
      #### Vockley _REF
      
      BED<- data.frame(matrix(vector(), dim(MPRA_Real_tile_carried_variants_sel_wide)[1], 7,
                              dimnames=list(c(),
                                            c("chrom","start","end","name","Vockley_REF","strand","type"))),
                       stringsAsFactors=F)
      
      
      BED$chrom<-as.character(MPRA_Real_tile_carried_variants_sel_wide$chr)
      BED$start<-MPRA_Real_tile_carried_variants_sel_wide$start
      BED$end<-MPRA_Real_tile_carried_variants_sel_wide$stop
      BED$name<-paste(MPRA_Real_tile_carried_variants_sel_wide$REAL_TILE_Plus_carried_variants,MPRA_Real_tile_carried_variants_sel_wide$Cell_Type,MPRA_Real_tile_carried_variants_sel_wide$DEF_CLASS,sep="|")
      BED$Vockley_REF<-MPRA_Real_tile_carried_variants_sel_wide$Vockley_REF
      BED$strand<-'1'
      BED$type<-'utr'
      
      cat("BED\n")
      cat(str(BED))
      cat("\n")
      
      setwd(path4)
      
      svglite(paste('SUSHI','_Vockley_REF_','.svg',sep=''), width = 8, height = 8)
      
      pg = plotGenes(BED[order(BED$start, decreasing=F),],chrom2,chromstart2,chromend2,
                     types = BED$type,
                     colorby=round(BED$Vockley_REF,2),
                     colorbycol= SushiColors(7),colorbyrange=c(labels_Vockley_REF[1],labels_Vockley_REF[length(labels_Vockley_REF)]),
                     labeltext=TRUE,maxrows=50,height=0.4,plotgenetype="box")
      
      
      labelgenome( chrom2, chromstart2,chromend2,n=10,scale="bp")
      addlegend(pg[[1]],palette=pg[[2]],title=paste("FC"),side="right",
                bottominset=0.4,topinset=0,xoffset=-.035,labelside="left",
                width=0.025,title.offset=0.055)
      
      dev.off()
      
      list_BED_Vockley_REF[[i]]<-BED
     
    }#k carried_variants_array[k]
  }# i VARS VARS[i]
  
  
  if(length(list_BED_Vockley_REF) >0)
  {
    
    
    BED_Vockley_REF = unique(as.data.frame(data.table::rbindlist(list_BED_Vockley_REF, fill = T)))
    
    
    cat("BED_Vockley_REF_0\n")
    cat(str(BED_Vockley_REF))
    cat("\n")
    
    
    BED_Vockley_REF$chr<-factor(BED_Vockley_REF$chr,
                   levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                            "chr17","chr18","chr19","chr20","chr21","chr22","chr23","chrX","chrY"),
                   ordered=T)
    
    BED_Vockley_REF<-BED_Vockley_REF[order(BED_Vockley_REF$chr,BED_Vockley_REF$start, decreasing = F),]
    
    ### export bed of tiles GRCh37----
    

    cat("BED_Vockley_REF_1\n")
    cat(str(BED_Vockley_REF))
    cat("\n")
    
    
    gr_bed_files <- GRanges(
      seqnames = BED_Vockley_REF$chr,
      ranges=IRanges(
        start=BED_Vockley_REF$start,
        end=BED_Vockley_REF$end,
        name=BED_Vockley_REF$name,
        score=BED_Vockley_REF$Vockley_REF
        ))
    
    gr_bed_files_unique<-unique(gr_bed_files)
    
    setwd(path2)
    
    filename_20<-paste("Vockley_REF_MPRA","_GRCh37",".bed", sep='')
    
    
    export.bed(gr_bed_files_unique,con=filename_20)
    
    
  }
  
  if(length(list_BED_LogFC) >0)
  {
    
    
    BED_LogFC = unique(as.data.frame(data.table::rbindlist(list_BED_LogFC, fill = T)))
    
    
    cat("BED_LogFC_0\n")
    cat(str(BED_LogFC))
    cat("\n")
    
    
    BED_LogFC$chr<-factor(BED_LogFC$chr,
                                levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                                         "chr17","chr18","chr19","chr20","chr21","chr22","chr23","chrX","chrY"),
                                ordered=T)
    
    BED_LogFC<-BED_LogFC[order(BED_LogFC$chr,BED_LogFC$start, decreasing = F),]
    
    cat("BED_LogFC_1\n")
    cat(str(BED_LogFC))
    cat("\n")
    
    
    ### export bed of tiles GRCh37----
    
  
    
    gr_bed_files <- GRanges(
      seqnames = BED_LogFC$chr,
      ranges=IRanges(
        start=BED_LogFC$start,
        end=BED_LogFC$end,
        name=BED_LogFC$name,
        score=BED_LogFC$LogFC
      ))
    
    gr_bed_files_unique<-unique(gr_bed_files)
    
    setwd(path2)
    
    filename_20<-paste("LogFC_MPRA","_GRCh37",".bed", sep='')
    
    
    export.bed(gr_bed_files_unique,con=filename_20)
    
    
  }
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
