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
suppressMessages(library("reshape2", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("tzdb", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("cli", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("tidyverse", lib.loc="/nfs/team151/software/manuel_R_libs_4_1//"))
library("svglite", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("cowplot",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("digest",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("farver",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("labeling",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("ggeasy",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("R.methodsS3", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("R.oo", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("R.utils", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("ggridges", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("ggtext", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("glue", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("markdown", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")

highlight = function(x, pat, color="black", family="") {
  ifelse(grepl(pat, x), glue("<b style='font-family:{family}; color:{color}'>{x}</b>"), x)
}

opt = NULL

options(warn = 1)

file_reader = function(option_list)
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
  
 
  #### READ and transform indir ----
  
  indir = opt$indir
  
  cat("indir_\n")
  cat(sprintf(as.character(indir)))
  cat("\n")
  
  setwd(indir)
  
  file_list <- list.files(path=indir, include.dirs = FALSE)
  
  cat("file_list\n")
  cat(str(file_list))
  cat("\n")
  
  indexes_sel <- grep("_gDNA_deduplicated_sorted_unique\\.tsv\\.gz$|_cDNA_deduplicated_sorted_unique\\.tsv\\.gz$",file_list)
  
  cat("indexes_sel\n")
  cat(sprintf(as.character(indexes_sel)))
  cat("\n")
  
  file_list_sel <- as.data.frame(file_list[indexes_sel], stringsAsFactors=F)
  
  colnames(file_list_sel)<-"file"
  
  file_list_sel$file_extraction<-file_list_sel$file
  
  file_list_sel$file_extraction<-gsub("^R1[0-1]_6","K562_6",file_list_sel$file_extraction)
  file_list_sel$file_extraction<-gsub("^R1[0-1]_7","K562_7",file_list_sel$file_extraction)
  file_list_sel$file_extraction<-gsub("_gfpp","",file_list_sel$file_extraction)
  file_list_sel$file_extraction<-gsub("ALL","K562",file_list_sel$file_extraction)
  file_list_sel$file_extraction<-gsub("-GFP","R0minus",file_list_sel$file_extraction)
  file_list_sel$file_extraction<-gsub("pGFP","R0plus",file_list_sel$file_extraction)
  file_list_sel$file_extraction<-gsub("THP-1","THP1",file_list_sel$file_extraction)
  file_list_sel$file_extraction<-gsub("^R","K562_",file_list_sel$file_extraction)
  
  
  
  
  cat("file_list_sel_1\n")
  cat(str(file_list_sel))
  cat("\n")
  
  file_list_sel$Cell_Type<-gsub("_.+$","",file_list_sel$file_extraction)
  
  
  cat("file_list_sel_2\n")
  cat(str(file_list_sel))
  cat("\n")
  
  
  file_list_sel$Replicate<-gsub("^[^_]+_","",file_list_sel$file_extraction)
  file_list_sel$Replicate<-gsub("_.+$","",file_list_sel$Replicate)
  
  #file_list_sel$Replicate[which(file_list_sel$Cell_Type == "ALL")]<-paste("ALL",file_list_sel$Replicate[which(file_list_sel$Cell_Type == "ALL")],sep='_')
  
  
  cat("file_list_sel_3\n")
  cat(str(file_list_sel))
  cat("\n")
  
  
  file_list_sel$type<-gsub("^[^_]+_[^_]+_","",file_list_sel$file_extraction)
  file_list_sel$type<-gsub("_.+$","",file_list_sel$type)
  
  cat("file_list_sel_4\n")
  cat(str(file_list_sel))
  cat("\n")
  
  #file_list_sel$Cell_Type[which(file_list_sel$Cell_Type == "ALL")]<-"K562"
  
  file_list_sel$sample<-paste(file_list_sel$Cell_Type,
                              file_list_sel$Replicate,
                              file_list_sel$type,
                              sep='_')
  
  file_list_sel<-file_list_sel[,-which(colnames(file_list_sel) == "file_extraction")]
 
  
  
  file_list_sel$type<-factor(file_list_sel$type,
                               c("gDNA","cDNA"),
                               ordered=T)
  file_list_sel$Cell_Type<-factor(file_list_sel$Cell_Type,
                             c("K562","CHRF","HL60","THP1"),
                             ordered=T)
  file_list_sel$Replicate<-factor(file_list_sel$Replicate,
                             c("6","7","14","15","16","17","18","19","R0minus","R0plus",
                               "R1","R2","R3","R4","R5","R6","R7","R8","R9","R10","R11","R12","R13","R14","R15"),
                             ordered=T)
  
  file_list_sel[order(file_list_sel$Cell_Type,file_list_sel$Replicate,file_list_sel$type),]
  
  
  
    
  cat("file_list_sel\n")
  cat(str(file_list_sel))
  cat("\n")
  cat(sprintf(as.character(names(summary(file_list_sel$Cell_Type)))))
  cat("\n")
  cat(sprintf(as.character(summary(file_list_sel$Cell_Type))))
  cat("\n")
  cat(sprintf(as.character(names(summary(file_list_sel$type)))))
  cat("\n")
  cat(sprintf(as.character(summary(file_list_sel$type))))
  cat("\n")
  cat(sprintf(as.character(names(summary(file_list_sel$Replicate)))))
  cat("\n")
  cat(sprintf(as.character(summary(file_list_sel$Replicate))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(file_list_sel$sample))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(file_list_sel$sample)))))
  cat("\n")
 
  setwd(out)
  write.table(file_list_sel,file="test_file_table.tsv",sep="\t",quote=F,row.names = F)
  
  List_RESULTS<-list()
  
  
  CT_array<-unique(file_list_sel$Cell_Type)
  
  cat("CT_array\n")
  cat(str(CT_array))
  cat("\n")
  
  list_CT<-list()
  
  
  for(iteration_Cell_Type in 1:length(CT_array))
  {
    Cell_Type_sel<-CT_array[iteration_Cell_Type]
    
    cat("--------------------------------------->\t")
    cat(sprintf(as.character(Cell_Type_sel)))
    cat("\t")
    
    file_list_CT_sel<-droplevels(file_list_sel[which(file_list_sel$Cell_Type == Cell_Type_sel),])
    
    # cat("file_list_CT_sel\n")
    # cat(str(file_list_CT_sel))
    # cat("\n")
    
    Replicates_array<-unique(file_list_CT_sel$Replicate)
    
    # cat("Replicates_array\n")
    # cat(str(Replicates_array))
    # cat("\n")
    
    list_replicate<-list()
    
    
    for(iteration_Replicate in 1:length(Replicates_array))
    {
      Replicates_sel<-Replicates_array[iteration_Replicate]
      
      cat("--------------->\t")
      cat(sprintf(as.character(Replicates_sel)))
      cat("\t")
      
      file_list_CT_sel_Replicates_sel<-droplevels(file_list_CT_sel[which(file_list_CT_sel$Replicate == Replicates_sel),])
      
      # cat("file_list_CT_sel_Replicates_sel\n")
      # cat(str(file_list_CT_sel_Replicates_sel))
      # cat("\n")
      
      
      
      
        ## gDNA counts need to be added they were sequenced with different barcodes for this sample
        
        types_array<-unique(file_list_CT_sel_Replicates_sel$type)
        
        # cat("types_array\n")
        # cat(str(types_array))
        # cat("\n")
        
        list_type<-list()
        
        for(iteration_type in 1:length(types_array))
        {
          types_sel<-types_array[iteration_type]
          
          cat("------------------------>\t")
          cat(sprintf(as.character(types_sel)))
          cat("\n")
          
          file_list_CT_Replicates_types_sel<-droplevels(file_list_CT_sel_Replicates_sel[which(file_list_CT_sel_Replicates_sel$type == types_sel),])
          
          cat("file_list_CT_Replicates_types_sel\n")
          cat(str(file_list_CT_Replicates_types_sel))
          cat("\n")
          
          sample_sel<-unique(file_list_CT_Replicates_types_sel$sample)
          
          cat("--->\t")
          cat(sprintf(as.character(sample_sel)))
          cat("\n")
          
          
          
          
          list_K562_Rep6<-list()
          
          if(Cell_Type_sel == "K562" & Replicates_sel == 6 & types_sel == "gDNA")
          {
            FILES_array<-unique(file_list_CT_Replicates_types_sel$file)
            
            # cat("FILES_array\n")
            # cat(str(FILES_array))
            # cat("\n")
            
            for(iteration_FILES_array in 1:length(FILES_array))
            {
              setwd(indir)
              
              sel_file<-file_list_CT_Replicates_types_sel$file[iteration_FILES_array]
              
              # cat("sel_file\n")
              # cat(str(sel_file))
              # cat("\n")
              
              SIZE_gate<-file.info(sel_file)$size
              
              # cat("SIZE_gate\n")
              # cat(str(SIZE_gate))
              # cat("\n")
              
              if(SIZE_gate> 0)
              {
                LINE_gate<-length(readLines(sel_file))
                
                
                # cat("LINE_gate\n")
                # cat(str(LINE_gate))
                # cat("\n")
                
                if(LINE_gate> 0)
                {
                  df<-as.data.frame(fread(file=sel_file, sep=" ", header = T, fill=TRUE), stringsAsFactors = F)
                  colnames(df)<-c("counts","seq_name_bc")
                  
                  # cat("df\n")
                  # cat(str(df))
                  # cat("\n")
                  
                  
                  matrix<-as.matrix(df[,-which(colnames(df) == "seq_name_bc")])
                  row.names(matrix)<-df$seq_name_bc
                  colnames(matrix)<-"counts"
                  # 
                  # cat("matrix\n")
                  # cat(str(matrix))
                  # cat("\n")
                  
                                   
                  list_K562_Rep6[[iteration_FILES_array]]<-matrix
                  
                  # quit(status=1)
                }#LINE_gate
                else{
                  
                  cat(sprintf("empty_file\n"))
                  cat(sprintf(as.character(sel_file)))
                  cat("\n")
                  
                  
                }
              }#SIZE_gates
            }# iteration_FILES_array
            
            
            # cat("list_K562_Rep6\n")
            # cat(str(list_K562_Rep6))
            # cat("\n")
            
            
            m1<-list_K562_Rep6[[1]]
            
            # cat("m1\n")
            # cat(str(m1))
            # cat("\n")
            # 
            # cat(str(row.names(m1)))
            # cat("\n")
            
            m2<-list_K562_Rep6[[2]]
            
            # cat("m2\n")
            # cat(str(m2))
            # cat("\n")
            
            cols <- unique(c(colnames(m1), colnames(m2)))
            
            # cat("cols\n")
            # cat(str(cols))
            # cat("\n")
            
            rows <- unique(c(rownames(m1), rownames(m2)))
            
            # cat("rows\n")
            # cat(str(rows))
            # cat("\n")
            
            dummy_m1 <- matrix(0, nrow = length(rows), ncol = length(cols), 
                               dimnames = list(rows,cols))
            
            # cat("dummy_m1\n")
            # cat(str(dummy_m1))
            # cat("\n")
            
            dummy_m2 <- dummy_m1
            
           
            
            dummy_m1[rownames(m1), colnames(m1)] <- m1
            
            # cat("dummy_m1\n")
            # cat(str(dummy_m1))
            # cat("\n")
            
            dummy_m2[rownames(m2), colnames(m2)] <- m2
            
            combined<-dummy_m1 + dummy_m2
           
           cat("combined\n")
           cat(str(combined))
           cat("\n")
           
           df<-as.data.frame(combined, stringsAsFactors = F)
           colnames(df)<-c("counts")
           df$seq_name_bc<-row.names(df)
           df$type<-types_sel
           df$Cell_Type<-Cell_Type_sel
           df$Replicate<-Replicates_sel
           df$sample<-sample_sel
           
           cat("df\n")
           cat(str(df))
           cat("\n")
          
           
           list_type[[iteration_type]]<-df
           
            
           # quit(status = 1)
            
            
            
          }
          else{
            
            list_K562_Rep7<-list()
            
            
            if(Cell_Type_sel == "K562" & Replicates_sel == 7 & types_sel == "gDNA")
            {
              FILES_array<-unique(file_list_CT_Replicates_types_sel$file)
              
              # cat("FILES_array\n")
              # cat(str(FILES_array))
              # cat("\n")
              
              for(iteration_FILES_array in 1:length(FILES_array))
              {
                setwd(indir)
                
                sel_file<-file_list_CT_Replicates_types_sel$file[iteration_FILES_array]
                
                # cat("sel_file\n")
                # cat(str(sel_file))
                # cat("\n")
                
                SIZE_gate<-file.info(sel_file)$size
                
                # cat("SIZE_gate\n")
                # cat(str(SIZE_gate))
                # cat("\n")
                
                if(SIZE_gate> 0)
                {
                  LINE_gate<-length(readLines(sel_file))
                  
                  
                  # cat("LINE_gate\n")
                  # cat(str(LINE_gate))
                  # cat("\n")
                  
                  if(LINE_gate> 0)
                  {
                    df<-as.data.frame(fread(file=sel_file, sep=" ", header = T, fill=TRUE), stringsAsFactors = F)
                    colnames(df)<-c("counts","seq_name_bc")
                    
                    # cat("df\n")
                    # cat(str(df))
                    # cat("\n")
                    
                    
                    matrix<-as.matrix(df[,-which(colnames(df) == "seq_name_bc")])
                    row.names(matrix)<-df$seq_name_bc
                    colnames(matrix)<-"counts"
                    # 
                    # cat("matrix\n")
                    # cat(str(matrix))
                    # cat("\n")
                    
                    
                    list_K562_Rep7[[iteration_FILES_array]]<-matrix
                    
                    # quit(status=1)
                  }#LINE_gate
                  else{
                    
                    cat(sprintf("empty_file\n"))
                    cat(sprintf(as.character(sel_file)))
                    cat("\n")
                    
                    
                  }
                }#SIZE_gates
              }# iteration_FILES_array
              
              
              # cat("list_K562_Rep7\n")
              # cat(str(list_K562_Rep7))
              # cat("\n")
              
              
              m1<-list_K562_Rep7[[1]]
              
              # cat("m1\n")
              # cat(str(m1))
              # cat("\n")
              # 
              # cat(str(row.names(m1)))
              # cat("\n")
              
              m2<-list_K562_Rep7[[2]]
              
              # cat("m2\n")
              # cat(str(m2))
              # cat("\n")
              
              cols <- unique(c(colnames(m1), colnames(m2)))
              
              # cat("cols\n")
              # cat(str(cols))
              # cat("\n")
              
              rows <- unique(c(rownames(m1), rownames(m2)))
              
              # cat("rows\n")
              # cat(str(rows))
              # cat("\n")
              
              dummy_m1 <- matrix(0, nrow = length(rows), ncol = length(cols), 
                                 dimnames = list(rows,cols))
              
              # cat("dummy_m1\n")
              # cat(str(dummy_m1))
              # cat("\n")
              
              dummy_m2 <- dummy_m1
              
              
              
              dummy_m1[rownames(m1), colnames(m1)] <- m1
              
              # cat("dummy_m1\n")
              # cat(str(dummy_m1))
              # cat("\n")
              
              dummy_m2[rownames(m2), colnames(m2)] <- m2
              
              combined<-dummy_m1 + dummy_m2
              
              # cat("combined\n")
              # cat(str(combined))
              # cat("\n")
              
              df<-as.data.frame(combined, stringsAsFactors = F)
              colnames(df)<-c("counts")
              df$seq_name_bc<-row.names(df)
              df$type<-types_sel
              df$Cell_Type<-Cell_Type_sel
              df$Replicate<-Replicates_sel
              df$sample<-sample_sel
              
              # cat("df\n")
              # cat(str(df))
              # cat("\n")
              
              
              list_type[[iteration_type]]<-df
              
              
              # quit(status = 1)
              
            }
            else{
              
              setwd(indir)
              
              sel_file<-file_list_CT_Replicates_types_sel$file
              
              # cat("sel_file\n")
              # cat(str(sel_file))
              # cat("\n")
              
              SIZE_gate<-file.info(sel_file)$size
              
              # cat("SIZE_gate\n")
              # cat(str(SIZE_gate))
              # cat("\n")
              
              if(SIZE_gate> 0)
              {
                LINE_gate<-length(readLines(sel_file))
                
                
                # cat("LINE_gate\n")
                # cat(str(LINE_gate))
                # cat("\n")
                
                if(LINE_gate> 0)
                {
                  df<-as.data.frame(fread(file=sel_file, sep=" ", header = T, fill=TRUE), stringsAsFactors = F)
                  colnames(df)<-c("counts","seq_name_bc")
                  
                  df$type<-types_sel
                  df$Cell_Type<-Cell_Type_sel
                  df$Replicate<-Replicates_sel
                  df$sample<-sample_sel
                  
                  # cat("df\n")
                  # cat(str(df))
                  # cat("\n")
                  
                  list_type[[iteration_type]]<-df
                  
                  
                  # quit(status = 1)
                }# LINE_gate
                else{
                  
                  cat(sprintf("empty_file\n"))
                  cat(sprintf(as.character(sel_file)))
                  cat("\n")
                  
                  
                }
              }# SIZE_gate
            }# Cell_Type_sel == "K562" & Replicates_sel == 7 & Types_sel == "gDNA"
          }# Cell_Type_sel == "K562" & Replicates_sel == 6 & Types_sel == "gDNA"
        }# iteration_Types_array
        
        if(length(list_type) >0)
        {
          TABLE_type = unique(as.data.frame(data.table::rbindlist(list_type, fill = T)))
          
          cat("TABLE_type_0\n")
          cat(str(TABLE_type))
          cat("\n")
          
          list_replicate[[iteration_Replicate]]<-TABLE_type
        }
        
    }# iteration replicate
    
    if(length(list_type) >0)
    {
      TABLE_replicate = unique(as.data.frame(data.table::rbindlist(list_replicate, fill = T)))
      
      cat("TABLE_replicate_0\n")
      cat(str(TABLE_replicate))
      cat("\n")
      
      list_CT[[iteration_Cell_Type]]<-TABLE_replicate
    }
  }# iteration_Cell_Type
  
  TABLE_MPRA_RESULTS = unique(as.data.frame(data.table::rbindlist(list_CT, fill = T)))

  cat("TABLE_MPRA_RESULTS_0\n")
  cat(str(TABLE_MPRA_RESULTS))
  cat("\n")
  cat(sprintf(as.character(names(summary(TABLE_MPRA_RESULTS$type)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(TABLE_MPRA_RESULTS$Replicate)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(TABLE_MPRA_RESULTS$Cell_Type)))))
  cat("\n")
  cat(sprintf(as.character(summary(TABLE_MPRA_RESULTS$Cell_Type))))
  cat("\n")
  # quit(status=1)
  
 
  
  TABLE_MPRA_RESULTS$type<-factor(TABLE_MPRA_RESULTS$type,
                             c("gDNA","cDNA"),
                             ordered=T)
  TABLE_MPRA_RESULTS$Cell_Type<-factor(TABLE_MPRA_RESULTS$Cell_Type,
                                  c("K562","CHRF","HL60","THP1"),
                                  ordered=T)
  TABLE_MPRA_RESULTS$Replicate<-factor(TABLE_MPRA_RESULTS$Replicate,
                                       c("6","7","14","15","16","17","18","19","R0minus","R0plus",
                                         "R1","R2","R5","R6","R7","R8","R9","R10","R11","R12","R13","R14","R15"),
                                  ordered=T)
  
  TABLE_MPRA_RESULTS[order(TABLE_MPRA_RESULTS$Cell_Type,TABLE_MPRA_RESULTS$Replicate,TABLE_MPRA_RESULTS$type),]
  
  cat("TABLE_MPRA_RESULTS_1\n")
  cat(str(TABLE_MPRA_RESULTS))
  cat("\n")
  
  
 
  TABLE_MPRA_RESULTS$batch<-as.numeric(as.factor(TABLE_MPRA_RESULTS$sample))
  
  TABLE_MPRA_RESULTS$batch<-paste('b',TABLE_MPRA_RESULTS$batch,sep='')
  
  cat("TABLE_MPRA_RESULTS_2\n")
  cat(str(TABLE_MPRA_RESULTS))
  cat("\n")
  
  TABLE_MPRA_RESULTS$ALLELE<-gsub("^[^;]+;[^;]+;[^;]+;[^;]+;","",TABLE_MPRA_RESULTS$seq_name)
  TABLE_MPRA_RESULTS$ALLELE<-gsub(";.+$","",TABLE_MPRA_RESULTS$ALLELE)
  # TABLE_MPRA_RESULTS$ALLELE<-gsub("NO_FLAG.*$","",TABLE_MPRA_RESULTS$ALLELE)
  
  
  cat("TABLE_MPRA_RESULTS_3\n")
  cat(str(TABLE_MPRA_RESULTS))
  cat("\n")
  

  
  #### SAVE FILES ----
  
  filename=paste('RAW_RESULTS','.rds',sep='')
  
  setwd(out)
  saveRDS(TABLE_MPRA_RESULTS,file=filename)
}

###################################### REACTIVATE file_reader #####################################################

regularize_matrixes = function(option_list)
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
  
  ### Read initial selection----
  
  
  Initial_Selection<-as.data.frame(readRDS(file=opt$Initial_Selection) , stringsAsFactors=F)
  
  cat("Initial_Selection\n")
  cat(str(Initial_Selection))
  cat("\n")
  
  Initial_Selection<-unique(Initial_Selection[,-which(colnames(Initial_Selection) == "maf_origin")])
  
  cat("Initial_Selection_2\n")
  cat(str(Initial_Selection))
  cat("\n")
  
  Screened_variants<-Initial_Selection[!is.na(Initial_Selection$Fig2_Annot_Category),]
  
  cat("Screened_variants\n")
  cat(str(Screened_variants))
  cat("\n")
  cat(str(unique(Screened_variants$VAR)))
  cat("\n")
  
  # data <- simulateMPRA(tr = rep(2,10), da=NULL, nbatch=2, nbc=15)
  # 
  # 
  # 
  # cat("data\n")
  # cat(str(data))
  # cat("\n")
  
  #### READ and transform LONG_MATRIX ----
  
  LONG_MATRIX = read.table(opt$LONG_MATRIX, sep="\t", stringsAsFactors = F, header = T)
  
  
  colnames(LONG_MATRIX)[which(colnames(LONG_MATRIX) == "Real_Tile")]<-"REAL_TILE"
  
  cat("----->LONG_MATRIX_0\n")
  cat(str(LONG_MATRIX))
  cat("\n")
  cat(str(unique(LONG_MATRIX$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(LONG_MATRIX$Label))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(LONG_MATRIX$Label)))))
  cat("\n")
  
 CTRL_Labels<-c("Kousik_variant","Negative_Control_Genomic_Regions","POSITIVE_CTRL","UNDETERMINED_CTRL")
 
  LONG_MATRIX_CTRL<-LONG_MATRIX[which(LONG_MATRIX$Label%in%CTRL_Labels),]
  
  cat("LONG_MATRIX_CTRL\n")
  cat(str(LONG_MATRIX_CTRL))
  cat("\n")
  cat(str(unique(LONG_MATRIX_CTRL$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(LONG_MATRIX_CTRL$Label))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(LONG_MATRIX_CTRL$Label)))))
  cat("\n")
 
  LONG_MATRIX_ASSAYED<-LONG_MATRIX[-which(LONG_MATRIX$Label%in%CTRL_Labels),]
  
  cat("LONG_MATRIX_ASSAYED\n")
  cat(str(LONG_MATRIX_ASSAYED))
  cat("\n")
  cat(str(unique(LONG_MATRIX_ASSAYED$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(LONG_MATRIX_ASSAYED$Label))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(LONG_MATRIX_ASSAYED$Label)))))
  cat("\n")
  
  LONG_MATRIX_ASSAYED_filtered<-LONG_MATRIX_ASSAYED[which(LONG_MATRIX_ASSAYED$VAR%in%Screened_variants$VAR),]
  
  cat("LONG_MATRIX_ASSAYED_filtered\n")
  cat(str(LONG_MATRIX_ASSAYED_filtered))
  cat("\n")
  cat(str(unique(LONG_MATRIX_ASSAYED_filtered$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(LONG_MATRIX_ASSAYED_filtered$Label))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(LONG_MATRIX_ASSAYED_filtered$Label)))))
  cat("\n")
  
  LONG_MATRIX<-rbind(LONG_MATRIX_CTRL,LONG_MATRIX_ASSAYED_filtered)
  
  cat("----->LONG_MATRIX_1\n")
  cat(str(LONG_MATRIX))
  cat("\n")
  cat(str(unique(LONG_MATRIX$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(LONG_MATRIX$Label))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(LONG_MATRIX$Label)))))
  cat("\n")
  
  
  # ###################################################
  # quit(status = 1)
  # 
  
  
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
 
  
  
  indx.int<-c(which(colnames(LONG_MATRIX) == "REAL_TILE_Plus_carried_variants"),
              which(colnames(LONG_MATRIX) == "VAR"),
              which(colnames(LONG_MATRIX) == "Label"),
              which(colnames(LONG_MATRIX) == "factor4"))
  
  LONG_MATRIX_subset<-unique(LONG_MATRIX[,indx.int])
  
  cat("LONG_MATRIX_subset\n")
  cat(str(LONG_MATRIX_subset))
  cat("\n")
  
  
 
  
  #### deduplicate the factor 4
  
  
  
  LONG_MATRIX_subset.dt<-data.table(LONG_MATRIX_subset, key=c("VAR","REAL_TILE_Plus_carried_variants","Label"))
  
  
  LONG_MATRIX_subset_deduplicated<-as.data.frame(LONG_MATRIX_subset.dt[,.(factor4_string=paste(factor4, collapse = "|")), 
                                                              by =key(LONG_MATRIX_subset.dt)],)

  
  
  cat("LONG_MATRIX_subset_deduplicated_1\n")
  cat(str(LONG_MATRIX_subset_deduplicated))
  cat("\n")
  
  LONG_MATRIX_subset_deduplicated$factor4<-"NA"
  
  indx<-grep("High_AT",LONG_MATRIX_subset_deduplicated$factor4_string)
  
  cat("indx_AT\n")
  cat(str(indx))
  cat("\n")
  
  
  LONG_MATRIX_subset_deduplicated$factor4[indx]<-"High_AT"
  
  indx<-grep("High_GC",LONG_MATRIX_subset_deduplicated$factor4_string)
  
  
  cat("indx_GC\n")
  cat(str(indx))
  cat("\n")
  
  LONG_MATRIX_subset_deduplicated$factor4[indx]<-"High_GC"
  
  indx<-grep("Medium",LONG_MATRIX_subset_deduplicated$factor4_string)
  
  cat("indx_Medium\n")
  cat(str(indx))
  cat("\n")
  
  
  LONG_MATRIX_subset_deduplicated$factor4[indx]<-"Medium"
  
  LONG_MATRIX_subset_deduplicated<-LONG_MATRIX_subset_deduplicated[,-which(colnames(LONG_MATRIX_subset_deduplicated) == "factor4_string")]
  
  
  LONG_MATRIX_subset_deduplicated$factor4<-factor(LONG_MATRIX_subset_deduplicated$factor4,
                                                  levels=c("Medium","High_GC","High_AT"),
                                                  ordered=T)
  
  cat("LONG_MATRIX_subset_deduplicated_Medium\n")
  cat(str(LONG_MATRIX_subset_deduplicated))
  cat("\n")
 
  
  REAL_TILE_DF<-LONG_MATRIX_subset_deduplicated
  
  
  cat("REAL_TILE_DF_0\n")
  cat(str(REAL_TILE_DF))
  cat("\n")
  cat(sprintf(as.character(names(summary(REAL_TILE_DF$factor4)))))
  cat("\n")
  cat(sprintf(as.character(summary(REAL_TILE_DF$factor4))))
  cat("\n")
 
  
  
  
 
  
  ### check compound variants
  
  
  check<-REAL_TILE_DF[grep("\\|",REAL_TILE_DF$REAL_TILE),]
  
  cat("check\n")
  cat(str(check))
  cat("\n")
  
 
  
  
  REAL_TILE_DF$chr<-gsub("_.+$","",REAL_TILE_DF$VAR)
  REAL_TILE_DF$pos<-gsub("^[^_]+_","",REAL_TILE_DF$VAR)
  REAL_TILE_DF$pos<-as.integer(gsub("_.+$","",REAL_TILE_DF$pos))
  
  
  
  cat("REAL_TILE_DF_POST_MERGE\n")
  cat(str(REAL_TILE_DF))
  cat("\n")
  
  
 
  ## REVERT NCGR that are not ctrls ----
  
  REVERT_chr<-c("chr2")
  REVERT_POS<-c(171580283,
                175191238)
  
  
  
  
  check<-REAL_TILE_DF[which(REAL_TILE_DF$pos%in%REVERT_POS & REAL_TILE_DF$chr%in%REVERT_chr),]
  
  cat("check_POST_MERGE\n")
  cat(str(check))
  cat("\n")
  
  VARS_CTRL<-unique(check$VAR)
  
  cat("VARS_CTRL_1\n")
  cat(str(VARS_CTRL))
  cat("\n")
  
  REAL_TILE_DF$Label[which(REAL_TILE_DF$VAR%in%VARS_CTRL)]<-"UNDETERMINED_CTRL"
  
  
  
  #### CTRLS check ----
  
  
  REAL_TILE_DF_NEG_CTRLS<-REAL_TILE_DF[which(REAL_TILE_DF$Label == 'Negative_Control_Genomic_Regions'),]
  
  cat("REAL_TILE_DF_NEG_CTRLS_POST_MERGE\n")
  cat(str(REAL_TILE_DF_NEG_CTRLS))
  cat("\n")
  
  VARS_CTRL<-unique(REAL_TILE_DF_NEG_CTRLS$VAR)
  
  cat("VARS_CTRL_2\n")
  cat(str(VARS_CTRL))
  cat("\n")
  
  
  
 
  #### READ RAW data ----
  
  
  setwd(out)
  
  filename=paste('RAW_RESULTS','.rds',sep='')
  TABLE_MPRA_RESULTS <-readRDS(file=filename)
  
  colnames(TABLE_MPRA_RESULTS)[which(colnames(TABLE_MPRA_RESULTS) == "ALLELE")]<-"carried_variants"
  
  cat("TABLE_MPRA_RESULTS_0\n")
  cat(str(TABLE_MPRA_RESULTS))
  cat("\n")
  
  check<-TABLE_MPRA_RESULTS[grep("\\|",TABLE_MPRA_RESULTS$carried_variants),]
  
  cat("check_0\n")
  cat(str(check))
  cat("\n")
  
  
  
  
  
  
  
  TABLE_MPRA_RESULTS$KEY<-gsub(";.+","",TABLE_MPRA_RESULTS$seq_name_bc)
  TABLE_MPRA_RESULTS$TILE<-gsub("^[^;]+;","",TABLE_MPRA_RESULTS$seq_name_bc)
  TABLE_MPRA_RESULTS$TILE<-gsub(";.+$","",TABLE_MPRA_RESULTS$TILE)
  
  # TABLE_MPRA_RESULTS$REAL_TILE<-paste(paste(TABLE_MPRA_RESULTS$KEY,TABLE_MPRA_RESULTS$TILE,sep="__"),TABLE_MPRA_RESULTS$carried_variants,sep=";")
  
  
  TABLE_MPRA_RESULTS$condition<-gsub("^[^;]+;[^;]+;","",TABLE_MPRA_RESULTS$seq_name_bc)
  TABLE_MPRA_RESULTS$condition<-gsub(";.+$","",TABLE_MPRA_RESULTS$condition)
  
  
  
  TABLE_MPRA_RESULTS$bc<-gsub("[^;]+;[^;]+;[^;]+;[^;]+;[^;]+;","",TABLE_MPRA_RESULTS$seq_name_bc)
  
  
 
  
  TABLE_MPRA_RESULTS$batch<-"NA"
  
  TABLE_MPRA_RESULTS$master_sample<-paste(TABLE_MPRA_RESULTS$Cell_Type,TABLE_MPRA_RESULTS$Replicate,sep="_")
  
  
  cat("TABLE_MPRA_RESULTS_2\n")
  cat(str(TABLE_MPRA_RESULTS))
  cat("\n")
  
  
 
  
  
  check<-TABLE_MPRA_RESULTS[grep("\\|",TABLE_MPRA_RESULTS$carried_variants),]
  
  
  
  
  cat("check_0\n")
  cat(str(check))
  cat("\n")
  
  check.carried_variants<-unique(check$carried_variants)
  
  cat("check.carried_variants\n")
  cat(str(check.carried_variants))
  cat("\n")
  
 
  
  #### ordered samples ----
  
  
  TABLE_MPRA_RESULTS<-TABLE_MPRA_RESULTS[order(TABLE_MPRA_RESULTS$Cell_Type,TABLE_MPRA_RESULTS$Replicate,TABLE_MPRA_RESULTS$type),]
  
  levels_samples<-unique(as.character(TABLE_MPRA_RESULTS$sample))
  
  TABLE_MPRA_RESULTS$sample<-factor(TABLE_MPRA_RESULTS$sample,
                                    levels=levels_samples,
                                    ordered=T)
  
  levels_master_sample<-unique(as.character(TABLE_MPRA_RESULTS$master_sample))
  
  TABLE_MPRA_RESULTS$master_sample<-factor(TABLE_MPRA_RESULTS$master_sample,
                                    levels=levels_master_sample,
                                    ordered=T)
  
  TABLE_MPRA_RESULTS$TILE<-factor(TABLE_MPRA_RESULTS$TILE,
                                           levels=c("NINE","TWO_THIRDS","HALF","ONE_THIRD","A_TENTH"),
                                           ordered=T)
  
  TABLE_MPRA_RESULTS$condition<-factor(TABLE_MPRA_RESULTS$condition,
                                           levels=c("REF","ALT"),
                                           ordered=T)
  
  TABLE_MPRA_RESULTS<-TABLE_MPRA_RESULTS[order(TABLE_MPRA_RESULTS$master_sample,TABLE_MPRA_RESULTS$type,TABLE_MPRA_RESULTS$condition,TABLE_MPRA_RESULTS$KEY,TABLE_MPRA_RESULTS$TILE),]
  
  
  
  cat("TABLE_MPRA_RESULTS_1\n")
  cat(str(TABLE_MPRA_RESULTS))
  cat("\n")
  
  #TABLE_MPRA_RESULTS$batch<-paste("b",as.numeric(TABLE_MPRA_RESULTS$sample),sep='')
  TABLE_MPRA_RESULTS$batch<-paste("b",as.numeric(TABLE_MPRA_RESULTS$master_sample),sep='')
  
  
  cat("TABLE_MPRA_RESULTS_2\n")
  cat(str(TABLE_MPRA_RESULTS))
  cat("\n")
  cat(sprintf(as.character(names(summary(TABLE_MPRA_RESULTS$condition)))))
  cat("\n")
  cat(sprintf(as.character(summary(TABLE_MPRA_RESULTS$condition))))
  cat("\n")
 
 
  
  Deconvolve_table<-unique(TABLE_MPRA_RESULTS[,c(which(colnames(TABLE_MPRA_RESULTS) == "batch"),
                                          which(colnames(TABLE_MPRA_RESULTS) == "master_sample"))])
  
  cat("Deconvolve_table_2\n")
  cat(str(Deconvolve_table))
  cat("\n")
  
  levels_batch<-unique(as.character(TABLE_MPRA_RESULTS$batch))
  
  cat("levels_batch\n")
  cat(str(levels_batch))
  cat("\n")
  
  
  
  TABLE_MPRA_RESULTS$batch<-factor(TABLE_MPRA_RESULTS$batch,
                                    levels=levels_batch,
                                    ordered=T)
  
  #### DT ----
  
  
  TABLE_MPRA_RESULTS.dt<-data.table(TABLE_MPRA_RESULTS,
                                    key=c("batch","type","KEY","TILE","condition","carried_variants"))
  
  
  AGGREGATION_MPRA_RESULTS<-as.data.frame(TABLE_MPRA_RESULTS.dt[, .(aggregated_counts=sum(counts),
                                                                    string_barcodes=paste(bc, collapse = "|")),
                                                                by=key(TABLE_MPRA_RESULTS.dt)]
                                          ,stringsAsFactors=F)
  
  
  cat("AGGREGATION_MPRA_RESULTS\n")
  cat(str(AGGREGATION_MPRA_RESULTS))
  cat("\n")
  
    
  AGGREGATION_MPRA_RESULTS.dt<-data.table(AGGREGATION_MPRA_RESULTS,
                                          key=c("batch","type","KEY","TILE","condition","carried_variants","aggregated_counts","string_barcodes"))
  
  AGGREGATION_MPRA_RESULTS_2<-as.data.frame(AGGREGATION_MPRA_RESULTS.dt[, .(n_barcodes=length(unlist(strsplit(string_barcodes, split = "\\|")))),
                                                                        by=key(AGGREGATION_MPRA_RESULTS.dt)]
                                            ,stringsAsFactors=F)
  
  
  cat("AGGREGATION_MPRA_RESULTS_2\n")
  cat(str(AGGREGATION_MPRA_RESULTS_2))
  cat("\n")
  cat(sprintf(as.character(names(summary(AGGREGATION_MPRA_RESULTS_2$n_barcodes)))))
  cat("\n")
  cat(sprintf(as.character((summary(AGGREGATION_MPRA_RESULTS_2$n_barcodes)))))
  cat("\n")
  
  
  
  #### LOOP TO CREATE MATRIXES of 15 bc ----
  
  cat("TABLE_MPRA_RESULTS_3\n")
  cat(str(TABLE_MPRA_RESULTS))
  cat("\n")
  
  check<-TABLE_MPRA_RESULTS[grep("\\|", TABLE_MPRA_RESULTS$carried_variants),]
 # check<-TABLE_MPRA_RESULTS[which(TABLE_MPRA_RESULTS$KEY == "Element_12"),]
  
  
  # setwd(out)
  # 
  # write.table(check, file="test.tsv",sep="\t",quote=F, row.names=F)
  # 
  # quit(status=1)
  
  
  KEYS_array<-unique(TABLE_MPRA_RESULTS$KEY)
  
  cat("KEYS_array\n")
  cat(str(KEYS_array))
  cat("\n")
  
  list_gDNA_DEF<-list()
  list_cDNA_DEF<-list()
  
  for(i in 1:length(KEYS_array))
  {
    KEYS_array_sel<-KEYS_array[i]
    
    cat("------------------------------->")
    cat(sprintf(as.character(KEYS_array_sel)))
    cat("\n")
    
    TABLE_MPRA_RESULTS_sel<-TABLE_MPRA_RESULTS[which(TABLE_MPRA_RESULTS$KEY == KEYS_array_sel),]
    
    # cat("TABLE_MPRA_RESULTS_sel\n")
    # cat(str(TABLE_MPRA_RESULTS_sel))
    # cat("\n")
    
    
    
    TILES_array<-unique(TABLE_MPRA_RESULTS_sel$TILE)
    
    
    # cat("TILES_array\n")
    # cat(str(TILES_array))
    # cat("\n")
    
    for(l in 1:length(TILES_array))
    {
      TILES_array_sel<-TILES_array[l]
      
      cat("------>")
      cat(sprintf(as.character(TILES_array_sel)))
      cat("\n")
      
      
      TABLE_MPRA_RESULTS_sel_TILE_sel<-TABLE_MPRA_RESULTS_sel[which(TABLE_MPRA_RESULTS_sel$TILE == TILES_array_sel),]
      
      # cat("TABLE_MPRA_RESULTS_sel_TILE_sel\n")
      # cat(str(TABLE_MPRA_RESULTS_sel_TILE_sel))
      # cat("\n")
      
      carried_variants_array<-unique(TABLE_MPRA_RESULTS_sel_TILE_sel$carried_variants)
      
      
      # cat("carried_variants_array\n")
      # cat(str(carried_variants_array))
      # cat("\n")
      
      
      Tile<-"NA"
      
      if(TILES_array_sel == "NINE")
      {
        Tile<-"TILE_1"
      }
      if(TILES_array_sel == "TWO_THIRDS")
      {
        Tile<-"TILE_2"
      }
      if(TILES_array_sel == "HALF")
      {
        Tile<-"TILE_3"
      }
      if(TILES_array_sel == "ONE_THIRD")
      {
        Tile<-"TILE_4"
      }
      if(TILES_array_sel == "A_TENTH")
      {
        Tile<-"TILE_5"
      }
      
      
    
      
      
      for(p in 1:length(carried_variants_array))
      {
        carried_variants_array_sel<-carried_variants_array[p]
        
        # cat("--->")
        # cat(sprintf(as.character(carried_variants_array_sel)))
        # cat("\n")
        
        TABLE_MPRA_RESULTS_sel_TILE_sel_carried_variants_sel<-TABLE_MPRA_RESULTS_sel_TILE_sel[which(TABLE_MPRA_RESULTS_sel_TILE_sel$carried_variants == carried_variants_array_sel),]
        
        # cat("TABLE_MPRA_RESULTS_sel_TILE_sel_carried_variants_sel\n")
        # cat(str(TABLE_MPRA_RESULTS_sel_TILE_sel_carried_variants_sel))
        # cat("\n")
        
       
        ### REAL_TILE_Plus_carried_variants<-definition
        
        REAL_TILE_Plus_carried_variants<-paste(paste(KEYS_array_sel,Tile, sep="__"),carried_variants_array_sel, sep=";")
        
        
        cat("--------------------------------------------------------------------------------------------------------->")
        cat(sprintf(as.character(REAL_TILE_Plus_carried_variants)))
        cat("\n")
        
        # cat("levels_batch\n")
        # cat(str(levels_batch))
        # cat("\n")
        
        list_gDNA_batch<-list()
        list_cDNA_batch<-list()
        
        for(n in 1: length(levels_batch))
        {
          levels_batch_sel<-levels_batch[n]
          
          
          # cat(">")
          # cat(sprintf(as.character(levels_batch_sel)))
          # cat("\n")
          
          
          TABLE_MPRA_RESULTS_sel_TILE_sel_carried_variants_sel_batch_sel<-TABLE_MPRA_RESULTS_sel_TILE_sel_carried_variants_sel[which(TABLE_MPRA_RESULTS_sel_TILE_sel_carried_variants_sel$batch == levels_batch_sel),]
          
          # cat("TABLE_MPRA_RESULTS_sel_TILE_sel_carried_variants_sel_batch_sel\n")
          # cat(str(TABLE_MPRA_RESULTS_sel_TILE_sel_carried_variants_sel_batch_sel))
          # cat("\n")
          
          
     
          
          # cat("condition_array\n")
          # cat(str(condition_array))
          # cat("\n")
          
          #### PROBLEM NOT ALL TILES HAVE A REF (compound tiles of two variants) ## se line #948
          
          
          REF_definition<-TABLE_MPRA_RESULTS_sel_TILE_sel[which(TABLE_MPRA_RESULTS_sel_TILE_sel$condition == "REF" &
                                                                  TABLE_MPRA_RESULTS_sel_TILE_sel$batch == levels_batch_sel),]
          
          # cat("REF_definition\n")
          # cat(str(REF_definition))
          # cat("\n")
          
          condition_array<-unique(TABLE_MPRA_RESULTS_sel_TILE_sel_carried_variants_sel_batch_sel$condition)
          
          
          if(length(condition_array) != 2)
          {
            cat("--------->WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! length(condition_array) != 2\t")
            cat(sprintf(as.character(levels_batch_sel)))
            cat("\n")
            
            condition_array<-c("REF","ALT")
            # quit(status = 1)
            
          }
          
          if(REAL_TILE_Plus_carried_variants == 'Element_115__TILE_1;9_5079224_C_G|9_5079248_T_C')
          {
            
            # cat("REF_definition\n")
            # cat(str(REF_definition))
            # cat("\n")
            
            # setwd(out)
            # 
            # write.table(REF_definition, file="test.tsv",sep="\t",quote=F, row.names=F)
            # 
            # quit(status = 1)
            
          }
          
          
          ## DEFINE CT, sample, master sample
          
          
          Cell_Type_DUMMY<-unique(as.character(TABLE_MPRA_RESULTS_sel_TILE_sel_carried_variants_sel_batch_sel$Cell_Type))
          master_sample_DUMMY<-unique(as.character(TABLE_MPRA_RESULTS_sel_TILE_sel_carried_variants_sel_batch_sel$master_sample))
          Replicate_DUMMY<-unique(as.character(TABLE_MPRA_RESULTS_sel_TILE_sel_carried_variants_sel_batch_sel$Replicate))
          
          
          
          list_gDNA<-list()
          list_cDNA<-list()
        
          for(h in 1:length(condition_array))
          {
            condition_array_sel<-condition_array[h]
            
            
            # cat("----->")
            # cat(sprintf(as.character(condition_array_sel)))
            # cat("\n")
            
            
            
            if(condition_array_sel == "REF")
            {
              
              condition_sel<-REF_definition
            }else{
              
              
              condition_sel<-TABLE_MPRA_RESULTS_sel_TILE_sel_carried_variants_sel_batch_sel[which(TABLE_MPRA_RESULTS_sel_TILE_sel_carried_variants_sel_batch_sel$condition == condition_array_sel),]
              
            }
            
            if(REAL_TILE_Plus_carried_variants == 'Element_115__TILE_1;9_5079224_C_G|9_5079248_T_C')
            {
            
              # cat("condition_sel\n")
              # cat(str(condition_sel))
              # cat("\n")
            }
            
            # 
            type_array<-levels(condition_sel$type)
            
            
            # cat("type_array\n")
            # cat(str(type_array))
            # cat("\n")
            # 
            
            if(dim(condition_sel)[1] == 0)
            {
              cat("NO_ALLELE_DIVERSITY_IN_gDNA_or_cDNA>")
              cat("\t")
              
              if(carried_variants_array_sel == "NCGR")
              {
                cat("BUT_NCGR")
                cat("\n")
                
                
                condition_sel<-REF_definition
                
                condition_sel$condition<-"ALT"
                
              }else{
                
                cat("NO_NCGR")
                cat("\n")
                
                setwd(out)
                
                write.table(TABLE_MPRA_RESULTS_sel_TILE_sel_carried_variants_sel_batch_sel, file="test.tsv",sep="\t",quote=F, row.names=F)
                
                quit(status = 1)
                
              }
            }
        
         
            
            
              
              list_type<-list()
              
              for(k in 1:length(type_array))
              {
                
                type_array_sel<-type_array[k]
                
                
                
                
                condition_type_sel<-condition_sel[which(condition_sel$type == type_array_sel),]
                
                if(REAL_TILE_Plus_carried_variants == 'Element_115__TILE_1;9_5079224_C_G|9_5079248_T_C')
                {
                  # cat("condition_type_sel\n")
                  # cat(str(condition_type_sel))
                  # cat("\n")
                }
                
                missing_barcodes<-NULL
                
                if(dim(condition_type_sel)[1] >0)
                {
                  
                  missing_barcodes<-15-dim(condition_type_sel)[1]
                  
                  # cat("missing_barcodes\n")
                  # cat(str(missing_barcodes))
                  # cat("\n")
                }
                else{
                  
                  cat("NO_ALLELE_DIVERSITY_IN:\t")
                  cat(sprintf(as.character(type_array_sel)))
                  cat("\n")
                  
                  # setwd(out)
                  # 
                  # write.table(condition_sel, file="test.tsv",sep="\t",quote=F, row.names=F)
                  # 
                  # quit(status = 1)
                  
                  missing_barcodes<-15
                  
                  cat("missing_barcodes\n")
                  cat(str(missing_barcodes))
                  cat("\n")
                }
                
                if(missing_barcodes > 0)
                {
                  missing_barcodes_array<-seq(1,missing_barcodes,by=1)
                  
                  # cat("missing_barcodes_array\n")
                  # cat(str(missing_barcodes_array))
                  # cat("\n")
                  
                  temp<- data.frame(matrix(vector(), 0, 
                                           dim(condition_type_sel)[2],
                                           dimnames=list(c(),
                                                         colnames(condition_type_sel)
                                           )),
                                    stringsAsFactors=F)
                  
                  # colnames(temp)<-colnames(DEF_TABLE_DESIGN_DROPOUTS)
                  
                  # cat("temp_PRE\n")
                  # cat(str(temp))
                  # cat("\n")
                  
                  
                  
                  for(missing_barcodes_array_iteration in 1:length(missing_barcodes_array))
                  {
                    iteration_sel<-missing_barcodes_array[missing_barcodes_array_iteration]
                    
                    # cat("iteration_sel\n")
                    # cat(str(iteration_sel))
                    # cat("\n")
                    
                    DUMMY_token<-paste("DUMMY",iteration_sel,sep="_")
                    sample_DUMMY<-paste(master_sample_DUMMY,as.character(type_array_sel),sep="_")
                    seq_name_bc_DUMMY<-paste(KEYS_array_sel,TILES_array_sel,as.character(condition_array_sel),"NA",carried_variants_array_sel,DUMMY_token,sep=";")
                    
                    
                    
                    temp_df<-as.data.frame(cbind(0,seq_name_bc_DUMMY,sample_DUMMY,as.character(type_array_sel),Cell_Type_DUMMY,Replicate_DUMMY,levels_batch_sel,carried_variants_array_sel,
                                                 KEYS_array_sel,TILES_array_sel,as.character(condition_array_sel),DUMMY_token,master_sample_DUMMY), stringsAsFactor=F)
                    
                    
                    
                    colnames(temp_df)<-colnames(temp)
                    
                    
                    # cat("temp_df\n")
                    # cat(str(temp_df))
                    # cat("\n")
                    
                    temp<-rbind(temp,temp_df)
                    
                    #condition_type_sel[position_to_fill,]<-c()
                    
                    
                    
                  }# missing_barcodes_array_iteration
                  
                  
                  # cat("temp_POST\n")
                  # cat(str(temp))
                  # cat("\n")
                  
                  if(dim(condition_type_sel)[1] >0)
                  {
                    condition_type_sel<-rbind(condition_type_sel,temp)
                    
                    # cat("condition_type_sel_POST\n")
                    # cat(str(condition_type_sel))
                    # cat("\n")
                  }else{
                    
                    condition_type_sel<-temp
                    
                    # cat("condition_type_sel_POST\n")
                    # cat(str(condition_type_sel))
                    # cat("\n")
                    
                  }
                  
                }else{
                  
                  condition_type_sel<-condition_type_sel
                }# missing_barcodes > 0
                  
                if(REAL_TILE_Plus_carried_variants == 'Element_115__TILE_1;9_5079224_C_G|9_5079248_T_C')
                {
                  # cat("condition_type_sel_POST\n")
                  # cat(str(condition_type_sel))
                  # cat("\n")
                }
                  
                  
                  vector_bc_numbers<-0
                  
                  if(condition_array_sel == "REF")
                  {
                    
                    vector_bc_numbers<-seq(1,15,by=1)
                  }else{
                    
                    vector_bc_numbers<-seq(16,30,by=1)
                    
                  }
                  
                  condition_type_sel$bc_name<-paste("bc",vector_bc_numbers,sep='')
                  
                  condition_type_sel$REAL_TILE_Plus_carried_variants<-REAL_TILE_Plus_carried_variants
                  
                  condition_type_sel$counts<-as.numeric(condition_type_sel$counts)
                  
                  
                  condition_type_sel_subset<-unique(condition_type_sel[,c(which(colnames(condition_type_sel) == "REAL_TILE_Plus_carried_variants"),which(colnames(condition_type_sel) == "batch"),
                                                                   which(colnames(condition_type_sel) == "bc_name"),which(colnames(condition_type_sel) == "counts"))])
                  
                  condition_type_sel_subset$batch_bc_name<-paste(condition_type_sel_subset$batch,condition_type_sel_subset$bc_name,sep="_")
                  
                  # cat("condition_type_sel_subset_0\n")
                  # cat(str(condition_type_sel_subset))
                  # cat("\n")
                  
                  condition_type_sel_subset<-condition_type_sel_subset[,-c(which(colnames(condition_type_sel_subset) == "batch"),
                                                                           which(colnames(condition_type_sel_subset) == "bc_name"))]
                  
                  # cat("condition_type_sel_subset_1\n")
                  # cat(str(condition_type_sel_subset))
                  # cat("\n")
                  
                  condition_type_sel_subset_wide<-as.data.frame(pivot_wider(condition_type_sel_subset,
                                                             id_cols=REAL_TILE_Plus_carried_variants,
                                                             names_from=batch_bc_name,
                                                             values_from=counts), stringsAsFactors=F)
                  
                  if(REAL_TILE_Plus_carried_variants == 'Element_115__TILE_1;9_5079224_C_G|9_5079248_T_C')
                  {
                  
                    # cat("condition_type_sel_subset_wide\n")
                    # str(condition_type_sel_subset_wide)
                    # cat("\n")
                    # 
                    # cat("type_array_sel\n")
                    # str(as.character(type_array_sel))
                    # cat("\n")
                  }
                  
                  
                  list_type[[as.character(type_array_sel)]]<-condition_type_sel_subset_wide
              }# k type_array
              
              names_RESULTS_1_array<-names(list_type)
              
              if(REAL_TILE_Plus_carried_variants == 'Element_115__TILE_1;9_5079224_C_G|9_5079248_T_C')
              {
                
                # cat("names_RESULTS_1_array\n")
                # cat(str(names_RESULTS_1_array))
                # cat("\n")
              }
              
              for(iteration_names_RESULTS_1_array in 1:length(names_RESULTS_1_array))
              {
                names_RESULTS_1_array_sel<-names_RESULTS_1_array[iteration_names_RESULTS_1_array]
                
                # cat(">")
                # cat(sprintf(as.character(names_RESULTS_1_array_sel)))
                # cat("\n")
                
                RESULTS_1 = list_type[[names_RESULTS_1_array_sel]]
                  # unique(as.data.frame(data.table::rbindlist(, fill = T)))
                
                if(REAL_TILE_Plus_carried_variants == 'Element_115__TILE_1;9_5079224_C_G|9_5079248_T_C')
                {
              
                  # cat("RESULTS_1\n")
                  # cat(str(RESULTS_1))
                  # cat("\n")
                }
                
               
                
                if(names_RESULTS_1_array_sel == "gDNA")
                {
                  list_gDNA[[as.character(condition_array_sel)]]<-RESULTS_1
                  
                  
                }else{
                  
                  if(names_RESULTS_1_array_sel == "cDNA")
                  {
                    list_cDNA[[as.character(condition_array_sel)]]<-RESULTS_1
                    
                  }
                }
              }# names_RESULTS_1_array
          }# h condition array
          
          ## gDNA
          
          gDNA_REF<-list_gDNA[["REF"]]
          
          if(REAL_TILE_Plus_carried_variants == 'Element_115__TILE_1;9_5079224_C_G|9_5079248_T_C')
          {
            
            # cat("gDNA_REF\n")
            # cat(str(gDNA_REF))
            # cat("\n")
          }
          
          gDNA_ALT<-list_gDNA[["ALT"]]
          
          if(REAL_TILE_Plus_carried_variants == 'Element_115__TILE_1;9_5079224_C_G|9_5079248_T_C')
          {
            
            # cat("gDNA_ALT\n")
            # cat(str(gDNA_ALT))
            # cat("\n")
          }
          
         
          
          gDNA_DEF<-merge(gDNA_REF,
                          gDNA_ALT,
                          by="REAL_TILE_Plus_carried_variants",
                          all=T)
          
          if(REAL_TILE_Plus_carried_variants == 'Element_115__TILE_1;9_5079224_C_G|9_5079248_T_C')
          {
          
            # cat("gDNA_DEF\n")
            # cat(str(gDNA_DEF))
            # cat("\n")
          }
          
          if(dim(gDNA_DEF)[2] != 31)
          {
            
            cat("WARNING IN BATCH not 31 columns in regularized matrix gDNA>")
            cat("\n")
            setwd(out)
            
            write.table(gDNA_DEF, file="test.tsv",sep="\t",quote=F, row.names=F)
            
            quit(status = 1)
            
          }
          
          list_gDNA_batch[[as.character(levels_batch_sel)]]<-gDNA_DEF
          
          ## cDNA
          
          cDNA_REF<-list_cDNA[["REF"]]
          
          if(REAL_TILE_Plus_carried_variants == 'Element_115__TILE_1;9_5079224_C_G|9_5079248_T_C')
          {
            # cat("cDNA_REF\n")
            # cat(str(cDNA_REF))
            # cat("\n")
          }
          
          cDNA_ALT<-list_cDNA[["ALT"]]
          
          if(REAL_TILE_Plus_carried_variants == 'Element_115__TILE_1;9_5079224_C_G|9_5079248_T_C')
          {
            # cat("cDNA_ALT\n")
            # cat(str(cDNA_ALT))
            # cat("\n")
          }
          
          cDNA_DEF<-merge(cDNA_REF,
                          cDNA_ALT,
                          by="REAL_TILE_Plus_carried_variants",
                          all=T)
          
          if(REAL_TILE_Plus_carried_variants == 'Element_115__TILE_1;9_5079224_C_G|9_5079248_T_C')
          {
            # cat("cDNA_DEF\n")
            # cat(str(cDNA_DEF))
            # cat("\n")
          }
          
          if(dim(cDNA_DEF)[2] != 31)
          {
            
            cat("WARNING IN BATCH not 31 columns in regularized matrix cDNA>")
            cat("\n")
            setwd(out)
            
            write.table(cDNA_DEF, file="test.tsv",sep="\t",quote=F, row.names=F)
            
            quit(status = 1)
            
          }
          
          list_cDNA_batch[[as.character(levels_batch_sel)]]<-cDNA_DEF
          
        }#n levels_batches
        
        ### gDNA
        
        names_gDNA_batch_array<-names(list_gDNA_batch)
        
        # cat("names_gDNA_batch_array\n")
        # cat(str(names_gDNA_batch_array))
        # cat("\n")
        
        
        attractor_gDNA<-NULL
        
        for(iteration_names_gDNA_batch_array in 1:length(names_gDNA_batch_array))
        {
          
          names_gDNA_batch_array_sel<-names_gDNA_batch_array[iteration_names_gDNA_batch_array]
          
          # cat(sprintf(as.character(names_gDNA_batch_array_sel)))
          # cat("\n")
          
          
          gDNA_object_after_batch = list_gDNA_batch[[names_gDNA_batch_array_sel]]
          
          # cat("gDNA_object_after_batch\n")
          # cat(str(gDNA_object_after_batch))
          # cat("\n")
          
          
          
          
          
          if(iteration_names_gDNA_batch_array == 1)
          {
            
            attractor_gDNA<-gDNA_object_after_batch
            
          }else{
            
            attractor_gDNA<-merge(attractor_gDNA,
                                  gDNA_object_after_batch,
                                  by="REAL_TILE_Plus_carried_variants",
                                  all=T)
            
          }
          
          
        }
        
        attractor_gDNA.matrix<-as.matrix(attractor_gDNA[,-1])
        row.names(attractor_gDNA.matrix)<-attractor_gDNA$REAL_TILE_Plus_carried_variants
        
        # cat("attractor_gDNA.matrix\n")
        # cat(str(attractor_gDNA.matrix))
        # cat("\n")
        
        
        if(dim(attractor_gDNA.matrix)[2] != 1050)
        {
          
          cat("WARNING IN  REAL_TILE_Plus_carried_variants not 1050 columns in regularized matrix gDNA>")
          cat("\n")
          cat(sprintf(as.character(dim(attractor_gDNA.matrix)[2])))
          cat("\n")
          
          setwd(out)
          
          write.table(attractor_gDNA, file="test.tsv",sep="\t",quote=F, row.names=F)
          
          quit(status = 1)
          
        }
        
        
        list_gDNA_DEF[[REAL_TILE_Plus_carried_variants]]<-attractor_gDNA.matrix
        
        ### cDNA
        
        names_cDNA_batch_array<-names(list_cDNA_batch)
        
        # cat("names_cDNA_batch_array\n")
        # cat(str(names_cDNA_batch_array))
        # cat("\n")
        
        
        attractor_cDNA<-NULL
        
        for(iteration_names_cDNA_batch_array in 1:length(names_cDNA_batch_array))
        {
          
          names_cDNA_batch_array_sel<-names_cDNA_batch_array[iteration_names_cDNA_batch_array]
          
          # cat(sprintf(as.character(names_cDNA_batch_array_sel)))
          # cat("\n")
          
          
          cDNA_object_after_batch = list_cDNA_batch[[names_cDNA_batch_array_sel]]
          
          # cat("cDNA_object_after_batch\n")
          # cat(str(cDNA_object_after_batch))
          # cat("\n")
          
          
          
          
          
          if(iteration_names_cDNA_batch_array == 1)
          {
            
            attractor_cDNA<-cDNA_object_after_batch
            
          }else{
            
            attractor_cDNA<-merge(attractor_cDNA,
                                  cDNA_object_after_batch,
                                  by="REAL_TILE_Plus_carried_variants",
                                  all=T)
            
          }
          
          
        }
        
        attractor_cDNA.matrix<-as.matrix(attractor_cDNA[,-1])
        row.names(attractor_cDNA.matrix)<-attractor_cDNA$REAL_TILE_Plus_carried_variants
        
        # cat("attractor_cDNA.matrix\n")
        # cat(str(attractor_cDNA.matrix))
        # cat("\n")
        
        
        if(dim(attractor_cDNA.matrix)[2] != 1050)
        {
          
          cat("WARNING IN  REAL_TILE_Plus_carried_variants not 1050 columns in regularized matrix cDNA>")
          cat("\n")
          setwd(out)
          
          write.table(attractor_cDNA, file="test.tsv",sep="\t",quote=F, row.names=F)
          
          quit(status = 1)
          
        }
        
        
        list_cDNA_DEF[[REAL_TILE_Plus_carried_variants]]<-attractor_cDNA.matrix
      }#p carried variants
    }#l TILE
  }#i KEY
  
  
  ### gDNA
  
  names_gDNA_DEF_array<-names(list_gDNA_DEF)
  
  cat("names_gDNA_DEF_array\n")
  cat(str(names_gDNA_DEF_array))
  cat("\n")
  
  
  appender_gDNA<-NULL
  
  for(iteration_names_gDNA_DEF_array in 1:length(names_gDNA_DEF_array))
  {
    
    names_gDNA_DEF_array_sel<-names_gDNA_DEF_array[iteration_names_gDNA_DEF_array]
    
    # cat(sprintf(as.character(names_gDNA_DEF_array_sel)))
    # cat("\n")
    
    
    gDNA_object_after_DEF = list_gDNA_DEF[[names_gDNA_DEF_array_sel]]
    
    # cat("gDNA_object_after_DEF\n")
    # cat(str(gDNA_object_after_DEF))
    # cat("\n")
    
    
    
    
    
    if(iteration_names_gDNA_DEF_array == 1)
    {
      
      appender_gDNA<-gDNA_object_after_DEF
      
    }else{
      
      appender_gDNA<-rbind(appender_gDNA,
                            gDNA_object_after_DEF)
      
    }
    
    
  }
  
  cat("appender_gDNA\n")
  cat(str(appender_gDNA))
  cat("\n") 
  
  
  
  gDNA.df <- data.frame(REAL_TILE_Plus_carried_variants = row.names(appender_gDNA), appender_gDNA)
  
  cat("gDNA.df\n")
  cat(str(gDNA.df))
  cat("\n")
  
  gDNA.df.m<-melt(gDNA.df, id.vars = "REAL_TILE_Plus_carried_variants")
  
  cat("gDNA.df.m_0\n")
  cat(str(gDNA.df.m))
  cat("\n")
  
  gDNA.df.m$type<-"gDNA"
  gDNA.df.m$batch<-gsub("_.+$","",gDNA.df.m$variable)
  gDNA.df.m$bc<-gsub("^[^_]+_","",gDNA.df.m$variable)
  
  REF_bc<-paste("bc",seq(1,15,by=1),sep='')
  ALT_bc<-paste("bc",seq(16,30,by=1),sep='')
  
  
  gDNA.df.m$condition<-"NA"
  
  gDNA.df.m$condition[gDNA.df.m$bc%in%REF_bc]<-"REF"
  gDNA.df.m$condition[gDNA.df.m$bc%in%ALT_bc]<-"ALT"
  
  gDNA.df.m<-gDNA.df.m[,-which(colnames(gDNA.df.m) == "variable")]
  
  colnames(gDNA.df.m)[which(colnames(gDNA.df.m) == "value")]<-"counts"
  
  cat("gDNA.df.m_1\n")
  cat(str(gDNA.df.m))
  cat("\n")
  
  
  ### cDNA
  
  names_cDNA_DEF_array<-names(list_cDNA_DEF)
  
  cat("names_cDNA_DEF_array\n")
  cat(str(names_cDNA_DEF_array))
  cat("\n")
  
  
  appender_cDNA<-NULL
  
  for(iteration_names_cDNA_DEF_array in 1:length(names_cDNA_DEF_array))
  {
    
    names_cDNA_DEF_array_sel<-names_cDNA_DEF_array[iteration_names_cDNA_DEF_array]
    
    # cat(sprintf(as.character(names_cDNA_DEF_array_sel)))
    # cat("\n")
    
    
    cDNA_object_after_DEF = list_cDNA_DEF[[names_cDNA_DEF_array_sel]]
    
    # cat("cDNA_object_after_DEF\n")
    # cat(str(cDNA_object_after_DEF))
    # cat("\n")
    
    
    
    
    
    if(iteration_names_cDNA_DEF_array == 1)
    {
      
      appender_cDNA<-cDNA_object_after_DEF
      
    }else{
      
      appender_cDNA<-rbind(appender_cDNA,
                           cDNA_object_after_DEF)
      
    }
    
    
  }
  
  cat("appender_cDNA\n")
  cat(str(appender_cDNA))
  cat("\n") 
  
  
  
  cDNA.df <- data.frame(REAL_TILE_Plus_carried_variants = row.names(appender_cDNA), appender_cDNA)
  
  cat("cDNA.df\n")
  cat(str(cDNA.df))
  cat("\n")
  
  cDNA.df.m<-melt(cDNA.df, id.vars = "REAL_TILE_Plus_carried_variants")
  
  cat("cDNA.df.m_0\n")
  cat(str(cDNA.df.m))
  cat("\n")
  
  cDNA.df.m$type<-"cDNA"
  cDNA.df.m$batch<-gsub("_.+$","",cDNA.df.m$variable)
  cDNA.df.m$bc<-gsub("^[^_]+_","",cDNA.df.m$variable)
  
  REF_bc<-paste("bc",seq(1,15,by=1),sep='')
  ALT_bc<-paste("bc",seq(16,30,by=1),sep='')
  
  
  cDNA.df.m$condition<-"NA"
  
  cDNA.df.m$condition[cDNA.df.m$bc%in%REF_bc]<-"REF"
  cDNA.df.m$condition[cDNA.df.m$bc%in%ALT_bc]<-"ALT"
  
  cDNA.df.m<-cDNA.df.m[,-which(colnames(cDNA.df.m) == "variable")]
  
  colnames(cDNA.df.m)[which(colnames(cDNA.df.m) == "value")]<-"counts"
  
  cat("cDNA.df.m_1\n")
  cat(str(cDNA.df.m))
  cat("\n")
  
  
  #### extended Rosetta ----
  
  
  Rosetta_extended<-rbind(gDNA.df.m,
                          cDNA.df.m)
  
  
  cat("Rosetta_extended_1\n")
  cat(str(Rosetta_extended))
  cat("\n")
  
  Rosetta_extended<-merge(Rosetta_extended,
                          Deconvolve_table,
                          by="batch",
                          all=T)
  
  cat("Rosetta_extended_2\n")
  cat(str(Rosetta_extended))
  cat("\n")
  
  Rosetta_extended$Cell_Type<-gsub("_.+$","",Rosetta_extended$master_sample)
  Rosetta_extended$Replicate<-gsub("^[^_]+_","",Rosetta_extended$master_sample)
  
  cat("Rosetta_extended_3\n")
  cat(str(Rosetta_extended))
  cat("\n")
  
  
  Rosetta_extended$KEY<-gsub("__.+$","",Rosetta_extended$REAL_TILE_Plus_carried_variants)
  Rosetta_extended$TILE<-gsub("^[^_]+_","",Rosetta_extended$REAL_TILE_Plus_carried_variants)
  Rosetta_extended$TILE<-gsub("^[^_]+__","",Rosetta_extended$TILE)
  Rosetta_extended$TILE<-gsub(";.+$","",Rosetta_extended$TILE)
  
  Rosetta_extended$carried_variants<-gsub("^[^_]+_","",Rosetta_extended$REAL_TILE_Plus_carried_variants)
  Rosetta_extended$carried_variants<-gsub("^[^_]+__","",Rosetta_extended$carried_variants)
  Rosetta_extended$carried_variants<-gsub("[^;]+;","",Rosetta_extended$carried_variants)
  
  
  
  
  
  cat("Rosetta_extended_4\n")
  cat(str(Rosetta_extended))
  cat("\n")
  
 
  
  
  Rosetta_extended<-merge(Rosetta_extended,
                          REAL_TILE_DF,
                          by="REAL_TILE_Plus_carried_variants",
                          all=T)
  
  
  cat("Rosetta_extended_5\n")
  cat(str(Rosetta_extended))
  cat("\n")
  
  
  Rosetta_extended$type<-factor(Rosetta_extended$type,
                                levels=c("gDNA","cDNA"),
                                ordered=T)
  Rosetta_extended$Cell_Type<-factor(Rosetta_extended$Cell_Type,
                                levels=c("K562","CHRF","HL60","THP1"),
                                ordered=T)
  levels_batch<-unique(as.character(TABLE_MPRA_RESULTS$batch))
  
  cat("levels_batch\n")
  cat(str(levels_batch))
  cat("\n")
  
  
  Rosetta_extended$batch<-factor(Rosetta_extended$batch,
                                levels=levels_batch,
                                ordered=T)
  
  Rosetta_extended$condition<-factor(Rosetta_extended$condition,
                                levels=c("REF","ALT"),
                                ordered=T)
  
  Rosetta_extended$TILE<-factor(Rosetta_extended$TILE,
                                  levels=c("TILE_1","TILE_2","TILE_3","TILE_4","TILE_5"),
                                  ordered=T)
  
  levels_master_sample<-unique(as.character(TABLE_MPRA_RESULTS$master_sample))
  
  Rosetta_extended$master_sample<-factor(Rosetta_extended$master_sample,
                                           levels=levels_master_sample,
                                           ordered=T)
  Rosetta_extended$bc<-factor(Rosetta_extended$bc,
                                     levels=c(REF_bc,ALT_bc),
                                     ordered=T)
  
  
  cat("Rosetta_extended_6_ORDERED_FACTORS\n")
  cat(str(Rosetta_extended))
  cat("\n")
  
  Rosetta_extended<-Rosetta_extended[order(Rosetta_extended$master_sample,Rosetta_extended$type,Rosetta_extended$condition,Rosetta_extended$KEY,Rosetta_extended$TILE),]
  
  CTRL_Labels<-c("Kousik_variant","Negative_Control_Genomic_Regions","POSITIVE_CTRL","UNDETERMINED_CTRL")
  
  Rosetta_extended_CTRL<-Rosetta_extended[which(Rosetta_extended$Label%in%CTRL_Labels),]
  
  cat("Rosetta_extended_CTRL\n")
  cat(str(Rosetta_extended_CTRL))
  cat("\n")
  cat(str(unique(Rosetta_extended_CTRL$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Rosetta_extended_CTRL$Label))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Rosetta_extended_CTRL$Label)))))
  cat("\n")
  
  Rosetta_extended_ASSAYED<-Rosetta_extended[-which(Rosetta_extended$Label%in%CTRL_Labels),]
  
  cat("Rosetta_extended_ASSAYED\n")
  cat(str(Rosetta_extended_ASSAYED))
  cat("\n")
  cat(str(unique(Rosetta_extended_ASSAYED$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Rosetta_extended_ASSAYED$Label))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Rosetta_extended_ASSAYED$Label)))))
  cat("\n")
  
  Rosetta_extended_ASSAYED_filtered<-Rosetta_extended_ASSAYED[which(Rosetta_extended_ASSAYED$VAR%in%Screened_variants$VAR),]
  
  cat("Rosetta_extended_ASSAYED_filtered\n")
  cat(str(Rosetta_extended_ASSAYED_filtered))
  cat("\n")
  cat(str(unique(Rosetta_extended_ASSAYED_filtered$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Rosetta_extended_ASSAYED_filtered$Label))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Rosetta_extended_ASSAYED_filtered$Label)))))
  cat("\n")
  
  Rosetta_extended<-rbind(Rosetta_extended_CTRL,Rosetta_extended_ASSAYED_filtered)
  
  cat("----->Rosetta_extended_1\n")
  cat(str(Rosetta_extended))
  cat("\n")
  cat(str(unique(Rosetta_extended$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Rosetta_extended$Label))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Rosetta_extended$Label)))))
  cat("\n")
  
  # ########################################################################################################################################################################################################################################################
  # quit(status = 1)
  
  
  indx.int<-c(which(colnames(Rosetta_extended) == "REAL_TILE_Plus_carried_variants"),
              which(colnames(Rosetta_extended) == "KEY"),
              which(colnames(Rosetta_extended) == "TILE"),
              which(colnames(Rosetta_extended) == "carried_variants"),
              which(colnames(Rosetta_extended) == "VAR"),
              which(colnames(Rosetta_extended) == "Label"),
              which(colnames(Rosetta_extended) == "factor4"))
              
              
  # which(colnames(Rosetta_extended) == "Cell_Type"),
  # which(colnames(Rosetta_extended) == "Replicate"),
  # which(colnames(Rosetta_extended) == "master_sample"),
              
  Rosetta_subset<-unique(Rosetta_extended[,indx.int])
  
  cat("Rosetta_subset_0\n")
  cat(str(Rosetta_subset))
  cat("\n")
  
  
  
  
  
              
  
  
  
 
 
  
 
  #### SAVE ----
  
  setwd(out)

  
  filename_1<-paste("matrix_gDNA_PRE",".rds", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  appender_gDNA<-as.data.frame(appender_gDNA)
  
  appender_gDNA$REAL_TILE_Plus_carried_variants<-row.names(appender_gDNA)
  
  cat("appender_gDNA_PRE\n")
  cat(str(appender_gDNA))
  cat("\n")
  cat(str(unique(appender_gDNA$REAL_TILE_Plus_carried_variants)))
  cat("\n")
  
  appender_gDNA_filtered<-appender_gDNA[which(appender_gDNA$REAL_TILE_Plus_carried_variants%in%Rosetta_subset$REAL_TILE_Plus_carried_variants),]
  
  cat("appender_gDNA_filtered_PRE\n")
  cat(str(appender_gDNA_filtered))
  cat("\n")
  cat(str(unique(appender_gDNA_filtered$REAL_TILE_Plus_carried_variants)))
  cat("\n")
  
  
  saveRDS(appender_gDNA_filtered,file=filename_1)#,sep="\t",row.names = F,quote=F)
  
  filename_1<-paste("matrix_cDNA_PRE",".rds", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  appender_cDNA<-as.data.frame(appender_cDNA)
  
  appender_cDNA$REAL_TILE_Plus_carried_variants<-row.names(appender_cDNA)
  
  cat("appender_cDNA_PRE\n")
  cat(str(appender_cDNA))
  cat("\n")
  cat(str(unique(appender_cDNA$REAL_TILE_Plus_carried_variants)))
  cat("\n")
  
  appender_cDNA_filtered<-appender_cDNA[which(appender_cDNA$REAL_TILE_Plus_carried_variants%in%Rosetta_subset$REAL_TILE_Plus_carried_variants),]
  
  cat("appender_cDNA_filtered_PRE\n")
  cat(str(appender_cDNA_filtered))
  cat("\n")
  cat(str(unique(appender_cDNA_filtered$REAL_TILE_Plus_carried_variants)))
  cat("\n")
  
  
  saveRDS(appender_cDNA_filtered,file=filename_1)#,sep="\t",row.names = F,quote=F)
  
  
  filename_1<-paste("Rosetta_extended",".rds", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  saveRDS(Rosetta_extended, file=filename_1)
  
 
  
  
  filename_1<-paste("Rosetta_df",".rds", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  saveRDS(Rosetta_subset, file=filename_1)
  
  
  
  
  
  
}

NOR_and_MPRAnalyze_files = function(option_list)
{
  
  library("MPRAnalyze", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
  library("gtools", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/")
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
  library("matrixStats", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
  library("MatrixGenerics", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
  library("SummarizedExperiment", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
  
  
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
  
  #### READ RAW data ----
  
  
  setwd(out)
  
  filename=paste('Rosetta_extended','.rds',sep='')
  Rosetta_extended <-readRDS(file=filename)
  
  
  cat("Rosetta_extended_0\n")
  cat(str(Rosetta_extended))
  cat("\n")
  cat(sprintf(as.character(names(summary(Rosetta_extended$batch)))))
  cat("\n")
  cat(sprintf(as.character(summary(Rosetta_extended$batch))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Rosetta_extended$type)))))
  cat("\n")
  cat(sprintf(as.character(summary(Rosetta_extended$type))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Rosetta_extended$condition)))))
  cat("\n")
  cat(sprintf(as.character(summary(Rosetta_extended$condition))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Rosetta_extended$master_sample)))))
  cat("\n")
  cat(sprintf(as.character(summary(Rosetta_extended$master_sample))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Rosetta_extended$bc)))))
  cat("\n")
  cat(sprintf(as.character(summary(Rosetta_extended$bc))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Rosetta_extended$Cell_Type)))))
  cat("\n")
  cat(sprintf(as.character(summary(Rosetta_extended$Cell_Type))))
  cat("\n")
  
  Deconvolve_table<-unique(Rosetta_extended[,c(which(colnames(Rosetta_extended) == "master_sample"),
                                               which(colnames(Rosetta_extended) == "batch"))])
  
  cat("Deconvolve_table_0\n")
  cat(str(Deconvolve_table))
  cat("\n")
  rm(Rosetta_extended)
  
  filename=paste('Rosetta_df','.rds',sep='')
  Rosetta <-readRDS(file=filename)
  
  
  cat("Rosetta_0\n")
  cat(str(Rosetta))
  cat("\n")
  
  filename=paste('matrix_gDNA_PRE','.rds',sep='')
  gDNA_PRE <-readRDS(file=filename)
  
  # cat("gDNA_PRE\n")
  # cat(str(gDNA_PRE))
  # cat("\n")
  
  matrix_gDNA_PRE<-as.matrix(gDNA_PRE[,-which(colnames(gDNA_PRE) == "REAL_TILE_Plus_carried_variants")])
  
  row.names(matrix_gDNA_PRE)<-gDNA_PRE$REAL_TILE_Plus_carried_variants
  
  cat("matrix_gDNA_PRE\n")
  cat(str(matrix_gDNA_PRE))
  cat("\n")
  
  rm(gDNA_PRE)
  
  filename=paste('matrix_cDNA_PRE','.rds',sep='')
  #cDNA_PRE <-as.data.frame(fread(file=filename, sep="\t",header=T), stringsAsFactors = F)
  cDNA_PRE <-readRDS(file=filename)
  

  # cat("cDNA_PRE\n")
  # cat(str(cDNA_PRE))
  # cat("\n")

  matrix_cDNA_PRE<-as.matrix(cDNA_PRE[,-which(colnames(cDNA_PRE) == "REAL_TILE_Plus_carried_variants")])

  row.names(matrix_cDNA_PRE)<-cDNA_PRE$REAL_TILE_Plus_carried_variants

  cat("matrix_cDNA_PRE\n")
  cat(str(matrix_cDNA_PRE))
  cat("\n")

  rm(cDNA_PRE)
  
  
  #### ANNOT dfs ----
  
  colnames_counts<-colnames(matrix_gDNA_PRE)
  
  
  
  cat("colnames_counts_0\n")
  cat(str(colnames_counts))
  cat("\n")
  
  
  df_bc_batch<-gsub("_.+$","",colnames_counts)
  df_bc_bc<-gsub("^[^_]+_","",colnames_counts)
  
  
  
  df_bc_condition<-as.data.frame(cbind(df_bc_batch,df_bc_bc),
                                 stringsAsFactors=F)
  
  colnames(df_bc_condition)<-c("batch","bc")
  
  df_bc_condition$condition<-"NA"
  
  cat("df_bc_condition_0\n")
  cat(str(df_bc_condition))
  cat("\n")
  
  REF_bc<-paste("bc",seq(1,15,by=1),sep='')
  ALT_bc<-paste("bc",seq(16,30,by=1),sep='')
  
  df_bc_condition$condition[which(df_bc_condition$bc%in%REF_bc)]<-"REF"
  df_bc_condition$condition[which(df_bc_condition$bc%in%ALT_bc)]<-"ALT"
  
  
  Annot_gDNA<-merge(Deconvolve_table,
                    df_bc_condition,
                    by="batch")
  
  
  cat("Annot_gDNA_0\n")
  cat(str(Annot_gDNA))
  cat("\n")
  
  
 
  
  ###
  
  colnames_counts<-colnames(matrix_cDNA_PRE)
  
  
  
  cat("colnames_counts_0\n")
  cat(str(colnames_counts))
  cat("\n")
  
  
  df_bc_batch<-gsub("_.+$","",colnames_counts)
  df_bc_bc<-gsub("^[^_]+_","",colnames_counts)
  
  
  
  df_bc_condition<-as.data.frame(cbind(df_bc_batch,df_bc_bc),
                                 stringsAsFactors=F)
  
  colnames(df_bc_condition)<-c("batch","bc")
  
  df_bc_condition$condition<-"NA"
  
  cat("df_bc_condition_0\n")
  cat(str(df_bc_condition))
  cat("\n")
  
  REF_bc<-paste("bc",seq(1,15,by=1),sep='')
  ALT_bc<-paste("bc",seq(16,30,by=1),sep='')
  
  df_bc_condition$condition[which(df_bc_condition$bc%in%REF_bc)]<-"REF"
  df_bc_condition$condition[which(df_bc_condition$bc%in%ALT_bc)]<-"ALT"
  
  
  Annot_cDNA<-merge(Deconvolve_table,
                    df_bc_condition,
                    by="batch")
  
  
  cat("Annot_cDNA_0\n")
  cat(str(Annot_cDNA))
  cat("\n")
  
  
  #### indexes ctrls ----
  
  
  Rosetta_NCGR<-Rosetta[which(Rosetta$Label == 'Negative_Control_Genomic_Regions'),]
  
  cat("Rosetta_NCGR\n")
  cat(str(Rosetta_NCGR))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Rosetta_NCGR$KEY))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Rosetta_NCGR$KEY)))))
  cat("\n")
  
  indx.ctrls<-which(row.names(matrix_gDNA_PRE)%in%Rosetta_NCGR$REAL_TILE_Plus_carried_variants)
  
  
  
  check.matrix_gDNA_PRE<-matrix_gDNA_PRE[indx.ctrls,]
  
  cat("check.matrix_gDNA_PRE_0\n")
  cat(str(check.matrix_gDNA_PRE))
  cat("\n")
  
  
  check.matrix_cDNA_PRE<-matrix_cDNA_PRE[indx.ctrls,]
  
  cat("check.matrix_cDNA_PRE_0\n")
  cat(str(check.matrix_cDNA_PRE))
  cat("\n")
  
  
  FLAG_check_ctrls<-sum(row.names(check.matrix_gDNA_PRE) == row.names(check.matrix_cDNA_PRE))
  
  cat("FLAG_check_ctrls_0\n")
  cat(str(FLAG_check_ctrls))
  cat("\n")
  
  
  if(FLAG_check_ctrls == length(row.names(check.matrix_gDNA_PRE)))
  {
    
    
  }else{
    
    cat(sprintf("Inconcordant ctrls\n"))
    
    quit(status=1)
  }
  
  
  
  # quit(status = 1)
  
  #### MPRA normalization by UQ ----
  
  obj <- MpraObject(dnaCounts = matrix_gDNA_PRE, 
                    rnaCounts = matrix_cDNA_PRE, 
                    dnaAnnot = Annot_gDNA,
                    rnaAnnot = Annot_cDNA,
                    controls = indx.ctrls)
  
  cat("obj_1\n")
  cat(str(obj))
  cat("\n")
  
  obj <- estimateDepthFactors(obj, 
                              lib.factor = "batch", 
                              which.lib = "both",
                              depth.estimator = "uq")
  cat("obj_2\n")
  cat(str(obj))
  cat("\n")
  
  dnaDepth<-as.numeric(obj@dnaDepth)
  
  cat("dnaDepth\n")
  cat(str(dnaDepth))
  cat("\n")
  
  dnaDepth_df<-as.data.frame(cbind(colnames(matrix_gDNA_PRE),dnaDepth), stringsAsFactors=F)
  
  colnames(dnaDepth_df)<-c("matrix_colnames","dnaDepth")
  
  cat("dnaDepth_df\n")
  cat(str(dnaDepth_df))
  cat("\n")
  
  
  
  cat("dnaDepth\n")
  cat(sprintf(as.character(names(summary(as.factor(dnaDepth))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(dnaDepth)))))
  cat("\n")
  
  rnaDepth<-as.numeric(obj@rnaDepth)
  
  cat("rnaDepth\n")
  cat(str(rnaDepth))
  cat("\n")
  
  rnaDepth_df<-as.data.frame(cbind(colnames(matrix_cDNA_PRE),rnaDepth), stringsAsFactors=F)
  
  colnames(rnaDepth_df)<-c("matrix_colnames","rnaDepth")
  
  cat("rnaDepth_df\n")
  cat(str(rnaDepth_df))
  cat("\n")
  
  cat("rnaDepth\n")
  cat(sprintf(as.character(names(summary(as.factor(rnaDepth))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(rnaDepth)))))
  cat("\n")
  
  # quit(status=1)
  
  
  
  
  ####  matrix NOR ----
  
  
  colnames_gDNA<-colnames(matrix_gDNA_PRE)
  
  cat("colnames_gDNA\n")
  cat(str(colnames_gDNA))
  cat("\n")
  
  
  matrix_gDNA_NOR<-t(t(matrix_gDNA_PRE) / dnaDepth)
  
  cat("matrix_gDNA_NOR\n")
  cat(str(matrix_gDNA_NOR))
  cat("\n")
  
  matrix_gDNA_NOR<-round(matrix_gDNA_NOR,1)
  
  cat("matrix_gDNA_NOR\n")
  cat(str(matrix_gDNA_NOR))
  cat("\n")
  
  colnames_cDNA<-colnames(matrix_cDNA_PRE)
  
  cat("colnames_cDNA\n")
  cat(str(colnames_cDNA))
  cat("\n")
  
  matrix_cDNA_NOR<-t(t(matrix_cDNA_PRE) / rnaDepth)
  
  cat("matrix_cDNA_NOR\n")
  cat(str(matrix_cDNA_NOR))
  cat("\n")
  
  matrix_cDNA_NOR<-round(matrix_cDNA_NOR,1)
  
  cat("matrix_cDNA_NOR\n")
  cat(str(matrix_cDNA_NOR))
  cat("\n")
  
  
  #### SAVE ----
  
  setwd(out)
  
  filename_1<-paste("matrix_gDNA_NOR",".rds", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  
  matrix_gDNA_NOR<-as.data.frame(matrix_gDNA_NOR)
  
  matrix_gDNA_NOR$REAL_TILE_Plus_carried_variants<-row.names(matrix_gDNA_NOR)
  
  
  
  #write.table(matrix_gDNA_NOR,file=filename_1,sep="\t",row.names = F,quote=F)
  
  saveRDS(matrix_gDNA_NOR,file=filename_1)#,sep="\t",row.names = F,quote=F)
  
  
  
  filename_1<-paste("matrix_cDNA_NOR",".rds", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  
  matrix_cDNA_NOR<-as.data.frame(matrix_cDNA_NOR)
  
  matrix_cDNA_NOR$REAL_TILE_Plus_carried_variants<-row.names(matrix_cDNA_NOR)
  
  
  
  #write.table(matrix_cDNA_NOR,file=filename_1,sep="\t",row.names = F,quote=F)
  
  saveRDS(matrix_cDNA_NOR,file=filename_1)#,sep="\t",row.names = F,quote=F)
  
  
 
  filename_1<-paste("Annot_gDNA_df",".rds", sep='')

  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  
  saveRDS(Annot_gDNA,file=filename_1)#,sep="\t",row.names = F,quote=F)
  
  
  filename_1<-paste("Annot_cDNA_df",".rds", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  
  saveRDS(Annot_cDNA,file=filename_1)#,sep="\t",row.names = F,quote=F)
  
  filename_1<-paste("dnaDepth_df",".tsv", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  
  write.table(dnaDepth_df,file=filename_1,sep="\t",row.names = F,quote=F)
  
  filename_1<-paste("rnaDepth_df",".tsv", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  
  write.table(rnaDepth_df,file=filename_1,sep="\t",row.names = F,quote=F)
  
  
}
parameters_FC_Vockley_calculations = function(option_list)
{
  
  library("gtools", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/")
  
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
  
  #### READ data ----
  
  
  setwd(out)
  
  filename=paste('Rosetta_extended','.rds',sep='')
  Rosetta_extended <-readRDS(file=filename)
  
  
  cat("Rosetta_extended_0\n")
  cat(str(Rosetta_extended))
  cat("\n")
  cat(sprintf(as.character(names(summary(Rosetta_extended$batch)))))
  cat("\n")
  cat(sprintf(as.character(summary(Rosetta_extended$batch))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Rosetta_extended$type)))))
  cat("\n")
  cat(sprintf(as.character(summary(Rosetta_extended$type))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Rosetta_extended$condition)))))
  cat("\n")
  cat(sprintf(as.character(summary(Rosetta_extended$condition))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Rosetta_extended$master_sample)))))
  cat("\n")
  cat(sprintf(as.character(summary(Rosetta_extended$master_sample))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Rosetta_extended$bc)))))
  cat("\n")
  cat(sprintf(as.character(summary(Rosetta_extended$bc))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Rosetta_extended$Cell_Type)))))
  cat("\n")
  cat(sprintf(as.character(summary(Rosetta_extended$Cell_Type))))
  cat("\n")
  
  Deconvolve_table<-unique(Rosetta_extended[,c(which(colnames(Rosetta_extended) == "master_sample"),
                                               which(colnames(Rosetta_extended) == "batch"))])
  
  cat("Deconvolve_table_0\n")
  cat(str(Deconvolve_table))
  cat("\n")
  
  
  filename=paste('Rosetta_df','.rds',sep='')
  Rosetta <-readRDS(file=filename)
  
  
  cat("Rosetta_0\n")
  cat(str(Rosetta))
  cat("\n")
  
  rm(Rosetta_extended)
  
  
  filename=paste('matrix_gDNA_PRE','.rds',sep='')
  gDNA_PRE <-readRDS(file=filename)
  
  cat("gDNA_PRE\n")
  cat(str(gDNA_PRE))
  cat("\n")
  
  # matrix_gDNA_PRE<-as.matrix(gDNA_PRE[,-which(colnames(gDNA_PRE) == "REAL_TILE_Plus_carried_variants")])
  # 
  # row.names(matrix_gDNA_PRE)<-gDNA_PRE$REAL_TILE_Plus_carried_variants
  # 
  # cat("matrix_gDNA_PRE\n")
  # cat(str(matrix_gDNA_PRE))
  # cat("\n")
  # 
  # rm(gDNA_PRE)
  
  filename=paste('matrix_cDNA_PRE','.rds',sep='')
  #cDNA_PRE <-as.data.frame(fread(file=filename, sep="\t",header=T), stringsAsFactors = F)
  cDNA_PRE <-readRDS(file=filename)
  
  
  cat("cDNA_PRE\n")
  cat(str(cDNA_PRE))
  cat("\n")
  
  # matrix_cDNA_PRE<-as.matrix(cDNA_PRE[,-which(colnames(cDNA_PRE) == "REAL_TILE_Plus_carried_variants")])
  # 
  # row.names(matrix_cDNA_PRE)<-cDNA_PRE$REAL_TILE_Plus_carried_variants
  # 
  # cat("matrix_cDNA_PRE\n")
  # cat(str(matrix_cDNA_PRE))
  # cat("\n")
  # 
  # rm(cDNA_PRE)
  
  
  
  
  filename=paste('matrix_gDNA_NOR','.rds',sep='')
  gDNA_NOR <-readRDS(file=filename)
  
  cat("gDNA_NOR\n")
  cat(str(gDNA_NOR))
  cat("\n")
  
  # matrix_gDNA_NOR<-as.matrix(gDNA_NOR[,-which(colnames(gDNA_NOR) == "REAL_TILE_Plus_carried_variants")])
  # 
  # row.names(matrix_gDNA_NOR)<-gDNA_NOR$REAL_TILE_Plus_carried_variants
  # 
  # cat("matrix_gDNA_NOR\n")
  # cat(str(matrix_gDNA_NOR))
  # cat("\n")
  # 
  # rm(gDNA_NOR)
  
  filename=paste('matrix_cDNA_NOR','.rds',sep='')
  #cDNA_NOR <-as.data.frame(fread(file=filename, sep="\t",header=T), stringsAsFactors = F)
  cDNA_NOR <-readRDS(file=filename)
  
  
  cat("cDNA_NOR\n")
  cat(str(cDNA_NOR))
  cat("\n")
  
  # matrix_cDNA_NOR<-as.matrix(cDNA_NOR[,-which(colnames(cDNA_NOR) == "REAL_TILE_Plus_carried_variants")])
  # 
  # row.names(matrix_cDNA_NOR)<-cDNA_NOR$REAL_TILE_Plus_carried_variants
  # 
  # cat("matrix_cDNA_NOR\n")
  # cat(str(matrix_cDNA_NOR))
  # cat("\n")
  # 
  # rm(cDNA_NOR)
  
 
  
  
  
  
  filename_1<-paste("Annot_gDNA_df",".rds", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  Annot_gDNA_df<-readRDS(file=filename_1)#,sep="\t")#,header=T), stringsAsFactors=F)
  
  cat("Annot_gDNA_df\n")
  cat(str(Annot_gDNA_df))
  cat("\n")
  
  filename_1<-paste("Annot_cDNA_df",".rds", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  Annot_cDNA_df<-readRDS(file=filename_1)#,sep="\t",header=T), stringsAsFactors=F)
  
  cat("Annot_cDNA_df\n")
  cat(str(Annot_cDNA_df))
  cat("\n")
  
  
  filename_1<-paste("dnaDepth_df",".tsv", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  dnaDepth_df<-as.data.frame(fread(file=filename_1,sep="\t",header=T), stringsAsFactors=F)
  
  cat("dnaDepth_df\n")
  cat(str(dnaDepth_df))
  cat("\n")
  
  filename_1<-paste("rnaDepth_df",".tsv", sep='')
  
  # cat("filename_1\n")
  # cat(sprintf(as.character(filename_1)))
  # cat("\n")
  
  rnaDepth_df<-as.data.frame(fread(file=filename_1,sep="\t",header=T), stringsAsFactors=F)
  
  cat("rnaDepth_df\n")
  cat(str(rnaDepth_df))
  cat("\n")
  
  
  
  
  
  
 
  
  
  
  ##### Annot_gDNA_df and Annot_cDNA_df MERGE and ordered levels ----
  
 
  
 
  Annot_gDNA_df$Cell_Type<-gsub("_.+$","",Annot_gDNA_df$master_sample)
  Annot_gDNA_df$Replicate<-gsub("^[^_]+_","",Annot_gDNA_df$master_sample)
  Annot_gDNA_df$type<-"gDNA"
  
  
  Annot_cDNA_df$Cell_Type<-gsub("_.+$","",Annot_cDNA_df$master_sample)
  Annot_cDNA_df$Replicate<-gsub("^[^_]+_","",Annot_cDNA_df$master_sample)
  Annot_cDNA_df$type<-"cDNA"
  
  Annot_DEF<-rbind(Annot_gDNA_df,Annot_cDNA_df)
  
  
  Annot_DEF$type<-factor(Annot_DEF$type,
                             c("gDNA","cDNA"),
                             ordered=T)
  Annot_DEF$Cell_Type<-factor(Annot_DEF$Cell_Type,
                                  c("K562","CHRF","HL60","THP1"),
                                  ordered=T)
  # Annot_DEF$Replicate<-factor(Annot_DEF$Replicate,
  #                                 c("Rep1","Rep2","Rep3","Rep4","Rep5","Rep6","Rep7","Rep8","Rep9"),
  #                                 ordered=T)
  
  
  Annot_DEF<-Annot_DEF[order(Annot_DEF$Cell_Type,Annot_DEF$master_sample,Annot_DEF$type),]
  
  
  
  levels_master_samples<-levels(Annot_DEF$master_sample)
  # 
  # Annot_DEF$master_sample<-factor(Annot_DEF$master_sample,
  #                          levels=levels_master_samples,
  #                          ordered=T)
  # 
  # levels_batchs<-unique(as.character(Annot_DEF$batch))
  # 
  # Annot_DEF$batch<-factor(Annot_DEF$batch,
  #                                 levels=levels_batchs,
  #                                 ordered=T)
  
  cat("Annot_DEF_2\n")
  cat(str(Annot_DEF))
  cat("\n")
  
  Annot_DEF_subset<-unique(Annot_DEF[,c(which(colnames(Annot_DEF) == "Cell_Type"),
                                        which(colnames(Annot_DEF) == "master_sample"),
                                        which(colnames(Annot_DEF) == "batch"),
                                              which(colnames(Annot_DEF) == "type"),
                                                    which(colnames(Annot_DEF) == "Replicate"))])
  
  cat("Annot_DEF_subset\n")
  cat(str(Annot_DEF_subset))
  cat("\n")
  
  levels_batchs<-levels(Annot_DEF_subset$batch)
  
  cat("levels_batchs\n")
  cat(str(levels_batchs))
  cat("\n")
  
  #### melt of matrixes ----
  
  gDNA_NOR.m<-melt(gDNA_NOR, id.variables="REAL_TILE_Plus_carried_variants")
  
  cat("gDNA_NOR.m_0\n")
  cat(str(gDNA_NOR.m))
  cat("\n")
  
  #gDNA_NOR.m$variable<-gsub("^aggregated_counts_","",gDNA_NOR.m$variable)
  
  
  
  
  
  gDNA_NOR.m$batch<-gsub("_.+$","",gDNA_NOR.m$variable)
  gDNA_NOR.m$bc<-gsub("^[^_]+_","",gDNA_NOR.m$variable)
  
  gDNA_NOR.m$GROUP<-"NOR"
  gDNA_NOR.m$type<-"gDNA"
  
  
  
  gDNA_NOR.m<-gDNA_NOR.m[,-which(colnames(gDNA_NOR.m) == "variable")]
  
  cat("gDNA_NOR.m_1\n")
  cat(str(gDNA_NOR.m))
  cat("\n")
  
  cDNA_NOR.m<-melt(cDNA_NOR, id.variables="REAL_TILE_Plus_carried_variants")
  
  cat("cDNA_NOR.m_0\n")
  cat(str(cDNA_NOR.m))
  cat("\n")
  
  #cDNA_NOR.m$variable<-gsub("^aggregated_counts_","",cDNA_NOR.m$variable)
  
  
  cDNA_NOR.m$batch<-gsub("_.+$","",cDNA_NOR.m$variable)
  cDNA_NOR.m$bc<-gsub("^[^_]+_","",cDNA_NOR.m$variable)
  
  cDNA_NOR.m$GROUP<-"NOR"
  cDNA_NOR.m$type<-"cDNA"
  
  
  
  cDNA_NOR.m<-cDNA_NOR.m[,-which(colnames(cDNA_NOR.m) == "variable")]
  
  cat("cDNA_NOR.m_1\n")
  cat(str(cDNA_NOR.m))
  cat("\n")
  
  ###
  
  gDNA_PRE.m<-melt(gDNA_PRE, id.variables="REAL_TILE_Plus_carried_variants")
  
  cat("gDNA_PRE.m_0\n")
  cat(str(gDNA_PRE.m))
  cat("\n")
  
 # gDNA_PRE.m$variable<-gsub("^aggregated_counts_","",gDNA_PRE.m$variable)
  
  
  gDNA_PRE.m$batch<-gsub("_.+$","",gDNA_PRE.m$variable)
  gDNA_PRE.m$bc<-gsub("^[^_]+_","",gDNA_PRE.m$variable)
  
  gDNA_PRE.m$GROUP<-"PRE"
  gDNA_PRE.m$type<-"gDNA"
  
  
  
  gDNA_PRE.m<-gDNA_PRE.m[,-which(colnames(gDNA_PRE.m) == "variable")]
  
  cat("gDNA_PRE.m_1\n")
  cat(str(gDNA_PRE.m))
  cat("\n")
  
  cDNA_PRE.m<-melt(cDNA_PRE, id.variables="REAL_TILE_Plus_carried_variants")
  
  cat("cDNA_PRE.m_0\n")
  cat(str(cDNA_PRE.m))
  cat("\n")
  
 # cDNA_PRE.m$variable<-gsub("^aggregated_counts_","",cDNA_PRE.m$variable)
  
  
  cDNA_PRE.m$batch<-gsub("_.+$","",cDNA_PRE.m$variable)
  cDNA_PRE.m$bc<-gsub("^[^_]+_","",cDNA_PRE.m$variable)
  
  cDNA_PRE.m$GROUP<-"PRE"
  cDNA_PRE.m$type<-"cDNA"
  
  
  
  cDNA_PRE.m<-cDNA_PRE.m[,-which(colnames(cDNA_PRE.m) == "variable")]
  
  cat("cDNA_PRE.m_1\n")
  cat(str(cDNA_PRE.m))
  cat("\n")
  
  
  
  
  #### rbind of matrixes ----
  
  
  REP_df<-rbind(gDNA_NOR.m,cDNA_NOR.m,gDNA_PRE.m,cDNA_PRE.m)
  
  
  REP_df$type<-factor(REP_df$type,
                       levels=c("gDNA","cDNA"),
                       ordered=T)
  
  REP_df$GROUP<-factor(REP_df$GROUP,
                             levels=c("PRE","NOR"),
                             ordered=T)
  
  
  REP_df$batch<-factor(REP_df$batch,
                       levels=levels_batchs,
                       ordered=T)
  
  cat("REP_df_1\n")
  cat(str(REP_df))
  cat("\n")
  
  #### merge REP with Annot_DEF ----
  
  
  REP_df<-merge(REP_df,
                Annot_DEF,
                by=c("batch","bc","type"),
                all.x=T)
  
  cat("REP_df_2\n")
  cat(str(REP_df))
  cat("\n")
  
  REP_df<-REP_df[order(REP_df$master_sample,REP_df$type),]
  
  
  REP_df$sample<-paste(REP_df$master_sample,REP_df$type,sep="_")
  
  
  levels_sample<-unique(as.character(REP_df$sample))
  
  cat("levels_sample\n")
  cat(str(levels_sample))
  cat("\n")
  
  REP_df$sample<-factor(REP_df$sample,
                        levels=levels_sample,
                        ordered=T)
  
  cat("REP_df_3\n")
  cat(str(REP_df))
  cat("\n")
  
  
  #### interaction sample and group ----
  
  REP_df$interaction_1<-interaction(REP_df$GROUP,REP_df$sample, sep="|",lex.order = T)
  
  
  cat("REP_df_3\n")
  cat(str(REP_df))
  cat("\n")
  cat("interaction_1\n")
  cat(sprintf(as.character(names(summary(REP_df$interaction_1)))))
  cat("\n")
  cat(sprintf(as.character(summary(REP_df$interaction_1))))
  cat("\n")
  
  
  REP_df<-REP_df[order(REP_df$interaction_1),]
  
  cat("REP_df_4\n")
  cat(str(REP_df))
  cat("\n")
  cat("interaction_1\n")
  cat(sprintf(as.character(names(summary(REP_df$interaction_1)))))
  cat("\n")
  cat(sprintf(as.character(summary(REP_df$interaction_1))))
  cat("\n")
  
  #### LISTS DEFINITION ----
  
  List_values<-list()
  List_medians<-list()
  
  #### POOL for cDNA and gDNA calculation -----
  
 
  REP_df.dt<-data.table(REP_df,
                        key=c("REAL_TILE_Plus_carried_variants","interaction_1"))
  
  
  cat("----------------------------------------------------------------------------------------------------------------->REP_df.dt\n")
  cat(str(REP_df.dt))
  cat("\n")
  
  
  dna_rna_pool_df<-unique(as.data.frame(REP_df.dt[,.(Median_value=median(value),
                                                     type=type,
                                                   Cell_Type=Cell_Type,
                                                   GROUP=GROUP,
                                                   batch=batch,
                                                   master_sample=master_sample,
                                                   Replicate=Replicate),by=key(REP_df.dt)],stringsAsFactors=F))
  
  
  cat("dna_rna_pool_df_1\n")
  cat(str(dna_rna_pool_df))
  cat("\n")
  
  dna_rna_pool_df_NOR<-dna_rna_pool_df[which(dna_rna_pool_df$GROUP =="NOR"),-which(colnames(dna_rna_pool_df) == "interaction_1")]
  
  
  # cat("dna_rna_pool_df_NOR_1\n")
  # cat(str(dna_rna_pool_df_NOR))
  # cat("\n")
  
  dna_rna_pool_df_NOR<-droplevels(dna_rna_pool_df_NOR)
  
  # cat("dna_rna_pool_df_NOR_2\n")
  # cat(str(dna_rna_pool_df_NOR))
  # cat("\n")
  
  
  dna_pool_df_NOR<-dna_rna_pool_df_NOR[which(dna_rna_pool_df_NOR$type == "gDNA"),]
  
  # cat("dna_pool_df_NOR_1\n")
  # cat(str(dna_pool_df_NOR))
  # cat("\n")
  
  dna_pool_df_NOR<-droplevels(dna_pool_df_NOR)
  
  cat("dna_pool_df_NOR_2\n")
  cat(str(dna_pool_df_NOR))
  cat("\n")
  
  List_values[['dna_pool']]<-dna_pool_df_NOR
  
  
  dna_pool_df_NOR.dt<-data.table(dna_pool_df_NOR,
                        key=c("REAL_TILE_Plus_carried_variants","Cell_Type"))
  
  dna_pool_median<-unique(as.data.frame(dna_pool_df_NOR.dt[,.(median_Median_value=median(Median_value)),
                                                           by=key(dna_pool_df_NOR.dt)],stringsAsFactors=F))
  
  
  cat("dna_pool_median_0\n")
  cat(str(dna_pool_median))
  cat("\n")
  
  List_medians[['dna_pool']]<-dna_pool_median
  
  
  dna_pool_df_NOR_wide<-pivot_wider(dna_pool_df_NOR,
                                          id_cols=REAL_TILE_Plus_carried_variants,
                                          names_from=master_sample,
                                          values_from=Median_value)
  
  cat("dna_pool_df_NOR_wide\n")
  cat(str(dna_pool_df_NOR_wide))
  cat("\n")
  
  dna_pool_df_NOR_wide.matrix<-as.matrix(dna_pool_df_NOR_wide[,-1])
  
  row.names(dna_pool_df_NOR_wide.matrix)<-dna_pool_df_NOR_wide$REAL_TILE_Plus_carried_variants
  
  cat("dna_pool_df_NOR_wide.matrix\n")
  cat(str(dna_pool_df_NOR_wide.matrix))
  cat("\n")
  
  
  rna_pool_df_NOR<-dna_rna_pool_df_NOR[which(dna_rna_pool_df_NOR$type == "cDNA"),]
  
  # cat("rna_pool_df_NOR_1\n")
  # cat(str(rna_pool_df_NOR))
  # cat("\n")
  
  rna_pool_df_NOR<-droplevels(rna_pool_df_NOR)
  
  cat("rna_pool_df_NOR_2\n")
  cat(str(rna_pool_df_NOR))
  cat("\n")
  
  List_values[['rna_pool']]<-rna_pool_df_NOR
  
  rna_pool_df_NOR.dt<-data.table(rna_pool_df_NOR,
                                 key=c("REAL_TILE_Plus_carried_variants","Cell_Type"))
  
  rna_pool_median<-unique(as.data.frame(rna_pool_df_NOR.dt[,.(median_Median_value=median(Median_value)),
                                                           by=key(rna_pool_df_NOR.dt)],stringsAsFactors=F))
  
  
  cat("rna_pool_median_0\n")
  cat(str(rna_pool_median))
  cat("\n")
  
  List_medians[['rna_pool']]<-rna_pool_median
  
  
  rna_pool_df_NOR_wide<-pivot_wider(rna_pool_df_NOR,
                                    id_cols=REAL_TILE_Plus_carried_variants,
                                    names_from=master_sample,
                                    values_from=Median_value)
  
  cat("rna_pool_df_NOR_wide\n")
  cat(str(rna_pool_df_NOR_wide))
  cat("\n")
  
  rna_pool_df_NOR_wide.matrix<-as.matrix(rna_pool_df_NOR_wide[,-1])
  
  row.names(rna_pool_df_NOR_wide.matrix)<-rna_pool_df_NOR_wide$REAL_TILE_Plus_carried_variants
  
  cat("rna_pool_df_NOR_wide.matrix\n")
  cat(str(rna_pool_df_NOR_wide.matrix))
  cat("\n")
  
  #### LogFC ----
  
  
  LogFC_matrix<- foldchange2logratio(foldchange(rna_pool_df_NOR_wide.matrix+0.0000001,dna_pool_df_NOR_wide.matrix+0.0000001))
  
  cat("LogFC_matrix\n")
  cat(str(LogFC_matrix))
  cat("\n")
  
  LogFC_df<-as.data.frame(LogFC_matrix, stringsAsFactors=F)
  
  
  LogFC_df$REAL_TILE_Plus_carried_variants<-row.names(LogFC_df)
  
  cat("LogFC_df\n")
  cat(str(LogFC_df))
  cat("\n")
  
  LogFC_df.m<-melt(LogFC_df, id.cols="REAL_TILE_Plus_carried_variants")
  
  cat("LogFC_df.m_1\n")
  cat(str(LogFC_df.m))
  cat("\n")
  
  colnames(LogFC_df.m)[which(colnames(LogFC_df.m) == "variable")]<-"master_sample"
  colnames(LogFC_df.m)[which(colnames(LogFC_df.m) == "value")]<-"LogFC"
  
  
  LogFC_df.m$master_sample<-factor(LogFC_df.m$master_sample,
                                   levels=levels_master_samples,ordered=T)
  
  
  cat("LogFC_df.m_2\n")
  cat(str(LogFC_df.m))
  cat("\n")
  
  LogFC_df.m<-merge(LogFC_df.m,
                    Annot_DEF_subset,
                    by="master_sample")
  
  
  cat("LogFC_df.m_3\n")
  cat(str(LogFC_df.m))
  cat("\n")
  
  
  ### NO NA
  
  LogFC_df.m_NO_NA<-LogFC_df.m[!is.na(LogFC_df.m$LogFC),]
  
  
  cat("LogFC_df.m_NO_NA_1\n")
  cat(str(LogFC_df.m_NO_NA))
  cat("\n")
  
  ### Infinite reverted
  
  LogFC_df.m_NO_NA$LogFC[(is.infinite(LogFC_df.m_NO_NA$LogFC) &  LogFC_df.m_NO_NA$LogFC > 0)]<-60
  # LogFC_df.m_NO_NA<-LogFC_df.m_NO_NA[-(is.infinite(LogFC_df.m_NO_NA$LogFC) &  LogFC_df.m_NO_NA$LogFC < 0),]
  #LogFC_df.m_NO_NA$LogFC[(is.infinite(LogFC_df.m_NO_NA$LogFC) &  LogFC_df.m_NO_NA$LogFC < 0)]<--1*10000
  
  LogFC_df.m_NO_NA$GROUP<-"NOR"
  
  cat("LogFC_df.m_NO_NA_2_infinite_reverted\n")
  cat(str(LogFC_df.m_NO_NA))
  cat("\n")
  
  
  


  A<-summary(LogFC_df.m_NO_NA$LogFC)
  
  cat("summary_LogFC\n")
  cat(sprintf(as.character(names(A))))
  cat("\n")
  cat(sprintf(as.character(A)))
  cat("\n")
  
  
  check<-LogFC_df.m_NO_NA[which(LogFC_df.m_NO_NA$LogFC >0),]
  
  cat("check_LogFC\n")
  cat(str(check))
  cat("\n")
  cat(sprintf(as.character(names(summary(check$master_sample)))))
  cat("\n")
  cat(sprintf(as.character(summary(check$master_sample))))
  cat("\n")
  
  A<-summary(check$LogFC)
  
  cat("summary_LogFC\n")
  cat(sprintf(as.character(names(A))))
  cat("\n")
  cat(sprintf(as.character(A)))
  cat("\n")
  
  
  check<-LogFC_df.m_NO_NA[is.na(LogFC_df.m_NO_NA$LogFC),]
  
  cat("check_2\n")
  cat(str(check))
  cat("\n")
  
  
  List_values[['LogFC']]<-LogFC_df.m_NO_NA
  
  LogFC_df.m_NO_NA.dt<-data.table(LogFC_df.m_NO_NA,
                                 key=c("REAL_TILE_Plus_carried_variants","Cell_Type"))
  
  LogFC_median<-unique(as.data.frame(LogFC_df.m_NO_NA.dt[,.(median_LogFC=median(LogFC)),
                                                           by=key(LogFC_df.m_NO_NA.dt)],stringsAsFactors=F))
  
  
  cat("LogFC_median_0\n")
  cat(str(LogFC_median))
  cat("\n")
  
  List_medians[['LogFC']]<-LogFC_median
  
  
  check_LogFC<-LogFC_median[which(LogFC_median$median_LogFC >0),]
  
  cat("check_LogFC_LogFC\n")
  cat(str(check_LogFC))
  cat("\n")
  cat(sprintf(as.character(names(summary(check_LogFC$Cell_Type)))))
  cat("\n")
  cat(sprintf(as.character(summary(check_LogFC$Cell_Type))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(check_LogFC$REAL_TILE_Plus_carried_variants))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(check_LogFC$REAL_TILE_Plus_carried_variants)))))
  cat("\n")
  
  A<-summary(check_LogFC$median_LogFC)
  
  cat("summary_median_LogFC\n")
  cat(sprintf(as.character(names(A))))
  cat("\n")
  cat(sprintf(as.character(A)))
  cat("\n")
 
 
  #### REF and ALT dna and rna ----
  
  REP_df.dt<-data.table(REP_df,
                        key=c("REAL_TILE_Plus_carried_variants","interaction_1","condition"))
  
  
  cat("-------------------------------------------------------------------------------------------------------------------------->condition\n")
  cat(str(REP_df.dt))
  cat("\n")
  
  
  dna_rna_condition_df<-unique(as.data.frame(REP_df.dt[,.(Median_value=median(value),
                                                          type=type,
                                                          Cell_Type=Cell_Type,
                                                          GROUP=GROUP,
                                                          batch=batch,
                                                          master_sample=master_sample,
                                              Replicate=Replicate),by=key(REP_df.dt)],stringsAsFactors=F))
  
  
  cat("dna_rna_condition_df_1\n")
  cat(str(dna_rna_condition_df))
  cat("\n")
  
  ### REF
  
  
  dna_rna_condition_df_NOR_REF<-unique(dna_rna_condition_df[which(dna_rna_condition_df$GROUP =="NOR" &
                                                             dna_rna_condition_df$condition == "REF"),-c(which(colnames(dna_rna_condition_df) == "interaction_1"),
                                                                                                         which(colnames(dna_rna_condition_df) == "condition"))])
  
  
  # cat("dna_rna_condition_df_NOR_REF_1\n")
  # cat(str(dna_rna_condition_df_NOR_REF))
  # cat("\n")
  
  dna_rna_condition_df_NOR_REF<-droplevels(dna_rna_condition_df_NOR_REF)
  
  cat("dna_rna_condition_df_NOR_REF_2\n")
  cat(str(dna_rna_condition_df_NOR_REF))
  cat("\n")
  
  
  
  dna_condition_df_NOR_REF<-dna_rna_condition_df_NOR_REF[which(dna_rna_condition_df_NOR_REF$type == "gDNA"),]
  
  # cat("dna_condition_df_NOR_REF_1\n")
  # cat(str(dna_condition_df_NOR_REF))
  # cat("\n")
  
  dna_condition_df_NOR_REF<-droplevels(dna_condition_df_NOR_REF)
  
  cat("dna_condition_df_NOR_REF_2\n")
  cat(str(dna_condition_df_NOR_REF))
  cat("\n")
  
  List_values[['dna_REF']]<-dna_condition_df_NOR_REF
  
  dna_condition_df_NOR_REF.dt<-data.table(dna_condition_df_NOR_REF,
                                 key=c("REAL_TILE_Plus_carried_variants","Cell_Type"))
  
  dna_REF_median<-unique(as.data.frame(dna_condition_df_NOR_REF.dt[,.(median_Median_value=median(Median_value)),
                                                           by=key(dna_condition_df_NOR_REF.dt)],stringsAsFactors=F))
  
  
  cat("dna_REF_median_0\n")
  cat(str(dna_REF_median))
  cat("\n")
  
  List_medians[['dna_REF']]<-dna_REF_median
  
  
  dna_condition_df_NOR_REF_wide<-pivot_wider(dna_condition_df_NOR_REF,
                                    id_cols=REAL_TILE_Plus_carried_variants,
                                    names_from=master_sample,
                                    values_from=Median_value)
  
  cat("dna_condition_df_NOR_REF_wide\n")
  cat(str(dna_condition_df_NOR_REF_wide))
  cat("\n")
  
  dna_condition_df_NOR_REF_wide.matrix<-as.matrix(dna_condition_df_NOR_REF_wide[,-1])
  
  row.names(dna_condition_df_NOR_REF_wide.matrix)<-dna_condition_df_NOR_REF_wide$REAL_TILE_Plus_carried_variants
  
  cat("dna_condition_df_NOR_REF_wide.matrix\n")
  cat(str(dna_condition_df_NOR_REF_wide.matrix))
  cat("\n")
  
  
  
  
  rna_condition_df_NOR_REF<-dna_rna_condition_df_NOR_REF[which(dna_rna_condition_df_NOR_REF$type == "cDNA"),]
  
  # cat("rna_condition_df_NOR_REF_1\n")
  # cat(str(rna_condition_df_NOR_REF))
  # cat("\n")
  
  rna_condition_df_NOR_REF<-droplevels(rna_condition_df_NOR_REF)
  
  cat("rna_condition_df_NOR_REF_2\n")
  cat(str(rna_condition_df_NOR_REF))
  cat("\n")
  
  List_values[['rna_REF']]<-rna_condition_df_NOR_REF
  
  rna_condition_df_NOR_REF.dt<-data.table(rna_condition_df_NOR_REF,
                                          key=c("REAL_TILE_Plus_carried_variants","Cell_Type"))
  
  rna_REF_median<-unique(as.data.frame(rna_condition_df_NOR_REF.dt[,.(median_Median_value=median(Median_value)),
                                                                   by=key(rna_condition_df_NOR_REF.dt)],stringsAsFactors=F))
  
  
  cat("rna_REF_median_0\n")
  cat(str(rna_REF_median))
  cat("\n")
  
  List_medians[['rna_REF']]<-rna_REF_median
  
  
  rna_condition_df_NOR_REF_wide<-pivot_wider(rna_condition_df_NOR_REF,
                                             id_cols=REAL_TILE_Plus_carried_variants,
                                             names_from=master_sample,
                                             values_from=Median_value)
  
  cat("rna_condition_df_NOR_REF_wide\n")
  cat(str(rna_condition_df_NOR_REF_wide))
  cat("\n")
  
  rna_condition_df_NOR_REF_wide.matrix<-as.matrix(rna_condition_df_NOR_REF_wide[,-1])
  
  row.names(rna_condition_df_NOR_REF_wide.matrix)<-rna_condition_df_NOR_REF_wide$REAL_TILE_Plus_carried_variants
  
  cat("rna_condition_df_NOR_REF_wide.matrix\n")
  cat(str(rna_condition_df_NOR_REF_wide.matrix))
  cat("\n")
  
  
  
  ### ALT
  
  
  dna_rna_condition_df_NOR_ALT<-unique(dna_rna_condition_df[which(dna_rna_condition_df$GROUP =="NOR" &
                                                                    dna_rna_condition_df$condition == "ALT"),-c(which(colnames(dna_rna_condition_df) == "interaction_1"),
                                                                                                                which(colnames(dna_rna_condition_df) == "condition"))])
  
  
  # cat("dna_rna_condition_df_NOR_ALT_1\n")
  # cat(str(dna_rna_condition_df_NOR_ALT))
  # cat("\n")
  
  dna_rna_condition_df_NOR_ALT<-droplevels(dna_rna_condition_df_NOR_ALT)
  
  cat("dna_rna_condition_df_NOR_ALT_2\n")
  cat(str(dna_rna_condition_df_NOR_ALT))
  cat("\n")
  
  
  
  dna_condition_df_NOR_ALT<-dna_rna_condition_df_NOR_ALT[which(dna_rna_condition_df_NOR_ALT$type == "gDNA"),]
  
  # cat("dna_condition_df_NOR_ALT_1\n")
  # cat(str(dna_condition_df_NOR_ALT))
  # cat("\n")
  
  dna_condition_df_NOR_ALT<-droplevels(dna_condition_df_NOR_ALT)
  
  cat("dna_condition_df_NOR_ALT_2\n")
  cat(str(dna_condition_df_NOR_ALT))
  cat("\n")
  
  List_values[['dna_ALT']]<-dna_condition_df_NOR_ALT
  
  dna_condition_df_NOR_ALT.dt<-data.table(dna_condition_df_NOR_ALT,
                                          key=c("REAL_TILE_Plus_carried_variants","Cell_Type"))
  
  dna_ALT_median<-unique(as.data.frame(dna_condition_df_NOR_ALT.dt[,.(median_Median_value=median(Median_value)),
                                                                   by=key(dna_condition_df_NOR_ALT.dt)],stringsAsFactors=F))
  
  
  cat("dna_ALT_median_0\n")
  cat(str(dna_ALT_median))
  cat("\n")
  
  List_medians[['dna_ALT']]<-dna_ALT_median
  
  
  dna_condition_df_NOR_ALT_wide<-pivot_wider(dna_condition_df_NOR_ALT,
                                             id_cols=REAL_TILE_Plus_carried_variants,
                                             names_from=master_sample,
                                             values_from=Median_value)
  
  cat("dna_condition_df_NOR_ALT_wide\n")
  cat(str(dna_condition_df_NOR_ALT_wide))
  cat("\n")
  
  dna_condition_df_NOR_ALT_wide.matrix<-as.matrix(dna_condition_df_NOR_ALT_wide[,-1])
  
  row.names(dna_condition_df_NOR_ALT_wide.matrix)<-dna_condition_df_NOR_ALT_wide$REAL_TILE_Plus_carried_variants
  
  cat("dna_condition_df_NOR_ALT_wide.matrix\n")
  cat(str(dna_condition_df_NOR_ALT_wide.matrix))
  cat("\n")
  
  
  
  
  rna_condition_df_NOR_ALT<-dna_rna_condition_df_NOR_ALT[which(dna_rna_condition_df_NOR_ALT$type == "cDNA"),]
  
  # cat("rna_condition_df_NOR_ALT_1\n")
  # cat(str(rna_condition_df_NOR_ALT))
  # cat("\n")
  
  rna_condition_df_NOR_ALT<-droplevels(rna_condition_df_NOR_ALT)
  
  cat("rna_condition_df_NOR_ALT_2\n")
  cat(str(rna_condition_df_NOR_ALT))
  cat("\n")
  
  List_values[['rna_ALT']]<-rna_condition_df_NOR_ALT
  
  rna_condition_df_NOR_ALT.dt<-data.table(rna_condition_df_NOR_ALT,
                                          key=c("REAL_TILE_Plus_carried_variants","Cell_Type"))
  
  rna_ALT_median<-unique(as.data.frame(rna_condition_df_NOR_ALT.dt[,.(median_Median_value=median(Median_value)),
                                                                   by=key(rna_condition_df_NOR_ALT.dt)],stringsAsFactors=F))
  
  
  cat("rna_ALT_median_0\n")
  cat(str(rna_ALT_median))
  cat("\n")
  
  List_medians[['rna_ALT']]<-rna_ALT_median
  
  
  rna_condition_df_NOR_ALT_wide<-pivot_wider(rna_condition_df_NOR_ALT,
                                             id_cols=REAL_TILE_Plus_carried_variants,
                                             names_from=master_sample,
                                             values_from=Median_value)
  
  cat("rna_condition_df_NOR_ALT_wide\n")
  cat(str(rna_condition_df_NOR_ALT_wide))
  cat("\n")
  
  rna_condition_df_NOR_ALT_wide.matrix<-as.matrix(rna_condition_df_NOR_ALT_wide[,-1])
  
  row.names(rna_condition_df_NOR_ALT_wide.matrix)<-rna_condition_df_NOR_ALT_wide$REAL_TILE_Plus_carried_variants
  
  cat("rna_condition_df_NOR_ALT_wide.matrix\n")
  cat(str(rna_condition_df_NOR_ALT_wide.matrix))
  cat("\n")
  
  
  #### Vockley_REF ----
  

  
  Vockley_REF_matrix<-(rna_condition_df_NOR_REF_wide.matrix+0.0000001/dna_condition_df_NOR_REF_wide.matrix+0.0000001)/(rna_condition_df_NOR_ALT_wide.matrix+0.0000001/dna_condition_df_NOR_ALT_wide.matrix+0.0000001)
  
  cat("--------------------------------------------------------------------------------------------->Vockley_REF_matrix\n")
  cat(str(Vockley_REF_matrix))
  cat("\n")
  
  Vockley_REF_df<-as.data.frame(Vockley_REF_matrix, stringsAsFactors=F)
  
  
  Vockley_REF_df$REAL_TILE_Plus_carried_variants<-row.names(Vockley_REF_df)
  
  cat("Vockley_REF_df\n")
  cat(str(Vockley_REF_df))
  cat("\n")
  
  Vockley_REF_df.m<-melt(Vockley_REF_df, id.cols="REAL_TILE_Plus_carried_variants")
  
  cat("Vockley_REF_df.m_1\n")
  cat(str(Vockley_REF_df.m))
  cat("\n")
  
  colnames(Vockley_REF_df.m)[which(colnames(Vockley_REF_df.m) == "variable")]<-"master_sample"
  colnames(Vockley_REF_df.m)[which(colnames(Vockley_REF_df.m) == "value")]<-"Vockley_REF"
  
  
  Vockley_REF_df.m$master_sample<-factor(Vockley_REF_df.m$master_sample,
                                   levels=levels_master_samples,ordered=T)
  
  
  cat("Vockley_REF_df.m_2\n")
  cat(str(Vockley_REF_df.m))
  cat("\n")
  
  Vockley_REF_df.m<-merge(Vockley_REF_df.m,
                    Annot_DEF_subset,
                    by="master_sample")
  
  
  cat("Vockley_REF_df.m_3\n")
  cat(str(Vockley_REF_df.m))
  cat("\n")
  
  
  ### NO NA
  
  Vockley_REF_df.m_NO_NA<-Vockley_REF_df.m[!is.na(Vockley_REF_df.m$Vockley_REF),]
  
  
  cat("Vockley_REF_df.m_NO_NA_1\n")
  cat(str(Vockley_REF_df.m_NO_NA))
  cat("\n")
  
  ### Infinite reverted
  
  Vockley_REF_df.m_NO_NA$Vockley_REF[(is.infinite(Vockley_REF_df.m_NO_NA$Vockley_REF) &  Vockley_REF_df.m_NO_NA$Vockley_REF > 0)]<-60
  
  # Vockley_REF_df.m_NO_NA<-Vockley_REF_df.m_NO_NA[-(is.infinite(Vockley_REF_df.m_NO_NA$Vockley_REF) &  Vockley_REF_df.m_NO_NA$Vockley_REF < 0),]
  
  
  
  # Vockley_REF_df.m_NO_NA$Vockley_REF[(is.infinite(Vockley_REF_df.m_NO_NA$Vockley_REF) &  Vockley_REF_df.m_NO_NA$Vockley_REF < 0)]<--1*10000
  
  Vockley_REF_df.m_NO_NA$GROUP<-"NOR"
  
  cat("Vockley_REF_df.m_NO_NA_2_infinite_reverted\n")
  cat(str(Vockley_REF_df.m_NO_NA))
  cat("\n")
  
  
  
  
  
  A<-summary(Vockley_REF_df.m_NO_NA$Vockley_REF)
  
  cat("summary_Vockley_REF\n")
  cat(sprintf(as.character(names(A))))
  cat("\n")
  cat(sprintf(as.character(A)))
  cat("\n")
  
  
  check<-Vockley_REF_df.m_NO_NA[which(Vockley_REF_df.m_NO_NA$Vockley_REF >1.2 | Vockley_REF_df.m_NO_NA$Vockley_REF <0.8),]
  
  cat("check_Vockley_REF\n")
  cat(str(check))
  cat("\n")
  cat(sprintf(as.character(names(summary(check$master_sample)))))
  cat("\n")
  cat(sprintf(as.character(summary(check$master_sample))))
  cat("\n")
  
  A<-summary(check$Vockley_REF)
  
  cat("summary_Vockley_REF\n")
  cat(sprintf(as.character(names(A))))
  cat("\n")
  cat(sprintf(as.character(A)))
  cat("\n")
  
  
  check<-Vockley_REF_df.m_NO_NA[is.na(Vockley_REF_df.m_NO_NA$Vockley_REF),]
  
  cat("check_2\n")
  cat(str(check))
  cat("\n")
  
  
  List_values[['Vockley_REF']]<-Vockley_REF_df.m_NO_NA
  
  Vockley_REF_df.m_NO_NA.dt<-data.table(Vockley_REF_df.m_NO_NA,
                                  key=c("REAL_TILE_Plus_carried_variants","Cell_Type"))
  
  Vockley_REF_median<-unique(as.data.frame(Vockley_REF_df.m_NO_NA.dt[,.(median_Vockley_REF=median(Vockley_REF)),
                                                         by=key(Vockley_REF_df.m_NO_NA.dt)],stringsAsFactors=F))
  
  
  cat("Vockley_REF_median_0\n")
  cat(str(Vockley_REF_median))
  cat("\n")
  
  List_medians[['Vockley_REF']]<-Vockley_REF_median
  
  
  check_Vockley_REF<-Vockley_REF_median[which(Vockley_REF_median$median_Vockley_REF >1.2 | Vockley_REF_median$median_Vockley_REF <0.8),]
  
  cat("check_Vockley_REF_Vockley_REF\n")
  cat(str(check_Vockley_REF))
  cat("\n")
  cat(sprintf(as.character(names(summary(check_Vockley_REF$Cell_Type)))))
  cat("\n")
  cat(sprintf(as.character(summary(check_Vockley_REF$Cell_Type))))
  cat("\n")
  # cat(sprintf(as.character(names(summary(as.factor(check_Vockley_REF$REAL_TILE_Plus_carried_variants))))))
  # cat("\n")
  # cat(sprintf(as.character(summary(as.factor(check_Vockley_REF$REAL_TILE_Plus_carried_variants)))))
  # cat("\n")
  
  A<-summary(check_Vockley_REF$median_Vockley_REF)
  
  cat("summary_median_Vockley_REF\n")
  cat(sprintf(as.character(names(A))))
  cat("\n")
  cat(sprintf(as.character(A)))
  cat("\n")
 
  
  
  
  check_DEF<-check_LogFC[which(check_LogFC$REAL_TILE_Plus_carried_variants%in%check_Vockley_REF$REAL_TILE_Plus_carried_variants),]
  
  
  
  cat("check_DEF_Vockley_REF\n")
  cat(str(check_DEF))
  cat("\n")
  
  
  
  #### dna.REF.vs.dna.ALT ----
  
  
  
  dna.REF.vs.dna.ALT_matrix<-(dna_condition_df_NOR_REF_wide.matrix+0.0000001)/(dna_condition_df_NOR_ALT_wide.matrix+0.0000001)
  
  cat("--------------------------------------------------------------------------------------------->dna.REF.vs.dna.ALT_matrix\n")
  cat(str(dna.REF.vs.dna.ALT_matrix))
  cat("\n")
  
  dna.REF.vs.dna.ALT_df<-as.data.frame(dna.REF.vs.dna.ALT_matrix, stringsAsFactors=F)
  
  
  dna.REF.vs.dna.ALT_df$REAL_TILE_Plus_carried_variants<-row.names(dna.REF.vs.dna.ALT_df)
  
  cat("dna.REF.vs.dna.ALT_df\n")
  cat(str(dna.REF.vs.dna.ALT_df))
  cat("\n")
  
  dna.REF.vs.dna.ALT_df.m<-melt(dna.REF.vs.dna.ALT_df, id.cols="REAL_TILE_Plus_carried_variants")
  
  cat("dna.REF.vs.dna.ALT_df.m_1\n")
  cat(str(dna.REF.vs.dna.ALT_df.m))
  cat("\n")
  
  colnames(dna.REF.vs.dna.ALT_df.m)[which(colnames(dna.REF.vs.dna.ALT_df.m) == "variable")]<-"master_sample"
  colnames(dna.REF.vs.dna.ALT_df.m)[which(colnames(dna.REF.vs.dna.ALT_df.m) == "value")]<-"dna.REF.vs.dna.ALT"
  
  
  dna.REF.vs.dna.ALT_df.m$master_sample<-factor(dna.REF.vs.dna.ALT_df.m$master_sample,
                                         levels=levels_master_samples,ordered=T)
  
  
  cat("dna.REF.vs.dna.ALT_df.m_2\n")
  cat(str(dna.REF.vs.dna.ALT_df.m))
  cat("\n")
  
  dna.REF.vs.dna.ALT_df.m<-merge(dna.REF.vs.dna.ALT_df.m,
                          Annot_DEF_subset,
                          by="master_sample")
  
  
  cat("dna.REF.vs.dna.ALT_df.m_3\n")
  cat(str(dna.REF.vs.dna.ALT_df.m))
  cat("\n")
  
  
  ### NO NA
  
  dna.REF.vs.dna.ALT_df.m_NO_NA<-dna.REF.vs.dna.ALT_df.m[!is.na(dna.REF.vs.dna.ALT_df.m$dna.REF.vs.dna.ALT),]
  
  
  cat("dna.REF.vs.dna.ALT_df.m_NO_NA_1\n")
  cat(str(dna.REF.vs.dna.ALT_df.m_NO_NA))
  cat("\n")
  
  ### Infinite reverted
  
  dna.REF.vs.dna.ALT_df.m_NO_NA$dna.REF.vs.dna.ALT[(is.infinite(dna.REF.vs.dna.ALT_df.m_NO_NA$dna.REF.vs.dna.ALT) &  dna.REF.vs.dna.ALT_df.m_NO_NA$dna.REF.vs.dna.ALT > 0)]<-60
  #dna.REF.vs.dna.ALT_df.m_NO_NA$dna.REF.vs.dna.ALT[(is.infinite(dna.REF.vs.dna.ALT_df.m_NO_NA$dna.REF.vs.dna.ALT) &  dna.REF.vs.dna.ALT_df.m_NO_NA$dna.REF.vs.dna.ALT < 0)]<--1*10000
  
  # dna.REF.vs.dna.ALT_df.m_NO_NA<-dna.REF.vs.dna.ALT_df.m_NO_NA[-(is.infinite(dna.REF.vs.dna.ALT_df.m_NO_NA$dna.REF.vs.dna.ALT) &  dna.REF.vs.dna.ALT_df.m_NO_NA$dna.REF.vs.dna.ALT < 0),]
  
  
  
  
  dna.REF.vs.dna.ALT_df.m_NO_NA$GROUP<-"NOR"
  
  cat("dna.REF.vs.dna.ALT_df.m_NO_NA_2_infinite_reverted\n")
  cat(str(dna.REF.vs.dna.ALT_df.m_NO_NA))
  cat("\n")
  
  
  
  
  
  A<-summary(dna.REF.vs.dna.ALT_df.m_NO_NA$dna.REF.vs.dna.ALT)
  
  cat("summary_dna.REF.vs.dna.ALT\n")
  cat(sprintf(as.character(names(A))))
  cat("\n")
  cat(sprintf(as.character(A)))
  cat("\n")
  
  
  # check<-dna.REF.vs.dna.ALT_df.m_NO_NA[which(dna.REF.vs.dna.ALT_df.m_NO_NA$dna.REF.vs.dna.ALT >1.2 | dna.REF.vs.dna.ALT_df.m_NO_NA$dna.REF.vs.dna.ALT <0.8),]
  # 
  # cat("check_dna.REF.vs.dna.ALT\n")
  # cat(str(check))
  # cat("\n")
  # cat(sprintf(as.character(names(summary(check$master_sample)))))
  # cat("\n")
  # cat(sprintf(as.character(summary(check$master_sample))))
  # cat("\n")
  # 
  # A<-summary(check$dna.REF.vs.dna.ALT)
  # 
  # cat("summary_dna.REF.vs.dna.ALT\n")
  # cat(sprintf(as.character(names(A))))
  # cat("\n")
  # cat(sprintf(as.character(A)))
  # cat("\n")
  
  
  check<-dna.REF.vs.dna.ALT_df.m_NO_NA[is.na(dna.REF.vs.dna.ALT_df.m_NO_NA$dna.REF.vs.dna.ALT),]
  
  cat("check_2\n")
  cat(str(check))
  cat("\n")
  
  
  List_values[['dna.REF.vs.dna.ALT']]<-dna.REF.vs.dna.ALT_df.m_NO_NA
  
  dna.REF.vs.dna.ALT_df.m_NO_NA.dt<-data.table(dna.REF.vs.dna.ALT_df.m_NO_NA,
                                        key=c("REAL_TILE_Plus_carried_variants","Cell_Type"))
  
  dna.REF.vs.dna.ALT_median<-unique(as.data.frame(dna.REF.vs.dna.ALT_df.m_NO_NA.dt[,.(median_dna.REF.vs.dna.ALT=median(dna.REF.vs.dna.ALT)),
                                                                     by=key(dna.REF.vs.dna.ALT_df.m_NO_NA.dt)],stringsAsFactors=F))
  
  
  cat("dna.REF.vs.dna.ALT_median_0\n")
  cat(str(dna.REF.vs.dna.ALT_median))
  cat("\n")
  
  List_medians[['dna.REF.vs.dna.ALT']]<-dna.REF.vs.dna.ALT_median
  
  
  # check_dna.REF.vs.dna.ALT<-dna.REF.vs.dna.ALT_median[which(dna.REF.vs.dna.ALT_median$median_dna.REF.vs.dna.ALT >1.2 | dna.REF.vs.dna.ALT_median$median_dna.REF.vs.dna.ALT <0.8),]
  # 
  # cat("check_dna.REF.vs.dna.ALT_dna.REF.vs.dna.ALT\n")
  # cat(str(check_dna.REF.vs.dna.ALT))
  # cat("\n")
  # cat(sprintf(as.character(names(summary(check_dna.REF.vs.dna.ALT$Cell_Type)))))
  # cat("\n")
  # cat(sprintf(as.character(summary(check_dna.REF.vs.dna.ALT$Cell_Type))))
  # cat("\n")
  # # cat(sprintf(as.character(names(summary(as.factor(check_dna.REF.vs.dna.ALT$REAL_TILE_Plus_carried_variants))))))
  # # cat("\n")
  # # cat(sprintf(as.character(summary(as.factor(check_dna.REF.vs.dna.ALT$REAL_TILE_Plus_carried_variants)))))
  # # cat("\n")
  # 
  # A<-summary(check_dna.REF.vs.dna.ALT$median_dna.REF.vs.dna.ALT)
  # 
  # cat("summary_median_dna.REF.vs.dna.ALT\n")
  # cat(sprintf(as.character(names(A))))
  # cat("\n")
  # cat(sprintf(as.character(A)))
  # cat("\n")
  
  
  #### rna.REF.vs.rna.ALT ----
  
  
  
  rna.REF.vs.rna.ALT_matrix<-(rna_condition_df_NOR_REF_wide.matrix+0.0000001)/(rna_condition_df_NOR_ALT_wide.matrix+0.0000001)
  
  cat("--------------------------------------------------------------------------------------------->rna.REF.vs.rna.ALT_matrix\n")
  cat(str(rna.REF.vs.rna.ALT_matrix))
  cat("\n")
  
  rna.REF.vs.rna.ALT_df<-as.data.frame(rna.REF.vs.rna.ALT_matrix, stringsAsFactors=F)
  
  
  rna.REF.vs.rna.ALT_df$REAL_TILE_Plus_carried_variants<-row.names(rna.REF.vs.rna.ALT_df)
  
  cat("rna.REF.vs.rna.ALT_df\n")
  cat(str(rna.REF.vs.rna.ALT_df))
  cat("\n")
  
  rna.REF.vs.rna.ALT_df.m<-melt(rna.REF.vs.rna.ALT_df, id.cols="REAL_TILE_Plus_carried_variants")
  
  cat("rna.REF.vs.rna.ALT_df.m_1\n")
  cat(str(rna.REF.vs.rna.ALT_df.m))
  cat("\n")
  
  colnames(rna.REF.vs.rna.ALT_df.m)[which(colnames(rna.REF.vs.rna.ALT_df.m) == "variable")]<-"master_sample"
  colnames(rna.REF.vs.rna.ALT_df.m)[which(colnames(rna.REF.vs.rna.ALT_df.m) == "value")]<-"rna.REF.vs.rna.ALT"
  
  
  rna.REF.vs.rna.ALT_df.m$master_sample<-factor(rna.REF.vs.rna.ALT_df.m$master_sample,
                                                levels=levels_master_samples,ordered=T)
  
  
  cat("rna.REF.vs.rna.ALT_df.m_2\n")
  cat(str(rna.REF.vs.rna.ALT_df.m))
  cat("\n")
  
  rna.REF.vs.rna.ALT_df.m<-merge(rna.REF.vs.rna.ALT_df.m,
                                 Annot_DEF_subset,
                                 by="master_sample")
  
  
  cat("rna.REF.vs.rna.ALT_df.m_3\n")
  cat(str(rna.REF.vs.rna.ALT_df.m))
  cat("\n")
  
  
  ### NO NA
  
  rna.REF.vs.rna.ALT_df.m_NO_NA<-rna.REF.vs.rna.ALT_df.m[!is.na(rna.REF.vs.rna.ALT_df.m$rna.REF.vs.rna.ALT),]
  
  
  cat("rna.REF.vs.rna.ALT_df.m_NO_NA_1\n")
  cat(str(rna.REF.vs.rna.ALT_df.m_NO_NA))
  cat("\n")
  
  ### Infinite reverted
  
  rna.REF.vs.rna.ALT_df.m_NO_NA$rna.REF.vs.rna.ALT[(is.infinite(rna.REF.vs.rna.ALT_df.m_NO_NA$rna.REF.vs.rna.ALT) &  rna.REF.vs.rna.ALT_df.m_NO_NA$rna.REF.vs.rna.ALT > 0)]<-60
  # rna.REF.vs.rna.ALT_df.m_NO_NA<-rna.REF.vs.rna.ALT_df.m_NO_NA[-(is.infinite(rna.REF.vs.rna.ALT_df.m_NO_NA$rna.REF.vs.rna.ALT) &  rna.REF.vs.rna.ALT_df.m_NO_NA$rna.REF.vs.rna.ALT < 0),]
  
#  rna.REF.vs.rna.ALT_df.m_NO_NA$rna.REF.vs.rna.ALT[(is.infinite(rna.REF.vs.rna.ALT_df.m_NO_NA$rna.REF.vs.rna.ALT) &  rna.REF.vs.rna.ALT_df.m_NO_NA$rna.REF.vs.rna.ALT < 0)]<--1*10000
  
  rna.REF.vs.rna.ALT_df.m_NO_NA$GROUP<-"NOR"
  
  cat("rna.REF.vs.rna.ALT_df.m_NO_NA_2_infinite_reverted\n")
  cat(str(rna.REF.vs.rna.ALT_df.m_NO_NA))
  cat("\n")
  
  
  
  
  
  A<-summary(rna.REF.vs.rna.ALT_df.m_NO_NA$rna.REF.vs.rna.ALT)
  
  cat("summary_rna.REF.vs.rna.ALT\n")
  cat(sprintf(as.character(names(A))))
  cat("\n")
  cat(sprintf(as.character(A)))
  cat("\n")
  
  
  # check<-rna.REF.vs.rna.ALT_df.m_NO_NA[which(rna.REF.vs.rna.ALT_df.m_NO_NA$rna.REF.vs.rna.ALT >1.2 | rna.REF.vs.rna.ALT_df.m_NO_NA$rna.REF.vs.rna.ALT <0.8),]
  # 
  # cat("check_rna.REF.vs.rna.ALT\n")
  # cat(str(check))
  # cat("\n")
  # cat(sprintf(as.character(names(summary(check$master_sample)))))
  # cat("\n")
  # cat(sprintf(as.character(summary(check$master_sample))))
  # cat("\n")
  # 
  # A<-summary(check$rna.REF.vs.rna.ALT)
  # 
  # cat("summary_rna.REF.vs.rna.ALT\n")
  # cat(sprintf(as.character(names(A))))
  # cat("\n")
  # cat(sprintf(as.character(A)))
  # cat("\n")
  
  
  check<-rna.REF.vs.rna.ALT_df.m_NO_NA[is.na(rna.REF.vs.rna.ALT_df.m_NO_NA$rna.REF.vs.rna.ALT),]
  
  cat("check_2\n")
  cat(str(check))
  cat("\n")
  
  
  List_values[['rna.REF.vs.rna.ALT']]<-rna.REF.vs.rna.ALT_df.m_NO_NA
  
  rna.REF.vs.rna.ALT_df.m_NO_NA.dt<-data.table(rna.REF.vs.rna.ALT_df.m_NO_NA,
                                               key=c("REAL_TILE_Plus_carried_variants","Cell_Type"))
  
  rna.REF.vs.rna.ALT_median<-unique(as.data.frame(rna.REF.vs.rna.ALT_df.m_NO_NA.dt[,.(median_rna.REF.vs.rna.ALT=median(rna.REF.vs.rna.ALT)),
                                                                                   by=key(rna.REF.vs.rna.ALT_df.m_NO_NA.dt)],stringsAsFactors=F))
  
  
  cat("rna.REF.vs.rna.ALT_median_0\n")
  cat(str(rna.REF.vs.rna.ALT_median))
  cat("\n")
  
  List_medians[['rna.REF.vs.rna.ALT']]<-rna.REF.vs.rna.ALT_median
  
  
  # check_rna.REF.vs.rna.ALT<-rna.REF.vs.rna.ALT_median[which(rna.REF.vs.rna.ALT_median$median_rna.REF.vs.rna.ALT >1.2 | rna.REF.vs.rna.ALT_median$median_rna.REF.vs.rna.ALT <0.8),]
  # 
  # cat("check_rna.REF.vs.rna.ALT_rna.REF.vs.rna.ALT\n")
  # cat(str(check_rna.REF.vs.rna.ALT))
  # cat("\n")
  # cat(sprintf(as.character(names(summary(check_rna.REF.vs.rna.ALT$Cell_Type)))))
  # cat("\n")
  # cat(sprintf(as.character(summary(check_rna.REF.vs.rna.ALT$Cell_Type))))
  # cat("\n")
  # # cat(sprintf(as.character(names(summary(as.factor(check_rna.REF.vs.rna.ALT$REAL_TILE_Plus_carried_variants))))))
  # # cat("\n")
  # # cat(sprintf(as.character(summary(as.factor(check_rna.REF.vs.rna.ALT$REAL_TILE_Plus_carried_variants)))))
  # # cat("\n")
  # 
  # A<-summary(check_rna.REF.vs.rna.ALT$median_rna.REF.vs.rna.ALT)
  # 
  # cat("summary_median_rna.REF.vs.rna.ALT\n")
  # cat(sprintf(as.character(names(A))))
  # cat("\n")
  # cat(sprintf(as.character(A)))
  # cat("\n")
  
 
  
  
  #### Graph PRE and POST CELL TyPES aggregate_counts ----
  
  A<-summary(log10(dna_rna_pool_df$Median_value+ 0.00001))
  
  
  
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
  
  
  #### path2 ----
  
  path2<-paste(out,'NORMALIZATION_plots','/', sep='')
  
  cat("path2\n")
  cat(sprintf(as.character(path2)))
  cat("\n")
  
  if (file.exists(path2)){
    
    unlink(path2)
    dir.create(file.path(path2))
    
    
  } else {
    dir.create(file.path(path2))
    
  }
  
  
  
  
  for(i in 1:length(levels(dna_rna_pool_df$Cell_Type)))
  {
    cell_type_sel<-levels(dna_rna_pool_df$Cell_Type)[i]
    
    cat("cell_type_sel\n")
    cat(sprintf(as.character(cell_type_sel)))
    cat("\n")
    
    dna_rna_pool_df_subset<-dna_rna_pool_df[which(dna_rna_pool_df$Cell_Type == cell_type_sel),]
    
    cat("dna_rna_pool_df_subset_0\n")
    cat(str(dna_rna_pool_df_subset))
    cat("\n")
    
    
    ggridge_plot<-ggplot(data =dna_rna_pool_df_subset[order(dna_rna_pool_df_subset$interaction_1),],
                         aes(y=interaction_1,
                             x=log10(dna_rna_pool_df_subset$Median_value+ 0.00001),
                             fill=type,
                             alpha=GROUP)) +
      stat_density_ridges(quantile_lines = TRUE, quantiles = 2, geom="density_ridges",scale = 3)+
      theme_bw()+
      theme(axis.title.y=element_text(size=24, family="sans"),
            axis.title.x=element_text(size=24, family="sans"),
            axis.text.y=element_text(angle=0,size=10, color="black", family="sans"),
            axis.text.x=element_markdown(angle=45, vjust=1,hjust=1,size=18, color="black", family="sans"),
            legend.title=element_text(size=16,color="black", family="sans"),
            legend.text=element_text(size=12,color="black", family="sans"))+
      scale_y_discrete(name=NULL, drop=T)+
      scale_x_continuous(name="Counts (log10(value + 0.00001))",breaks=breaks.log_value,labels=labels_value, limits=c(breaks.log_value[1],breaks.log_value[length(breaks.log_value)]))+
      scale_fill_manual(values = c("red","green"),
                        drop=F)+
      scale_alpha_manual(values=c(0.5,2)) +
      theme(legend.position = "none")+
      ggeasy::easy_center_title()
  
    
    setwd(path2)
    
    svgname<-paste("ggridge_NORMALIZATION_NAMED_",cell_type_sel,".svg",sep='')
    makesvg = TRUE
    
    if (makesvg == TRUE)
    {
      ggsave(svgname, plot= ggridge_plot,
             device="svg",
             height=10, width=12)
    }
    
    cat("ggridge_NORMALIZATION_NAMED:\n")
    
    
    dna_rna_condition_df_subset_REF<-dna_rna_condition_df[which(dna_rna_condition_df$Cell_Type == cell_type_sel &
                                                              dna_rna_condition_df$condition == "REF"),]
    
    cat("dna_rna_condition_df_subset_REF_REF_0\n")
    cat(str(dna_rna_condition_df_subset_REF))
    cat("\n")
    
    
    ggridge_plot<-ggplot(data =dna_rna_condition_df_subset_REF[order(dna_rna_condition_df_subset_REF$interaction_1),],
                         aes(y=interaction_1,
                             x=log10(dna_rna_condition_df_subset_REF$Median_value+ 0.00001),
                             fill=type,
                             alpha=GROUP)) +
      stat_density_ridges(quantile_lines = TRUE, quantiles = 2, geom="density_ridges",scale = 3)+
      theme_bw()+
      theme(axis.title.y=element_text(size=24, family="sans"),
            axis.title.x=element_text(size=24, family="sans"),
            axis.text.y=element_text(angle=0,size=10, color="black", family="sans"),
            axis.text.x=element_markdown(angle=45, vjust=1,hjust=1,size=18, color="black", family="sans"),
            legend.title=element_text(size=16,color="black", family="sans"),
            legend.text=element_text(size=12,color="black", family="sans"))+
      scale_y_discrete(name=NULL, drop=T)+
      scale_x_continuous(name="Counts (log10(value + 0.00001))",breaks=breaks.log_value,labels=labels_value, limits=c(breaks.log_value[1],breaks.log_value[length(breaks.log_value)]))+
      scale_fill_manual(values = c("red","green"),
                        drop=F)+
      scale_alpha_manual(values=c(0.5,2)) +
      theme(legend.position = "none")+
      ggeasy::easy_center_title()
    
    
    setwd(path2)
    
    svgname<-paste("ggridge_NORMALIZATION_REF_NAMED_",cell_type_sel,".svg",sep='')
    makesvg = TRUE
    
    if (makesvg == TRUE)
    {
      ggsave(svgname, plot= ggridge_plot,
             device="svg",
             height=10, width=12)
    }
    
    cat("ggridge_NORMALIZATION_REF_NAMED:\n")
    
    
    dna_rna_condition_df_subset_ALT<-dna_rna_condition_df[which(dna_rna_condition_df$Cell_Type == cell_type_sel &
                                                                  dna_rna_condition_df$condition == "ALT"),]
    
    cat("dna_rna_condition_df_subset_ALT_ALT_0\n")
    cat(str(dna_rna_condition_df_subset_ALT))
    cat("\n")
    
    
    ggridge_plot<-ggplot(data =dna_rna_condition_df_subset_ALT[order(dna_rna_condition_df_subset_ALT$interaction_1),],
                         aes(y=interaction_1,
                             x=log10(dna_rna_condition_df_subset_ALT$Median_value+ 0.00001),
                             fill=type,
                             alpha=GROUP)) +
      stat_density_ridges(quantile_lines = TRUE, quantiles = 2, geom="density_ridges",scale = 3)+
      theme_bw()+
      theme(axis.title.y=element_text(size=24, family="sans"),
            axis.title.x=element_text(size=24, family="sans"),
            axis.text.y=element_text(angle=0,size=10, color="black", family="sans"),
            axis.text.x=element_markdown(angle=45, vjust=1,hjust=1,size=18, color="black", family="sans"),
            legend.title=element_text(size=16,color="black", family="sans"),
            legend.text=element_text(size=12,color="black", family="sans"))+
      scale_y_discrete(name=NULL, drop=T)+
      scale_x_continuous(name="Counts (log10(value + 0.00001))",breaks=breaks.log_value,labels=labels_value, limits=c(breaks.log_value[1],breaks.log_value[length(breaks.log_value)]))+
      scale_fill_manual(values = c("red","green"),
                        drop=F)+
      scale_alpha_manual(values=c(0.5,2)) +
      theme(legend.position = "none")+
      ggeasy::easy_center_title()
    
    
    setwd(path2)
    
    svgname<-paste("ggridge_NORMALIZATION_ALT_NAMED_",cell_type_sel,".svg",sep='')
    makesvg = TRUE
    
    if (makesvg == TRUE)
    {
      ggsave(svgname, plot= ggridge_plot,
             device="svg",
             height=10, width=12)
    }
    
    cat("ggridge_NORMALIZATION_ALT_NAMED:\n")
    
  }#i Cell_Type
  
 ##### SAVE LISTS with parameters ----
  
  setwd(out)
  filename_1<-paste("list_values",".rds", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  saveRDS(List_values, file=filename_1)
  
  
  filename_1<-paste("list_medians",".rds", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  saveRDS(List_medians, file=filename_1)
  

}

# ###############################################################################################################################################################################################################################################
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
    make_option(c("--LONG_MATRIX"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Initial_Selection"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--barcodes_per_tile"), type="character", default=NULL, 
                metavar="type1", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Threshold_EXP"), type="numeric", default=NULL, 
                metavar="Threshold_EXP", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--indir"), type="character", default=NULL, 
                metavar="indir", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--swapped_NEG.CTRLS"), type="character", default=NULL, 
                metavar="NEG.CTRLS", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--K562_replicates"), type="character", default=NULL, 
                metavar="NEG.CTRLS", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--CHRF_replicates"), type="character", default=NULL, 
                metavar="NEG.CTRLS", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--HL60_replicates"), type="character", default=NULL, 
                metavar="NEG.CTRLS", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--THP1_replicates"), type="character", default=NULL, 
                metavar="NEG.CTRLS", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="out", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "137_MPRA_normalization_and_filtering_Rscript_v2.R
                        --REAL_TILE_DF FILE.txt
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
  
  file_reader(opt)
  regularize_matrixes(opt)
  NOR_and_MPRAnalyze_files(opt)
  parameters_FC_Vockley_calculations(opt)
 
  
  
}


###########################################################################

system.time( main() )
