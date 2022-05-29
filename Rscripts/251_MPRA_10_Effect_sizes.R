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


opt = NULL

options(warn=1)

data_wrangling = function(option_list)
{
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  
  
  #### READ and transform out ----
  
  finemap_prob_Threshold = opt$finemap_prob_Threshold
  
  cat("finemap_prob_Threshold\n")
  cat(sprintf(as.character(finemap_prob_Threshold)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### Read TOME_correspondence ----
  
  TOME_correspondence = read.table(opt$TOME_correspondence, sep="\t", stringsAsFactors = F, header = T)
  
  cat("TOME_correspondence_\n")
  str(TOME_correspondence)
  cat("\n")
  
  TOME_correspondence$Lineage[TOME_correspondence$phenotype == "wbc"]<-"ALL_wbc_lineage"
  
  indx.int<-c(which(colnames(TOME_correspondence) == "phenotype"),
              which(colnames(TOME_correspondence) == "Lineage"))
  
  
  TOME_correspondence_subset<-unique(TOME_correspondence[,indx.int])
  
  cat("TOME_correspondence_subset_\n")
  str(TOME_correspondence_subset)
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(TOME_correspondence_subset$Lineage))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(TOME_correspondence_subset$Lineage)))))
  cat("\n")
  
  
  #### KEY_collpase_Plus_Variant_lineage_CLASSIFICATION ----
  
  
  KEY_collapse<-readRDS(file=opt$KEY_collpase_Plus_Variant_lineage_CLASSIFICATION)
  
  cat("KEY_collapse_0\n")
  cat(str(KEY_collapse))
  cat("\n")
  
  
  check<-KEY_collapse[is.na(KEY_collapse$CLASS_erythroid),]
  
  # cat("check_0\n")
  # cat(str(check))
  # cat("\n")
  # cat(sprintf(as.character(names(summary(check$VEP_DEF_LABELS)))))
  # cat("\n")
  # cat(sprintf(as.character(summary(check$VEP_DEF_LABELS))))
  # cat("\n")
  
  
  ##### exclude CTRLS -----
  
  KEY_collapse<-KEY_collapse[which(KEY_collapse$Label_2 == "ASSAYED_VARIANT"),]
  
  KEY_collapse<-droplevels(KEY_collapse)
  
  cat("KEY_collapse_excluded_CTRLS\n")
  cat(str(KEY_collapse))
  cat("\n")
  
  
  check<-KEY_collapse[is.na(KEY_collapse$CLASS_erythroid),]
  
  # cat("check_0\n")
  # cat(str(check))
  # cat("\n")
  
  
  ### Read dB file ----
  
  dB = as.data.frame(fread(file=opt$dB, sep="\t", stringsAsFactors = F, header = T), stringsAsFactors =F)
  
  dB$abs_finemap_z<-abs(dB$finemap_z)
  
  cat("dB_\n")
  cat(str(dB))
  cat("\n")
  

  dB_subset<-dB[which(dB$VAR%in%KEY_collapse$VAR),]
  
  # cat("dB_subset_\n")
  # cat(str(dB_subset))
  # cat("\n")
  
  dB_subset_thresholded<-dB_subset[which(dB_subset$finemap_prob >= finemap_prob_Threshold),]
  
  cat("dB_subset_thresholded_\n")
  cat(str(dB_subset_thresholded))
  cat("\n")
  
  ##### LOOP classify variants -----
  
  
  VARS<-unique(as.character(KEY_collapse$VAR))
  
  
  cat("VARS_\n")
  cat(str(VARS))
  cat("\n")
  
 
  temp<-data.frame(matrix(vector(), 0, 8,
                          dimnames=list(c(),
                                        c("phenotype","VAR","abs_finemap_z","LINEAGE",
                                          "KEY","enhancer_CLASS","E_Plus_ASE_CLASS","Cell_Type"))
  ))
  

  
  
  
  for(i in 1:length(VARS))
  {
    VAR_sel<-VARS[i]
    
    cat("--------------------->\t")
    cat(sprintf(as.character(VAR_sel)))
    cat("\t")
    
    carried_variants_sel<-gsub("chr","",VAR_sel)
    
    cat(sprintf(as.character(carried_variants_sel)))
    cat("\t")
  
    
    KEY_collapse_sel<-KEY_collapse[which(KEY_collapse$VAR == VAR_sel &
                                           KEY_collapse$carried_variants  == carried_variants_sel),]
    
    # cat("KEY_collapse_sel_\n")
    # cat(str(KEY_collapse_sel))
    # cat("\n")
    # 
    Element_sel<-unique(as.character(KEY_collapse_sel$KEY))
    
    
    cat(sprintf(as.character(Element_sel)))
    cat("\n")
    
    if(dim(KEY_collapse_sel)[1] > 0)
    {
    
        Cell_Types_array<-levels(KEY_collapse_sel$Cell_Type)
        
        # cat("Cell_Types_array_\n")
        # cat(str(Cell_Types_array))
        # cat("\n")
        
        for(l in 1:length(Cell_Types_array))
        {
          Cell_Types_array_sel<-Cell_Types_array[l]
          
          
          # cat("--------------------->\t")
          # cat(sprintf(as.character(Cell_Types_array_sel)))
          # cat("\t")
          
          
          KEY_collapse_CT<-KEY_collapse_sel[which(KEY_collapse_sel$Cell_Type == Cell_Types_array_sel),]
          
          # cat("KEY_collapse_CT_\n")
          # cat(str(KEY_collapse_CT))
          # cat("\n")
          
          if(dim(KEY_collapse_CT)[1] > 0)
          {
            enhancer_CLASS_string<-paste(KEY_collapse_CT$enhancer_CLASS,collapse=";")
            
            # cat("enhancer_CLASS_string>\t")
            # cat(sprintf(as.character(enhancer_CLASS_string)))
            # cat("\n")
            
            Element_enhancer_CLASS<-NULL
            
            indx_0<-grep("0", enhancer_CLASS_string)
            
            # cat("indx_0_\n")
            # cat(str(indx_0))
            # cat("\n")
            
            if(length(indx_0) > 0)
            {
              
              Element_enhancer_CLASS<-"0"
            }
            
            indx_1<-grep("1", enhancer_CLASS_string)
            
            # cat("indx_1_\n")
            # cat(str(indx_1))
            # cat("\n")
            
            if(length(indx_1) > 0)
            {
              
              Element_enhancer_CLASS<-"1"
            }
            
            indx_2<-grep("2", enhancer_CLASS_string)
            
            # cat("indx_2_\n")
            # cat(str(indx_2))
            # cat("\n")
            
            if(length(indx_2) > 0)
            {
              
              Element_enhancer_CLASS<-"2"
            }
            
            indx_3<-grep("3 or >", enhancer_CLASS_string)
            
            # cat("indx_3_\n")
            # cat(str(indx_3))
            # cat("\n")
            
            if(length(indx_3) > 0)
            {
              
              Element_enhancer_CLASS<-"3 or >"
            }
            
            E_Plus_ASE_CLASS_string<-paste(KEY_collapse_CT$E_Plus_ASE_CLASS,collapse=";")
            
            # cat("E_Plus_ASE_CLASS_string>\t")
            # cat(sprintf(as.character(E_Plus_ASE_CLASS_string)))
            # cat("\n")
            
            Element_E_Plus_ASE_CLASS<-NULL
            
            indx_0<-grep("0", E_Plus_ASE_CLASS_string)
            
            # cat("indx_0_\n")
            # cat(str(indx_0))
            # cat("\n")
            
            if(length(indx_0) > 0)
            {
              
              Element_E_Plus_ASE_CLASS<-"0"
            }
            
            indx_1<-grep("1", E_Plus_ASE_CLASS_string)
            
            # cat("indx_1_\n")
            # cat(str(indx_1))
            # cat("\n")
            
            if(length(indx_1) > 0)
            {
              
              Element_E_Plus_ASE_CLASS<-"1"
            }
            
            indx_2<-grep("2", E_Plus_ASE_CLASS_string)
            
            # cat("indx_2_\n")
            # cat(str(indx_2))
            # cat("\n")
            
            if(length(indx_2) > 0)
            {
              
              Element_E_Plus_ASE_CLASS<-"2"
            }
            
            indx_3<-grep("3 or >", E_Plus_ASE_CLASS_string)
            
            # cat("indx_3_\n")
            # cat(str(indx_3))
            # cat("\n")
            
            if(length(indx_3) > 0)
            {
              
              Element_E_Plus_ASE_CLASS<-"3 or >"
            }
            # cat(sprintf(as.character(Element_sel)))
            # cat("\t")
            # cat(sprintf(as.character(Element_enhancer_CLASS)))
            # cat("\t")
            # cat(sprintf(as.character(Element_E_Plus_ASE_CLASS)))
            # cat("\n")
            
            
            #### phenotypes and absolute finemap_z ----
            
          
            dB_subset_thresholded_sel<-dB_subset_thresholded[which(dB_subset_thresholded$VAR == VAR_sel),]
            
            # cat("dB_subset_thresholded_sel_\n")
            # cat(str(dB_subset_thresholded_sel))
            # cat("\n")
            
            indx_phenotypes<-c(which(colnames(dB_subset_thresholded_sel) == "VAR"),
                                      which(colnames(dB_subset_thresholded_sel) == "phenotype"),which(colnames(dB_subset_thresholded_sel) == "abs_finemap_z"))
            
            
            subset_phenotypes<-unique(dB_subset_thresholded_sel[,indx_phenotypes])
            
            
            # cat("subset_phenotypes_0\n")
            # cat(str(subset_phenotypes))
            # cat("\n")
            
            subset_phenotypes<-merge(subset_phenotypes,
                                     TOME_correspondence_subset,
                                     by="phenotype")
            
            
            # cat("subset_phenotypes_1\n")
            # cat(str(subset_phenotypes))
            # cat("\n")
            
            subset_phenotypes$KEY<-Element_sel
            subset_phenotypes$enhancer_CLASS<-Element_enhancer_CLASS
            subset_phenotypes$E_Plus_ASE_CLASS<-Element_E_Plus_ASE_CLASS
            subset_phenotypes$Cell_Type<-Cell_Types_array_sel
            
            # cat("subset_phenotypes_2\n")
            # cat(str(subset_phenotypes))
            # cat("\n")
            
            temp<-rbind(temp,subset_phenotypes)
            
            
          }# dim(KEY_collapse_CT)[1] > 0
        }# l Cell_Types_array
    }# dim(KEY_collapse_sel)[1] > 0
  }#i
  
  cat("temp_FINAL\n")
  cat(str(temp))
  cat("\n")
  
  temp$enhancer_CLASS<-factor(temp$enhancer_CLASS,
                                                    levels = c("0","1","2","3 or >"),
                                                    ordered=T)
  
  
  temp$E_Plus_ASE_CLASS<-factor(temp$E_Plus_ASE_CLASS,
                              levels = c("0","1","2","3 or >"),
                              ordered=T)
  
  temp$Cell_Type<-factor(temp$Cell_Type,
                                levels = c("K562","CHRF","HL60","THP1"),
                                ordered=T)
  
  temp$Lineage<-factor(temp$Lineage,
                         levels = c("erythroid_lineage","mega_lineage","gran_mono_lineage","lymph_lineage","ALL_wbc_lineage"),
                         ordered=T)
  
  temp$phenotype<-factor(as.character(temp$phenotype),
                                 levels=c("plt","mpv","pdw","pct","H-IPF","IPF_perc","IPF","P-LCR",
                                          "rbc","mcv","hct","mch","mchc","hgb","rdw_cv","MacrorR","MicroR","RBC-He","Delta-He","RET-RBC-FSC","Hyer-He","RDW-SD","RPI","HFR","LFR","MFR","IG","IG_perc",
                                          "ret","ret_p","irf","hlr","hlr_p","mrv","mscv","RET-FSC","IRF-FSC","RET-He","RET-UPP",
                                          "lymph","lymph_p","LY-FSC","LY-FSC-DW","LY-SSC","LY-SSC-DW","LY-SFL","LY-SFL-DW",
                                          "mono","neut","eo","baso","wbc","mono_p","neut_p","eo_p","baso_p",
                                          "MO-FSC","MO-FSC-DW","MO-SSC","MO-SSC-DW","MO-SFL","MO-SFL-DW",
                                          "NE-FSC","NE-FSC-DW","NE-SSC","NE-SSC-DW","NE-SFL","NE-SFL-DW"
                                 )
                                 , ordered =T)
  
  #### Nicole's new phenotype labels ----
  
  temp$phenotype_DEF<-revalue(temp$phenotype, 
                                      c("plt"="PLT#",
                                        "mpv"="MPV",
                                        "pdw"="PDW",
                                        "pct"="PCT",
                                        "H-IPF"="H-IPF",
                                        "IPF_perc"="IPF%",
                                        "IPF"="IPF#",
                                        "P-LCR"="P-LCR",
                                        "rbc"="RBC#",
                                        "mcv"="MCV",
                                        "hct"="HCT",
                                        "mch"="MCH",
                                        "mchc"="MCHC",
                                        "hgb"="HGB",
                                        "rdw_cv"='RDW',
                                        "MacrorR"='MacroR',
                                        "MicroR"='MicroR',
                                        "RBC-He"='RBC-He',
                                        "Delta-He"='Delta-He',
                                        "Hyer-He"='Hyper-He',
                                        "RDW-SD"='RDW-SD',
                                        "RPI"='RPI',
                                        "HFR"='HFR',
                                        "MFR"='MFR',
                                        "LFR"='LFR',
                                        "IG"='IG#',
                                        "IG_perc"='IG%',
                                        "ret"='RET#',
                                        "ret_p"='RET%',
                                        "irf"='IRF',
                                        "hlr"='HLSR#',
                                        "hlr_p"='HLSR%',
                                        "mrv"='MRV',
                                        "mscv"='MSCV',
                                        "RET-FSC"='RET-FSC',
                                        "RET-RBC-FSC"='RET-FSC',
                                        "IRF-FSC"='IRF-FSC',
                                        "RET-He"='RET-He',
                                        "RET-UPP"='RET-UPP',
                                        "lymph"='LYMPH#',
                                        "lymph_p"='LYMPH%',
                                        "LY-FSC"="LY-FSC",
                                        "LY-FSC-DW"="LY-FSC-DW",
                                        "LY-SSC"="LY-SSC",
                                        "LY-SSC-DW"="LY-SSC-DW",
                                        "LY-SFL"="LY-SFL",
                                        "LY-SFL-DW"="LY-SFL-DW",
                                        "baso"='BASO#',
                                        "baso_p"='BASO%',
                                        "eo"='EO#',
                                        "eo_p"='EO%',
                                        "mono"='MONO#',
                                        "mono_p"='MONO%',
                                        "MO-FSC"="MO-FSC",
                                        "MO-FSC-DW"="MO-FSC-DW",
                                        "MO-SSC"="MO-SSC",
                                        "MO-SSC-DW"="MO-SSC-DW",
                                        "MO-SFL"="MO-SFL",
                                        "MO-SFL-DW"="MO-SFL-DW",
                                        "neut"='NEUT#',
                                        "neut_p"='NEUT%',
                                        "NE-FSC"="NE-FSC",
                                        "NE-FSC-DW"="NE-FSC-DW",
                                        "NE-SSC"="NE-SSC",
                                        "NE-SSC-DW"="NE-SSC-DW",
                                        "NE-SFL"="NE-SFL",
                                        "NE-SFL-DW"="NE-SFL-DW",
                                        "wbc"='WBC#'))
  
  
  temp$phenotype_DEF<-factor(as.character(temp$phenotype_DEF),
                                     levels=c("PLT#","MPV","PDW","PCT","H-IPF","IPF%","IPF#","P-LCR",
                                              "RBC#","MCV","HCT","MCH","MCHC","HGB",'RDW','MacroR','MicroR','RBC-He','Delta-He','Hyper-He','RDW-SD','RPI','HFR','MFR','LFR','IG#','IG%',
                                              'RET#','RET%','IRF','HLSR#','HLSR%','MRV','MSCV',"RET-FSC","RET-RBC-FSC","IRF-FSC","RET-He","RET-UPP",
                                              'LYMPH#','LYMPH%',"LY-FSC","LY-FSC-DW","LY-SSC","LY-SSC-DW","LY-SFL","LY-SFL-DW",
                                              'BASO#','BASO%','EO#','EO%',
                                              'MONO#','MONO%',"MO-FSC","MO-FSC-DW","MO-SSC","MO-SSC-DW","MO-SFL","MO-SFL-DW",
                                              'NEUT#','NEUT%',"NE-FSC","NE-FSC-DW","NE-SSC","NE-SSC-DW","NE-SFL","NE-SFL-DW",
                                              'WBC#'), ordered = T)
  
        
  temp<-droplevels(temp)
  
  cat("temp_FINAL_2\n")
  cat(str(temp))
  cat("\n")
  cat(sprintf(as.character(names(summary(temp$enhancer_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(temp$enhancer_CLASS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(temp$E_Plus_ASE_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(temp$E_Plus_ASE_CLASS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(temp$Cell_Type)))))
  cat("\n")
  cat(sprintf(as.character(summary(temp$Cell_Type))))
  cat("\n")
  cat(sprintf(as.character(names(summary(temp$Lineage)))))
  cat("\n")
  cat(sprintf(as.character(summary(temp$Lineage))))
  cat("\n")
  cat(sprintf(as.character(names(summary(temp$phenotype)))))
  cat("\n")
  cat(sprintf(as.character(summary(temp$phenotype))))
  cat("\n")
  
 

  
  ##### LINEAGE_plots  ----
  
  
  path5<-paste(out,'GWAS_plots','/', sep='')
  
  cat("path5\n")
  cat(sprintf(as.character(path5)))
  cat("\n")
  
  
  if (file.exists(path5)){
    
    
    
    
  } else {
    dir.create(file.path(path5))
    
  }
  
  
  ####    SAVE ----
  
  setwd(path5)
  
  # quit(status = 1)
  

  filename_1<-paste("Element_collapse_Plus_Variant_Phenotypes",".rds", sep='')
  saveRDS(file=filename_1,temp)
  
  
}

data_wrangling_FC_Vockley_REF = function(option_list)
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
 
  ### Read RAW data ----
  
  
  
  path5<-paste(out,'GWAS_plots','/', sep='')
  
  cat("path5\n")
  cat(sprintf(as.character(path5)))
  cat("\n")
  
  
  if (file.exists(path5)){
    
    
    
    
  } else {
    dir.create(file.path(path5))
    
  }
  
  
  ####    SAVE ----
  
  setwd(path5)
  
  # quit(status = 1)
  
  
  filename_1<-paste("Element_collapse_Plus_Variant_Phenotypes",".rds", sep='')
  
  
  Element_collapse<-readRDS(file=filename_1)
  
  cat("Element_collapse\n")
  cat(str(Element_collapse))
  cat("\n")
  cat(sprintf(as.character(names(summary(Element_collapse$enhancer_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Element_collapse$enhancer_CLASS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Element_collapse$E_Plus_ASE_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Element_collapse$E_Plus_ASE_CLASS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Element_collapse$Cell_Type)))))
  cat("\n")
  cat(sprintf(as.character(summary(Element_collapse$Cell_Type))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Element_collapse$Lineage)))))
  cat("\n")
  cat(sprintf(as.character(summary(Element_collapse$Lineage))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Element_collapse$phenotype)))))
  cat("\n")
  cat(sprintf(as.character(summary(Element_collapse$phenotype))))
  cat("\n")
  
  #### LOOP ----
  
  VARS<-unique(as.character(Element_collapse$VAR))
  
  
  cat("VARS_\n")
  cat(str(VARS))
  cat("\n")
  
  Cell_Types_array<-levels(Element_collapse$Cell_Type)
  
  cat("Cell_Types_array_\n")
  cat(str(Cell_Types_array))
  cat("\n")
  
  
  
  temp<-data.frame(matrix(vector(), 0, 7,
                          dimnames=list(c(),
                                        c("VAR","carried_variants","KEY","LogFC",
                                          "Vockley_REF","abs_Vockley_REF","Cell_Type"))
  ))
  
  
  
  
  for(i in 1:length(VARS))
  {
    VAR_sel<-VARS[i]
    
    cat("--------------------->\t")
    cat(sprintf(as.character(VAR_sel)))
    cat("\n")
    
    carried_variants_sel<-gsub("chr","",VAR_sel)
    
    # cat(sprintf(as.character(carried_variants_sel)))
    # cat("\t")
    
    
    MPRA_Real_tile_sel<-MPRA_Real_tile[which(MPRA_Real_tile$VAR == VAR_sel &
                                           MPRA_Real_tile$carried_variants  == carried_variants_sel),]
    
    # cat("MPRA_Real_tile_sel_\n")
    # cat(str(MPRA_Real_tile_sel))
    # cat("\n")
    
       
    
    if(dim(MPRA_Real_tile_sel)[1] > 0)
    {
           
      for(l in 1:length(Cell_Types_array))
      {
        Cell_Types_array_sel<-Cell_Types_array[l]
        
        
        # cat("--------------------->\t")
        # cat(sprintf(as.character(Cell_Types_array_sel)))
        # cat("\t")
        
        MPRA_Real_tile_sel_CT_sel<-MPRA_Real_tile_sel[which(MPRA_Real_tile_sel$Cell_Type == Cell_Types_array_sel),]
        
        # cat("MPRA_Real_tile_sel_CT_sel_\n")
        # cat(str(MPRA_Real_tile_sel_CT_sel))
        # cat("\n")
        
       
        
        if(dim(MPRA_Real_tile_sel_CT_sel)[1] > 0)
        {
          MPRA_Real_tile_LogFC_Vockley_REF<-MPRA_Real_tile_sel_CT_sel[c(which(MPRA_Real_tile_sel_CT_sel$variable == "LogFC"),
                                                                        which(MPRA_Real_tile_sel_CT_sel$variable == "Vockley_REF")),]
          
          # cat("MPRA_Real_tile_LogFC_Vockley_REF_\n")
          # cat(str(MPRA_Real_tile_LogFC_Vockley_REF))
          # cat("\n")
          
         
          MPRA_Real_tile_LogFC_Vockley_REF_wide<-as.data.frame(pivot_wider(MPRA_Real_tile_LogFC_Vockley_REF,
                                                                  id_cols=c("KEY","TILE","VAR","carried_variants"),
                                                                  names_from=variable,
                                                                  values_from=value),stringsAsFactors=F)
          
          MPRA_Real_tile_LogFC_Vockley_REF_wide$abs_Vockley_REF<-abs(MPRA_Real_tile_LogFC_Vockley_REF_wide$Vockley_REF)
          
          
          # cat("MPRA_Real_tile_LogFC_Vockley_REF_wide_0\n")
          # cat(str(MPRA_Real_tile_LogFC_Vockley_REF_wide))
          # cat("\n")
          
          MPRA_Real_tile_LogFC_max<-setDT(MPRA_Real_tile_LogFC_Vockley_REF_wide)[, .SD[which.max(LogFC)], by=.(KEY)]
          
          # cat("MPRA_Real_tile_LogFC_max_0\n")
          # cat(str(MPRA_Real_tile_LogFC_max))
          # cat("\n")
          
          
          
          MPRA_Real_tile_abs_Vockley_REF_max<-setDT(MPRA_Real_tile_LogFC_Vockley_REF_wide)[, .SD[which.max(abs_Vockley_REF)], by=.(KEY)]
          
          # cat("MPRA_Real_tile_abs_Vockley_REF_max_0\n")
          # cat(str(MPRA_Real_tile_abs_Vockley_REF_max))
          # cat("\n")
          
          
          A.df<-as.data.frame(cbind(VAR_sel,carried_variants_sel,MPRA_Real_tile_abs_Vockley_REF_max$KEY,MPRA_Real_tile_LogFC_max$LogFC,
                              MPRA_Real_tile_abs_Vockley_REF_max$Vockley_REF,MPRA_Real_tile_abs_Vockley_REF_max$abs_Vockley_REF), stringsAsFactors=F)
          
          colnames(A.df)<-c("VAR","carried_variants","KEY","LogFC","Vockley_REF","abs_Vockley_REF")
          
          # cat("A.df_0\n")
          # cat(str(A.df))
          # cat("\n")
        }# dim(MPRA_Real_tile_sel_CT_sel)[1] > 0
        
        A.df$Cell_Type<-Cell_Types_array_sel
        
        temp<-rbind(temp,
                    A.df)
        # cat("temp_0\n")
        # cat(str(temp))
        # cat("\n")
      }#l
    }#dim(MPRA_Real_tile_sel)[1] > 0
  }#i
  
  
  temp$LogFC<-as.numeric(temp$LogFC)
  temp$Vockley_REF<-as.numeric(temp$Vockley_REF)
  temp$abs_Vockley_REF<-as.numeric(temp$abs_Vockley_REF)
  
  temp$Cell_Type<-factor(temp$Cell_Type,
                         levels = c("K562","CHRF","HL60","THP1"),
                         ordered=T)
  
  cat("temp_FINAL\n")
  cat(str(temp))
  cat("\n")
  
  
  cat("Element_collapse_0\n")
  cat(str(Element_collapse))
  cat("\n")
  
  
  Element_collapse<-merge(Element_collapse,
                          temp,
                          by=c("VAR","KEY","Cell_Type"))
  
  
  
  cat("Element_collapse_1\n")
  cat(str(Element_collapse))
  cat("\n")
  
  ##### LINEAGE_plots  ----
  
  
  path5<-paste(out,'GWAS_plots','/', sep='')
  
  cat("path5\n")
  cat(sprintf(as.character(path5)))
  cat("\n")
  
  
  if (file.exists(path5)){
    
    
    
    
  } else {
    dir.create(file.path(path5))
    
  }
  
  
  ####    SAVE ----
  
  setwd(path5)
  
  # quit(status = 1)
  
  
  filename_1<-paste("Element_collapse_Plus_Variant_Phenotypes_Plus_LogFC_Vockley_REF",".rds", sep='')
  saveRDS(file=filename_1,Element_collapse)
}

ggridge_enhancer_CLASS = function(option_list)
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
  
  
  
  path5<-paste(out,'GWAS_plots','/', sep='')
  
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
  
  
  filename_1<-paste("Element_collapse_Plus_Variant_Phenotypes_Plus_LogFC_Vockley_REF",".rds", sep='')
  
  
  Element_collapse<-readRDS(file=filename_1)
  
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
  
  
  
  #### LOOP ggridge ----
  
  
  Cell_Types_array<-levels(Element_collapse$Cell_Type)
  
  cat("Cell_Types_array_2\n")
  cat(str(Cell_Types_array))
  cat("\n")
  
  
  list_STATS_DEF<-list()
  list_LM_FC_DEF<-list()
  list_LM_abs_Vockley_REF_DEF<-list()
  
  
  for(i in 1:length(Cell_Types_array))
  {
    
    Cell_Types_array_sel<-Cell_Types_array[i]
    
    cat("--->\t")
    cat(sprintf(as.character(Cell_Types_array_sel)))
    cat("\n")
    
    
    Element_collapse_CT<-Element_collapse[which(Element_collapse$Cell_Type == Cell_Types_array_sel),]
    
    Element_collapse_CT<-droplevels(Element_collapse_CT)
    
    # cat("Element_collapse_CT\n")
    # cat(str(Element_collapse_CT))
    # cat("\n")
    
    
    AT_LEAST_CLASSES<-c("AT_LEAST_CLASS_1","AT_LEAST_CLASS_2","AT_LEAST_CLASS_3")
    
    list_STATS_CLASSES<-list()
    list_LM_FC_CLASSES<-list()
    list_LM_abs_Vockley_REF_CLASSES<-list()
    
    
    for(k in 1:length(AT_LEAST_CLASSES))
    {
      
      AT_LEAST_CLASSES_sel<-AT_LEAST_CLASSES[k]
      
      # cat("-CLASS COMPARISON-->\t")
      # cat(sprintf(as.character(AT_LEAST_CLASSES_sel)))
      # cat("\t")
      
      indx.CLASS<-which(colnames(Element_collapse_CT) == AT_LEAST_CLASSES_sel)
      
      # cat(sprintf(as.character(indx.CLASS)))
      # cat("\n")
      
      Element_collapse_CT_CLASS_sel<-Element_collapse_CT[!is.na(Element_collapse_CT[,indx.CLASS]),]
      
      cat("Element_collapse_CT_CLASS_sel\n")
      cat(str(Element_collapse_CT_CLASS_sel))
      cat("\n")
      
      if(dim(Element_collapse_CT_CLASS_sel)[1] >0)
      {
        
        Lineage_array<-levels(Element_collapse_CT$Lineage)
        
        cat("Lineage_array\n")
        cat(str(Lineage_array))
        cat("\n")
        
        list_STATS_Lineage<-list()
        list_LM_FC_Lineage<-list()
        list_LM_abs_Vockley_REF_Lineage<-list()
        
        
        for(l in 1:length(Lineage_array))
        {
          Lineage_array_sel<-Lineage_array[l]
          
          cat("--->\t")
          cat(sprintf(as.character(Lineage_array_sel)))
          cat("\n")
          
          Element_collapse_CT_CLASS_Lineage_sel<-Element_collapse_CT_CLASS_sel[which(Element_collapse_CT_CLASS_sel$Lineage == Lineage_array_sel),]
          
          Element_collapse_CT_CLASS_Lineage_sel<-droplevels(Element_collapse_CT_CLASS_Lineage_sel)
          
          cat("Element_collapse_CT_CLASS_Lineage_sel\n")
          cat(str(Element_collapse_CT_CLASS_Lineage_sel))
          cat("\n")
          
          
          phenotype_array<-levels(Element_collapse_CT_CLASS_Lineage_sel$phenotype)
          
          cat("phenotype_array\n")
          cat(str(phenotype_array))
          cat("\n")
          
          list_STATS_phenotype<-list()
          list_LM_FC_phenotype<-list()
          list_LM_abs_Vockley_REF_phenotype<-list()
          
          
          for(m in 1:length(phenotype_array))
          {
            phenotype_array_sel<-phenotype_array[m]
            
            cat("--->\t")
            cat(sprintf(as.character(phenotype_array_sel)))
            cat("\n")
            
            Element_collapse_CT_CLASS_Lineage_phenotype_sel<-Element_collapse_CT_CLASS_Lineage_sel[which(Element_collapse_CT_CLASS_Lineage_sel$phenotype == phenotype_array_sel),]
            
            
            
            if(dim(Element_collapse_CT_CLASS_Lineage_phenotype_sel)[1] >0)
            {
              cat("Element_collapse_CT_CLASS_Lineage_phenotype_sel\n")
              cat(str(Element_collapse_CT_CLASS_Lineage_phenotype_sel))
              cat("\n")
              
              indx.CLASS_2<-which(colnames(Element_collapse_CT_CLASS_Lineage_phenotype_sel) == AT_LEAST_CLASSES_sel)
              
              
              Element_collapse_CT_CLASS_Lineage_phenotype_sel_ACTIVE<-Element_collapse_CT_CLASS_Lineage_phenotype_sel[which(as.numeric(Element_collapse_CT_CLASS_Lineage_phenotype_sel[,indx.CLASS_2]) == 2),]
              
              cat("Element_collapse_CT_CLASS_Lineage_phenotype_sel_ACTIVE\n")
              cat(str(Element_collapse_CT_CLASS_Lineage_phenotype_sel_ACTIVE))
              cat("\n")
              
              
              
              if(dim(Element_collapse_CT_CLASS_Lineage_phenotype_sel_ACTIVE)[1] >= 2)
              {
                
                Element_collapse_CT_CLASS_Lineage_phenotype_sel_ACTIVE$FC<-logratio2foldchange(Element_collapse_CT_CLASS_Lineage_phenotype_sel_ACTIVE$LogFC, base=2)
                
                
                # cat("Element_collapse_CT_CLASS_Lineage_phenotype_sel_ACTIVE_2\n")
                # cat(str(Element_collapse_CT_CLASS_Lineage_phenotype_sel_ACTIVE))
                # cat("\n")
                # quit(status = 1)
                
                ####models of relation effect size GWAS and effect size MPRA ----
                
                linearModel_abs_Vockley_REF <- summary(lm(abs_finemap_z ~ abs_Vockley_REF, 
                                                          data=Element_collapse_CT_CLASS_Lineage_phenotype_sel_ACTIVE))
                
                # cat("linearModel_abs_Vockley_REF\n")
                # cat(str(linearModel_abs_Vockley_REF))
                # cat("\n")
                # 
                
                r.squared_linearModel_abs_Vockley_REF<-round(linearModel_abs_Vockley_REF$r.squared,3)
                adj.r.squared_linearModel_abs_Vockley_REF<-round(linearModel_abs_Vockley_REF$adj.r.squared,3)
                
                
                
                linearModel_abs_Vockley_REF_coeffcient_df.m<-melt(linearModel_abs_Vockley_REF$coefficients)
                
                colnames(linearModel_abs_Vockley_REF_coeffcient_df.m)[which(colnames(linearModel_abs_Vockley_REF_coeffcient_df.m)=="Var1")]<-"Terms"
                colnames(linearModel_abs_Vockley_REF_coeffcient_df.m)[which(colnames(linearModel_abs_Vockley_REF_coeffcient_df.m)=="Var2")]<-"Parameters"
                
                linearModel_abs_Vockley_REF_coeffcient_df.m$Terms<-as.character(linearModel_abs_Vockley_REF_coeffcient_df.m$Terms)
                linearModel_abs_Vockley_REF_coeffcient_df.m$Parameters<-as.character(linearModel_abs_Vockley_REF_coeffcient_df.m$Parameters)
                linearModel_abs_Vockley_REF_coeffcient_df.m$Cell_Type<-Cell_Types_array_sel
                #  linearModel_abs_Vockley_REF_coeffcient_df.m$LINEAGE_CLASSIF_DEF<-LINEAGE_CLASSIF_DEF_array_sel
                linearModel_abs_Vockley_REF_coeffcient_df.m$adj.r.squared_linearModel_abs_Vockley_REF<-adj.r.squared_linearModel_abs_Vockley_REF
                
                
                
                # cat("linearModel_abs_Vockley_REF_coeffcient_df.m\n")
                # cat(str(linearModel_abs_Vockley_REF_coeffcient_df.m))
                # cat("\n")
                # 
                
                FC_slope<-round(linearModel_abs_Vockley_REF_coeffcient_df.m$value[which(linearModel_abs_Vockley_REF_coeffcient_df.m$Parameters== "Estimate" &
                                                                                          linearModel_abs_Vockley_REF_coeffcient_df.m$Terms == "abs_Vockley_REF")],2)
                
                
                # cat("FC_slope\n")
                # cat(str(FC_slope))
                # cat("\n")
                
                
                FC_intercept<-round(linearModel_abs_Vockley_REF_coeffcient_df.m$value[which(linearModel_abs_Vockley_REF_coeffcient_df.m$Parameters== "Estimate" &
                                                                                              linearModel_abs_Vockley_REF_coeffcient_df.m$Terms == "(Intercept)")],2)
                
                
                # cat("FC_intercept\n")
                # cat(str(FC_intercept))
                # cat("\n")
                
                A.df<-as.data.frame(cbind(phenotype_array_sel,"abs_Vockley_REF",FC_slope,FC_intercept,dim(Element_collapse_CT_CLASS_Lineage_phenotype_sel_ACTIVE)[1],r.squared_linearModel_abs_Vockley_REF))
                
                # cat("A.df\n")
                # cat(str(A.df))
                # cat("\n")
                
                colnames(A.df)<-c("phenotype","MPRA_parameter","LM_slope","LM_intercept","n","r.squared")
                
                # cat("A.df\n")
                # cat(str(A.df))
                # cat("\n")
                
                list_LM_abs_Vockley_REF_phenotype[[m]]<-A.df
                
                
                
                # quit(status = 1)
                
                
                ####models of relation effect size GWAS and effect size MPRA ----
                
                
                
                
                linearModel_FC <- summary(lm(abs_finemap_z ~ FC, 
                                             data=Element_collapse_CT_CLASS_Lineage_phenotype_sel_ACTIVE))
                
                # cat("linearModel_FC\n")
                # cat(str(linearModel_FC))
                # cat("\n")
                
                
                r.squared_linearModel_FC<-round(linearModel_FC$r.squared,3)
                adj.r.squared_linearModel_FC<-round(linearModel_FC$adj.r.squared,3)
                
                
                
                linearModel_FC_coeffcient_df.m<-melt(linearModel_FC$coefficients)
                
                colnames(linearModel_FC_coeffcient_df.m)[which(colnames(linearModel_FC_coeffcient_df.m)=="Var1")]<-"Terms"
                colnames(linearModel_FC_coeffcient_df.m)[which(colnames(linearModel_FC_coeffcient_df.m)=="Var2")]<-"Parameters"
                
                linearModel_FC_coeffcient_df.m$Terms<-as.character(linearModel_FC_coeffcient_df.m$Terms)
                linearModel_FC_coeffcient_df.m$Parameters<-as.character(linearModel_FC_coeffcient_df.m$Parameters)
                linearModel_FC_coeffcient_df.m$Cell_Type<-Cell_Types_array_sel
                #  linearModel_FC_coeffcient_df.m$LINEAGE_CLASSIF_DEF<-LINEAGE_CLASSIF_DEF_array_sel
                linearModel_FC_coeffcient_df.m$adj.r.squared_linearModel_FC<-adj.r.squared_linearModel_FC
                
                
                
                # cat("linearModel_FC_coeffcient_df.m\n")
                # cat(str(linearModel_FC_coeffcient_df.m))
                # cat("\n")
                
                
                FC_slope<-round(linearModel_FC_coeffcient_df.m$value[which(linearModel_FC_coeffcient_df.m$Parameters== "Estimate" &
                                                                             linearModel_FC_coeffcient_df.m$Terms == "FC")],2)
                
                
                # cat("FC_slope\n")
                # cat(str(FC_slope))
                # cat("\n")
                
                
                FC_intercept<-round(linearModel_FC_coeffcient_df.m$value[which(linearModel_FC_coeffcient_df.m$Parameters== "Estimate" &
                                                                                 linearModel_FC_coeffcient_df.m$Terms == "(Intercept)")],2)
                
                
                # cat("FC_intercept\n")
                # cat(str(FC_intercept))
                # cat("\n")
                # 
                A.df<-as.data.frame(cbind(phenotype_array_sel,"FC",FC_slope,FC_intercept,dim(Element_collapse_CT_CLASS_Lineage_phenotype_sel_ACTIVE)[1],r.squared_linearModel_FC))
                
                # cat("A.df\n")
                # cat(str(A.df))
                # cat("\n")
                # 
                colnames(A.df)<-c("phenotype","MPRA_parameter","LM_slope","LM_intercept","n","r.squared")
                
                # cat("A.df\n")
                # cat(str(A.df))
                # cat("\n")
                
                list_LM_FC_phenotype[[m]]<-A.df
                
                
                # quit(status = 1)
                
              }#dim(Element_collapse_CT_CLASS_Lineage_phenotype_sel_ACTIVE)[1] >0
              
              
              
              
              ### STATS ----
              
              
              CLASS_diversity<-unique(as.character(Element_collapse_CT_CLASS_Lineage_phenotype_sel[,indx.CLASS]))
              
              # cat("CLASS_diversity\n")
              # cat(str(CLASS_diversity))
              # cat("\n")
              
              if(length(CLASS_diversity) >1)
              {
                Element_collapse_CT_CLASS_Lineage_phenotype_sel.dt<-data.table(Element_collapse_CT_CLASS_Lineage_phenotype_sel,
                                                                               key=c(AT_LEAST_CLASSES_sel))
                
                
                # cat("Element_collapse_CT_CLASS_Lineage_phenotype_sel.dt_2\n")
                # cat(str(Element_collapse_CT_CLASS_Lineage_phenotype_sel.dt))
                # cat("\n")
                
                Summary_table_Quartiles<-as.data.frame(Element_collapse_CT_CLASS_Lineage_phenotype_sel.dt[, .(minAbsolute=round(as.numeric(summary(abs_finemap_z)[1]),3),
                                                                                                              FirstQuartile=round(as.numeric(summary(abs_finemap_z)[2]),3),
                                                                                                              median=round(as.numeric(summary(abs_finemap_z)[3]),3),
                                                                                                              ThirdQuartile=round(as.numeric(summary(abs_finemap_z)[5]),3),
                                                                                                              maxAbsolute=round(as.numeric(summary(abs_finemap_z)[6]),3)),
                                                                                                          by=key(Element_collapse_CT_CLASS_Lineage_phenotype_sel.dt)]
                                                       ,stringsAsFactors=F)
                
                
                # cat("Summary_table_Quartiles\n")
                # cat(str(Summary_table_Quartiles))
                # cat("\n")
                
                Summary_table_Quartiles_wide<-as.data.frame(pivot_wider(Summary_table_Quartiles,
                                                                        names_from=colnames(Summary_table_Quartiles)[which(colnames(Summary_table_Quartiles) == AT_LEAST_CLASSES_sel)],
                                                                        values_from=colnames(Summary_table_Quartiles)[-which(colnames(Summary_table_Quartiles) == AT_LEAST_CLASSES_sel)]),stringsAsFactors=F)
                
                # cat("Summary_table_Quartiles_wide\n")
                # cat(str(Summary_table_Quartiles_wide))
                # cat("\n")
                
                #### STATS Kruskall-Wallis data not normally distributed ----
                
                
                Kruskal.Wallis.category<-kruskal.test(abs_finemap_z ~ Element_collapse_CT_CLASS_Lineage_phenotype_sel[,indx.CLASS], data = Element_collapse_CT_CLASS_Lineage_phenotype_sel)
                
                
                # cat("Kruskal.Wallis.category\n")
                # cat(str(Kruskal.Wallis.category))
                # cat("\n")
                
                logpvalue.Kruskal.Wallis.category<-round(-1*log10(Kruskal.Wallis.category$p.value),2)
                #
                # cat("logpvalue.Kruskal.Wallis.category\n")
                # cat(sprintf(as.character(logpvalue.Kruskal.Wallis.category)))
                # cat("\n")
                
                
                colnames(Summary_table_Quartiles_wide)<-gsub("_0","_CTRL",colnames(Summary_table_Quartiles_wide))
                
                # cat("Summary_table_Quartiles_1.5\n")
                # cat(str(Summary_table_Quartiles))
                # cat("\n")
                
                colnames(Summary_table_Quartiles_wide)<-gsub("_[^CTRL]","_COMPARISON",colnames(Summary_table_Quartiles_wide))
                colnames(Summary_table_Quartiles_wide)<-gsub(" or >","",colnames(Summary_table_Quartiles_wide))
                colnames(Summary_table_Quartiles_wide)<-gsub(" ","",colnames(Summary_table_Quartiles_wide))
                
                
                Summary_table_Quartiles_wide$log_pval<-logpvalue.Kruskal.Wallis.category
                
                Summary_table_Quartiles_wide$phenotype<-phenotype_array_sel
                
                
                # cat("Summary_table_Quartiles_wide_2\n")
                # cat(str(Summary_table_Quartiles_wide))
                # cat("\n")
                
                list_STATS_phenotype[[m]]<-Summary_table_Quartiles_wide
                
                # if(logpvalue.Kruskal.Wallis.category >= 1.3)
                # {
                #   quit(status=1)
                #   
                # }
                
              }# CLASS_diversity
            }# dim(Element_collapse_CT_CLASS_Lineage_phenotype_sel)[1] >0
          }#m phenotype_array_sel
          
          
          if(length(list_LM_FC_phenotype) >0)
          {
            LM_FC_phenotype = unique(as.data.frame(data.table::rbindlist(list_LM_FC_phenotype, fill = T)))
            
            LM_FC_phenotype$Lineage<-Lineage_array_sel
            
            cat("LM_FC_phenotype\n")
            cat(str(LM_FC_phenotype))
            cat("\n")
            
            list_LM_FC_Lineage[[l]]<-LM_FC_phenotype
            
          }
          
          if(length(list_LM_abs_Vockley_REF_phenotype) >0)
          {
            LM_abs_Vockley_REF_phenotype = unique(as.data.frame(data.table::rbindlist(list_LM_abs_Vockley_REF_phenotype, fill = T)))
            
            LM_abs_Vockley_REF_phenotype$Lineage<-Lineage_array_sel
            
            cat("LM_abs_Vockley_REF_phenotype\n")
            cat(str(LM_abs_Vockley_REF_phenotype))
            cat("\n")
            
            list_LM_abs_Vockley_REF_Lineage[[l]]<-LM_abs_Vockley_REF_phenotype
            
          }
          
          if(length(list_STATS_phenotype) >0)
          {
            STATS_phenotype = unique(as.data.frame(data.table::rbindlist(list_STATS_phenotype, fill = T)))
            
            STATS_phenotype$Lineage<-Lineage_array_sel
            
            # cat("STATS_phenotype\n")
            # cat(str(STATS_phenotype))
            # cat("\n")
            
            list_STATS_Lineage[[l]]<-STATS_phenotype
          }
          
        }#l Lineage_array
        
        
        if(length(list_LM_FC_Lineage) >0)
        {
          LM_FC_LINEAGE = unique(as.data.frame(data.table::rbindlist(list_LM_FC_Lineage, fill = T)))
          
          LM_FC_LINEAGE$CLASS_COMPARISON<-AT_LEAST_CLASSES_sel
          
          cat("LM_FC_LINEAGE\n")
          cat(str(LM_FC_LINEAGE))
          cat("\n")
          
          list_LM_FC_CLASSES[[k]]<-LM_FC_LINEAGE
        }
        
        if(length(list_LM_abs_Vockley_REF_Lineage) >0)
        {
          LM_abs_Vockley_REF_LINEAGE = unique(as.data.frame(data.table::rbindlist(list_LM_abs_Vockley_REF_Lineage, fill = T)))
          
          LM_abs_Vockley_REF_LINEAGE$CLASS_COMPARISON<-AT_LEAST_CLASSES_sel
          
          cat("LM_abs_Vockley_REF_LINEAGE\n")
          cat(str(LM_abs_Vockley_REF_LINEAGE))
          cat("\n")
          
          list_LM_abs_Vockley_REF_CLASSES[[k]]<-LM_abs_Vockley_REF_LINEAGE
        }
        
        if(length(list_STATS_Lineage) >0)
        {
          STATS_LINEAGE = unique(as.data.frame(data.table::rbindlist(list_STATS_Lineage, fill = T)))
          
          STATS_LINEAGE$CLASS_COMPARISON<-AT_LEAST_CLASSES_sel
          
          # cat("STATS_LINEAGE\n")
          # cat(str(STATS_LINEAGE))
          # cat("\n")
          
          list_STATS_CLASSES[[k]]<-STATS_LINEAGE
        }
      }# dim(Element_collapse_CT_CLASS_sel)[1] >0
      # #############################################################################################################
      # quit(status = 1)
    }# k AT_LEAST_CLASSES_sel
    
    if(length(list_LM_abs_Vockley_REF_CLASSES) >0)
    {
      LM_abs_Vockley_REF_CLASSES = unique(as.data.frame(data.table::rbindlist(list_LM_abs_Vockley_REF_CLASSES, fill = T)))
      
      LM_abs_Vockley_REF_CLASSES$Cell_Type<-Cell_Types_array_sel
      
      cat("LM_abs_Vockley_REF_CLASSES\n")
      cat(str(LM_abs_Vockley_REF_CLASSES))
      cat("\n")
      
      list_LM_abs_Vockley_REF_DEF[[i]]<-LM_abs_Vockley_REF_CLASSES
    }
    
    if(length(list_LM_FC_CLASSES) >0)
    {
      LM_FC_CLASSES = unique(as.data.frame(data.table::rbindlist(list_LM_FC_CLASSES, fill = T)))
      
      LM_FC_CLASSES$Cell_Type<-Cell_Types_array_sel
      
      # cat("LM_FC_CLASSES\n")
      # cat(str(LM_FC_CLASSES))
      # cat("\n")
      
      list_LM_FC_DEF[[i]]<-LM_FC_CLASSES
    }
    
    
    if(length(list_STATS_CLASSES) >0)
    {
      STATS_CLASSES = unique(as.data.frame(data.table::rbindlist(list_STATS_CLASSES, fill = T)))
      
      STATS_CLASSES$Cell_Type<-Cell_Types_array_sel
      
      # cat("STATS_CLASSES\n")
      # cat(str(STATS_CLASSES))
      # cat("\n")
      
      list_STATS_DEF[[i]]<-STATS_CLASSES
    }
    
  }#i Cell_Types_array[i]
  
  if(length(list_LM_abs_Vockley_REF_DEF) >0)
  {
    LM_abs_Vockley_REF_CT = unique(as.data.frame(data.table::rbindlist(list_LM_abs_Vockley_REF_DEF, fill = T)))
    
    
    cat("LM_abs_Vockley_REF_CT\n")
    cat(str(LM_abs_Vockley_REF_CT))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(LM_abs_Vockley_REF_CT$Cell_Type))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(LM_abs_Vockley_REF_CT$Cell_Type)))))
    cat("\n")
    
  }
  
  if(length(list_LM_FC_DEF) >0)
  {
    LM_FC_CT = unique(as.data.frame(data.table::rbindlist(list_LM_FC_DEF, fill = T)))
    
    
    cat("LM_FC_CT\n")
    cat(str(LM_FC_CT))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(LM_FC_CT$Cell_Type))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(LM_FC_CT$Cell_Type)))))
    cat("\n")
    
  }
  
  LM_DEF<-rbind(LM_FC_CT,LM_abs_Vockley_REF_CT)
  cat("LM_DEF\n")
  cat(str(LM_DEF))
  cat("\n")
  
  if(length(list_STATS_DEF) >0)
  {
    STATS_CT = unique(as.data.frame(data.table::rbindlist(list_STATS_DEF, fill = T)))
    
    
    cat("STATS_CT\n")
    cat(str(STATS_CT))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(STATS_CT$Cell_Type))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(STATS_CT$Cell_Type)))))
    cat("\n")
    
  }
  
  
  #### SAVE ----
  
  path5<-paste(out,'GWAS_plots','/', sep='')
  
  cat("path5\n")
  cat(sprintf(as.character(path5)))
  cat("\n")
  
  
  if (file.exists(path5)){
    
    
    
    
  } else {
    dir.create(file.path(path5))
    
  }
  
  
  setwd(path5)
  
  saveRDS(STATS_CT,file="phenotypes_enhancer_CLASS_stats.rds")
  
  saveRDS(LM_DEF,file="phenotypes_enhancer_CLASS_LM.rds")
  
}

ggridge_E_Plus_ASE_CLASS = function(option_list)
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
  
  
  
  path5<-paste(out,'GWAS_plots','/', sep='')
  
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
  
  
  filename_1<-paste("Element_collapse_Plus_Variant_Phenotypes_Plus_LogFC_Vockley_REF",".rds", sep='')
  
  
  Element_collapse<-readRDS(file=filename_1)
  
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
  
  
  
  #### LOOP ggridge ----
  
  
  Cell_Types_array<-levels(Element_collapse$Cell_Type)
  
  cat("Cell_Types_array_2\n")
  cat(str(Cell_Types_array))
  cat("\n")
  
  
  list_STATS_DEF<-list()
  list_LM_FC_DEF<-list()
  list_LM_abs_Vockley_REF_DEF<-list()
  
  
  for(i in 1:length(Cell_Types_array))
  {
    
    Cell_Types_array_sel<-Cell_Types_array[i]
    
    cat("--->\t")
    cat(sprintf(as.character(Cell_Types_array_sel)))
    cat("\n")
    
    
    Element_collapse_CT<-Element_collapse[which(Element_collapse$Cell_Type == Cell_Types_array_sel),]
    
    Element_collapse_CT<-droplevels(Element_collapse_CT)
    
    # cat("Element_collapse_CT\n")
    # cat(str(Element_collapse_CT))
    # cat("\n")
    
    
    AT_LEAST_CLASSES<-c("AT_LEAST_CLASS_1","AT_LEAST_CLASS_2","AT_LEAST_CLASS_3")
    
    list_STATS_CLASSES<-list()
    list_LM_FC_CLASSES<-list()
    list_LM_abs_Vockley_REF_CLASSES<-list()
    
    
    for(k in 1:length(AT_LEAST_CLASSES))
    {
      
      AT_LEAST_CLASSES_sel<-AT_LEAST_CLASSES[k]
      
      # cat("-CLASS COMPARISON-->\t")
      # cat(sprintf(as.character(AT_LEAST_CLASSES_sel)))
      # cat("\t")
      
      indx.CLASS<-which(colnames(Element_collapse_CT) == AT_LEAST_CLASSES_sel)
      
      # cat(sprintf(as.character(indx.CLASS)))
      # cat("\n")
      
      Element_collapse_CT_CLASS_sel<-Element_collapse_CT[!is.na(Element_collapse_CT[,indx.CLASS]),]
      
      cat("Element_collapse_CT_CLASS_sel\n")
      cat(str(Element_collapse_CT_CLASS_sel))
      cat("\n")
      
      if(dim(Element_collapse_CT_CLASS_sel)[1] >0)
      {
        
        Lineage_array<-levels(Element_collapse_CT$Lineage)
        
        cat("Lineage_array\n")
        cat(str(Lineage_array))
        cat("\n")
        
        list_STATS_Lineage<-list()
        list_LM_FC_Lineage<-list()
        list_LM_abs_Vockley_REF_Lineage<-list()
        
        
        for(l in 1:length(Lineage_array))
        {
          Lineage_array_sel<-Lineage_array[l]
          
          cat("--->\t")
          cat(sprintf(as.character(Lineage_array_sel)))
          cat("\n")
          
          Element_collapse_CT_CLASS_Lineage_sel<-Element_collapse_CT_CLASS_sel[which(Element_collapse_CT_CLASS_sel$Lineage == Lineage_array_sel),]
          
          Element_collapse_CT_CLASS_Lineage_sel<-droplevels(Element_collapse_CT_CLASS_Lineage_sel)
          
          cat("Element_collapse_CT_CLASS_Lineage_sel\n")
          cat(str(Element_collapse_CT_CLASS_Lineage_sel))
          cat("\n")
          
          
          phenotype_array<-levels(Element_collapse_CT_CLASS_Lineage_sel$phenotype)
          
          cat("phenotype_array\n")
          cat(str(phenotype_array))
          cat("\n")
          
          list_STATS_phenotype<-list()
          list_LM_FC_phenotype<-list()
          list_LM_abs_Vockley_REF_phenotype<-list()
          
          
          for(m in 1:length(phenotype_array))
          {
            phenotype_array_sel<-phenotype_array[m]
            
            cat("--->\t")
            cat(sprintf(as.character(phenotype_array_sel)))
            cat("\n")
            
            Element_collapse_CT_CLASS_Lineage_phenotype_sel<-Element_collapse_CT_CLASS_Lineage_sel[which(Element_collapse_CT_CLASS_Lineage_sel$phenotype == phenotype_array_sel),]
            
            
            
            if(dim(Element_collapse_CT_CLASS_Lineage_phenotype_sel)[1] >0)
            {
              cat("Element_collapse_CT_CLASS_Lineage_phenotype_sel\n")
              cat(str(Element_collapse_CT_CLASS_Lineage_phenotype_sel))
              cat("\n")
              
              indx.CLASS_2<-which(colnames(Element_collapse_CT_CLASS_Lineage_phenotype_sel) == AT_LEAST_CLASSES_sel)
              
              
              Element_collapse_CT_CLASS_Lineage_phenotype_sel_ACTIVE<-Element_collapse_CT_CLASS_Lineage_phenotype_sel[which(as.numeric(Element_collapse_CT_CLASS_Lineage_phenotype_sel[,indx.CLASS_2]) == 2),]
              
              cat("Element_collapse_CT_CLASS_Lineage_phenotype_sel_ACTIVE\n")
              cat(str(Element_collapse_CT_CLASS_Lineage_phenotype_sel_ACTIVE))
              cat("\n")
              
              
              
              if(dim(Element_collapse_CT_CLASS_Lineage_phenotype_sel_ACTIVE)[1] >= 2)
              {
                
                Element_collapse_CT_CLASS_Lineage_phenotype_sel_ACTIVE$FC<-logratio2foldchange(Element_collapse_CT_CLASS_Lineage_phenotype_sel_ACTIVE$LogFC, base=2)
                
                
                # cat("Element_collapse_CT_CLASS_Lineage_phenotype_sel_ACTIVE_2\n")
                # cat(str(Element_collapse_CT_CLASS_Lineage_phenotype_sel_ACTIVE))
                # cat("\n")
                # quit(status = 1)
                
                ####models of relation effect size GWAS and effect size MPRA ----
                
                linearModel_abs_Vockley_REF <- summary(lm(abs_finemap_z ~ abs_Vockley_REF, 
                                                          data=Element_collapse_CT_CLASS_Lineage_phenotype_sel_ACTIVE))
                
                # cat("linearModel_abs_Vockley_REF\n")
                # cat(str(linearModel_abs_Vockley_REF))
                # cat("\n")
                # 
                
                r.squared_linearModel_abs_Vockley_REF<-round(linearModel_abs_Vockley_REF$r.squared,3)
                adj.r.squared_linearModel_abs_Vockley_REF<-round(linearModel_abs_Vockley_REF$adj.r.squared,3)
                
                
                
                linearModel_abs_Vockley_REF_coeffcient_df.m<-melt(linearModel_abs_Vockley_REF$coefficients)
                
                colnames(linearModel_abs_Vockley_REF_coeffcient_df.m)[which(colnames(linearModel_abs_Vockley_REF_coeffcient_df.m)=="Var1")]<-"Terms"
                colnames(linearModel_abs_Vockley_REF_coeffcient_df.m)[which(colnames(linearModel_abs_Vockley_REF_coeffcient_df.m)=="Var2")]<-"Parameters"
                
                linearModel_abs_Vockley_REF_coeffcient_df.m$Terms<-as.character(linearModel_abs_Vockley_REF_coeffcient_df.m$Terms)
                linearModel_abs_Vockley_REF_coeffcient_df.m$Parameters<-as.character(linearModel_abs_Vockley_REF_coeffcient_df.m$Parameters)
                linearModel_abs_Vockley_REF_coeffcient_df.m$Cell_Type<-Cell_Types_array_sel
                #  linearModel_abs_Vockley_REF_coeffcient_df.m$LINEAGE_CLASSIF_DEF<-LINEAGE_CLASSIF_DEF_array_sel
                linearModel_abs_Vockley_REF_coeffcient_df.m$adj.r.squared_linearModel_abs_Vockley_REF<-adj.r.squared_linearModel_abs_Vockley_REF
                
                
                
                # cat("linearModel_abs_Vockley_REF_coeffcient_df.m\n")
                # cat(str(linearModel_abs_Vockley_REF_coeffcient_df.m))
                # cat("\n")
                # 
                
                FC_slope<-round(linearModel_abs_Vockley_REF_coeffcient_df.m$value[which(linearModel_abs_Vockley_REF_coeffcient_df.m$Parameters== "Estimate" &
                                                                                          linearModel_abs_Vockley_REF_coeffcient_df.m$Terms == "abs_Vockley_REF")],2)
                
                
                # cat("FC_slope\n")
                # cat(str(FC_slope))
                # cat("\n")
                
                
                FC_intercept<-round(linearModel_abs_Vockley_REF_coeffcient_df.m$value[which(linearModel_abs_Vockley_REF_coeffcient_df.m$Parameters== "Estimate" &
                                                                                              linearModel_abs_Vockley_REF_coeffcient_df.m$Terms == "(Intercept)")],2)
                
                
                # cat("FC_intercept\n")
                # cat(str(FC_intercept))
                # cat("\n")
                
                A.df<-as.data.frame(cbind(phenotype_array_sel,"abs_Vockley_REF",FC_slope,FC_intercept,dim(Element_collapse_CT_CLASS_Lineage_phenotype_sel_ACTIVE)[1],r.squared_linearModel_abs_Vockley_REF))
                
                # cat("A.df\n")
                # cat(str(A.df))
                # cat("\n")
                
                colnames(A.df)<-c("phenotype","MPRA_parameter","LM_slope","LM_intercept","n","r.squared")
                
                # cat("A.df\n")
                # cat(str(A.df))
                # cat("\n")
                
                list_LM_abs_Vockley_REF_phenotype[[m]]<-A.df
                
                
                
                # quit(status = 1)
                
                
                ####models of relation effect size GWAS and effect size MPRA ----
                
                
                
                
                linearModel_FC <- summary(lm(abs_finemap_z ~ FC, 
                                             data=Element_collapse_CT_CLASS_Lineage_phenotype_sel_ACTIVE))
                
                # cat("linearModel_FC\n")
                # cat(str(linearModel_FC))
                # cat("\n")
                
                
                r.squared_linearModel_FC<-round(linearModel_FC$r.squared,3)
                adj.r.squared_linearModel_FC<-round(linearModel_FC$adj.r.squared,3)
                
                
                
                linearModel_FC_coeffcient_df.m<-melt(linearModel_FC$coefficients)
                
                colnames(linearModel_FC_coeffcient_df.m)[which(colnames(linearModel_FC_coeffcient_df.m)=="Var1")]<-"Terms"
                colnames(linearModel_FC_coeffcient_df.m)[which(colnames(linearModel_FC_coeffcient_df.m)=="Var2")]<-"Parameters"
                
                linearModel_FC_coeffcient_df.m$Terms<-as.character(linearModel_FC_coeffcient_df.m$Terms)
                linearModel_FC_coeffcient_df.m$Parameters<-as.character(linearModel_FC_coeffcient_df.m$Parameters)
                linearModel_FC_coeffcient_df.m$Cell_Type<-Cell_Types_array_sel
                #  linearModel_FC_coeffcient_df.m$LINEAGE_CLASSIF_DEF<-LINEAGE_CLASSIF_DEF_array_sel
                linearModel_FC_coeffcient_df.m$adj.r.squared_linearModel_FC<-adj.r.squared_linearModel_FC
                
                
                
                # cat("linearModel_FC_coeffcient_df.m\n")
                # cat(str(linearModel_FC_coeffcient_df.m))
                # cat("\n")
                
                
                FC_slope<-round(linearModel_FC_coeffcient_df.m$value[which(linearModel_FC_coeffcient_df.m$Parameters== "Estimate" &
                                                                             linearModel_FC_coeffcient_df.m$Terms == "FC")],2)
                
                
                # cat("FC_slope\n")
                # cat(str(FC_slope))
                # cat("\n")
                
                
                FC_intercept<-round(linearModel_FC_coeffcient_df.m$value[which(linearModel_FC_coeffcient_df.m$Parameters== "Estimate" &
                                                                                 linearModel_FC_coeffcient_df.m$Terms == "(Intercept)")],2)
                
                
                # cat("FC_intercept\n")
                # cat(str(FC_intercept))
                # cat("\n")
                # 
                A.df<-as.data.frame(cbind(phenotype_array_sel,"FC",FC_slope,FC_intercept,dim(Element_collapse_CT_CLASS_Lineage_phenotype_sel_ACTIVE)[1],r.squared_linearModel_FC))
                
                # cat("A.df\n")
                # cat(str(A.df))
                # cat("\n")
                # 
                colnames(A.df)<-c("phenotype","MPRA_parameter","LM_slope","LM_intercept","n","r.squared")
                
                # cat("A.df\n")
                # cat(str(A.df))
                # cat("\n")
                
                list_LM_FC_phenotype[[m]]<-A.df
                
                
                # quit(status = 1)
                
              }#dim(Element_collapse_CT_CLASS_Lineage_phenotype_sel_ACTIVE)[1] >0
              
              
              
              
              ### STATS ----
              
              
              CLASS_diversity<-unique(as.character(Element_collapse_CT_CLASS_Lineage_phenotype_sel[,indx.CLASS]))
              
              # cat("CLASS_diversity\n")
              # cat(str(CLASS_diversity))
              # cat("\n")
              
              if(length(CLASS_diversity) >1)
              {
                Element_collapse_CT_CLASS_Lineage_phenotype_sel.dt<-data.table(Element_collapse_CT_CLASS_Lineage_phenotype_sel,
                                                                               key=c(AT_LEAST_CLASSES_sel))
                
                
                # cat("Element_collapse_CT_CLASS_Lineage_phenotype_sel.dt_2\n")
                # cat(str(Element_collapse_CT_CLASS_Lineage_phenotype_sel.dt))
                # cat("\n")
                
                Summary_table_Quartiles<-as.data.frame(Element_collapse_CT_CLASS_Lineage_phenotype_sel.dt[, .(minAbsolute=round(as.numeric(summary(abs_finemap_z)[1]),3),
                                                                                                              FirstQuartile=round(as.numeric(summary(abs_finemap_z)[2]),3),
                                                                                                              median=round(as.numeric(summary(abs_finemap_z)[3]),3),
                                                                                                              ThirdQuartile=round(as.numeric(summary(abs_finemap_z)[5]),3),
                                                                                                              maxAbsolute=round(as.numeric(summary(abs_finemap_z)[6]),3)),
                                                                                                          by=key(Element_collapse_CT_CLASS_Lineage_phenotype_sel.dt)]
                                                       ,stringsAsFactors=F)
                
                
                # cat("Summary_table_Quartiles\n")
                # cat(str(Summary_table_Quartiles))
                # cat("\n")
                
                Summary_table_Quartiles_wide<-as.data.frame(pivot_wider(Summary_table_Quartiles,
                                                                        names_from=colnames(Summary_table_Quartiles)[which(colnames(Summary_table_Quartiles) == AT_LEAST_CLASSES_sel)],
                                                                        values_from=colnames(Summary_table_Quartiles)[-which(colnames(Summary_table_Quartiles) == AT_LEAST_CLASSES_sel)]),stringsAsFactors=F)
                
                # cat("Summary_table_Quartiles_wide\n")
                # cat(str(Summary_table_Quartiles_wide))
                # cat("\n")
                
                #### STATS Kruskall-Wallis data not normally distributed ----
                
                
                Kruskal.Wallis.category<-kruskal.test(abs_finemap_z ~ Element_collapse_CT_CLASS_Lineage_phenotype_sel[,indx.CLASS], data = Element_collapse_CT_CLASS_Lineage_phenotype_sel)
                
                
                # cat("Kruskal.Wallis.category\n")
                # cat(str(Kruskal.Wallis.category))
                # cat("\n")
                
                logpvalue.Kruskal.Wallis.category<-round(-1*log10(Kruskal.Wallis.category$p.value),2)
                #
                # cat("logpvalue.Kruskal.Wallis.category\n")
                # cat(sprintf(as.character(logpvalue.Kruskal.Wallis.category)))
                # cat("\n")
                
                
                colnames(Summary_table_Quartiles_wide)<-gsub("_0","_CTRL",colnames(Summary_table_Quartiles_wide))
                
                # cat("Summary_table_Quartiles_1.5\n")
                # cat(str(Summary_table_Quartiles))
                # cat("\n")
                
                colnames(Summary_table_Quartiles_wide)<-gsub("_[^CTRL]","_COMPARISON",colnames(Summary_table_Quartiles_wide))
                colnames(Summary_table_Quartiles_wide)<-gsub(" or >","",colnames(Summary_table_Quartiles_wide))
                colnames(Summary_table_Quartiles_wide)<-gsub(" ","",colnames(Summary_table_Quartiles_wide))
                
                
                Summary_table_Quartiles_wide$log_pval<-logpvalue.Kruskal.Wallis.category
                
                Summary_table_Quartiles_wide$phenotype<-phenotype_array_sel
                
                
                # cat("Summary_table_Quartiles_wide_2\n")
                # cat(str(Summary_table_Quartiles_wide))
                # cat("\n")
                
                list_STATS_phenotype[[m]]<-Summary_table_Quartiles_wide
                
                # if(logpvalue.Kruskal.Wallis.category >= 1.3)
                # {
                #   quit(status=1)
                #   
                # }
                
              }# CLASS_diversity
            }# dim(Element_collapse_CT_CLASS_Lineage_phenotype_sel)[1] >0
          }#m phenotype_array_sel
          
          
          if(length(list_LM_FC_phenotype) >0)
          {
            LM_FC_phenotype = unique(as.data.frame(data.table::rbindlist(list_LM_FC_phenotype, fill = T)))
            
            LM_FC_phenotype$Lineage<-Lineage_array_sel
            
            cat("LM_FC_phenotype\n")
            cat(str(LM_FC_phenotype))
            cat("\n")
            
            list_LM_FC_Lineage[[l]]<-LM_FC_phenotype
            
          }
          
          if(length(list_LM_abs_Vockley_REF_phenotype) >0)
          {
            LM_abs_Vockley_REF_phenotype = unique(as.data.frame(data.table::rbindlist(list_LM_abs_Vockley_REF_phenotype, fill = T)))
            
            LM_abs_Vockley_REF_phenotype$Lineage<-Lineage_array_sel
            
            cat("LM_abs_Vockley_REF_phenotype\n")
            cat(str(LM_abs_Vockley_REF_phenotype))
            cat("\n")
            
            list_LM_abs_Vockley_REF_Lineage[[l]]<-LM_abs_Vockley_REF_phenotype
            
          }
          
          if(length(list_STATS_phenotype) >0)
          {
            STATS_phenotype = unique(as.data.frame(data.table::rbindlist(list_STATS_phenotype, fill = T)))
            
            STATS_phenotype$Lineage<-Lineage_array_sel
            
            # cat("STATS_phenotype\n")
            # cat(str(STATS_phenotype))
            # cat("\n")
            
            list_STATS_Lineage[[l]]<-STATS_phenotype
          }
          
        }#l Lineage_array
        
        
        if(length(list_LM_FC_Lineage) >0)
        {
          LM_FC_LINEAGE = unique(as.data.frame(data.table::rbindlist(list_LM_FC_Lineage, fill = T)))
          
          LM_FC_LINEAGE$CLASS_COMPARISON<-AT_LEAST_CLASSES_sel
          
          cat("LM_FC_LINEAGE\n")
          cat(str(LM_FC_LINEAGE))
          cat("\n")
          
          list_LM_FC_CLASSES[[k]]<-LM_FC_LINEAGE
        }
        
        if(length(list_LM_abs_Vockley_REF_Lineage) >0)
        {
          LM_abs_Vockley_REF_LINEAGE = unique(as.data.frame(data.table::rbindlist(list_LM_abs_Vockley_REF_Lineage, fill = T)))
          
          LM_abs_Vockley_REF_LINEAGE$CLASS_COMPARISON<-AT_LEAST_CLASSES_sel
          
          cat("LM_abs_Vockley_REF_LINEAGE\n")
          cat(str(LM_abs_Vockley_REF_LINEAGE))
          cat("\n")
          
          list_LM_abs_Vockley_REF_CLASSES[[k]]<-LM_abs_Vockley_REF_LINEAGE
        }
        
        if(length(list_STATS_Lineage) >0)
        {
          STATS_LINEAGE = unique(as.data.frame(data.table::rbindlist(list_STATS_Lineage, fill = T)))
          
          STATS_LINEAGE$CLASS_COMPARISON<-AT_LEAST_CLASSES_sel
          
          # cat("STATS_LINEAGE\n")
          # cat(str(STATS_LINEAGE))
          # cat("\n")
          
          list_STATS_CLASSES[[k]]<-STATS_LINEAGE
        }
      }# dim(Element_collapse_CT_CLASS_sel)[1] >0
      # #############################################################################################################
      # quit(status = 1)
    }# k AT_LEAST_CLASSES_sel
    
    if(length(list_LM_abs_Vockley_REF_CLASSES) >0)
    {
      LM_abs_Vockley_REF_CLASSES = unique(as.data.frame(data.table::rbindlist(list_LM_abs_Vockley_REF_CLASSES, fill = T)))
      
      LM_abs_Vockley_REF_CLASSES$Cell_Type<-Cell_Types_array_sel
      
      cat("LM_abs_Vockley_REF_CLASSES\n")
      cat(str(LM_abs_Vockley_REF_CLASSES))
      cat("\n")
      
      list_LM_abs_Vockley_REF_DEF[[i]]<-LM_abs_Vockley_REF_CLASSES
    }
    
    if(length(list_LM_FC_CLASSES) >0)
    {
      LM_FC_CLASSES = unique(as.data.frame(data.table::rbindlist(list_LM_FC_CLASSES, fill = T)))
      
      LM_FC_CLASSES$Cell_Type<-Cell_Types_array_sel
      
      # cat("LM_FC_CLASSES\n")
      # cat(str(LM_FC_CLASSES))
      # cat("\n")
      
      list_LM_FC_DEF[[i]]<-LM_FC_CLASSES
    }
    
    
    if(length(list_STATS_CLASSES) >0)
    {
      STATS_CLASSES = unique(as.data.frame(data.table::rbindlist(list_STATS_CLASSES, fill = T)))
      
      STATS_CLASSES$Cell_Type<-Cell_Types_array_sel
      
      # cat("STATS_CLASSES\n")
      # cat(str(STATS_CLASSES))
      # cat("\n")
      
      list_STATS_DEF[[i]]<-STATS_CLASSES
    }
    
  }#i Cell_Types_array[i]
  
  if(length(list_LM_abs_Vockley_REF_DEF) >0)
  {
    LM_abs_Vockley_REF_CT = unique(as.data.frame(data.table::rbindlist(list_LM_abs_Vockley_REF_DEF, fill = T)))
    
    
    cat("LM_abs_Vockley_REF_CT\n")
    cat(str(LM_abs_Vockley_REF_CT))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(LM_abs_Vockley_REF_CT$Cell_Type))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(LM_abs_Vockley_REF_CT$Cell_Type)))))
    cat("\n")
    
  }
  
  if(length(list_LM_FC_DEF) >0)
  {
    LM_FC_CT = unique(as.data.frame(data.table::rbindlist(list_LM_FC_DEF, fill = T)))
    
    
    cat("LM_FC_CT\n")
    cat(str(LM_FC_CT))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(LM_FC_CT$Cell_Type))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(LM_FC_CT$Cell_Type)))))
    cat("\n")
    
  }
  
  LM_DEF<-rbind(LM_FC_CT,LM_abs_Vockley_REF_CT)
  cat("LM_DEF\n")
  cat(str(LM_DEF))
  cat("\n")
  
  if(length(list_STATS_DEF) >0)
  {
    STATS_CT = unique(as.data.frame(data.table::rbindlist(list_STATS_DEF, fill = T)))
    
    
    cat("STATS_CT\n")
    cat(str(STATS_CT))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(STATS_CT$Cell_Type))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(STATS_CT$Cell_Type)))))
    cat("\n")
    
  }
  
  
  #### SAVE ----
  
  path5<-paste(out,'GWAS_plots','/', sep='')
  
  cat("path5\n")
  cat(sprintf(as.character(path5)))
  cat("\n")
  
  
  if (file.exists(path5)){
    
    
    
    
  } else {
    dir.create(file.path(path5))
    
  }
  
  
  setwd(path5)
  
  saveRDS(STATS_CT,file="phenotypes_E_Plus_ASE_CLASS_stats.rds")
  
  saveRDS(LM_DEF,file="phenotypes_E_Plus_ASE_CLASS_LM.rds")
  
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
  
  
  path5<-paste(out,'GWAS_plots','/', sep='')
  
  cat("path5\n")
  cat(sprintf(as.character(path5)))
  cat("\n")
  
  
  if (file.exists(path5)){
    
    
    
    
  } else {
    dir.create(file.path(path5))
    
  }
  
  
  
  setwd(path5)
  
  
  #### READ DATA ----
  
  
  
  enhancer_CLASS_stats<-readRDS(file="phenotypes_enhancer_CLASS_stats.rds")
  enhancer_CLASS_stats$ER_Type<-"enhancer_CLASS"
  
  cat("enhancer_CLASS_stats\n")
  cat(str(enhancer_CLASS_stats))
  cat("\n")
  
  
  E_Plus_ASE_CLASS_stats<-readRDS(file="phenotypes_E_Plus_ASE_CLASS_stats.rds")
  E_Plus_ASE_CLASS_stats$ER_Type<-"E_Plus_ASE_CLASS"
  
  cat("E_Plus_ASE_CLASS_stats\n")
  cat(str(E_Plus_ASE_CLASS_stats))
  cat("\n")
  
  
  DEF_stats<-rbind(enhancer_CLASS_stats,
                   E_Plus_ASE_CLASS_stats)
  
  cat("DEF_stats\n")
  cat(str(DEF_stats))
  cat("\n")
  
  #### LM
  
 
  setwd(path5)
  

  
  enhancer_CLASS_LM<-readRDS(file="phenotypes_enhancer_CLASS_LM.rds")
  enhancer_CLASS_LM$ER_Type<-"enhancer_CLASS"
  
  cat("enhancer_CLASS_LM\n")
  cat(str(enhancer_CLASS_LM))
  cat("\n")
  
  
  E_Plus_ASE_CLASS_LM<-readRDS(file="phenotypes_E_Plus_ASE_CLASS_LM.rds")
  E_Plus_ASE_CLASS_LM$ER_Type<-"E_Plus_ASE_CLASS"
  
  cat("E_Plus_ASE_CLASS_LM\n")
  cat(str(E_Plus_ASE_CLASS_LM))
  cat("\n")
  
  
  DEF_LM<-rbind(enhancer_CLASS_LM,
                   E_Plus_ASE_CLASS_LM)
  
  cat("DEF_LM\n")
  cat(str(DEF_LM))
  cat("\n")
  
  
  setwd(path5)
  
  write.table(DEF_LM,file="LM_object.tsv",sep="\t",quote=F,row.names = F)
  
  
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
  
  
  data_wrangling(opt)
  data_wrangling_FC_Vockley_REF(opt)
  ggridge_enhancer_CLASS(opt)
  ggridge_E_Plus_ASE_CLASS(opt)
  print_join(opt)
 
}
  
  
  
 

###########################################################################

system.time( main() )
