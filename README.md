# MPRA_Programmed_pipeline


############# SAMPLE PROCESSING


command line:
 
$ bash Sample_processing.sh sample_file_GLOBAL.txt \<path to leave output files\> \<batch\> \<path to the MPRA_programmed_pipeline\> \<path to the original fastq files produced by bcl2fastq\> \<mem\> \<processors\> \<bsub queue\>

example:

$ nohup bash ~/Scripts/MPRA_programmed_pipeline/Sample_processing.sh /nfs/team151_data03/MPRA_Programmed/sample_file_GLOBAL.txt /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST/ HT_TEST /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/ /nfs/team151_data03/MPRA_Programmed/ 16000 4 normal &

Memory requirements: 16-48 GB  (for R16_gDNA)

## sample_file_GLOBAL.txt carries the information of how to merge the fastqs and where are the UMIs and the R1 and R2

head -4 sample_file_GLOBAL.txt

/nfs/team151_data03/MPRA_Programmed/R10_6_gDNA_S1_L001_R1_001.fastq.gz  /nfs/team151_data03/MPRA_Programmed/R10_6_gDNA_S1_L002_R1_001.fastq.gz  /nfs/team151_data03/MPRA_Programmed/R10_6_gDNA_R1.fastq.gz
/nfs/team151_data03/MPRA_Programmed/R10_6_gDNA_S1_L001_R2_001.fastq.gz  /nfs/team151_data03/MPRA_Programmed/R10_6_gDNA_S1_L002_R2_001.fastq.gz  /nfs/team151_data03/MPRA_Programmed/R10_6_gDNA_R2.fastq.gz
/nfs/team151_data03/MPRA_Programmed/R10_6_gDNA_S1_L001_R3_001.fastq.gz  /nfs/team151_data03/MPRA_Programmed/R10_6_gDNA_S1_L002_R3_001.fastq.gz  /nfs/team151_data03/MPRA_Programmed/R10_6_gDNA_R3.fastq.gz
/nfs/team151_data03/MPRA_Programmed/R10_6_gDNA_S1_L001_I1_001.fastq.gz  /nfs/team151_data03/MPRA_Programmed/R10_6_gDNA_S1_L002_I1_001.fastq.gz  /nfs/team151_data03/MPRA_Programmed/R10_6_gDNA_I1.fastq.gz

for R10_6_gDNA R10_6_gDNA_R2.fastq.gz is the file with UMIs . Umi-tools extract would take UMIs from this file and add it to /nfs/team151_data03/MPRA_Programmed/R10_6_gDNA_R1.fastq.gz and to /nfs/team151_data03/MPRA_Programmed/R10_6_gDNA_R3.fastq.gz. /nfs/team151_data03/MPRA_Programmed/R10_6_gDNA_I1.fastq.gz has the sample indexes but the files have already been demultiplexed

- see an example of sample_file in examples folder


## SOFTWARE DEPENDENCIES FOR SAMPLE PROCESSING

-Trimmomatic-0.39 (change the path in the script)
-umi_tools version: 1.1.1 (change the path in the script)
-flash2 v2.2.00
-bwa
-samtools
-bamtools 2.5.1 (change the path in the script, in this case it is a conda environment)


## The barcode reference to align the fastq "Library_TRIMMED_15bp_Carried_Variants.fasta" is in the Dependencies folder with the indexes needed by bwa.

############# MPRA analysis

$ bash MPRA_pipeline.sh <path to leave output files> <path to the MPRA_programmed_pipeline> <path to the FINAL output files from sample_processing> <mem> <processors> <bsub queue>

$ bash MPRA_pipeline.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/ /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/ /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST/ 16000 4 normal

## FILE DEPENDENCIES

untar the ALL_db.tsv.tar.gz in the Dependencies folder

## SOFTWARE DEPENDENCIES

-R-4.1.0, you have to set the path to Rscript at the top of the bash script: Rscript=/software/R-4.1.0/bin/Rscript
- Rpackages. Important! change the default lib.loc in every script as it direct to my folder in the Sanger farm. See per script but all told they are:

optparse
MPRAnalyze
reshape2
labeling
ggrepel
ggeasy
farver
cowplot
zoo
withr
viridisLite
tzdb
TxDb.Hsapiens.UCSC.hg19.knownGene
ttutils
tidyverse
sysfonts
svglite
Sushi
SummarizedExperiment
showtextdb
showtext
S4Vectors
R.utils
rtracklayer
rstudioapi
R.oo
R.methodsS3
plyr
org.Hs.eg.db
OrganismDbi
memoise
matrixStats
MatrixGenerics
markdown
liftOverlib.loc
lemon
labeling
jsonlite
IRanges
Homo.sapiens
gwascat
gtools
gridExtra
gridBase
grid
GO.db
glue
ggtext
ggridges
ggplot2
ggforce
ggeasy
GenomicRanges
GenomicFeatures
GenomeInfoDb
farver
extrafontdb
extrafont
dplyr
digest
digest
data.table
curl
crayon
cowplot
cli
broom
biomaRt
BiocGenerics
Biobase
backports
AnnotationDbi
 
