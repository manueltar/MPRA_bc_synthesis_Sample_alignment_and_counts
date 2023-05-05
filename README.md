# MPRA_bc_synthesis_Sample_alignment_and_counts Sample Preparation


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
-bamtools 2.5.1


## The barcode reference to align the fastq "Library_TRIMMED_15bp_Carried_Variants.fasta" is in the Dependencies folder with the indexes needed by bwa.