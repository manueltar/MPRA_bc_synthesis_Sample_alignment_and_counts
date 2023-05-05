#!/bin/bash
bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing//Lane_merger_CHRF_R1_gDNA_CHRF_gDNA_R1.out -q normal -n5 -J 141_CHRF_R1_gDNA_1 -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"cat /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_gDNA_S4_L001_R1_001.fastq.gz /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_gDNA_S4_L002_R1_001.fastq.gz  > /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_gDNA_R1.fastq.gz"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing//Lane_merger_CHRF_R1_gDNA_CHRF_gDNA_R2.out -q normal -n5 -J 141_CHRF_R1_gDNA_2 -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"cat /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_gDNA_S4_L001_R2_001.fastq.gz /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_gDNA_S4_L002_R2_001.fastq.gz  > /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_gDNA_R2.fastq.gz"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing//Lane_merger_CHRF_R1_gDNA_CHRF_gDNA_R3.out -q normal -n5 -J 141_CHRF_R1_gDNA_3 -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"cat /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_gDNA_S4_L001_R3_001.fastq.gz /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_gDNA_S4_L002_R3_001.fastq.gz  > /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_gDNA_R3.fastq.gz"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing//Lane_merger_CHRF_R1_gDNA_CHRF_gDNA_I1.out -q normal -n5 -J 141_CHRF_R1_gDNA_4 -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"cat /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_gDNA_S4_L001_I1_001.fastq.gz /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_gDNA_S4_L002_I1_001.fastq.gz  > /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_gDNA_I1.fastq.gz"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing//CHRF_R1_gDNA_Processing.out -q normal -n5 -w"done(141_CHRF_R1_gDNA_1) && done(141_CHRF_R1_gDNA_2) && done(141_CHRF_R1_gDNA_3) && done(141_CHRF_R1_gDNA_4)" -J 141_Trimmomatic_CHRF_R1_gDNA_job -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"java -jar /nfs/users/nfs_m/mt19/sOFTWARE/Trimmomatic-0.39/trimmomatic-0.39.jar SE /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_gDNA_R3.fastq.gz /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/Trimmed_CHRF_R1_gDNA_R3.fastq.gz HEADCROP:6"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_Extracting_UMI_R1_CHRF_R1_gDNA.out -M 20000 -w"done(141_Trimmomatic_CHRF_R1_gDNA_job)" -J 141_Extracting_UMI_R1_CHRF_R1_gDNA_job -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -n5 -q normal -- \
"/nfs/team151/software/umi-tools/bin/umi_tools extract --bc-pattern=NNNNNNNNNN --stdin=/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_gDNA_R2.fastq.gz --read2-in=/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_gDNA_R1.fastq.gz --stdout=Output_folder/CHRF_R1_gDNA_R1_Plus_UMIS.fastq.gz --read2-stdout"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_Extracting_UMI_R3_CHRF_R1_gDNA.out -M 20000 -w"done(141_Trimmomatic_CHRF_R1_gDNA_job)" -J 141_Extracting_UMI_R3_CHRF_R1_gDNA_job -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -n5 -q normal -- \
"/nfs/team151/software/umi-tools/bin/umi_tools extract --bc-pattern=NNNNNNNNNN --stdin=/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_gDNA_R2.fastq.gz --read2-in=/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/Trimmed_CHRF_R1_gDNA_R3.fastq.gz --stdout=Output_folder/CHRF_R1_gDNA_R3_Plus_UMIS.fastq.gz --read2-stdout"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_R1_R3_merging_CHRF_R1_gDNA.out -q normal -n5 -w"done(141_Extracting_UMI_R1_CHRF_R1_gDNA_job) && done(141_Extracting_UMI_R3_CHRF_R1_gDNA_job)" -J 141_R1_R3_merging_CHRF_R1_gDNA_job -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"flash2 -z --output-prefix=CHRF_R1_gDNA -d Output_folder/ -t 1 -r 11 -f 11 -s 2 -m 10 -x 0.122 Output_folder/CHRF_R1_gDNA_R1_Plus_UMIS.fastq.gz Output_folder/CHRF_R1_gDNA_R3_Plus_UMIS.fastq.gz"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_R1_R3_merging_CHRF_R1_gDNA.out -M 4000  -w"done(141_R1_R3_merging_CHRF_R1_gDNA_job)" -J 141_CLEANING_1_CHRF_R1_gDNA_job  -R"select[mem>=4000] rusage[mem=4000] span[hosts=1]" -n1 -q normal -- \
"rm /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_gDNA_R1.fastq.gz /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/Trimmed_CHRF_R1_gDNA_R3.fastq.gz /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_gDNA_R3.fastq.gz"


#########################################################################################################################################################################
#########################################################################################################################################################################
bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_aligning_CHRF_R1_gDNA.out -q normal -n5 -w"done(141_R1_R3_merging_CHRF_R1_gDNA_job)" -J 141_aligning_CHRF_R1_gDNA_job -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"bwa aln -l 10 -O 100 -E 100 /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Dependencies/Library_TRIMMED_15bp_Carried_Variants.fasta Output_folder/CHRF_R1_gDNA.extendedFrags.fastq.gz > Output_folder/CHRF_R1_gDNA.sai"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_converting_and_filtering_CHRF_R1_gDNA.out -q normal -n5 -w"done(141_aligning_CHRF_R1_gDNA_job)" -J 141_converting_and_filtering_CHRF_R1_gDNA_job -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"bwa samse /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Dependencies/Library_TRIMMED_15bp_Carried_Variants.fasta Output_folder/CHRF_R1_gDNA.sai Output_folder/CHRF_R1_gDNA.extendedFrags.fastq.gz > Output_folder/CHRF_R1_gDNA.sam"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_converting_and_filtering_CHRF_R1_gDNA.out -q normal -n5 -w"done(141_converting_and_filtering_CHRF_R1_gDNA_job)" -J 141_sam_to_bam_CHRF_R1_gDNA_job -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"samtools view -@5 -h -F 256 Output_folder/CHRF_R1_gDNA.sam|samtools sort -@5 -o Output_folder/CHRF_R1_gDNA.bam"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_converting_and_filtering_CHRF_R1_gDNA.out -q normal -n5 -w"done(141_sam_to_bam_CHRF_R1_gDNA_job)" -J 141_filter_and_rest_CHRF_R1_gDNA_job -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"bash /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Bash_subscript/268_MPRA_Sample_procesing_UMI_TOOLS_v2_sub_bash.sh Output_folder/ /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/ 141 20000 5 normal CHRF_R1_gDNA CHRF_R1_gDNA"


#########################################################################################################################################################################
#########################################################################################################################################################################
#########################################################################################################################################################################
#########################################################################################################################################################################




bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing//Lane_merger_CHRF_R1_cDNA_CHRF_cDNA_R1.out -q normal -n5 -J 141_CHRF_R1_cDNA_1 -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"cat /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_cDNA_S9_L001_R1_001.fastq.gz /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_cDNA_S9_L002_R1_001.fastq.gz  > /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_cDNA_R1.fastq.gz"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing//Lane_merger_CHRF_R1_cDNA_CHRF_cDNA_R2.out -q normal -n5 -J 141_CHRF_R1_cDNA_2 -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"cat /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_cDNA_S9_L001_R2_001.fastq.gz /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_cDNA_S9_L002_R2_001.fastq.gz  > /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_cDNA_R2.fastq.gz"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing//Lane_merger_CHRF_R1_cDNA_CHRF_cDNA_R3.out -q normal -n5 -J 141_CHRF_R1_cDNA_3 -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"cat /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_cDNA_S9_L001_R3_001.fastq.gz /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_cDNA_S9_L002_R3_001.fastq.gz  > /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_cDNA_R3.fastq.gz"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing//Lane_merger_CHRF_R1_cDNA_CHRF_cDNA_I1.out -q normal -n5 -J 141_CHRF_R1_cDNA_4 -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"cat /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_cDNA_S9_L001_I1_001.fastq.gz /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_cDNA_S9_L002_I1_001.fastq.gz  > /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_cDNA_I1.fastq.gz"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing//CHRF_R1_cDNA_Processing.out -q normal -n5 -w"done(141_CHRF_R1_cDNA_1) && done(141_CHRF_R1_cDNA_2) && done(141_CHRF_R1_cDNA_3) && done(141_CHRF_R1_cDNA_4)" -J 141_Trimmomatic_CHRF_R1_cDNA_job -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"java -jar /nfs/users/nfs_m/mt19/sOFTWARE/Trimmomatic-0.39/trimmomatic-0.39.jar SE /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_cDNA_R3.fastq.gz /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/Trimmed_CHRF_R1_cDNA_R3.fastq.gz HEADCROP:6"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_Extracting_UMI_R1_CHRF_R1_cDNA.out -M 20000 -w"done(141_Trimmomatic_CHRF_R1_cDNA_job)" -J 141_Extracting_UMI_R1_CHRF_R1_cDNA_job -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -n5 -q normal -- \
"/nfs/team151/software/umi-tools/bin/umi_tools extract --bc-pattern=NNNNNNNNNN --stdin=/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_cDNA_R2.fastq.gz --read2-in=/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_cDNA_R1.fastq.gz --stdout=Output_folder/CHRF_R1_cDNA_R1_Plus_UMIS.fastq.gz --read2-stdout"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_Extracting_UMI_R3_CHRF_R1_cDNA.out -M 20000 -w"done(141_Trimmomatic_CHRF_R1_cDNA_job)" -J 141_Extracting_UMI_R3_CHRF_R1_cDNA_job -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -n5 -q normal -- \
"/nfs/team151/software/umi-tools/bin/umi_tools extract --bc-pattern=NNNNNNNNNN --stdin=/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_cDNA_R2.fastq.gz --read2-in=/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/Trimmed_CHRF_R1_cDNA_R3.fastq.gz --stdout=Output_folder/CHRF_R1_cDNA_R3_Plus_UMIS.fastq.gz --read2-stdout"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_R1_R3_merging_CHRF_R1_cDNA.out -q normal -n5 -w"done(141_Extracting_UMI_R1_CHRF_R1_cDNA_job) && done(141_Extracting_UMI_R3_CHRF_R1_cDNA_job)" -J 141_R1_R3_merging_CHRF_R1_cDNA_job -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"flash2 -z --output-prefix=CHRF_R1_cDNA -d Output_folder/ -t 1 -r 11 -f 11 -s 2 -m 10 -x 0.122 Output_folder/CHRF_R1_cDNA_R1_Plus_UMIS.fastq.gz Output_folder/CHRF_R1_cDNA_R3_Plus_UMIS.fastq.gz"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_R1_R3_merging_CHRF_R1_cDNA.out -M 4000  -w"done(141_R1_R3_merging_CHRF_R1_cDNA_job)" -J 141_CLEANING_1_CHRF_R1_cDNA_job  -R"select[mem>=4000] rusage[mem=4000] span[hosts=1]" -n1 -q normal -- \
"rm /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_cDNA_R1.fastq.gz /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/Trimmed_CHRF_R1_cDNA_R3.fastq.gz /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R1_cDNA_R3.fastq.gz"


#########################################################################################################################################################################
#########################################################################################################################################################################
bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_aligning_CHRF_R1_cDNA.out -q normal -n5 -w"done(141_R1_R3_merging_CHRF_R1_cDNA_job)" -J 141_aligning_CHRF_R1_cDNA_job -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"bwa aln -l 10 -O 100 -E 100 /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Dependencies/Library_TRIMMED_15bp_Carried_Variants.fasta Output_folder/CHRF_R1_cDNA.extendedFrags.fastq.gz > Output_folder/CHRF_R1_cDNA.sai"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_converting_and_filtering_CHRF_R1_cDNA.out -q normal -n5 -w"done(141_aligning_CHRF_R1_cDNA_job)" -J 141_converting_and_filtering_CHRF_R1_cDNA_job -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"bwa samse /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Dependencies/Library_TRIMMED_15bp_Carried_Variants.fasta Output_folder/CHRF_R1_cDNA.sai Output_folder/CHRF_R1_cDNA.extendedFrags.fastq.gz > Output_folder/CHRF_R1_cDNA.sam"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_converting_and_filtering_CHRF_R1_cDNA.out -q normal -n5 -w"done(141_converting_and_filtering_CHRF_R1_cDNA_job)" -J 141_sam_to_bam_CHRF_R1_cDNA_job -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"samtools view -@5 -h -F 256 Output_folder/CHRF_R1_cDNA.sam|samtools sort -@5 -o Output_folder/CHRF_R1_cDNA.bam"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_converting_and_filtering_CHRF_R1_cDNA.out -q normal -n5 -w"done(141_sam_to_bam_CHRF_R1_cDNA_job)" -J 141_filter_and_rest_CHRF_R1_cDNA_job -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"bash /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Bash_subscript/268_MPRA_Sample_procesing_UMI_TOOLS_v2_sub_bash.sh Output_folder/ /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/ 141 20000 5 normal CHRF_R1_cDNA CHRF_R1_cDNA"


#########################################################################################################################################################################
#########################################################################################################################################################################
#########################################################################################################################################################################
#########################################################################################################################################################################




bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing//Lane_merger_CHRF_R2_gDNA_CHRF_gDNA_R1.out -q normal -n5 -J 141_CHRF_R2_gDNA_1 -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"cat /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_gDNA_S5_L001_R1_001.fastq.gz /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_gDNA_S5_L002_R1_001.fastq.gz  > /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_gDNA_R1.fastq.gz"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing//Lane_merger_CHRF_R2_gDNA_CHRF_gDNA_R2.out -q normal -n5 -J 141_CHRF_R2_gDNA_2 -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"cat /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_gDNA_S5_L001_R2_001.fastq.gz /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_gDNA_S5_L002_R2_001.fastq.gz  > /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_gDNA_R2.fastq.gz"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing//Lane_merger_CHRF_R2_gDNA_CHRF_gDNA_R3.out -q normal -n5 -J 141_CHRF_R2_gDNA_3 -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"cat /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_gDNA_S5_L001_R3_001.fastq.gz /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_gDNA_S5_L002_R3_001.fastq.gz  > /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_gDNA_R3.fastq.gz"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing//Lane_merger_CHRF_R2_gDNA_CHRF_gDNA_I1.out -q normal -n5 -J 141_CHRF_R2_gDNA_4 -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"cat /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_gDNA_S5_L001_I1_001.fastq.gz /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_gDNA_S5_L002_I1_001.fastq.gz  > /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_gDNA_I1.fastq.gz"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing//CHRF_R2_gDNA_Processing.out -q normal -n5 -w"done(141_CHRF_R2_gDNA_1) && done(141_CHRF_R2_gDNA_2) && done(141_CHRF_R2_gDNA_3) && done(141_CHRF_R2_gDNA_4)" -J 141_Trimmomatic_CHRF_R2_gDNA_job -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"java -jar /nfs/users/nfs_m/mt19/sOFTWARE/Trimmomatic-0.39/trimmomatic-0.39.jar SE /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_gDNA_R3.fastq.gz /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/Trimmed_CHRF_R2_gDNA_R3.fastq.gz HEADCROP:6"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_Extracting_UMI_R1_CHRF_R2_gDNA.out -M 20000 -w"done(141_Trimmomatic_CHRF_R2_gDNA_job)" -J 141_Extracting_UMI_R1_CHRF_R2_gDNA_job -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -n5 -q normal -- \
"/nfs/team151/software/umi-tools/bin/umi_tools extract --bc-pattern=NNNNNNNNNN --stdin=/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_gDNA_R2.fastq.gz --read2-in=/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_gDNA_R1.fastq.gz --stdout=Output_folder/CHRF_R2_gDNA_R1_Plus_UMIS.fastq.gz --read2-stdout"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_Extracting_UMI_R3_CHRF_R2_gDNA.out -M 20000 -w"done(141_Trimmomatic_CHRF_R2_gDNA_job)" -J 141_Extracting_UMI_R3_CHRF_R2_gDNA_job -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -n5 -q normal -- \
"/nfs/team151/software/umi-tools/bin/umi_tools extract --bc-pattern=NNNNNNNNNN --stdin=/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_gDNA_R2.fastq.gz --read2-in=/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/Trimmed_CHRF_R2_gDNA_R3.fastq.gz --stdout=Output_folder/CHRF_R2_gDNA_R3_Plus_UMIS.fastq.gz --read2-stdout"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_R1_R3_merging_CHRF_R2_gDNA.out -q normal -n5 -w"done(141_Extracting_UMI_R1_CHRF_R2_gDNA_job) && done(141_Extracting_UMI_R3_CHRF_R2_gDNA_job)" -J 141_R1_R3_merging_CHRF_R2_gDNA_job -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"flash2 -z --output-prefix=CHRF_R2_gDNA -d Output_folder/ -t 1 -r 11 -f 11 -s 2 -m 10 -x 0.122 Output_folder/CHRF_R2_gDNA_R1_Plus_UMIS.fastq.gz Output_folder/CHRF_R2_gDNA_R3_Plus_UMIS.fastq.gz"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_R1_R3_merging_CHRF_R2_gDNA.out -M 4000  -w"done(141_R1_R3_merging_CHRF_R2_gDNA_job)" -J 141_CLEANING_1_CHRF_R2_gDNA_job  -R"select[mem>=4000] rusage[mem=4000] span[hosts=1]" -n1 -q normal -- \
"rm /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_gDNA_R1.fastq.gz /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/Trimmed_CHRF_R2_gDNA_R3.fastq.gz /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_gDNA_R3.fastq.gz"


#########################################################################################################################################################################
#########################################################################################################################################################################
bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_aligning_CHRF_R2_gDNA.out -q normal -n5 -w"done(141_R1_R3_merging_CHRF_R2_gDNA_job)" -J 141_aligning_CHRF_R2_gDNA_job -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"bwa aln -l 10 -O 100 -E 100 /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Dependencies/Library_TRIMMED_15bp_Carried_Variants.fasta Output_folder/CHRF_R2_gDNA.extendedFrags.fastq.gz > Output_folder/CHRF_R2_gDNA.sai"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_converting_and_filtering_CHRF_R2_gDNA.out -q normal -n5 -w"done(141_aligning_CHRF_R2_gDNA_job)" -J 141_converting_and_filtering_CHRF_R2_gDNA_job -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"bwa samse /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Dependencies/Library_TRIMMED_15bp_Carried_Variants.fasta Output_folder/CHRF_R2_gDNA.sai Output_folder/CHRF_R2_gDNA.extendedFrags.fastq.gz > Output_folder/CHRF_R2_gDNA.sam"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_converting_and_filtering_CHRF_R2_gDNA.out -q normal -n5 -w"done(141_converting_and_filtering_CHRF_R2_gDNA_job)" -J 141_sam_to_bam_CHRF_R2_gDNA_job -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"samtools view -@5 -h -F 256 Output_folder/CHRF_R2_gDNA.sam|samtools sort -@5 -o Output_folder/CHRF_R2_gDNA.bam"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_converting_and_filtering_CHRF_R2_gDNA.out -q normal -n5 -w"done(141_sam_to_bam_CHRF_R2_gDNA_job)" -J 141_filter_and_rest_CHRF_R2_gDNA_job -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"bash /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Bash_subscript/268_MPRA_Sample_procesing_UMI_TOOLS_v2_sub_bash.sh Output_folder/ /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/ 141 20000 5 normal CHRF_R2_gDNA CHRF_R2_gDNA"


#########################################################################################################################################################################
#########################################################################################################################################################################
#########################################################################################################################################################################
#########################################################################################################################################################################




bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing//Lane_merger_CHRF_R2_cDNA_CHRF_cDNA_R1.out -q normal -n5 -J 141_CHRF_R2_cDNA_1 -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"cat /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_cDNA_S10_L001_R1_001.fastq.gz /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_cDNA_S10_L002_R1_001.fastq.gz  > /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_cDNA_R1.fastq.gz"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing//Lane_merger_CHRF_R2_cDNA_CHRF_cDNA_R2.out -q normal -n5 -J 141_CHRF_R2_cDNA_2 -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"cat /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_cDNA_S10_L001_R2_001.fastq.gz /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_cDNA_S10_L002_R2_001.fastq.gz  > /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_cDNA_R2.fastq.gz"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing//Lane_merger_CHRF_R2_cDNA_CHRF_cDNA_R3.out -q normal -n5 -J 141_CHRF_R2_cDNA_3 -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"cat /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_cDNA_S10_L001_R3_001.fastq.gz /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_cDNA_S10_L002_R3_001.fastq.gz  > /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_cDNA_R3.fastq.gz"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing//Lane_merger_CHRF_R2_cDNA_CHRF_cDNA_I1.out -q normal -n5 -J 141_CHRF_R2_cDNA_4 -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"cat /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_cDNA_S10_L001_I1_001.fastq.gz /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_cDNA_S10_L002_I1_001.fastq.gz  > /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_cDNA_I1.fastq.gz"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing//CHRF_R2_cDNA_Processing.out -q normal -n5 -w"done(141_CHRF_R2_cDNA_1) && done(141_CHRF_R2_cDNA_2) && done(141_CHRF_R2_cDNA_3) && done(141_CHRF_R2_cDNA_4)" -J 141_Trimmomatic_CHRF_R2_cDNA_job -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"java -jar /nfs/users/nfs_m/mt19/sOFTWARE/Trimmomatic-0.39/trimmomatic-0.39.jar SE /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_cDNA_R3.fastq.gz /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/Trimmed_CHRF_R2_cDNA_R3.fastq.gz HEADCROP:6"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_Extracting_UMI_R1_CHRF_R2_cDNA.out -M 20000 -w"done(141_Trimmomatic_CHRF_R2_cDNA_job)" -J 141_Extracting_UMI_R1_CHRF_R2_cDNA_job -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -n5 -q normal -- \
"/nfs/team151/software/umi-tools/bin/umi_tools extract --bc-pattern=NNNNNNNNNN --stdin=/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_cDNA_R2.fastq.gz --read2-in=/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_cDNA_R1.fastq.gz --stdout=Output_folder/CHRF_R2_cDNA_R1_Plus_UMIS.fastq.gz --read2-stdout"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_Extracting_UMI_R3_CHRF_R2_cDNA.out -M 20000 -w"done(141_Trimmomatic_CHRF_R2_cDNA_job)" -J 141_Extracting_UMI_R3_CHRF_R2_cDNA_job -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -n5 -q normal -- \
"/nfs/team151/software/umi-tools/bin/umi_tools extract --bc-pattern=NNNNNNNNNN --stdin=/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_cDNA_R2.fastq.gz --read2-in=/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/Trimmed_CHRF_R2_cDNA_R3.fastq.gz --stdout=Output_folder/CHRF_R2_cDNA_R3_Plus_UMIS.fastq.gz --read2-stdout"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_R1_R3_merging_CHRF_R2_cDNA.out -q normal -n5 -w"done(141_Extracting_UMI_R1_CHRF_R2_cDNA_job) && done(141_Extracting_UMI_R3_CHRF_R2_cDNA_job)" -J 141_R1_R3_merging_CHRF_R2_cDNA_job -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"flash2 -z --output-prefix=CHRF_R2_cDNA -d Output_folder/ -t 1 -r 11 -f 11 -s 2 -m 10 -x 0.122 Output_folder/CHRF_R2_cDNA_R1_Plus_UMIS.fastq.gz Output_folder/CHRF_R2_cDNA_R3_Plus_UMIS.fastq.gz"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_R1_R3_merging_CHRF_R2_cDNA.out -M 4000  -w"done(141_R1_R3_merging_CHRF_R2_cDNA_job)" -J 141_CLEANING_1_CHRF_R2_cDNA_job  -R"select[mem>=4000] rusage[mem=4000] span[hosts=1]" -n1 -q normal -- \
"rm /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_cDNA_R1.fastq.gz /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/Trimmed_CHRF_R2_cDNA_R3.fastq.gz /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/CHRF_R2_cDNA_R3.fastq.gz"


#########################################################################################################################################################################
#########################################################################################################################################################################
bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_aligning_CHRF_R2_cDNA.out -q normal -n5 -w"done(141_R1_R3_merging_CHRF_R2_cDNA_job)" -J 141_aligning_CHRF_R2_cDNA_job -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"bwa aln -l 10 -O 100 -E 100 /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Dependencies/Library_TRIMMED_15bp_Carried_Variants.fasta Output_folder/CHRF_R2_cDNA.extendedFrags.fastq.gz > Output_folder/CHRF_R2_cDNA.sai"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_converting_and_filtering_CHRF_R2_cDNA.out -q normal -n5 -w"done(141_aligning_CHRF_R2_cDNA_job)" -J 141_converting_and_filtering_CHRF_R2_cDNA_job -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"bwa samse /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Dependencies/Library_TRIMMED_15bp_Carried_Variants.fasta Output_folder/CHRF_R2_cDNA.sai Output_folder/CHRF_R2_cDNA.extendedFrags.fastq.gz > Output_folder/CHRF_R2_cDNA.sam"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_converting_and_filtering_CHRF_R2_cDNA.out -q normal -n5 -w"done(141_converting_and_filtering_CHRF_R2_cDNA_job)" -J 141_sam_to_bam_CHRF_R2_cDNA_job -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"samtools view -@5 -h -F 256 Output_folder/CHRF_R2_cDNA.sam|samtools sort -@5 -o Output_folder/CHRF_R2_cDNA.bam"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/outfile_converting_and_filtering_CHRF_R2_cDNA.out -q normal -n5 -w"done(141_sam_to_bam_CHRF_R2_cDNA_job)" -J 141_filter_and_rest_CHRF_R2_cDNA_job -M 20000 -R"select[mem>=20000] rusage[mem=20000] span[hosts=1]" -- \
"bash /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Bash_subscript/268_MPRA_Sample_procesing_UMI_TOOLS_v2_sub_bash.sh Output_folder/ /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Log_files_Sample_Processing/ 141 20000 5 normal CHRF_R2_cDNA CHRF_R2_cDNA"


#########################################################################################################################################################################
#########################################################################################################################################################################
#########################################################################################################################################################################
#########################################################################################################################################################################




