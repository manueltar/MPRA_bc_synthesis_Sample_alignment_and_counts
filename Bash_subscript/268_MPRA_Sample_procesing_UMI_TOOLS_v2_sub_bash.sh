#!/usr/bin/env bash
 
 
MASTER_ROUTE=$1
batch=$2
mem=$3
pc=$4
queue=$5
master_prefix=$6
name_sample=$7

source /software/hgi/installs/anaconda3/etc/profile.d/conda.sh


conda_bamtools=$(echo "/nfs/team151/software/bamtools/")

conda deactivate

conda activate $conda_bamtools

output_dir=$(echo "$MASTER_ROUTE")



echo "-------------------------------->$batch-------------------->""$master_prefix""-------------------->""$name_sample"


#------------> bamtools filter

        type=$(echo "filtering_2")
        name_filtering_2=$(echo "$batch""_""$type""_""$name_sample""_job")

        output_bwa_bam_filtered_DEF=$(echo "$output_dir""$master_prefix""_filtered_DEF.bam")

	conda activate $conda_bamtools
	
        bsub -G team151 -o $outfile_converting -q normal -n$pc -w"done($name_filtering_1)"  -J $name_filtering_2 -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
        "bamtools filter -tag NM:0 -in $output_bwa_bam_filtered -out $output_bwa_bam_filtered_DEF"



        type=$(echo "sorting")
        name_sorting=$(echo "$batch""_""$type""_""$name_sample""_job")

        output_bwa_bam_filtered_sorted=$(echo "$output_dir""$master_prefix""_filtered_sorted.bam")


        bsub -G team151 -o $outfile_converting -q normal -n$pc -w"done($name_filtering_2)" -J $name_sorting -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
        "samtools sort -@$pc $output_bwa_bam_filtered_DEF -o $output_bwa_bam_filtered_sorted"




        type=$(echo "indexing")
        name_indexing=$(echo "$batch""_""$type""_""$name_sample""_job")


        bsub -G team151 -o $outfile_converting -q normal -n$pc -w"done($name_sorting)" -J $name_indexing -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \




        type=$(echo "CLEANING_2")
        name_CLEANING=$(echo "$batch""_""$type""_""$name_sample""_job")

        bsub -G team151 -o $outfile_converting -M 4000  -w"done($name_indexing)" -J $name_CLEANING  -R"select[mem>=4000] rusage[mem=4000] span[hosts=1]" -n1 -q $queue -- \
        "rm $sam_output $sai_output" #$output_bwa_bam_filtered $output_bwa_bam_filtered_DEF"


        type=$(echo "umi_tools_grouping")
        outfile_umi_tools_grouping=$(echo "$output_dir""outfile""_""$type""_""$name_sample"".out")
        touch $outfile_umi_tools_grouping
        echo -n "" > $outfile_umi_tools_grouping
        name_umi_tools_grouping=$(echo "$batch""_""$type""_""$name_sample""_job")



        output_umi_tools_group=$(echo "$output_dir""$master_prefix""_group.tsv")
        output_umi_tools_log=$(echo "$output_dir""$master_prefix""_group.log")

	type=$(echo "BASHING")
        outfile_BASHING=$(echo "$output_dir""outfile""_""$type""_""$name_sample"".out")
        touch $outfile_BASHING
        echo -n "" > $outfile_BASHING


        bsub -G team151 -o $outfile_umi_tools_grouping -M $mem -J $name_umi_tools_grouping -w"done($name_indexing)" -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
        "umi_tools group -I $output_bwa_bam_filtered_sorted --group-out=$output_umi_tools_group --log=$output_umi_tools_log"

	type=$(echo "UMI_sorting")
        name_UMI_sorting=$(echo  "$batch""_""$type""_""$name_sample""_job")
        output_umi_summary=$(echo "$output_dir""$master_prefix""_UMI.tsv")


        bsub -G team151 -o $outfile_BASHING -M $mem  -w"done($name_umi_tools_grouping)" -J $name_UMI_sorting  -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
        "cut -f2,7 $output_umi_tools_group|sort --parallel=$pc -S 20G > $output_umi_summary"


        type=$(echo "UMI_unique")
        name_UMI_unique=$(echo "$batch""_""$type""_""$name_sample""_job")
        output_umi_summary_unique=$(echo "$output_dir""$master_prefix""_UMI_unique.tsv")

        bsub -G team151 -o $outfile_BASHING -M $mem  -w"done($name_UMI_sorting)" -J $name_UMI_unique  -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
        "uniq -c $output_umi_summary > $output_umi_summary_unique"


        type=$(echo "Deduplicated_sorting")
        name_Deduplicated_sorting=$(echo  "$batch""_""$type""_""$name_sample""_job")
        output_umi_tools_deduplicated_sorted=$(echo "$output_dir""$master_prefix""_deduplicated_sorted.tsv")


        bsub -G team151 -o $outfile_BASHING -M $mem  -w"done($name_UMI_unique)" -J $name_Deduplicated_sorting  -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
        "cat $output_umi_summary_unique|tr -s ' '|cut -d ' ' -f 3|cut -f1|sort --parallel=$pc -S 20G > $output_umi_tools_deduplicated_sorted"


        type=$(echo "Deduplicated_unique")
        name_Deduplicated_unique=$(echo "$batch""_""$type""_""$name_sample""_job")
        output_umi_tools_deduplicated_sorted_unique=$(echo "$output_dir""$master_prefix""_deduplicated_sorted_unique.tsv")

        bsub -G team151 -o $outfile_BASHING -M $mem  -w"done($name_Deduplicated_sorting)" -J $name_Deduplicated_unique  -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
        "uniq -c $output_umi_tools_deduplicated_sorted > $output_umi_tools_deduplicated_sorted_unique"


        type=$(echo "Deduplicated_compress")
        name_Deduplicated_compress=$(echo "$batch""_""$type""_""$name_sample""_job")
        output_umi_tools_group_collapsed=$(echo "$output_umi_tools_group_unique"".gz")

        bsub -G team151 -o $outfile_BASHING -M $mem  -w"done($name_Deduplicated_unique)" -J $name_Deduplicated_compress  -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
        "gzip -f $output_umi_tools_deduplicated_sorted_unique"



        type=$(echo "UMI_compress")
        name_UMI_compress=$(echo "$batch""_""$type""_""$name_sample""_job")
        output_umi_summary_collapsed=$(echo "$output_umi_summary_unique"".gz")

        bsub -G team151 -o $outfile_BASHING -M $mem  -w"done($name_Deduplicated_unique)" -J $name_UMI_compress  -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
        "gzip -f $output_umi_summary_unique"



        type=$(echo "CLEANING_3")
        name_CLEANING=$(echo "$batch""_""$type""_""$name_sample""_job")

        bsub -G team151 -o $outfile_BASHING -M 4000  -w"done($name_Deduplicated_compress)" -J $name_CLEANING  -R"select[mem>=4000] rusage[mem=4000] span[hosts=1]" -n1 -q $queue -- \
        "rm $output_umi_tools_group $output_umi_summary $output_umi_tools_deduplicated_sorted"


