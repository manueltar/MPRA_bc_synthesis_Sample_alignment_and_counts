#!/usr/bin/env bash
 
 
MASTER_ROUTE=$1
Log_files_path=$2
batch=$3
mem=$4
pc=$5
queue=$6
master_prefix=$7
name_sample=$8

source /software/hgi/installs/anaconda3/etc/profile.d/conda.sh


conda_bamtools=$(echo "/nfs/team151/software/bamtools/")

conda deactivate

conda activate $conda_bamtools

output_dir=$(echo "$MASTER_ROUTE")

outfile_subscript=$(echo "$Log_files_path""outfile_subscript""_""$name_sample"".out")
touch $outfile_subscript
echo -n "" > $outfile_subscript

outfile_umi_tools_grouping=$(echo "$Log_files_path""outfile""_""umi_tools_grouping""_""$name_sample"".out")
touch $outfile_umi_tools_grouping
echo -n "" > $outfile_umi_tools_grouping

outfile_BASHING=$(echo "$Log_files_path""outfile""_""BASHING""_""$name_sample"".out")
touch $outfile_BASHING
echo -n "" > $outfile_BASHING

echo "-------------------------------->$batch-------------------->""$master_prefix""-------------------->""$name_sample"

    
# This step filters out all reads that have any changes comprared to the 11bp barcode in the reference to which they align
# https://sourceforge.net/p/bio-bwa/mailman/message/26514149/
# > NM      Number of nucleotide differences (i.e. edit distance to the 
# > reference sequence). Edit distance equals the sum of mismatches, gap 
# > opens and gap extensions. Mismatches or gaps in clipped sequences are 
# > not counted.

        type=$(echo "filtering_2")
        name_filtering_2=$(echo "$batch""_""$type""_""$name_sample""_job")

        output_bwa_bam_filtered_DEF=$(echo "$output_dir""$master_prefix""_filtered_DEF.bam")
        output_bwa_bam=$(echo "$output_dir""$master_prefix"".bam")

	conda activate $conda_bamtools
	
        bsub -G team151 -o $outfile_subscript -q normal -n$pc -J $name_filtering_2 -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
        "bamtools filter -tag NM:0 -in $output_bwa_bam -out $output_bwa_bam_filtered_DEF"


	# This step sorts the filtered bam
	
        type=$(echo "sorting")
        name_sorting=$(echo "$batch""_""$type""_""$name_sample""_job")

        output_bwa_bam_filtered_sorted=$(echo "$output_dir""$master_prefix""_filtered_sorted.bam")


        bsub -G team151 -o $outfile_subscript -q normal -n$pc -w"done($name_filtering_2)" -J $name_sorting -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
        "samtools sort -@$pc $output_bwa_bam_filtered_DEF -o $output_bwa_bam_filtered_sorted"


	# This step indexes the sorted and filtered bam

        type=$(echo "indexing")
        name_indexing=$(echo "$batch""_""$type""_""$name_sample""_job")


        bsub -G team151 -o $outfile_subscript -q normal -n$pc -w"done($name_sorting)" -J $name_indexing -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
	"samtools index $output_bwa_bam_filtered_sorted"


	# This step eliminates the sam and sai files for every sample
	
        type=$(echo "CLEANING_2")
        name_CLEANING=$(echo "$batch""_""$type""_""$name_sample""_job")
	sam_output=$(echo "$output_dir""$master_prefix"".sam")
	sai_output=$(echo "$output_dir""$master_prefix"".sai")

        bsub -G team151 -o $outfile_subscript -M 4000  -w"done($name_indexing)" -J $name_CLEANING  -R"select[mem>=4000] rusage[mem=4000] span[hosts=1]" -n1 -q $queue -- \
        "rm $sam_output $sai_output" #$output_bwa_bam_filtered $output_bwa_bam_filtered_DEF"

	# This step groups the UMIs

        type=$(echo "umi_tools_grouping")
        name_umi_tools_grouping=$(echo "$batch""_""$type""_""$name_sample""_job")



        output_umi_tools_group=$(echo "$output_dir""$master_prefix""_group.tsv")
        output_umi_tools_log=$(echo "$output_dir""$master_prefix""_group.log")

	type=$(echo "BASHING")



        bsub -G team151 -o $outfile_umi_tools_grouping -M $mem -J $name_umi_tools_grouping -w"done($name_indexing)" -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
        "/nfs/team151/software/umi-tools/bin/umi_tools group -I $output_bwa_bam_filtered_sorted --group-out=$output_umi_tools_group --log=$output_umi_tools_log"



	# This step takes the contig and the final_umi field of the group.tsv file and sorts them
	# the final output has the element with all the UMIs from all the reads mapped to it:

	# Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA     AAAAATTTTA
	# Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA     AAAAATTTTA
	# Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA     AACACCCTTA
	# Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA     ACGCCCCATA
	# Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA     ACGCCCCATA
	# Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA     AGAGATGACT
	# Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA     AGTTTAGCCT
	# Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA     CAAAGGTGGC
	# Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA     CACTAAGCCC
	# Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA     CAGGGTCAAT
	# Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA     CGAAAGTCAA
	
	type=$(echo "UMI_sorting")
        name_UMI_sorting=$(echo  "$batch""_""$type""_""$name_sample""_job")
        output_umi_summary=$(echo "$output_dir""$master_prefix""_UMI.tsv")


        bsub -G team151 -o $outfile_BASHING -M $mem  -w"done($name_umi_tools_grouping)" -J $name_UMI_sorting  -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
        "cut -f2,7 $output_umi_tools_group|sort --parallel=$pc -S 20G > $output_umi_summary"

	# This step collapses and counts the unique UMIS per contig

	# 2 Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA     AAAAATTTTA
	# 1 Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA     AACACCCTTA
	# 2 Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA     ACGCCCCATA
	# 1 Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA     AGAGATGACT

	
        type=$(echo "UMI_unique")
        name_UMI_unique=$(echo "$batch""_""$type""_""$name_sample""_job")
        output_umi_summary_unique=$(echo "$output_dir""$master_prefix""_UMI_unique.tsv")

        bsub -G team151 -o $outfile_BASHING -M $mem  -w"done($name_UMI_sorting)" -J $name_UMI_unique  -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
        "uniq -c $output_umi_summary > $output_umi_summary_unique"

	# This step extracts the contig and removes the counts of UMIS

	# Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA
	# Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA
	# Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA
	# Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA
	# Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA
	# Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA
	# Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA
	# Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA
	# Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA
	# Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA
	# Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA
	# Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA
	# Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA
	# Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA
	# Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA
	# Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA

	
        type=$(echo "Deduplicated_sorting")
        name_Deduplicated_sorting=$(echo  "$batch""_""$type""_""$name_sample""_job")
        output_umi_tools_deduplicated_sorted=$(echo "$output_dir""$master_prefix""_deduplicated_sorted.tsv")


        bsub -G team151 -o $outfile_BASHING -M $mem  -w"done($name_UMI_unique)" -J $name_Deduplicated_sorting  -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
        "cat $output_umi_summary_unique|tr -s ' '|cut -d ' ' -f 3|cut -f1|sort --parallel=$pc -S 20G > $output_umi_tools_deduplicated_sorted"

	# This step counts the contig appearances, the final counts

	# 16 Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA

        type=$(echo "Deduplicated_unique")
        name_Deduplicated_unique=$(echo "$batch""_""$type""_""$name_sample""_job")
        output_umi_tools_deduplicated_sorted_unique=$(echo "$output_dir""$master_prefix""_deduplicated_sorted_unique.tsv")

        bsub -G team151 -o $outfile_BASHING -M $mem  -w"done($name_Deduplicated_sorting)" -J $name_Deduplicated_unique  -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
        "uniq -c $output_umi_tools_deduplicated_sorted > $output_umi_tools_deduplicated_sorted_unique"

	# This step compresses the table of the results of the contig counts

        type=$(echo "Deduplicated_compress")
        name_Deduplicated_compress=$(echo "$batch""_""$type""_""$name_sample""_job")
        output_umi_tools_group_collapsed=$(echo "$output_umi_tools_group_unique"".gz")

        bsub -G team151 -o $outfile_BASHING -M $mem  -w"done($name_Deduplicated_unique)" -J $name_Deduplicated_compress  -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
        "gzip -f $output_umi_tools_deduplicated_sorted_unique"

	# This step compresses the table of the UMIS found in case I need to revise it

	# 2 Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA     AAAAATTTTA
	# 1 Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA     AACACCCTTA
	# 2 Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA     ACGCCCCATA
	# 1 Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA     AGAGATGACT
	# 1 Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA     AGTTTAGCCT
	# 1 Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA     CAAAGGTGGC
	# 1 Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA     CACTAAGCCC
	# 1 Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA     CAGGGTCAAT
	# 1 Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA     CGAAAGTCAA
	# 1 Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA     CTCCTGCTGG
	# 1 Element_100;A_TENTH;ALT;1;5_1093511_G_A;AAACCAGGCAA     CTCGAAGTGT


        type=$(echo "UMI_compress")
        name_UMI_compress=$(echo "$batch""_""$type""_""$name_sample""_job")
        output_umi_summary_collapsed=$(echo "$output_umi_summary_unique"".gz")

        bsub -G team151 -o $outfile_BASHING -M $mem  -w"done($name_Deduplicated_unique)" -J $name_UMI_compress  -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
        "gzip -f $output_umi_summary_unique"

	# This step cleans intermediate files

        type=$(echo "CLEANING_3")
        name_CLEANING=$(echo "$batch""_""$type""_""$name_sample""_job")

        bsub -G team151 -o $outfile_BASHING -M 4000  -w"done($name_Deduplicated_compress)" -J $name_CLEANING  -R"select[mem>=4000] rusage[mem=4000] span[hosts=1]" -n1 -q $queue -- \
        "rm $output_umi_tools_group $output_umi_summary $output_umi_tools_deduplicated_sorted"


