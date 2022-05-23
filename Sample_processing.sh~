#!/bin/bash>

 
#### Rscript

input_file=$1
MASTER_ROUTE=$2
batch=$3
path_bash_scripts=$4

path_Bash_sub_scripts=$(echo "$path_bash_scripts""Bash_sub_scripts/")
path_Dependencies=$(echo "$path_bash_scripts""Dependencies/")

path_to_fastq_file=$5

mem=$6
pc=$7
queue=$8


output_dir=$(echo "$MASTER_ROUTE")
type_for_R=$(echo "100_percent")
declare -a arr

 
output_dir=$MASTER_ROUTE


end_of_file=0
counter=0
array=()
array_Merge=()                  #

cd $output_dir

while [[ $end_of_file == 0 ]]
do
  read -r line
  end_of_file=$?

#echo "LINE:$line"

    if [[ "$line" =~ [^[:space:]] ]]; then

    ((counter++))
    a=($(echo "$line" | tr "\t" '\n'))
    Lane_1=${a[0]}
    Lane_2=${a[1]}
    Merge=${a[2]}
  
    echo "PROCESSED: $Lane_1 $Lane_2 > $Merge COUNTER $counter"

    

    name_sample=$(echo "$Merge"|sed -r 's/_[RI][0-9]+\.fastq\.gz//g')


    
    name_sample=$(echo "$name_sample"|sed -r "s|$path_to_fastq_file||g")
    name_job=$(echo "$batch""_""$name_sample""_""$counter")

    
    Replicate=$(echo "$name_sample"|sed 's/_gDNA//g')
    Replicate=$(echo "$Replicate"|sed 's/_cDNA//g')


    Type=$(echo "$name_sample"|sed "s|$Replicate||g")
    Type=$(echo "$Type"|sed -r 's/^_//g')
    
    echo "------------------------>$name_sample $Type $Replicate"    
    
    echo "------------------------>name_job: $name_job"

    echo "------------------------>$Lane_1 $Lane_2 $Merge"

#exit

    stderrfile2=$(echo "$output_dir""/""Lane_merger_"$name_sample"_"$Type".out")
 
    touch $stderrfile2

    echo -n "" > $stderrfile2

     

#    exit    
    bsub -G team151 -o $stderrfile2 -q $queue -n$pc -J $name_job -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
   "cat $Lane_1 $Lane_2  > $Merge"



    wait_done=$(echo "done(""$name_job"") &&")
    array+=($wait_done)
    array_Merge+=($Merge)
    echo "---------------------------------------------------------------------------$Merge------------------------------------------>$counter\n"


       if [ $counter == "4" ]; then

        printf -- "%s\n" "${array_Merge[@]}"
        string=$(printf " %s" "${array[@]}")
#       echo "string:->$string<-"
        string2=$(echo "$string"|sed -r 's/^\s+//g')
        string3=$(echo "$string2"|sed -r 's/\s+&&$//g')
#       echo "string3:->$string3<-"
#       string4=$(echo "-w\"""$string3""\"")
#        echo "string4:->$string4<-"

        R1=${array_Merge[0]}
        UMI=${array_Merge[1]}
        R3=${array_Merge[2]}
        master_prefix=$name_sample

        echo "----------$name_sample--------->$R1 $UMI $R3 $master_prefix<-"



	#	---------------------> Dependency Trimmomatic-0.39
	
        outfile_Trimmomatic=$(echo "$output_dir""/""$name_sample""_Processing.out")
        touch $outfile_Trimmomatic
        echo -n "" > $outfile_Trimmomatic

	
        name_trimmomatic=$(echo "$batch""_""Trimmomatic""_""$name_sample""_job")
        output_trimmomatic_r3=$(echo "$output_dir""$Type""_""Trimmed_""$R3")
	output_trimmomatic_r3=$(echo "$output_trimmomatic_r3"|sed -r "s|$path_to_fastq_file||g")
	
        echo "->$name_trimmomatic $output_trimmomatic_r3<-"

        bsub -G team151 -o $outfile_Trimmomatic -q $queue -n$pc -w"$string3" -J $name_trimmomatic -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
        "java -jar /nfs/users/nfs_m/mt19/sOFTWARE/Trimmomatic-0.39/trimmomatic-0.39.jar SE $R3 $output_trimmomatic_r3 HEADCROP:6"


	#       ---------------------> umi_tools

	type=$(echo "Extracting_UMI_R1")
	outfile_Extracting_UMI_R1=$(echo "$output_dir""outfile""_""$type""_""$name_sample"".out")
	touch $outfile_Extracting_UMI_R1
	echo -n "" > $outfile_Extracting_UMI_R1
	name_Extracting_UMI_R1=$(echo "$batch""_""$type""_""$name_sample""_job")
	
	output_processed_R1=$(echo "$output_dir""$master_prefix""_R1_Plus_UMIS"".fastq.gz")


	bsub -G team151 -o $outfile_Extracting_UMI_R1 -M $mem -w"done($name_trimmomatic)" -J $name_Extracting_UMI_R1 -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
	"/nfs/team151/software/umi-tools/bin/umi_tools extract --bc-pattern=NNNNNNNNNN --stdin=$UMI --read2-in=$R1 --stdout=$output_processed_R1 --read2-stdout"


	type=$(echo "Extracting_UMI_R3")
	outfile_Extracting_UMI_R3=$(echo "$output_dir""outfile""_""$type""_""$name_sample"".out")
	touch $outfile_Extracting_UMI_R3
	echo -n "" > $outfile_Extracting_UMI_R3
	name_Extracting_UMI_R3=$(echo "$batch""_""$type""_""$name_sample""_job")

	
	output_processed_R3=$(echo "$output_dir""$master_prefix""_R3_Plus_UMIS"".fastq.gz")


	bsub -G team151 -o $outfile_Extracting_UMI_R3 -M $mem -w"done($name_trimmomatic)" -J $name_Extracting_UMI_R3 -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
	"/nfs/team151/software/umi-tools/bin/umi_tools extract --bc-pattern=NNNNNNNNNN --stdin=$UMI --read2-in=$output_trimmomatic_r3 --stdout=$output_processed_R3 --read2-stdout"


        #       ---------------------> flash


	 type=$(echo "R1_R3_merging")
        outfile_merging=$(echo "$output_dir""outfile""_""$type""_""$name_sample"".out")
        touch $outfile_merging
        echo -n "" > $outfile_merging
         name_merging=$(echo "$batch""_""$type""_""$name_sample""_job")
	 flash2_output=$(echo "$output_dir""$master_prefix"".extendedFrags.fastq.gz")
	
	
	bsub -G team151 -o $outfile_merging -q $queue -n$pc -w"done($name_Extracting_UMI_R1) && done($name_Extracting_UMI_R3)" -J $name_merging -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
        "flash2 -z --output-prefix=$master_prefix -t 1 -r 11 -f 11 -s 2 -m 10 -x 0.122 $output_processed_R1 $output_processed_R3"
	

        type=$(echo "CLEANING_1")
        name_CLEANING=$(echo "$batch""_""$type""_""$name_sample""_job")

        bsub -G team151 -o $outfile_merging -M 4000  -w"done($name_merging)" -J $name_CLEANING  -R"select[mem>=4000] rusage[mem=4000] span[hosts=1]" -n1 -q $queue -- \
        "rm $R1 $R2 $output_trimmomatic_r3 $R3"


#######################################################################################################################################################

	#	----> bwa
	
	type=$(echo "aligning")
        outfile_aligning=$(echo "$output_dir""outfile""_""$type""_""$name_sample"".out")
        touch $outfile_aligning
        echo -n "" > $outfile_aligning
        name_aligning=$(echo "$batch""_""$type""_""$name_sample""_job")
	REFERENCE=$(echo "$path_Dependencies""Library_TRIMMED_15bp_Carried_Variants.fasta")

	sai_output=$(echo "$output_dir""$master_prefix"".sai")	
	
	bsub -G team151 -o $outfile_aligning -q $queue -n$pc -w"done($name_merging)" -J $name_aligning -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
        "bwa aln -l 10 -O 100 -E 100 $REFERENCE $flash2_output > $sai_output"

	 
	type=$(echo "converting_and_filtering")
        outfile_converting=$(echo "$output_dir""outfile""_""$type""_""$name_sample"".out")
        touch $outfile_converting
        echo -n "" > $outfile_converting
        name_converting=$(echo "$batch""_""$type""_""$name_sample""_job")
        REFERENCE=$(echo "$path_Dependencies""Library_TRIMMED_15bp_Carried_Variants.fasta")

	flash2_output=$(echo "$output_dir""$master_prefix"".extendedFrags.fastq.gz")
	sam_output=$(echo "$output_dir""$master_prefix"".sam")

	
	bsub -G team151 -o $outfile_converting -q $queue -n$pc -w"done($name_aligning)" -J $name_converting -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
        "bwa samse $REFERENCE $sai_output $flash2_output > $sam_output"


	  type=$(echo "clean_2""_""$master_prefix")
	  outfile_clean_2=$(echo "$output_dir""outfile_""$type"".out")
	  touch $outfile_clean_2
	  echo -n "" > $outfile_clean_2
	  name_clean_2=$(echo "$type""_job")


	  bsub -G team151 -o $outfile_clean_2 -q $queue -n$pc -w"done($name_converting)" -J $name_clean_2 -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
       "rm $sai_output"



	  type=$(echo "sam_to_bam_and_filter_first_alignment""_""$master_prefix")
	  outfile_sam_to_bam_and_filter_first_alignment=$(echo "$output_dir""outfile_""$type"".out")
	  touch $outfile_sam_to_bam_and_filter_first_alignment
	  echo -n "" > $outfile_sam_to_bam_and_filter_first_alignment
	  name_sam_to_bam_and_filter_first_alignment=$(echo "$type""_job")


	  sam_output=$(echo "$output_dir""$master_prefix"".sam")
	  output_bwa_bam=$(echo "$output_dir""$master_prefix"".bam")


	  bsub -G team151 -o $outfile_sam_to_bam_and_filter_first_alignment -q $queue -n$pc -w"done($name_converting)" -J $name_sam_to_bam_and_filter_first_alignment -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
	   "samtools view -@$pc -h -F 20,4,256 $sam_output|samtools sort -@$pc -o $output_bwa_bam"

	   ####### Get unalagined reads

    type=$(echo "unaligned_reads""_""$master_prefix")
    outfile_unaligned_reads=$(echo "$output_dir""outfile_""$type"".out")
    touch $outfile_unaligned_reads
    echo -n "" > $outfile_unaligned_reads
    name_unaligned_reads=$(echo "$type""_job")


    sam_output=$(echo "$output_dir""$master_prefix"".sam")
    unaligned_bwa_bam=$(echo "$output_dir""Unaligned_reads_""$master_prefix"".bam")


    bsub -G team151 -o $outfile_unaligned_reads -q $queue -n$pc -w"done($name_converting)" -J $name_unaligned_reads -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
    "samtools view -@$pc -h -f 20,4 -F 256 $sam_output|samtools sort -@$pc -o $unaligned_bwa_bam"

    ####### clean_2_5



    type=$(echo "clean_2_5""_""$master_prefix")
    name_clean_2=$(echo "$type""_job")


    bsub -G team151 -o $outfile_clean_2 -q $queue -n$pc -w"done($name_unaligned_reads)" -J $name_clean_2 -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
    "rm $sam_output"

	
##### indexing_1 #####


    type=$(echo "indexing_1""_""$master_prefix")
    outfile_indexing_1=$(echo "$output_dir""outfile_""$type"".out")
    touch $outfile_indexing_1
    echo -n "" > $outfile_indexing_1
    name_indexing_1=$(echo "$type""_job")


    sam_output=$(echo "$output_dir""$master_prefix"".sam")
    output_bwa_bam=$(echo "$output_dir""$master_prefix"".bam")


    bsub -G team151 -o $outfile_indexing_1 -q $queue -n$pc -w"done($name_sam_to_bam_and_filter_first_alignment)" -J $name_indexing_1 -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
	 "samtools index $output_bwa_bam"

      ################################################################## bamtools_filter ####################################################################################################

    source /software/hgi/installs/anaconda3/etc/profile.d/conda.sh


    conda_bamtools=$(echo "/nfs/team151/software/bamtools/")

    conda deactivate

    conda activate $conda_bamtools


    type=$(echo "bamtools_filter""_""$master_prefix")
    outfile_bamtools_filter=$(echo "$output_dir""outfile_""$type"".out")
    touch $outfile_bamtools_filter
    echo -n "" > $outfile_bamtools_filter
    name_bamtools_filter=$(echo "$type""_job")


    output_bwa_bam=$(echo "$output_dir""$master_prefix"".bam")
    output_bwa_bam_filtered=$(echo "$output_dir""$master_prefix""_filtered"".bam")


     bsub -G team151 -o $outfile_bamtools_filter -q $queue -n$pc -w"done($name_indexing_1)" -J $name_bamtools_filter -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
	  "bamtools filter -tag NM:0 -in $output_bwa_bam -out $output_bwa_bam_filtered"

    type=$(echo "indexing_2""_""$master_prefix")
    outfile_indexing_2=$(echo "$output_dir""outfile_""$type"".out")
    touch $outfile_indexing_2
    echo -n "" > $outfile_indexing_2
    name_indexing_2=$(echo "$type""_job")

    output_bwa_bam_filtered=$(echo "$output_dir""$master_prefix""_filtered"".bam")


    bsub -G team151 -o $outfile_indexing_2 -q $queue -n$pc -w"done($name_bamtools_filter)" -J $name_indexing_2 -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
	 "samtools index $output_bwa_bam_filtered"

    ########################################## umi_tools ##########################################
    ###################################################################################
    ##################################################################################

        type=$(echo "umi_tools_grouping""_""$master_prefix")
    outfile_umi_tools_grouping=$(echo "$output_dir""outfile""_""$type"".out")
    touch $outfile_umi_tools_grouping
    echo -n "" > $outfile_umi_tools_grouping
    name_umi_tools_grouping=$(echo "$type""_job")


    output_bwa_bam_filtered=$(echo "$output_dir""$master_prefix""_filtered"".bam")
    output_umi_tools_group=$(echo "$output_dir""$master_prefix""_group.tsv")
    output_umi_tools_log=$(echo "$output_dir""$master_prefix""_group.log")
    output_bwa_bam=$(echo "$output_dir""$master_prefix"".bam")

    bsub -G team151 -o $outfile_umi_tools_grouping -M $mem -J $name_umi_tools_grouping -w"done($name_indexing_2)" -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
    "/nfs/team151/software/umi-tools/bin/umi_tools group -I $output_bwa_bam_filtered --group-out=$output_umi_tools_group --log=$output_umi_tools_log"


    type=$(echo "bash_cut_and_sort""_""$master_prefix")
    outfile_bash_cut_and_sort=$(echo "$output_dir""outfile""_""$type"".out")
    touch $outfile_bash_cut_and_sort
    echo -n "" > $outfile_bash_cut_and_sort
    name_bash_cut_and_sort=$(echo "$type""_job")


    output_umi_tools_group=$(echo "$output_dir""$master_prefix""_group.tsv")
    output_umi_summary=$(echo "$output_dir""$master_prefix""_UMI.tsv")
    max_memory_G=$(expr $mem / 1000)
    max_memory_G=$(echo "$max_memory_G""G")

    echo "max_memory_G_1: $max_memory_G"


    bsub -G team151 -o $outfile_bash_cut_and_sort -M $mem -J $name_bash_cut_and_sort -w"done($name_umi_tools_grouping)" -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
    "cut -f2,7 $output_umi_tools_group|sort --parallel=$pc -S $max_memory_G > $output_umi_summary"


    type=$(echo "bash_unique""_""$master_prefix")
    outfile_bash_unique=$(echo "$output_dir""outfile""_""$type"".out")
    touch $outfile_bash_unique
    echo -n "" > $outfile_bash_unique
    name_bash_unique=$(echo "$type""_job")


    output_umi_summary=$(echo "$output_dir""$master_prefix""_UMI.tsv")
    output_umi_summary_unique=$(echo "$output_dir""$master_prefix""_UMI_unique.tsv")


    bsub -G team151 -o $outfile_bash_unique -M $mem -J $name_bash_unique -w"done($name_bash_cut_and_sort)" -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
	 "uniq -c $output_umi_summary > $output_umi_summary_unique"

    type=$(echo "deduplicate""_""$master_prefix")
    outfile_deduplicate=$(echo "$output_dir""outfile""_""$type"".out")
    touch $outfile_deduplicate
    echo -n "" > $outfile_deduplicate
    name_deduplicate=$(echo "$type""_job")

    output_umi_summary_unique=$(echo "$output_dir""$master_prefix""_UMI_unique.tsv")

    output_umi_tools_deduplicated_sorted=$(echo "$output_dir""$master_prefix""_deduplicated_sorted.tsv")

    max_memory_G=$(expr $mem / 1000)
    max_memory_G=$(echo "$max_memory_G""G")

    echo "max_memory_G_2: $max_memory_G"


    bsub -G team151 -o $outfile_deduplicate -M $mem -J $name_deduplicate -w"done($name_bash_unique)" -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
	 "cat $output_umi_summary_unique|tr -s ' '|cut -d ' ' -f 3|cut -f1|sort --parallel=$pc -S $max_memory_G > $output_umi_tools_deduplicated_sorted"

    type=$(echo "bash_unique_COMPRESS""_""$master_prefix")
    name_bash_unique_COMPRESS=$(echo "$type""_job")

    bsub -G team151 -o $outfile_bash_unique -q $queue -M $mem -J $name_bash_unique_COMPRESS -w"done($name_deduplicate)" -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -- \
	 "gzip -f $output_umi_summary_unique"

    type=$(echo "bash_unique_DEDUPLICATED""_""$master_prefix")
    outfile_bash_unique_DEDUPLICATED=$(echo "$output_dir""outfile""_""$type"".out")
    touch $outfile_bash_unique_DEDUPLICATED
    echo -n "" > $outfile_bash_unique_DEDUPLICATED
    name_bash_unique_DEDUPLICATED=$(echo "$type""_job")

    output_umi_tools_deduplicated_sorted=$(echo "$output_dir""$master_prefix""_deduplicated_sorted.tsv")
    output_umi_tools_deduplicated_sorted_unique=$(echo "$output_dir""$master_prefix""_deduplicated_sorted_unique.tsv")


    bsub -G team151 -o $outfile_bash_unique_DEDUPLICATED -M $mem -J $name_bash_unique_DEDUPLICATED -w"done($name_deduplicate)" -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
    "uniq -c $output_umi_tools_deduplicated_sorted > $output_umi_tools_deduplicated_sorted_unique"

    type=$(echo "bash_unique_DEDUPLICATED_COMPRESS""_""$master_prefix")
    name_bash_unique_DEDUPLICATED_COMPRESS=$(echo "$type""_job")

    bsub -G team151 -o $outfile_bash_unique_DEDUPLICATED -M $mem -J $name_bash_unique_DEDUPLICATED_COMPRESS -w"done($name_bash_unique_DEDUPLICATED)" -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
    "gzip -f $output_umi_tools_deduplicated_sorted_unique"



    type=$(echo "clean_3""_""$master_prefix")
    outfile_clean_3=$(echo "$output_dir""outfile_""$type"".out")
    touch $outfile_clean_3
    echo -n "" > $outfile_clean_3
    name_clean_3=$(echo "$type""_job")


    bsub -G team151 -o $outfile_clean_3 -q $queue -n$pc -w"done($name_bash_unique_DEDUPLICATED_COMPRESS) && done($name_bash_unique_COMPRESS)" -J $name_clean_3 -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
    "rm $output_umi_tools_deduplicated_sorted $output_umi_summary"



	echo "------------------------------------------->$counter\n"
	counter=0
        array=()
        array_Merge=()


       fi



  fi #no spaces
done < "$input_file"

















