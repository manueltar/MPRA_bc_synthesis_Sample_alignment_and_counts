#!/bin/bash>

 
#### Rscript

sample_file=$1

echo "$sample_file"

path_to_leave_output_files=$2
batch=$3
path_to_original_fastq_files=$4
path_to_scripts_folder=$5
mem=$6
pc=$7
queue=$8


####

dependencies_folder=$(echo "$path_to_scripts_folder""Dependencies""/")
bash_subscript_folder=$(echo "$path_to_scripts_folder""Bash_subscript""/")
Rscripts_folder=$(echo "$path_to_scripts_folder""Rscripts""/")
Log_files_path=$(echo "$path_to_scripts_folder""Log_files_Sample_Processing""/")

rm -rf $Log_files_path

mkdir -p $Log_files_path

output=$(echo "$path_to_scripts_folder""Sample_Processing_printed_script.sh")

touch $output
echo -n "" > $output

echo "#!/bin/bash"  >> $output




output_dir=$(echo "$path_to_leave_output_files")

declare -a arr



end_of_file=0
counter=0
array=()
array_Merge=()                  #

while [[ $end_of_file == 0 ]]
do
  read -r line
  end_of_file=$?

echo "LINE:$line"

    if [[ "$line" =~ [^[:space:]] ]]; then

    ((counter++))
    a=($(echo "$line" | tr "\t" '\n'))
    Lane_1=${a[0]}
    Lane_2=${a[1]}
    Merge=${a[2]}

  
#    echo "PROCESSED: $Lane_1 $Lane_2 > $Merge COUNTER $counter"


    name_sample=$(echo "$Merge"|sed -r 's/_[RI][0-9]+\.fastq\.gz//g')
    name_job=$(echo "$batch""_""$name_sample""_""$counter")
    Type=$(echo "$name_sample"|sed -r 's/R[0-9]+_//g')
    Replicate=$(echo "$name_sample"|sed -r 's/_.+//g')
    read_file=$(echo "$Merge"|sed -r 's/\.fastq\.gz//g')
    read_file=$(echo "$read_file"|sed -r 's/^[^_]+_[^_]+_[^_]+_//g')


    
    echo "------------------------>name_job: $name_job"
    echo "------------------------>$name_sample $Type $Replicate"
    echo "------------------------>$Lane_1 $Lane_2 $Merge"
    echo "------------------------>COUNTER $counter"
    echo "------------------------>read_file2 $read_file"

    

    Log_file=$(echo "$Log_files_path""/""Lane_merger_"$name_sample"_"$Type"""_""$read_file"".out")

    touch $Log_file

    echo -n "" > $Log_file

#    echo "------------------------>COUNTER $Log_file"

    # Merge files for the same sample and read (R1-R4) that are split in 2 lanes

    Lane_1=$(echo "$path_to_original_fastq_files""$Lane_1")    
    Lane_2=$(echo "$path_to_original_fastq_files""$Lane_2")
    Merge=$(echo "$path_to_original_fastq_files""$Merge")
    
    echo "bsub -G team151 -o $Log_file -q normal -n$pc -J $name_job -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\" >> $output
    echo "\"cat $Lane_1 $Lane_2  > $Merge\"" >> $output

    echo -e "\n"  >> $output

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



        outfile_Trimmomatic=$(echo "$Log_files_path""/""$name_sample""_Processing.out")
        touch $outfile_Trimmomatic
        echo -n "" > $outfile_Trimmomatic

        name_trimmomatic=$(echo "$batch""_""Trimmomatic""_""$name_sample""_job")
	output_trimmomatic_r3=$(echo "$path_to_original_fastq_files""Trimmed""_""$master_prefix""_""R3"".fastq.gz")

	echo "$output_trimmomatic_r3"

	echo "->$name_trimmomatic $output_trimmomatic_r3<-"


	# This step removes the 6 bases of the BamHI restriction site from the read2 of the sequencing (the R3 file in our series). This is not necessary for the 15 first bases of R1 because it was run in Dark cycling
	
        echo "bsub -G team151 -o $outfile_Trimmomatic -q normal -n$pc -w\"$string3\" -J $name_trimmomatic -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\" >> $output
#        echo "bsub -G team151 -o $outfile_Trimmomatic -q normal -n$pc -J $name_trimmomatic -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\" >> $output
        echo "\"java -jar /nfs/users/nfs_m/mt19/sOFTWARE/Trimmomatic-0.39/trimmomatic-0.39.jar SE $R3 $output_trimmomatic_r3 HEADCROP:6\"" >> $output

	echo -e "\n"  >> $output

	#### This step extracts the UMI from the UMI read file (R2) and adds it to R1 file

	type=$(echo "Extracting_UMI_R1")
	outfile_Extracting_UMI_R1=$(echo "$Log_files_path""outfile""_""$type""_""$name_sample"".out")
	touch $outfile_Extracting_UMI_R1
	echo -n "" > $outfile_Extracting_UMI_R1
	name_Extracting_UMI_R1=$(echo "$batch""_""$type""_""$name_sample""_job")
	
	output_processed_R1=$(echo "$output_dir""$master_prefix""_R1_Plus_UMIS"".fastq.gz")


	echo "bsub -G team151 -o $outfile_Extracting_UMI_R1 -M $mem -w\"done($name_trimmomatic)\" -J $name_Extracting_UMI_R1 -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
	echo "\"/nfs/team151/software/umi-tools/bin/umi_tools extract --bc-pattern=NNNNNNNNNN --stdin=$UMI --read2-in=$R1 --stdout=$output_processed_R1 --read2-stdout\"" >> $output

	#### This step extracts the UMI from the UMI read file (R2) and adds it to the processed R3 file
	
	echo -e "\n"  >> $output
	
	type=$(echo "Extracting_UMI_R3")
	outfile_Extracting_UMI_R3=$(echo "$Log_files_path""outfile""_""$type""_""$name_sample"".out")
	touch $outfile_Extracting_UMI_R3
	echo -n "" > $outfile_Extracting_UMI_R3
	name_Extracting_UMI_R3=$(echo "$batch""_""$type""_""$name_sample""_job")

	
	output_processed_R3=$(echo "$output_dir""$master_prefix""_R3_Plus_UMIS"".fastq.gz")


	echo "bsub -G team151 -o $outfile_Extracting_UMI_R3 -M $mem -w\"done($name_trimmomatic)\" -J $name_Extracting_UMI_R3 -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
	echo "\"/nfs/team151/software/umi-tools/bin/umi_tools extract --bc-pattern=NNNNNNNNNN --stdin=$UMI --read2-in=$output_trimmomatic_r3 --stdout=$output_processed_R3 --read2-stdout\"" >> $output


	echo -e "\n"  >> $output

	# This step merges the UMI added processed_R1 and the UMI added processed_R3, both of them should contain the 11 bp barcode of the MPRA each one in one orientation

	 type=$(echo "R1_R3_merging")
        outfile_merging=$(echo "$Log_files_path""outfile""_""$type""_""$name_sample"".out")
        touch $outfile_merging
        echo -n "" > $outfile_merging
         name_merging=$(echo "$batch""_""$type""_""$name_sample""_job")
	 flash2_output=$(echo "$output_dir""$master_prefix"".extendedFrags.fastq.gz")
	
	
	echo "bsub -G team151 -o $outfile_merging -q normal -n$pc -w\"done($name_Extracting_UMI_R1) && done($name_Extracting_UMI_R3)\" -J $name_merging -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\"  >> $output
        echo "\"flash2 -z --output-prefix=$master_prefix -d $output_dir -t 1 -r 11 -f 11 -s 2 -m 10 -x 0.122 $output_processed_R1 $output_processed_R3\"" >> $output
	
	echo -e "\n"  >> $output

#	exit
	# This step removes the first set of intermediate files
	
        type=$(echo "CLEANING_1")
        name_CLEANING=$(echo "$batch""_""$type""_""$name_sample""_job")

        echo "bsub -G team151 -o $outfile_merging -M 4000  -w\"done($name_merging)\" -J $name_CLEANING  -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q $queue -- \\" >> $output
        echo "\"rm $R1 $output_trimmomatic_r3 $R3\"" >> $output

	echo -e "\n"  >> $output


	echo "#########################################################################################################################################################################"  >> $output
	echo "#########################################################################################################################################################################"  >> $output

	# This is the aligning step with bwa aln to match the 11 bp barcode to the respective regulatory sequence

	type=$(echo "aligning")
        outfile_aligning=$(echo "$Log_files_path""outfile""_""$type""_""$name_sample"".out")
        touch $outfile_aligning
        echo -n "" > $outfile_aligning
        name_aligning=$(echo "$batch""_""$type""_""$name_sample""_job")
	REFERENCE=$(echo "$dependencies_folder""Library_TRIMMED_15bp_Carried_Variants.fasta")

	sai_output=$(echo "$output_dir""$master_prefix"".sai")	
	
	echo "bsub -G team151 -o $outfile_aligning -q normal -n$pc -w\"done($name_merging)\" -J $name_aligning -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\"  >> $output
#	echo "bsub -G team151 -o $outfile_aligning -q normal -n$pc -J $name_aligning -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\"  >> $output
        echo "\"bwa aln -l 10 -O 100 -E 100 $REFERENCE $flash2_output > $sai_output\"" >> $output

        echo -e "\n"  >> $output

	# This step converts the sai format to sam
	
	type=$(echo "converting_and_filtering")
        outfile_converting=$(echo "$Log_files_path""outfile""_""$type""_""$name_sample"".out")
        touch $outfile_converting
        echo -n "" > $outfile_converting
        name_converting=$(echo "$batch""_""$type""_""$name_sample""_job")
	REFERENCE=$(echo $dependencies_folder"Library_TRIMMED_15bp_Carried_Variants.fasta")
	flash2_output=$(echo "$output_dir""$master_prefix"".extendedFrags.fastq.gz")
	sam_output=$(echo "$output_dir""$master_prefix"".sam")

	
	echo "bsub -G team151 -o $outfile_converting -q normal -n$pc -w\"done($name_aligning)\" -J $name_converting -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\"  >> $output
        echo "\"bwa samse $REFERENCE $sai_output $flash2_output > $sam_output\"" >> $output

        echo -e "\n"  >> $output

	# This step converts sam to bam and removes not primary alignment reads
	
	type=$(echo "sam_to_bam")
        name_sam_to_bam=$(echo "$batch""_""$type""_""$name_sample""_job")

	output_bwa_bam=$(echo "$output_dir""$master_prefix"".bam")
#	MAPQ=$(echo "10")
	
	echo "bsub -G team151 -o $outfile_converting -q normal -n$pc -w\"done($name_converting)\" -J $name_sam_to_bam -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\"  >> $output
        echo "\"samtools view -@$pc -h -F 256 $sam_output|samtools sort -@$pc -o $output_bwa_bam\"" >> $output

        echo -e "\n"  >> $output

	# This step initiates the sub bash script because it uses a different conda environment

	type=$(echo "filter_and_rest")
        name_filter_and_rest=$(echo "$batch""_""$type""_""$name_sample""_job")
	bash_subscript=$(echo "$bash_subscript_folder""268_MPRA_Sample_procesing_UMI_TOOLS_v2_sub_bash.sh")


	
	echo "bsub -G team151 -o $outfile_converting -q normal -n$pc -w\"done($name_sam_to_bam)\" -J $name_filter_and_rest -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\"  >> $output
        echo "\"bash $bash_subscript $path_to_leave_output_files $Log_files_path $batch $mem $pc $queue $master_prefix $name_sample\"" >> $output

        echo -e "\n"  >> $output



	echo "------------------------------------------->$counter\n"
	counter=0
        array=()
        array_Merge=()


	echo "#########################################################################################################################################################################"  >> $output
	echo "#########################################################################################################################################################################"  >> $output
	echo "#########################################################################################################################################################################"  >> $output
	echo "#########################################################################################################################################################################"  >> $output

	echo -e "\n"  >> $output
	echo -e "\n"  >> $output
    
       fi



  fi #no spaces
done < "$sample_file"


bash $output














