#!/bin/bash>

 
#### Rscript

input_file=$1
MASTER_ROUTE=$2
batch=$3
mem=$4
pc=$5
queue=$6
output=$7


touch $output
echo -n "" > $output

echo "#!/bin/bash"  >> $output



output_dir=$(echo "$MASTER_ROUTE")
type_for_R=$(echo "100_percent")
declare -a arr


output_dir2=$MASTER_ROUTE
output_dir=$MASTER_ROUTE


end_of_file=0
counter=0
array=()
array_Merge=()                  #

cd $output_dir2

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
    name_job=$(echo "$batch""_""$name_sample""_""$counter")
    Type=$(echo "$name_sample"|sed -r 's/R[0-9]+_//g')
    Replicate=$(echo "$name_sample"|sed -r 's/_.+//g')

    echo "------------------------>name_job: $name_job"
    echo "------------------------>$name_sample $Type $Replicate"
    echo "------------------------>$Lane_1 $Lane_2 $Merge"

    stderrfile2=$(echo "$output_dir2""/""Lane_merger_"$name_sample"_"$Type".out")

    touch $stderrfile2

    echo -n "" > $stderrfile2

    

    
    echo "bsub -G team151 -o $stderrfile2 -q normal -n$pc -J $name_job -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\" >> $output
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



        outfile_Trimmomatic=$(echo "$output_dir2""/""$name_sample""_Processing.out")
        touch $outfile_Trimmomatic
        echo -n "" > $outfile_Trimmomatic

        name_trimmomatic=$(echo "$batch""_""Trimmomatic""_""$name_sample""_job")
        output_trimmomatic_r3=$(echo "$Type""_""Trimmed_""$R3")
        echo "->$name_trimmomatic $output_trimmomatic_r3<-"

        echo "bsub -G team151 -o $outfile_Trimmomatic -q normal -n$pc -w\"$string3\" -J $name_trimmomatic -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\" >> $output
#        echo "bsub -G team151 -o $outfile_Trimmomatic -q normal -n$pc -J $name_trimmomatic -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\" >> $output
        echo "\"java -jar /nfs/users/nfs_m/mt19/sOFTWARE/Trimmomatic-0.39/trimmomatic-0.39.jar SE $R3 $output_trimmomatic_r3 HEADCROP:6\"" >> $output

	echo -e "\n"  >> $output

	type=$(echo "Extracting_UMI_R1")
	outfile_Extracting_UMI_R1=$(echo "$output_dir""outfile""_""$type""_""$name_sample"".out")
	touch $outfile_Extracting_UMI_R1
	echo -n "" > $outfile_Extracting_UMI_R1
	name_Extracting_UMI_R1=$(echo "$batch""_""$type""_""$name_sample""_job")
	
	output_processed_R1=$(echo "$output_dir""$master_prefix""_R1_Plus_UMIS"".fastq.gz")


	echo "bsub -G team151 -o $outfile_Extracting_UMI_R1 -M $mem -w\"done($name_trimmomatic)\" -J $name_Extracting_UMI_R1 -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
	echo "\"umi_tools extract --bc-pattern=NNNNNNNNNN --stdin=$UMI --read2-in=$R1 --stdout=$output_processed_R1 --read2-stdout\"" >> $output


	echo -e "\n"  >> $output
	
	type=$(echo "Extracting_UMI_R3")
	outfile_Extracting_UMI_R3=$(echo "$output_dir""outfile""_""$type""_""$name_sample"".out")
	touch $outfile_Extracting_UMI_R3
	echo -n "" > $outfile_Extracting_UMI_R3
	name_Extracting_UMI_R3=$(echo "$batch""_""$type""_""$name_sample""_job")

	
	output_processed_R3=$(echo "$output_dir""$master_prefix""_R3_Plus_UMIS"".fastq.gz")


	echo "bsub -G team151 -o $outfile_Extracting_UMI_R3 -M $mem -w\"done($name_trimmomatic)\" -J $name_Extracting_UMI_R3 -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
	echo "\"umi_tools extract --bc-pattern=NNNNNNNNNN --stdin=$UMI --read2-in=$output_trimmomatic_r3 --stdout=$output_processed_R3 --read2-stdout\"" >> $output


	echo -e "\n"  >> $output

	 type=$(echo "R1_R3_merging")
        outfile_merging=$(echo "$output_dir""outfile""_""$type""_""$name_sample"".out")
        touch $outfile_merging
        echo -n "" > $outfile_merging
         name_merging=$(echo "$batch""_""$type""_""$name_sample""_job")
	 flash2_output=$(echo "$output_dir""$master_prefix"".extendedFrags.fastq.gz")
	
	
	echo "bsub -G team151 -o $outfile_merging -q normal -n$pc -w\"done($name_Extracting_UMI_R1) && done($name_Extracting_UMI_R3)\" -J $name_merging -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\"  >> $output
        echo "\"flash2 -z --output-prefix=$master_prefix -t 1 -r 11 -f 11 -s 2 -m 10 -x 0.122 $output_processed_R1 $output_processed_R3\"" >> $output
	
	echo -e "\n"  >> $output

        type=$(echo "CLEANING_1")
        name_CLEANING=$(echo "$batch""_""$type""_""$name_sample""_job")

        echo "bsub -G team151 -o $outfile_merging -M 4000  -w\"done($name_merging)\" -J $name_CLEANING  -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q $queue -- \\" >> $output
        echo "\"rm $R1 $R2 $output_trimmomatic_r3 $R3\"" >> $output

	echo -e "\n"  >> $output


	echo "#########################################################################################################################################################################"  >> $output
	echo "#########################################################################################################################################################################"  >> $output


	type=$(echo "aligning")
        outfile_aligning=$(echo "$output_dir""outfile""_""$type""_""$name_sample"".out")
        touch $outfile_aligning
        echo -n "" > $outfile_aligning
        name_aligning=$(echo "$batch""_""$type""_""$name_sample""_job")
	REFERENCE=$(echo "/lustre/scratch115/teams/soranzo/projects/NEW_MPRA/reference_files/Library_TRIMMED_15bp_Carried_Variants.fasta")

	sai_output=$(echo "$output_dir""$master_prefix"".sai")	
	
	echo "bsub -G team151 -o $outfile_aligning -q normal -n$pc -w\"done($name_merging)\" -J $name_aligning -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\"  >> $output
#	echo "bsub -G team151 -o $outfile_aligning -q normal -n$pc -J $name_aligning -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\"  >> $output
        echo "\"bwa aln -l 10 -O 100 -E 100 $REFERENCE $flash2_output > $sai_output\"" >> $output

        echo -e "\n"  >> $output

	 
	type=$(echo "converting_and_filtering")
        outfile_converting=$(echo "$output_dir""outfile""_""$type""_""$name_sample"".out")
        touch $outfile_converting
        echo -n "" > $outfile_converting
        name_converting=$(echo "$batch""_""$type""_""$name_sample""_job")
	REFERENCE=$(echo "/lustre/scratch115/teams/soranzo/projects/NEW_MPRA/reference_files/Library_TRIMMED_15bp_Carried_Variants.fasta")
	flash2_output=$(echo "$output_dir""$master_prefix"".extendedFrags.fastq.gz")
	sam_output=$(echo "$output_dir""$master_prefix"".sam")

	
	echo "bsub -G team151 -o $outfile_converting -q normal -n$pc -w\"done($name_aligning)\" -J $name_converting -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\"  >> $output
        echo "\"bwa samse $REFERENCE $sai_output $flash2_output > $sam_output\"" >> $output

        echo -e "\n"  >> $output


	type=$(echo "sam_to_bam")
        name_sam_to_bam=$(echo "$batch""_""$type""_""$name_sample""_job")

	output_bwa_bam=$(echo "$output_dir""$master_prefix"".bam")
#	MAPQ=$(echo "10")
	
	echo "bsub -G team151 -o $outfile_converting -q normal -n$pc -w\"done($name_converting)\" -J $name_sam_to_bam -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\"  >> $output
        echo "\"samtools view -@$pc -h -F 256 $sam_output|samtools sort -@$pc -o $output_bwa_bam\"" >> $output

        echo -e "\n"  >> $output

	type=$(echo "filter_and_rest")
        name_filter_and_rest=$(echo "$batch""_""$type""_""$name_sample""_job")


	
	echo "bsub -G team151 -o $outfile_converting -q normal -n$pc -w\"done($name_sam_to_bam)\" -J $name_filter_and_rest -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\"  >> $output
        echo "\"bash /nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/268_MPRA_Sample_procesing_UMI_TOOLS_v2_sub_bash.sh $MASTER_ROUTE $batch $mem $pc $queue $master_prefix $name_sample\"" >> $output

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
done < "$input_file"


echo "bash $output"














