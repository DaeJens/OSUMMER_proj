#!/bin/bash

#############################################################################
#
#   This bash script is intended to perform an analysis pipeline of genomic
#   data. It automatically detects the input files from the toAnalyze
#   directory. It expects the input to be formatted into two files, a
#   forward and a reverse read file (ex: sample_R1.fastq, sample_R2.fastq).
#   It then uses a variety of analysis toolkits including bwa, samtools,
#   picard, and GATK to manipulate and analyze the samples. It references
#   files in the refFiles directory. All required reference files are shown
#   below. The output files are put in the output directory (ex. output/sample)
#   The script will create the output directory if it does not already exist.
#   Before each step, the script checks if there is already a file in the 
#   output directory which matches the output name. If so it assumes this step
#   has already been performed and skips it. The script outputs information
#   into a logfile that shows how long each step takes and what was ran or
#   skipped.
#
#   Command to run:  ./Eagle
#
#   Expected reference Files:
#
#   refFiles_GRCm39/GRCm39.dna_sm.fa                    Reference  file
#   refFiles/GRCm39_structural_variations.vcf	        Known structural variations
#   refFiles/GRCm39_snp.vcf         			        Known snp sites
#   refFiles/GRC39_bait_interval.bed		            Bait intervals
#
#############################################################################


#This function confirms with the user that they want to run the analysis
#on the samples that it detected. Also serves as a manual check to make sure
#the script detected the sample names correctly
function userConfirmation(){
    while []
    do
        echo -n "Start Analysis on above samples? y/n: ">>/dev/tty
        read -n 1 answer    #Grabs user input

        #If answer is yes it goes to comparisonConfirmation, if the answer
        #is no then it closes the program. Otherwise it reprompts the user.
        if [ "$answer" = "y" ]; then
            echo ""
            comparisonConfirmation
            return 0
        elif [ "$answer" = "n" ]; then
            echo "Execution aborted.">>$logfile
            echo -e "\nExiting Program.">>/dev/tty
            exit
        else
            echo -e "\nInvalid Selection. Try Again">>/dev/tty
        fi
    done
}

#This function asks the user whether they want to do a comparison between the
#normal sample and the tumor samples. If no is selected it skips the last
#four steps of the analysis and ends with the mpileup command. 
function comparisonConfirmation(){
    while []
    do
        echo -n "Do you wish to perform tumor/normal comparisons? y/n: ">>/dev/tty
        read -n 1 answer    #Grabs user input

        #If the answer is yes it then displays the choices of the samples
        #and prompts the user to choose the normal sample from the list.
        if [ "$answer" = "y" ]; then
            echo ""
            performComparisons=1    #Used later to check whether to skip

            #Displays the choices for the normal sample
            for ((i=0;i<${#sampleNames[@]};i+=1))
            do
                num="$(($i+1))"
                echo "$num) ${sampleNames[i]}">>/dev/tty
            done
            echo "$((${#sampleNames[@]}+1))) I change my mind">>/dev/tty

            #Asks the user to choose which sample
            while []
            do
                echo -n "Select which sample is the normal sample: ">>/dev/tty
                read -n 1 choice
                echo ""
                if [[ $choice != [0-9]* ]]; then
                    echo -e "Invalid Selection. Try Again">>/dev/tty
                elif [ $choice -le ${#sampleNames[@]} -a $choice -gt 0 ]; then
                    normalSample=${sampleNames["$(($choice-1))"]}
                    echo "Tumor/normal comparison will be performed.">>$logfile
                    echo "$normalSample selected as normal sample.">>$logfile
                    echo "$normalSample selected. Nice choice!">>/dev/tty
                    return 0
                elif [ $choice == 4 ]; then
                    echo "You flaked out of the comparison.">>$logfile
                    echo "Tumor/Normal comparison skipped.">>/dev/tty
                    performComparisons=0
                    return 0
                else
                    echo -e "Invalid Selection. Try Again">>/dev/tty
                fi
            done
        elif [ "$answer" == "n" ]; then
            echo "Comparison skipped.">>$logfile
            echo -e "\nComparison skipped.">>/dev/tty
            performComparisons=0
            return 0
        else
            echo -e "\nInvalid Selection. Try Again">>/dev/tty
        fi
    done
}

#Checks for -help tag and prints out the file header
if [ "$#" -ne "0" ] && [ $1 == "-help" ]; then
    head -36 Eagle
    exit
fi

###############################################################################
#
# Error Checking: This section checks that everything is in place so that the
# script will be able to run smoothly. It makes sure that all of the reference
# files are there and that the genome file is indexed. It also check to make
# sure that there are files in the toAnalyze/ directory. 
#
###############################################################################

#Checks to see whether all of the referenceFiles are there.
refFiles=(./refFiles_GRC39/GRCm39.dna_sm.fa ./refFiles_GRC39/GRCm39_snp.vcf ./refFiles_GRC39/GRCm39_structural_variations.vcf ./refFiles_GRC39/GRC39_bait_interval.bed)
for i in ${refFiles[@]}; 
do
    if [ ! -f $i ]; then
        echo "Missing one of the reference files. $i">>/dev/tty
        echo "Execution Aborted.">>/dev/tty
        exit
    fi
done

#Checks to make sure the genome file has been indexed
if [ ! -f "./refFiles_GRC39/GRCm39.bwt" ]; then
    echo "Reference genome file has not been indexed.">>/dev/tty
    echo "Please run ""'""bwa index refFiles/genome.fa""'"" to index it.">>/dev/tty
    echo "Execution Aborted.">>/dev/tty
    exit
fi

#Extracts the sample names from the toAnalyze folder grabs every other
#file since each sample should have an R1 and R2 file.
shopt -s nullglob
array=(toAnalyze/*)
sampleNames=()
arraylen=${#array[@]}
for (( i=0; i<${arraylen}; i+=2));
do
    filename=${array[i]%_*}
    filename=${filename#*/}
    sampleNames[${#sampleNames[@]}]=${filename}
done

if [ ${#sampleNames[@]} == 0 ]; then
    echo "No files to Analyze. Please put input files in toAnalyze/"
    exit
fi

###############################################################################
#
# Beginning of Execution:
# This is where the script starts to execute the pipeline. It first checks
# to with the user to confirm the samples to analyze and then to check if
# they want to do a tumor/normal comparison. It then performs the countdown
# and blasts off!
#
###############################################################################

#Creates the logfile directory if it does not already exist
if [ ! -d "logfiles" ]; then
    mkdir logfiles
fi

#Creates the logfile name by appending the date
logfile="logfiles/log_GRC39$(date +"%m%d%Y").txt"
echo -e "----------Starting new Execution:$(date +"%T")----------">>$logfile

#Redirects stderr to the commandoutput text file
commandoutput="logfiles/commandoutput$(date +"%m%d%Y%T").txt"
exec 2>>$commandoutput

#Gets user confirmation to begin running the program
echo "Samples to analyze:">>/dev/tty
echo ${sampleNames[@]}>>/dev/tty
echo "Samples to analyze:">>$logfile
echo ${sampleNames[@]}>>$logfile
normalSample=""
performComparisons=0
userConfirmation

#Countdown and Blast Off!
printf "Starting Execution in: 3\r">>/dev/tty
sleep 1
printf "Starting Execution in: 2\r">>/dev/tty
sleep 1
printf "Starting Execution in: 1\r">>/dev/tty
sleep 1
printf "Starting Execution in: Blast off!!\r">>/dev/tty
sleep 1
printf " >>[0000])>                                            \r">>/dev/tty
trail="~"
for ((i=0;i<68;i+=1));
do
    printf "$trail >>[0000])>\r">>/dev/tty
    trail+="~"
    sleep 0.05
done

#Checks to see if the folders in output already exist and creates them if not
echo -e "\nCreating output directories:">>$logfile
if [ ! -d "output_GRCm39" ]; then
    mkdir output_GRCm39
fi
for i in ${sampleNames[@]};
do
    i+=/
    foldername=output_GRCm39/$i
    if [ -d  $foldername ]; then
        echo "Temp folder $i already exists.">>$logfile
    else
        mkdir $foldername
        echo "Created output folder $i.">>$logfile
    fi
done

printf "                                         (0%%) bwa alignment                    \r">>/dev/tty
#Starting bwa alignment
echo -e "\nBeginning bwa alignment ($(date +"%T")):">>$logfile
for i in ${sampleNames[@]};
do
    
    R1="toAnalyze/null_R1.fastq"
    R2="toAnalyze/null_R2.fastq"
    header=''"'"'@RG\tID:null\tLB:null\tSM:null\tPL:ILLUMINA'"'"
    output="output_GRCm39/null/null.sam"
    check="output_GRCm39/null/null.bam"
    command="bwa mem -M -t 12 -R ${header//null/$i} ./refFiles_GRC39/GRCm39 ${R1/null/$i} ${R2/null/$i} > ${output//null/$i} &"
    #Checks if the bam file already exists since sam is deleted later
    if [ -f ${check//null/$i} ]; then
        echo "Skipping $i, ${output//null/$i} already exists.">>$logfile
        continue
    fi
    echo "Running bwa alignment on ${R1//null/$i}.">>$logfile
    echo $command
    eval $command
done
wait
echo -e "DONE WITH ALIGNMENT: $(date +"%T")\n" >> $logfile


printf "####                                     (9%%) samtools view            \r">>/dev/tty
#starting samtools sam to bam conversion
echo -e "\nBeginning samtools view ($(date +"%T")):">>$logfile
for i in ${sampleNames[@]};
do
    bam="output_GRCm39/null/null.bam"
    sam="output_GRCm39/null/null.sam"
    command="samtools view -S -b -o ${bam//null/$i} ${sam//null/$i} &"
    #Checks if the bam file already exists
    if [ -f ${bam//null/$i} ]; then
        echo "Skipping $i, ${output//null/$i} already exists.">>$logfile
        continue
    fi
    echo "Running samtools view on ${sam//null/$i}.">>$logfile
    eval $command
done
wait
echo -e "DONE WITH CONVERSIONS: $(date +"%T")\n" >>$logfile


printf "#####                                    (12%%) samtools sort           \r">>/dev/tty
#Starting samtools sort
echo -e "\nBeginning samtools sort ($(date +"%T")):">>$logfile
for i in ${sampleNames[@]};
do
    #TEMPORARY!!!#############################
    todelete="output_GRCm39/null/null.sam"   #
    rm ${todelete//null/$i}                  #
    #TEMPORARY!!!#############################   
    
    bam="output_GRCm39/null/null.bam"
    output="output_GRCm39/null/null_sorted.bam"
    command="samtools sort ${bam//null/$i} -o ${output//null/$i} &"
    filename=$output
    if [ -f ${filename//null/$i} ]; then
        echo "Skipping $i, ${output//null/$i} already exists.">>$logfile
        continue
    fi
    echo "Running samtools sort on ${bam//null/$i}.">>$logfile
    eval $command
done
wait
echo -e "DONE WITH SORTING: $(date +"%T")\n" >>$logfile


printf "######                                   (15%%) samtools index          \r">>/dev/tty
#Starting samtools index
echo -e "\nBeginning samtools index ($(date +"%T")):">>$logfile
for i in ${sampleNames[@]};
do 
    bam="output_GRCm39/null/null_sorted.bam"
    output="output_GRCm39/null/null_sorted.bai"
    command="samtools index ${bam//null/$i} ${output//null/$i} &"
    if [ -f ${output//null/$i} ]; then
        echo "Skipping $i, ${output//null/$i} already exists.">>$logfile
        continue
    fi  
    echo "Running samtools index on ${bam//null/$i}.">>$logfile
    eval $command
done
wait
echo -e "DONE WITH INDEXING: $(date +"%T")\n" >>$logfile


printf "######                                   (16%%) MarkDuplicates          \r">>/dev/tty
#Starting picard MarkDuplicates
echo -e "\nBeginning picard MarkDuplicates ($(date +"%T")):">>$logfile
for i in ${sampleNames[@]};
do
    input="output_GRCm39/null/null_sorted.bam"
    output="output_GRCm39/null/null_sorted.dup.bam"
    metrics="output_GRCm39/null/null.duplicate_metrics"
    if [ -f ${output//null/$i} ]; then
        echo "Skipping $i, ${output//null/$i} already exists.">>$logfile
        continue
    fi
    command="picard-tools MarkDuplicates I=${input//null/$i} O=${output//null/$i} M=${metrics//null/$i} VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true &"
    echo -e "Running picard MarkDuplicates on ${input//null/$i}.">>$logfile
    eval $command
done
wait
echo -e "DONE WITH MARKING DUPLICATES: $(date +"%T")\n" >>$logfile


printf "########                                 (20%%) BuildBamIndex           \r">>/dev/tty
#Starting picard BuildBamIndex
echo -e "\nBeginning picard BuildBamIndex ($(date +"%T")):">>$logfile
for i in ${sampleNames[@]};
do
    input="output_GRCm39/null/null_sorted.dup.bam"
    output="output_GRCm39/null/null_sorted.dup.bai"
    if [ -f ${output//null/$i} ]; then
        echo "Skipping $i, ${output//null/$i} already exists.">>$logfile
        continue
    fi
    command="picard-tools BuildBamIndex INPUT=${input//null/$i} VALIDATION_STRINGENCY=LENIENT &"
    echo -e "Running picard BuildBamIndex on ${input//null/$i}.">>$logfile
    eval $command
done
wait
echo -e "DONE WITH BUILDBAMINDEX: $(date +"%T")\n" >>$logfile


printf "##########                               (25%%) BaseRecalibrator        \r">>/dev/tty
#Starting GATK BaseRecalibrator
echo -e "\nBeginning GATK BaseRecalibrator ($(date +"%T")):">>$logfile
for i in ${sampleNames[@]};
do
    input="output_GRCm39/null/null_sorted.dup.bam"
    output="output_GRCm39/null/null.recal_data.table"
    if [ -f ${output//null/$i} ]; then
        echo "Skipping $i, ${output//null/$i} already exists." >> $logfile
        continue
    fi
    command="gatk BaseRecalibrator -R ${refFiles[0]} -I ${input//null/$i} --known-sites ${refFiles[1]} --known-sites ${refFiles[2]} -O ${output//null/$i} >> $commandoutput &"
    eval $command
    echo "Running BaseRecalibrator on ${input//null/$i}." >> $logfile
done
wait
echo -e "DONE WITH BASERECALIBRATOR: $(date +"%T")\n" >> $logfile


printf "################                         (40%%) PrintReads              \r">>/dev/tty
#Starting GATK PrintReads
echo -e "\nBeginning GATK PrintReads ($(date +"%T")):">>$logfile
for i in ${sampleNames[@]};
do
    input="output_GRCm39/null/null_sorted.dup.bam"
    BQSR="output_GRCm39/null/null.recal_data.table"
    output="output_GRCm39/null/null.realign.recal.bam"
    if [ -f ${output//null/$i} ]; then
        echo "Skipping $i, ${output//null/$i} already exists." >> $logfile
        continue
    fi
    command="gatk ApplyBQSR -R ${refFiles[0]} -I ${input//null/$i} --bqsr-recal-file ${BQSR//null/$i} -O ${output//null/$i} >> $commandoutput &"
    eval $command
    echo "Running PrintReads on ${input//null/$i}." >> $logfile
done
wait
echo -e "DONE WITH PRINTREADS: $(date +"%T")\n" >> $logfile



printf "#########################                (63%%) CollectInsertSizeMetrics\r">>/dev/tty
#Starting picard CollectInsertSizeMetrics
echo "Beginning CollectgInsertSizeMetrics ($(date +"%T")):" >> $logfile
for i in ${sampleNames[@]};
do
    input="output_GRCm39/null/null.realign.recal.bam"
    output="output_GRCm39/null/null.insert_size_metrics.txt"
    hist="output_GRCm39/null/null.insert_size_histogram.pdf"
    if [ -f ${output//null/$i} ]; then
        echo "Skipping $i, ${output//null/$i} already exists." >> $logfile
        continue
    fi
    command="picard-tools CollectInsertSizeMetrics I=${input//null/$i} O=${output//null/$i} H=${hist//null/$i} M=0.5 VALIDATION_STRINGENCY=LENIENT &"
    eval $command
    echo "Running CollectInsertSizeMetrics on ${input//null/$i}." >> $logfile
done
wait
echo -e "DONE WITH COLLECTINSERTSIZEMETRICS: $(date +"%T")\n" >> $logfile


printf "##########################               (64%%) CollectHsMetrics        \r">>/dev/tty
#Starting picard CollectHsMetrics
echo "Beginning CollectHsMetris ($(date +"%T")):" >> $logfile
for i in ${sampleNames[@]};
do
    input="output_GRCm39/null/null.realign.recal.bam"
    output="output_GRCm39/null/null.hs_metrics.txt"
    if [ -f ${output//null/$i} ]; then
        echo "Skipping $i, ${output//null/$i} already exists." >> $logfile
        continue
    fi
    command="picard-tools CollectHsMetrics I=${input//null/$i}  O=${output//null/$i} R=${refFiles[0]} BAIT_INTERVALS=${refFiles[3]} TARGET_INTERVALS=${refFiles[3]} VALIDATION_STRINGENCY=LENIENT &"
    eval $command
    echo "Running CollectHsMetrics on ${input//null/$i}." >> $logfile
done
wait
echo -e "DONE WITH COLLECTHSMETRICS: $(date +"%T")\n" >> $logfile


printf "###########################              (66%%) mpileup                 \r">>/dev/tty
#Starting samtools mpileup
echo "Beginning samtools mpileup ($(date +"%T")):" >> $logfile
for i in ${sampleNames[@]};
do
    input="output_GRCm39/null/null.realign.recal.bam"
    output="output_GRCm39/null/null.mpileup"
    if [ -f ${output//null/$i} ]; then
        echo "Skipping $i, ${output//null/$i} already exists." >> $logfile
        continue
    fi
    command="samtools mpileup -q 1 -f ${refFiles[0]} ${input//null/$i} > ${output//null/$i} &"
    eval $command
    echo "Running mpileup on ${input//null/$i}." >> $logfile
done
wait
echo -e "DONE WITH MPILEUP: $(date +"%T")\n" >> $logfile

if [ "$performComparisons" == 0 ]; then
    echo "Tumor/Normal comparison skipped." >> $logfile
    printf "Done! Back to work!                                   \n">>/dev/tty
    exit
fi



printf "###################################      (88%%) somatic                 \r">>/dev/tty
#Starting VarScan somatic
echo -e "\nBeginning VarScan somatic ($(date +"%T")):">>$logfile
for i in ${sampleNames[@]};
do
    input="output_GRCm39/null/null.mpileup"
    output="output_GRCm39/null/normal_null"
    output=${output//normal/$normalSample}
    if [ $i = $normalSample ]; then
        continue
    elif [ -f ${output//null/$i}.snp.vcf ]; then
        echo "Skipping $i, ${output//null/$i}.snp.vcf already exists.">>$logfile
        continue
    fi
    echo "Running somatic on ${input//null/$i}." >> $logfile
    command="varscan somatic ${input//null/$normalSample} ${input//null/$i} ${output//null/$i} --output-vcf 1 &"
    eval $command
done
wait
echo -e "DONE WITH VARSCAN SOMATIC: $(date +"%T")\n" >> $logfile

#Starting VarScan processSomatic
echo -e "\nBeginning VarScan processSomatic ($(date +"%T")):">>$logfile
for i in ${sampleNames[@]};
do
    input="output_GRCm39/null/normal_null"
    input=${input//normal/$normalSample}
    if [ $i = $normalSample ]; then
        continue
    elif [ -f ${input//null/$i}.snp.Somatic.hc.vcf ]; then
        echo "Skipping $i, output files already exist.">>$logfile
        continue
    fi
    echo "Running VarScan processSommatic on ${input//null/$i}." >> $logfile
    command1="varscan processSomatic ${input//null/$i}.snp.vcf &"
    command2="varscan processSomatic ${input//null/$i}.indel.vcf &"
    eval $command1
    eval $command2
done
wait
echo -e "DONE WITH VARSCAN PROCESSSOMATIC: $(date +"%T")\n" >> $logfile

#Starting grep/awk
echo -e "\nBeginning grep/awk ($(date +"%T")):">>$logfile
for i in ${sampleNames[@]};
do
    input="output_GRCm39/null/normal_null"
    input=${input//normal/$normalSample}
    output="output_GRCm39/null/MOD_normal_null"
    output=${output//normal/$normalSample}
    if [ $i = $normalSample ]; then
        continue
    elif [ -f ${output//null/$i}.snp.Somatic.hc.vcf ]; then
        echo "Skipping $i, output files already exist.">>$logfile
        continue
    fi
    echo "Using grep/awk to reformat ${input//null/$i}." >> $logfile
    command1="grep -v ""'""#""'"" ${input//null/$i}.snp.Somatic.hc.vcf | awk ""'""{print \$1,\$2,\$1,\"@\",\$2,\"@\",\$4,\"@\",\$5,\$4,\$5,\$6,\$7}""'"" - > ${output//null/$i}.snp.Somatic.hc.vcf &"
    command2="grep -v ""'""#""'"" ${input//null/$i}.indel.Somatic.hc.vcf | awk ""'""{print \$1,\$2,\$1,\"@\",\$2,\"@\",\$4,\"@\",\$5,\$4,\$5,\$6,\$7}""'"" - > ${output//null/$i}.indel.Somatic.hc.vcf &"
    eval $command1
    eval $command2
done
wait
echo -e "DONE WITH GREP/AWK: $(date +"%T")\n" >> $logfile

#Starting sed
echo -e "\nBeginning sed ($(date +"%T")):">>$logfile
for i in ${sampleNames[@]};
do
    if [ $i = $normalSample ]; then
        continue
    fi
    echo Using sed to reformat ${input//null/$i}.snp.Somatic.hc.vcf>>$logfile
    input="output_GRCm39/null/MOD_normal_null"
    input=${input//normal/$normalSample}
    sedat="sed -i ""'""s/ @ /_/g""'"" ${input//null/$i}.snp.Somatic.hc.vcf"
    sedtab="sed -i ""'""s/ /\\t/g""'"" ${input//null/$i}.snp.Somatic.hc.vcf"
    sedchr="sed -i ""'""s/chr//g""'"" ${input//null/$i}.snp.Somatic.hc.vcf"
    eval $sedat
    eval $sedtab
    eval $sedchr
    echo Using sed to reformat ${input//null/$i}.indel.Somatic.hc.vcf>>$logfile
    sedat2="sed -i ""'""s/ @ /_/g""'"" ${input//null/$i}.indel.Somatic.hc.vcf"
    sedtab2="sed -i ""'""s/ /\\t/g""'"" ${input//null/$i}.indel.Somatic.hc.vcf"
    sedchr2="sed -i ""'""s/chr//g""'"" ${input//null/$i}.indel.Somatic.hc.vcf"
    eval $sedat2
    eval $sedtab2
    eval $sedchr2
done
echo -e "DONE WITH SED: $(date +"%T")\n" >> $logfile

printf "Done! Back to work!                                       \n">>/dev/tty
