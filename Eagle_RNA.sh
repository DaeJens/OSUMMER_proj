#!/bin/bash

#############################################################################
#
#   This bash script is intended to perform an analysis pipeline of RNA-seq
#   data. It automatically detects the input files from the toAnalyze
#   directory. It expects the input to be formatted into two files, a
#   forward and a reverse read file (ex: sample_R1.fastq.gz, sample_R2.fastq.gz).
#   It then aligns the sequences and marks duplicates. It references
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
#   refFiles/genome.fa                      Reference  file
#   refFiles/indelintervals.list            Target intervals (from GATK RTC)
#   refFiles/dbSNP137.indels.sort.vcf       Known indel sites
#   refFiles/dbSNP137.snps.sort.vcf         Known snp sites
#   refFiles/chr.00-All.sort.vcf            Known mutation sites
#   refFiles/mm10baits.interval_list        Bait intervals
#
#############################################################################


# This function confirms with the user that they want to run the analysis
# on the samples that it detected. Also serves as a manual check to make sure
# the script detected the sample names correctly
function userConfirmation(){
    while []
    do
        echo -n "Start Analysis on above samples? y/n: ">>/dev/tty
        read -n 1 answer    # Grabs user input

        # If answer is yes it goes to comparisonConfirmation, if the answer
        # is no then it closes the program. Otherwise it reprompts the user.
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
if [ ! -f "./refFiles/GRCm39/starIndex/SA" ]; then
    echo "Reference genome file has not been indexed.">>/dev/tty
    echo "Please run ""'""STAR index refFiles/genome.fa""'"" to index it.">>/dev/tty
    echo "Execution Aborted.">>/dev/tty
    exit
fi

#Extracts the sample names from the toAnalyze folder grabs every other
#file since each sample should have an R1 and R2 file.
shopt -s nullglob
array=(Analyze_RNA/*)
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
logfile="logfiles/log_GRC39_RNA$(date +"%m%d%Y%H").txt"
echo -e "----------Starting new Execution:$(date +"%T")----------">>$logfile

#Redirects stderr to the commandoutput text file
commandoutput="logfiles/commandoutput$(date +"%m%d%Y%H").txt"
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
if [ ! -d "output_GRCm39_RNA" ]; then
    mkdir output_GRCm39_RNA
fi
for i in ${sampleNames[@]};
do
    i+=/
    foldername=output_GRCm39_RNA/$i
    if [ -d  $foldername ]; then
        echo "Temp folder $i already exists.">>$logfile
    else
        mkdir $foldername
        echo "Created output folder $i.">>$logfile
    fi
done

printf "STAR alignment  $(date)  \n">>/dev/tty
#Starting STAR alignment
echo -e "\nBeginning STAR alignment ($(date +"%T")):">>$logfile
for i in ${sampleNames[@]};
do
    R1="Analyze_RNA/null_R1.fq.gz"
    R2="Analyze_RNA/null_R2.fq.gz"
    header=''"'"'@RG\tID:null\tLB:null\tSM:null\tPL:ILLUMINA'"'"
    output="output_GRCm39_RNA/null/null."
    check="output_GRCm39_RNA/null/null.bam"
    command="STAR --genomeDir /home/OSUMC.EDU/jens45/00_refFiles/GRCm39/starIndex --genomeLoad LoadAndKeep --readFilesCommand zcat --readFilesIn ${R1/null/$i} ${R2/null/$i} --outFileNamePrefix ${output//null/$i} --runThreadN 4 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMattrIHstart 0 --alignSoftClipAtReferenceEnds No &"
    #Checks if the bam file already exists since sam is deleted later
    if [ -f ${check//null/$i} ]; then
        echo "Skipping $i, ${output//null/$i} already exists.">>$logfile
        continue
    fi
    echo "Running STAR alignment on ${R1//null/$i}.">>$logfile
    eval $command
done
wait
echo -e "DONE WITH ALIGNMENT: $(date +"%T")\n" >> $logfile

STAR --genomeLoad Remove  
printf "samtools view- converting to bam   $(date)  \n">>/dev/tty
#starting samtools sam to bam conversion
echo -e "\nBeginning samtools view ($(date +"%T")):">>$logfile
for i in ${sampleNames[@]};
do
    bam="output_GRCm39_RNA/null/null.bam"
    sam="output_GRCm39_RNA/null/null.Aligned.out.sam"
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


printf "samtools sort   $(date)   \n">>/dev/tty
#Starting samtools sort
echo -e "\nBeginning samtools sort ($(date +"%T")):">>$logfile
for i in ${sampleNames[@]};
do
    #TEMPORARY!!!#############################################
    todelete="output_GRCm39_RNA/null/null.Aligned.out.sam"   #
    rm ${todelete//null/$i}                                  #
    #TEMPORARY!!!#############################################   
    
    bam="output_GRCm39_RNA/null/null.bam"
    output="output_GRCm39_RNA/null/null_sorted.bam"
    command="samtools sort -@ 4 -o ${output//null/$i} ${bam//null/$i}  &"
    filename=$output".bam"
    if [ -f ${output//null/$i}]; then
        echo "Skipping $i, ${output//null/$i} already exists.">>$logfile
        continue
    fi
    echo "Running samtools sort on ${bam//null/$i}.">>$logfile
    eval $command
done
wait
echo -e "DONE WITH SORTING: $(date +"%T")\n" >>$logfile


printf "samtools index   $(date)     \n">>/dev/tty
#Starting samtools index
echo -e "\nBeginning samtools index ($(date +"%T")):">>$logfile
for i in ${sampleNames[@]};
do 
    bam="output_GRCm39_RNA/null/null_sorted.bam"
    output="output_GRCm39_RNA/null/null_sorted.bai"
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

printf "AddReadGroups   $(date)    \n">>/dev/tty 
#add read groups
echo -e "\nBeginning Add Read Groups ($(date +"%T")):">>$logfile
for i in ${sampleNames[@]};
do
    input="output_GRCm39_RNA)/null/null_sorted.bam"
    output="output_GRCm39_RNA/null/null_sortedr.bam"
        if [ -f ${output//null/$i} ]; then
        echo "Skipping $i, ${output//null/$i} already exists.">>$logfile
        continue
    fi
    command="picard-tools AddOrReplaceReadGroups I=${input//null/$i} O=${output//null/$i} SO=coordinate RGID={$i} RGLB=lib1 RGPL=ILLUMINA RGPU=machine RGSM={$i}} &"
    echo -e "Running picard Add Read Groups on ${input//null/$i}.">>$logfile
    eval $command
done
wait
echo -e "DONE adding Read Groups: $(date +"%T")\n" >>$logfile


printf "MarkDuplicates   $(date)       \n">>/dev/tty
#Starting picard MarkDuplicates
echo -e "\nBeginning picard MarkDuplicates ($(date +"%T")):">>$logfile
for i in ${sampleNames[@]};
do
    input="output_GRCm39_RNA/null/null_sortedr.bam"
    output="output_GRCm39_RNA/null/null_sorted.dup.bam"
    metrics="output/null/null.duplicate_metrics"
    if [ -f ${output//null/$i} ]; then
        echo "Skipping $i, ${output//null/$i} already exists.">>$logfile
        continue
    fi
    command="picard-tools MarkDuplicates INPUT=${input//null/$i} OUTPUT=${output//null/$i} METRICS_FILE=${metrics//null/$i} VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true &"
    echo -e "Running picard MarkDuplicates on ${input//null/$i}.">>$logfile
    eval $command
done
wait
echo -e "DONE WITH MARKING DUPLICATES: $(date +"%T")\n" >>$logfile

printf "Done! Back to work!                                       \n">>/dev/tty

