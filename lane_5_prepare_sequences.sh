#!/bin/bash
# Enable history expantion in order to save the comands to the log-file.
set -o history -o histexpand

# Options
# fastx_trimmer
	# First base to keep
	F=14

# cutadapt
	# Adaptors to trim
	A1="GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG"	# TruSeq Adapter, Index 1
	A2="GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG"	# TruSeq Adapter, Index 2
	A6="GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG"	# TruSeq Adapter, Index 6
	A7="GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG"	# TruSeq Adapter, Index 7
	A8="GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG"	# TruSeq Adapter, Index 8
	A9="AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"		# Illumina Single End PCR Primer 1


	# quality threshold for the 3' trimming.
	Q=20	

	# Minimum number of bases overlaping with the adaptor sequence.
	O=10

	# Maximum error rate when identifying adaptors. 2> cutadapt.log
	E="0.1"	
	N=1
	M=50

# fastq_quality_filter
	K=20
	P=95

# Log-files
LOGFILE=${PWD##*/}.log
ERROR=${PWD##*/}.err
CU_LOG=${PWD##*/}.cutadapt.log
PAIR_LOG=${PWD##*/}.pairSeq.log
date >> $LOGFILE
printf "\n" >> $LOGFILE
# Log versions of programs used
echo "# Versions of software used:" >> $LOGFILE
echo $(fastx_trimmer -h) >> $LOGFILE
printf "\n" >> $LOGFILE
echo "cutadapt: " $(cutadapt --version) >> $LOGFILE
printf "\n" >> $LOGFILE

# Print some informative error meassages
err() {
    echo "$1 exited unexpectedly";
        exit 1;
	}

# Function for checking the exit code of a child process
ckeckExit() {
if [ "$1" == "0" ]; then
	printf "[Done] $2 `date`\n\n" >> $LOGFILE;
	else
		err "[Error] $2 returned non-0 exit code in $PWD" >> $LOGFILE;
		fi
		}

#for dir in /nobackup/data11/Perca_fluviatilis/data/J.Sturve_15_01/*
#do
#	cd $dir/*
	echo $PWD
   	fastqc *.fastq.gz
	# Gunzip the fastq files
	gunzip $(ls 5*_1.fastq.gz) 2> $ERROR
	gunzip $(ls 5*_2.fastq.gz) 2> $ERROR

	file1=$PWD/$(ls 5*_1.fastq)
	file2=$PWD/$(ls 5*_2.fastq)

	# Trim 5' end of reads
	printf "# fastx_trimmer\n" >> $LOGFILE
	printf "[ `date` ]\n" >> $LOGFILE
	fastx_trimmer -Q33 -i $file1 -f $F -o "${file1%.fastq}.FXT.fastq" 2>> $ERROR
		echo !! >> $LOGFILE
		ckeckExit $? "fastx_trimmer on file 1"
	printf "[ `date` ]\n" >> $LOGFILE
	fastx_trimmer -Q33 -i $file2 -f $F -o "${file2%.fastq}.FXT.fastq" 2>> $ERROR
		echo !! >> $LOGFILE
		ckeckExit $? "fastx_trimmer on file 2"
	
	# Remove adaptors
	printf "\n# cutadapt\n" >> $LOGFILE
	printf "[ `date` ]\n" >> $LOGFILE
	cutadapt -b $A1 -b $A2 -b $A6 -b $A7 -b $A8 -b $A9 -q $Q -O $O -e $E -n $N -m $M -o "${file1%.fastq}.FXT.CA.fastq" "${file1%.fastq}.FXT.fastq" >> $CU_LOG 2>> $ERROR
		echo !! >> $LOGFILE
		ckeckExit $? "cutadapt on file 1"
	printf "[ `date` ]\n" >> $LOGFILE
	cutadapt -b $A1 -b $A2 -b $A6 -b $A7 -b $A8 -b $A9 -q $Q -O $O -e $E -n $N -m $M -o "${file2%.fastq}.FXT.CA.fastq" "${file2%.fastq}.FXT.fastq" >> $CU_LOG 2>> $ERROR
		echo !! >> $LOGFILE
		ckeckExit $? "cutadapt on file 2"

	# Quality filtering
	printf "\n# fastq_quality_filter\n" >> $LOGFILE
	printf "[ `date` ]\n" >> $LOGFILE
	fastq_quality_filter -Q33 -q $K -p $P -i "${file1%.fastq}.FXT.CA.fastq" -o "${file1%.fastq}.FXT.CA.FQF.fastq" 2>> $ERROR
		echo !! >> $LOGFILE
		ckeckExit $? "fastq_quality_filter on file 1"
	printf "[ `date` ]\n" >> $LOGFILE
	fastq_quality_filter -Q33 -q $K -p $P -i "${file2%.fastq}.FXT.CA.fastq" -o "${file2%.fastq}.FXT.CA.FQF.fastq" 2>> $ERROR
		echo !! >> $LOGFILE
		ckeckExit $? "fastq_quality_filter on file 2"

	# Sort paired and singlet reads
	printf "\n# pairSeq.py\n" >> $LOGFILE
	printf "[ `date` ]\n" >> $LOGFILE
	pairSeq.py "${file1%.fastq}.FXT.CA.FQF.fastq" "${file2%.fastq}.FXT.CA.FQF.fastq" " " 2>> $ERROR >> $PAIR_LOG
		echo !! >> $LOGFILE
		ckeckExit $? "fastq_quality_filter on file 1"
	
	# Quality check
	fastqc "${file1%.fastq}.FXT.CA.FQF.fastq"
	fastqc "${file2%.fastq}.FXT.CA.FQF.fastq"

#done
