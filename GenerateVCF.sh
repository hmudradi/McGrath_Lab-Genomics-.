#!/bin/bash

reads1=0 #a-input reads file 1
reads2=0 #b-input reads file 2
ref=0 #r-Reference genome file (path or file name if in current dir)
realign=0 #e-perform read re-alignment (flag)
output=0 #o-Output VCF file name
millsFile=0 #f-Mills file location
gunzip=0 #z-Output VCF file should be in gunzipped (*.vcf.gz --> *.vcf)
v=0 #v-verbose mode (flag)
index=0 #i-Index the output BAM file (use 'samtools index' command)
help=0 #h-print how to use info (what arguments it takes in) if this is present, dont do any other options
answer=0 #variable for checking whether the input file needs to be overwritten
picard=0 #variable to assign location of picard.jar

while getopts "a:b:r:eo:f:zvip:h" option
do
	case $option in
	a) reads1=$OPTARG;;
	b) reads2=$OPTARG;;
	r) ref=$OPTARG;;
	e) realign=1;;
	o) output=$OPTARG;;
	f) millsFile=$OPTARG;;
	z) gunzip=1;;
	v) v=1;echo "Verbose Mode is on";;
	i) index=1;;
	p) picard=$OPTARG;;
	h) help=1;echo "Welcome to the help Page";;
	esac
done

#file check for if the input files exist
if [ -f *1.fq -a -f *2.fq ]
then
	echo "both of the paired reads files exist"
elif [ ! -f *1.fq -a ! -f *2.fq ]
then
	echo "neither of the paired reads files exist"
	exit
elif [ ! -f *1.fq ]
then
	echo "first paired reads file does not exist"
	exit
elif [ ! -f *2.fq ]
then
	echo "second paired reads file does not exist"
	exit
fi

#file check to see if reference genome exists
if [ -f *.fa ]
then
	echo "the reference genome exists"
elif [ ! -f *.fa ]
then
	echo "the reference genome does not exist or it is compressed"
	exit
fi

if [ $help = 1 ]
then
	echo "
	NOTE: You must have bwa, samtools, bcftools, and GATK version 3.7.0 installed to use this script
	GenomeAnalysisTK.jar must be present in the directory in where you call this script

	Options:
	-a input the name of the first reads file. if the file is in another directory, specify the file path.
	-b input the name of the second reads file. if the file is in another directory, specify the file path.
	-r input the name of the reference genome file or the path to the file. The file must be uncompressed.
	-e this flag toggles on the realignment function of the script.
	-o this option allows the user to specify the output file name of the VCF file that will be produced
	-f input the name of the used Mills file or the path to the file
	-z this flag will automatically unzip the *.vcf.gz file produced into a *.vcf file
	-v verbose mode will be turned on
	-i this option will index the output BAM file
	-h this flag will open this help page"
else
	#prepare the reference for mapping by giving bwa the path to the reference file
	if [ $v = 1 ]
        then
                echo "The reference genome is being prepared for mapping"
        fi
	bwa index $ref
	echo "ref seq prepped for mapping"

	#map the reads to the reference genome (time consuming)
	if [ $v = 1 ]
        then
                echo "The paired end reads are being mapped to the reference genome"
        fi
	bwa mem -R '@RG\tID:foo\tSM:bar\tLB:library1' $ref $reads1 $reads2 > lane.sam
	echo "paired reads mapped to ref seq"

	#clean up the read pairing info and flags
	if [ $v = 1 ]
	then
		echo "the read pairing information is being cleaned up"
	fi
	samtools fixmate -O bam lane.sam lane_fixmate.bam
	echo "the read paring info has been cleaned"

	#sort into coordinate order
	if [ $v = 1 ]
	then
		echo "the the mapped reads are being sorted into coordinate order"
	fi
	samtools sort -O bam -o lane_sorted.bam -T /tmp/lane_temp lane_fixmate.bam
	echo "the mapped reads have been sorted into coordinate order"

	#Start Improvements:
	#first generate the ref seq .fai and .dict files
	if [ $v = 1 ]
	then
		echo "the *.fai and *.dict files are being generated"
	fi
	samtools faidx $ref
	java -jar $picard CreateSequenceDictionary -R $ref
	echo "the .fai and .dict files were generated"

	#index the sorted bam file for the next 2 java commads (otherwise, they throw errors saying that the BAM file is not indexed)
	if [ $v = 1 ]
	then
		echo "the sorted BAM file is being indexed for realignment"
	fi
	samtools index lane_sorted.bam
	echo "the sorted BAM file has been indexed and is ready for realignment"

	#reduce the number of miscalls of INDELS in the data by realigning with the next 2 commands
	if [ $v = 1 ]
	then
		echo "the reads are being realigned and imrpoved using GATK"
	fi
	#if statement to prompt realignment if the -e flag is called
	if [ $realign = 1 ]
	then
		java -Xmx2g -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R $ref -I lane_sorted.bam -o lane.intervals --known $millsFile
		java -Xmx4g -jar GenomeAnalysisTK.jar -T IndelRealigner -R $ref -I lane_sorted.bam -targetIntervals lane.intervals -known $millsFile -o lane_realigned.bam
		echo "The realignment occured since -e was called"

		if [ $index = 1 ]
		then
			#index the realigned bam file
			if [ $v = 1 ]
			then
				echo "the realigned BAM file is being indexed"
			fi
                	samtools index lane_realigned.bam
                	echo "the indexing and realignment was done"

			#lastly, call the variants and generate the .vcf file
			#if statement to check if the vcf.gz file already exists
			if [ -f $output.vcf* ]
			then
        			echo "The vcf file already exists"
        			read -p "How would you like to proceed? Enter overwrite or exit: " option

        			case $option in
                			overwrite)if [ $v = 1 ]
                        			then
                                			echo "the variants are being called and the $output.vcf file is beign overwritten"
                        			fi
						bcftools mpileup -Ou -f $ref lane_realigned.bam | bcftools call -vmO z -o $output.vcf.gz
                			echo "the file was overwritten"
					exit;;
                			exit) exit;;
        			esac
			fi

			#verbose mode check for creating .vcf file
			if [ $v = 1 ]
			then
				echo "the variants are being called and the $output.vcf.gz file is geing generated"
			fi
                	bcftools mpileup -Ou -f $ref lane_realigned.bam | bcftools call -vmO z -o $output.vcf.gz
                	echo "the vcf generation and realignment was done"
		else
			echo "The -i flag was not called so indexing was not performed. Script cannot proceed without indexing the output BAM"
		fi
	else
		echo "The -e flag was not called and realignment was not performed. Script cannot proceed without realignment"
	fi
fi

#z flag conditional to unzip the produced *.vcf.gz file
if [ $gunzip = 1 ]
then
	if [ $v = 1 ]
	then
        	echo "the $output.vcf.gz file is being unzipped"
	fi
	gunzip *.vcf.gz
	echo "the *.vcf.gz file was unzipped"
fi