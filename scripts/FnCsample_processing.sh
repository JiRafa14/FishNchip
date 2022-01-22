#! /bin/bash

#Reading parameters

## sampledir
## ind (para que sepa el script si es el 1, el 2...). Debera haber para ctrl=TRUE : IND=1, IND=2... ;; para ctrl=FALSE : IND=1 ..
## paired
## ctrl = TRUE/FALSE (para que sepa si tiene que tirar el script para chip o para ctrl

SAMPLEDIR=$1
IND=$2
PAIRED=$3
CTRL=$4
TOTALSAMP=$5
EXPNUMBER=$6
USESUMMITS=$7
WORDLENGTH=$8
SIZE=$9
PVAL=${10}
TFNAME=${11}
UPLEN=${12}
DOWNLEN=${13}

#Procesing samples: quality analysis and reads mapping

cd $SAMPLEDIR

if [ $PAIRED = FALSE ]
then
	if [ $CTRL = FALSE ]
	then
		echo "Processing unpaired chip $IND"

	        cd chip
		fastqc chip${IND}.fastq.gz

                bowtie2 -x ../../genome/index -U chip${IND}.fastq.gz -S chip${IND}.sam
                samtools sort -o chip${IND}.bam chip${IND}.sam
                samtools index chip${IND}.bam
		rm chip${IND}.sam
		bamCoverage -bs 5 --normalizeUsing CPM --bam chip${IND}.bam -o chip${IND}.bw


		cd ..
	elif [ $CTRL = TRUE ]
	then
		echo "Processing unpaired control $IND"

                cd ctrl
		fastqc ctrl${IND}.fastq.gz

                bowtie2 -x ../../genome/index -U ctrl${IND}.fastq.gz -S ctrl${IND}.sam
                samtools sort -o ctrl${IND}.bam ctrl${IND}.sam
                samtools index ctrl${IND}.bam
		rm ctrl${IND}.sam
		bamCoverage -bs 5 --normalizeUsing CPM --bam ctrl${IND}.bam -o ctrl${IND}.bw

		cd ..
	fi


elif [ $PAIRED = TRUE ]
then

	if [ $CTRL = FALSE ]
	then
		echo "Processing paired chip $IND"

		cd chip
		fastqc chip${IND}_A.fastq.gz
		fastqc chip${IND}_B.fastq.gz

                bowtie2 -x ../../genome/index -1 chip${IND}_A.fastq.gz -2 chip${IND}_B.fastq.gz -S chip${IND}.sam
                samtools sort -o chip${IND}.bam chip${IND}.sam
                samtools index chip${IND}.bam
		rm chip${IND}.sam
		bamCoverage -bs 5 --normalizeUsing CPM --bam chip${IND}.bam -o chip${IND}.bw


		cd ..

	elif [ $CTRL = TRUE ]
	then
		echo "Processing paired control $IND"

                cd ctrl
		fastqc ctrl${IND}_A.fastq.gz
		fastqc ctrl${IND}_B.fastq.gz

                bowtie2 -x ../../genome/index -1 ctrl${IND}_A.fastq.gz -2 ctrl${IND}_B.fastq.gz -S ctrl${IND}.sam
                samtools sort -o ctrl${IND}.bam ctrl${IND}.sam
                samtools index ctrl${IND}.bam

		rm ctrl${IND}.sam

		bamCoverage -bs 5 --normalizeUsing CPM --bam ctrl${IND}.bam -o ctrl${IND}.bw


		cd ..
	fi

fi

echo ""
echo ".sam files removed"
echo "Generated .bam .bai and .bw files"
echo ""



## Write in blackboard

if [ $CTRL = FALSE ]
then

	echo "chip ${IND} is finished" >> ../logs/everyone_present.txt

elif [ $CTRL = TRUE ]
then
	echo "control ${IND} is finished" >> ../logs/everyone_present.txt

fi

## Read Blackboard

NSAMPLE=$( wc -l ../logs/everyone_present.txt | awk '{ print $1 }' )

if [ $TOTALSAMP -eq $NSAMPLE ]
then

	echo ""
	echo "all samples finished"
	echo "proceeding to peak analysis"

	for IND in `seq 1 ${EXPNUMBER}`
	do
        	cd ../logs

	        sbatch ../scripts/peakAnalysis.sh $IND $SAMPLEDIR $USESUMMITS $EXPNUMBER $WORDLENGTH $SIZE $PVAL $TFNAME $UPLEN $DOWNLEN

		cd ../results
	done



fi


