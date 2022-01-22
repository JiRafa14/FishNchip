#! /bin/bash

## This script obtains the transcription factor binding sites (peaks calling)

IND=$1
SAMPLEDIR=$2
USESUMMITS=$3
EXPNUMBER=$4
WORDLENGTH=$5
SIZE=$6
PVAL=$7
TFNAME=$8
UPLEN=$9
DOWNLEN=${10}

cd $SAMPLEDIR


# Generating [name]_peaks.narrowPeak and [name]_summits.bed for each replica

cd ../results

macs2 callpeak -t ../samples/chip/chip${IND}.bam -c ../samples/ctrl/ctrl${IND}.bam -f BAM --outdir . -n sample${IND}

echo "Peaks calling for sample${IND} completed"
echo "Results saved in results folder"



## Write in blackboard


echo "replica ${IND} is finished" >> ../logs/secondbb.txt

## Read Blackboard

FINISHEDREPS=$( wc -l ../logs/secondbb.txt | awk '{ print $1 }' )

IND=1
NEXT=$(($IND +1))

if [ $FINISHEDREPS = $EXPNUMBER ]
then
	# Merging replica data

	if [ $USESUMMITS = TRUE ]
	then

		

	        if [ $EXPNUMBER -gt 1 ]
        	then
                bedtools intersect -a sample${IND}_summits.bed -b sample${NEXT}_summits.bed > merged_beds.bed
		
			if [ $EXPNUMBER -gt 2 ]
			then

        	        	for OND in `seq 2 $EXPNUMBER`
	                	do
          	        		bedtools intersect -a sample${OND}_summits.bed -b merged_beds.bed > merged_beds2.bed
					rm merged_beds.bed
					mv merged_beds2.bed merged_beds.bed
				
		                done
			fi
		else

		mv sample${IND}_summits.bed merged_beds.bed

		fi

		DNABEDS="../results/merged_beds.bed"

		## Finding enriched words in the DNA peaks obtained
		findMotifsGenome.pl merged_beds.bed ../genome/genome.fa enriched_words -len $WORDLENGTH -size $SIZE

        	

	elif [ $USESUMMITS = FALSE ]
	then


	        if [ $EXPNUMBER -gt 1 ]
        	then

		 bedtools intersect -a sample${IND}_peaks.narrowPeak -b sample${NEXT}_peaks.narrowPeak > merged_peaks.narrowPeak
            		if [ $EXPNUMBER -gt 2 ]
			then

			    	for OND in `seq 2 $EXPNUMBER`
	                	do
        	                	bedtools intersect -a sample${OND}_peaks.narrowPeak -b merged_peaks.narrowPeak > merged_peaks2.narrowPeak
					rm merged_peaks.narrowPeak
					mv merged_peaks2.narrowPeak merged_peaks.narrowPeak

                		done
			fi
		else
		mv sample${IND}_peaks.narrowPeak merged_peaks.narrowPeak

		fi


		DNABEDS="../results/merged_peaks.narrowPeak"

		## Finding enriched words in the DNA peaks obtained
                findMotifsGenome.pl merged_peaks.narrowPeak ../genome/genome.fa enriched_words -len $WORDLENGTH -size $SIZE

	

	fi


	echo ""
	echo "Either summits.bed or peaks.narrowPeak merged succesfully"
	echo "Results found in results folder as merged_beds.bed or merged_narrow_narrowPeak"
	echo ""

	## R script call for regulome analysis

	Rscript  ../scripts/peaks.R $DNABEDS $TFNAME $UPLEN $DOWNLEN $PVAL


fi
