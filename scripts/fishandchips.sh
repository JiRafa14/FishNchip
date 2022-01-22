#! /bin/bash

## Help message

if [ $# -ne 1 ]
then
        echo ""
        echo "usage: chipsnfish <param-file> "
        echo ""
        echo "param-file: file where the parametres are specified"
        echo "see README for further information or take given param.txt as an example"
        echo ""
        exit

fi
echo ""
echo "Thaks for using  FISHANDCHIPS"
echo ""
echo ""

echo "--READING PARAMETERS--"

## Param reading

PARAM=$1


WD=$(grep working_directory: $PARAM | awk '{ print $2 }')
echo "Working directory = $WD"
FD=$(grep folder_name: $PARAM | awk '{ print $2 }')
echo "Folder name = $FD"
CHIPNUM=$(grep number_chip: $PARAM | awk  '{ print $2 }')
echo "Number of CHiP replicas = $CHIPNUM"
CTRLNUM=$(grep number_ctrl: $PARAM | awk  '{ print $2 }')
echo "Number of control replicas = $CTRLNUM"


GEN=$(grep path_to_genome: $PARAM | awk '{ print $2 }')
echo "Reference genome = $GEN"
ANN=$(grep path_to_annotation: $PARAM | awk '{ print $2 }')
echo "Genome Annotation = $ANN"
PAIRED=$(grep is_paired: $PARAM | awk '{ print $2 }')
echo "Paired files = $PAIRED"

SCRIPTS=$(grep path_to_scripts: $PARAM | awk '{ print $2 }')
echo "Scripts folder = $SCRIPTS"
SAMPLEFORMAT=$(grep sample_format: $PARAM | awk '{ print $2 }')
echo "Sample format = $SAMPLEFORMAT"
EXPNUMBER=$(grep number_of_experiments: $PARAM | awk '{ print $2 }')
echo "Number of replicas = $EXPNUMBER"
USESUMMITS=$(grep use_summits: $PARAM | awk '{ print $2 }')
echo "Use summits = $USESUMMITS"
WORDLENGTH=$(grep word_length: $PARAM | awk '{ print $2 }')
echo "Word length = $WORDLENGTH"
SIZE=$(grep size: $PARAM | awk '{ print $2 }')
echo "Size for word search: $SIZE"
PVAL=$(grep pvalue_cutoff: $PARAM | awk '{ print $2 }')
echo "Pvalue cutoff for GO enrichment: $PVAL"
TFNAME=$(grep transcription_factor_name: $PARAM | awk '{ print $2 }')
echo "Transcription factor name: $TFNAME"
UPLEN=$(grep upstream_promoter_length: $PARAM | awk '{ print $2 }')
echo "Upstream promoter length: $UPLEN"
DOWNLEN=$(grep downstream_promoter_length: $PARAM | awk '{ print $2 }')
echo "Downstream promoter length: $DOWNLEN"


echo ""
echo "Parameters have been loaded"

## Creating Workspace

cd $WD
mkdir $FD
cd $FD
mkdir genome annotation samples scripts results logs
cd samples
mkdir ctrl chip


if [ $SAMPLEFORMAT == "FILE" ]
then
	if [ $PAIRED = FALSE ]
	then
		IND=1
		while [ $IND -le $CHIPNUM ]
		do

        		echo ""
		        echo "copying CHiP $IND"
	
			cd chip
	
		        ROUTE=$( grep "chip${IND}:" $PARAM | awk '{ print $2 }' )
        		cp $ROUTE chip${IND}.fastq.gz
	
	
			IND=$(($IND + 1))
			cd ../

		done
	
		IND=1
		while [ $IND -le $CTRLNUM ]
		do
			echo ""
			echo "copying control $IND"
	
			cd ctrl
	
        		ROUTE=$( grep "ctrl${IND}:" $PARAM | awk '{ print $2 }' )
	        	cp $ROUTE ctrl${IND}.fastq.gz

		        IND=$(($IND + 1))
		        cd ../


		done
	
	elif [ $PAIRED = TRUE ]
	then
		IND=1
		while [ $IND -le $CHIPNUM ]
		do
			echo ""
			echo "copying CHiP $IND "
	
			cd chip
	
			ROUTE=$( grep "chip${IND}_A:" $PARAM | awk '{ print $2 }' )
	                cp $ROUTE chip${IND}_A.fastq.gz
	
        	        ROUTE=$( grep "chip${IND}_B_:" $PARAM | awk '{ print $2 }' )
                	cp $ROUTE chip${IND}_B.fastq.gz
	
	
			IND=$(($IND + 1))
                	cd ..
	
	
	
		done
	
		IND=1
		while [ $IND -le $CTRLNUM ]
		do
			echo ""
			echo "copying control $IND"
	
			cd ctrl
                	ROUTE=$( grep "ctrl${IND}_A:" $PARAM | awk '{ print $2 }' )
	                cp $ROUTE ctrl${IND}_A.fastq.gz
	
        	        ROUTE=$( grep "ctrl${IND}_B:" $PARAM | awk '{ print $2 }' )
                	cp $ROUTE ctrl${IND}_B.fastq.gz
	
        	        IND=$(($IND + 1))
                	cd ../
		done
	fi

elif [ $SAMPLEFORMAT == "LINK" ]
then

        if [ $PAIRED = FALSE ]
	then
                IND=1
                while [ $IND -le $CHIPNUM ]
                do

                        echo ""
                        echo "copying CHiP $IND"

                        cd chip

                        ROUTE=$( grep "chip${IND}:" $PARAM | awk '{ print $2 }' )
                        fastq-dump --gzip --split-files $ROUTE
			mv ${ROUTE}_1.fastq.gz chip${IND}.fastq.gz

                        IND=$(($IND + 1))
                        cd ../

                done

                IND=1
                while [ $IND -le $CTRLNUM ]
                do
                        echo ""
                        echo "copying control $IND"

                        cd ctrl

                        ROUTE=$( grep "ctrl${IND}:" $PARAM | awk '{ print $2 }' )
			fastq-dump --gzip --split-files ${ROUTE}
			mv ${ROUTE}_1.fastq.gz ctrl${IND}.fastq.gz

                        IND=$(($IND + 1))
                        cd ../


                done
        elif [ $PAIRED = TRUE ]
	then

                IND=1
                while [ $IND -le $CHIPNUM ]
                do
                        echo ""
                        echo "copying CHiP $IND "

                        cd chip

                        ROUTE=$( grep "chip${IND}:" $PARAM | awk '{ print $2 }' )
                        fastq-dump --gzip --split-files ${ROUTE}


                        mv ${ROUTE}_1.fastq.gz chip${IND}_A.fastq.gz
			mv ${ROUTE}_2.fastq.gz chip${IND}_B.fastq.gz

                        IND=$(($IND + 1))
                        cd ..



                done

                IND=1
                while [ $IND -le $CTRLNUM ]
                do
                        echo ""
                        echo "copying control $IND"

                        cd ctrl
                        ROUTE=$( grep "ctrl${IND}:" $PARAM | awk '{ print $2 }' )
                        fastq-dump --gzip --split-files ${ROUTE}


                        mv ${ROUTE}_1.fastq.gz ctrl${IND}_A.fastq.gz
                        mv ${ROUTE}_2.fastq.gz ctrl${IND}_B.fastq.gz

			cd ..
                done
        fi

fi

## Placing scripts in scripts folder

cd ../scripts

cp $SCRIPTS/FnCsample_processing.sh .
cp $SCRIPTS/peaks.R .
cp $SCRIPTS/peakAnalysis.sh .

cd ../genome
cp $GEN genome.fa.gz
gunzip genome.fa.gz


cd ../annotation
cp $ANN annotation.gtf.gz
gunzip annotation.gtf.gz

cd ../

## Creating index

cd genome
bowtie2-build genome.fa index


cd ../



TOTALSAMP=$(($CHIPNUM + $CTRLNUM))
IND=1
CTRL=FALSE
while [ $IND -le $CHIPNUM ]
do
        cd logs

        sbatch ../scripts/FnCsample_processing.sh ${WD}/${FD}/samples $IND $PAIRED $CTRL $TOTALSAMP $EXPNUMBER $USESUMMITS $WORDLENGTH $SIZE $PVAL $TFNAME $UPLEN $DOWNLEN

	echo "CHiP $IND submitted to processing"

        IND=$(($IND + 1))
	cd ..

done

IND=1
CTRL=TRUE
while [ $IND -le $CTRLNUM ]
do
	cd logs

	sbatch ../scripts/FnCsample_processing.sh ${WD}/${FD}/samples $IND $PAIRED $CTRL $TOTALSAMP $EXPNUMBER $USESUMMITS $WORDLENGTH $SIZE $PVAL $TFNAME $UPLEN $DOWNLEN

	echo "Control $IND submitted to processing"

	IND=$(($IND + 1))
	cd ..
	
done
