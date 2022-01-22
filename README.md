# FISHANDCHIPS: An new automatic pipeline to process and analyse ChIP-seq data of Arabidopsis Thaliana

Authors: 

Rodrigo Bedera Garcia and Rafael Morales Marquez

University of Seville

## 1. Summary

## 2. Installation

## 3. Dependencies

## 4. Input

## 5. Usage

### a. Parameters read

### b. Workspace generation

### c. Index building

### d. Sample load

### e. Sample processing

#### - Quality control

#### - Reads alignment

### f. Peaks calling


## 6. Example


## 1. Summary 

This repository is made up of 3 bash scripts (fishandchips, FnC_sample_processing and peakAnalysis) and one R script (peaks.R), to analyse unlimited ChIP-seq samples of transcription factors from Arabidopsis thaliana, comparing all of them with the same control samples. Actually, by running the main script, fishandchips, the whole analysis is done, since this script is programmed to run the other three scripts when required.

This way, with only one command, the whole ChIP-seq samples analysis is done and organized in intuitive directories. Briefly, this analysis consists of the quality control of the samples, the reads mapping, the peaks calling, the regulome determination and Gene Set Enrichment Analysis.


## 2. Installation

To install the FISHANDCHIPS software, follow these steps:

1. Download scripts folder and parameters_file.txt (located at example)

2. Place the scripts folder and parameters_file.txt wherever is desired

3. In the parameteres_file.txt, write the path_to_scripts path. Final dash ("/") must not be writtend at the end.

4. Add the scripts folder to the bash.rc file. This enables fishandchips.sh script to be launched from any folder.

5. Installation done!


## 3. Dependencies

Tools needed to run the bash scripts: bowtie2 (https://howtoinstall.co/es/bowtie2), fastqc (https://howtoinstall.co/es/fastqc), samtools (https://howtoinstall.co/es/samtools) , bedtools (https://bedtools.readthedocs.io/en/latest/content/installation.html) , SRA-toolkit (https://hpc.nih.gov/apps/sratoolkit.html) ,  macs2 (https://command-not-found.com/macs2) and homer (http://homer.ucsd.edu/homer/introduction/install.html). These tools can be installed running the following command: sudo apt-get install <tool_name>

Packages needed to run the R script: ChipSeeker (https://www.bioconductor.org/packages/release/bioc/html/ChIPseeker.html), DO.db (https://www.bioconductor.org/packages/release/data/annotation/html/DO.db.html), GO.db (https://www.bioconductor.org/packages/release/data/annotation/html/GO.db.html), clusterProfiler (https://www.bioconductor.org/packages/release/bioc/html/clusterProfiler.html),  TxDb.Athaliana.Biomart.plantsmart28 (https://www.bioconductor.org/packages/release/data/annotation/html/TxDb.Athaliana.BioMart.plantsmart28.html) and org.At.tair.db (https://www.bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html). All packages can be downloaded from Bioconductor, running the following command in R: BiocManager::install(“<package_name>”). 


## 4. Input

The pipeline FISHANDCHIPS requires only one parameter. This parameter, the so-called “parameters file” is a TXT file with the parameters needed by the pipeline. To know more about the file content, go to section 5.a. 

## 5. Usage

To run the pipeline, write on the command line (having previously written the custom parameters):

fishandchips.sh /full/path/to/parameters_file.txt


### a. Parameters read 

The parameters passed to the function need to be written on a TXT file, the “parameters file”. An example of a parameters file can be found on the directory test. When the pipeline is run with this file, a case study is generated. To see the example, go to section 6.

The parameters file includes the following information:


-working_directory: The path to the directory where will be generated the experiment directory.

-folder_name: The name of the experiment directory, which will be located in the working directory. Every file created by the pipeline will be saved in a subdirectory of the experiment directory.

-number_chip: Number of samples for every transcription factor sample. 

-number_ctrl: Number of samples for every control sample.

-path_to_genome: The path to the file to be used as reference genome, in FASTA format. It must be in .fa.gz format. .gff3.gz was tested and seems to work well too.

-path_to_annotation: The path to the file to be used as reference genome annotation, in GTF format. It must be in gtf.gz format.

-paired: Write FALSE for single-end sequenced samples and TRUE for paired-end sequenced samples.

-path_to_scripts: The path to every scripts that the pipeline uses. The final slash of the path mustn´t be written, as shown in parameters_file.txt

-sample_format: The format in which the samples used are found fot the pipeline. Write FILE for file format or LINK for download link format. If a link is provided, it must be compatible with the fastq-dump function.

-number_of_experiments: Number of replicas that have been obtained. This pipeline is designed to analyse all the samples provided, returning a quality control and .bam , .bam.bai files. However, when performing peaks calling with macs2, a pair of control and CHiP data must be introduced. This way, if 3 CHiP samples and 4 control samples are obtained, number_of_experiments: 3 (3 full pairs).

-use_summits: Depending on the restriction level desired, use summit (TRUE) (more restrictive) or narrowpeaks (FALSE) (less restrictive).

-word_length: Word length enrichment searched by homer. Default is 8,10,12

-size: Size of the peak in nucleotides where homer searches for word enrichment. Default is 200

-pvalue_cutoff: Threshold used to study the enrichment of gene ontology terms. By default, it has a value of 0.05.

-transcription_factor_name: Name of the transcription factor in the study.

-upstream_promoter_length: Upstream length of the promoter from the TSS. By default, it has a value of 1000.


-downstream_promoter_length: Downstream length of the promoter from the TSS. By default, it has a value of 1000.


When a parameters file is passed to the pipeline FISHANDCHIPS, the pipeline stores the lines of the parameters file in different variables, and prints out on the screen every variable for each parameter.


### b. Workspace generation

To create the workspace, the route of the working directory is run. Once there, the experiment directory is created with the name given on the parameter file, then accessed, and there, five directories are created: genome, annotation, samples, results and scripts.

After that, using their paths, written on the parameters file, the genome and annotation files are copied to the directories genome and annotation, and extracted under the names “genome.fa” and “annotation.gtf”, respectively. 

### c. Index building

To build the index, the pipeline leaves the annotation directory and accesses the genome directory.

Using the function bowtie2-build, which receives the genome.fa and a word to use as a suffix of the index files, the pipeline builds an index of the reference genome, and saves it on the genome directory, under the suffix .index.

### d. Sample load 

We have 4 possible combinations of paired and sample format. Paired can be true or false. In addition, the sample format can be file or link.
Normally, they will be read as sample1, sample2. However, when samples are paired, and the sample format is FILE, 2 files are needed. These are read as sample1_A, sample1_B, ctrl1_A, ctrl1_B, etc.
When the samples are in link format and paired, only the download link is needed, as both files are created in the process.

### e. Sample processing 

FISHANDCHIPS throws FnCsample_proccesing.sh for each chip and control sample. As mentioned above the number of total samples can be odd. All these samples are analysed.
However, in the following steps of the pipeline we will stick with the paired samples.
The processing of the samples is done, according to whether the samples are paired or whether they are chip or control samples.
Depending on the format of the samples (link or file) a different reading is performed. As the samples are being processed, the blackboard is filled in. When this step is finished, the next script, peakAnalysis.sh, is launched.


#### - Quality control

Once the pipeline runs the sample_processing bash script, a quality control is made for the sample passed to the script, using the function fastqc, which receives the name of the sample (SAMPLE1) or samples -if paired-end, SAMPLE1_A and SAMPLE1_B- .


#### - Reads alignment

The next step is to align the reads to the index of the reference genome, previously built. This way, the different zones of the genome where reads accumulate, which probably correspond to binding regions of the transcription factor, can be identified. To make the alignment, the bowtie2 function is used. The parameters passed to this function are the following: 

-x    /path/to/index_of_reference_genome

-U or -1 -2   , depending if the samples are single or paired-end, respectively. 

-S   output name of the SAM file generated. 


After the SAM file is created, since it is a very heavy file, it is converted to a BAM file under the same name. To make this conversion, the function samtools sort is used. The parameters passed to this function are the output BAM file name (SAMPLENAME.bam) and the SAM file. Then, the SAM file is removed.


### f. Peaks calling

Even though the reads are aligned on the genome, the user cannot know, only with that information, which are the possible target genes for the transcription factors analysed. To check if the result of the mapping is significant, a statistical analysis is made. Using a contrast of hypothesis, the error probability can be determined. 

To call peaks, a sliding window goes through the genome, counting reads on the window defined, for chips and control samples. Then, the reads of each pair of samples are compared, studying the fold-change and p-value. These values represent the difference in the number of reads for the chip and control sample, and how significant that difference is, respectively. It is important to know that, since this process is a multiple testing, a correction of the p-value, called q-value or FDR, needs to be made. Therefore, peaks are regions in which there is a significant difference of reads between the chip and the control samples.

To call peaks, the function macs2 with the command callpeak is used. The arguments passed to this function are the following:

-t   path to the BAM files for the replicas of the chip sample.

-c   path to the BAM files for the replicas of the control sample.

-f   BAM format, which is the format of the input files.

-n   name of the output file that are being compared.

-- outdir to specify the output directory, which is still the results directory.

PeakAnalysis.sh is launched for each pair. For example, if we have 3 chip samples and 4 control samples, peakAnalysis.sh would run 3 times. The sample that is left loose does not continue in the analysis. As the samples are being processed, the second blackboard is filled in. When this step is finished, it is checked if the number of experiments is greater than 1. Depending on how comprehensive the study wants to be, the .bed or .NarrowPeak file is generated.
On the contrary, if we have more than one experiment, a loop is generated in which the bed or NarrowPeak files generated for each experiment are compared, generating a merged_bed or merged_NarrowPeak file.


### g. Regulome determination and Gene Set Enrichment Analysis


The first thing that the R script, called peaks_script.R does, is to store the arguments passed on the bash script on variables. Then, the script loads the necessary packages for the analysis. To know the packages needed, go to section 3. 

Then, the script gets the genetic information from the organism, in this case Arabidopsis thaliana, and gets the universe, in this case, the genes id of the whole organism. After this, the script reads the different files.

The next step is to define the promoters for narrow peaks and summits, using the base pair up and downstream given in the parameters file, and to annotate those peaks.

Then, the annotation is converted to a data frame. Now, from each data frame, the regulome is determined, simply by extracting the elements whose annotation is “Promoter”. The genes ID of those elements is the regulome, which the pipeline saves on a TXT file called <transcription factor name>_target_genes.txt

Once the regulome is determined, the next step is the Gene Set Enrichment Analysis (GSEA) using GO terms.

The enrichGO function receives the regulome, the universe defined, an OrgDb object -in this case, that for Arabidopsis thaliana-, the ontology terms to analyse -all, in this case-, the p-value cutoff and the keyType -in this case, TAIR-. 

After this, the analysys of the peaks is over. 2 plots are generated in the Rplots.pdf file: the transcription factor binding sites and the distribution of genomic loci relative to the TSS. Also, the gene set enrichment analysis results is saved in GO_enrichment.csv, and promoters in promoters.csv. As said above, the regulome is saved in <transcription factor name>_target_genes.txt.


## 6. Example

We use an study about the role Arabidopsis thaliana homeobox gene 1 in the control of the plant architecture (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8092594/). In this study the scientists do a Chip-seq (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157332). The pipeline can use this samples and the parameters_file.txt has the default values set according to this study . 

For this study, we used the reference genome and the annotation file provided in the data folder.
