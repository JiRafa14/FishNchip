## Script created to obtain the regulome. 
## Authors: Rodrigo Bedera Garcia
##          Rafael Morales Marquez

## Needed inputs : 
##  DNA beds file (.bed or .narrowPeak)
##  Name of TF analyzed
##  Upstream promoter length
##  Downstream promoter length
##  Pvalue cutoff for GO enrichment

arguments <- commandArgs(trailingOnly = TRUE)
print(arguments)

dnabeds <- arguments[1] #Path to dnabeds file


tfname <- as.character(arguments[2]) #Name of TF
upprom <- as.numeric(arguments[3]) #Promoter length upstream
downprom <- as.numeric(arguments[4]) #Promoter length downstream
pval <- as.numeric(arguments[5])

print("..")
print(c(dnabeds,tfname,upprom,downprom,pval))

## Install needed packages

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("ChIPseeker")
library(ChIPseeker)

#BiocManager::install("DO.db") 
#BiocManager::install("GO.db") 

library("DO.db")
library("GO.db")


#BiocManager::install("TxDb.Athaliana.BioMart.plantsmart28")
library(TxDb.Athaliana.BioMart.plantsmart28)


txdb <- TxDb.Athaliana.BioMart.plantsmart28

as.data.frame(head(genes(txdb)))

## Read dna peaks file

peaks <- readPeakFile(peakfile = dnabeds, header=FALSE)

## Obtaining promoters

print("Obtaining promoters")

promoter <- getPromoters(TxDb=txdb, 
                         upstream=upprom, 
                         downstream=downprom)

peakAnno <- annotatePeak(peak = peaks, 
                              tssRegion=c(- upprom,downprom),
                              TxDb=txdb)

plotAnnoPie(peakAnno,main = paste(c(tfname," binding sites"),collapse = ""))

plotDistToTSS(peakAnno,
              title="Distribution of genomic loci relative to TSS",
              ylab = "Genomic Loci (%) (5' -> 3')")

## Writing out results
print("Writing out results")

annot <- as.data.frame(peakAnno)

target.printing <- function(a,b)
{
	if (a >= 2000 | b >= 2000)
	{
		annot$geneId[annot$annotation == "Promoter (<=1kb)" | annot$annotation == "Promoter (1-2kb)" |
                                         annot$annotation == "Promoter (2-3kb)"]	
	}

	else if (a < 2000 & b < 2000)
	{
		 annot$geneId[annot$annotation == "Promoter"]
	}
}

target.genes <- target.printing(upprom,downprom)

write(x = target.genes,file = paste(c(tfname,"_target_genes.txt"),collapse = ""))

## GO terms enrichment

#BiocManager::install("clusterProfiler")

library("clusterProfiler")

#BiocManager::install("org.At.tair.db")
library(org.At.tair.db)

goenrichment <- clusterProfiler::enrichGO(gene = target.genes,OrgDb = org.At.tair.db,ont = "ALL",keyType = "TAIR", pvalueCutoff = pval)
head(goenrichment)

goenrichment <- as.data.frame(goenrichment)
write.csv(x = goenrichment,file = "GO_enrichment.csv")

promoter <- as.data.frame(promoter)
write.csv(x = promoter,file = "promoters.csv")
