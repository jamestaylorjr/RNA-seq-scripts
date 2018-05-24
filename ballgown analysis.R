#set up envrionment
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("alyssafrazee/RSkittleBrewer","ballgown","genefilter","dplyr","devtools"))

library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)
library(readxl)

#set working directory to wherever your data is
#setwd(C:/Users/Jim/Desktop)

#load phenotype data file. The second line avoids an error saying the names are in the wrong order. 
pheno_data = read.csv("6hr/pheno_data.csv")
pheno_data = pheno_data[order(pheno_data$id),]

#create ballgown object with data from stringtie

bg = ballgown(dataDir = "6hr", samplePattern = "W", pData= pheno_data)

#filter low abundance genes
bg_filt = subset(bg, "rowVars(texpr(bg)) >1",genomesubset = TRUE)

#statistical tests
results_transcripts = stattest(bg_filt, feature="transcript",covariate = "treatment",getFC = TRUE, meas = "FPKM")
results_genes = stattest(bg_filt, feature="gene", covariate = "treatment", getFC = TRUE, meas = "FPKM")

#add gene names and id to results transcripts
results_transcripts = data.frame(transcriptNames=ballgown::transcriptNames(bg_filt),geneNames=ballgown::geneNames(bg_filt), geneIDs=ballgown::geneIDs(bg_filt),results_transcripts)

results_transcripts = arrange(results_transcripts,pval)
results_genes = arrange(results_genes,pval)

#get list of jgi gene names from ensembl biomart and import here
#remove the hyperlinks it provides and use concat to add "transcript:" before every EHK##
#also I changed the title of the transcript column to transcriptNames to make life easier and so it matched the ballgown jargon
#mart_export <-read_excel("/path/to/excel/file.xls")

results_transcripts$JGI <- mart_export$GeneID[match(results_transcripts$transcriptNames,mart_export$transcriptNames)]

#import signalp data from jgi
signalp <- read.delim("path/to/file.tab")

#fix gene numbers so they match current format
signalp$proteinid <- paste("TRIVIDRAFT_",signalp$proteinid, sep="")

#convert FC to log2FC
results_transcripts$log2FC <- log2(results_transcripts$fc)

#filter genes below 1.5 log2FC
filtered_results <- results_transcripts[which(results_transcripts$log2FC < -1.5 & results_transcripts$qval < 0.05),]

#add signalp data to results file
filtered_results$signalp <- signalp$hmm_signalpep_probability[match(filtered_results$JGI,signalp$proteinid)]

#get only putatively secreted proteins
secreted = filtered_results[which(filtered_results$signalp > 0.5),]

write.csv(results_transcripts,"total_results_transcripts.csv",row.names=FALSE)
write.csv(secreted,"secreted_results_transcripts.csv",row.names = FALSE)


fpkm= texpr(bg,meas = "FPKM")
fpkm = log2(fpkm+1)

# plot fpkm
#plot(fpkm[id,]~pheno_data$variable, main=ballgown::transcriptNames(bg)[id],pch=19, xlab=variable, ylab = 'log2(FPKM+1)')

