#set up envrionment
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("alyssafrazee/RSkittleBrewer","ballgown","genefilter","dplyr","devtools"))

library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)

#set working directory to wherever your data is
#setwd(C:/Users/Jim/Desktop)

#load phenotype data file. The second line avoids an error saying the names are in the wrong order. 
pheno_data = read.csv("phenodata.csv")
pheno_data = pheno_data[order(pheno_data$id),]

#create ballgown object with data from stringtie

bg = ballgown(dataDir = "ballgown", samplePattern = "W", pData= pheno_data)

#filter low abundance genes
bg_filt = subset(bg, "rowVars(texpr(bg)) >1",genomesubset = TRUE)

#statistical tests
results_transcripts = stattest(bg_filt, feature="transcript",covariate = "time",getFC = TRUE, meas = "FPKM")
results_genes = stattest(bg_filt, feature="gene", covariate = "time", getFC = TRUE, meas = "FPKM")

#add gene names and id to results transcripts
results_transcripts = data.frame(transcriptNames=ballgown::transcriptNames(bg_filt),geneNames=ballgown::geneNames(bg_filt), geneIDs=ballgown::geneIDs(bg_filt),results_transcripts)

results_transcripts = arrange(results_transcripts,pval)
results_genes = arrange(results_genes,pval)

#get list of jgi gene names from ensembl biomart and import here
#remove the hyperlinks it provides and use concat to add "transcript:" before every EHK##
#also I changed the title of the transcript column to transcriptNames to make life easier and so it matched the ballgown jargon
#mart_export <-read_excel("/path/to/excel/file.xls")

results_transcripts$JGI <- mart_export$GeneID[match(results_transcripts$transcriptNames,mart_export$transcriptNames)]

write.csv(results_transcripts,"results_transcripts-time.csv",row.names=FALSE)
write.csv(results_genes,"results_genes-time.csv",row.names = FALSE)

fpkm= texpr(bg,meas = "FPKM")
fpkm = log2(fpkm+1)

# plot fpkm
#plot(fpkm[id,]~pheno_data$variable, main=ballgown::transcriptNames(bg)[id],pch=19, xlab=variable, ylab = 'log2(FPKM+1)')

