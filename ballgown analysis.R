#set up envrionment
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("alyssafrazee/RSkittleBrewer","ballgown","genefilter","dplyr","devtools"))

library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)
library(readxl)
library(tidyr)
library(data.table)
library(ggplot2)

#set working directory to wherever your data is
setwd("C:/Users/jimta/Desktop")

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
full_transcripts = stattest(bg,feature = "transcript",covariate = "treatment",meas = "FPKM")
#add gene names and id to results transcripts
results_transcripts = data.frame(transcriptNames=ballgown::transcriptNames(bg_filt),geneNames=ballgown::geneNames(bg_filt), geneIDs=ballgown::geneIDs(bg_filt),results_transcripts)
full_transcripts = data.frame(transcriptNames=ballgown::transcriptNames(bg),geneNames=ballgown::geneNames(bg), geneIDs=ballgown::geneIDs(bg),full_transcripts)

results_transcripts = arrange(results_transcripts,pval)
results_genes = arrange(results_genes,pval)

#get list of jgi gene names from ensembl biomart and import here
#remove the hyperlinks it provides and use concat to add "transcript:" before every EHK##
#also I changed the title of the transcript column to transcriptNames to make life easier and so it matched the ballgown jargon
mart_export <-read_excel("C:/Users/jimta/Downloads/mart_export.xls")

results_transcripts$JGI <- mart_export$GeneID[match(results_transcripts$transcriptNames,mart_export$transcriptNames)]

full_transcripts$JGI <- mart_export$GeneID[match(full_transcripts$transcriptNames,mart_export$transcriptNames)]

#import signalp data from jgi
signalp <- read.delim("C:/Users/jimta/Downloads/Tvirens_v2.signalp_FrozenGeneCatalog_20100318.tab")

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

#output files
#write.csv(results_transcripts,"total_results_transcripts.csv",row.names=FALSE)
#write.csv(secreted,"secreted_results_transcripts.csv",row.names = FALSE)

#generate dataset over time for graphing

#new pheno_data to select only B73+WT at 6, 30, 54 hours
pheno = read.csv("pheno.csv")
pheno = pheno[order(pheno$id),]

bg_time = ballgown(dataDir = "time", samplePattern = "W", pData= pheno)

fpkm= texpr(bg_time,meas = "FPKM")
fpkm = log2(fpkm+1)

id <- secreted$id

#creating and rearranging tables for generation of plots
selected <- as.data.frame(fpkm[id,])

cnames <- as.data.frame(fpkm[id,])

cnames <- as.data.frame(colnames(cnames))

cnames <- separate(cnames,"colnames(cnames)",into = c("fpkm","index"),sep = "FPKM.")

cnum <- length(cnames$index)

c <- data.frame(matrix(ncol = cnum, nrow = 0))

names(c) <- cnames$index

names(selected) <- cnames$index

c[1,]<-pheno$time

comb <- rbind(selected,c)

comb <- t(comb)

names(comb) <- c("a","b","c","d","e","f","g","h","i","j","l","time")

comb <- as.data.frame(comb)


# # plot fpkm
#ggplot(comb) + geom_point(aes(x=comb$time,y=comb$a))

grex <- function(jgi){
  
  cgi = class(jgi)
  
  if (cgi == "list" | cgi == "character") {
    
    for (a in jgi) {
      a=paste("TRIVIDRAFT_",a,sep = "")
      
      b <- full_transcripts$id[match(a,full_transcripts$JGI)]
      
      jpeg(paste(a,".jpeg",sep = ""))
      plot(fpkm[b,]~pheno$time, main = a, xlab="Time",ylab=paste("FPKM",a))
      dev.off()
      
    }
    
  }
  
  
}
