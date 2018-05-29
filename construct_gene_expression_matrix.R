#########################################################################################################
#	> File Name: construct_gene_expression_matrix.R
#	> This program constructs the RNA-seq gene expression matrix from htseq_count result files 
#   > and cufflinks result files, and normalize the obtained RNA-seq expression matrix
#	> Author: Hua Yu
#	> Heavily Edited by James Taylor
#	> Mail: huayu@genetics.ac.cn 
#	Created Time: 2014-08-10
#########################################################################################################



# normalizing RNA-seq expression data
normalizationData <- function (path,suffix,pathadd){
	if(suffix == "htseq_count.txt"){
		rawCount <- as.matrix(read.delim(file=paste(path,"htseq_count.txt",sep="/"),header=TRUE,row.names="GENE_ID",stringsAsFactors=FALSE))
		condition <- read.csv(file = "phenodata.csv",header = TRUE,row.names = 1)
		if(all(rownames(condition) %in% colnames(rawCount)) == FALSE){
		  return(10)
		}
		countdata <- rawCount[,rownames(condition)]
		if(all(rownames(condition) == colnames(countdata)) == FALSE){
		  return(11)
		  
		}
		
		DESeqCountData <- DESeqDataSetFromMatrix(countData = countdata,colData = condition, design = ~ treatment)
		DESeqCountData <- estimateSizeFactors(DESeqCountData)
		cat("The normalization factors for the different sample runs are",sizeFactors(DESeqCountData))
		normalizedMedianData <- counts(DESeqCountData,normalized = TRUE)
		write.table(normalizedMedianData,file=paste(path,"Median_htseq_DESeq_Normalized_Count.txt",sep="/"),row.names = TRUE,col.names = NA,quote=FALSE)
		DispersionsBlind <- estimateDispersions(DESeqCountData)
		normalizedVSTData <- getVarianceStabilizedData(DispersionsBlind)
		write.table(normalizedVSTData,file=paste(path,"VST_htseq_DESeq_Normalized_Count.txt",sep="/"),row.names = TRUE,col.names = NA,quote=FALSE)
		edgeRobjectTMM <- DGEList(counts=rawCount,genes=rownames(rawCount))
		edgeRobjectTMM <- calcNormFactors(edgeRobjectTMM)
		for(colIndex in seq(1,ncol(edgeRobjectTMM$counts),by=1)){
			edgeRobjectTMM$counts[,colIndex]=(edgeRobjectTMM$counts[,colIndex])/(edgeRobjectTMM$samples$norm.factors[colIndex])
		}
		write.table(edgeRobjectTMM$counts,file=paste(path,"TMM_htseq_EdgeR_Normalized_Count.txt",sep="/"),row.names = TRUE,col.names = NA,quote=FALSE)
		edgeRobjectUQ <- DGEList(counts=rawCount,genes=rownames(rawCount))
		edgeRobjectUQ <- calcNormFactors(edgeRobjectUQ,method="upperquartile",p=0.9)
		for(colIndex in seq(1,ncol(edgeRobjectUQ$counts),by=1)){
			edgeRobjectUQ$counts[,colIndex]=(edgeRobjectUQ$counts[,colIndex])/(edgeRobjectUQ$samples$norm.factors[colIndex])
		}
		write.table(edgeRobjectUQ$counts,file=paste(path,"UQ_htseq_EdgeR_Normalized_Count.txt",sep="/"),row.names = TRUE,col.names = NA,quote=FALSE)
	}
	if(suffix == "genes.fpkm_tracking"){
		rawCount <- read.delim(file=paste(pathadd,"htseqGenerawCountForSample.txt",sep="/"),header=TRUE,row.names=1,stringsAsFactors=FALSE)
		FPKMCount <- read.delim(file=paste(path,"cufflinksGeneFPKMCountforSample.txt",sep="/"),header=TRUE,row.names=1,stringsAsFactor=FALSE)
		condition <- colnames(rawCount)
		for(i in seq(1,length(condition),by = 1)){
			if((length(which(rawCount[,names(rawCount)==condition[i]] == 0))/length(rawCount[,names(rawCount)==condition[i]])) > 0.9){
				rawCount[,names(rawCount)==condition[i]] <- NULL
			}
		}
		geneNames <- rownames(rawCount)
		delMarks <- c()
		for(i in seq(1,length(geneNames),by = 1)){
			if((length(which((rawCount[i,]) < 10))/length(rawCount[i,])) > 0.8){
				delMarks <- append(delMarks,i)
			}else if(sd(as.numeric(rawCount[i,]),na.rm=TRUE)/mean(as.numeric(rawCount[i,]),na.rm=TRUE) <= 0.5){
				delMarks <- append(delMarks,i)
			}
		}
		if(!is.null(delMarks)){
			rawCount <- rawCount[-delMarks,]
		}
		commonGenes <- intersect(rownames(rawCount),rownames(FPKMCount))
		commonSamples <- intersect(colnames(rawCount),colnames(FPKMCount))
		FPKMCount <- FPKMCount[commonGenes,commonSamples]
		write.table(FPKMCount,file=paste(path,"FPKM_Cufflinks_Normalized_Count.txt",sep="/"),row.names = TRUE,col.names = NA,quote=FALSE)
	}
}

# loading necessary R libraries
library(DESeq2)
library(edgeR)

### invoking functions to conduct below tasks
normalizationData("C:/Users/jimta/Desktop","htseq_count.txt")
