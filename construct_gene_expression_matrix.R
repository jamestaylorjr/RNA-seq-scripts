#########################################################################################################
#	> File Name: construct_gene_expression_matrix.R
#	> This program constructs the RNA-seq gene expression matrix from htseq_count result files 
#   > and cufflinks result files, and normalize the obtained RNA-seq expression matrix
#	> Author: Hua Yu
#	> Mail: huayu@genetics.ac.cn 
#	Created Time: 2014-08-10
#########################################################################################################

# extracting the gene expression value from htseq-count result files or cufflinks result files 
extractGeneCount <- function(path,suffix){
	cat(path,"\n",suffix,"\n")
	if(suffix == "htseq_count.txt"){
		files <- file.path(path,suffix, fsep = .Platform$file.sep)
		expressionMatrix <- read.delim(file=files[1],header=FALSE,stringsAsFactors=FALSE, sep = "\t")
		files <- files[2:length(files)]
		uniqueIndex <- !duplicated(expressionMatrix[,1])
		expressionMatrix <- expressionMatrix[uniqueIndex,]
		expressionMatrix <- as.character(expressionMatrix[,1])
		sampleNames <- expressionMatrix
		columnNames <- c("GENE_ID")
		for (filename in files){
			if(file.info(filename)[1] > 0 & is.na(filename) == FALSE){
				matrixStep <- read.delim(file=filename,header=FALSE,stringsAsFactors=FALSE, sep = "\t")
				uniqueIndex <- !duplicated(matrixStep[,1])
				matrixStep <- matrixStep[uniqueIndex,]
				rownames(matrixStep) <- matrixStep[,1]
				if(class(matrixStep[,2]) == "character"){
					matrixStep[,2] = as.numeric(matrixStep[,2])
					matrixStep[which(is.na(matrixStep[,2])),2] <- 0
				}
				expressionMatrix <- as.matrix(cbind(expressionMatrix,matrixStep[sampleNames,2]))
				posIndex <- regexpr('[D|E|S]RR\\d+',filename)
				runid <- substring(filename,posIndex,posIndex+attr(posIndex,'match.length')-1)
				cat(runid,"\n")
				columnNames <- append(columnNames,runid)
			}
		}
		expressionMatrix <- as.data.frame(expressionMatrix)
		names(expressionMatrix) <- columnNames
		write.table(expressionMatrix,file=paste(path,"htseqGeneRawCount.txt",sep="/"), quote=FALSE, row.names=F, sep='\t')
	}
	if(suffix == "genes.fpkm_tracking"){
		files <- Sys.glob(file.path(path,"*","*",suffix))
		expressionMatrix <- read.delim(file=files[1],header=TRUE,stringsAsFactors=FALSE)
		uniqueIndex <- !duplicated(expressionMatrix[,1])
		expressionMatrix <- expressionMatrix[uniqueIndex,]
		expressionMatrix <- as.character(expressionMatrix[,1])
		sampleNames <- expressionMatrix
		columnNames <- c("GENE_ID")
		for (filename in files){
			htfilename <- filename
			htfilename <- sub("cufflinksfile","htseqfile",htfilename)
			htfilename <- sub("\\/genes.fpkm_tracking",".htseq_count.txt",htfilename)
			cat(htfilename,"\n")
			if(file.info(filename)[1] > 0){
				matrixStep <- read.delim(file=filename,header=TRUE,stringsAsFactors=FALSE)
				uniqueIndex <- !duplicated(matrixStep[,1])
				matrixStep <- matrixStep[uniqueIndex,]
				rownames(matrixStep) <- matrixStep[,1]
				matrixStep[,1] <- NULL
				matrixStep <- matrixStep[sampleNames,]
				if(class(matrixStep$FPKM) == "character"){
					matrixStep$FPKM = as.numeric(matrixStep$FPKM)
					matrixStep$FPKM[which(is.na(matrixStep$FPKM))] <- 0
				}
				expressionMatrix <- cbind(expressionMatrix,matrixStep$FPKM)
				posIndex <- regexpr('[D|E|S]RR\\d+',filename)
				runid <- substring(filename,posIndex,posIndex+attr(posIndex,'match.length')-1)
				cat(runid,"\n");
				columnNames <- append(columnNames,runid)
			}
		}
		expressionMatrix <- as.data.frame(expressionMatrix)
		names(expressionMatrix) <- columnNames
		write.table(expressionMatrix,file=paste(path,"cufflinksFPKM.txt",sep="/"),row.names=FALSE,quote=FALSE, sep='\t')
	}	
}

# combining expression value in each run for each gene according to the NCBI run information
combineGeneCount <- function(path,suffix,relatRunAndSample){
	if(suffix == "htseq_count.txt"){
		rawCount = read.delim(file=paste(path,"htseqGeneRawCount.txt",sep="/"),header=TRUE,stringsAsFactors=FALSE)
		rawCount <- rawCount[1:(nrow(rawCount)-5),]
		geneNames = rawCount$GENE_ID
		rawCount$GENE_ID <- NULL
		rawCountTrans <- as.data.frame(t(rawCount))
		colnames(rawCountTrans) <- geneNames
		runID <- rownames(rawCountTrans)
		matchPosition <- match(runID,relatRunAndSample[,1])
		RemainRunID <- runID[which(!is.na(matchPosition))]
		rawCountTrans <- rawCountTrans[RemainRunID,]
		matchPosition <- matchPosition[!is.na(matchPosition)]
		rawCountTrans[which(is.na(rawCountTrans),arr.ind=TRUE)] <- 0
		samples <- unique(relatRunAndSample[matchPosition,2])
		rawCountForSample <- samples
		for(geneName in geneNames){
			cat("Gene ID:",geneName,"\n")
			geneExpressionValue <- tapply(rawCountTrans[,match(geneName,colnames(rawCountTrans))],relatRunAndSample[matchPosition,2],sum)
			cat("The corresponding sample name:",names(geneExpressionValue),"\n")
			cat("The corresponding expression value:",as.numeric(geneExpressionValue),"\n")
			rawCountForSample <- cbind(rawCountForSample,as.numeric(geneExpressionValue[samples]))
		}
		rawCountForSample <- as.data.frame(rawCountForSample)
		rawCountForSample[,1] <- NULL
		rawCountForSample = t(rawCountForSample)
		colnames(rawCountForSample) <- samples
		rownames(rawCountForSample) <- geneNames
		write.table(rawCountForSample,file=paste(path,"htseqGenerawCountForSample.txt",sep="/"),sep = "\t",row.names = TRUE,col.names = NA,quote=FALSE)
	}
	if(suffix == "genes.fpkm_tracking"){
		if(sysinfo[1] == "Linux"){
			FPKMCount = read.delim(file=paste(path,"cufflinksFPKM.txt",sep="/"),header=TRUE,stringsAsFactors=FALSE)
		}
		if(sysinfo[1] == "Windows"){
			FPKMCount = read.delim(file=paste(path,"cufflinksFPKM.txt",sep="//"),header=TRUE,stringsAsFactors=FALSE)
		}
		geneNames = FPKMCount$GENE_ID
		FPKMCount$GENE_ID <- NULL
		FPKMCountTrans <- as.data.frame(t(FPKMCount))
		colnames(FPKMCountTrans) <- geneNames
		runID <- rownames(FPKMCountTrans)
		#debug 2: save(FPKMCountTrans,file="FPKMCountTrans.RData")
		matchPosition <- match(runID,relatRunAndSample[,1])
		RemainRunID <- runID[which(!is.na(matchPosition))]
		FPKMCountTrans <- FPKMCountTrans[RemainRunID,]
		matchPosition <- matchPosition[!is.na(matchPosition)]
		FPKMCountTrans[which(is.na(FPKMCountTrans),arr.ind=TRUE)] <- 0
		samples <- unique(relatRunAndSample[matchPosition,2])
		FPKMCountForSample <- samples
		for(geneName in geneNames){
			cat("Gene ID:",geneName,"\n")
			geneExpressionValue <- tapply(FPKMCountTrans[,match(geneName,colnames(FPKMCountTrans))],relatRunAndSample[matchPosition,2],sum)
			cat("The corresponding sample name:",names(geneExpressionValue),"\n")
			cat("The corresponding expression value:",as.numeric(geneExpressionValue),"\n")
			FPKMCountForSample <- cbind(FPKMCountForSample,as.numeric(geneExpressionValue[samples]))
		}
		FPKMCountForSample <- as.data.frame(FPKMCountForSample)
		FPKMCountForSample[,1] <- NULL
		FPKMCountForSample = t(FPKMCountForSample)
		colnames(FPKMCountForSample) <- samples
		rownames(FPKMCountForSample) <- geneNames
		FPKMCountForSample[which(is.na(FPKMCountForSample),arr.ind=TRUE)] <- 0
		write.table(FPKMCountForSample,file=paste(path,"cufflinksGeneFPKMCountforSample.txt",sep="/"),sep = "\t",row.names = TRUE,col.names = NA,quote=FALSE)
	}
}

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

### extracting the relationships between runs and samples from the NCBI run information file
extractRelationBetweenRunAndSample <- function(path,filename){
	runInfo <- read.csv(file=paste(path,filename,sep="/"),header=TRUE,stringsAsFactors=FALSE)
	results <- cbind(runInfo$time,runInfo$id)
	return(results)
}

# args <- commandArgs(trailingOnly=T)
# options(warn = -1)
# cat(args,"\n");
# if (length(args) < 3){
# 	print("incorrect input parameters\n")
# 	quit(state=1)
# }

# loading necessary R libraries
library(DESeq2)
library(edgeR)

### invoking functions to conduct below tasks
extractGeneCount("C:/Users/jimta/Desktop","htseq_count.txt")
relatRunAndSample <- extractRelationBetweenRunAndSample("C:/Users/jimta/Desktop","phenodata.csv")
combineGeneCount("C:/Users/jimta/Desktop","htseq_count.txt",relatRunAndSample)
normalizationData("C:/Users/jimta/Desktop","htseq_count.txt")