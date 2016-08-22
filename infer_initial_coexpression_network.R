consCoexpressionNetwork_wgcna <- function(path,fileName,species,nThreads,commonGenesAndSamples=NULL){
	library(WGCNA)
	library(corpcor)
	options(stringsAsFactors = FALSE)
	cat(nThreads, "\n")
	normalizedData <- read.table(file=fileName, sep=" ", header=TRUE, row.names=1, stringsAsFactors=FALSE)
	normalizedData <- normalizedData[match(commonGenesAndSamples$Genes,rownames(normalizedData)),match(commonGenesAndSamples$Samples,colnames(normalizedData))]
	normalizedData[1,1] <- as.numeric(normalizedData[1,1])
	normalizedDataMatrix <- data.matrix(t(normalizedData))
	randomMatrix <- matrix(runif((ncol(normalizedDataMatrix)*nrow(normalizedDataMatrix)),min=1e-20,max=2e-20),nrow=nrow(normalizedDataMatrix))
	normalizedDataMatrix <- normalizedDataMatrix + randomMatrix
	geneNames <- colnames(normalizedDataMatrix)
	
	enableWGCNAThreads(nThreads = as.integer(nThreads))
	powers = seq(from=2,to=10,by=2)
	sft <- pickSoftThreshold(normalizedDataMatrix, powerVector = powers, verbose = 5)
	R2value <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
	maxIndex <- which.max(R2value)
	optimalPower <- sft$fitIndices[maxIndex,1]
	adjacencyMatrix <- adjacency(normalizedDataMatrix,power=optimalPower)
	TOMMatrix <- TOMsimilarity(adjacencyMatrix)
	TOMVector <- sm2vec(TOMMatrix)
	indexes = sm.index(TOMMatrix)
	edgeList = cbind(indexes,TOMVector)
	sortIndex <- order(edgeList[,3],decreasing=TRUE)
	edgeList <- edgeList[sortIndex,]
	
	# permutation
	randomRowIndexes <- sample(nrow(normalizedData),10000,replace=FALSE)
	randomMatrix <- normalizedData[randomRowIndexes,]
	for(i in seq(1,length(randomRowIndexes))){
		randomColIndexes <- sample(ncol(normalizedData),ncol(normalizedData),replace=FALSE)
		randomMatrix[i,] <- randomMatrix[i,randomColIndexes] 
	}
	randomMatrix[1,1] <- as.numeric(randomMatrix[1,1])
	randomMatrix <- data.matrix(t(randomMatrix))
	randomsft <- pickSoftThreshold(randomMatrix, powerVector = powers, verbose = 5)
	R2value <- -sign(randomsft$fitIndices[,3])*randomsft$fitIndices[,2]
	maxIndex <- which.max(R2value)
	optimalPower <- randomsft$fitIndices[maxIndex,1]
	randomAdjacencyMatrix <- adjacency(randomMatrix,power=optimalPower)
	randomTOMMatrix <- TOMsimilarity(randomAdjacencyMatrix)
	randomTOMVector <- sm2vec(randomTOMMatrix)
	sortIndex <- order(randomTOMVector,decreasing=TRUE)
	randomTOMVector <- randomTOMVector[sortIndex]
	threValue <- randomTOMVector[round(0.00001*length(randomTOMVector))]
	cat(fileName,threValue,"\n")
	
 	aboveIndexes <- which(edgeList[,3] > threValue)
	edgeList <- edgeList[aboveIndexes,]
	save(edgeList,file=paste(path,paste(fileName,"_wgcna_network_edgeIndex.RData", sep=""),sep="/"))
}

consCoexpressionNetwork_bc3net <- function(path,fileName,species,commonGenesAndSamples=NULL){
	library(bc3net)
	library(igraph)
	normalizedData <- read.table(file=fileName,sep=" ", header=TRUE, row.names=1, stringsAsFactors=FALSE)
	normalizedData <- normalizedData[match(commonGenesAndSamples$Genes,rownames(normalizedData)),match(commonGenesAndSamples$Samples,colnames(normalizedData))]
	normalizedData[1,1] <- as.numeric(normalizedData[1,1])
	normalizedDataMatrix <- data.matrix(normalizedData)
	randomMatrix <- matrix(runif((ncol(normalizedDataMatrix)*nrow(normalizedDataMatrix)),min=1e-20,max=2e-20),nrow=nrow(normalizedDataMatrix))
	normalizedDataMatrix <- normalizedDataMatrix + randomMatrix
	bc3netNetwork <- bc3net(normalizedDataMatrix,boot=100)
	edgeList <- get.data.frame(bc3netNetwork)
	sortIndex <- order(edgeList[,3],decreasing=TRUE)
	edgeList <- edgeList[sortIndex,]
	save(edgeList,file=paste(path,paste(fileName,"_bc3net_network_final_edgeIndex.RData",sep=""),sep="/"))
}

consCoexpressionNetwork_ggm <- function(path,fileName,species,commonGenesAndSamples=NULL){
	library(GeneNet)
	library(corpcor)
	normalizedData <- read.table(file=fileName, sep=" ", header=TRUE, row.names=1, stringsAsFactors=FALSE)
	normalizedData <- normalizedData[match(commonGenesAndSamples$Genes,rownames(normalizedData)),match(commonGenesAndSamples$Samples,colnames(normalizedData))]
	normalizedData[1,1] <- as.numeric(normalizedData[1,1])
	geneNames <- rownames(normalizedData)
	
	### Estimate partial correlation matrix
	normalizedDataMatrix <- data.matrix(normalizedData)
	normalizedDataMatrix <- t(normalizedDataMatrix)
	randomMatrix <- matrix(runif((ncol(normalizedDataMatrix)*nrow(normalizedDataMatrix)),min=1e-20,max=2e-20),nrow=nrow(normalizedDataMatrix))
	normalizedDataMatrix <- normalizedDataMatrix + randomMatrix
	partialCorMatrix <- ggm.estimate.pcor(normalizedDataMatrix)
	pcor = sm2vec(partialCorMatrix)
	indexes = sm.index(partialCorMatrix)
	edgeList = cbind(indexes,pcor)
	sortIndex <- order(edgeList[,3],decreasing=TRUE)
	edgeList <- edgeList[sortIndex,]
	fdr.out = fdrtool(edgeList[,3],statistic="correlation",plot=FALSE)
	pval = fdr.out$pval
	qval = fdr.out$qval
	prob = 1 - fdr.out$lfdr
	pcor <- edgeList[,3]
	
	### permutation
	randomRowIndexes <- sample(nrow(normalizedData),10000,replace=FALSE)
	randomMatrix <- normalizedData[randomRowIndexes,]
	for(i in seq(1,length(randomRowIndexes))){
		randomColIndexes <- sample(ncol(normalizedData),ncol(normalizedData),replace=FALSE)
		randomMatrix[i,] <- randomMatrix[i,randomColIndexes] 
	}
	randomMatrix[1,1] <- as.numeric(randomMatrix[1,1])
	randomMatrix <- data.matrix(t(randomMatrix))
	randomCorMatrix <- ggm.estimate.pcor(randomMatrix)
	
	randomcor = sm2vec(randomCorMatrix)
	posrandomcor <- randomcor[which(randomcor >= 0)]
	negrandomcor <- abs(randomcor[which(randomcor < 0)])
	sortIndex1 <- order(posrandomcor,decreasing=TRUE)
	posrandomcor <- posrandomcor[sortIndex1]
	sortIndex2 <- order(negrandomcor,decreasing=TRUE)
	negrandomcor <- negrandomcor[sortIndex2]
	
	threValue1 <- posrandomcor[round(0.00001*length(posrandomcor))]
	threValue2 <- negrandomcor[round(0.00001*length(negrandomcor))]
	cat(fileName,threValue1,threValue2,"\n")
	
	aboveIndexes <- which(edgeList[,3] > threValue1 | edgeList[,3] < -threValue2)
	edgeList <- edgeList[aboveIndexes,]
	save(edgeList,file=paste(path,paste(fileName,"_ggm_network_final_edgeIndex.RData",sep=""),sep="/"))
}

getInterSet <- function(filelists){
	expData <- read.table(file=filelists[1],header=TRUE,row.names=1,sep=" ",stringsAsFactors=FALSE)
	commonGenes <- rownames(expData)
	commonSamples <- colnames(expData)
	for(fileName in filelists[2:length(filelists)]){
		expData <- read.table(file=fileName,header=TRUE,row.names=1,sep=" ",stringsAsFactors=FALSE)	
		commonGenes <- intersect(commonGenes,rownames(expData))
		commonSamples <- intersect(commonSamples,colnames(expData))
	}
	cat("The common genes between htseq-count and cufflinks are:", commonGenes, "\n")
	cat("The common samples between htseq-count and cufflinks are:", commonSamples, "\n")
	return(list(Genes=commonGenes,Samples=commonSamples))
}

###############################################################################
args <- commandArgs(trailingOnly=T)
if(length(args) != 4){
	print("incorrect input parameters\n")
	quit(status=1)
}

filelists <- list.files(args[2],pattern="_Count\\.txt$",full.name=TRUE)
cat("The gene expression abundance files are:",filelists,"\n")
commonGenesAndSamples <- getInterSet(filelists)

for (fileName in filelists) { 
	filebasename <- gsub("\\.txt$","",fileName)
	filebasename <- gsub(".*\\/", "", filebasename)
	cat(filebasename,"\n")
	
	### Weighted Gene Co-expression Network Inference
	if(args[1] == "WGCNA"){
		consCoexpressionNetwork_wgcna(args[2],fileName,args[3],commonGenesAndSamples)
	}
	
	### Bagging Based Inference
	if(args[1] == "BC3NET"){
		consCoexpressionNetwork_bc3net(args[2],fileName,commonGenesAndSamples)
	}

	### Graphical Gaussian Model for Network Inference
	if(args[1] == "GGM"){
		consCoexpressionNetwork_ggm(args[2],fileName,commonGenesAndSamples)
	}
}