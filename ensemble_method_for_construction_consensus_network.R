######################################################################################################################################################################
#	> File Name: ensemble_method_for_construction_consensus_network.R
#	> Inferring genome-wide gene co-expression network using massive-scale RNA-seq samples
#   > by integrating the predictions of three different network inference algorithms
#	> Author: Hua Yu
#	> Mail: huayu@genetics.ac.cn 
#	Created Time: 2017-10-17
######################################################################################################################################################################

consCoexpressionNetwork_wgcna <- function(path,fileName,nThreads,percThreshold,commonGenesAndSamples=NULL){
	library(WGCNA)
	library(corpcor)
	options(stringsAsFactors = FALSE)
	# cat(nThreads, "\n")
	normalizedData <- read.table(file=fileName, sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)
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
	randomRowIndexes <- sample(nrow(normalizedData),nrow(normalizedData),replace=FALSE)
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
	threValue <- randomTOMVector[round(as.numeric(percThreshold)*length(randomTOMVector))]
	cat(fileName,threValue,"\n")
	
 	aboveIndexes <- which(edgeList[,3] > threValue)
	edgeList <- edgeList[aboveIndexes,]
	fileName <- gsub("\\.txt$","",fileName)
	fileName <- gsub(".*\\/", "", fileName)
	save(edgeList,file=paste(path,paste(fileName,"_wgcna_network_final_edgeIndex.RData", sep=""),sep="/"))
}

consCoexpressionNetwork_bc3net <- function(path,fileName,commonGenesAndSamples=NULL){
	library(bc3net)
	library(igraph)
	normalizedData <- read.table(file=fileName,sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)
	normalizedData <- normalizedData[match(commonGenesAndSamples$Genes,rownames(normalizedData)),match(commonGenesAndSamples$Samples,colnames(normalizedData))]
	normalizedData[1,1] <- as.numeric(normalizedData[1,1])
	normalizedDataMatrix <- data.matrix(normalizedData)
	randomMatrix <- matrix(runif((ncol(normalizedDataMatrix)*nrow(normalizedDataMatrix)),min=1e-20,max=2e-20),nrow=nrow(normalizedDataMatrix))
	normalizedDataMatrix <- normalizedDataMatrix + randomMatrix
	bc3netNetwork <- bc3net(normalizedDataMatrix,boot=100)
	edgeList <- get.data.frame(bc3netNetwork)
	sortIndex <- order(edgeList[,3],decreasing=TRUE)
	edgeList <- edgeList[sortIndex,]
	fileName <- gsub("\\.txt$","",fileName)
	fileName <- gsub(".*\\/", "", fileName)
	save(edgeList,file=paste(path,paste(fileName,"_bc3net_network_final_edgeIndex.RData",sep=""),sep="/"))
}

consCoexpressionNetwork_ggm <- function(path,fileName,percThreshold,commonGenesAndSamples=NULL){
	library(GeneNet)
	library(corpcor)
	normalizedData <- read.table(file=fileName, sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)
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
	randomRowIndexes <- sample(nrow(normalizedData),nrow(normalizedData),replace=FALSE)
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
	
	threValue1 <- posrandomcor[round(as.numeric(percThreshold)*length(posrandomcor))]
	threValue2 <- negrandomcor[round(as.numeric(percThreshold)*length(negrandomcor))]
	cat(fileName,threValue1,threValue2,"\n")
	
	aboveIndexes <- which(edgeList[,3] > threValue1 | edgeList[,3] < -threValue2)
	edgeList <- edgeList[aboveIndexes,]
	fileName <- gsub("\\.txt$","",fileName)
	fileName <- gsub(".*\\/", "", fileName)
	save(edgeList,file=paste(path,paste(fileName,"_ggm_network_final_edgeIndex.RData",sep=""),sep="/"))
}

integrateNetworks <- function(dataPath,geneNames,S1N,S2N){
	library(corpcor)
	wgcnaNets <- list.files(dataPath,pattern="Count_wgcna_network_final_edgeIndex\\.RData$",full.name=TRUE)
	bc3netNets <- list.files(dataPath,pattern="Count_bc3net_network_final_edgeIndex\\.RData$",full.name=TRUE)
	ggmNets <- list.files(dataPath,pattern="Count_ggm_network_final_edgeIndex\\.RData$",full.name=TRUE)
	wgcnaAdjacencyMatrix <- matrix(0,nrow=length(geneNames),ncol=length(geneNames))
	wgcnaBoolMatrix <- matrix(1,nrow=length(geneNames),ncol=length(geneNames))
	wgcnaNumMatrix <- matrix(0,nrow=length(geneNames),ncol=length(geneNames))
	for(wgcnaNetIndex in seq(1,length(wgcnaNets))){
		tempMatrix <- matrix(0,nrow=length(geneNames),ncol=length(geneNames))
		load(wgcnaNets[wgcnaNetIndex])
		edgeList[,3] <- (edgeList[,3]-min(edgeList[,3]))/(max(edgeList[,3])-min(edgeList[,3]))
		tempMatrix[as.matrix(edgeList[,1:2])] <- edgeList[,3]
		wgcnaAdjacencyMatrix <- wgcnaAdjacencyMatrix + tempMatrix
		wgcnaNumMatrix[as.matrix(edgeList[,1:2])] <- wgcnaNumMatrix[as.matrix(edgeList[,1:2])] + 1
	}
	belowIndexes <- which(wgcnaNumMatrix < S1N,arr.ind=TRUE)
	wgcnaAdjacencyMatrix[belowIndexes] <- 0
	wgcnaBoolMatrix[belowIndexes] <- 0
	wgcnaAdjacencyMatrix <- wgcnaAdjacencyMatrix/(length(wgcnaNets)-1)
	
	bc3netAdjacencyMatrix <- matrix(0,nrow=length(geneNames),ncol=length(geneNames))
	bc3netBoolMatrix <- matrix(1,nrow=length(geneNames),ncol=length(geneNames))
	bc3netNumMatrix <- matrix(0,nrow=length(geneNames),ncol=length(geneNames))
	for(bc3netNetIndex in seq(1,length(bc3netNets))){
		tempMatrix <- matrix(0,nrow=length(geneNames),ncol=length(geneNames))
		load(bc3netNets[bc3netNetIndex])
		edgeNodes1 <- match(edgeList[,1],geneNames)
		edgeNodes2 <- match(edgeList[,2],geneNames)
		edgeList <- cbind(edgeNodes1,edgeNodes2,edgeList[,3])
		for(index in seq(1,nrow(edgeList))){
			if(edgeList[index,1] < edgeList[index,2]){
				tmp <- edgeList[index,1]
				edgeList[index,1] <- edgeList[index,2]
				edgeList[index,2] <- tmp
			}
		}
		edgeList[,3] <- (edgeList[,3]-min(edgeList[,3]))/(max(edgeList[,3])-min(edgeList[,3]))
		tempMatrix[as.matrix(edgeList[,1:2])] <- edgeList[,3]
		bc3netAdjacencyMatrix <- bc3netAdjacencyMatrix + tempMatrix
		bc3netNumMatrix[as.matrix(edgeList[,1:2])] <- bc3netNumMatrix[as.matrix(edgeList[,1:2])] + 1 
	}
	belowIndexes <- which(bc3netNumMatrix < S1N,arr.ind=TRUE)
	bc3netAdjacencyMatrix[belowIndexes] <- 0
	bc3netBoolMatrix[belowIndexes] <- 0
	bc3netAdjacencyMatrix <- bc3netAdjacencyMatrix/(length(bc3netNets))
	
	ggmAdjacencyMatrix <- matrix(0,nrow=length(geneNames),ncol=length(geneNames))
	ggmBoolMatrix <- matrix(1,nrow=length(geneNames),ncol=length(geneNames))
	ggmNumMatrix <- matrix(0,nrow=length(geneNames),ncol=length(geneNames))
	for(ggmNetIndex in seq(1,length(ggmNets))){
		tempMatrix <- matrix(0,nrow=length(geneNames),ncol=length(geneNames))
		load(ggmNets[ggmNetIndex])
		edgeList[,3] <- (edgeList[,3]-min(edgeList[,3]))/(max(edgeList[,3])-min(edgeList[,3]))
		tempMatrix[as.matrix(edgeList[,1:2])] <- edgeList[,3]
		ggmAdjacencyMatrix <- ggmAdjacencyMatrix + tempMatrix
		ggmNumMatrix[as.matrix(edgeList[,1:2])] <- ggmNumMatrix[as.matrix(edgeList[,1:2])] + 1
	}
	belowIndexes <- which(ggmNumMatrix < S1N,arr.ind=TRUE)
	ggmAdjacencyMatrix[belowIndexes] <- 0
	ggmBoolMatrix[belowIndexes] <- 0
	ggmAdjacencyMatrix <- ggmAdjacencyMatrix/(length(ggmNets)-1)
	
	adjacencyMatrix <- (wgcnaAdjacencyMatrix + bc3netAdjacencyMatrix + ggmAdjacencyMatrix)/3
	boolMatrix <- wgcnaBoolMatrix + bc3netBoolMatrix + ggmBoolMatrix
	candiIndexes <- which(boolMatrix >= S2N,arr.ind=TRUE)
	cat("candiIndexes:",candiIndexes[1,],"\n")
	edgeList <- cbind(candiIndexes,adjacencyMatrix[candiIndexes])
	cat("edgeList:",edgeList[1:5,3],"\n")
	sortIndex <- order(edgeList[,3],decreasing=TRUE)
	edgeList <- edgeList[sortIndex,]
	edgeList[,3] <- (edgeList[,3]-min(edgeList[,3]))/(max(edgeList[,3])-min(edgeList[,3]))
	# save(edgeList,file=paste(dataPath,"final_geneNet_edgeIndex.RData",sep="/"))
	netMatrix <- cbind(geneNames[edgeList[,1]],geneNames[edgeList[,2]],edgeList[,3])
	write.table(netMatrix,file=paste(dataPath,"final_geneNet.tab",sep="/"),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
}

getInterSet <- function(filelists){
	cat(filelists[1],"\n")
	expData <- read.table(file=filelists[1],header=TRUE,row.names=1,sep="\t",stringsAsFactors=FALSE)
	commonGenes <- rownames(expData)
	commonSamples <- colnames(expData)
	for(fileName in filelists[2:length(filelists)]){
		expData <- read.table(file=fileName,header=TRUE,row.names=1,sep="\t",stringsAsFactors=FALSE)	
		commonGenes <- intersect(commonGenes,rownames(expData))
		commonSamples <- intersect(commonSamples,colnames(expData))
	}
	cat("The common genes between htseq-count and cufflinks are:", commonGenes, "\n")
	cat("The common samples between htseq-count and cufflinks are:", commonSamples, "\n")
	return(list(Genes=commonGenes,Samples=commonSamples))
}


args <- commandArgs(trailingOnly=T)
cat(length(args),"\n")
if(length(args)!=5){
	cat("The number of input parameters are incorrect\n")
}

filelists <- list.files(args[1],pattern="_Count\\.txt$",full.name=TRUE)
commonGenesAndSamples <- getInterSet(filelists)
for (fileName in filelists) { 
	### Weighted gene co-expression network inference
	consCoexpressionNetwork_wgcna(args[1],fileName,args[2],args[3],commonGenesAndSamples)
	
	### Bagging based network inference
	consCoexpressionNetwork_bc3net(args[1],fileName,commonGenesAndSamples)

	### Graphical gaussian model network inference
	consCoexpressionNetwork_ggm(args[1],fileName,args[3],commonGenesAndSamples)
}
integrateNetworks(args[1],commonGenesAndSamples$Genes,args[4],args[5])

# Rscript ensemble_method_for_construction_consensus_network.R /public/Project/Project_YuHua/CoexpressNetwork/Data/RiceData/networkfile/NetMiner-master/testdata 4 0.001 3 2