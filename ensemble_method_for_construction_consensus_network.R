integrateNetworks <- function(dataPath,geneNames){
	library(corpcor)
	wgcnaNets <- list.files(dataPath,pattern="Normalized_Count_wgcna_network_final_edgeIndex\\.RData$",full.name=TRUE)
	bc3netNets <- list.files(dataPath,pattern="Normalized_Count_bc3net_network_final_edgeIndex\\.RData$",full.name=TRUE)
	ggmNets <- list.files(dataPath,pattern="Normalized_Count_ggm_network_final_edgeIndex\\.RData$",full.name=TRUE)
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
	belowIndexes <- which(wgcnaNumMatrix < 4,arr.ind=TRUE)
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
	belowIndexes <- which(bc3netNumMatrix < 3,arr.ind=TRUE)
	bc3netAdjacencyMatrix[belowIndexes] <- 0
	bc3netBoolMatrix[belowIndexes] <- 0
	bc3netAdjacencyMatrix <- bc3netAdjacencyMatrix/(length(bc3netNets)-1)
	
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
	belowIndexes <- which(ggmNumMatrix < 4,arr.ind=TRUE)
	ggmAdjacencyMatrix[belowIndexes] <- 0
	ggmBoolMatrix[belowIndexes] <- 0
	ggmAdjacencyMatrix <- ggmAdjacencyMatrix/(length(ggmNets)-1)
	
	adjacencyMatrix <- (wgcnaAdjacencyMatrix + bc3netAdjacencyMatrix + ggmAdjacencyMatrix)/3
	boolMatrix <- wgcnaBoolMatrix + bc3netBoolMatrix + ggmBoolMatrix
	candiIndexes <- which(boolMatrix >= 2,arr.ind=TRUE)
	cat("candiIndexes:",candiIndexes[1,],"\n")
	edgeList <- cbind(candiIndexes,adjacencyMatrix[candiIndexes])
	cat("edgeList:",edgeList[1:5,3],"\n")
	sortIndex <- order(edgeList[,3],decreasing=TRUE)
	edgeList <- edgeList[sortIndex,]
	edgeList[,3] <- (edgeList[,3]-min(edgeList[,3]))/(max(edgeList[,3])-min(edgeList[,3]))
	save(edgeList,file=paste(dataPath,"final_geneNet_edgeIndex.RData",sep="/"))
	netMatrix <- cbind(geneNames[edgeList[,1]],geneNames[edgeList[,2]],edgeList[,3])
	write.table(netMatrix,file=paste(dataPath,"final_geneNet.tab",sep="/"),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
}


getInterSet <- function(filelists){
	expData <- read.table(file=filelists[1],header=TRUE,row.names=1,sep=" ",stringsAsFactors=FALSE)
	commonGenes <- rownames(expData)
	commonSamples <- colnames(expData)
	for(filename in filelists[2:length(filelists)]){
		expData <- read.table(file=filename,header=TRUE,row.names=1,sep=" ",stringsAsFactors=FALSE)	
		commonGenes <- intersect(commonGenes,rownames(expData))
		commonSamples <- intersect(commonSamples,colnames(expData))
	}
	cat("The common genes between htseq-count and cufflinks are:", commonGenes, "\n")
	cat("The common samples between htseq-count and cufflinks are:", commonSamples, "\n")
	return(commonGenes)
}

args <- commandArgs(TRUE)
filelists <- list.files(args[1],pattern="_Count\\.txt$",full.name=TRUE)
commonGenes <- getInterSet(filelists)
integrateNetworks(args[1],commonGenes)

# Rscript --vanilla /public/Project/Project_YuHua/CoexpressNetwork/Program/networkAnalysis/integrateNetworks.R /public/Project/Project_YuHua/CoexpressNetwork/Data/RiceData/networkfile