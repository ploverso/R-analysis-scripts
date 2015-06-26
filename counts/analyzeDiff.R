
source("headTail.R")
setwd("diffExp")

comps <- c("rn.astVSrest.txt", "rn.neuVSrest.txt", "rn.opcVSrest.txt",
	"mm.astVSrest.txt", "mm.neuVSrest.txt", "mm.opcVSrest.txt", "hg.astVSneu.txt", "")

topGenes <- data.frame(matrix(ncol=8, nrow=40))
colnames(topGenes) <- comps
i <- 1
for(comp in comps){
	
	
	if(i == 8){
		mat <- read.delim("hg.astVSneu.txt")
		mat$logFC <- mat$logFC * -1
		mat <- mat[order(mat$logFC, decreasing = T),]
		mat <- mat[mat$logCPM > 1.3,]
		
		genes1 <- mat[,1][1:40]
		genes2 <- c()
		for(gene in genes1){
			geneShr <- strsplit(gene, ".", fixed=T)[[1]][1]
			genes2 <- c(genes2, geneShr)
		}
		topGenes[,i] <- genes2
	}else{
		mat <- read.delim(comp)
		mat <- mat[order(mat$logFC, decreasing = T),]
		mat <- mat[mat$logCPM > 1.3,]
		
		genes1 <- mat[,1][1:40]
		genes2 <- c()
		for(gene in genes1){
			geneShr <- strsplit(gene, ".", fixed=T)[[1]][1]
			genes2 <- c(genes2, geneShr)
		}
		topGenes[,comp] <- genes2
		i <- i + 1
	}
}
setwd("..")

colnames(topGenes) <- c("rn.Astrocyte", "rn.Neuron", "rn.OPC", "mm.Astrocyte", "mm.Neuron", "mm.OPC", "hg.Astrocyte", "hg.Neuron")

library(biomaRt)
mart<- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
for(cell in colnames(topGenes)){
	topGenes[,cell] <- as.character(getBMlist(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "external_gene_name"),
		values=topGenes[,cell],mart= mart)$external_gene_name)
	
}
