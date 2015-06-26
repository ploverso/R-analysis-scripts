


spiaRun <- function(compType){
	lrt <- fullLRT[[compType]]
	top <- data.frame(topTags(lrt, n=nrow(lrt$table)))

	library(org.Mm.eg.db)
	egMap <- as.list(org.Mm.egENSEMBL2EG)
	entrez <- c()
	for(gNameInd in 1:nrow(top)){
		egName <- egMap[[strsplit(rownames(top)[gNameInd], ".", fixed=T)[[1]][1]]]
		if(!is.null(egName)){
			for(gName in egName){
				if(!gName %in% entrez){
					entrez <- c(entrez, gName)
					break
				}
			}
			entrez <- c(entrez, -1)
		}else{
			entrez <- c(entrez, -1)
		}
	}

	top <- top[entrez != -1,]
	entrez <- entrez[entrez != -1]

	top <- cbind(top, entrez)

	sig_genes <- subset(top, FDR<0.05)$logFC
	names(sig_genes) <- subset(top, FDR<0.05)$entrez
	all_genes <- top$entrez
	 
	# run SPIA. 
	spia_result <- spia(de=sig_genes, all=all_genes, organism="mmu", plots=F,
		data.dir="/home/peter/Seafile/myFiles/Research/cui/rnaSeq/new/SPIA/inst/extdata/")
	png(paste0("./spiaData/", compType, "-spia.png"))
	diffNames <- plotP(spia_result, threshold=0.1, pert=0.1, nde=0.1)
	dev.off()
	write.table(spia_result[spia_result$ID %in% diffNames,], file=paste0("./spiaData/", compType, "-spia.txt"), 
		sep="\t", quote=F, row.names=F)
	
	return(diffNames)
}
