
keggMap <- function(pairID, fcTable){
	
	print(paste0("pairID: ", pairID))
	
	#Map the gene lists and fold change values to KEGG using the M. musculus kegg set from gageData
	
	library(pathview)
	library(gage)
	library(gageData)
	library(org.Mm.eg.db)
	egMap <- as.list(org.Mm.egENSEMBL2EG)
	exp.fc <- c()
	for(gNameInd in 1:length(fcTable)){
		egName <- egMap[[strsplit(names(fcTable)[gNameInd], ".", fixed=T)[[1]][1]]]
		if(!is.null(egName)){
			for(gName in egName){
				if(!gName %in% names(exp.fc)){
					exp.fc <- c(exp.fc, fcTable[gNameInd])
					names(exp.fc)[length(exp.fc)] <- gName
				}
			}
		}
	}
	
	#Species-specific
	data(kegg.sets.mm)
	fc.kegg.p <- gage(exp.fc, gsets = kegg.sets.mm, ref = NULL, samp = NULL)
	sel <- fc.kegg.p$greater[, "p.val"] < 0.05 & !is.na(fc.kegg.p$greater[, "p.val"])
	path.ids <- rownames(fc.kegg.p$greater)[sel]
	sel.l <- fc.kegg.p$less[, "p.val"] < 0.05 & !is.na(fc.kegg.p$less[,"p.val"])
	path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
	path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)
	write.table(fc.kegg.p, file=paste0("gageOut/", pairID, ".txt"), quote=F, sep="\t")
	if(length(path.ids2) == 0){
		path.ids2 <- ""
	}
	return(path.ids2)
}

keggGraph <- function(pairID, path.ids, exp.fc){
	library(pathview)
	library(org.Mm.eg.db)
	egMap <- as.list(org.Mm.egENSEMBL2EG)
	egCounts  <- data.frame(matrix(nrow=0, ncol=ncol(exp.fc)))
	for(gNameInd in 1:nrow(exp.fc)){
		egName <- egMap[[strsplit(rownames(exp.fc)[gNameInd], ".", fixed=T)[[1]][1]]]
		if(!is.null(egName)){
			for(gName in egName){
				if(!gName %in% rownames(egCounts)){
					egCounts <- rbind(egCounts, exp.fc[gNameInd,])
					rownames(egCounts)[nrow(egCounts)] <- gName
				}
			}
		}
	}
	
	#Use the KEGG R. norvegicus pathway charts to chart our own expression patterns
	pway <- paste0(pairID, "pathways")
	out.suffix="edger"
	dir.create(pway)
	setwd(pway)
	#Species-specific
	pv.out.list <- sapply(path.ids, function(pid) pathview(gene.data = egCounts, pathway.id = pid,
		species = "mmu", out.suffix=out.suffix, limit=list(gene=5, cpd=5), expand.node=T,
		bins = list(gene = 20, cpd= 20), low = list(gene = "blue", cpd = "blue")))
	
	setwd("..")
}




