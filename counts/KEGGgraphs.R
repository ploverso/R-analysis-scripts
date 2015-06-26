


#Map the gene lists and fold change values to KEGG using the M. musculus kegg set from gageData
library(org.Mm.eg.db)
library(pathview)
library(gage)
library(gageData)

#~ counts <- read.delim("countsEG.txt")
sel.rn <- rowSums(cnts) != 0
counts <- counts[sel.rn,]

egMap <- as.list(org.Mm.egENSEMBL2EG)
egCounts  <- data.frame(matrix(nrow=0, ncol=ncol(counts)))
for(gNameInd in 1:nrow(counts)){
	egName <- egMap[[strsplit(rownames(counts)[gNameInd], ".", fixed=T)[[1]][1]]]
	if(!is.null(egName)){
		for(gName in egName){
			if(!gName %in% rownames(egCounts)){
				egCounts <- rbind(egCounts, counts[gNameInd,])
				rownames(egCounts)[nrow(egCounts)] <- gName
			}
		}
	}
}

#Species-specific
data(kegg.sets.mm)

analyzePair <- function(pairID, cnts.norm, ref.idx, samp.idx, compType="paired"){
	
	out.suffix="edger"
	fc.kegg.p <- gage(cnts.norm, gsets = kegg.sets.mm, ref = NULL, samp = NULL)
	sel <- fc.kegg.p$greater[, "q.val"] < 0.1 & !is.na(fc.kegg.p$greater[, "q.val"])
	path.ids <- rownames(fc.kegg.p$greater)[sel]
	sel.l <- fc.kegg.p$less[, "q.val"] < 0.1 & !is.na(fc.kegg.p$less[,"q.val"])
	path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
	path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)

	#Use the KEGG R. norvegicus pathway charts to chart our own expression patterns
	pway <- paste0(pairID, "pathways")
	dir.create(pway)
	setwd(pway)
	#Species-specific
	pv.out.list <- sapply(path.ids2, function(pid) pathview(gene.data = cnts.norm, pathway.id = pid, species = "mmu", out.suffix=out.suffix))
	write.csv(path.ids2, file="ids.csv", quote=F)
	write.table(fc.kegg.p, file="geneData.txt", quote=F, sep="\t")
	setwd("..")
}

analyzePair("ast.vtVSvv", cnts.norm, c(1,2), c(7, 8))
