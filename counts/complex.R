
ensToSym <- function(ensNames){
	library(biomaRt)
	
	genes2 <- c()
	for(gene in ensNames){
		geneShr <- strsplit(gene, ".", fixed=T)[[1]][1]
		genes2 <- c(genes2, geneShr)
	}
	
	mart<- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
	symNames <- as.character(getBMlist(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "external_gene_name"),
		values=genes2,mart= mart, giveWarning=F)$external_gene_name)
	return(symNames)
}


library(edgeR)
counts <- read.delim("counts.txt")
Group <- factor(substring(colnames(counts), 1, 5))
keep <- rowSums(cpm(counts)) >= 3
counts <- counts[keep,]
procCounts <- DGEList(counts=counts, group=Group)
procCounts <- calcNormFactors(procCounts)
procCounts <- estimateCommonDisp(procCounts)
procCounts <- estimateTagwiseDisp(procCounts)
design <- model.matrix(~0+Group)
colnames(design) <- levels(Group)
procCounts <- estimateGLMCommonDisp(procCounts,design)
procCounts <- estimateGLMTrendedDisp(procCounts,design)
procCounts <- estimateGLMTagwiseDisp(procCounts,design)
fit <- glmFit(procCounts,design)


my.contrasts <- makeContrasts(
	rn.astVSrest = rnast - (rnneu + rnopc)/2,
	rn.neuVSrest = rnneu - (rnast + rnopc)/2,
	rn.opcVSrest = rnopc - (rnneu + rnast)/2,
	mm.astVSrest = mmast - (mmneu + mmopc)/2,
	mm.neuVSrest = mmneu - (mmast + mmopc)/2,
	mm.opcVSrest = mmopc - (mmneu + mmast)/2,
	ast.vtVSvv = (hgast + rnast)/2 - mmast,
	neu.vtVSvv = (hgneu + rnneu)/2 - mmneu,
	opc.vtVSvv = rnopc - mmopc,
	all.vtVSvv = (hgast + hgneu + rnast + rnneu + rnopc)/5 - (mmast + mmneu + mmopc)/3,
	all.rnVSmm = (rnast + rnneu + rnopc)/3 - (mmast + mmneu + mmopc)/3,
	ast.rnVShg = rnast - hgast,
	ast.rnVSmm = rnast - mmast,
	ast.hgVSmm = hgast - mmast,
	neu.rnVShg = rnneu - hgneu,
	neu.rnVSmm = rnneu - mmneu,
	neu.hgVSmm = hgneu - mmneu,
	rn.astVSneu = rnast - rnneu,
	mm.astVSneu = mmast - mmneu,
	hg.astVSneu = hgast - hgneu,
	rn.astVSopc = rnast - rnopc,
	mm.astVSopc = mmast - mmopc,
	rn.neuVSopc = rnneu - rnopc,
	mm.neuVSopc = mmneu - mmopc,
	vt.astVSneu = (rnast + hgast)/2 - (rnneu + hgneu)/2,
	vv.astVSneu = mmast - mmneu,
	vt.astVSopc = (rnast + hgast)/2 - rnopc,
	vv.astVSopc = mmast - mmopc,
	vt.neuVSopc = (rnneu + hgneu)/2 - rnopc,
	vv.neuVSopc = mmneu - mmopc,
	astneu.vtVSvv = ((rnast - rnneu) + (hgast - hgneu))/2 - (mmast - mmneu),
	astopc.vtVSvv = ((rnast + hgast)/2 - rnopc) - (mmast - mmopc),
	neuopc.vtVSvv = ((rnneu + hgneu)/2 - rnopc) - (mmneu - mmopc),
	astRest.vtVSvv = ((rnast + hgast)/2 - (rnneu + hgneu + rnopc)/3) - (mmast - (mmneu + mmopc)/2),
	neuRest.vtVSvv = ((rnneu + hgneu)/2 - (rnast + hgast + rnopc)/3) - (mmneu - (mmast + mmopc)/2),
	opcRest.vtVSvv = (rnopc - (rnneu + hgneu + rnast + hgast)/4) - (mmopc - (mmneu + mmast)/2),
	astneu.rnVSvv = (rnast - rnneu) - (mmast - mmneu),
	astopc.rnVSvv = (rnast - rnopc) - (mmast - mmopc),
	neuopc.rnVSvv = (rnneu - rnopc) - (mmneu - mmopc),
	astRest.rnVSvv = (rnast - (rnneu + rnopc)/2) - (mmast - (mmneu + mmopc)/2),
	neuRest.rnVSvv = (rnneu  - (rnast + rnopc)/2) - (mmneu - (mmast + mmopc)/2),
	opcRest.rnVSvv = (rnopc - (rnneu + rnast)/2) - (mmopc - (mmneu + mmast)/2),
	vt.astVSrest = (rnast + hgast)/2 - (rnneu + hgneu + rnopc)/3,
	vt.neuVSrest = (rnneu + hgneu)/2 - (rnast + hgast + rnopc)/3,
	vt.opcVSrest = rnopc - (rnneu + hgneu + rnast + hgast)/4,
	all.astVSneu = (rnast + mmast + hgast)/3 - (rnneu + mmneu + hgneu)/3,
	all.astVSrest = (rnast + mmast + hgast)/3 - (rnneu + mmneu + hgneu + rnopc + mmopc)/5,
	all.neuVSrest = (rnneu + mmneu + hgneu)/3 - (rnast + mmast + hgast + rnopc + mmopc)/5,
	all.opcVSrest = (rnopc + mmopc)/2 - (rnneu + mmneu + hgneu + rnast + mmast + hgast)/6,
	
	ast.rnmm2 = (rnast + mmneu + mmopc)/3 - (mmast + rnneu + rnopc)/3,
	neu.rnmm2 = (mmast + rnneu + mmopc)/3 - (rnast + mmneu + rnopc)/3,
	opc.rnmm2 = (mmast + mmneu + rnopc)/3 - (rnast + rnneu + mmopc)/3,
	
	ast.rnmm3 = (rnast - (rnast + rnneu + rnopc)/3) - (mmast - (mmast + mmneu + mmopc)/3),
	neu.rnmm3 = (rnneu - (rnast + rnneu + rnopc)/3) - (mmneu - (mmast + mmneu + mmopc)/3),
	opc.rnmm3 = (rnopc - (rnast + rnneu + rnopc)/3) - (mmopc - (mmast + mmneu + mmopc)/3),
	
levels=design)

#~ keggDiff <- list()
#~ geneDiff <- list()
#~ fullLRT <- list()

new <- c("ast.rnmm3", "neu.rnmm3", "opc.rnmm3")
source("diffPath.R")
#~ for(id in colnames(my.contrasts)){
for(id in new){
	
	lrt <- glmLRT(fit, contrast=my.contrasts[,id])
	fcTable <- data.frame(topTags(lrt, n=nrow(lrt$table)))
	geneDiff[[id]] <- data.frame(topTags(lrt, n=nrow(lrt$table)))
	fullLRT[[id]] <- lrt
	logFC <- fcTable[,1]
	names(logFC) <- rownames(fcTable)
	if(length(logFC) > 0){
		
		keggDiff[[id]] <- keggMap(id, logFC)
	}else{
		keggDiff[[id]] <- ""
	}
}


logCPM <- cpm(procCounts)
logCPM <- data.frame(log2(logCPM + 1))
source("createHeatmaps.R")

#~ for(compName in names(geneDiff)){
for(compName in new){
	print(compName)
	thisLRT <- fullLRT[[compName]]
	topGenesUp <- geneDiff[[compName]][geneDiff[[compName]]$logFC > 0,]
	topGenesDown <- geneDiff[[compName]][geneDiff[[compName]]$logFC < 0,]
	topGenes <- rbind(topGenesUp[1:20,], topGenesDown[1:20,])
	sampCols <- c()
	comps <- thisLRT$comparison
	comps <- strsplit(gsub("[^[:alpha:] ]", "", comps), " ")[[1]]
	for(comp in comps){
		sampCols <- c(sampCols, which(thisLRT$design[,comp] == 1))
	}
	sampCols <- sort(sampCols)
	names(sampCols) <- NULL
	sampCols <- logCPM[,sampCols]
	sampCols <- sampCols[rownames(topGenes),]
	rownames(sampCols) <- ensToSym(rownames(sampCols))
	rownames(sampCols)[which(rownames(sampCols) == "NA")] <- rownames(topGenes)[which(rownames(sampCols) == "NA")]
#~ 	dendroHeatmap(sampCols, paste0("40 most significant DEGs: ", compName), paste0("./outputs/", compName, ".png"), doSort=F)
}


#~ for(compName in names(geneDiff)){
#~ 	write.table(cbind(rownames(geneDiff[[compName]]), geneDiff[[compName]]),
#~ 		file=paste0("./diffExp/", compName, ".txt"), quote=F, row.names=F, sep="\t")
#~ }

#~ spiaDIR <-"../../SPIA/R/"
#~ for(spiaFile in dir(spiaDIR)){
#~ 	source(paste0(spiaDIR, spiaFile))
#~ }
#~ source("spiaRun.R")

#~ spiaPways <- list()


pwayMake <- function(compName){
	spiaPways[[compName]] <- spiaRun(compName)
	expData <- geneDiff[[compName]]
	expData <- expData[expData$FDR < 0.05,]
	expData <- data.frame(expData$logFC, row.names=rownames(expData))
	setwd("./pathways/spia/")
	keggGraph(paste0(compName, "-spia"), spiaPways[[compName]], expData)
	setwd("../gage")
	keggGraph(paste0(compName, "-gage"), keggDiff[[compName]], expData)
	setwd("../..")
}
#~ mclapply(names(geneDiff), pwayMake, mc.cores = getOption("mc.cores", 5L))
#~ mclapply(new, pwayMake, mc.cores = getOption("mc.cores", 5L))


#~ keggGraph("./pathways/gage/vtVSvv.byCell.combined", c(keggDiff[["ast.vtVSvv"]], keggDiff[["neu.vtVSvv"]], keggDiff[["opc.vtVSvv"]]),
#~ 	data.frame(geneDiff[["ast.vtVSvv"]]$logFC, geneDiff[["neu.vtVSvv"]]$logFC, geneDiff[["opc.vtVSvv"]]$logFC,
#~ 	row.names=rownames(fit$counts)))
#~ 
#~ keggGraph("./pathways/gage/astVSneu.bySpecies.combined", c(keggDiff[["rn.astVSneu"]], keggDiff[["hg.astVSneu"]], keggDiff[["mm.astVSneu"]]),
#~ 	data.frame(geneDiff[["rn.astVSneu"]]$logFC, geneDiff[["hg.astVSneu"]]$logFC, geneDiff[["mm.astVSneu"]]$logFC,
#~ 	row.names=rownames(fit$counts)))

#~ 
#~ logCPM <- cpm(procCounts)
#~ logCPM <- data.frame(log2(logCPM + 1))
#~ source("createHeatmaps.R")
#~ dendroHeatmap(logCPM, "200 Top Expressed Genes By Sample (log2(cpm))", "top200.bySample.png")
#~ 
#~ logCPM <- cpm(procCounts)
#~ logCPM <- data.frame(logCPM[,1], rowMeans(logCPM[,2:3]), rowMeans(logCPM[,4:5]), rowMeans(logCPM[,6:7]), 
#~ 	rowMeans(logCPM[,8:9]), rowMeans(logCPM[,10:11]), rowMeans(logCPM[,12:13]), logCPM[,14])
#~ colnames(logCPM) <- c("Rn_ast", "Rn_neu", "Rn_opc", "Mm_ast", "Mm_neu", "Mm_opc", "Hg_ast", "Hg_neu")
#~ logCPM <- log2(logCPM + 1)
#~ dendroHeatmap(logCPM, "200 Top Expressed Genes By Celltype (log2(cpm))", "top200.byCelltype.png")

#~ cellCPM <- data.frame(matrix(nrow=0, ncol=8))
#~ for(cellType in colnames(logCPM)){
#~ 	logCPM <- logCPM[rev(order(logCPM[,cellType])),]
#~ 	cellCPM <- rbind(cellCPM, head(logCPM, n=25))
#~ }
#~ dendroHeatmap(cellCPM, "25 Top Expressed Genes Per Celltype (log2(cpm))", "top25.perCelltype.png", doSort = F)
