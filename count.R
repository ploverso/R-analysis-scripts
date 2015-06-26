
library(Rsubread)
vitro <- c("rnor5", "hg38", "mm10")
origDir <- getwd()
for(species in vitro){
	
	cells <- c("ast", "neu", "opc")
	if(species == "hg38"){
		cells <- c("ast", "neu")
	}
	
	setwd(origDir)
	setwd(species)
	annotLoc <- paste0("../annotations/", species, ".gtf")
	
	for(cell in cells){
		cellData <- paste0(cell, ".bam")
		
		isPair <- F
		if(species == "mm10"){
			isPair <- T
		}
		
		featureCounts(files,annot.ext=annotLoc,isGTFAnnotationFile=T,
		GTF.featureType="exon",GTF.attrType="gene_id",useMetaFeatures=TRUE,
		allowMultiOverlap=FALSE,isPairedEnd=isPair,nthreads=4L,reportReads=T)
	}
	
}
