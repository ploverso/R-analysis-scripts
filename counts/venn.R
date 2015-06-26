
#~ source("headTail.R")
#~ setwd("diffExp")
#~ library("VennDiagram")
#~ 
#~ fls <- dir()
#~ geneLists <- list()
#~ for(fl in fls){
#~ 	
#~ 	diffData <- read.delim(fl)
#~ 	fl <- strhead(fl, -4)
#~ 	diffData <- diffData[diffData$FDR < 0.05,]
#~ 	
#~ 	geneLists[[fl]] <- diffData[,1]
#~ 	geneLists[[paste0(fl, ".up")]] <- diffData[,1][diffData$logFC > 0]
#~ 	geneLists[[paste0(fl, ".down")]] <- diffData[,1][diffData$logFC < 0]
#~ 	
#~ }
#~ 
#~ setwd("..")

#~ venn.diagram(geneLists[c(44,89,98)], filename="./vennDiagrams/venn_vtVSvv_up.png", imagetype="png", 
#~ 	width=1500, height=1500, resolution=400, cat.dist=-0.05, category.names=c("Astrocyte", "Neuron", "OPC"))
#~ venn.diagram(geneLists[c(45,90,99)], filename="./vennDiagrams/venn_vtVSvv_down.png", imagetype="png", 
#~ 	width=1500, height=1500, resolution=400, cat.dist=-0.05, category.names=c("Astrocyte", "Neuron", "OPC"))

venn.diagram(geneLists[c(44,95,116)], filename="./vennDiagrams/venn_rnVSvvCorr_up.png", imagetype="png", 
	width=1500, height=1500, resolution=400, cat.dist=-0.05, category.names=c("Astrocyte", "Neuron", "OPC"))
venn.diagram(geneLists[c(45,96,117)], filename="./vennDiagrams/venn_rnVSvvCorr_down.png", imagetype="png", 
	width=1500, height=1500, resolution=400, cat.dist=-0.05, category.names=c("Astrocyte", "Neuron", "OPC"))


#~ venn.diagram(geneLists[c(23,29,74)], filename="./vennDiagrams/venn_vtVSvvCellcomp_up.png", imagetype="png",
#~ 	width=1500, height=1500, resolution=400, cat.dist=c(0.03, 0.06, 0.06),
#~ 	category.names=c("Ast/Neu", "Ast/OPC", "Neu/OPC"), cat.pos=c(180, 180, 0))
#~ venn.diagram(geneLists[c(27,33,78)], filename="./vennDiagrams/venn_vtVSvvCellcomp_down.png", imagetype="png", 
#~ 	width=1500, height=1500, resolution=400, cat.dist=c(0.03, 0.06, 0.06),
#~ 	category.names=c("Ast/Neu", "Ast/OPC", "Neu/OPC"), cat.pos=c(180, 180, 0))

#~ venn.diagram(geneLists[c(20,26,71)], filename="./vennDiagrams/venn_rnVSvvCellcomp_up.png", imagetype="png",
#~ 	width=1500, height=1500, resolution=400, cat.dist=c(0.03, 0.06, 0.06),
#~ 	category.names=c("Ast/Neu", "Ast/OPC", "Neu/OPC"), cat.pos=c(180, 180, 0))
#~ venn.diagram(geneLists[c(24,30,75)], filename="./vennDiagrams/venn_rnVSvvCellcomp_down.png", imagetype="png", 
#~ 	width=1500, height=1500, resolution=400, cat.dist=c(0.03, 0.06, 0.06),
#~ 	category.names=c("Ast/Neu", "Ast/OPC", "Neu/OPC"), cat.pos=c(180, 180, 0))

#~ j <- 1
#~ for(i in names(geneLists)){
#~ 	cat(i, "\t")
#~ 	cat(length(geneLists[[i]]), "\t")
#~ 	if(j == 3){
#~ 		j <- 1
#~ 		cat("\n")
#~ 	}else{j <- j + 1}
#~ 	
#~ }
