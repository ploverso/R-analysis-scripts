vivoUp <- scan("vivoUp.txt", what="character")
vitroUp <- scan("vitroUp.txt", what="character")

library(biomaRt)
mart<- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
vivoUpM <- as.character(getBMlist(filters= "external_gene_name", attributes= c("ensembl_gene_id", "external_gene_name"),
		values=vivoUp,mart= mart)$ensembl_gene_id)

vitroUpM <- as.character(getBMlist(filters= "external_gene_name", attributes= c("ensembl_gene_id", "external_gene_name"),
		values=vitroUp,mart= mart)$ensembl_gene_id)

source("headTail.R")

vivoUpM <- vivoUpM[vivoUpM != "NA"]
vitroUpM <- vitroUpM[vitroUpM != "NA"]

for(myComp in c("ast.vtVSvv", "ast.rnVSmm", "ast.rnmm3")){
	
	gd <- geneDiff[[myComp]]
	gdVv <- gd[strhead(rownames(gd), 18) %in% vivoUpM,]
	gdVt <- gd[strhead(rownames(gd), 18) %in% vitroUpM,]
	
	write.table(gdVv, file=paste0(myComp, ".vivo.txt"), quote=F, sep="\t")
	write.table(gdVt, file=paste0(myComp, ".vitro.txt"), quote=F, sep="\t")
	
}
