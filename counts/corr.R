
myCPM <- cpm(procCounts)
myCPM <- data.frame(log2(myCPM + 1))

corResults <- data.frame(matrix(0, nrow=ncol(myCPM), ncol=ncol(myCPM)))
rownames(corResults) <- colnames(myCPM)
colnames(corResults) <- colnames(myCPM)

for(id1 in 1:length(colnames(myCPM))){
	for(id2 in 1:length(colnames(myCPM))){
		
		corResults[id1, id2] <- cor(myCPM[,id1], myCPM[,id2], method="spearman")
		
	}
}
