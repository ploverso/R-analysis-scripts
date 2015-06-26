

load("trueExp_mm.RData")

final <- countedData
counts <-final[[1]][[1]]

for(samp in names(final)[2:length(names(final))]){
	
	newData <- final[[samp]][[1]]
	oldRows <- rownames(counts)
	newRows <- rownames(newData)
	allRows <- intersect(oldRows, newRows)
	
	counts <- data.frame(counts[rownames(counts) %in% allRows,])
	newData <- data.frame(newData[rownames(newData) %in% allRows,])
	
	counts <- cbind(counts, newData)
	
}

colnames(counts) <- names(final)
counts <- counts[order(rownames(counts)),]

write.table(counts, file="counts_mm.txt", quote=F, sep="\t")

