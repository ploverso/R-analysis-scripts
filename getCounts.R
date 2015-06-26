
load("trueExp_hg.RData")

final <- countedData
counts <-final[[1]][[1]]
counts <- data.frame(counts[order(rownames(counts)),])

for(samp in names(final)[2:length(names(final))]){
	newData <- final[[samp]][[1]]
	newData <- data.frame(newData[order(rownames(newData)),])
	newData <- data.frame(newData[rownames(newData) %in% rownames(counts),])
	rows <- rownames(counts)[rownames(counts) %in% rownames(newData)]
	counts <- data.frame(counts[rownames(counts) %in% rownames(newData),])
	
	rownames(counts) <- rows
	
	newData <- newData[names(newData)]
	
	counts <- cbind(counts, newData)
}

colnames(counts) <- names(final)

write.table(counts, file="counts_hg.txt", quote=F, sep="\t")
