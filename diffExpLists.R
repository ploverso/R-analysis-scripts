#Usage: place in the same directory as your gene_exp.diff file. Run using:


#Read the cuffdiff output
cuff.read=read.delim(file="gene_exp.diff", sep="\t")

#---------------------------------
#Everything between these lines is used to detect whether there is more than one comparison being made
#And if so, to separate each comparison from one another
y <- levels(cuff.read$sample_1)
for(x in levels(cuff.read$sample_2)){
	if(x %in% y){}else{y <- c(y, x)
	}
}
f <- function (snames, pairlist){
	if(length(snames) == 1){
		return(pairlist)
	}
	lHead <- head(snames, n=1)
	lTail <- tail(snames, n=-1)
	for(i in seq(length(lTail))){
		#print(i)
		newPair <- c(lHead, lTail[i])
		#print(newPair)
		pairlist <- c(pairlist, list(newPair))
	}
	#print("recursing")
	pairlist <- f(lTail, pairlist)
	return(pairlist)
}
pairs <- f(y, list())

comps <- list()

for(pair in pairs){
	p1 <- pair[1]
	p2 <- pair[2]
	theComp <- cuff.read[(cuff.read$sample_1==p1 & cuff.read$sample_2==p2),]
	if(length(theComp[,1]) == 0){
		theComp <- cuff.read[(cuff.read$sample_1==p2 & cuff.read$sample_2==p1),]
	}
	comps <- c(comps, list(theComp))
}
#---------------------------------




for(retVal in comps){
	p1 <- retVal$sample_1[1]
	p2 <- retVal$sample_2[1]
	up1 <- retVal$significant == "yes" & retVal$log2.fold_change > 0
	up1 <- retVal$gene[up1]
	up2 <- retVal$significant == "yes" & retVal$log2.fold_change < 0
	up2 <- retVal$gene[up2]
	pn1 <- paste0(p1, "-", p2, "-up-", substr(p1,1,3), ".diff")
	pn2 <- paste0(p1, "-", p2, "-up-", substr(p2,1,3), ".diff")
	write.table(up1, file = pn1, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
	write.table(up2, file = pn2, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
}
