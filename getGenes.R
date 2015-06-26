#This script will extract all genes into a table making it easy to compare across all samples

strtail <- function(s, n) {
   	if(n < 0){
   		return(substring(s, 1 - n))
   	}
   	else{
    	return(substring(s, nchar(s) - n + 1))
   	}
}

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

retTab <- data.frame(comps[[1]]$gene)
colnames(retTab) <- c("geneNames")
used <- c()
for(retVal in comps){
	s1 <- retVal$sample_1[1]
	s1 <- levels(factor(s1))
	s2 <- retVal$sample_2[1]
	s2 <- levels(factor(s2))
	if(s1 %in% used){}else{
		used <- c(used, s1)
		retTab <- cbind(retTab, retVal$value_1)
		colnames(retTab)[length(colnames(retTab))] <- s1
	}
	if(s2 %in% used){}else{
		used <- c(used, s2)
		retTab <- cbind(retTab, retVal$value_2)
		colnames(retTab)[length(colnames(retTab))] <- s2
	}
}

write.table(retTab, file = "expByGene.diff", sep = "\t", row.names = FALSE, quote = FALSE, col.names = T)


