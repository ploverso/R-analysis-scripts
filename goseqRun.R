#Usage: place in the same directory as your gene_exp.diff file. Run using:
#> source("goseqRun.R")
#NOTE: this script only works as-is for R. norvegicus data
#If you have another species, change the lines marked "#Species-specific"

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
		newPair <- c(lHead, lTail[i])
		pairlist <- c(pairlist, list(newPair))
	}
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
require(goseq)
require(KEGGREST)

#Species-specific:
require(org.Rn.eg.db)

for(retVal in comps){
	#Identify the pairing, using the first letter of each sample name.
	#If you have similarly names samples feel free to change the naming sceme
	cuff.read <- retVal
	p1 <- cuff.read$sample_1[1]
	p2 <- cuff.read$sample_2[1]
	pairID <- paste0(substr(p1,1,1), substr(p2,1,1))
	print(paste0("pairID: ", pairID))

	#Generate a named vector where the contents are a 0/1 boolean representing whether the gene is significantly expressed
	#and where the names of the elements are the ENSEMBL gene names
	isdif <- c()
	dif1 <- c()
	dif2 <- c()
	i <- 0
	for (test in cuff.read$significant){
		i <- i + 1
		if (test == "yes"){
			isdif <- c(isdif, 1)
			if(cuff.read$log2fold_change[i] > 0){
				dif1 <- c(dif1, 1)
				dif2 <- c(dif2, 0)
			}else{
				dif1 <- c(dif1, 0)
				dif2 <- c(dif2, 1)
			}
		}else{
			isdif <- c(isdif, 0)
			dif1 <- c(dif1, 0)
			dif2 <- c(dif2, 0)
		}
	}
	names(isdif) <- cuff.read$gene_id
	names(dif1) <- cuff.read$gene_id
	names(dif2) <- cuff.read$gene_id

	#Generate a pwf for the gene set, and apply the mapping to it.
	#Species-specific
	pwf <- nullp(isdif, "rn5", "ensGene", plot.fit=F)
	#Species-specific
	gos <- goseq(pwf, "rn5", "ensGene", test.cats=c("KEGG"))
	
	pwf1 <- nullp(dif1, "rn5", "ensGene", plot.fit=F)
	#Species-specific
	gos1 <- goseq(pwf1, "rn5", "ensGene", test.cats=c("KEGG"))
	
	pwf2 <- nullp(dif2, "rn5", "ensGene", plot.fit=F)
	#Species-specific
	gos2 <- goseq(pwf2, "rn5", "ensGene", test.cats=c("KEGG"))
	
	des <- c()
	des1 <- c()
	des2 <- c()
	for(kcat in gos1$category){
		print(kcat)
		kcat <- paste0("rno", kcat)
		x <- keggGet(kcat)
		#Get KEGG description
		kdes <- x[[1]][2]
		des1 <- c(des1, kdes)
	}
	for(kcat in gos2$category){
		print(kcat)
		kcat <- paste0("rno", kcat)
		x <- keggGet(kcat)
		#Get KEGG description
		kdes <- x[[1]][2]
		des2 <- c(des2, kdes)
	}
	for(kcat in gos$category){
		print(kcat)
		kcat <- paste0("rno", kcat)
		x <- keggGet(kcat)
		#Get KEGG description
		kdes <- x[[1]][2]
		des <- c(des, kdes)
	}
	desc <- as.vector(unlist(des1))
	gos1 <- cbind(gos1, desc)
	desc <- as.vector(unlist(des2))
	gos2 <- cbind(gos2, desc)
	desc <- as.vector(unlist(des))
	gos <- cbind(gos, desc)

	#Output the results
	pway <- paste0(pairID, "goseq")
	dir.create(pway)
	setwd(pway)
	pway <- paste0(pairID, "-goseq.csv")
	pway1 <- paste0(p1, "-goseq.csv")
	pway2 <- paste0(p2, "-goseq.csv")
	write.table(gos, file=pway, sep="\t")
	write.table(gos1, file=pway1, sep="\t")
	write.table(gos2, file=pway2, sep="\t")
	setwd("..")
}
