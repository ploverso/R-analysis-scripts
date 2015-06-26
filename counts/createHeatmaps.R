
dendroHeatmap <- function(logCPM, chartTitle, fname, doSort = T){
	library(gplots)
	library(RColorBrewer)
	
	my <- logCPM
	if(doSort){
		rowv <- T
		dendro <- "both"
		sortcmd <- "my[with(my, order(-pmax("
		for(colnm in colnames(my)){
			sortcmd <- paste0(sortcmd, colnm, ", ")
		}
		sortcmd <- substr(sortcmd, 1, nchar(sortcmd) - 2)
		sortcmd <- paste0(sortcmd, "))),]")
		my <- eval(parse(text=sortcmd))
	}else{
		rowv <- T
		dendro <- "both"
	}
	#Plot the data using the gplot package
	colfunc <- colorRampPalette(c("red", "red", "red", "white", "white", "blue"))
	colors <- colfunc(200)
	labs <- matrix(row.names(my))
	for(i in 2:ncol(my)){
		labs <- cbind(labs, row.names(my))
	}
	
	my <- as.matrix(my)
	png(fname, width = 600, height = 700)
	heatmap.2(my, col=colors, srtCol = 40, margins = c(5, 5), labRow = labs, keysize = 1.5, lwid=c(1, 8,1),
		lmat=rbind( c(0,0,0), c(0, 3,0), c(2,1,0), c(0, 4,0) ), lhei = c(1,0.5,9,1.5), cexRow=1.2, cexCol=1.55,
		tracecol = "black", dendrogram = dendro, Colv=rowv, Rowv= rowv)
	title(chartTitle, cex.main = 2)
	dev.off()
}
