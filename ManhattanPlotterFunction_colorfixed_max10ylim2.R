#############################Begin Function#############################

#First we define the function with required parameters and optional parameters with their defaults

ManhattanPlot_AJS_cut <- function(chr, pos, pvalue, snp, genomesig = 0, genomesug = 0,photoname = "Manhattan Plot", outname = "ManhattanPlot_AJS", colors = c("navyblue","gray69"), sigcolor = "darkred", sugcolor = "indianred", ncex = .5, plotsymbol = 18, sugcex = 1, sigcex = 1.5, nonautosome = c(23,24,25,26),xlabel = "Chromosomal Positions",ylabel = "-log10(p-value)", pdf = "TRUE",topsnplabels = "FALSE", pvalue_miss = "NA") {
	
	

	#bind the input data into one table and prune out missing data based on the pvalue_miss option and then 
	#redefines the input variables chr, pos, pvalue, snp to reflect the cleaned data
	data_cl <- data.frame(snp,chr,pos,pvalue)
	data_cl <- subset(data_cl,data_cl$pvalue != pvalue_miss)
	chr <- data_cl$chr
	pos <- data_cl$pos
	pvalue <- data_cl$pvalue
	snp <- data_cl$snp
	
	#convert genomesig and genomsug pvalue options into the –log10 form for use on the plot
	if (genomesug != 0) {
		genomesug <- -log10(genomesug)
	}
	if (genomesig != 0) {
		genomesig <- -log10(genomesig)
	}
	
	#create a list of unique chromosome codes
	chroms <- as.numeric(levels(as.factor(chr)))
	
	#create a table with important chromosome data for each unique chromosome code.
	chrominfo <- data.frame()
	colflag <- 0
	
	#loop through each unique chromosome code and store in a table the code [1], the max bp position [2], 
	#the cumulative bp position to mark the beginning of the chromosome on the plot [3], the cumulative 
	#midpoint of each chromosome to draw the label at [4], the label to be used on the plot for each 
	#chromosome [5], the color code for the chromosome [6]
	for (i in 1:length(chroms)) {
		chrominfo[i,1] <- chroms[i]
		chrominfo[i,2] <- max(pos[chr == chroms[i]], na.rm = T)
		ifelse(i == 1, chrominfo[i,3] <- 0, chrominfo[i,3] <- chrominfo[i-1,2] + chrominfo[i-1,3])
		chrominfo[i,4] <- chrominfo[i,3] + (.5*chrominfo[i,2])
		if (chroms[i] <= 22) {
			chrominfo[i,5] <- chroms[i]
		} else if (chroms[i] == nonautosome[1]) {
			chrominfo[i,5] <- "X"
		} else if (chroms[i] == nonautosome[2]) {
			chrominfo[i,5] <- "Y"
		} else if (chroms[i] == nonautosome[3]) {
			chrominfo[i,5] <- "XY"
		} else if (chroms[i] == nonautosome[4]) {
			chrominfo[i,5] <- "MT"
		}
		if (colflag == 0) {
			chrominfo[i,6] <- colors[1]
			colflag <- 1
		} else {
			chrominfo[i,6] <- colors[2]
			colflag <- 0
		}
		#5/2/13 AM edit, this makes it so that the unplaced stuff is ALWAYS grey, so the color choice is uniform (meta analysis files do not always have unplaced SNPs)
		if (chroms[i] == 0) 
		{
			chrominfo[i,6] <- colors[2]
			colflag <- 0
		}
	}
	
	#define a list of the sizes of each point on the graph based on its significance level and the 
	#defined cex for ncex, sugcex, and sigcex and the defined genomesig and genomesug options
	ptcex <- pvalue
	for (i in 1:length(pvalue)){ 
		if (genomesig == 0) {
			if (genomesug == 0) {
				ptcex <- ncex
			}else {
				if (-log10(pvalue[i]) > genomesug) {
					ptcex[i] <- sugcex
				}else (ptcex[i] <- ncex)
			}
		}else {
			if (-log10(pvalue[i]) > genomesig) {
				ptcex[i] <- sigcex
			}else {
				if (genomesug == 0) {
					ptcex[i] <- ncex
				}else {
					if (-log10(pvalue[i]) > genomesug) {
						ptcex[i] <- sugcex
					} else (ptcex[i] <- ncex)
				}
			}
		}
	}
	
	#create an output file, either pdf or jpeg based on specified pdf option
	if (pdf == "TRUE") {
		pdf(file = paste(outname, ".pdf", sep = ""),11,8.5)
	} else {
		jpeg(filename = paste(outname, ".jpeg", sep = ""),width = 11, height = 8.5, units = "in", res = 500,quality = 100)
	}
	
	#define the maximum value along the x chromosome to be used for defining plot and axis sizes 
	xmax <- chrominfo[length(chrominfo[,2]),2] + chrominfo[length(chrominfo[,3]),3]
	
	#defining the plot window size within the output file by specifying margin sizes and mute automatic 
	#labels to allow custom labels
	#bty set to l to remove boundaries on upper right
	par(mar=c(5,5,5,5),xaxs = "i", yaxs = "i", bty="l", las=1)
	
	#plot the points and define size of axes
	plot(pos + chrominfo[,3][match(chr,chrominfo[,1])],-log10(pvalue), col = chrominfo[,6][match(chr,chrominfo[,1])],cex = ptcex, pch = plotsymbol, ann = F, xaxt = "n", xlim = c(0,1.015*xmax), ylim = c(2, max(10)),cex.axis=1.45, cex.lab=1.5)
	
	#label the plot and axes
	title(main = photoname, col.main = colors[1], cex.main = 2, font.main = 2, xlab = xlabel,ylab = ylabel,cex.lab=1.5, cex.axis=1.45 )
	
	#draw and label axis tick marks on x axis
	axis(1, at = chrominfo[,4],labels = chrominfo[,5],cex.axis=1.3, cex.lab=1.5)
	
	
	#draw significant and suggestive lines according to genomesig and genomesug
	genomesigcheck <- 1
	genomesugcheck <- 1
	if (genomesig != 0) {
		if (TRUE | (max(-log10(pvalue)) > genomesig)) {
			abline(h=genomesig,lty = 1, col = sigcolor)
			text(1.0075*xmax,genomesig,label = "",col = sigcolor, xpd = T, pos = 4)
			genomesigcheck <- 1
		}
	}
	if (genomesug != 0) {
		if (TRUE | (max(-log10(pvalue)) > genomesug)) {
			abline(h=genomesug,lty = 2, col = sugcolor)
			text(1.0075*xmax,genomesug,label = "",col = sugcolor, xpd = T, pos = 4)
			genomesugcheck <- 1
		}
	}
	
	#label snps according to topsnplabels
	if (topsnplabels == "ALL") {
		for (i in chrominfo[,1]) {
			textxcord <- pos[pvalue == min(pvalue[chr == i])] + chrominfo[,3][chrominfo[,1] == i]
			textycord <- -log10(min(pvalue[chr == i]))
			textlabel <- snp[pvalue == min(pvalue[chr == i])]
			textpvalue <- min(pvalue[chr == i])
			if (genomesigcheck == 1) {
				if (textycord > genomesig) {
					text(textxcord,textycord,textlabel,pos = 3,cex = sigcex*(2/3), col = sigcolor, offset = sigcex, xpd = T)
					text(textxcord,textycord,paste("p=",textpvalue),pos = 3,cex = sigcex*(2/3), col = sigcolor, xpd = T)
				}
			}
			if (genomesugcheck == 1) {
				if (genomesig == 0) {
					if (textycord > genomesug) {
						text(textxcord,textycord,textlabel,pos = 3,cex = sugcex*(2/3), col = sugcolor, offset = sugcex, xpd = T)
						text(textxcord,textycord,paste("p=",textpvalue),pos = 3,cex = sugcex*(2/3), col = sugcolor, xpd = T)
					}
				}
				else {
					if (textycord > genomesug) {
						if (textycord < genomesig) {
							text(textxcord,textycord,textlabel,pos = 3,cex = sugcex*(2/3), col = sugcolor, offset = sugcex, xpd = T)
							text(textxcord,textycord,paste("p=",textpvalue),pos = 3,cex = sugcex*(2/3), col = sugcolor, xpd = T)

						}
					}
				}
			}
		}
	}
	else {
		if (topsnplabels == "TOP") {
			textxcord <- pos[pvalue == min(pvalue)] + chrominfo[,3][chrominfo[,1] == chr[pvalue == min(pvalue)]]
			textycord <- -log10(min(pvalue))
			textlabel <- snp[pvalue == min(pvalue)]
			textpvalue <- min(pvalue)
			if (genomesigcheck == 1) {
				if (textycord > genomesig) {
					text(textxcord,textycord,textlabel,pos = 3,cex = sigcex*(2/3), col = sigcolor, offset = sigcex, xpd = T)
					text(textxcord,textycord,paste("p=",textpvalue),pos = 3,cex = sigcex*(2/3), col = sigcolor, xpd = T)
				}
			}
			if (genomesugcheck == 1) {
				if (genomesig == 0) {
					if (textycord > genomesug) {
						text(textxcord,textycord,textlabel,pos = 3,cex = sugcex*(2/3), col = sugcolor, offset = sugcex, xpd = T)
						text(textxcord,textycord,paste("p=",textpvalue),pos = 3,cex = sugcex*(2/3), col = sugcolor, xpd = T)
					}
				}
				else {
					if (textycord > genomesug) {
						if (textycord < genomesig) {
							text(textxcord,textycord,textlabel,pos = 3,cex = sugcex*(2/3), col = sugcolor, offset = sugcex, xpd = T)
							text(textxcord,textycord,paste("p=",textpvalue),pos = 3,cex = sugcex*(2/3), col = sugcolor, xpd = T)
						}
					}
				}
			}
		}
		else {
			if (topsnplabels == "SIG") {
				for (i in chrominfo[,1]) {
					textxcord <- pos[pvalue == min(pvalue[chr == i])] + chrominfo[,3][chrominfo[,1] == i]
					textycord <- -log10(min(pvalue[chr == i]))
					textlabel <- snp[pvalue == min(pvalue[chr == i])]
					textpvalue <- min(pvalue[chr == i])
					if (genomesigcheck == 1) {
						if (textycord > genomesig) {
							text(textxcord,textycord,textlabel,pos = 3,cex = sigcex*(2/3), col = sigcolor, offset = sigcex, xpd = T)
							text(textxcord,textycord,paste("p=",textpvalue),pos = 3,cex = sigcex*(2/3), col = sigcolor, xpd = T)
						}
					}
				}
			}
		}
	}
	#close output file
	dev.off()
}

#############################End Function#############################
