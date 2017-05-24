args <- commandArgs(trailingOnly = TRUE)
 scriptloc <- args[1]
 results <- args[2]
 outfile <- args[3]
 colorchoice <- args[4]

 blue=rgb(153,191,254,max=255)
 red=rgb(230,185,184,max=255)
 green=rgb(147,205, 221,max=255)
 purple=rgb(179,162,199,max=255)

##Read in data
 print("Loading meta analysis results")
#get number of rows in file (reads file faster)
 nr <- as.numeric(system(paste("wc -l", results, " | awk \'{print $1}\' "),intern=T))
#get list of numeric column classes (reads file faster)
 dat_temp <- read.table(results, header=F,nrows=50,stringsAsFactors=F)
 classes_to_use <- sapply(dat_temp , class)

#read file
 dat1 <- read.table (results, header=T,stringsAsFactors=F,colClasses=classes_to_use ,nrow=nr)
 names(dat1) <- c("SNP","CHR","BP","P")
#Load manhattanplot script
 print("Plotting data")
 source(scriptloc)

 attach(dat1)

    ManhattanPlot_AJS_cut(CHR, BP, P, SNP, genomesig = 5*(10^-8), genomesug = 0,photoname = '', 
    outname = outfile, colors = c(rgb(78,78,77,max=255),eval(parse(text = colorchoice))), sigcolor = 'darkred', sugcolor = 'indianred', ncex = 1, 
    sugcex = 1, sigcex = 1, nonautosome = c(23,24,25,26),xlabel = 'Chromosomal Positions',ylabel = '-log10(p-value)', 
    pdf = 'FALSE',topsnplabels = 'FALSE', pvalue_miss = 'NA', sigsnpcolor = "red")
