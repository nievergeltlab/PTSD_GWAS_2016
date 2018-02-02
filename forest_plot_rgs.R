args <- commandArgs(trailingOnly = TRUE)
 descriptor <-args[1] #Genetic correlation data, as a .csv
 sigval <- args[2] #Cutoff FDR Q value for which traits to include in plo
 traitlabel <- args[3] #Label for trait that will appear on plot
 groups <- args[4] #Column defining subgroups
 ukbb <- args[5]
 outfile <- args[6] #Output name

 #Sample usage:
 #Rscript forest_plot_rgs.R rgv2.csv 0.05 trait2 Category no rg_overall
  #Rscript forest_plot_rgs.R rgfemv2.csv 0.05 trait2 Category no rg_fem
 
 
 source('rma_hidden.r')
 source("forest.rma2.r")
sigval <- as.numeric(sigval)
 print(is.numeric(sigval))
 
# snpresults="results_cat/rs73154700.use"
# descriptor="pgc_ptsd_study_order_v7.csv"
# anc="eur"
# methodz="FE"
# outfile="rs73154700.plot_RE"

#Load metafor library and plyr (for mapped values)
library(metafor)

#Load rg file
descriptors <- read.csv(descriptor,header=T,stringsAsFactors=F,na.strings=c("NA","#N/A","#VALUE!"))


#Only plot FDR significant results
dat <- subset(descriptors, descriptors$redundant == 0 & descriptors$p_bonf <= sigval)

if(ukbb == "no")
{
 dat <- subset(dat, dat$Category != "ukbb")
}

#Order the data by category
dat <- dat[order(dat$Category),]

#Which studys are NOT in the original data 

print("N significant traits being plotted:")
print(dim(dat)[1])
print(head(dat))

#For each category, assign a color
descriptors$colors  <- NA
catcolors <- scan(what=character())
psychiatric blue
personality lightblue
reproductive red
sleeping green
anthropometric lightgreen
aging orange
education darkorange
smoking_behaviour brown
cardiometabolic grey
autoimmune navy

catcolors <- data.frame(matrix(catcolors,ncol=2,byrow=T),stringsAsFactors=F)
names(catcolors) <- c("category","color")

print(unique(dat$Category))
head(catcolors)
dat$colors <-NA 

for (categories in unique(dat$Category))
{
 #Match color

 dat[which(dat$Category == categories),]$colors <- catcolors[which(catcolors$category == categories),]$color

 }
 print(data.frame(dat$trait2, dat$colors))
#These will match if the number of inputs in the starting directory 
#matches the number of inputs in the study description file

meta_res <- rma(yi=dat$rg,sei=dat$se,slab=paste(dat[,traitlabel]," (p=",formatC(dat$p, format = "e", digits = 0),")",sep=""))

print(names(meta_res))
print(meta_res$slab)
#save(meta_res,file=paste('results_cat/',outfile,'.R',sep=''))
pdf(paste('plots/',outfile,"_.pdf",sep=''),7,9)
forest.rma2(meta_res, col=rev(dat$colors),pch=16,psize=1,addfit=F,annotate=F,xlab="Genetic Correlation (rg)",at=c(-1,-0.75,-0.5,-0.25,0,.25,.5,.7,1),ilab.pos=4) #Notice I revoerse ordered hte colors, because there is somehow a problem with it
dev.off()

