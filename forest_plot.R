args <- commandArgs(trailingOnly = TRUE)
 snpresults <- args[1]
 descriptor <-args[2]
 anc <- args[3]
 outfile <- args[4]

#Load metafor library
library(metafor)

#load study descriptors
descriptors <- read.csv(descriptor,header=T,stringsAsFactors=F)

#Use only subset of data of given ancestry
if (anc != "all") 
{ 
 if (anc =="aam") 
 { 
  descriptors <- subset(descriptors,ancestry %in% c("aam","saf-afr" )) 
 }else {
   descriptors <- subset(descriptors,ancestry == anc) 
  } 
}



print(dim(descriptors)[1])

#Load data (should be in described format of file name followed by 
dat0 <- read.table(snpresults,header=F,stringsAsFactors=F)
names(dat0) <- c("study","snp","A1","A2","OR","SE","P")

#Put ORs on log scale (betas)
dat0$B <- log(dat0$OR)


sink(file=paste('plots/',outfile,"_",anc,'.log',sep=''))

#Merge in study descriptors for the plot
dat <- merge(dat0,descriptors,by="study")

dat <- dat[order(dat$cases,decreasing=T),]
#Which studys are NOT in the original data 
print("Not included (missing)")
print(descriptors[-which(descriptors$study %in% dat$study),]$study)

print("N Studies analyzed:")
print(dim(dat)[1])


#These will match if the number of inputs in the starting directory 
#matches the number of inputs in the study description file



meta_res <- rma(yi=dat$B,sei=dat$SE,method="FE",slab=dat$abbr)
print("P value is")
meta_res$pval
pdf(paste('plots/',outfile,"_",anc,'.pdf',sep=''),7,9)
forest(meta_res)
dev.off()
