 args <- commandArgs(trailingOnly = TRUE)
 snpresults <- args[1]
 descriptor <-args[2]
 anc <- args[3]

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
dat0 <- read.table(paste('results_cat/',snpresults,'.use',sep=''),header=F,stringsAsFactors=F)
names(dat0) <- c("study","snp","A1","A2","OR","SE","P")

#Put ORs on log scale (betas)
dat0$B <- log(dat0$OR)

#Merge in study descriptors for the plot
dat <- merge(dat0,descriptors,by="study")

#Which studys are NOT in the original data 
print("Not included (missing)")
print(descriptors[-which(descriptors$study %in% dat$study),]$study)

print("N Studies analyzed:")
print(dim(dat)[1])

#Collapse data for metasoft
ddim <- dim(dat[,c("B","SE")])
testmat <- rep(NA,ddim[1]*ddim[2])

testmat[seq(1,length(testmat),by=2)] <- dat$B
testmat[seq(2,length(testmat),by=2)] <- dat$SE

testmat <- matrix(testmat,nrow=1)
print(testmat)
#row.names(testmat) <- descriptors$abbr
write.table(testmat,paste(snpresults,'.metasoft_in',sep=''),row.names=T,quote=F,col.names=F)

system(paste('java -jar Metasoft.jar -input ', paste(snpresults,'.metasoft_in',sep=''),  '-output', paste(snpresults,'.metasoft_out')))




