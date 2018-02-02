args <- commandArgs(trailingOnly = TRUE)
 snpresults <- args[1]
 gender <- args[2]
 descriptor <-args[3]
 anc <- args[4]
 methodz <- args[5]
 outfile <- args[6]

# snpresults="results_cat/rs73154700.use"
# descriptor="pgc_ptsd_study_order_v7.csv"
# anc="eur"
# methodz="FE"
# outfile="rs73154700.plot_RE"

#Load metafor library and plyr (for mapped values)
library(metafor)
library(plyr)

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

if (gender == "males")
{
 descriptors$study <- descriptors$study_males
 descriptors <- subset(descriptors, included_male == 1)
}
if (gender == "females")
{
 descriptors$study <- descriptors$study_females
 descriptors <- subset(descriptors, included_fem == 1)
}

print(dim(descriptors)[1])

#Load data (should be in described format of file name followed by 
dat0 <- read.table(snpresults,header=F,stringsAsFactors=F,colClasses = c("character", "character", "character", "character", 
"numeric", "numeric",  "numeric",  "numeric",  "numeric",  "numeric",  "numeric",  "numeric"))
names(dat0) <- c("study","snp","A1","A2","Maf","Info","Ncases","Ncontrols","Neff","OR","SE","P")

#Put ORs on log scale (betas)
dat0$B <- log(dat0$OR)
dat0$B2 <- NA


##Align alleles

#Set a reference for pivoting (Just choose row 1 in the data)
pivot=c(dat0[1,"A1"],dat0[1,"A2"])




#For every row, check allele alignment
for ( i in 1:dim(dat0)[1])
{
 betaval=dat0[i,]$B
 snp=c(dat0[i,"A1"],dat0[i,"A2"])
 flipped=mapvalues(snp, c("A", "C", "G", "T"),  c("T", "G", "C", "A"))

 if (pivot[1] == snp[1] & pivot[2] == snp[2])
 {
  #No changes need to be made, alleles already aligned
 } else if (pivot[1] == snp[2] & pivot[2] == snp[1]) 
 {
  #A1 and A2 need to be flipped, change sign of beta value
   betaval=-1*betaval
   print(paste("Direction flip for", dat0[i,],sep=" "))
 } else if ( !(flipped[1] %in% pivot) |  !(flipped[2] %in% pivot)) 
 {
  #Check if allele might be flippable, if not, give NA value 
  betaval=NA
 } else if (flipped[1] == snp[1] &  flipped[2] == snp[2]) 
 {
   #It's on the wrong strand but nothing actually needs to be done
 } else if (flipped[1] == snp[2] &  flipped[2] == snp[1]) 
 {
  #It's on the wrong strand and needs to be flipped
   print(paste("Direction and allele flip for", dat0[i,]$study,sep=" "))
   betaval=-1*betaval
 } else {
 betaval=NA
 }
 
 #Now reassign adjusted beta
 dat0[i,]$B2 <- betaval

}
                        


sink(file=paste('plots/',outfile,"_",anc,'.log',sep=''))

#Merge in study descriptors for the plot
dat <- merge(dat0,descriptors,by="study")

dat <- dat[order(dat$cases,decreasing=c(T)),]
#Which studys are NOT in the original data 
print("Not included (missing)")
print(descriptors[-which(descriptors$study %in% dat$study),]$study)

print("N Studies analyzed:")
print(dim(dat)[1])


#These will match if the number of inputs in the starting directory 
#matches the number of inputs in the study description file

sum(dat$B * 1/dat$SE^2) / sum(1/dat$SE^2)


meta_res <- rma(yi=dat$B2,sei=dat$SE,method=methodz,slab=dat$abbr)
print("P value is")
meta_res$pval
print("B value is")
meta_res$estimate

save(meta_res,file=paste('results_cat/',outfile,"_",anc,'.R',sep=''))
pdf(paste('plots/',outfile,"_",anc,'_',gender,'.pdf',sep=''),7,9)
forest(meta_res)
dev.off()

#Write output file for metasoft
dout <- unlist(matrix(c(dat$B, dat$SE),ncol=2,byrow=F))

write.table(t(c(outfile,t(dout + 0.0000001))),paste("results_cat/", outfile,"_",anc,'.msin',sep=''),quote=F,row.names=F,col.names=F)
write.table(paste(dat$abbr,sep="" ),paste("results_cat/", outfile,"_",anc,'.msinstudynames',sep=''),quote=F,row.names=F,col.names=F) # Label by Study Number and abbr
write.table(dat$caroline_order2,paste("results_cat/", outfile,"_",anc,'.msinstudyorder',sep=''),quote=F,row.names=F,col.names=F) #Order by number of cases

print(warnings())

