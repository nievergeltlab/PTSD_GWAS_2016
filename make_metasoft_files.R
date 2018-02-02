#Read .csv file containing an ordered list of studies to be included in the meta analysis
 args <- commandArgs(trailingOnly = TRUE)
 dataset <- args[1]
 gender <- args[2]
 
 if (gender == "mf")
 {
  gender=""
 }
dat <- read.csv(dataset,header=T,na.strings=c("NA","#N/A"),stringsAsFactors=F)

paste(dat[,4])
for (chrom in c(1:22))
{
 ioutdat <- paste(na.omit(unlist(dat[,4:14])),chrom,sep="_")
 eaoutdat <- paste(na.omit(unlist(dat[,c(4,9,10,11)])),chrom,sep="_")
 aaoutdat <- paste(na.omit(unlist(dat[,c(5,12,13,14)])),chrom,sep="_")
 hnaoutdat <-  paste(na.omit(unlist(dat[,c(6)])),chrom,sep="_")

 write.table(t(ioutdat),  paste('metasoft_scripts/',gender,'/metasoft_all_chr',chrom,'.msoftin',sep=""),row.names=F,col.names=F,quote=F)
 write.table(t(eaoutdat), paste('metasoft_scripts/',gender,'/metasoft_eur_chr',chrom,'.msoftin',sep=""),row.names=F,col.names=F,quote=F)
 write.table(t(aaoutdat), paste('metasoft_scripts/',gender,'/metasoft_aam_chr',chrom,'.msoftin',sep=""),row.names=F,col.names=F,quote=F)
 write.table(t(hnaoutdat),paste('metasoft_scripts/',gender,'/metasoft_lat_chr',chrom,'.msoftin',sep=""),row.names=F,col.names=F,quote=F)

}