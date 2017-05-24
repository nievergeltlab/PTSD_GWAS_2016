#Read .csv file containing an ordered list of studies to be included in the meta analysis
dat <- read.csv('study_order.csv',header=T,na.strings=c("NA","#N/A"),stringsAsFactors=F)
mi_header <- read.csv('meta_header.mi',header=F,stringsAsFactors=F)$V1
mi_footer <- read.csv('meta_footer.mi',header=F,stringsAsFactors=F)$V1

#Example data head
#Make sure that COMPLETE file names are there

#Group	ORder	Study	EA	AA	HNA	EA2	AA2
#1	1	MRS	daner_mrsc_eur_analysis_run2.gz 	daner_mrsc_aam_analysis_run2.gz	daner_mrsc_lat_analysis_run2.gz	#N/A	#N/A
#1	2	ONG (More military)	daner_onga_eur_analysis_run2.gz	#N/A	#N/A	#N/A	#N/A
#1	3	SAFR (Adding African diversity)	#N/A	daner_safr_oth_analysis_run2.gz	#N/A	#N/A	daner_safr_afr_analysis_run3.gz
#1	6	DNHS (more AA)	#N/A	daner_dnhs_aam_analysis_run3.gz	#N/A	#N/A	#N/A


#note:

#i is the incrementor for overall progress
#eai is the incrementor for EA meta
#aai is hte incrementor for AA meta
#lati is the incrementor for LAT meta

#Note I use AA2 and EA2 because sometimes I added two datasets at once, e.g. if one study has 2 components.

#Feb 16 tbd: add in the chromosome number using hte paste ocmmand


for (chrom in c(1:22))
{
 ioutdat <- c()
 eaoutdat <- c()
 aaoutdat <- c()
 hnaoutdat <- c()



 eaoutdat  <- paste(na.omit(unlist(dat[,c("EA","EA2","EA3","EA4")])),chrom,sep="_")
 eaout <- paste("PROCESS ",  "metal_inputs/",t(eaoutdat),sep="")
 write.table(eaoutdat,'ea_studies_list.txt')

 
 for (test in 1:length(eaoutdat))  
 {

  eaoutp <- eaout[-test] #, paste('OUTFILE results/LOOeur_',test,'_',chrom,  '_ .tbl',sep=''))
  eaoutp2 <- c(eaoutp, paste('OUTFILE results/LOOeur_',test,'_',chrom,  '_ .tbl',sep=''))

  eaoutfilename <- paste('metal_scripts/LOOeur_',test,'_',chrom, '.mi',sep='')
  #write.table(c(mi_header,eaoutp2,mi_footer), eaoutfilename ,quote=F,row.names=F,col.names=F)
 }




}
