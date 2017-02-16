
#Read .csv file containing an ordered list of studies to be included in the meta analysis
dat <- read.csv('study_order.csv',header=T,na.strings=c("NA","#N/A"),stringsAsFactors=F)

#Example data head
#Group	ORder	Study	EA	AA	HNA	EA2	AA2
#1	1	MRS	daner_mrsc_eur_analysis_run2.gz.metal 	daner_mrsc_aam_analysis_run2.gz.metal 	daner_mrsc_lat_analysis_run2.gz.metal	#N/A	#N/A
#1	2	ONG (More military)	daner_onga_eur_analysis_run2.gz.metal	#N/A	#N/A	#N/A	#N/A
#1	3	SAFR (Adding African diversity)	#N/A	daner_safr_oth_analysis_run2.gz.metal	#N/A	#N/A	daner_safr_afr_analysis_run3.gz.metal
#1	6	DNHS (more AA)	#N/A	daner_dnhs_aam_analysis_run3.gz.metal	#N/A	#N/A	#N/A


#note:

#i is the incrementor for overall progress
#eai is the incrementor for EA meta
#aai is hte incrementor for AA meta
#lati is the incrementor for LAT meta

#Note I use AA2 and EA2 because sometimes I added two datasets at once, e.g. if one study has 2 components.

ioutdat <- c()
eaoutdat <- c()
aaoutdat <- c()
hnaoutdat <- c()

for (i in 1:dim(dat)[1])
{
 ioutdat  <- c(ioutdat, na.omit(unlist(dat[i,4:8])))
 iout <- paste("PROCESS", t(ioutdat))
 iout <- c(iout, paste('OUTFILE results/all_',i,'_', ' .tbl',sep=''))

 ioutfilename <- paste('temporary_files/all',i,'.mi',sep='')
 write.table(iout, ioutfilename ,quote=F,row.names=F,col.names=F)
 system(paste('cat metal_inputs/meta_header.mi ', ioutfilename, ' metal_inputs/meta_footer.mi >', ioutfilename,'.mif',sep=""))
 
 if(!is.na(dat[i,]$EA) | !is.na(dat[i,]$EA2))
 {
   eaoutdat  <- c(eaoutdat, na.omit(unlist(dat[i,c("EA","EA2")])))
   eaout <- paste("PROCESS", t(eaoutdat))

   eaout <- c(eaout, paste('OUTFILE results/eur_',i,'_', ' .tbl',sep=''))

   eaoutfilename <- paste('temporary_files/eur',i,'.mi',sep='')
   write.table(eaout, eaoutfilename ,quote=F,row.names=F,col.names=F)
   system(paste('cat metal_inputs/meta_header.mi ', eaoutfilename, ' metal_inputs/meta_footer.mi >', eaoutfilename,'.mif',sep=""))
 }

 if(!is.na(dat[i,]$AA) | !is.na(dat[i,]$AA2))
 {
   aaoutdat  <- c(aaoutdat, na.omit(unlist(dat[i,c("AA","AA2")])))
   aaout <- paste("PROCESS", t(aaoutdat))

   aaout <- c(aaout, paste('OUTFILE results/aam_',i,'_', ' .tbl',sep=''))

   aaoutfilename <- paste('temporary_files/aam',i,'.mi',sep='')
   write.table(aaout, aaoutfilename ,quote=F,row.names=F,col.names=F)
   system(paste('cat metal_inputs/meta_header.mi ', aaoutfilename, ' metal_inputs/meta_footer.mi >', aaoutfilename,'.mif',sep=""))
 }

 if(!is.na(dat[i,]$HNA))
 {
   hnaoutdat  <- c(hnaoutdat, na.omit(unlist(dat[i,c("HNA")])))
   hnaout <- paste("PROCESS", t(hnaoutdat))

   hnaout <- c(hnaout, paste('OUTFILE results/lat_',i,'_', ' .tbl',sep=''))

   hnaoutfilename <- paste('temporary_files/lat',i,'.mi',sep='')
   write.table(hnaout, hnaoutfilename ,quote=F,row.names=F,col.names=F)
   system(paste('cat metal_inputs/meta_header.mi ', hnaoutfilename, ' metal_inputs/meta_footer.mi >', hnaoutfilename,'.mif',sep=""))
 }

}

 
