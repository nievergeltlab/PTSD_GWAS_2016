
#Read .csv file containing an ordered list of studies to be included in the meta analysis
dat <- read.csv('study_order.csv',header=T,na.strings=c("NA","#N/A"),stringsAsFactors=F)
mi_header <- read.csv('meta_header.mi',header=F)
mi_footer <- read.csv('meta_footer.mi',header=F)

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

#Feb 16 tbd: add in the chromosome number using hte paste ocmmand

ioutdat <- c()
eaoutdat <- c()
aaoutdat <- c()
hnaoutdat <- c()

for (chrom in c(1:22))
 {
 for (i in 1:dim(dat)[1]) #If just doing ALL data and not interested in intermediates, take the last line of the data only (the dimension)
 {
  ioutdat  <- c(ioutdat, paste(na.omit(unlist(dat[i,4:8])),chrom,sep="_")) #Here 4:8 are the relevant columns. May have to change if adding new populations or working with reduced column set
  iout <- paste("PROCESS", t(ioutdat))
  iout <- c(iout, paste('OUTFILE results/all_',i,'_',chrom, ' .tbl',sep='')) #Results file name
  ioutfilename <- paste('temporary_files/all',i,chrom,'.mi',sep='') #Metal script input file name
  write.table(c(mi_header,iout,mi_footer), ioutfilename ,quote=F,row.names=F,col.names=F)

  if(!is.na(dat[i,]$EA) | !is.na(dat[i,]$EA2))
  {
    eaoutdat  <- c(eaoutdat, paste(na.omit(unlist(dat[i,c("EA","EA2")])),chrom,sep="_"))
    eaout <- paste("PROCESS", t(eaoutdat))

    eaout <- c(eaout, paste('OUTFILE results/eur_',i,'_',chrom,  ' .tbl',sep=''))

    eaoutfilename <- paste('temporary_files/eur',i,chrom, '.mi',sep='')
    write.table(c(mi_header,eaout,mi_footer), eaoutfilename ,quote=F,row.names=F,col.names=F)
  }

  if(!is.na(dat[i,]$AA) | !is.na(dat[i,]$AA2))
  {
    aaoutdat  <- c(aaoutdat, paste(na.omit(unlist(dat[i,c("AA","AA2")])),chrom,sep="_"))
    aaout <- paste("PROCESS", t(aaoutdat))

    aaout <- c(aaout, paste('OUTFILE results/aam_',i,'_',chrom,  ' .tbl',sep=''))

    aaoutfilename <- paste('temporary_files/aam',i,chrom, '.mi',sep='')
    write.table(c(mi_header,aaout,mi_footer), aaoutfilename ,quote=F,row.names=F,col.names=F)
  }

  if(!is.na(dat[i,]$HNA))
  {
    hnaoutdat  <- c(hnaoutdat, paste(na.omit(unlist(dat[i,c("HNA")])),chrom,sep="_"))
    hnaout <- paste("PROCESS", t(hnaoutdat))

    hnaout <- c(hnaout, paste('OUTFILE results/lat_',i,'_',chrom,  ' .tbl',sep=''))

    hnaoutfilename <- paste('temporary_files/lat',i,chrom, '.mi',sep='')
    write.table(c(mi_header,hnaout,mi_footer), hnaoutfilename ,quote=F,row.names=F,col.names=F)
  }

 }
}
 
