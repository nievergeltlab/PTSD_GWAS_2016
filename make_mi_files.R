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

 for (i in 1:dim(dat)[1]) #If just doing ALL data and not interested in intermediates, take the last line of the data only (the dimension), otherwise 1:dim(dat)[1]
 {
  ioutdat  <- c(ioutdat, paste(na.omit(unlist(dat[i,4:14])),chrom,sep="_")) #Here 4:14 are the relevant columns. May have to change if adding new populations or working with reduced column set
  iout <- paste("PROCESS ", "metal_inputs/", t(ioutdat),sep="")
  iout <- c(iout, paste('OUTFILE results/all_',i,'_',chrom, '_ .tbl',sep='')) #Results file name
  ioutfilename <- paste('metal_scripts/all_',i,'_',chrom,'.mi',sep='') #Metal script input file name
  write.table(c(mi_header,iout,mi_footer), ioutfilename ,quote=F,row.names=F,col.names=F)

  if(!is.na(dat[i,]$EA) | !is.na(dat[i,]$EA2) | !is.na(dat[i,]$EA3) | !is.na(dat[i,]$EA4))
  {
    eaoutdat  <- c(eaoutdat, paste(na.omit(unlist(dat[i,c("EA","EA2","EA3","EA4")])),chrom,sep="_"))
    eaout <- paste("PROCESS ",  "metal_inputs/",t(eaoutdat),sep="")

    eaout <- c(eaout, paste('OUTFILE results/eur_',i,'_',chrom,  '_ .tbl',sep=''))

    eaoutfilename <- paste('metal_scripts/eur_',i,'_',chrom, '.mi',sep='')
    write.table(c(mi_header,eaout,mi_footer), eaoutfilename ,quote=F,row.names=F,col.names=F)
  }

  if(!is.na(dat[i,]$AA) | !is.na(dat[i,]$AA2) | !is.na(dat[i,]$AA3) | !is.na(dat[i,]$AA4))
  {
    aaoutdat  <- c(aaoutdat, paste(na.omit(unlist(dat[i,c("AA","AA2","AA3","AA4")])),chrom,sep="_"))
    aaout <- paste("PROCESS ",  "metal_inputs/", t(aaoutdat),sep="")

    aaout <- c(aaout, paste('OUTFILE results/aam_',i,'_',chrom,  '_ .tbl',sep=''))

    aaoutfilename <- paste('metal_scripts/aam_',i,'_',chrom, '.mi',sep='')
    write.table(c(mi_header,aaout,mi_footer), aaoutfilename ,quote=F,row.names=F,col.names=F)
  }

  if(!is.na(dat[i,]$HNA))
  {
    hnaoutdat  <- c(hnaoutdat, paste(na.omit(unlist(dat[i,c("HNA")])),chrom,sep="_"))
    hnaout <- paste("PROCESS ",  "metal_inputs/",t(hnaoutdat),sep="")

    hnaout <- c(hnaout, paste('OUTFILE results/lat_',i,'_',chrom,  '_ .tbl',sep=''))

    hnaoutfilename <- paste('metal_scripts/lat_',i,'_',chrom, '.mi',sep='')
    write.table(c(mi_header,hnaout,mi_footer), hnaoutfilename ,quote=F,row.names=F,col.names=F)
  }

 }
}
