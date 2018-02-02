#Read .csv file containing an ordered list of studies to be included in the meta analysis
 args <- commandArgs(trailingOnly = TRUE)
 dataset <- args[1]
 mscriptdir <- args[2]
 indir <- args[3]
 outdir <- args[4]
 suff <- args[5]

 dat <- read.csv(dataset,header=T,na.strings=c("NA","#N/A"),stringsAsFactors=F)
 mi_header <- read.csv('meta_header.mi',header=F,stringsAsFactors=F)$V1
 mi_footer <- read.csv('meta_footer.mi',header=F,stringsAsFactors=F)$V1

#Example data head
#Make sure that COMPLETE file names are there (i.e. include file extensions in the names!)

#Group	ORder	Study	EA	AA	HNA	EA2	AA2
#1	1	MRS	daner_mrsc_eur_analysis_run2.gz 	daner_mrsc_aam_analysis_run2.gz	daner_mrsc_lat_analysis_run2.gz	#N/A	#N/A
#1	2	ONG (More military)	daner_onga_eur_analysis_run2.gz	#N/A	#N/A	#N/A	#N/A
#1	3	SAFR (Adding African diversity)	#N/A	daner_safr_oth_analysis_run2.gz	#N/A	#N/A	daner_safr_afr_analysis_run3.gz
#1	6	DNHS (more AA)	#N/A	daner_dnhs_aam_analysis_run3.gz	#N/A	#N/A	#N/A

#Note I use AA2 and EA2 because sometimes I added two datasets at once, e.g. if one study has 2 components.

for (chrom in c(1:22))
{

 for (anc in c("all","eur","aam","lat"))
 {
  if (anc == "all")
  {
   anccols=c(4:14)
  } else if (anc == "eur")
  {
   anccols=c("EA","EA2","EA3","EA4")
  } else if (anc == "aam")
  {
   anccols=c("AA","AA2","AA3","AA4")
  } else if (anc == "lat")
  {
   anccols=c("HNA")
  } else {
  print( "no valid ancestry set!")
  }
  
  ioutdat  <- c(paste(na.omit(unlist(dat[,anccols])),chrom,sep="_")) 
  iout <- paste("PROCESS ", indir,  "/", t(ioutdat),sep="")
  #Results file name will be called this
   iout <- c(iout, paste('OUTFILE ', outdir, '/',anc,"_",suff,'_',chrom, '_ .tbl',sep='')) 
  #Metal script input file name will be called this:
   ioutfilename <- paste(mscriptdir, "/", anc, "_",suff,'_',chrom,'.mi',sep='')        
  #Write results
   write.table(c(mi_header,iout,mi_footer), ioutfilename ,quote=F,row.names=F,col.names=F)
    

 }  
}

