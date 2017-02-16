#Make a folder for all studies to be included.
mkdir starting_data

#To the User: Move all Ricopili results files into starting_data
#To the User: Make a list of all studies of unrelateds called studylist.txt
#To the User: Make a list of all studies of relateds called studylist_related.txt

#Make a folder of where filtered data will go
mkdir metal_inputs

#From one of the datasets where all SNP positions are given, list them all SNP positions into a file
zcat starting_data/file_name_here.gz  | awk '{print  $2,$1,$3}' | LC_ALL=C sort -k1,1b > phase3.locations
LC_ALL=C sort -k1,1b phase3.locations > phase3.locations2
#Save these as an R object

#Module load R
R
 nr <- as.numeric(system(paste("wc -l phase3.locations | awk \'{print $1}\' "),intern=T))
 #get list of numeric column classes (reads file faster)
 dat_temp <- read.table("phase3.locations", header=T,nrows=50,stringsAsFactors=F)
 classes_to_use <- sapply(dat_temp , class)
 #read file
 dat1 <- read.table ("phase3.locations", header=T,stringsAsFactors=F,nrow=nr,colClasses=classes_to_use)

 snp_locations <- dat1
 
 save(snp_locations, file="phase3_snplocations_v2_oct31_2016.R")
 q("no")
 
#Export METAL columns only, with MAF and info filters
#Info filter : Filter to markers with info > 0.6. Power is approx 80% at info > 0.6, according to Zheng & Scheet 2011
#MAF filter: Filter at MAF 0.5%
#NA Fitler: Filter out markers without results
maf=0.005
info=0.6

for study in  $(cat studylist.txt )
do
 echo "Filtering $study"
 for chr in {1..22}
 do
  gzip -d -c starting_data/"$study" | awk -v chr=$chr -v maf=$maf -v info=$info'{if  (NR == 1) print "SNP","A1","A2","OR","SE","P" ;  if (($1 == chr) && ($8 > info) &&  ($6 > maf) && ($7 > maf) &&  ($6 < 1-maf) && ($7 < 1-maf))  print $2,$4,$5,$9,$10,$11}' | grep -v NA > metal/"$study"_"$chr"
 done
done

#Not complete as of Feb 16. 
#Filter GEMMA outputs (related subjects)
#Assumes that GEMMA was run in a way that prints out all p-values (wald, lrt, mle, reml, score)
for study in  $(cat studylist_related.txt )
do
 echo "Filtering $study"
 for chr in {1..22}
 do
  gzip -d -c starting_data/"$study" | awk -v chr=$chr -v maf=$maf  '{if  (NR == 1) print "SNP","A1","A2","OR","SE","P" ;  if (($1 ==chr) && ($7 > maf) && ($7 < 1-maf))  print $2,$5,$6,exp($8),$9,$14}' > metal/"$study"_"$chr"
 done
done

#To the user: prepare an excel file that is a sequence of when studies will be entered.
#Sample file:

 #Group	Order	Study	EA	AA	HNA	EA2	AA2
 #1	1	MRS	daner_mrsc_eur_analysis_run2.gz 	daner_mrsc_aam_analysis_run2.gz	daner_mrsc_lat_analysis_run2.gz #N/A	#N/A
 #1	3	SAFR (Adding African diversity)	#N/A	daner_safr_oth_analysis_run2.gz #N/A	#N/A	daner_safr_afr_analysis_run3.gz

#Now that the metal inputs have been made, we can run metal

#Now that metal has been run, we can make the sequence of manhattan plots
