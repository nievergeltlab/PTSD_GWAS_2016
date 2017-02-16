
#Make a folder for all studies to be included.
mkdir starting_data

#Move all Ricopili results files into here
#...

#Make a list of all studies of unrelateds called studylist.txt
#Make a list of all studies of relateds called studylist_related.txt

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
  gzip -d -c starting_data/"$study" | awk -v chr=$chr -v maf=$maf -v info=$info'{if  (NR == 1) print "SNP","A1","A2","OR","SE","P" ;  if (($1 == chr) && ($8 > info) &&  ($6 > maf) && ($7 > maf) &&  ($6 < 1-maf) && ($7 < 1-maf))  print $2,$4,$5,$9,$10,$11}' | grep -v NA > metal/"$study"_"$chr".metal
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
  gzip -d -c starting_data/"$study" | awk -v chr=$chr -v maf=$maf  '{if  (NR == 1) print "SNP","A1","A2","OR","SE","P" ;  if (($1 ==chr) && ($7 > maf) && ($7 < 1-maf))  print $2,$5,$6,exp($8),$9,$14}' > metal/"$study"_"$chr".metal
 done
done

