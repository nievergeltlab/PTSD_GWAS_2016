workingdir=/home/cnieverg/report_data/testrun

#To the user: Install metal, write the path to the binary here
metalpath=/home/cnieverg/meta_istss/generic-metal/metal


cd $workingdir

#Make a folder for all studies to be included.
mkdir starting_data

#Make a folder of where filtered data will go
mkdir metal_inputs

#Make a folder of where the metal scripts will go
mkdir metal_scripts

mkdir results
mkdir plots
mkdir temporary_files
mkdir errandout

#To the User: Move all Ricopili results files into starting_data
#To the User: Make a list of all studies of unrelateds called studylist.txt
#To the User: Make a list of all studies of relateds called studylist_related.txt


#Unzip the datasets, e..g
#ls *.tar > tarfiles.txt

#for file in $(cat tarfiles.txt)
#do
# tar xvf "$file" --exclude '*males*' --exclude '*females*' --wildcards 
#done

#Get MF
ls *.tar > tarfiles.txt

for file in $(cat tarfiles.txt)
do
 tar xvf "$file" --strip-components=4   --wildcards '*males*' '*females*'
 rm $file
done

#From one of the datasets where all SNP positions are given, list them all SNP positions into a file
#To the user: give a file name
zcat starting_data/daner_gtpc_aam_analysis_run3.gz  | awk '{print  $2,$1,$3}' | LC_ALL=C sort -k1,1b > phase3.locations
LC_ALL=C sort -k1,1b phase3.locations > phase3.locations2
#Save these as an R object

#module load R
R
 nr <- as.numeric(system(paste("wc -l phase3.locations | awk \'{print $1}\' "),intern=T))
 #get list of numeric column classes (reads file faster)
 dat_temp <- read.table("phase3.locations", header=T,nrows=50,stringsAsFactors=F)
 classes_to_use <- sapply(dat_temp , class)
 #read file
 dat1 <- read.table ("phase3.locations", header=T,stringsAsFactors=F,nrow=nr,colClasses=classes_to_use)

 snp_locations <- dat1
workingdir=/home/cnieverg/report_data/testrun

#To the user: Install metal, write the path to the binary here
metalpath=/home/cnieverg/meta_istss/generic-metal/metal


cd $workingdir

#Make a folder for all studies to be included.
mkdir starting_data

#Make a folder of where filtered data will go
mkdir metal_inputs

#Make a folder of where the metal scripts will go
mkdir metal_scripts

mkdir results
mkdir plots
mkdir temporary_files
mkdir errandout

#To the User: Move all Ricopili results files into starting_data
#To the User: Make a list of all studies of unrelateds called studylist.txt
#To the User: Make a list of all studies of relateds called studylist_related.txt


#Unzip the datasets, e..g
#ls *.tar > tarfiles.txt

#for file in $(cat tarfiles.txt)
#do
# tar xvf "$file" --exclude '*males*' --exclude '*females*' --wildcards 
#done

#Get MF
ls *.tar > tarfiles.txt

for file in $(cat tarfiles.txt)
do
 tar xvf "$file" --strip-components=4   --wildcards '*males*' '*females*'
 rm $file
done

#From one of the datasets where all SNP positions are given, list them all SNP positions into a file
#To the user: give a file name
zcat starting_data/daner_gtpc_aam_analysis_run3.gz  | awk '{print  $2,$1,$3}' | LC_ALL=C sort -k1,1b > phase3.locations
LC_ALL=C sort -k1,1b phase3.locations > phase3.locations2
#Save these as an R object

#module load R
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
maf=0.01
info=0.6


ls starting_data | grep -v qimr*  | grep -v vetsa* | grep -v -w old | grep -v Metadata_ > studylist.txt


#If you missed a study, create a list like this:
#echo daner_comc_aam_analysis_run3.gz  daner_meg2_af3_analysis_run5.gz | sed 's/ /\n/g' > studylist.txt

njobs=$(wc -l studylist.txt | awk '{print $1}')

qsub -t2-$njobs -l walltime=00:05:00 split_files.qsub -d $workingdir -e errandout/ -o errandout/ -F "-f studylist.txt -s starting_data -o metal_inputs -m 0.01 -i 0.6"


#Filter GEMMA outputs (related subjects)
#Assumes that GEMMA was run in a way that prints out all p-values (wald, lrt, mle, reml, score)
#Also need to fix the beta values to put them on a logit scale,recalculate SE

#awk '{print $1,$18}' starting_data/Metadata_allchr.txt | LC_ALL=C sort -k1,1b > starting_data/Metadata_allchr_qimr.txt

#LC_ALL=C sort -k2,2b <(zcat starting_data/qimr_rp_v1_oct6_2016.gz) > starting_data/qimr_rp_v1_oct6_2016.sorted

#LC_ALL=C join -1 1 -2 2 starting_data/Metadata_allchr_qimr.txt starting_data/qimr_rp_v1_oct6_2016.sorted | cut -d " " -f 2- | awk '{if ($1 != ".") print}' | sort -g -k 10 | gzip > starting_data/qimr_rp_v1_oct6_2016.fixed.gz 


LC_ALL=C sort -k2,2b <(zcat PGC_PTSD_allchrfemale_lmm_with5PCs.assoc.txt.gz) > PGC_PTSD_allchrfemale_lmm_with5PCs.assoc.txt.sorted
LC_ALL=C join -1 1 -2 2 Metadata_allchr_qimr.txt  PGC_PTSD_allchrfemale_lmm_with5PCs.assoc.txt.sorted | cut -d " " -f 2- | awk '{if ($1 != ".") print}' | sort -g -k 10 | gzip > qimr_females_rp_v1_oct6_2016.fixed.gz 


LC_ALL=C sort -k2,2b <(zcat PGC_PTSD_allchrmale_lmm_with5PCs.assoc.txt.gz) > PGC_PTSD_allchrmale_lmm_with5PCs.assoc.txt.sorted
LC_ALL=C join -1 1 -2 2 Metadata_allchr_qimr.txt  PGC_PTSD_allchrmale_lmm_with5PCs.assoc.txt.sorted | cut -d " " -f 2- | awk '{if ($1 != ".") print}' | sort -g -k 10 | gzip > qimr_males_rp_v1_oct6_2016.fixed.gz 



#20.6% ptsd prevalence



#I went on to an external server and put the QIMR on the logit scale

#Fix the queensland data
#columns 15 and 16 are the LOG SCALE BETA and SE. columns 8 and 9 are the typical scale 
for study in qimr_males_rp_v1_oct6_2016.fixed.logit.gz # $(cat studylist_related.txt )
do
 echo "Filtering $study"
 gzip -d -c starting_data/"$study" >  metal_inputs/"$study".unzip
 for chr in {1..22}
 do
  awk -v chr=$chr -v maf=$maf  '{if  (NR == 1) print "SNP","A1","A2","OR","SE","P" ;  if (($2 ==chr) && ($7 > maf) && ($7 < 1-maf) && ($15 > -5) && ($15 < 5))  print $1,$5,$6,exp($15),$16,$12}'  metal_inputs/"$study".unzip > metal_inputs/"$study"_"$chr"
 done
 rm -f  metal_inputs/"$study".unzip
done

#Fix the VETSA data, it
#columns 13 and 14 are the LOG SCALE BETA and SE. 
#added a special thing in case of horrible estimates ( do not PRINT lines where beta on hte log scale is abs (b) > 5
for study in  vets_rp_v3_april21_2017.gemma.logit.gz #  $(cat studylist.txt )
do
 echo "Filtering $study"
 gzip -d -c starting_data/"$study" > metal_inputs/"$study".unzip
 for chr in {1..22}
 do
  awk -v chr=$chr -v maf=$maf -v info=$info '{if  (NR == 1) print "SNP","A1","A2","OR","SE","P" ;  if (($2 == chr) && ($6 > info) &&  ($9 > maf) && ($9 < 1-maf) && ($12 != "NA") && ($13 < 5) && ($13 > -5) )  print $1,$7,$8,exp($13),$14,$12}' metal_inputs/"$study".unzip  > metal_inputs/"$study"_"$chr"
 done
 rm -f metal_inputs/"$study".unzip
done



#Run the script on the data to make the sequence of metal input files.
#Requires as an input an ordered list of studies called study_order.csv. Read the .R file for a description of how the input should look
#Note: The very last file generated by the script is the one you want to use if you want to do a total meta analysis instead of a growing sequence of them.
#e.g. if you have an input of 3 files, there is a file called all_

Rscript make_mi_files.R
Rscript make_mi_files_MIL.R
Rscript make_mi_files_CIV.R
Rscript make_mi_files_FEM.R
Rscript make_mi_files_FEM2.R
Rscript make_mi_files_FEM3.R

Rscript make_mi_files_MAL.R
Rscript make_mi_files_NOAS.R

Rscript make_mi_files_NOVS.R

Rscript make_mi_files_LOO.R
Rscript make_mi_files_LOOF.R
Rscript make_mi_files_LOOM.R

Rscript make_mi_files_SOBP.R

pgc_ptsd_study_order_v4_sobp.csv

###Run metal

#First list all .mi inputs, then write a loop sending them to metal 
ls metal_scripts | grep CIV | grep CIVeur_18 > metal_jobs3.txt #Put in a grep for the last row of the study input data if you ONLY want the total data meta analysis

ls metal_scripts | grep FEM2 | grep FEMeur_47 > metal_jobs4.txt #Put in a grep for the last row of the study input data if you ONLY want the total data meta analysis
ls metal_scripts | grep MAL | grep MALeur_39 > metal_jobs5.txt #Put in a grep for the last row of the study input data if you ONLY want the total data meta analysis

ls metal_scripts | grep FEM2 | grep FEM2eur_47 > metal_jobs6.txt #Put in a grep for the last row of the study input data if you ONLY want the total data meta analysis

ls metal_scripts | grep FEM3 | grep FEM3eur_47 > metal_jobs7.txt #Put in a grep for the last row of the study input data if you ONLY want the total data meta analysis
ls metal_scripts | grep NOAS | grep NOASeur_39 > metal_jobs8.txt #Put in a grep for the last row of the study input data if you ONLY want the total data meta analysis
ls metal_scripts | grep SOBP  > metal_jobs10.txt #Put in a grep for the last row of the study input data if you ONLY want the total data meta analysis

ls metal_scripts | grep LOOM  > metal_jobs11.txt #Put in a grep for the last row of the study input data if you ONLY want the total data meta analysis



qsub -t1 -l walltime=00:05:00 run_meta_v2.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m metal_jobs11.txt -p $metalpath"

ls metal_scripts | grep LOOM  > metal_jobs12.txt #Put in a grep for the last row of the study input data if you ONLY want the total data meta analysis

ncommands=$(wc -l metal_jobs12.txt | awk '{print $1}')
nodesize=16
nodeuse=$(($nodesize ))
totjobs=$(( ($ncommands + $nodeuse - 1 ) / $nodeuse ))

qsub -t1-$totjobs -l walltime=00:05:00 run_meta_v2_loo.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m metal_jobs12.txt -p $metalpath -n $nodeuse"

#Concatenate over all chromosomes
 cat results/MILeur_23_*.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/MILeur_23.results
  cat results/CIVeur_18_*.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/CIVeur_18.results
 
 cat results/FEMeur_47_*.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/FEMeur_47.results

 cat results/FEM2eur_47_*.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/FEM2eur_47.results
 cat results/FEM3eur_47_*.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/FEM3eur_47.results
  
 cat results/MALeur_39_*.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/MALeur_39.results
 cat results/NOASeur_39_*.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/NOASeur_39.results
 cat results/NOVSeur_39_*.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/NOVSeur_39.results

 for i in {1..42}
 do
  cat results/LOOeur_"$i"_*_1.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/LOOeur_"$i"_.results
 done

 cat results/eur_39_*.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/eur_39.results
 cat results/all_39_*.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/all_39.results
 cat results/aam_39_*.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/aam_39.results
 cat results/lat_38_*.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/lat_38.results

sort -g -k 6 results_cat/eur_39.results | head -n 30 > results_cat/eur_39.results.top30
sort -g -k 6 results_cat/aam_39.results | head -n 30 > results_cat/aam_39.results.top30
sort -g -k 6 results_cat/lat_38.results | head -n 30 > results_cat/lat_38.results.top30
sort -g -k 6 results_cat/all_39.results | head -n 30 > results_cat/all_39.results.top30




 #Manhattan plot of results (list of things to plot is the input)


 ls results_cat | grep -E eur | grep SOBP > metal_outputsA.txt
 ls results_cat | grep -E aam | grep  SOBP > metal_outputsB.txt
 ls results_cat | grep -E all | grep SOBP > metal_outputsC.txt
 ls results_cat | grep -E lat | grep SOBP > metal_outputsD.txt


 njobs=$(wc -l metal_outputsC.txt | awk '{print $1}')
 njobs=$(wc -l metal_outputsD.txt | awk '{print $1}')

 qsub -t1 -l walltime=00:20:00 mh_plot.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m mh_plot_pgc.R -s ManhattanPlotterFunction_colorfixed_max10ylim2_pgc.R -l phase3.locations2 -d metal_outputsA.txt -c green -p 0.05"
  njobs=$(wc -l metal_outputsA.txt | awk '{print $1}')
 qsub -t1 -l walltime=00:20:00 mh_plot.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m mh_plot_pgc.R -s ManhattanPlotterFunction_colorfixed_max10ylim2_pgc.R -l phase3.locations2 -d metal_outputsB.txt -c red -p 0.05"
  njobs=$(wc -l metal_outputsB.txt | awk '{print $1}')


 
 qsub -t1 -l walltime=00:20:00 mh_plot.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m mh_plot_pgc.R -s ManhattanPlotterFunction_colorfixed_max10ylim2_pgc.R -l phase3.locations2 -d metal_outputsC.txt -c blue -p 0.05"
 qsub -t1 -l walltime=00:20:00 mh_plot.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m mh_plot_pgc.R -s ManhattanPlotterFunction_colorfixed_max10ylim2_pgc.R -l phase3.locations2 -d metal_outputsD.txt -c purple -p 0.05"


ncommands=$(wc -l metal_outputsD.txt | awk '{print $1}')
nodesize=16
nodeuse=$(($nodesize ))
totjobs=$(( ($ncommands + $nodeuse - 1 ) / $nodeuse ))

 qsub -t1-$totjobs -l walltime=00:20:00 mh_plot_rapid.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m mh_plot_pgc.R -s ManhattanPlotterFunction_colorfixed_max10ylim2_pgc.R -l phase3.locations2 -d metal_outputsA.txt -c green -p 0.05"
 
 qsub -t1-$totjobs -l walltime=00:20:00 mh_plot_rapid.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m mh_plot_pgc.R -s ManhattanPlotterFunction_colorfixed_max10ylim2_pgc.R -l phase3.locations2 -d metal_outputsB.txt -c red -p 0.05"
 qsub -t1-$totjobs -l walltime=00:20:00 mh_plot_rapid.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m mh_plot_pgc.R -s ManhattanPlotterFunction_colorfixed_max10ylim2_pgc.R -l phase3.locations2 -d metal_outputsC.txt -c blue -p 0.05"
 qsub -t1-$totjobs -l walltime=00:20:00 mh_plot_rapid.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m mh_plot_pgc.R -s ManhattanPlotterFunction_colorfixed_max10ylim2_pgc.R -l phase3.locations2 -d metal_outputsD.txt -c purple -p 0.05"
 

 echo SORTED_PTSD_EA9_ALL_study_specific_PCs1.txt > metal_outputsL.txt


 qsub -t1 -l walltime=00:20:00 mh_plot.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m mh_plot_pgc.R -s ManhattanPlotterFunction_colorfixed_max10ylim2_pgc.R -l phase3.locations2 -d metal_outputsL.txt -c green -p 0.05"




ls results_cat | grep -E lat |  grep 38 > metal_outputsD.txt
 qsub -t1 -l walltime=00:20:00 make_locuszoom.qsub -d $workingdir -e errandout/ -o errandout/ -F "-l phase3.locations2 -d metal_outputsD.txt"


#Produce QQ plots for results
 awk '{print $6}'  results_cat/eur_39.results | grep -v NA > results_cat/eur_39.results.p
 awk '{print $6}'  results_cat/all_39.results | grep -v NA > results_cat/all_39.results.p
 awk '{print $6}'  results_cat/aam_39.results | grep -v NA > results_cat/aam_39.results.p
 awk '{print $6}'  results_cat/lat_38.results | grep -v NA > results_cat/lat_38.results.p

 ls results_cat | grep .p$ > qq_plot.files
 qsub -l walltime=00:20:00 qq_plot.qsub -d $workingdir -e errandout/ -o errandout/ -F "-q qq_plot.R -e 1 -p qq_plot.files "




#Make a forest plot
forestsnp=rs57753395
forestsnpchr=13

#Input file is a .csv of rows with: chr,SNP,ancestry


#Get genotype data for forest plot SNPS
for fp in $(tail -n1 fucker.csv)
do
 forestsnp=$(echo $fp | awk 'BEGIN {FS=","}{print $2}')
 forestsnpchr=$(echo $fp | awk 'BEGIN {FS=","}{print $1}') 
 
 echo $forestsnp $forestsnpchr
 grep -w -m1 $forestsnp metal_inputs/*_"$forestsnpchr" | sed 's/:/ /g' | sed 's/metal_inputs\///g' | sed 's/\.gz_[0-9]*/.gz /g' > results_cat/"$forestsnp".use
done

#Plot forest plots
for fp in $(cat fucker.csv | awk 'NR==5{print}')
do
 forestsnp=$(echo $fp | awk 'BEGIN {FS=","}{print $2}')
 ancgroup=$(echo $fp | awk 'BEGIN {FS=","}{print $3}') 
 echo $forestsnp $ancgroup
 Rscript forest_plot.R results_cat/"$forestsnp".use pgc_ptsd_study_order_v4_descriptors.csv $ancgroup "$forestsnp".plot
done



#Look up genotyping quality of top hits
13 rs57753395 eur  # eur tophit 2. maf 5%, near LOC105370330. Can't find proxy directly (phase 3). The second hit under this was rs73553233. Proxy is rs11618293, Is genotyped on the psych array, most illumina, but not psycchip. Results for this SNP are p ~ 0.0002,error bars are looser for it, effects similar
6 rs34582172 eur # eur tophit 1. maf 40%, is intergenic. Proxy is likely rs4612196, which is still GWS. 
9 rs9410042  lat # lat tophit 1: overall maf 13% (5% in ea and lat, 40% in afr) , on ehmt 1. proxy is rs3125782, its near GWS
18 rs3862732 lat # lat tophit 2: maf is 25% in ea/lat. is likely genotyped
13 rs115539978 aam # aam tophit. intergenic, closest is mir 5007 . rs9536997 is i think I supporting it, however this marker is af~40%,  yet hte top hit is af~=2%. rs12862673 and rs9569233 rs4383040 are near by genotyped. theyre p values are totally insig

ls starting_data > all_input_data.txt
for i in $(cat all_input_data.txt)
do
 echo "Checking file $i"
 zgrep -w -m1 rs4612196 starting_data/"$i" >> check/rs4612196.assoc.dosage
done
#Conclusion: genotyped proxies will differ by study. I think I should use the info metric on the SNP itself.

#look up genotyped proxy marker and plot that

#ls starting_data | grep -v qimr*  | grep -v vetsa* | grep -v -w old | grep -v Metadata_ > studylist.txt
#njobs=$(wc -l studylist.txt | awk '{print $1}')
#qsub -t1-$njobs -l walltime=00:05:00 split_files_infoscores.qsub -d $workingdir -e errandout/ -o errandout/ -F "-f studylist.txt -s starting_data -o snp_info -m 0.01 -i 0.6"
#I don't have necessary columns for vetsa and qimr, leaving out for now...

echo 6,rs4612196,eur > rs4612196.csv
for fp in $(cat rs4612196.csv )
do
 forestsnp=$(echo $fp | awk 'BEGIN {FS=","}{print $2}')
 forestsnpchr=$(echo $fp | awk 'BEGIN {FS=","}{print $1}') 
 
 echo $forestsnp $forestsnpchr
 zgrep -w -m1 $forestsnp snp_info/*_"$forestsnpchr".gz | sed 's/:/ /g' | sed 's/snp_info\///g' | sed 's/\.gz_[0-9]*//g' > results_cat/"$forestsnp".info
done


zgrep -w -m1 rs11618293 snp_info/*_"$forestsnpchr".gz | sed 's/:/ /g' | sed 's/snp_info\///g' | sed 's/\.gz_[0-9]*//g' > results_cat/rs11618293.info



#Add civilian, military, etc. subgrouping


#Make an LDSC file
grep -v ??????? results_cat/eur_39.results  | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"12813","35640"}' | zip  > eur_39_moststudies.results.zip

grep -v ??? results_cat/MILeur_23.results  | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"6086","21545"}' | zip  > MILeur_23_moststudies.results.zip
grep -v ??? results_cat/CIVeur_18.results  | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"6208","12612"}' | zip  > CIVeur_18_moststudies.results.zip


munge_sumstats.py --sumstats MILeur_23_moststudies.results --out MILeur_23_moststudies_munge
munge_sumstats.py --sumstats CIVeur_18_moststudies.results --out CIVeur_23_moststudies_munge

ldsc.py \
--rg MILeur_23_moststudies_munge.sumstats.gz,CIVeur_23_moststudies_munge.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out mil_civ
less mil_civ.log


grep -v ? results_cat/eur_39.results  | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"12813","35640"}' | zip  > eur_39_allstudies.results.zip

grep -v ? results_cat/eur_39.results  | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"12813","35640"}' | zip  > eur_39_allstudies.results.zip
grep -v ?? results_cat/FEMeur_47.results | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"5914","8759"}' | zip  > FEMeur_47_allstudies.results.zip
grep -v ?? results_cat/MALeur_39.results | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"6381","23580"}' | zip  > MALeur_39_allstudies.results.zip


grep -v ?? results_cat/FEM2eur_47.results | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"5315","8008"}'  > FEM2eur_47_allstudies.results

grep -v ?? results_cat/FEM3eur_47.results | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"3488","7011"}' | zip > FEM3eur_47_allstudies.results.zip


awk 'BEGIN {FS="?"; var=0 } {if (NF>0){ var=var + (NF-1) } } END{print var}' testfile1.txt


grep -v ?? results_cat/NOASeur_39.results | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"9963","26810"}'   > NOASeur_39_allstudies.results.zip

grep -v ?? results_cat/NOVSeur_39.results | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"12000","34000"}'   > NOVSeur_39_allstudies.results.zip

results_cat/NOASeur_39.results
 
unzip MALeur_39_allstudies.results.zip
unzip FEMeur_47_allstudies.results.zip

#All studies

awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"12813","35640"}'  results_cat/eur_39.results  > eur_39_allstudies.results
awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"5914","8759"}' results_cat/FEMeur_47.results  > FEMeur_47_allstudies.results
awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"6381","23580"}' results_cat/MALeur_39.results > MALeur_39_allstudies.results

munge_sumstats.py --sumstats eur_39_allstudies.results --out results_cat/eur_39_allstudies_munge
munge_sumstats.py --sumstats FEMeur_47_allstudies.results --out results_cat/FEMeur_47_allstudies_munge
munge_sumstats.py --sumstats MALeur_39_allstudies.results --out results_cat/MALeur_39_allstudies_munge


cd results_cat/

ldsc.py \
--h2 eur_39_allstudies_munge.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--samp-prev 0.264 \
--pop-prev 0.09 \
--out eur38_h2

ldsc.py \
--h2 FEMeur_47_allstudies_munge.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--samp-prev 0.40 \
--pop-prev 0.11 \
--out femeur47_h2

ldsc.py \
--h2 MALeur_39_allstudies_munge.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--samp-prev 0.21 \
--pop-prev 0.05 \
--out maleur38_h2

#Most studies ( to match LOO


grep -v ??  results_cat/eur_39.results  | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"12813","35640"}'  > eur_39_moststudies.results
grep -v ??  results_cat/FEMeur_47.results | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"5914","8759"}'   > FEMeur_47_moststudies.results
grep -v ?? results_cat/MALeur_39.results | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"6381","23580"}'   > MALeur_39_moststudies.results

munge_sumstats.py --sumstats eur_39_moststudies.results --out results_cat/eur_39_moststudies_munge
munge_sumstats.py --sumstats FEMeur_47_moststudies.results --out results_cat/FEMeur_47_moststudies_munge
munge_sumstats.py --sumstats MALeur_39_moststudies.results --out results_cat/MALeur_39_moststudies_munge


cd results_cat/

ldsc.py \
--h2 eur_39_moststudies_munge.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--samp-prev 0.264 \
--pop-prev 0.09 \
--out eur39_h2_most

ldsc.py \
--h2 FEMeur_47_moststudies_munge.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--samp-prev 0.40 \
--pop-prev 0.11 \
--out femeur47_h2_most

ldsc.py \
--h2 MALeur_39_moststudies_munge.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--samp-prev 0.21 \
--pop-prev 0.05 \
--out maleur39_h2_most



#h2 is lower if removing small studies

##LOO - fix so that it doesnt do 50 of these
ls results_cat | grep LOOM > ldsc_jobs15.txt #Put in a grep for the last row of the study input data if you ONLY want the total data meta analysis
njobs=$(wc -l ldsc_jobs15.txt | awk '{print $1}')
qsub -t1-$njobs -l walltime=00:05:00 ldsc.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m ldsc_jobs15.txt -p results_cat"


grep -v ??? results_cat/eur_39.results | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6, "1556","6618" }'  | zip  > eur_9_allstudies.results.zip

grep -v ?.*?.* results_cat/all_39.results | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6, "19000","50000" }'  | zip  > all_39_allstudies.results.zip

 for i in {1..33}
 do
  cat results/LOOMeur_"$i"_*_1.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/LOOMeur_"$i"_.results

  #rm -f LOOeur_"$i"_allstudies.results
  #rm -f results_cat/LOOeur_"$i"_.results
 done


grep #Look up genotyping quality of top hits
grep -E "rs57753395|rs34582172|rs9410042|rs3862732|rs115539978" /home/cnieverg/report_data/testrun/results_cat/eur_39.results
grep -E "rs57753395|rs34582172|rs9410042|rs3862732|rs115539978" /home/cnieverg/report_data/testrun/results_cat/aam_39.results
grep -E "rs57753395|rs34582172|rs9410042|rs3862732|rs115539978" /home/cnieverg/report_data/testrun/results_cat/lat_38.results
 
#Add positions to data
 ls results_cat | grep -E eur_39\|lat_38\|aam_39 | grep -v NOVS | grep -v MAL | grep -v NOVS |  grep .results$  | grep lat > metal_outputsA.txt

 qsub -t1 -l walltime=00:20:00 add_positions.qsub -d $workingdir -e errandout/ -o errandout/ -F "-l phase3.locations2 -d metal_outputsA.txt "



bash add_positions.qsub -l phase3.locations2 -d results_cat/lat_38.results

#Filter to select snps and do ldsc

 maf=0.01
 info=0.9
 zcat  /home/cnieverg/report_data/testrun/starting_data/daner_mrsc_eur_analysis_run3_males.gz | awk  -v maf=$maf -v info=$info '{if  (NR == 1) print "SNP";  if (($8 > info) &&  ($7 > maf) && ($7 < 1-maf)  ) print $2}' | LC_ALL=C sort -k1b,1 > good_snps.snplist

awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"12813","35640"}'  results_cat/eur_39.results | LC_ALL=C sort -k1b,1 > eur_39.results.munge.sorted
LC_ALL=C join good_snps.snplist eur_39.results.munge.sorted | sort -g -k6 > eur_39.results.munge.sorted.fin
munge_sumstats.py --sumstats eur_39.results.munge.sorted.fin --out results_cat/eur_39_allstudies_goodsnpsonly_munge
cd results_cat/

ldsc.py \
--h2 eur_39_allstudies_goodsnpsonly_munge.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--samp-prev 0.264 \
--pop-prev 0.09 \
--out eur38_goodsnps_h2

#Make h2 results into a LZ plot
grep  Liability results_cat/LOOeur_*_.results_h2.log | awk 'BEGIN{FS=":"; OFS=" "}{print $1,$3}' | sed 's/results_cat\/LOOeur_//g' | sed 's/_.results_h2.log//g' | sed 's/(//g' | sed 's/)//g' |  > results_cat/looeur_h2.ldsc_res
grep  Liability results_cat/LOOFeur_*_.results_h2.log | awk 'BEGIN{FS=":"; OFS=" "}{print $1,$3}' | sed 's/results_cat\/LOOFeur_//g' | sed 's/_.results_h2.log//g' | sed 's/(//g' | sed 's/)//g' > results_cat/loofeur_h2.ldsc_res
grep  Liability results_cat/LOOMeur_*_.results_h2.log | awk 'BEGIN{FS=":"; OFS=" "}{print $1,$3}' | sed 's/results_cat\/LOOMeur_//g' | sed 's/_.results_h2.log//g' | sed 's/(//g' | sed 's/)//g' > results_cat/loomeur_h2.ldsc_res


R
 library(metafor)
 study_abbrs <- read.csv('pgc_ptsd_study_order_v6.csv',stringsAsFactors=F)
 
#All studies
 liability <- read.table('results_cat/looeur_h2.ldsc_res',header=F)
 names(liability) <- c("study_ord","h2","se")

 studies <- read.table('ea_studies_list.txt')
 studies$study <- c(1:dim(studies)[1])
 names(studies) <- c("study","study_ord")
 studies$study <- gsub("_22","",studies$study)


 named_studies <- merge(studies,study_abbrs,by="study",all.x=T)
 dat <- merge(liability,named_studies,by="study_ord",all.x=T)
 dat <- dat[order(dat$cases,decreasing=T),]
 meta_res <- rma(yi=dat$h2,sei=dat$se,method="FE",slab=dat$abbr)

 pdf("plots/LOO_forest.pdf",7,9)
 forest(meta_res)
 dev.off()

#Fem studies
 liability <- read.table('results_cat/loofeur_h2.ldsc_res',header=F)
 names(liability) <- c("study_ord","h2","se")

 studies <- read.table('FEM_ea_studies_list.txt')
 studies$study <- c(1:dim(studies)[1])
 names(studies) <- c("study_females","study_ord")
 studies$study_females <- gsub("_22","",studies$study_females)


 named_studies <- merge(studies,study_abbrs,by="study_females",all.x=T)
 dat <- merge(liability,named_studies,by="study_ord",all.x=T)
 dat <- dat[order(dat$cases,decreasing=T),]
 meta_res <- rma(yi=dat$h2,sei=dat$se,method="FE",slab=dat$abbr)

 pdf("plots/LOOF_forest.pdf",7,9)
 forest(meta_res)
 dev.off()




R
 library(metafor)
 study_abbrs <- read.csv('pgc_ptsd_study_order_v6.csv',stringsAsFactors=F)
 
#All studies
 liability <- read.table('results_cat/looeur_h2.ldsc_res',header=F)
 names(liability) <- c("study_ord","h2","se")

 studies <- read.table('ea_studies_list.txt')
 studies$study <- c(1:dim(studies)[1])
 names(studies) <- c("study","study_ord")
 studies$study <- gsub("_22","",studies$study)


 named_studies <- merge(studies,study_abbrs,by="study",all.x=T)
 dat <- merge(liability,named_studies,by="study_ord",all.x=T)
 dat <- dat[order(dat$cases,decreasing=T),]
 meta_res <- rma(yi=dat$h2,sei=dat$se,method="FE",slab=dat$abbr)

 pdf("plots/LOO_forest.pdf",7,9)
 forest(meta_res)
 dev.off()

#male studies
 liability <- read.table('results_cat/loomeur_h2.ldsc_res',header=F)
 names(liability) <- c("study_ord","h2","se")

 studies <- read.table('MALea_studies_list.txt')
 studies$study <- c(1:dim(studies)[1])
 names(studies) <- c("study_males","study_ord")
 studies$study_males <- gsub("_22","",studies$study_males)


 named_studies <- merge(studies,study_abbrs,by="study_males",all.x=T)
 dat <- merge(liability,named_studies,by="study_ord",all.x=T)
 dat <- dat[order(dat$cases,decreasing=T),]
 dat <- subset(dat,included_male == 1)
 meta_res <- rma(yi=dat$h2,sei=dat$se,method="FE",slab=dat$abbr)

 pdf("plots/LOOM_forest.pdf",7,9)
 forest(meta_res)
 dev.off()




 



done

#Fix the VETSA data, it
#columns 13 and 14 are the LOG SCALE BETA and SE. 
#added a special thing in case of horrible estimates ( do not PRINT lines where beta on hte log scale is abs (b) > 5
for study in  vets_rp_v3_april21_2017.gemma.logit.gz #  $(cat studylist.txt )
do
 echo "Filtering $study"
 gzip -d -c starting_data/"$study" > metal_inputs/"$study".unzip
 for chr in {1..22}
 do
  awk -v chr=$chr -v maf=$maf -v info=$info '{if  (NR == 1) print "SNP","A1","A2","OR","SE","P" ;  if (($2 == chr) && ($6 > info) &&  ($9 > maf) && ($9 < 1-maf) && ($12 != "NA") && ($13 < 5) && ($13 > -5) )  print $1,$7,$8,exp($13),$14,$12}' metal_inputs/"$study".unzip  > metal_inputs/"$study"_"$chr"
 done
 rm -f metal_inputs/"$study".unzip
done



#Run the script on the data to make the sequence of metal input files.
#Requires as an input an ordered list of studies called study_order.csv. Read the .R file for a description of how the input should look
#Note: The very last file generated by the script is the one you want to use if you want to do a total meta analysis instead of a growing sequence of them.
#e.g. if you have an input of 3 files, there is a file called all_

Rscript make_mi_files.R
Rscript make_mi_files_MIL.R
Rscript make_mi_files_CIV.R
Rscript make_mi_files_FEM.R
Rscript make_mi_files_FEM2.R
Rscript make_mi_files_FEM3.R

Rscript make_mi_files_MAL.R
Rscript make_mi_files_NOAS.R

Rscript make_mi_files_NOVS.R

Rscript make_mi_files_LOO.R
Rscript make_mi_files_LOOF.R
Rscript make_mi_files_LOOM.R

Rscript make_mi_files_SOBP.R

pgc_ptsd_study_order_v4_sobp.csv

###Run metal

#First list all .mi inputs, then write a loop sending them to metal 
ls metal_scripts | grep CIV | grep CIVeur_18 > metal_jobs3.txt #Put in a grep for the last row of the study input data if you ONLY want the total data meta analysis

ls metal_scripts | grep FEM2 | grep FEMeur_47 > metal_jobs4.txt #Put in a grep for the last row of the study input data if you ONLY want the total data meta analysis
ls metal_scripts | grep MAL | grep MALeur_39 > metal_jobs5.txt #Put in a grep for the last row of the study input data if you ONLY want the total data meta analysis

ls metal_scripts | grep FEM2 | grep FEM2eur_47 > metal_jobs6.txt #Put in a grep for the last row of the study input data if you ONLY want the total data meta analysis

ls metal_scripts | grep FEM3 | grep FEM3eur_47 > metal_jobs7.txt #Put in a grep for the last row of the study input data if you ONLY want the total data meta analysis
ls metal_scripts | grep NOAS | grep NOASeur_39 > metal_jobs8.txt #Put in a grep for the last row of the study input data if you ONLY want the total data meta analysis
ls metal_scripts | grep SOBP  > metal_jobs10.txt #Put in a grep for the last row of the study input data if you ONLY want the total data meta analysis

ls metal_scripts | grep LOOM  > metal_jobs11.txt #Put in a grep for the last row of the study input data if you ONLY want the total data meta analysis



qsub -t1 -l walltime=00:05:00 run_meta_v2.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m metal_jobs11.txt -p $metalpath"

ls metal_scripts | grep LOOM  > metal_jobs12.txt #Put in a grep for the last row of the study input data if you ONLY want the total data meta analysis

ncommands=$(wc -l metal_jobs12.txt | awk '{print $1}')
nodesize=16
nodeuse=$(($nodesize ))
totjobs=$(( ($ncommands + $nodeuse - 1 ) / $nodeuse ))

qsub -t1-$totjobs -l walltime=00:05:00 run_meta_v2_loo.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m metal_jobs12.txt -p $metalpath -n $nodeuse"

#Concatenate over all chromosomes
 cat results/MILeur_23_*.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/MILeur_23.results
  cat results/CIVeur_18_*.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/CIVeur_18.results
 
 cat results/FEMeur_47_*.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/FEMeur_47.results

 cat results/FEM2eur_47_*.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/FEM2eur_47.results
 cat results/FEM3eur_47_*.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/FEM3eur_47.results
  
 cat results/MALeur_39_*.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/MALeur_39.results
 cat results/NOASeur_39_*.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/NOASeur_39.results
 cat results/NOVSeur_39_*.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/NOVSeur_39.results

 for i in {1..42}
 do
  cat results/LOOeur_"$i"_*_1.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/LOOeur_"$i"_.results
 done

 cat results/eur_39_*.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/eur_39.results
 cat results/all_39_*.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/all_39.results
 cat results/aam_39_*.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/aam_39.results
 cat results/lat_38_*.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/lat_38.results

sort -g -k 6 results_cat/eur_39.results | head -n 30 > results_cat/eur_39.results.top30
sort -g -k 6 results_cat/aam_39.results | head -n 30 > results_cat/aam_39.results.top30
sort -g -k 6 results_cat/lat_38.results | head -n 30 > results_cat/lat_38.results.top30
sort -g -k 6 results_cat/all_39.results | head -n 30 > results_cat/all_39.results.top30




 #Manhattan plot of results (list of things to plot is the input)


 ls results_cat | grep -E eur | grep SOBP > metal_outputsA.txt
 ls results_cat | grep -E aam | grep  SOBP > metal_outputsB.txt
 ls results_cat | grep -E all | grep SOBP > metal_outputsC.txt
 ls results_cat | grep -E lat | grep SOBP > metal_outputsD.txt


 njobs=$(wc -l metal_outputsC.txt | awk '{print $1}')
 njobs=$(wc -l metal_outputsD.txt | awk '{print $1}')

 qsub -t1 -l walltime=00:20:00 mh_plot.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m mh_plot_pgc.R -s ManhattanPlotterFunction_colorfixed_max10ylim2_pgc.R -l phase3.locations2 -d metal_outputsA.txt -c green -p 0.05"
  njobs=$(wc -l metal_outputsA.txt | awk '{print $1}')
 qsub -t1 -l walltime=00:20:00 mh_plot.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m mh_plot_pgc.R -s ManhattanPlotterFunction_colorfixed_max10ylim2_pgc.R -l phase3.locations2 -d metal_outputsB.txt -c red -p 0.05"
  njobs=$(wc -l metal_outputsB.txt | awk '{print $1}')


 
 qsub -t1 -l walltime=00:20:00 mh_plot.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m mh_plot_pgc.R -s ManhattanPlotterFunction_colorfixed_max10ylim2_pgc.R -l phase3.locations2 -d metal_outputsC.txt -c blue -p 0.05"
 qsub -t1 -l walltime=00:20:00 mh_plot.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m mh_plot_pgc.R -s ManhattanPlotterFunction_colorfixed_max10ylim2_pgc.R -l phase3.locations2 -d metal_outputsD.txt -c purple -p 0.05"


ncommands=$(wc -l metal_outputsD.txt | awk '{print $1}')
nodesize=16
nodeuse=$(($nodesize ))
totjobs=$(( ($ncommands + $nodeuse - 1 ) / $nodeuse ))

 qsub -t1-$totjobs -l walltime=00:20:00 mh_plot_rapid.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m mh_plot_pgc.R -s ManhattanPlotterFunction_colorfixed_max10ylim2_pgc.R -l phase3.locations2 -d metal_outputsA.txt -c green -p 0.05"
 
 qsub -t1-$totjobs -l walltime=00:20:00 mh_plot_rapid.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m mh_plot_pgc.R -s ManhattanPlotterFunction_colorfixed_max10ylim2_pgc.R -l phase3.locations2 -d metal_outputsB.txt -c red -p 0.05"
 qsub -t1-$totjobs -l walltime=00:20:00 mh_plot_rapid.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m mh_plot_pgc.R -s ManhattanPlotterFunction_colorfixed_max10ylim2_pgc.R -l phase3.locations2 -d metal_outputsC.txt -c blue -p 0.05"
 qsub -t1-$totjobs -l walltime=00:20:00 mh_plot_rapid.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m mh_plot_pgc.R -s ManhattanPlotterFunction_colorfixed_max10ylim2_pgc.R -l phase3.locations2 -d metal_outputsD.txt -c purple -p 0.05"
 

 echo SORTED_PTSD_EA9_ALL_study_specific_PCs1.txt > metal_outputsL.txt


 qsub -t1 -l walltime=00:20:00 mh_plot.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m mh_plot_pgc.R -s ManhattanPlotterFunction_colorfixed_max10ylim2_pgc.R -l phase3.locations2 -d metal_outputsL.txt -c green -p 0.05"




ls results_cat | grep -E lat |  grep 38 > metal_outputsD.txt
 qsub -t1 -l walltime=00:20:00 make_locuszoom.qsub -d $workingdir -e errandout/ -o errandout/ -F "-l phase3.locations2 -d metal_outputsD.txt"


#Produce QQ plots for results
 awk '{print $6}'  results_cat/eur_39.results | grep -v NA > results_cat/eur_39.results.p
 awk '{print $6}'  results_cat/all_39.results | grep -v NA > results_cat/all_39.results.p
 awk '{print $6}'  results_cat/aam_39.results | grep -v NA > results_cat/aam_39.results.p
 awk '{print $6}'  results_cat/lat_38.results | grep -v NA > results_cat/lat_38.results.p

 ls results_cat | grep .p$ > qq_plot.files
 qsub -l walltime=00:20:00 qq_plot.qsub -d $workingdir -e errandout/ -o errandout/ -F "-q qq_plot.R -e 1 -p qq_plot.files "




#Make a forest plot
forestsnp=rs57753395
forestsnpchr=13

#Input file is a .csv of rows with: chr,SNP,ancestry


#Get genotype data for forest plot SNPS
for fp in $(tail -n1 fucker.csv)
do
 forestsnp=$(echo $fp | awk 'BEGIN {FS=","}{print $2}')
 forestsnpchr=$(echo $fp | awk 'BEGIN {FS=","}{print $1}') 
 
 echo $forestsnp $forestsnpchr
 grep -w -m1 $forestsnp metal_inputs/*_"$forestsnpchr" | sed 's/:/ /g' | sed 's/metal_inputs\///g' | sed 's/\.gz_[0-9]*/.gz /g' > results_cat/"$forestsnp".use
done

#Plot forest plots
for fp in $(cat fucker.csv | awk 'NR==5{print}')
do
 forestsnp=$(echo $fp | awk 'BEGIN {FS=","}{print $2}')
 ancgroup=$(echo $fp | awk 'BEGIN {FS=","}{print $3}') 
 echo $forestsnp $ancgroup
 Rscript forest_plot.R results_cat/"$forestsnp".use pgc_ptsd_study_order_v4_descriptors.csv $ancgroup "$forestsnp".plot
done



#Look up genotyping quality of top hits
13 rs57753395 eur  # eur tophit 2. maf 5%, near LOC105370330. Can't find proxy directly (phase 3). The second hit under this was rs73553233. Proxy is rs11618293, Is genotyped on the psych array, most illumina, but not psycchip. Results for this SNP are p ~ 0.0002,error bars are looser for it, effects similar
6 rs34582172 eur # eur tophit 1. maf 40%, is intergenic. Proxy is likely rs4612196, which is still GWS. 
9 rs9410042  lat # lat tophit 1: overall maf 13% (5% in ea and lat, 40% in afr) , on ehmt 1. proxy is rs3125782, its near GWS
18 rs3862732 lat # lat tophit 2: maf is 25% in ea/lat. is likely genotyped
13 rs115539978 aam # aam tophit. intergenic, closest is mir 5007 . rs9536997 is i think I supporting it, however this marker is af~40%,  yet hte top hit is af~=2%. rs12862673 and rs9569233 rs4383040 are near by genotyped. theyre p values are totally insig

ls starting_data > all_input_data.txt
for i in $(cat all_input_data.txt)
do
 echo "Checking file $i"
 zgrep -w -m1 rs4612196 starting_data/"$i" >> check/rs4612196.assoc.dosage
done
#Conclusion: genotyped proxies will differ by study. I think I should use the info metric on the SNP itself.

#look up genotyped proxy marker and plot that

#ls starting_data | grep -v qimr*  | grep -v vetsa* | grep -v -w old | grep -v Metadata_ > studylist.txt
#njobs=$(wc -l studylist.txt | awk '{print $1}')
#qsub -t1-$njobs -l walltime=00:05:00 split_files_infoscores.qsub -d $workingdir -e errandout/ -o errandout/ -F "-f studylist.txt -s starting_data -o snp_info -m 0.01 -i 0.6"
#I don't have necessary columns for vetsa and qimr, leaving out for now...

echo 6,rs4612196,eur > rs4612196.csv
for fp in $(cat rs4612196.csv )
do
 forestsnp=$(echo $fp | awk 'BEGIN {FS=","}{print $2}')
 forestsnpchr=$(echo $fp | awk 'BEGIN {FS=","}{print $1}') 
 
 echo $forestsnp $forestsnpchr
 zgrep -w -m1 $forestsnp snp_info/*_"$forestsnpchr".gz | sed 's/:/ /g' | sed 's/snp_info\///g' | sed 's/\.gz_[0-9]*//g' > results_cat/"$forestsnp".info
done


zgrep -w -m1 rs11618293 snp_info/*_"$forestsnpchr".gz | sed 's/:/ /g' | sed 's/snp_info\///g' | sed 's/\.gz_[0-9]*//g' > results_cat/rs11618293.info



#Add civilian, military, etc. subgrouping


#Make an LDSC file
grep -v ??????? results_cat/eur_39.results  | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"12813","35640"}' | zip  > eur_39_moststudies.results.zip

grep -v ??? results_cat/MILeur_23.results  | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"6086","21545"}' | zip  > MILeur_23_moststudies.results.zip
grep -v ??? results_cat/CIVeur_18.results  | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"6208","12612"}' | zip  > CIVeur_18_moststudies.results.zip


munge_sumstats.py --sumstats MILeur_23_moststudies.results --out MILeur_23_moststudies_munge
munge_sumstats.py --sumstats CIVeur_18_moststudies.results --out CIVeur_23_moststudies_munge

ldsc.py \
--rg MILeur_23_moststudies_munge.sumstats.gz,CIVeur_23_moststudies_munge.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out mil_civ
less mil_civ.log


grep -v ? results_cat/eur_39.results  | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"12813","35640"}' | zip  > eur_39_allstudies.results.zip

grep -v ? results_cat/eur_39.results  | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"12813","35640"}' | zip  > eur_39_allstudies.results.zip
grep -v ?? results_cat/FEMeur_47.results | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"5914","8759"}' | zip  > FEMeur_47_allstudies.results.zip
grep -v ?? results_cat/MALeur_39.results | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"6381","23580"}' | zip  > MALeur_39_allstudies.results.zip


grep -v ?? results_cat/FEM2eur_47.results | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"5315","8008"}'  > FEM2eur_47_allstudies.results

grep -v ?? results_cat/FEM3eur_47.results | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"3488","7011"}' | zip > FEM3eur_47_allstudies.results.zip


awk 'BEGIN {FS="?"; var=0 } {if (NF>0){ var=var + (NF-1) } } END{print var}' testfile1.txt


grep -v ?? results_cat/NOASeur_39.results | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"9963","26810"}'   > NOASeur_39_allstudies.results.zip

grep -v ?? results_cat/NOVSeur_39.results | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"12000","34000"}'   > NOVSeur_39_allstudies.results.zip

results_cat/NOASeur_39.results
 
unzip MALeur_39_allstudies.results.zip
unzip FEMeur_47_allstudies.results.zip

#All studies

awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"12813","35640"}'  results_cat/eur_39.results  > eur_39_allstudies.results
awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"5914","8759"}' results_cat/FEMeur_47.results  > FEMeur_47_allstudies.results
awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"6381","23580"}' results_cat/MALeur_39.results > MALeur_39_allstudies.results

munge_sumstats.py --sumstats eur_39_allstudies.results --out results_cat/eur_39_allstudies_munge
munge_sumstats.py --sumstats FEMeur_47_allstudies.results --out results_cat/FEMeur_47_allstudies_munge
munge_sumstats.py --sumstats MALeur_39_allstudies.results --out results_cat/MALeur_39_allstudies_munge


cd results_cat/

ldsc.py \
--h2 eur_39_allstudies_munge.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--samp-prev 0.264 \
--pop-prev 0.09 \
--out eur38_h2

ldsc.py \
--h2 FEMeur_47_allstudies_munge.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--samp-prev 0.40 \
--pop-prev 0.11 \
--out femeur47_h2

ldsc.py \
--h2 MALeur_39_allstudies_munge.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--samp-prev 0.21 \
--pop-prev 0.05 \
--out maleur38_h2

#Most studies ( to match LOO


grep -v ??  results_cat/eur_39.results  | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"12813","35640"}'  > eur_39_moststudies.results
grep -v ??  results_cat/FEMeur_47.results | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"5914","8759"}'   > FEMeur_47_moststudies.results
grep -v ?? results_cat/MALeur_39.results | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"6381","23580"}'   > MALeur_39_moststudies.results

munge_sumstats.py --sumstats eur_39_moststudies.results --out results_cat/eur_39_moststudies_munge
munge_sumstats.py --sumstats FEMeur_47_moststudies.results --out results_cat/FEMeur_47_moststudies_munge
munge_sumstats.py --sumstats MALeur_39_moststudies.results --out results_cat/MALeur_39_moststudies_munge


cd results_cat/

ldsc.py \
--h2 eur_39_moststudies_munge.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--samp-prev 0.264 \
--pop-prev 0.09 \
--out eur39_h2_most

ldsc.py \
--h2 FEMeur_47_moststudies_munge.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--samp-prev 0.40 \
--pop-prev 0.11 \
--out femeur47_h2_most

ldsc.py \
--h2 MALeur_39_moststudies_munge.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--samp-prev 0.21 \
--pop-prev 0.05 \
--out maleur39_h2_most



#h2 is lower if removing small studies

##LOO - fix so that it doesnt do 50 of these
ls results_cat | grep LOOM > ldsc_jobs15.txt #Put in a grep for the last row of the study input data if you ONLY want the total data meta analysis
njobs=$(wc -l ldsc_jobs15.txt | awk '{print $1}')
qsub -t1-$njobs -l walltime=00:05:00 ldsc.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m ldsc_jobs15.txt -p results_cat"


grep -v ??? results_cat/eur_39.results | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6, "1556","6618" }'  | zip  > eur_9_allstudies.results.zip

grep -v ?.*?.* results_cat/all_39.results | awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6, "19000","50000" }'  | zip  > all_39_allstudies.results.zip

 for i in {1..33}
 do
  cat results/LOOMeur_"$i"_*_1.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/LOOMeur_"$i"_.results

  #rm -f LOOeur_"$i"_allstudies.results
  #rm -f results_cat/LOOeur_"$i"_.results
 done


grep #Look up genotyping quality of top hits
grep -E "rs57753395|rs34582172|rs9410042|rs3862732|rs115539978" /home/cnieverg/report_data/testrun/results_cat/eur_39.results
grep -E "rs57753395|rs34582172|rs9410042|rs3862732|rs115539978" /home/cnieverg/report_data/testrun/results_cat/aam_39.results
grep -E "rs57753395|rs34582172|rs9410042|rs3862732|rs115539978" /home/cnieverg/report_data/testrun/results_cat/lat_38.results
 
#Add positions to data
 ls results_cat | grep -E eur_39\|lat_38\|aam_39 | grep -v NOVS | grep -v MAL | grep -v NOVS |  grep .results$  | grep lat > metal_outputsA.txt

 qsub -t1 -l walltime=00:20:00 add_positions.qsub -d $workingdir -e errandout/ -o errandout/ -F "-l phase3.locations2 -d metal_outputsA.txt "



bash add_positions.qsub -l phase3.locations2 -d results_cat/lat_38.results

#Filter to select snps and do ldsc

 maf=0.01
 info=0.9
 zcat  /home/cnieverg/report_data/testrun/starting_data/daner_mrsc_eur_analysis_run3_males.gz | awk  -v maf=$maf -v info=$info '{if  (NR == 1) print "SNP";  if (($8 > info) &&  ($7 > maf) && ($7 < 1-maf)  ) print $2}' | LC_ALL=C sort -k1b,1 > good_snps.snplist

awk '{if(NR==1) print "SNP","A1","A2","BETA","SE","P", "N_case","N_Control"; else print $1,$2,$3,$4,$5,$6,"12813","35640"}'  results_cat/eur_39.results | LC_ALL=C sort -k1b,1 > eur_39.results.munge.sorted
LC_ALL=C join good_snps.snplist eur_39.results.munge.sorted | sort -g -k6 > eur_39.results.munge.sorted.fin
munge_sumstats.py --sumstats eur_39.results.munge.sorted.fin --out results_cat/eur_39_allstudies_goodsnpsonly_munge
cd results_cat/

ldsc.py \
--h2 eur_39_allstudies_goodsnpsonly_munge.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--samp-prev 0.264 \
--pop-prev 0.09 \
--out eur38_goodsnps_h2

#Make h2 results into a LZ plot
grep  Liability results_cat/LOOeur_*_.results_h2.log | awk 'BEGIN{FS=":"; OFS=" "}{print $1,$3}' | sed 's/results_cat\/LOOeur_//g' | sed 's/_.results_h2.log//g' | sed 's/(//g' | sed 's/)//g' |  > results_cat/looeur_h2.ldsc_res
grep  Liability results_cat/LOOFeur_*_.results_h2.log | awk 'BEGIN{FS=":"; OFS=" "}{print $1,$3}' | sed 's/results_cat\/LOOFeur_//g' | sed 's/_.results_h2.log//g' | sed 's/(//g' | sed 's/)//g' > results_cat/loofeur_h2.ldsc_res
grep  Liability results_cat/LOOMeur_*_.results_h2.log | awk 'BEGIN{FS=":"; OFS=" "}{print $1,$3}' | sed 's/results_cat\/LOOMeur_//g' | sed 's/_.results_h2.log//g' | sed 's/(//g' | sed 's/)//g' > results_cat/loomeur_h2.ldsc_res


R
 library(metafor)
 study_abbrs <- read.csv('pgc_ptsd_study_order_v6.csv',stringsAsFactors=F)
 
#All studies
 liability <- read.table('results_cat/looeur_h2.ldsc_res',header=F)
 names(liability) <- c("study_ord","h2","se")

 studies <- read.table('ea_studies_list.txt')
 studies$study <- c(1:dim(studies)[1])
 names(studies) <- c("study","study_ord")
 studies$study <- gsub("_22","",studies$study)


 named_studies <- merge(studies,study_abbrs,by="study",all.x=T)
 dat <- merge(liability,named_studies,by="study_ord",all.x=T)
 dat <- dat[order(dat$cases,decreasing=T),]
 meta_res <- rma(yi=dat$h2,sei=dat$se,method="FE",slab=dat$abbr)

 pdf("plots/LOO_forest.pdf",7,9)
 forest(meta_res)
 dev.off()


