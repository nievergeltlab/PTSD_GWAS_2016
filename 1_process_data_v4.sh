#Set working dir
 workingdir=/home/cnieverg/report_data/testrun

#To the user: Install metal, write the path to the binary here
 metalpath=/home/cnieverg/meta_istss/generic-metal/metal

 cd $workingdir

#Make a folder for all studies to be included.
 mkdir starting_data

#Make a folder of where filtered data will go
 mkdir metal_inputs
 mkdir metal_inputs/males
 mkdir metal_inputs/females

#Make a folder of where the metal scripts will go
 mkdir metal_scripts
 mkdir metal_scripts/males
 mkdir metal_scripts/females

#Results location
 mkdir metal_results
 mkdir metal_results/males
 mkdir metal_results/females
 
 
#Make directory for metasoft file input, scripts for running
 mkdir metasoft_scripts
 mkdir metasoft_scripts/males
 mkdir metasoft_scripts/females
 
 
 mkdir metasoft_results
 mkdir metasoft_results/males
 mkdir metasoft_results/females
 
#RE meta analysis file locations
 mkdir metasoft_inputs
 mkdir metasoft_inputs/males
 mkdir metasoft_inputs/females

 
#PLot and output locations 
 mkdir plots
 mkdir temporary_files
 mkdir errandout


#To the User: Move all Ricopili results files into folder starting_data

#Unzip the summary datasets, e.g.
#This is for the combined gender sets
 ls *.tar > tarfiles.txt

 for file in $(cat tarfiles.txt)
 do
  tar xvf "$file" --exclude '*males*' --exclude '*females*' --wildcards 
 done

# #This is for the gender stratified sets
# ls *.tar > tarfiles.txt

# for file in $(cat tarfiles.txt)
# do
 # tar xvf "$file" --strip-components=4   --wildcards '*males*' '*females*'
 # rm $file
# done

#Using a datasets where all SNP positions are given, list all SNP positions into a file
#To the user: give a file name. Here I used grady positions
 zcat starting_data/daner_gtpc_aam_analysis_run3.gz  | awk '{print  $2,$1,$3}' | LC_ALL=C sort -k1,1b > phase3.locations
 LC_ALL=C sort -k1,1b phase3.locations > phase3.locations2

#These SNP positions will be saved as an R file also
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
#THis will be levle 1 filtering just to make the data small, then will put on higher filtering later
#NA Fitler: Filter out markers without results
 maf=0.005
 info=0.6

#To the User: Make a list of all studies of unrelateds called studylist.txt
ls starting_data | grep PTSD13_All_MAF1_INFO4_ricopili.txt.affix.logscale.gz > studylist.txt
echo daner_coga_eur_analysis_run3.meta.gz > studylist.txt

#To the User: Make a list of all studies of relateds called studylist_related.txt

ls starting_data  | grep PTSD13_All_MAF1_INFO4_ricopili.txt.affix.logscale.gz > studylist.txt
ls starting_data  | grep -P _male\|_Male > male_studylist.txt
ls starting_data  | grep -P _female\|_Female > female_studylist.txt

ls starting_data  | grep -P _male\|_Male | grep TRAC > male_studylist.txt
ls starting_data  | grep -P _female\|_Female | grep TRAC > female_studylist.txt

#Code for checking ns..

for femaleset in $(ls | grep 22)
do
 awk -v fsaz=$femaleset 'NR==2{print fsaz,$6,$7}' $femaleset >> ~/report_data/testrun/female_ns.txt
done

#If you missed a study, create a list like this:
#echo daner_comc_aam_analysis_run3.gz  daner_meg2_af3_analysis_run5.gz | sed 's/ /\n/g' > studylist.txt
 njobs=$(wc -l studylist.txt | awk '{print $1}')
 qsub -t1-$njobs -l walltime=00:05:00 split_files.qsub -d $workingdir -e errandout/ -o errandout/ -F "-f studylist.txt -s $workingdir -o metal_inputs -m $maf -i $info"

 njobs=$(wc -l male_studylist.txt | awk '{print $1}')
 qsub -t1-$njobs -l walltime=00:05:00 split_files.qsub -d $workingdir -e errandout/ -o errandout/ -F "-f male_studylist.txt -s $workingdir -o metal_inputs/males -m $maf -i $info"

 njobs=$(wc -l female_studylist.txt | awk '{print $1}')
 qsub -t1-$njobs -l walltime=00:05:00 split_files.qsub -d $workingdir -e errandout/ -o errandout/ -F "-f female_studylist.txt -s $workingdir -o metal_inputs/females -m $maf -i $info"
 

#UKBB had to be done separately, it was a format into itself, the gemma ones are now edited using the gemma pipeline (note to self, archive what I did with the qimr data)

###Fixed effects meta-analysis

#Make metal input files
#Requires as an input an ordered list of studies called e.g. study_order.csv. Read the .R file for a description of how the input should look
#Note: SOme settings for mETAL are HARD SET in metal_header.mi and metal_footer.mi. EDIT THESE FILES to your preferences for metal!!

 Rscript make_mi_files_v2.r study_order.csv metal_scripts/ metal_inputs/ metal_results/ mf
 Rscript make_mi_files_v2.r study_order_males.csv metal_scripts/males  metal_inputs/males metal_results/males males
 Rscript make_mi_files_v2.r study_order_females.csv metal_scripts/females  metal_inputs/females metal_results/females females

#Without UKBB (only the overall and european need to be redone!!!)
 Rscript make_mi_files_v2.r noukbb_study_order.csv metal_scripts_noukbb/ metal_inputs/ metal_results_noukbb/ mf
 Rscript make_mi_files_v2.r study_order_males_noukbb.csv metal_scripts_noukbb/ metal_inputs/males metal_results_noukbb/males males
 Rscript make_mi_files_v2.r study_order_females_noukbb.csv metal_scripts_noukbb/ metal_inputs/females metal_results_noukbb/females females

#Cross validated run
 mkdir metal_results_cv/
 mkdir metal_scripts_cv/
 
for groupnum in {1..10}
do
 Rscript make_mi_files_v2_cv.r crossval/groupa_"$groupnum".csv metal_scripts_cv/ metal_inputs/ metal_results_cv/ mf_cva_"$groupnum"
 Rscript make_mi_files_v2_cv.r crossval/groupb_"$groupnum".csv metal_scripts_cv/ metal_inputs/ metal_results_cv/ mf_cvb_"$groupnum"
done
 
 
###Run metal

 ls metal_scripts/*  > metal_jobs1.txt # | grep -P "eur|all"  > metal_jobs1.txt #Put in a grep for the last row of the study input data if you ONLY want the total data meta analysis
 ls metal_scripts/females/*  > metal_jobs2.txt # | grep -P "eur|all"  > metal_jobs1.txt #Put in a grep for the last row of the study input data if you ONLY want the total data meta analysis
 ls metal_scripts/males/*   > metal_jobs3.txt # | grep -P "eur|all"  > metal_jobs1.txt #Put in a grep for the last row of the study input data if you ONLY want the total data meta analysis

 ls metal_scripts_noukbb/* | grep males > metal_jobs4.txt 
 ls metal_scripts_noukbb/* | grep females > metal_jobs5.txt 

 cat metal_jobs4.txt  metal_jobs5.txt  > metal_jobs.txt
 
 ls metal_scripts_summary/* | grep eur  > metal_jobs6.txt 
 ls metal_scripts_genotype/* | grep eur >> metal_jobs6.txt 
 
 ls metal_scripts_smallstudies/* | grep eur > metal_jobs.txt 
 
  ls metal_scripts_cv/* | grep eur > metal_jobs.txt 
 
  cat metal_jobs6.txt  > metal_jobs.txt
 
 
 ncommands=$(wc -l metal_jobs.txt | awk '{print $1}')
 nodesize=16
 nodeuse=$(($nodesize ))
 totjobs=$(( ($ncommands + $nodeuse - 1 ) / $nodeuse ))

 qsub -t1-$totjobs -l walltime=00:10:00 run_meta_v2_loo.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m metal_jobs.txt -p $metalpath -n $nodeuse"

 
#Concatenate over all chromosomes (sample file shown)

#all data
 cat metal_results/eur_*_1.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/eur_dec28_2017_maf01_info6.results
 cat metal_results/aam_*_1.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/aam_dec28_2017_maf01_info6.results
 cat metal_results/lat_*_1.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/lat_dec28_2017_maf01_info6.results
 cat metal_results/all_*_1.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/all_dec28_2017_maf01_info6.results

 #Males
 cat metal_results/males/eur_*_1.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/eur_dec28_2017_maf01_info6_males.results
 cat metal_results/males/aam_*_1.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/aam_dec28_2017_maf01_info6_males.results
 cat metal_results/males/lat_*_1.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/lat_dec28_2017_maf01_info6_males.results
 cat metal_results/males/all_*_1.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/all_dec28_2017_maf01_info6_males.results

#Females
 cat metal_results/females/eur_*_1.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/eur_dec28_2017_maf01_info6_females.results
 cat metal_results/females/aam_*_1.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/aam_dec28_2017_maf01_info6_females.results
 cat metal_results/females/lat_*_1.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/lat_dec28_2017_maf01_info6_females.results
 cat metal_results/females/all_*_1.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/all_dec28_2017_maf01_info6_females.results
 
 #all data, NO UKBB
 cat metal_results_noukbb/eur_*_1.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/eur_dec28_2017_maf01_info6_noukbb.results
 cat metal_results_noukbb/all_*_1.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/all_dec28_2017_maf01_info6_noukbb.results

 cat metal_results_noukbb/males/eur_*_1.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' | gzip > results_cat/eur_dec28_2017_maf01_info6_males_noukbb.results.gz
 cat metal_results_noukbb/females/eur_*_1.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' | gzip > results_cat/eur_dec28_2017_maf01_info6_females_noukbb.results.gz

#Freeze 1
 cat metal_results_freeze1/eur_*_1.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/eur_dec28_2017_maf01_info6_freeze1.results

#Freeze 2
 cat metal_results_freeze2/eur_*_1.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/eur_dec28_2017_maf01_info6_freeze2.results

#Freeze 2 (no ukbb)
 cat metal_results_freeze2noukbb/eur_*_1.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/eur_dec28_2017_maf01_info6_freeze2noukbb.results

#Summary only v genotype studies
 cat metal_results_genotype/eur_*_1.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/eur_dec28_2017_maf01_info6_genotype.results
 cat metal_results_summary/eur_*_1.tbl  | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/eur_dec28_2017_maf01_info6_summary.results

#Smallstudies
 cat metal_results_smallstudies/eur_*_1.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/eur_dec28_2017_maf01_info6_smallstudies.results

#Cross validated
for group in a b
do
 for num in {1..10}
 do
  cat metal_results_cv/eur_mf_cv"$group"_"$num"_*_1.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/eur_dec28_2017_maf01_info6_cv"$group"_"$num".results
 done
done


##for LDSC/munging of Europeans, Calculate info score per study for all markers, take weighted avg info
# for files in $(ls metal_inputs)
 #do
#  awk '{if (NR==1) print "SNP","WInfo","Neff"print $1,$5*$8
 
#Filter to only markers with effective N > 25% of neff total
 awk '{qcol=$11; gsub("?","",qcol); if ($14 >= .25*86058 && length(qcol) >= 3 ) print }' results_cat/all_dec28_2017_maf01_info6.results | sort -g -k10 > results_cat/all_dec28_2017_maf01_info6.results_neff
 awk '{qcol=$11; gsub("?","",qcol); if ($14 >= .25*70237.5 && length(qcol) >= 3) print}' results_cat/eur_dec28_2017_maf01_info6.results | sort -g -k10 > results_cat/eur_dec28_2017_maf01_info6.results_neff
 awk '{qcol=$11; gsub("?","",qcol); if ($14 >= .25*11321.5 && length(qcol) >= 3) print}' results_cat/aam_dec28_2017_maf01_info6.results | sort -g -k10 > results_cat/aam_dec28_2017_maf01_info6.results_neff
 awk '{qcol=$11; gsub("?","",qcol); if ($14 >= .25*4498.96 && length(qcol) >= 3) print}' results_cat/lat_dec28_2017_maf01_info6.results | sort -g -k10 > results_cat/lat_dec28_2017_maf01_info6.results_neff
 
 awk '{qcol=$11; gsub("?","",qcol); if ($14 >= .25*36526.9 && length(qcol) >= 3) print}' results_cat/all_dec28_2017_maf01_info6_males.results | sort -g -k10 > results_cat/all_dec28_2017_maf01_info6_males.results_neff
 awk '{qcol=$11; gsub("?","",qcol); if ($14 >= .25*30565.9 && length(qcol) >= 3) print}' results_cat/eur_dec28_2017_maf01_info6_males.results | sort -g -k10 > results_cat/eur_dec28_2017_maf01_info6_males.results_neff
 awk '{qcol=$11; gsub("?","",qcol); if ($14 >= .25*4701.54 && length(qcol) >= 3) print}' results_cat/aam_dec28_2017_maf01_info6_males.results | sort -g -k10 > results_cat/aam_dec28_2017_maf01_info6_males.results_neff
 awk '{qcol=$11; gsub("?","",qcol); if ($14 >= .25*1259.45 && length(qcol) >= 3) print}' results_cat/lat_dec28_2017_maf01_info6_males.results | sort -g -k10 > results_cat/lat_dec28_2017_maf01_info6_males.results_neff
    
 awk '{qcol=$11; gsub("?","",qcol); if ($14 >= .25*47572.6 && length(qcol) >= 3) print}' results_cat/all_dec28_2017_maf01_info6_females.results | sort -g -k10 > results_cat/all_dec28_2017_maf01_info6_females.results_neff
 awk '{qcol=$11; gsub("?","",qcol); if ($14 >= .25*38426.1 && length(qcol) >= 3) print}' results_cat/eur_dec28_2017_maf01_info6_females.results | sort -g -k10 > results_cat/eur_dec28_2017_maf01_info6_females.results_neff
 awk '{qcol=$11; gsub("?","",qcol); if ($14 >= .25*6063.73 && length(qcol) >= 3) print}' results_cat/aam_dec28_2017_maf01_info6_females.results | sort -g -k10 > results_cat/aam_dec28_2017_maf01_info6_females.results_neff
 awk '{qcol=$11; gsub("?","",qcol); if ($14 >= .25*3082.71 && length(qcol) >= 1) print}' results_cat/lat_dec28_2017_maf01_info6_females.results | sort -g -k10 > results_cat/lat_dec28_2017_maf01_info6_females.results_neff # have to set to 1 for this group, N studies too small!
     
 #No ukBB
 awk '{qcol=$11; gsub("?","",qcol); if ($14 >= .25*32197.5 && length(qcol) >= 3) print}' results_cat/eur_dec28_2017_maf01_info6_noukbb.results | sort -g -k10 > results_cat/eur_dec28_2017_maf01_info6_noukbb.results_neff
 zcat results_cat/eur_dec28_2017_maf01_info6_males_noukbb.results.gz | awk '{qcol=$11; gsub("?","",qcol); if ($14 >= .25*17328.3 && length(qcol) >= 3) print}'  | sort -g -k10 | gzip > results_cat/eur_dec28_2017_maf01_info6_males_noukbb.results_neff.gz
 zcat results_cat/eur_dec28_2017_maf01_info6_females_noukbb.results.gz  | awk '{qcol=$11; gsub("?","",qcol); if ($14 >= .25*13752.9 && length(qcol) >= 3) print}' | sort -g -k10  | gzip  > results_cat/eur_dec28_2017_maf01_info6_females_noukbb.results_neff.gz

#Summary and not summary
 awk '{qcol=$11; gsub("?","",qcol); if ($14 >= .25*45067.4 && length(qcol) >= 3) print}' results_cat/eur_dec28_2017_maf01_info6_summary.results | sort -g -k10 | gzip > results_cat/eur_dec28_2017_maf01_info6_summary.results_neff.gz
 awk '{qcol=$11; gsub("?","",qcol); if ($14 >= .25*25170.1 && length(qcol) >= 3) print}' results_cat/eur_dec28_2017_maf01_info6_genotype.results | sort -g -k10 | gzip > results_cat/eur_dec28_2017_maf01_info6_genotype.results_neff.gz
  

 
 #Cross validated
for group in a b
do
 for num in {1..10}
 do
  samplesize=$(head -n2 results_cat/eur_dec28_2017_maf01_info6_cv"$group"_"$num".results | tail -n1 | awk '{print $14}' ) 
  echo $samplesize
  awk -v samples=$samplesize '{qcol=$11; gsub("?","",qcol); if ($14 >= .25*samples && length(qcol) >= 3) print}' results_cat/eur_dec28_2017_maf01_info6_cv"$group"_"$num".results | gzip > results_cat/eur_dec28_2017_maf01_info6_cv"$group"_"$num".results_neff.gz
 done
done

#Add positional info (for what??)
# LC_ALL=C join <(awk '{if(NR==1) print "MarkerName","CHR","BP"; else print $1,$2,$3}' phase3.locations2)  <(LC_ALL=C sort -k1b,1 results_cat/eur_dec18_2017_maf01_info6.results_neff) | sort -g -k 12 > results_cat/eur_dec18_2017_maf01_info6.results_neff_pos
# LC_ALL=C join <(awk '{if(NR==1) print "MarkerName","CHR","BP"; else print $1,$2,$3}' phase3.locations2)  <(LC_ALL=C sort -k1b,1 results_cat/all_dec18_2017_maf01_info6.results_neff) | sort -g -k 12 > results_cat/all_dec18_2017_maf01_info6.results_neff_pos
# LC_ALL=C join <(awk '{if(NR==1) print "MarkerName","CHR","BP"; else print $1,$2,$3}' phase3.locations2)  <(LC_ALL=C sort -k1b,1 results_cat/aam_dec18_2017_maf01_info6.results_neff) | sort -g -k 12 > results_cat/aam_dec18_2017_maf01_info6.results_neff_pos
# LC_ALL=C join <(awk '{if(NR==1) print "MarkerName","CHR","BP"; else print $1,$2,$3}' phase3.locations2)  <(LC_ALL=C sort -k1b,1 results_cat/lat_dec18_2017_maf01_info6.results_neff) | sort -g -k 12 > results_cat/lat_dec18_2017_maf01_info6.results_neff_pos
 
###Manhattan plot of results (list of things to plot is the input)
 ls results_cat | grep eur | grep neff$ > metal_outputsA.txt # echo eur_dec18_2017_maf01_info6.results_neff > metal_outputsA.txt
 ls results_cat | grep aam | grep neff$ > metal_outputsB.txt # echo aam_dec18_2017_maf01_info6.results_neff > metal_outputsB.txt
 ls results_cat | grep lat | grep neff$ > metal_outputsD.txt  # echo lat_dec18_2017_maf01_info6.results_neff > metal_outputsD.txt
 ls results_cat | grep all | grep neff$ > metal_outputsC.txt # echo all_dec18_2017_maf01_info6.results_neff > metal_outputsC.txt

 ncommands=$(wc -l metal_outputsA.txt | awk '{print $1}')
 nodesize=16
 nodeuse=$(($nodesize ))
 totjobs=$(( ($ncommands + $nodeuse - 1 ) / $nodeuse ))

 qsub -t1-$totjobs -l walltime=00:20:00 mh_plot_rapid_v2.qsub -d $workingdir -e errandout/ -o errandout/ -V -F "-m mh_plot_pgc.R -s ManhattanPlotterFunction_colorfixed_max10ylim2_pgc.R -l phase3.locations2 -d metal_outputsA.txt -c green -p 0.05 -n 10"
 qsub -t1-$totjobs -l walltime=00:20:00 mh_plot_rapid_v2.qsub -d $workingdir -e errandout/ -o errandout/ -V -F "-m mh_plot_pgc.R -s ManhattanPlotterFunction_colorfixed_max10ylim2_pgc.R -l phase3.locations2 -d metal_outputsB.txt -c red -p 0.05 -n 10"
 qsub -t1-$totjobs -l walltime=00:20:00 mh_plot_rapid_v2.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m mh_plot_pgc.R -s ManhattanPlotterFunction_colorfixed_max10ylim2_pgc.R -l phase3.locations2 -d metal_outputsC.txt -c blue -p 0.05 -n 10"
 qsub -t1-$totjobs -l walltime=00:20:00 mh_plot_rapid_v2.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m mh_plot_pgc.R -s ManhattanPlotterFunction_colorfixed_max10ylim2_pgc.R -l phase3.locations2 -d metal_outputsD.txt -c purple -p 0.05 -n 10" 
 
#QQ plots of results

#Print p-value column of all data
 for files in $(ls results_cat | grep neff )
 do
  awk '{print $10}' results_cat/"$files" | gzip > results_cat/"$files".p.gz
 done

#List all relevant files
 ls results_cat | grep .results_neff.p.gz$ | sed 's/.gz//g' > qq_plot.files

#Submit job
 qsub -l walltime=00:20:00 qq_plot.qsub -d $workingdir -e errandout/ -o errandout/ -F "-q qq_plot.R -e 1 -p qq_plot.files "

#Make locus-zoom ready files - for regional plots of top hits 

 qsub -t1 -l walltime=00:20:00 make_locuszoom.qsub -d $workingdir -e errandout/ -o errandout/ -F "-l phase3.locations2 -d metal_outputsA.txt -n 6"
 qsub -t1 -l walltime=00:20:00 make_locuszoom.qsub -d $workingdir -e errandout/ -o errandout/ -F "-l phase3.locations2 -d metal_outputsB.txt -n 6"
 qsub -t1 -l walltime=00:20:00 make_locuszoom.qsub -d $workingdir -e errandout/ -o errandout/ -F "-l phase3.locations2 -d metal_outputsD.txt -n 6"

 
###Forest plots

#User: Make a csv file with the chromosome, rs-id, and ancestry of each top hit (use 'all' if you want it plotted in all data)

echo "6,rs71563706,eur,mf
6,rs9364611,eur,mf
13,rs115539978,aam,mf
1,rs514370,eur,males
1,rs148757321,eur,males
19,chr19:53988841,eur,males
19,rs8112292,eur,males
9,rs3123508,lat,males
6,rs142174523,aam,males" > hits.csv #Latino hit is just the overall anc tophit, but a few markers away


9,rs9410042,lat,mf
#18,rs3862732,lat,mf" > mfhits.csv
13,rs57753395,eur,mf

echo "11,rs113290828,aam" > femhits.csv #ALmost GWS
 
echo "3,rs9822553,eur,mf
12,chr12:48914460,lat,males

14,chr14:26167750,aam,males
14,rs138943387,aam,males" > falsepos.csv


#Load summary GWAS data for forest plot SNPS

for fp in $(cat hits.csv  )
do
 forestsnpchr=$(echo $fp | awk 'BEGIN {FS=","}{print $1}') 
 forestsnp=$(echo $fp | awk 'BEGIN {FS=","}{print $2}')
 ancgroup=$(echo $fp | awk 'BEGIN {FS=","}{print $3}') 
 gender=$(echo $fp | awk 'BEGIN {FS=","}{print $4}') 

 if [ $gender == "mf" ]
 then
  grep -w -m1 $forestsnp metal_inputs/*_"$forestsnpchr" | sed 's/:/ /' | sed 's/metal_inputs\///g' | sed 's/\.gz_[0-9]*/.gz /g' > results_cat/"$forestsnp".use
  sex=""
  sex2="mf"
 fi
 if [ $gender ==  "males" ]
 then
  grep -w -m1 $forestsnp metal_inputs/males/*_"$forestsnpchr" | sed 's/:/ /' | sed 's/metal_inputs\///g' | sed 's/males\///g' | sed 's/\.gz_[0-9]*/.gz /g' > results_cat/"$forestsnp".use.males
  sex=".males"
  sex2="males"
 fi
 
 if [ $gender ==  "females" ]
 then
  grep -w -m1 $forestsnp metal_inputs/females/*_"$forestsnpchr" | sed 's/:/ /' | sed 's/metal_inputs\///g' | sed 's/females\///g' | sed 's/\.gz_[0-9]*/.gz /g' > results_cat/"$forestsnp".use.females
  sex=".females"
  sex2="females"
 fi

 echo "Extracting data for $forestsnp $forestsnpchr"

 Rscript forest_plot_v2.R results_cat/"$forestsnp".use"$sex" "$sex2" pgc_ptsd_study_order_v8.csv $ancgroup FE "$forestsnp".plot_RE
 java -jar Metasoft.jar -input results_cat/"$forestsnp".plot_RE_"$ancgroup".msin -output results_cat/"$forestsnp".plot_RE_"$ancgroup"_msoft -log results_cat/"$forestsnp".plot_RE.mslog -mvalue -mvalue_method mcmc 
 #If N studies small, mvalue calculation tractable. i.e. 10 or less. Otherwise use -mvalue_method mcmc 
 python pmplot.py results_cat/"$forestsnp".plot_RE_"$ancgroup".msin results_cat/"$forestsnp".plot_RE_"$ancgroup"_msoft  results_cat/"$forestsnp".plot_RE_"$ancgroup".msinstudynames   results_cat/"$forestsnp".plot_RE_"$ancgroup".msinstudyorder "$forestsnp".plot_RE "$forestsnp" plots/"$forestsnp".plot_RE_"$ancgroup".pmplot.pdf
 
done

args=(commandArgs(TRUE))
output_file = args[1]
input_file = args[2]
gene_file = args[3]
studynames_file = args[4]
studyorder_file = args[5]
pmplot_file = args[6]

Rscript forestpmplot.R  results_cat/"$forestsnp".plot_RE_"$ancgroup".msin results_cat/"$forestsnp".plot_RE_"$ancgroup"_msoft  results_cat/"$forestsnp".plot_RE_"$ancgroup".msinstudynames   results_cat/"$forestsnp".plot_RE_"$ancgroup".msinstudyorder "$forestsnp".plot_RE "$forestsnp" plots/"$forestsnp".plot_RE_"$ancgroup".pmplot.pdf

###Random effects meta-analysis (maybe left to its own script..)

#Make input file scripts
 Rscript make_metasoft_files.R study_order.csv mf
 Rscript make_metasoft_files.R study_order_males.csv males
 Rscript make_metasoft_files.R study_order_females.csv females
 
 Rscript make_metasoft_files.R noukbb_study_order.csv mf

#Make sure you check if results are with or without UKBB included..


##Reformat metasoft data:

 #Split the locations data into chunks
  #Needs to be done only one time
 #qsub -t1-22 -d $workingdir  -e errandout/ -o errandout/ make_phase3_locations.pbs -lwalltime=00:15:00 
 
 #Set frequency and info filters (should match FE GWAS)
 freq=0.01
 info=0.6
 
 ancgroup="all"

 #Reformat data for each chromosome (chromosomes are noted by job array thing) 
  for ancgroup in all aam  eur  lat 
 do
 for gender in mf # females #mf males females
 do
  qsub -t1-22 -d $workingdir  -e errandout/ -o errandout/ make_metasoft_bigchr_v2.pbs -lwalltime=1:20:00 -F "-r $ancgroup -f $freq -i $info -g $gender"
 done
done

#Run metasoft analysis (add -b bem for binary effects. take a long time, definitely longer than 20 mins!!)
for gender in mf # females # mf males females
do
 for ancgroup in eur # all # aam # eur  #lat aam #all
 do
  qsub -t1-22 -d $workingdir  -e errandout/ -o errandout/ run_metasoft.pbs -lwalltime=0:30:00  -F "-r $ancgroup -b xxx -g $gender"
 done
done


#Filter results to just the essential columns, concatenate them together 
 #note: this must be reformated for binary effects models
for gender in mf  # females males 
do
 if [ $gender ==  "mf" ]
 then
  sex=""
  sex2=""
 fi
 
 if [ $gender !=  "mf" ]
 then
  sex=$gender
  sex2=_"$gender"
 fi
  
 for ancgroup in eur # all #  aam  all  lat # 
 do
  for i in $(ls metasoft_results | grep "$ancgroup" | grep .msoftout$)
  do 
   
   awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$13,$15}' metasoft_results/"$sex"/"$i" | grep -v NA > temporary_files/"$i""$sex2".short
  done
  cat temporary_files/"$ancgroup"_chr*_p*.msoftout"$sex2".short | awk '{if (NR == 1 || $1 != "RSID") print}' > temporary_files/"$ancgroup"_allchr"$sex2".msoftoutpvs

 done

 #TBD: incorporate MAF and position and avg info into results

done

echo "SNP A1 A2 Direction Ncases NControls Neff NSTUDY PVALUE_FE BETA_FE STD_FE PVALUE_RE BETA_RE STD_RE PVALUE_RE2 ISQ Q " > reheader.txt
 #May want to not add in Neff until the end...
 
for gender in mf males females
do

 if [ $gender ==  "mf" ]
 then
  sex=""
  sex2=""
 fi
 
 if [ $gender !=  "mf" ]
 then
  sex=$gender
  sex2=_"$gender"
 fi
  
 for ancgroup in  eur #  lat eur aam all # eur  
 do
  #Add in A1/A2 alleles for each marker
  LC_ALL=C join  <( awk '{print $1,$2,$3,$11,$12,$13,$14}' results_cat/"$ancgroup"_dec28_2017_maf01_info6"$sex2"_noukbb.results_neff  | LC_ALL=C sort -g -k 1b,1 ) <(LC_ALL=C sort -k1b,1 temporary_files/"$ancgroup"_allchr"$sex2".msoftoutpvs) | cat reheader.txt  - > results_cat/"$ancgroup"_dec28_2017_maf01_info6"$sex2"_noukbb.msoftoutpvs_pos_neff
  awk '{if (NR==1) { pchoose="P"; bchoose="B" ;sechoose="SE"} else if($17 <= 0.05) { pchoose=$15 ;bchoose=$13; sechoose=$14} else {pchoose=$9 ; bchoose=$10; sechoose=$11}; print $1,$2,$3,bchoose,sechoose,pchoose,$17,$4,$5,$6,$7,$8}' results_cat/"$ancgroup"_dec28_2017_maf01_info6"$sex2"_noukbb.msoftoutpvs_pos_neff | sort -g -k 6 > results_cat/"$ancgroup"_dec28_2017_maf01_info6"$sex2"_noukbb.msoftoutpvs_pos_neff_re 
 done
done

##Transethnic random effects analysis
#Format data to metasoft format
#Am I sure that the A1 will be the same across analyses? It should, AS LONG AS, the first file used in each meta was one that I qced
mkdir metasoft_input_trans

for gender in males females # mf 
do

 if [ $gender ==  "mf" ]
 then
  sex=""
  sex2=""
 fi
 
 if [ $gender !=  "mf" ]
 then
  sex=$gender
  sex2=_"$gender"
 fi
 for ancgroup in eur aam lat 
 do
  zcat results_cat/"$ancgroup"_dec28_2017_maf01_info6"$sex2".results_neff.gz | awk '{print $1,$8,$9}'  | LC_ALL=C sort -k1b,1 > temporary_files/"$ancgroup"_dec28_2017_maf01_info6"$sex2".results_neff.msoft
 done
 
 LC_ALL=C join <(LC_ALL=C join temporary_files/eur_dec28_2017_maf01_info6"$sex2".results_neff.msoft temporary_files/aam_dec28_2017_maf01_info6"$sex2".results_neff.msoft ) temporary_files/lat_dec28_2017_maf01_info6"$sex2".results_neff.msoft | tail -n+2 | sed 's/0.0000/0.0001/g' > metasoft_input_trans/trans_dec28_2017_maf01_info6"$sex2".msoftdat.meta
 qsub -t1 -d $workingdir  -e errandout/ -o errandout/ run_metasoft.pbs -lwalltime=0:30:00  -F "-r trans -b xxx -g xxx"
done

#note: trans has never been done without UKBB
for gender in males females  mf 
do

 if [ $gender ==  "mf" ]
 then
  sex=""
  sex2=""
 fi
 
 if [ $gender !=  "mf" ]
 then
  sex=$gender
  sex2=_"$gender"
 fi
 
  awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$13,$15}' metasoft_results/trans_dec28_2017_maf01_info6"$sex2".msoftdat.meta.msoftout | grep -v NA > temporary_files/trans_dec28_2017_maf01_info6"$sex2""$sex2".short
  LC_ALL=C join  <( awk '{print $1,$2,$3,$11,$12,$13,$14}' results_cat/all_dec28_2017_maf01_info6"$sex2".results_neff  | LC_ALL=C sort -g -k 1b,1 ) <(LC_ALL=C sort -k1b,1 temporary_files/trans_dec28_2017_maf01_info6"$sex2""$sex2".short) | cat reheader.txt  - > results_cat/trans_dec28_2017_maf01_info6"$sex2".msoftoutpvs_pos_neff
  awk '{if (NR==1) { pchoose="P"; bchoose="B" ;sechoose="SE"} else if($17 <= 0.05) { pchoose=$15 ;bchoose=$13; sechoose=$14} else {pchoose=$9 ; bchoose=$10; sechoose=$11}; print $1,$2,$3,bchoose,sechoose,pchoose,$17,$4,$5,$6,$7,$8}' results_cat/trans_dec28_2017_maf01_info6"$sex2".msoftoutpvs_pos_neff | sort -g -k 6 > results_cat/trans_dec28_2017_maf01_info6"$sex2".msoftoutpvs_pos_neff_re 
done

#For final results, filter out the top SNPs that we found were not validated. If nothing needs to be changed, just make a shortcut..
zcat eur_dec28_2017_maf01_info6_females.msoftoutpvs_pos_neff.gz | grep -v rs9822553 | gzip > eur_dec28_2017_maf01_info6_females.msoftoutpvs_pos_nefff.gz 
zcat eur_dec28_2017_maf01_info6_males.msoftoutpvs_pos_neff.gz | grep -v rs9822553 | gzip > eur_dec28_2017_maf01_info6_males.msoftoutpvs_pos_nefff.gz 
zcat eur_dec28_2017_maf01_info6.msoftoutpvs_pos_neff.gz | grep -v rs9822553 | gzip > eur_dec28_2017_maf01_info6.msoftoutpvs_pos_nefff.gz 

zcat all_dec28_2017_maf01_info6_females.msoftoutpvs_pos_neff.gz | grep -v rs9822553 | gzip > all_dec28_2017_maf01_info6_females.msoftoutpvs_pos_nefff.gz 
zcat all_dec28_2017_maf01_info6_males.msoftoutpvs_pos_neff.gz | grep -v rs9822553 | gzip > all_dec28_2017_maf01_info6_males.msoftoutpvs_pos_nefff.gz 
zcat all_dec28_2017_maf01_info6.msoftoutpvs_pos_neff.gz | grep -v rs9822553 | gzip > all_dec28_2017_maf01_info6.msoftoutpvs_pos_nefff.gz 

zcat lat_dec28_2017_maf01_info6_females.msoftoutpvs_pos_neff.gz | grep -v chr12:48914460 | gzip > lat_dec28_2017_maf01_info6_females.msoftoutpvs_pos_nefff.gz 
zcat lat_dec28_2017_maf01_info6_males.msoftoutpvs_pos_neff.gz | grep -v chr12:48914460 | gzip > lat_dec28_2017_maf01_info6_males.msoftoutpvs_pos_nefff.gz 
zcat lat_dec28_2017_maf01_info6.msoftoutpvs_pos_neff.gz | grep -v chr12:48914460 | gzip > lat_dec28_2017_maf01_info6.msoftoutpvs_pos_nefff.gz 

ln -sf aam_dec28_2017_maf01_info6_females.msoftoutpvs_pos_neff.gz aam_dec28_2017_maf01_info6_females.msoftoutpvs_pos_nefff.gz 
ln -sf aam_dec28_2017_maf01_info6_males.msoftoutpvs_pos_neff.gz aam_dec28_2017_maf01_info6_males.msoftoutpvs_pos_nefff.gz 
ln -sf aam_dec28_2017_maf01_info6.msoftoutpvs_pos_neff.gz aam_dec28_2017_maf01_info6.msoftoutpvs_pos_nefff.gz 

#Manhattan and QQ plots

###Manhattan plot of results (list of things to plot is the input)
 ls results_cat | grep eur | grep nefff.gz$ > metal_outputsA.txt # echo eur_dec18_2017_maf01_info6.results_neff > metal_outputsA.txt
 ls results_cat | grep aam | grep nefff.gz$ > metal_outputsB.txt # echo aam_dec18_2017_maf01_info6.results_neff > metal_outputsB.txt
 ls results_cat | grep lat | grep nefff.gz$ > metal_outputsD.txt  # echo lat_dec18_2017_maf01_info6.results_neff > metal_outputsD.txt
 ls results_cat | grep all | grep nefff.gz$ > metal_outputsC.txt # echo all_dec18_2017_maf01_info6.results_neff > metal_outputsC.txt


 ncommands=$(wc -l metal_outputsA.txt | awk '{print $1}')
 nodesize=16
 nodeuse=$(($nodesize ))
 totjobs=$(( ($ncommands + $nodeuse - 1 ) / $nodeuse ))

 qsub -t1-$totjobs -l walltime=00:20:00 mh_plot_rapid_v2.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m mh_plot_pgc.R -s ManhattanPlotterFunction_colorfixed_max10ylim2_pgc.R -l phase3.locations2 -d metal_outputsA.txt -c green -p 0.05 -n 9"
 qsub -t1-$totjobs -l walltime=00:20:00 mh_plot_rapid_v2.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m mh_plot_pgc.R -s ManhattanPlotterFunction_colorfixed_max10ylim2_pgc.R -l phase3.locations2 -d metal_outputsB.txt -c red -p 0.05 -n 9"
 qsub -t1-$totjobs -l walltime=00:20:00 mh_plot_rapid_v2.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m mh_plot_pgc.R -s ManhattanPlotterFunction_colorfixed_max10ylim2_pgc.R -l phase3.locations2 -d metal_outputsC.txt -c blue -p 0.05 -n 9"
 qsub -t1-$totjobs -l walltime=00:20:00 mh_plot_rapid_v2.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m mh_plot_pgc.R -s ManhattanPlotterFunction_colorfixed_max10ylim2_pgc.R -l phase3.locations2 -d metal_outputsD.txt -c purple -p 0.05 -n 9" 


 
#Print p-value column of all data
 for files in $(ls results_cat | grep nefff.gz$ )
 do
  zcat results_cat/"$files" | awk '{print $15}'  | gzip > results_cat/"$files".p.gz
 done

#List all relevant files
 ls results_cat | grep nefff.gz.p.gz$ | sed 's/.gz$//g' > qq_plot.files

 qsub -l walltime=00:20:00 qq_plot.qsub -d $workingdir -e errandout/ -o errandout/ -F "-q qq_plot.R -e 1 -p qq_plot.files "

 #Make locus-zoom ready files - for regional plots of top hits
 
 qsub  -lwalltime=00:20:00 make_locuszoom.qsub -d $workingdir -e errandout/ -o errandout/ -F "-l phase3.locations2 -d metal_outputsA.txt -n 6"
 qsub  -lwalltime=00:20:00 make_locuszoom.qsub -d $workingdir -e errandout/ -o errandout/ -F "-l phase3.locations2 -d metal_outputsB.txt -n 6"
 #qsub  -lwalltime=00:20:00 make_locuszoom.qsub -d $workingdir -e errandout/ -o errandout/ -F "-l phase3.locations2 -d metal_outputsC.txt -n 6"
 qsub  -lwalltime=00:20:00 make_locuszoom.qsub -d $workingdir -e errandout/ -o errandout/ -F "-l phase3.locations2 -d metal_outputsD.txt -n 6"

 
#For Ben: Het test QQ plot - will do in eur and all
 for files in $(ls results_cat |  grep neff_ref.gz$ | grep lat )
 do
  zcat results_cat/"$files" | awk '{print $7}'  | gzip > results_cat/"$files"_het.p.gz
 done
 
  ls results_cat | grep het.p.gz$ | sed 's/.gz$//g' > qq_plot.files

 qsub -l walltime=00:20:00 qq_plot.qsub -d $workingdir -e errandout/ -o errandout/ -F "-q qq_plot.R -e 1 -p qq_plot.files "

 
###LDSC analyisis

 mkdir ldsc
 module load ldsc
 
##Filter markers to w_hm3 snplist, for LDSC (local
 
#FE Europeans
 LC_ALL=C join  <(awk '{print $1}' eur_w_ld_chr/w_hm3.noMHC.snplist.sorted) <(zcat eur_dec28_2017_maf01_info6.results_neff.gz  | awk '{if (NR == 1) $1="SNP"; print}'   | LC_ALL=C sort -k1b,1 ) | sort -g -k 10 | gzip > tempfiles/eur_dec28_2017_maf01_info6.results_neff.premunge.gz 
 LC_ALL=C join  <(awk '{print $1}' eur_w_ld_chr/w_hm3.noMHC.snplist.sorted) <(zcat eur_dec28_2017_maf01_info6_males.results_neff.gz    | awk '{if (NR == 1) $1="SNP"; print}'   | LC_ALL=C sort -k1b,1 ) | sort -g -k 10 | gzip > tempfiles/eur_dec28_2017_maf01_info6_males.results_neff.premunge.gz 
 LC_ALL=C join  <(awk '{print $1}' eur_w_ld_chr/w_hm3.noMHC.snplist.sorted) <(zcat eur_dec28_2017_maf01_info6_females.results_neff.gz   | awk '{if (NR == 1) $1="SNP";  print}'    | LC_ALL=C sort -k1b,1 ) | sort -g -k 10 | gzip > tempfiles/eur_dec28_2017_maf01_info6_females.results_neff.premunge.gz 

#FE Europeans, no UKBB
 LC_ALL=C join  <(awk '{print $1}' eur_w_ld_chr/w_hm3.noMHC.snplist.sorted) <(zcat eur_dec28_2017_maf01_info6_noukbb.results_neff.gz    | awk '{if (NR == 1) $1="SNP";  print}'   | LC_ALL=C sort -k1b,1 ) | sort -g -k 10 | gzip > tempfiles/eur_dec28_2017_maf01_info6_noukbb.results_neff.premunge.gz 
 LC_ALL=C join  <(awk '{print $1}' eur_w_ld_chr/w_hm3.noMHC.snplist.sorted) <(zcat eur_dec28_2017_maf01_info6_males_noukbb.results_neff.gz   | awk '{if (NR == 1) $1="SNP";  print}'    | LC_ALL=C sort -k1b,1 ) | sort -g -k 10 | gzip > tempfiles/eur_dec28_2017_maf01_info6_males_noukbb.results_neff.premunge.gz 
 LC_ALL=C join  <(awk '{print $1}' eur_w_ld_chr/w_hm3.noMHC.snplist.sorted) <(zcat eur_dec28_2017_maf01_info6_females_noukbb.results_neff.gz  | awk '{if (NR == 1) $1="SNP"; print}'     | LC_ALL=C sort -k1b,1 ) | sort -g -k 10 | gzip > tempfiles/eur_dec28_2017_maf01_info6_females_noukbb.results_neff.premunge.gz 

#UKBB
 LC_ALL=C join -1 1 -2 2 <(awk '{print $1}' eur_w_ld_chr/w_hm3.noMHC.snplist.sorted) <(zcat PTSD13_All_MAF1_INFO4_ricopili.txt.affix.logscale.gz  | LC_ALL=C sort -k2b,2 ) | sort -g -k 11 | gzip > tempfiles/PTSD13_All_MAF1_INFO4_ricopili.txt.affix.logscale.premunge.gz 
 LC_ALL=C join -1 1 -2 2 <(awk '{print $1}' eur_w_ld_chr/w_hm3.noMHC.snplist.sorted) <(zcat PTSD13_Males_MAF1_INFO4_ricopili.txt.affix.logscale.gz  | LC_ALL=C sort -k2b,2 ) | sort -g -k 11 | gzip > tempfiles/PTSD13_Males_MAF1_INFO4_ricopili.txt.affix.logscale.premunge.gz 
 LC_ALL=C join -1 1 -2 2 <(awk '{print $1}' eur_w_ld_chr/w_hm3.noMHC.snplist.sorted) <(zcat PTSD13_Females_MAF1_INFO4_ricopili.txt.affix.logscale.gz  | LC_ALL=C sort -k2b,2 ) | sort -g -k 11 | gzip > tempfiles/PTSD13_Females_MAF1_INFO4_ricopili.txt.affix.logscale.premunge.gz 

#Filter LDSC data for ldhub
zcat tempfiles/eur_dec28_2017_maf01_info6.results_neff.premunge.gz  | awk '{if(NR==1) print "SNP","A1","A2","MAF","BETA","P","N"; else print $1,$2,$3,$4,$8,$10,$12+$13}'   > eur_dec28_2017_maf01_info6.results_neff.ldhub
zcat tempfiles/eur_dec28_2017_maf01_info6_males.results_neff.premunge.gz  | awk '{if(NR==1) print "SNP","A1","A2","MAF","BETA","P","N"; else print $1,$2,$3,$4,$8,$10,$12+$13}'  > eur_dec28_2017_maf01_info6_males.results_neff.ldhub
zcat tempfiles/eur_dec28_2017_maf01_info6_females.results_neff.premunge.gz  | awk '{if(NR==1) print "SNP","A1","A2","MAF","BETA","P","N"; else print $1,$2,$3,$4,$8,$10,$12+$13}'   > eur_dec28_2017_maf01_info6_females.results_neff.ldhub

zcat tempfiles/eur_dec28_2017_maf01_info6_noukbb.results_neff.premunge.gz  | awk '{if(NR==1) print "SNP","A1","A2","MAF","BETA","P","N"; else print $1,$2,$3,$4,$8,$10,$12+$13}' > eur_dec28_2017_maf01_info6_noukbb.results_neff.ldhub
zcat tempfiles/eur_dec28_2017_maf01_info6_males_noukbb.results_neff.premunge.gz  | awk '{if(NR==1) print "SNP","A1","A2","MAF","BETA","P","N"; else print $1,$2,$3,$4,$8,$10,$12+$13}'  > eur_dec28_2017_maf01_info6_males_noukbb.results_neff.ldhub
zcat tempfiles/eur_dec28_2017_maf01_info6_females_noukbb.results_neff.premunge.gz  | awk '{if(NR==1) print "SNP","A1","A2","MAF","BETA","P","N"; else print $1,$2,$3,$4,$8,$10,$12+$13}'   > eur_dec28_2017_maf01_info6_females_noukbb.results_neff.ldhub

zip eur_dec28_2017_maf01_info6.results_neff.ldhub.zip eur_dec28_2017_maf01_info6.results_neff.ldhub
zip eur_dec28_2017_maf01_info6_males.results_neff.ldhub.zip eur_dec28_2017_maf01_info6_males.results_neff.ldhub
zip eur_dec28_2017_maf01_info6_females.results_neff.ldhub.zip eur_dec28_2017_maf01_info6_females.results_neff.ldhub

zip eur_dec28_2017_maf01_info6_noukbb.results_neff.ldhub.zip eur_dec28_2017_maf01_info6_noukbb.results_neff.ldhub
zip eur_dec28_2017_maf01_info6_males_noukbb.results_neff.ldhub.zip eur_dec28_2017_maf01_info6_males_noukbb.results_neff.ldhub
zip eur_dec28_2017_maf01_info6_females_noukbb.results_neff.ldhub.zip eur_dec28_2017_maf01_info6_females_noukbb.results_neff.ldhub

#Filter data for mrbase
 zcat results_cat/eur_dec28_2017_maf01_info6.results_neff.gz  | head -n 20 | gzip > results_cat/eur_dec28_2017_maf01_info6.msoftoutpvs_pos_nefff_mrbase.gz
 
#FE AAM
 LC_ALL=C join  <(awk '{print $1}' eur_w_ld_chr/w_hm3.noMHC.snplist.sorted) <(zcat aam_dec28_2017_maf01_info6.results_neff.gz     | LC_ALL=C sort -k1b,1 ) | sort -g -k 10 | gzip > tempfiles/aam_dec28_2017_maf01_info6.results_neff.premunge.gz 
 

#All  (Trans)
 LC_ALL=C join  <(awk '{print $1}' eur_w_ld_chr/w_hm3.noMHC.snplist.sorted) <(zcat trans_dec28_2017_maf01_info6.msoftoutpvs_pos_neff_re.gz     | LC_ALL=C sort -k1b,1 ) | sort -g -k 10 | gzip > tempfiles/trans_dec28_2017_maf01_info6.msoftoutpvs_pos_neff_re.premunge.gz 
 LC_ALL=C join  <(awk '{print $1}' eur_w_ld_chr/w_hm3.noMHC.snplist.sorted) <(zcat all_dec28_2017_maf01_info6.msoftoutpvs_pos_neff_re.gz     | LC_ALL=C sort -k1b,1 ) | sort -g -k 10 | gzip > tempfiles/all_dec28_2017_maf01_info6.msoftoutpvs_pos_neff_re.premunge.gz 
  
 
#Munge PGC FE data

 #PGC Eur + UKBB
 munge_sumstats.py --N-cas-col Ncases --N-con-col Ncontrols --sumstats  tempfiles/eur_dec28_2017_maf01_info6.results_neff.premunge.gz  --out tempfiles/eur_dec28_2017_maf01_info6.results_neff.munge.gz 
 munge_sumstats.py --N-cas-col Ncases --N-con-col Ncontrols --sumstats  tempfiles/eur_dec28_2017_maf01_info6_males.results_neff.premunge.gz  --out tempfiles/eur_dec28_2017_maf01_info6_males.results_neff.munge.gz 
 munge_sumstats.py --N-cas-col Ncases --N-con-col Ncontrols --sumstats  tempfiles/eur_dec28_2017_maf01_info6_females.results_neff.premunge.gz  --out tempfiles/eur_dec28_2017_maf01_info6_females.results_neff.munge.gz 
 
 #PGC Eur - UKBB
 munge_sumstats.py --N-cas-col Ncases --N-con-col Ncontrols --sumstats  tempfiles/eur_dec28_2017_maf01_info6_noukbb.results_neff.premunge.gz  --out tempfiles/eur_dec28_2017_maf01_info6_noukbb.results_neff.munge.gz 
 munge_sumstats.py --N-cas-col Ncases --N-con-col Ncontrols --sumstats  tempfiles/eur_dec28_2017_maf01_info6_males_noukbb.results_neff.premunge.gz  --out tempfiles/eur_dec28_2017_maf01_info6_males_noukbb.results_neff.munge.gz
 zcat tempfiles/eur_dec28_2017_maf01_info6_females_noukbb_MA.results_neff.munge.gz.sumstats.gz | awk '{if ($2 !="") print $1,$2,$3}' > tempfiles/eur_dec28_2017_maf01_info6_females_noukbb_MA.results_neff.alleles
 
 munge_sumstats.py --N-cas-col Ncases --N-con-col Ncontrols --sumstats  tempfiles/eur_dec28_2017_maf01_info6_males_noukbb.results_neff.premunge.gz   --merge-alleles tempfiles/eur_dec28_2017_maf01_info6_females_noukbb_MA.results_neff.alleles --out tempfiles/eur_dec28_2017_maf01_info6_males_noukbb_MA.results_neff.munge.gz
 
 munge_sumstats.py --N-cas-col Ncases --N-con-col Ncontrols --sumstats  tempfiles/eur_dec28_2017_maf01_info6_females_noukbb.results_neff.premunge.gz  --out tempfiles/eur_dec28_2017_maf01_info6_females_noukbb.results_neff.munge.gz 
 zcat tempfiles/eur_dec28_2017_maf01_info6_males_noukbb_MA.results_neff.munge.gz.sumstats.gz | awk '{if ($2 !="") print $1,$2,$3}' > tempfiles/eur_dec28_2017_maf01_info6_males_noukbb_MA.results_neff.alleles
 munge_sumstats.py --N-cas-col Ncases --N-con-col Ncontrols --sumstats  tempfiles/eur_dec28_2017_maf01_info6_females_noukbb.results_neff.premunge.gz  --merge-alleles tempfiles/eur_dec28_2017_maf01_info6_males_noukbb_MA.results_neff.alleles --out tempfiles/eur_dec28_2017_maf01_info6_females_noukbb_MA.results_neff.munge.gz 
 
#Munge UKBB
 munge_sumstats.py --N 126023 --sumstats  tempfiles/PTSD13_All_MAF1_INFO4_ricopili.txt.affix.logscale.premunge.gz  --out tempfiles/PTSD13_All_MAF1_INFO4_ricopili.txt.affix.logscale.munge.gz 

 munge_sumstats.py --N 55180  --sumstats  tempfiles/PTSD13_Males_MAF1_INFO4_ricopili.txt.affix.logscale.premunge.gz  --out tempfiles/PTSD13_Males_MAF1_INFO4_ricopili.txt.affix.logscale.munge.gz 
 #An error occurs calculating RG with women, so get alleles from PGC fem eur
 zcat tempfiles/eur_dec28_2017_maf01_info6_females_noukbb.results_neff.munge.gz.sumstats.gz | awk '{print $1,$2,$3}' > tempfiles/eur_dec28_2017_maf01_info6_females_noukbb.results_neff.alleles
 munge_sumstats.py --N 70843  --sumstats  tempfiles/PTSD13_Females_MAF1_INFO4_ricopili.txt.affix.logscale.premunge.gz --merge-alleles tempfiles/eur_dec28_2017_maf01_info6_females_noukbb.results_neff.alleles --out tempfiles/PTSD13_Females_MAF1_INFO4_ricopili.txt.affix.logscale.munge.gz 

#Munge AAM Fe
 munge_sumstats.py --N-col Neff --sumstats  tempfiles/aam_dec28_2017_maf01_info6.results_neff.premunge.gz   --out tempfiles/aam_dec28_2017_maf01_info6.results_neff.munge.gz 

 
#RG comparison:

#PGC Eur v UKBB (all)
 ldsc.py \
 --rg tempfiles/eur_dec28_2017_maf01_info6_noukbb.results_neff.munge.gz.sumstats.gz,tempfiles/PTSD13_All_MAF1_INFO4_ricopili.txt.affix.logscale.munge.gz.sumstats.gz \
 --ref-ld-chr eur_w_ld_chr/ \
 --w-ld-chr eur_w_ld_chr/ \
 --samp-prev 0.2645499,0.08222309 \
 --pop-prev 0.08,0.08 \
 --out ldsc/pgceur_ukbb_fe_rg

# (men)
 ldsc.py \
 --rg tempfiles/eur_dec28_2017_maf01_info6_males_noukbb.results_neff.munge.gz.sumstats.gz,tempfiles/PTSD13_Males_MAF1_INFO4_ricopili.txt.affix.logscale.munge.gz.sumstats.gz \
 --ref-ld-chr eur_w_ld_chr/ \
 --w-ld-chr eur_w_ld_chr/ \
 --samp-prev 0.2102481,0.06408119 \
 --pop-prev 0.05,0.05 \
 --out ldsc/pgceur_ukbb_males_fe_rg
 
# (women)
 ldsc.py \
 --rg tempfiles/eur_dec28_2017_maf01_info6_females_noukbb.results_neff.munge.gz.sumstats.gz,tempfiles/PTSD13_Females_MAF1_INFO4_ricopili.txt.affix.logscale.munge.gz.sumstats.gz \
 --ref-ld-chr eur_w_ld_chr/ \
 --w-ld-chr eur_w_ld_chr/ \
  --samp-prev  0.3914154,0.09635391 \
 --pop-prev 0.11,0.11 \
 --out ldsc/pgceur_ukbb_females_fe_rg
  
#RG comparison 2 men v women:

# with UKBB: men v women
 ldsc.py \
 --rg tempfiles/eur_dec28_2017_maf01_info6_males.results_neff.munge.gz.sumstats.gz,tempfiles/eur_dec28_2017_maf01_info6_females.results_neff.munge.gz.sumstats.gz \
 --ref-ld-chr eur_w_ld_chr/ \
 --w-ld-chr eur_w_ld_chr/ \
 --samp-prev 0.1158586,0.1761371 \
 --pop-prev 0.05,0.11 \
 --out ldsc/pgceurukbb_malesvfemales_fe_rg
 
#No ukbb: Men v women
 ldsc.py \
 --rg tempfiles/eur_dec28_2017_maf01_info6_males_noukbb_MA.results_neff.munge.gz.sumstats.gz,tempfiles/eur_dec28_2017_maf01_info6_females_noukbb_MA.results_neff.munge.gz.sumstats.gz \
 --ref-ld-chr eur_w_ld_chr/ \
 --w-ld-chr eur_w_ld_chr/ \
 --samp-prev 0.2102481,0.3914154 \
 --pop-prev 0.05,0.11 \
 --out ldsc/pgceur_ukbb_fe_rg

# UKBB only :men v women
 ldsc.py \
 --rg tempfiles/PTSD13_Males_MAF1_INFO4_ricopili.txt.affix.logscale.munge.gz.sumstats.gz,tempfiles/PTSD13_Females_MAF1_INFO4_ricopili.txt.affix.logscale.munge.gz.sumstats.gz \
 --ref-ld-chr eur_w_ld_chr/ \
 --w-ld-chr eur_w_ld_chr/ \
  --samp-prev  0.06408119,0.09635391 \
 --pop-prev 0.05,0.11 \
 --out ldsc/ukbb_malesvfemales_fe_rg
  

# UKBB men vs pgc women
 ldsc.py \
 --rg tempfiles/PTSD13_Males_MAF1_INFO4_ricopili.txt.affix.logscale.munge.gz.sumstats.gz,tempfiles/eur_dec28_2017_maf01_info6_females_noukbb_MA.results_neff.munge.gz.sumstats.gz  \
 --ref-ld-chr eur_w_ld_chr/ \
 --w-ld-chr eur_w_ld_chr/ \
  --samp-prev  0.06408119,0.3914154 \
 --pop-prev 0.05,0.11 \
 --out ldsc/ukbbmales_pgcfemales_fe_rg
  


### Cross validation analysis

 #Cross validated
 mkdir ldsc/cv
for num in {1..10}
do
 for group in a b
 do
  LC_ALL=C join  <(awk '{print $1}' eur_w_ld_chr/w_hm3.noMHC.snplist.sorted) <(zcat eur_dec28_2017_maf01_info6_cv"$group"_"$num".results_neff.gz     | LC_ALL=C sort -k1b,1 ) | sort -g -k 10 | gzip > tempfiles/eur_dec28_2017_maf01_info6_cv"$group"_"$num".results_neff.premunge.gz 
  munge_sumstats.py --N-col Neff --sumstats  tempfiles/eur_dec28_2017_maf01_info6_cv"$group"_"$num".results_neff.premunge.gz   --out tempfiles/eur_dec28_2017_maf01_info6_cv"$group"_"$num".results_neff.munge.gz 
 done
done
 
for num in {1..10}
do
 ldsc.py \
 --rg tempfiles/eur_dec28_2017_maf01_info6_cva_"$num".results_neff.munge.gz.sumstats.gz,tempfiles/eur_dec28_2017_maf01_info6_cvb_"$num".results_neff.munge.gz.sumstats.gz \
 --ref-ld-chr eur_w_ld_chr/ \
 --w-ld-chr eur_w_ld_chr/ \
 --out ldsc/cv_eur_"$num"_rg
done

#Get out summary info
rm ldsc/cv_eur_allrg.log
rm ldsc/cv_eur_allh2.log

for num in {1..10}
do
grep "Total Observed scale h2:" ldsc/cv_eur_"$num"_rg.log >> ldsc/cv_eur_allh2.log
grep "Genetic Correlation:"  ldsc/cv_eur_"$num"_rg.log >> ldsc/cv_eur_allrg.log
done

for group in a b
do
 for num in {1..10}
 do
  head -n2 eur_dec28_2017_maf01_info6_cv"$group"_"$num".results | tail -n1 | awk '{qcol=$11; gsub("?","",qcol); print length(qcol),$14}'  >> ldsc/cv_eur_all_ss"$group".log
 done
done


###Popcorn analysis

 #Mafs only necessary for ambiguous alleles
 
#Export python location
 export PYTHONPATH=$PYTHONPATH:/home/cnieverg/.local/lib/python2.7/site-packages
 
##Compare EA PGC to UKBB
 
#Compare PGC-PTSD to UKBB EAs only - to make this perfect, I should use the munged data. Beyond this, no reason to use hapmap3 munged data. Will not work for non white samples! 
 /home/cnieverg/.local/bin/popcorn fit -v 1 --cfile ../scores_euronly_allchr.txt --gen_effect --sfile1 tempfiles/eur_dec28_2017_maf01_info6_noukbb.results_neff.munge.gz.sumstats.gz --sfile2 tempfiles/PTSD13_All_MAF1_INFO4_ricopili.txt.affix.logscale.munge.gz.sumstats.gz --K1 .1 --P1 .5 --K2 .1 --P2 .5 pgcukbb_maf01_info6_fe_popcorn_matchldsc.txt
 #Results very similar to what I got without using hm3 anyway

 #popcorn reports .81 genetic correlation between datasets ( didnt check se). ldsc reports .71 (se .2). ref datasets are different but i think this about does it. h2 estimates for the two studies are alittle lower with popcorn (.027 and 0.08 comapred to something liek .037 and 0.10)
 
## Compare EA to AA
 
#Compare EA to AA, fixed effects results (Note: I think filtering to hapmap markers removed african topsnps! May not be good -other analysis had rg .54 or so. Indeed, rg goes down to ZERO due to filtering to these markers!
 /home/cnieverg/.local/bin/popcorn fit -v 1 --cfile ../scores_allchr.txt --gen_effect --sfile1 tempfiles/eur_dec28_2017_maf01_info6.results_neff.munge.gz.sumstats.gz  --sfile2 tempfiles/aam_dec28_2017_maf01_info6.results_neff.munge.gz.sumstats.gz  --K1 .1 --P1 .5 --K2 .1 --P2 .5 eur_aam_maf01_info6_neff_fe_popcorn.txt
 


 #Need mafs?
 zcat eur_dec28_2017_maf01_info6.msoftoutpvs_pos_neff_re.gz | awk '{OFS="\t"}{if(NR==1) print "SNP","a1","a2","af","N","beta","SE"; else if (  length($2) <=1 && length($3) <= 1) print $1,toupper($2),toupper($3),$4,$14,$8,$9}'  > eur_dec28_2017_maf01_info6.results_neff_popcorn
 zcat aam_dec28_2017_maf01_info6.msoftoutpvs_pos_neff_re.gz | awk '{OFS="\t"}{if(NR==1) print "SNP","a1","a2","af","N","beta","SE"; else if (  length($2) <=1 && length($3) <= 1) print $1,toupper($2),toupper($3),$4,$14,$8,$9}'  > aam_dec28_2017_maf01_info6.results_neff_popcorn
 /home/cnieverg/.local/bin/popcorn fit -v 1 --cfile ../scores_allchr.txt --gen_effect --sfile1 eur_dec28_2017_maf01_info6.results_neff_popcorn --sfile2 aam_dec28_2017_maf01_info6.results_neff_popcorn --K1 .1 --P1 .5 --K2 .1 --P2 .5 eur_aam_maf01_info6_fe_popcorn.txt
 
 
#LDSC

##With UKBB

#PGC Eur+Ukbb
 ldsc.py \
 --h2 tempfiles/eur_dec28_2017_maf01_info6.results_neff.munge.gz.sumstats.gz \
 --ref-ld-chr eur_w_ld_chr/ \
 --w-ld-chr eur_w_ld_chr/ \
 --samp-prev 0.1328699 \
 --pop-prev 0.08 \
 --out ldsc/pgceurukbb_fe_h2

# (men)
 ldsc.py \
 --h2  tempfiles/eur_dec28_2017_maf01_info6_males.results_neff.munge.gz.sumstats.gz \
 --ref-ld-chr eur_w_ld_chr/ \
 --w-ld-chr eur_w_ld_chr/ \
 --samp-prev 0.1158586  \
 --pop-prev 0.05 \
 --out ldsc/pgceurukbb_males_fe_h2
 
# (women)
 ldsc.py \
 --h2  tempfiles/eur_dec28_2017_maf01_info6_females.results_neff.munge.gz.sumstats.gz \
 --ref-ld-chr eur_w_ld_chr/ \
 --w-ld-chr eur_w_ld_chr/ \
 --samp-prev  0.1761371 \
 --pop-prev 0.11 \
 --out ldsc/pgceurukbb_females_fe_h2
  
 ## No ukbb
 
#men and women
 ldsc.py \
 --h2 tempfiles/eur_dec28_2017_maf01_info6_noukbb.results_neff.munge.gz.sumstats.gz \
 --ref-ld-chr eur_w_ld_chr/ \
 --w-ld-chr eur_w_ld_chr/ \
 --samp-prev 0.2645499 \
 --pop-prev 0.08 \
 --out ldsc/pgceur_fe_h2

# (men)
 ldsc.py \
 --h2  tempfiles/eur_dec28_2017_maf01_info6_males_noukbb.results_neff.munge.gz.sumstats.gz \
 --ref-ld-chr eur_w_ld_chr/ \
 --w-ld-chr eur_w_ld_chr/ \
 --samp-prev 0.2102481  \
 --pop-prev 0.05 \
 --out ldsc/pgceur_males_fe_h2
 
# (women)
 ldsc.py \
 --h2  tempfiles/eur_dec28_2017_maf01_info6_females_noukbb.results_neff.munge.gz.sumstats.gz \
 --ref-ld-chr eur_w_ld_chr/ \
 --w-ld-chr eur_w_ld_chr/ \
 --samp-prev  0.3914154 \
 --pop-prev 0.11 \
 --out ldsc/pgceur_females_fe_h2
  

  #UKBB alone
  #PGC Eur v UKBB (all)
 ldsc.py \
 --h2 tempfiles/PTSD13_All_MAF1_INFO4_ricopili.txt.affix.logscale.munge.gz.sumstats.gz \
 --ref-ld-chr eur_w_ld_chr/ \
 --w-ld-chr eur_w_ld_chr/ \
 --samp-prev 0.08222309 \
 --pop-prev 0.08 \
 --out ldsc/ukbb_fe_h2

# (men)
 ldsc.py \
 --h2 tempfiles/PTSD13_Males_MAF1_INFO4_ricopili.txt.affix.logscale.munge.gz.sumstats.gz \
 --ref-ld-chr eur_w_ld_chr/ \
 --w-ld-chr eur_w_ld_chr/ \
 --samp-prev 0.06408119 \
 --pop-prev 0.05 \
 --out ldsc/ukbb_males_fe_h2
 
# (women)
 ldsc.py \
 --h2 tempfiles/PTSD13_Females_MAF1_INFO4_ricopili.txt.affix.logscale.munge.gz.sumstats.gz \
 --ref-ld-chr eur_w_ld_chr/ \
 --w-ld-chr eur_w_ld_chr/ \
  --samp-prev  0.09635391 \
 --pop-prev 0.11 \
 --out ldsc/ukbb_females_fe_h2
  

 
 
 
 
#PGC Eur females
 ldsc.py \
 --h2 tempfiles/eur_dec28_2017_maf01_info6_females.results_neff.munge.gz.sumstats.gz \
 --ref-ld-chr eur_w_ld_chr/ \
 --w-ld-chr eur_w_ld_chr/ \
 --out ldsc/pgceur_women_fe_h2
 
  ldsc.py \
 --h2 tempfiles/eur_dec28_2017_maf01_info6_females.results_neff.munge.gz.sumstats.gz \
 --ref-ld-chr eur_ref_ld_chr/ \
 --w-ld-chr eur_ref_ld_chr/ \
 --out ldsc/pgceur_women_fe_h2_alt
 #Using w_ld (hapmap 3) vs ref_ld (1000g) results in slightly better h2. it was their recommendation to use w_ld
 
  ldsc.py \
 --h2 tempfiles/eur_dec28_2017_maf01_info6_females.results_neff.munge_ref.gz.sumstats.gz \
 --ref-ld-chr eur_ref_ld_chr/ \
 --w-ld-chr eur_ref_ld_chr/ \
 --out ldsc/pgceur_women_fe_h2_ref
 


  #Note: RE model gives lower genetic correlation, but higher H2 for PGC (with higher intercept)

 #just to see what happens with EA and AA


 
#Note: to use fuma output: Capitalize R2, remove space in genomic locus, get ridof markers with no associated p value

###LDSC (use kgp phase 3 refs)

cd /home/cnieverg/report_data/testrun/results_cat/
awk '{if(NR==1) print "SNP","A1","A2","MAF","BETA","P","N"; else print $1,$2,$3,$4,$8,$10,$12+$13}' eur_dec11_2017_maf01_info8.results > eur_dec11_2017_maf01_info8.results_ldsc
awk '{if(NR==1) print "SNP","A1","A2","MAF","BETA","P","N"; else print $1,$2,$3,$4,$8,$10,$12+$13}' eur_dec11_2017_maf01_info6.results > eur_dec11_2017_maf01_info6.results_ldsc

awk '{if(NR==1) print "SNP","A1","A2","MAF","BETA","P","N"; else print $1,toupper($2),toupper($3),$4,$8,$10,$14}' eur_dec28_2017_maf01_info6_females.results_neff > eur_dec28_2017_maf01_info6_females.results_neff_ldsc
zip eur_dec28_2017_maf01_info6_females.results_neff_ldsc.zip eur_dec28_2017_maf01_info6_females.results_neff_ldsc

awk '{if(NR==1) print "SNP","A1","A2","MAF","BETA","P","N"; else print $1,toupper($2),toupper($3),$4,$8,$10,$14}' eur_dec28_2017_maf01_info6_males.results_neff > eur_dec28_2017_maf01_info6_males.results_neff_ldsc
zip eur_dec28_2017_maf01_info6_males.results_neff_ldsc.zip eur_dec28_2017_maf01_info6_males.results_neff_ldsc


awk '{if(NR==1) print "SNP","A1","A2","Z","N"; else print $1,toupper($2),toupper($3),$4/$5,$17+$18}' results_cat/eur_dec18_2017.msoftoutpvs_pos_neff_re | grep -v chr* | grep -v "A T" | grep -v "T A" | grep -v "C G" | grep -v "G C" | gzip> results_cat/eur_dec18_2017.msoftoutpvs_pos_neff_re.results_ldsc.sumstats.gz
awk '{if(NR==1) print "SNP","A1","A2","Z","N"; else print $1,toupper($2),toupper($3),$4/$5,$17+$18}' results_cat/all_dec18_2017.msoftoutpvs_pos_neff_re | grep -v chr* | grep -v "A T" | grep -v "T A" | grep -v "C G" | grep -v "G C" | gzip> results_cat/all_dec18_2017.msoftoutpvs_pos_neff_re.results_ldsc.sumstats.gz


awk '{if(NR==1) print "SNP","A1","A2","Z","N"; else if($4 >= 0.01 ) print $1,toupper($2),toupper($3),($8/$9),$12+$13}' eur_dec11_2017_maf01_info9.results |  grep -v chr* | grep -v "A T" | grep -v "T A" | grep -v "C G" | grep -v "G C" | gzip > eur_dec11_2017_maf01_info9.results.sumstats.gz
awk '{if(NR==1) print "SNP","A1","A2","Z","N"; else if($4 >= 0.01 && $14 > 60000) print $1,toupper($2),toupper($3),($8/$9),$12+$13}' eur_dec11_2017_maf01_info9.results |  grep -v chr* | grep -v "A T" | grep -v "T A" | grep -v "C G" | grep -v "G C" | gzip > eur_dec11_2017_maf01_info9_nstudies.results.sumstats.gz

awk '{if(NR==1) print "SNP","A1","A2","Z","N"; else if($4 >= 0.01 ) print $1,toupper($2),toupper($3),($8/$9),$12+$13}' eur_dec11_2017_maf01_info8.results |  grep -v chr* | grep -v "A T" | grep -v "T A" | grep -v "C G" | grep -v "G C" | gzip > eur_dec11_2017_maf01_info8.results.sumstats.gz
awk '{if(NR==1) print "SNP","A1","A2","Z","N"; else if($4 >= 0.01 ) print $1,toupper($2),toupper($3),($8/$9),$12+$13}' eur_dec11_2017_maf01_info6.results |  grep -v chr* | grep -v "A T" | grep -v "T A" | grep -v "C G" | grep -v "G C" | gzip > eur_dec11_2017_maf01_info6.results.sumstats.gz
awk '{if(NR==1) print "SNP","A1","A2","Z","N"; else if($4 >= 0.01 && $14 > 60000 )print $1,toupper($2),toupper($3),($8/$9),$12+$13}' eur_dec11_2017_maf01_info6.results |  grep -v chr* | grep -v "A T" | grep -v "T A" | grep -v "C G" | grep -v "G C" | gzip > eur_dec11_2017_maf01_info6_nstudies.results.sumstats.gz


awk '{if(NR==1) print "SNP","A1","A2","B","SE","P"; else print $1,toupper($2),toupper($3),$8,$9,$10}' aam_dec28_2017_maf01_info6.results_neff | grep -v "A T" | grep -v "T A" | grep -v "C G" | grep -v "G C" | gzip > aam_dec28_2017_maf01_info6.results_neff.fuma.gz


awk '{if(NR==1) print "SNP","A1","A2","B","SE","P"; else print $1,toupper($2),toupper($3),$4,$5,$6}' results_cat/eur_dec18_2017.msoftoutpvs_pos_neff_re | grep -v "A T" | grep -v "T A" | grep -v "C G" | grep -v "G C" | gzip > results_cat/eur_dec18_2017.msoftoutpvs_pos_neff_re.fuma.gz
awk '{if(NR==1) print "SNP","A1","A2","Z","P-value","N"; else print $1,toupper($2),toupper($3),$4/$5,$6,$19}' results_cat/eur_dec18_2017.msoftoutpvs_pos_neff_re | grep -v "A T" | grep -v "T A" | grep -v "C G" | grep -v "G C" > results_cat/eur_dec18_2017_pneff_ldsc.txt
cd results_cat
zip eur_dec18_2017_pneff_ldsczip eur_dec18_2017_pneff_ldsc.txt

module load ldsc; munge_sumstats.py  --sumstats eur_dec11_2017_maf01_info8.results_ldsc  --out eur_dec11_2017_maf01_info8 
munge_sumstats.py  --sumstats eur_dec11_2017_maf01_info6.results_ldsc  --out eur_dec11_2017_maf01_info6 


ldsc.py \
--h2 results_cat/all_dec18_2017.msoftoutpvs_pos_neff_re.results_ldsc.sumstats.gz \
--ref-ld-chr results_cat/eur_w_ld_chr/ \
--w-ld-chr results_cat/eur_w_ld_chr/ \
--samp-prev 0.125342 \
--pop-prev 0.08 \
--out results_cat/all_dec18_2017.msoftoutpvs_pos_neff_re_bestpv



ldsc.py \
--h2 eur_dec11_2017_maf01_info8.results.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--samp-prev 0.125342 \
--pop-prev 0.08 \
--out eur_dec11_2017_maf01_info8

ldsc.py \
--h2 eur_dec11_2017_maf01_info6.results.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--samp-prev 0.125342 \
--pop-prev 0.08 \
--out eur_dec11_2017_maf01_info6

ldsc.py \
--h2 eur_dec11_2017_maf01_info6_nstudies.results.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--samp-prev 0.125342 \
--pop-prev 0.08 \
--out eur_dec11_2017_maf01_info6_nstudies

ldsc.py \
--h2 eur_dec11_2017_maf01_info9_nstudies.results.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--samp-prev 0.125342 \
--pop-prev 0.08 \
--out eur_dec11_2017_maf01_info9_nstudies

### Credible set analysis (get r2 from MAGMA)
Rscript gedcredible.R zdhhc14.txt
Rscript gedcredible.R park2.txt

###Look up genotyping quality of top hits
 ls starting_data > all_input_data.txt
 for i in $(cat all_input_data.txt)
 do
  echo "Checking file $i"
  zgrep -w -m1 $snp starting_data/"$i" >> check/rs4612196.assoc.dosage
 done


#Comparison of info filters:
#Did not markers filter subjects based on N:
## 0.0396 (0.0057) on info 9 data
##0.0371 (0.0057) for the info 8 data (on ldhub 0.0161 (0.0025))
##0.0362 (0.0061) for the info 6 data
#Marginal difference, but we preseve a hit at info .6
#Setting the info 6 data so that Neff must be > 60000:  0.0338 (0.0063)
#Setting the info 9 data so that Neff must be > 60000:  0.0341 (0.0065)
#In other words, this is like filtering on N good studies - removes differences between info filters

#
