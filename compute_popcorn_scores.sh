
#Convert all 1000g data to plink binary
for chr in {1..22}
do
 plink2 --vcf 1000g/ALL.chr"$chr".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep 1000g_eur.ref2 --make-bed --out 1000g/1000geur_chr"$chr"
 plink2 --vcf 1000g/ALL.chr"$chr".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep 1000g_afr.ref2 --make-bed --out 1000g/1000gafr_chr"$chr"
done



for chr in {1..22}
 do
 #Remove SNPs that appear multiple times!
 awk '{print $2}' 1000g/1000geur_chr"$chr".bim | sort -k1 | uniq -cd | awk '{print $2}' > 1000g/chr"$chr".duplicated 
 #Get allele freqs
 plink2 --bfile 1000g/1000geur_chr"$chr"_fixed --exclude 1000g/chr"$chr".duplicated --freq --make-bed --out 1000g/1000geur_chr"$chr"_fixed2
 plink2 --bfile 1000g/1000gafr_chr"$chr"_fixed --exclude 1000g/chr"$chr".duplicated --freq --make-bed --out 1000g/1000gafr_chr"$chr"_fixed2

 #If the data has been fixed, removed originals
  if [ -e 1000g/1000geur_chr"$chr"_fixed2.bed ]
  then
   rm -f 1000g/1000gafr_chr"$chr"_fixed.bed 1000g/1000geur_chr"$chr"_fixed.bed
  fi
done
 
#Compute scores for all markers in EA and AA samples
for chr in {1..22}
do
 popcorn compute -v 1 --SNPs_to_store 2000000 --bfile1 1000g/1000geur_chr"$chr"_fixed2 --bfile2 1000g/1000gafr_chr"$chr"_fixed2 --gen_effect scores_chr"$chr".txt 
done


#Concatenate scores
 cat  scores_chr*.txt | grep -P -v "C\tG"| grep -P -v "A\tT" | grep -P  -v "T\tA" | grep  -P -v "G\tC"  >  scores_allchr.txt

#Cat AFs for all chr
 cat 1000g/1000geur_chr*_fixed.frq | awk '{print $2,$5}' > 1000g/1000geur_allchr_fixed.frq
 cat 1000g/1000gafr_chr*_fixed.frq  | awk '{print $2,$5}'  > 1000g/1000gafr_allchr_fixed.frq

#Need to merge in allele freqs (no longer needed as I put these reusltso ut of meta analysis)
 LC_ALL=C join <(LC_ALL=C sort -k1b,1 eur_39.results)  <( cat 1000g/1000geur_allchr_fixed.frq | LC_ALL=C sort -k1b,1 ) > eur_39.fixed
 awk '{OFS="\t"}{if (NR==1) print "SNP",	"a1"	,"a2"	,"af",	"N"	,"beta",	"SE"; else print $1, toupper($2),toupper($3),$8,"48000",$4,$5}' eur_39.fixed > eur_39.popcorn
 eur_39.popcorn > eur_39.popcorn2
 

 LC_ALL=C join <(LC_ALL=C sort -k1b,1 aam_39.results)  <( cat  1000g/1000gafr_allchr_fixed.frq | LC_ALL=C sort -k1b,1 ) > aam_39.fixed
 awk '{OFS="\t"}{if (NR==1) print "SNP",	"a1"	,"a2"	,"af",	"N"	,"beta",	"SE"; else print $1, toupper($2),toupper($3),$8,"20000",$4,$5}' aam_39.fixed > aam_39.popcorn
 awk '{if(NR == 1 || (length($2) <=1 && length($3) <= 1)) print}' aam_39.popcorn > aam_39.popcorn2
 

popcorn fit -v 1 --cfile scores_allchr.txt --gen_effect --sfile1 eur_39.popcorn2 --sfile2 aam_39.popcorn2 eur39_aam39_popprev.ldsc --K1 .1 --P1 .25 --K2 .1 --P2 .25  
  

 
#12813","35640
#4387,15405

#LDSC of data

python ldsc-master/munge_sumstats.py --sumstats eur_39.results --N 48000 --out eur_39

python ldsc-master/ldsc.py \
--h2 eur_39.sumstats.gz \
--ref-ld-chr baseline_v1.1/baseline2. \
--w-ld-chr baseline_v1.1/baseline2. \
--samp-prev 0.25  \
--pop-prev 0.08 \
--out merged

python ldsc-master/ldsc.py \
--h2 eur_39.sumstats.gz \
--ref-ld-chr /mnt/sdb/genetics/ukbbmeta/eur_w_ld_chr/ \
--w-ld-chr /mnt/sdb/genetics/ukbbmeta/eur_w_ld_chr/ \

--out merged2



--overlap-annot \

--frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. \
--annot baseline_v1.1/baseline. \
baseline_v1.1/baseline.1.l2.ldscore
baseline_v1.1/baseline1.l2.ldscore


#Estimate LD scores
for chr in {1..22}
do
python ldsc-master/ldsc.py \
	--bfile 1000g/1000geur_chr"$chr"_fixed2 \
	--l2 \
	--ld-wind-cm 1 \ 
	--out 1000gphase3ld/1000geur_chr"$chr"
 done
 
for chr in {1..22}
do
zcat baseline_v1.1/baseline."$chr".l2.ldscore.gz | awk '{print $1,$2,$3,$4}' | gzip > baseline_v1.1/baseline2."$chr".l2.ldscore.gz 
zcat baseline_v1.1/baseline."$chr".annot.gz | awk '{print $1,$2,$3,$4}' | gzip > baseline_v1.1/baseline2."$chr".annot.gz 
awk '{print $1}' baseline_v1.1/baseline."$chr".l2.M_5_50 > baseline_v1.1/baseline2."$chr".l2.M_5_50
 done
 

