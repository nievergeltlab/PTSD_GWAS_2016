for snp in $(cat mega.snplist)
do
 grep  -m1 -w $snp results_cat/eur_39.results >> results_cat/cgd_mega.results
done

#Also get the marker info
mkdir mega_info

for study in $(cat european_studylist.txt)
do

 zgrep -f mega.snplist -w  snp_info/"$study"* >> mega_info/"$study".megainfoscores
done

cd mega_info

echo "SNP FRQ_A FRQ_U INFO ngt" > info_header.txt

for files in $(ls | grep .megainfoscores$)
do
 cat info_header.txt <(awk -F : '{$1=$1; print $2}' $files )  > "$files"_v2
done

#I dont have vetsa, qimr, or external study site info scores!

#Now we have to join all of the results files

R
datasheets=system('ls | grep .megainfoscores_v2',intern=T)
library(plyr)


#read each file into a data frame with the same name
for (i in datasheets)
{
	assign(
		i, read.table(paste(i), header=T, stringsAsFactors=F)[,c(1,4)]
 
		) 
}

#parse the text list of data frame names as a list of data frames
data_list <- eval( 
			parse( 
				text=paste(
					"list(", paste(datasheets, collapse=','), ")" 
					)
				)
			)
   
  datA <- join_all(data_list,by="SNP", type="left")
  outdat <- data.frame(datA$SNP,t(apply(datA[,-1],1,quantile,na.rm=T)))
  names(outdat) <- c("SNP","min","x25pct","median","x75pct","max")
  write.table(outdat,"mega_infoscores.txt",quote=F,row.names=F)

