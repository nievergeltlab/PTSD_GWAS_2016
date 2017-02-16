#!/bin/bash
while getopts m:s:l:d:o:c:p:d: option
do
  case "${option}"
    in
      m) mhplotscript=${OPTARG};;
      s) graphicsscript=${OPTARG};;
      l) snplocations=${OPTARG};;
      d) metaresults=${OPTARG};;
      o) plotname=${OPTARG};;
      c) color=${OPTARG};;
      p) pvalinclude=${OPTARG};;   
      z) workingdir=${OPTARG};;   
    esac
done

#Sample syntax
#bash mh_plot_v3.sh -m mh_plot_v2.r -s ManhattanPlotterFunction_colorfixed_max10ylim2.R -l phase3.locations2 -d /mnt/sdb/genetics/manhattan/lisa_nov4/results/istss_v11.tbl -o istss_v11.tbl.jpg -z /mnt/sdb/genetics/manhattan/ -c blue -p 0.05


cd $workingdir

echo "Filtering data to columns 1 and 6 (metal defaults). Only take results with p < $pvalinclude"
 cut -f 1,6 $metaresults | awk -v pvalinclude=$pvalinclude '{if ($2 < pvalinclude) print}' | LC_ALL=C sort -k1,1b > "$metaresults".temp.txt

echo "Joining to reference locations"
 LC_ALL=C join $snplocations "$metaresults".temp.txt > "$metaresults".temp2.txt

echo "Starting Manhattan plotting script"
Rscript $mhplotscript $graphicsscript "$metaresults".temp2.txt $plotname $color

echo "Cleanup"
 #rm "$metaresults".temp.txt
 #rm "$metaresults".temp2.txt
