##script R
fabrique_outlier.R
##script bash
mkdir mes_outliers
COUNTER=1
for monvcf in `ls mes_vcf_save/chr[1-9].vcf`;
do
for outl in `cat nom_84_outliers.txt`;
do
grep $outl $monvcf | awk '{ print $2}'  >> "mes_outliers/chr"$COUNTER"_outliers.list";
done
let COUNTER=COUNTER+1 
done



