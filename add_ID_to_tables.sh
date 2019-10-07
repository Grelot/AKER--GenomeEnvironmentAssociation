#mkdir tables/avec_id
#Get IA des SNP
COUNTER=1
for vcfFichier in `ls mes_vcf/chr*.vcf.gz`
do
	python get_id_snp.py -f <(zcat $vcfFichier) -f2 "tables/chr"$COUNTER"_stats.csv" -o "tables/avec_id/chr"$COUNTER"_stats_id.csv"
	let COUNTER=COUNTER+1
done

#merge all tables
cat tables/avec_id/chr*_stats_id.csv >> tables/avec_id/all_stats.csv 
#enlever les NA
grep -v 0,0,0 tables/avec_id/all_stats.csv > tables/avec_id/all_stats_propre.csv


