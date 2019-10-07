#creer des fichiers VCF lisible sur le package R PopGenome

##STEP1 fabriquer les VCFs et avoir les infos sur les chr
#sort par chrom et scaffold_ID et pos
sort -k 2 -V data/Positions.14409.txt > data/Positions_scaff_sorted.txt
#ne garder que 950 indv dans le fichier CSV
python get_col.py -scol data/Data.950.sauvages.txt -csv data/NoPool.14409.csv > data/950_14409.csv

#fabriquer les VCF
python convert_data2vcf.py -f data/950_14409.csv -f2 data/Positions_scaff_sorted.txt -vcf bon_exemple.vcf -o mes_vcf > infos_chr.txt

cd mes_vcf
cp * ../mes_vcf_save/

##STEP2 convertir les VCFs au bon format pour PopGenome
for file_vcf in `ls chr*.vcf`
do
	echo $file_vcf
	bgzip $file_vcf
	vcf-sort $file_vcf.gz > $file_vcf.sorted
	bgzip $file_vcf.sorted
	tabix -p vcf $file_vcf.sorted.gz
done

cd ..

##STEP3 lire les VCFs dans R PopGenome

lire_vcf_PopGenome.R

