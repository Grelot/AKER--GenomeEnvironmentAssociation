#================================================================================
#AUTHORS
#================================================================================
'''
GUERIN Pierre-Edouard
France, Montpellier
CNRS, CEFE
05/2017
'''
#================================================================================
#NOTICE
#================================================================================

'''
convert CSV file into list of individual allele A (SNPs) sequence
for homozygote get A or B allele
for heterozygote get [randomly] A or B allele
and put all the data into VCF FORMAT able to be read by Rpackage PopGenome

'''

#================================================================================
#MODULES
#================================================================================

import re
import sys
import subprocess
import argparse
import os.path
import random

#================================================================================
#CLASS
#================================================================================

class SNP:
	def __init__(self,nom,chrom,pos,alleleRef,alleleAlt):
		self.nom = nom
		self.dicOfIndv={}
		self.chrom = chrom
		self.pos = pos
		self.alleleRef=alleleRef
		self.alleleAlt=alleleAlt
	def rename_snp(self,nom):
		self.nom = nom
	def add_indv_GT(self,indv_name,genoGT):
		self.dicOfIndv[indv_name] = genoGT



class Chromosome:
	def __init__(self,numer):
		self.numero = numer
		self.dicOfSnp = {}
	def append_snp(self,name_snp, snp):
		self.dicOfSnp[name_snp] = snp
	def set_stats(self,minimum,maximum,length):
		self.minimum = minimum
		self.maximum = maximum
		self.length = length
	def remove_snp(self,snp_name):
		del self.dicOfSnp[snp_name]
	def update_snp(self,snp_name,new_snp):
		self.dicOfSnp[snp_name]=new_snp


	
	

#================================================================================
#FUNCTIONS
#================================================================================


#================================================================================
#ARGUMENTS
#================================================================================

parser = argparse.ArgumentParser(description='specie seq geo')
parser.add_argument("-o","--output", type=str)
parser.add_argument("-f","--inputFile",type=str)
parser.add_argument("-f2","--inputFile2",type=str)
parser.add_argument("-vcf","--vcfbase",type=str)
parser.add_argument("-f3","--inputFile3",type=str)

#================================================================================
#MAIN
#================================================================================
args = parser.parse_args()
resultsPath = args.output
inputFile = args.inputFile
inputFile2=args.inputFile2
vcfbase=args.vcfbase

previous_scaff=0
shift_pos = 0
previous_chrom=0
previous_scaff_pos=0
#recuperer la liste des coords CHR:POS des SNPs
with open(inputFile2,'r') as readFile2:	
	dicOfSNP = {}
	listOfChrom = []
	#initialiser la liste de chromosome
	for ch_numero in xrange(1,10):		
		listOfChrom.append(Chromosome(ch_numero))
	for ligne in readFile2.readlines()[1:]:
		ligneSplit=ligne.split()
		nom=ligneSplit[0]
		chrom=ligneSplit[1]
		scaff=ligneSplit[2]
		scaff_pos=ligneSplit[3]
		if chrom != "NA" and chrom != "Chromosome" and scaff_pos !="NA" and scaff != "NA" and scaff != "Scaffold":
			if previous_chrom != chrom:
				shift_pos=0
				previous_scaff=0
				previous_scaff_pos=0
			if previous_scaff != scaff:
				shift_pos+=int(previous_scaff_pos)+150000
			previous_scaff=scaff
			pos=int(scaff_pos)+shift_pos
			previous_scaff_pos=int(scaff_pos)
			previous_chrom=chrom
			dicOfSNP[nom] = SNP(nom,chrom,pos,"-","-")			
			listOfChrom[int(chrom)-1].append_snp(nom,dicOfSNP[nom])	
readFile2.close()




#A ce stade on a la liste des coords CHR:POS par nom pour chaque SNP
header=1
nb_chr=len(listOfChrom)
nb_chr=9
for ch_id in xrange(0,nb_chr):
	#on ouvre le fichier des genotypes par nom pour chaque indv pour chaque SNP
	with open(inputFile,'r') as readFile:
		for ligne in readFile.readlines():			
			ligneSplit=ligne.split(";")
			if header==1:
				official_name_indv=[]
				for chaque_indv in ligneSplit[10:-1]:
					official_name_indv.append(chaque_indv)
				header=0
			nom = ligneSplit[0]			
			#si le SNP est present sur le CHR
			if nom in listOfChrom[ch_id].dicOfSnp:
				alleleRef = str(ligneSplit[1])
				alleleAlt= str(ligneSplit[2])
				if alleleRef!="-" and alleleAlt !="-":
					#alleles = [alleleRef,alleleAlt]
					#genot de chaque indv
					chrom=listOfChrom[ch_id].dicOfSnp[nom].chrom
					pos=listOfChrom[ch_id].dicOfSnp[nom].pos
					new_snp=SNP(nom,chrom,int(pos),alleleRef,alleleAlt)
					listOfChrom[ch_id].dicOfSnp[nom]=new_snp								
					id_indv=0
					for indvGeno in ligneSplit[10:-1]:
						if indvGeno =="AA":
							genoGT="0|0"
						elif indvGeno == "BB":
							genoGT="1|1"
						else:
							genoGT = random.choice(["0|1","1|0"])
						indv_nom=str(official_name_indv[id_indv])
						listOfChrom[ch_id].dicOfSnp[nom].add_indv_GT(indv_nom, genoGT)
						id_indv+=1
				else:
					listOfChrom[ch_id].remove_snp(nom)
					#suprimer le SNP dont le genotype estmal defini								
	readFile.close()



#informations max/min/numberofSNP
numero_chr=1
for ch in listOfChrom:
	listOfSnp=[]
	for key,value in ch.dicOfSnp.items():
		listOfSnp.append(value.pos)
	sortedSnp = sorted(listOfSnp)
	chminimum =str(sortedSnp[0])
	chmaximum =str(sortedSnp[-1])
	chlength=str(len(sortedSnp))
	print "chr"+str(numero_chr), "minpos:"+chminimum+" |", "maxpos:"+chmaximum+" |","numberofSNPs:"+chlength
	numero_chr+=1



#ecrire les fichiers VCF

vcfheader=open(vcfbase, 'r').read()

numero_chr=1
indvheader =""
for indv in official_name_indv:
	indvheader+="{0}\t".format(indv)
indvheader=indvheader[0:-1]
indvheader+="\n"

vcftextbase=vcfheader[0:-1]+indvheader


for ch in listOfChrom:	
	#print key, value.chrom, value.pos, value.dicOfIndv["M03212"]
	outputFile=resultsPath+"/chr"+str(numero_chr)+".vcf"
	with open(outputFile,'w') as writeFile:
		vcftext=vcftextbase						
		for key, value in ch.dicOfSnp.items():			
			vcftext+="{0}\t{1}\t{2}\t".format(value.chrom,value.pos,value.nom)
			vcftext+="{0}\t{1}\t.\t.\t.\tGT".format(value.alleleRef,value.alleleAlt)
			for indv_nom in official_name_indv:
				vcftext+="\t{0}".format(value.dicOfIndv[indv_nom])
			vcftext+="\n"
		writeFile.write(vcftext)
	writeFile.close()
	numero_chr+=1
print "Done"





#for i in xrange(0,len(snps[0])):
#	print snps[0][0]
	



