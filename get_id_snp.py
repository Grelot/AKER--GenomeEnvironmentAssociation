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
get coordinates in a file of SNP and search his matching ID in a VCF file

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
#ARGUMENTS
#================================================================================

parser = argparse.ArgumentParser(description='specie seq geo')
parser.add_argument("-o","--output", type=str)
parser.add_argument("-f","--inputFile",type=str)
parser.add_argument("-f2","--inputFile2",type=str)



#================================================================================
#MAIN
#================================================================================
args = parser.parse_args()
outputFile = args.output
inputFile = args.inputFile
inputFile2=args.inputFile2


#get list of couple SNP coordinates/SNP ID
listOfSNP=[]
with open(inputFile,'r') as readFile:
	for ligne in readFile.readlines():
		if ligne[0] != "#":
			ligneSplit=ligne.split()
			snpCoords = int(ligneSplit[1])
			snpID = ligneSplit[2]
			snpInfo=(snpCoords,snpID)			
			listOfSNP.append(snpInfo)
readFile.close()



listOfMatchedSNP = [("x",-1)]


text=""
#get list of SNP coordinates that I want to know their ID
with open(inputFile2,'r') as readFile2:
	for ligne in readFile2.readlines()[1:]:
		ligneSplit=ligne.split(",")
		unknownSnpCoords = int(ligneSplit[2])
		for i in listOfSNP:			
			if i[0] == unknownSnpCoords:
				if i not in listOfMatchedSNP:
					text+="{0},{1},{2},{3},{4},{5},{6},{7}".format(i[1],ligneSplit[1],ligneSplit[2],ligneSplit[3],ligneSplit[4],ligneSplit[5],ligneSplit[6],ligneSplit[7])
					listOfMatchedSNP.append((i[0],i[1]))
					break
readFile2.close()


with open(outputFile,'w') as writeFile:
	writeFile.write(text)
writeFile.close()
	
	




			
		
		








