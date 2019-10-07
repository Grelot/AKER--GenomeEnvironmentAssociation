
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
parser.add_argument("-scol","--selectcol", type=str)
parser.add_argument("-csv","--csvFile",type=str)


#================================================================================
#MAIN
#================================================================================
args = parser.parse_args()
selectcol = args.selectcol
csvFile = args.csvFile

listOfScol = []
with open(selectcol,'r') as colFile:
	for ligne in colFile.readlines():
		listOfScol.append(ligne.split()[0])
colFile.close()

#prendre les ID des colonnes a selectionner
listOfID = []
id_h=0
with open(csvFile,'r') as readFile:
	headigne = readFile.readlines()[0]
	for h in headigne.split(";"):
		if h in listOfScol:
			listOfID.append(id_h)
		id_h+=1
readFile.close()

full_text=""
with open(csvFile,'r') as readFile:
	for ligne in readFile.readlines():
		ligneSplit = ligne.split(";")
		for h_id in listOfID:
			full_text+=ligneSplit[h_id]+";"
		full_text+="\n"

print full_text
		


			
