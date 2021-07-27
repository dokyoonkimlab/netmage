# Author: Vivek Sriram
# Last revised: 04/27/2021
# ------------------------------------------------
# Function: Script to prepare input data for Gephi.
#     The script will process summary PheWAS data and
#	  generate a dictionary of phenotypes mapped to
#     lists of lists. Each list contains lists of
#     SNP name, p-value, and MAF. If the user doesn't
#	  provide a column name for MAF or for p-value,
#	  then the corresponding information just won't
#	  be included
#

# Import statements
import os
import argparse
import sys
import glob

# Read in user arguments
parser = argparse.ArgumentParser(description='Take in user input for creation of node and edge maps for Gephi.')
parser.add_argument('--input-files', dest= 'inputfiles', default = None, required = True, metavar='input data (regex)', help='Regular expression for input data. Ex: "dir/*.csv"')
parser.add_argument('--datafile-output', dest='datafileout', default = "", metavar='data output file', help='Filepath for data output file (tsv)')
parser.add_argument('--phenotype-name', dest='phenotypename', metavar = 'phenotype colname', help='Name of column corresponding to phenotype')
parser.add_argument('--snp-name', dest='snpname', required=True, metavar = 'SNP ID colname', help='Name of column corresponding to SNP ID')
parser.add_argument('--maf-name', dest='mafname', metavar = 'MAF colname', help='Name of column corresponding to minor allele frequency')
parser.add_argument('--pvalue-name', dest='pvname', metavar = 'p-value colname', help='Name of column corresponding to p-value of association')
parser.add_argument('--delim', dest='delim', required=True, metavar = 'delimiter', help='Delimiter used in input PheWAS data file(s)')

args = parser.parse_args()

# Example usage:
# python makeIntermediateNodeMap.py --input-files "/Users/viveksrm/Documents/GitHub/peGBSSL/data/rawPheWASData/ldPrunedUKBB_noSymptoms_rawData/*.txt" 
# --datafile-output "./intermediateFile.csv" 
# --phenotype-name "phenotypeID" 
# --snp-name "ID" 
# --maf-name "af" 
# --pvalue-name "pval" 
# --delim " "

inputData = args.inputfiles
datafileOutputPath = args.datafileout

phenotypeColName = args.phenotypename
snpColName = args.snpname
mafColName = args.mafname
pValueColName = args.pvname

fileDelimiter = args.delim

if (fileDelimiter == "\\t"):
	fileDelimiter = '\t'
elif (fileDelimiter == "\\s"):
	fileDelimiter = ' '

if(fileDelimiter == None):
	sys.exit('ERROR: No file delimiter provided in input')

def getColumnNumbers(fileName, line):
	phenotypeColNumber = None
	mafColNumber = None
	pValueColNumber = None
	snpColNumber = None

	firstLine = line.strip()
	firstLineAsArray = firstLine.split(fileDelimiter)
	# Figure out the column indices corresponding to phenotype, MAF, case count, p-value, and SNP (if they exist)
	for i in range(0, len(firstLineAsArray)):
		if(firstLineAsArray[i] == phenotypeColName):
			phenotypeColNumber = i
		elif(mafColName != None and firstLineAsArray[i] == str(mafColName)):
			mafColNumber = i
		elif(pValueColName != None and firstLineAsArray[i] == str(pValueColName)):
			pValueColNumber = i
		elif(firstLineAsArray[i] == snpColName):
			snpColNumber = i

	# We need a column corresponding to snp name in our data. Force an exit if it can't be found
	if(snpColNumber == None):
		print("Header line is " + str(firstLine))
		sys.exit('ERROR: Unable to identify column corresponding to SNP name in file ' + fileName + '. Confirm that the right file has been uploaded and that the header is correct.')

	print("Identified column indices from " + str(fileName))
	return phenotypeColNumber, mafColNumber, pValueColNumber, snpColNumber



def updatePhenotypeSNPMapFromData(fileName, fileContent, phenotypeSnpMap, phenotypeColNumber, mafColNumber, pValueColNumber, fileDelimiter):
	# Save out the file name as the phenotype for now
	phenotype = os.path.splitext(os.path.basename(fileName))[0]
	print(phenotype)

	for line in fileContent:
		# Strip the line of any unnecessary characters
		line = line.strip()
		# Get the elements of each line
		lineAsArray = line.split(fileDelimiter)
		# If the line is not empty, process it
		if(lineAsArray != ['']):
			# If there's a phenotype column, get the name of the phenotype
			if(phenotypeColNumber != None):
				phenotype = lineAsArray[phenotypeColNumber]

			# If this is not the first line in the file
			if(lineAsArray[0] != 'ID'):
				# Get the name of the SNP
				currentSNP = lineAsArray[snpColNumber]
				currentValues = [currentSNP]

				if(pValueColNumber != None):
					currentPValue = lineAsArray[pValueColNumber]
					currentPValueFloat = float(currentPValue)
					currentValues.append(currentPValueFloat)
				if(mafColNumber != None):
					currentMAF = lineAsArray[mafColNumber]
					currentMAFFloat = float(currentMAF)
					currentValues.append(currentMAFFloat)

				# Set is used for the edge map, to find overlaps of associations between phenotypes
				if(phenotype not in phenotypeSnpMap):
					phenotypeSnpMap[phenotype] = []

				# Add the values to our list
				phenotypeSnpMap[phenotype].append(currentValues)

	return phenotypeSnpMap



print("Identifying phenotype to SNP mappings...")
# Initialize SNP Mappings
phenotypeSnpMap = {}
phenotypeSnpMap_sorted = {}

filesToBeParsed = glob.glob(inputData)
if(len(filesToBeParsed) == 0):
	sys.exit('ERROR: No input files have been provided. Check the accuracy of file names for the --input-files argument')

for name in filesToBeParsed:
	with open(name, 'r') as f:
		firstLine = next(f)
		# Figure out where each column is in the input data
		phenotypeColNumber, mafColNumber, pValueColNumber, snpColNumber = getColumnNumbers(name, firstLine)
		phenotypeSnpMap = updatePhenotypeSNPMapFromData(name, f, phenotypeSnpMap, phenotypeColNumber, mafColNumber, pValueColNumber, fileDelimiter)

print("Finished identifying phenotype to SNP mappings.")


with open(datafileOutputPath, 'w') as outfile:
	outfile.write("Phenotype\tAssociatedSNPs\n")
	for key, value in phenotypeSnpMap.items():
		outfile.write(str(key) + "\t" + str(value) + "\n")

print("Finished outputting file. All done!")
