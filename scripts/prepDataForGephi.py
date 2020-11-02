# Author: Vivek Sriram
# Last revised: 11/02/2020
# ------------------------------------------------
# Function: Script to prepare input data for Gephi.
#     The script will filter input data, create
#     maps of phenotypes to SNPs, generate a node
#     map, and generate an edge map.

# Import statements
import os
import math
import csv
import operator
import argparse
import sys

# Read in user arguments
parser = argparse.ArgumentParser(description='Take in user input for creation of node and edge maps for Gephi.')
parser.add_argument('--input-directory', dest= 'inputdirectory', default = None, metavar='path to input data directory', help='directory of input PheWAS data')
parser.add_argument('--input-file', dest= 'inputfile', default = None, metavar='path to input data file', help='filepath for input PheWAS data (csv)')
parser.add_argument('--edgemap-output', dest='edgemapout', default = "", metavar='path to edge map output file (csv)', help='filepath for edge map output (csv)')
parser.add_argument('--nodemap-input', dest='nodemapin', default = None, metavar='path to node map input file (csv)', help='filepath for node map input (csv)')
parser.add_argument('--nodemap-output', dest='nodemapout', default = "", metavar='path to node mapp output file (csv)', help='filepath for node map output (csv)')
parser.add_argument('--maf-threshold', dest='mafthresh', metavar = 'threshold for MAF filtration', help='threshold for minor allele frequency filtration')
parser.add_argument('--casecount-threshold', dest='ccthresh', metavar = 'threshold for case count filtration', help='threshold for case count filtration')
parser.add_argument('--pvalue-threshold', dest='pvalthresh', metavar = 'threshold for p-value', help='threshold for significance of association between SNP and phenotype')
parser.add_argument('--phenotype-name', dest='phenotypename', metavar = 'name of variable for phenotype ID', help='name of variable corresponding to phenotype ID')
parser.add_argument('--snp-name', dest='snpname', required=True, metavar = 'name of variable for SNP ID', help='name of variable corresponding to SNP ID')
parser.add_argument('--maf-name', dest='mafname', metavar = 'MAF variable name', help='name of variable corresponding to minor allele frequency')
parser.add_argument('--casecount-name', dest='ccname', metavar = 'case count variable name', help='name of variable corresponding to case count')
parser.add_argument('--pvalue-name', dest='pvname', metavar = 'name of variable for p-value association', help='name of variable corresponding to p-value of association')
parser.add_argument('--delim', dest='delim', required=True, metavar = 'delimiter used in input PheWAS datasets', help='delimiter used in input PheWAS data')

args = parser.parse_args()

# Example usage:
# python prepDataForGephi.py --input-data "./subset/" 
# --edgemap-output "./subsetEdgeMap.csv" 
# --nodemap-input "./phecodeDefinitions.csv" 
# --nodemap-output "./subsetNodeMap.csv" 
# --maf-threshold "0.01" 
# --casecount-threshold "200" 
# --pvalue-threshold "1e-8" 
# --phenotype-name "phenotypeID" 
# --snp-name "ID" 
# --maf-name "af" 
# --casecount-name "num_cases" 
# --pvalue-name "pval" 
# --delim " "

directoryPath = args.inputdirectory
inputFilePath = args.inputfile

if(directoryPath == None and inputFilePath == None):
	sys.exit('ERROR: User needs to supply the filepath to the input data file, or a directory containing the input data')

edgeMapOutputPath = args.edgemapout

nodeMapInputPath = args.nodemapin
nodeMapOutputPath = args.nodemapout

mafThreshold = args.mafthresh
caseCountThreshold = args.ccthresh
pValueThreshold = args.pvalthresh

phenotypeColName = args.phenotypename
snpColName = args.snpname
mafColName = args.mafname
caseCountColName = args.ccname
pValueColName = args.pvname

fileDelimiter = args.delim


# Figure out where each column is in the input data
phenotypeColNumber = None
mafColNumber = None
caseCountColNumber = None
pValueColNumber = None
snpColNumber = None

if(directoryPath != None):
	# Get the header row of the first file in our input directory
	firstFileName = os.listdir(directoryPath)[0]
	inputFilePath = os.path.join(directoryPath, firstFileName)

firstLine = None
with open(inputFilePath, 'r') as firstFile:
	firstLine = firstFile.readline()
	firstLine = firstLine.strip()
	firstLineAsArray = firstLine.split(fileDelimiter)
	# Figure out the column indices corresponding to phenotype, MAF, case count, p-value, and SNP (if they exist)
	for i in range(0, len(firstLineAsArray)):
		if(firstLineAsArray[i] == phenotypeColName):
			phenotypeColNumber = i
		elif(firstLineAsArray[i] == str(mafColName)):
			mafColNumber = i
		elif(firstLineAsArray[i] == str(caseCountColName)):
			caseCountColNumber = i
		elif(firstLineAsArray[i] == str(pValueColName)):
			pValueColNumber = i
		elif(firstLineAsArray[i] == snpColName):
			snpColNumber = i

# We need a column corresponding to snp name in our data. Force an exit if it can't be found
if(snpColNumber == None):
	sys.exit('ERROR: Cannot identify column corresponding to SNP name in input data')

if(nodeMapInputPath != None and phenotypeColNumber == None):
	sys.exit('ERROR: You have provided an input node map without specifying the phenotype ID column in your data. Either avoid including an input node map or give the name of the column for phenotype.')

# Given a line of data, indices for different variables, and threshold values, check if the line passes these thresholds
def linePassesFilter(line, mafColNumber, mafFilter, ccColNumber, ccFilter, pvalColNumber, pvalFilter):
	# If there is a line corresponding to case count, check if it passes input filter (should be at least the threshold)
	if(ccColNumber != None and ccFilter != None):
		if(float(line[ccColNumber]) < float(ccFilter)):
			return False
	# If there is a line corresponding to minor allele frequency, check if it passes input filter (should be at least the threshold)
	if(mafColNumber != None and mafFilter != None):
		if(float(line[mafColNumber]) < float(mafFilter)):
			return False
	# If there is a line corresponding to p-value, check if it passes input filter (should be at most the threshold)
	if(pvalColNumber != None and pvalFilter != None):
		if(float(line[pvalColNumber]) > float(pvalFilter)):
			return False
	return True


print("Identifying phenotype to SNP mappings...")
# Initialize SNP Mappings
phenotypeSnpMap_unsorted = {}
phenotypeSnpMap_sorted = {}
phenotypeSnpMap_snpAndPVal= {}

def updatePhenotypeSNPMapFromData(fileName, fileContent, phenotypeSnpMap_unsorted, phenotypeSnpMap_snpAndPVal, phenotypeColNumber, mafColNumber, mafThreshold, caseCountColNumber, caseCountThreshold, pValueColNumber, pValueThreshold, fileDelimiter):
	# Save out the file name as the phenotype for now
	phenotype = os.path.splitext(fileName)[0]

	# Get all lines in the file, except for the header line
	next(fileContent)
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
			
			# If the line passes all specified thresholds, we can use it
			if(linePassesFilter(lineAsArray, mafColNumber, mafThreshold, caseCountColNumber, caseCountThreshold, pValueColNumber, pValueThreshold)):
				# Get the name of the SNP
				snpName = lineAsArray[snpColNumber]

				# Set is used for the edge map, to find overlaps of associations between phenotypes
				currentSetOfSNPsForThisPhenotype = set()
				# Get the current set of SNPs corresponding to this phenotype, if one exists already
				if(phenotype in phenotypeSnpMap_unsorted.keys()):
					currentSetOfSNPsForThisPhenotype = phenotypeSnpMap_unsorted[phenotype]
				# Add the SNP to which this line corresponds to our set of associated SNPs for this phenotype
				currentSetOfSNPsForThisPhenotype.add(snpName)
				phenotypeSnpMap_unsorted[phenotype] = currentSetOfSNPsForThisPhenotype

				# If the input data includes p-value information, get the p-value for the SNP too
				if(pValueColNumber != None):
					snpPValue = lineAsArray[pValueColNumber]
					snpPValueFloat = float(snpPValue)
					
					# List is used for the node map to present a sorted list of nodes
					currentListOfSNPsForThisPhenotype = []
					if(phenotype in phenotypeSnpMap_snpAndPVal.keys()):
						currentListOfSNPsForThisPhenotype = phenotypeSnpMap_snpAndPVal[phenotype]
					# Add a tuple of (SNP ID, p-value) to our list
					currentListOfSNPsForThisPhenotype.append((snpName, snpPValueFloat))
					phenotypeSnpMap_snpAndPVal[phenotype] = currentListOfSNPsForThisPhenotype

	return (phenotypeSnpMap_sorted, phenotypeSnpMap_snpAndPVal)


if(directoryPath != None):
	# Create the SNP mapping from our files
	for filename in os.listdir(directoryPath):
		with open(os.path.join(directoryPath, filename), 'r') as f:
			phenotypeSnpMap_sorted, phenotypeSnpMap_snpAndPVal = updatePhenotypeSNPMapFromData(filename, f, phenotypeSnpMap_unsorted, phenotypeSnpMap_snpAndPVal, phenotypeColNumber, mafColNumber, mafThreshold, caseCountColNumber, caseCountThreshold, pValueColNumber, pValueThreshold, fileDelimiter)
else:
	with open(inputFilePath, 'r') as f:
		phenotypeSnpMap_sorted, phenotypeSnpMap_snpAndPVal = updatePhenotypeSNPMapFromData(inputFilePath, f, phenotypeSnpMap_unsorted, phenotypeSnpMap_snpAndPVal, phenotypeColNumber, mafColNumber, mafThreshold, caseCountColNumber, caseCountThreshold, pValueColNumber, pValueThreshold, fileDelimiter)

# Sort all the tuple lists by phenotype at the end
for key in phenotypeSnpMap_snpAndPVal:
	snpsForThisPhenotype = phenotypeSnpMap_snpAndPVal[key]
	# Sort the list by p-value (i.e., by the second element of each tuple)
	snpsForThisPhenotype.sort(key = lambda x: float(x[1]), reverse = False)
	# Get a list of the now sorted SNP IDs
	snpNamesSorted = [i[0] for i in snpsForThisPhenotype]
	# Update the list in our map
	phenotypeSnpMap_sorted[key] = snpNamesSorted

print("Finished identifying phenotype to SNP mappings.")



# If we want to create a node map file
if(nodeMapOutputPath != ""):
	print("Creating node map...")
	print("We will be processing " + str(len(phenotypeSnpMap_unsorted)) + " phenotypes in total for our node map.")

	with open(nodeMapOutputPath, 'a') as outf:
		nodeMapOutput = csv.writer(outf, delimiter = '\t')

		# If the user has not provided an input map off of which to work, we will come up with our own
		if(nodeMapInputPath == None):
			print("Initializing a new node map")
			header = ['phenotype'] + ['associatedSNPs']
			nodeMapOutput.writerow(header)

			# If p-value is involved, iterate over key/value pairs in List
			if(len(phenotypeSnpMap_sorted) > 0):
				for key1, value1 in phenotypeSnpMap_sorted.items():
					row = [str(key1)] + [str(list(value1))]
					nodeMapOutput.writerow(row)
					print("Processed phenotype " + key1 + " for node map")
			# Otherwise, iterate over key/value pairs in set
			else:
				for key2, value2 in phenotypeSnpMap_unsorted.items():
					row = [str(key2)] + [str(list(value2))]
					nodeMapOutput.writerow(row)
					print("Processed phenotype " + key2 + " for node map")

		# Otherwise, if the user has provided an input node map, where the first column is the phenotype,
		# then we can just append associated SNPs as a final column
		else:
			print("Appending to input node map")
			# Read in the node map input file
			with open(nodeMapInputPath, 'r') as inf:
				currentNodeMap = csv.reader(inf)

				# Add the header to our output node map
				header = next(currentNodeMap)
				header.append('associatedSNPs')
				nodeMapOutput.writerow(header)

				# Iterate through the lines of the current node map
				for row in currentNodeMap:
					# Get the phenotype of the current row
					currentPhenotype = row[0]

					# If there are elements in the phenotypeSnpMap_sorted, check there
					if(len(phenotypeSnpMap_sorted) > 0):
						if(currentPhenotype in phenotypeSnpMap_sorted.keys()):
							snpsForCurrPhenotype = str(list(phenotypeSnpMap_sorted[currentPhenotype]))
							row.append(snpsForCurrPhenotype)
							nodeMapOutput.writerow(row)
					# Otherwise, refer to the phenotypeSnpMap_unsorted
					else:
						if(currentPhenotype in phenotypeSnpMap_unsorted.keys()):
							snpsForCurrPhenotype = str(list(phenotypeSnpMap_unsorted[currentPhenotype]))
							row.append(snpsForCurrPhenotype)
							nodeMapOutput.writerow(row)

					print("Processed phenotype " + currentPhenotype + " for node map")

	print("Finished creating node map.")


# If we want to create the edge map
if(edgeMapOutputPath != ""):
	print("Creating edge map...")

	# Write to an output file that will include all pairs of phenotypes and their shared SNPs
	with open(edgeMapOutputPath, "a") as edgeMapOutput:
		edgeMapOutput.write("Source\tTarget\tWeight\tlistOfSharedSNPs\n")
		snpPairsAlreadySeen = set()

		# Do a nested for loop that will get the intersection of our lists for each pair of nodes
		i = 0
		for key1, value1 in phenotypeSnpMap_unsorted.items():
			i = i+1
			for key2, value2 in phenotypeSnpMap_unsorted.items():
				# If we've already seen 
				snpPairV1 = key1+"_"+key2
				snpPairV2 = key2+"_"+key1

				if((key1 != key2) and (snpPairV1 not in snpPairsAlreadySeen) and (snpPairV2 not in snpPairsAlreadySeen)):
					snpPairsAlreadySeen.add(snpPairV1)
					snpPairsAlreadySeen.add(snpPairV2)

					# Get the intersection of SNPs for the two phenotypes
					sharedSNPs = list(set.intersection(value1, value2))

					# Write this line to the output file
					if(len(sharedSNPs) > 0):
						outputLine = key1 + "\t" + key2 + "\t" + str(len(sharedSNPs)) + "\t" + str(sharedSNPs) + "\n"
						edgeMapOutput.write(outputLine)

			print("Processed phenotype " + key1 + ", " + str(i) + " out of " + str(len(phenotypeSnpMap_unsorted.keys())) + " phenotypes for the edge map")

	print("Finished creating edge map.")

print("All done!")

