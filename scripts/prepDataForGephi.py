# Author: Vivek Sriram
# Last revised: 10/25/2020
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
parser.add_argument('--inputdata', dest= 'inputdata', required=True, metavar='path to input data directory', help='directory of input PheWAS data')
parser.add_argument('--filteredoutput', dest='filteredout', default="", metavar='path to filtered data directory', help='directory to contain filtered PheWAS data')
parser.add_argument('--edgemapoutput', dest='edgemapout', default = "", metavar='path to edge map output file', help='filepath for edge map output')
parser.add_argument('--nodemapinput', dest='nodemapin', default = "", metavar='path to node map input file', help='filepath for node map input')
parser.add_argument('--nodemapoutput', dest='nodemapout', default = "", metavar='path to node mapp output file', help='filepath for node map output')
parser.add_argument('--mafthreshold', dest='mafthresh', default=0.01, metavar = 'threshold for MAF filtration', help='threshold for minor allele frequency filtration')
parser.add_argument('--casecountthreshold', dest='ccthresh', default=200, metavar = 'threshold for case count filtration', help='threshold for case count filtration')
parser.add_argument('--pvaluethreshold', dest='pvalthresh', default="1e-8", metavar = 'threshold for p-value', help='threshold for significance of association between SNP and phenotype')
parser.add_argument('--phecodename', dest='phecodename', required=True, metavar = 'name of variable for phenotype ID', help='name of variable corresponding to phenotype ID')
parser.add_argument('--snpname', dest='snpname', required=True, metavar = 'name of variable for SNP ID', help='name of variable corresponding to SNP ID')
parser.add_argument('--mafname', dest='mafname', metavar = 'MAF variable name', help='name of variable corresponding to minor allele frequency')
parser.add_argument('--casecountname', dest='ccname', metavar = 'case count variable name', help='name of variable corresponding to case count')
parser.add_argument('--pvaluename', dest='pvname', metavar = 'name of variable for p-value association', help='name of variable corresponding to p-value of association')
parser.add_argument('--delim', dest='delim', required=True, metavar = 'delimiter used in input PheWAS datasets', help='delimiter used in input PheWAS data')

args = parser.parse_args()

# Example usage:
# python prepDataForGephi.py --inputdata "./subset/" 
# --filteredoutput "./filteredSubset/" 
# --edgemapoutput "./subsetEdgeMap.csv" 
# --nodemapinput "./phecodeDefinitions.csv" 
# --nodemapoutput "./subsetNodeMap.csv" 
# --mafthreshold "0.01" 
# --casecountthreshold "200" 
# --pvaluethreshold "1e-8" 
# --phecodename "phenotypeID" 
# --snpname "ID" 
# --mafname "af" 
# --casecountname "num_cases" 
# --pvaluename "pval" 
# --delim " "

directoryPath = args.inputdata

filteredDirectoryPath = args.filteredout

edgeMapOutputPath = args.edgemapout

nodeMapInputPath = args.nodemapin
nodeMapOutputPath = args.nodemapout

mafThreshold = args.mafthresh
caseCountThreshold = args.ccthresh
pValueThreshold = args.pvalthresh

phecodeColName = args.phecodename
snpColName = args.snpname
mafColName = args.mafname
caseCountColName = args.ccname
pValueColName = args.pvname

fileDelimiter = args.delim


# Figure out where each column is in the input data
phecodeColNumber = None
mafColNumber = None
caseCountColNumber = None
pValueColNumber = None
snpColNumber = None

# Get the header row of the first file in our input directory
firstFileName = os.listdir(directoryPath)[0]
firstFile = open(directoryPath + firstFileName, 'r')
firstLine = firstFile.readlines()[0]
firstLineAsArray = firstLine.split(fileDelimiter)

for i in range(0, len(firstLineAsArray)):
    if(firstLineAsArray[i] == phecodeColName or firstLineAsArray[i] == phecodeColName + "\n"):
        phecodeColNumber = i
    elif(firstLineAsArray[i] == str(mafColName) or firstLineAsArray[i] == str(mafColName) + "\n"):
        mafColNumber = i
    elif(firstLineAsArray[i] == str(caseCountColName) or firstLineAsArray[i] == str(caseCountColName) + "\n"):
        caseCountColNumber = i
    elif(firstLineAsArray[i] == str(pValueColName) or firstLineAsArray[i] == str(pValueColName) + "\n"):
        pValueColNumber = i
    elif(firstLineAsArray[i] == snpColName or firstLineAsArray[i] == snpColName + "\n"):
        snpColNumber = i
firstFile.close()

# We need a column corresponding to phenotype name in our data. Force an exit if it can't be found
if(phecodeColNumber == None):
    sys.exit('ERROR: Cannot identify column corresponding to phenotype name in input data')
if(snpColNumber == None):
    sys.exit('ERROR: Cannot identify column corresponding to SNP name in input data')

# If we want to filter our input directory:
if(filteredDirectoryPath != ""):
    print("Filtering input data...")
    for filename in os.listdir(directoryPath):
        # Open a new file in the filter directory that corresponds to the filtered version of the current file
        currentFilterOutput = open(filteredDirectoryPath + "filtered_" + filename, 'a')

        # Read through the lines of the current file
        f = open(directoryPath + filename, 'r')
        lines = f.readlines()
        lineCounter = 0

        # Get the elements of each line (and make sure to remove the EOL character in the last column)
        for elem in lines:
            elemAsArray = elem.split(fileDelimiter)
            elemAsArray[len(elemAsArray)-1] = elemAsArray[len(elemAsArray)-1].rstrip()

            # Initialize the output line as well as a boolean that tells us if the line passes our desired filters
            outputLine = ""
            linePassesFilters = True

            # If there is a column corresponding to phenotype name, add it to the line
            if(phecodeColNumber != None):
                outputLine = outputLine + str(elemAsArray[phecodeColNumber])
            else:
                linePassesFilters = False

            # If we want to filter by case count
            if(caseCountColNumber != None):
                if(lineCounter == 0):
                    outputLine = outputLine + fileDelimiter + str(elemAsArray[caseCountColNumber])
                elif(float(elemAsArray[caseCountColNumber]) >= float(caseCountThreshold)):
                    outputLine = outputLine + fileDelimiter + str(elemAsArray[caseCountColNumber])
                else:
                    linePassesFilters = False

            # If we want to filter by MAF
            if(mafColNumber != None):
                if(lineCounter == 0):
                    outputLine = outputLine + fileDelimiter + str(elemAsArray[mafColNumber])
                elif(float(elemAsArray[mafColNumber]) >= float(mafThreshold)):
                    outputLine = outputLine + fileDelimiter + str(elemAsArray[mafColNumber])
                else:
                    linePassesFilters = False

            # If there is a column corresponding to SNP name, add it to the line
            if(snpColNumber != None):
                outputLine = outputLine + fileDelimiter + str(elemAsArray[snpColNumber])
            else:
                linePassesFilters = False

            # If we want to filter by p-value
            if(pValueColNumber != None):
                if(lineCounter == 0):
                    outputLine = outputLine + fileDelimiter + str(elemAsArray[pValueColNumber])
                elif(float(elemAsArray[pValueColNumber]) <= float(pValueThreshold)):
                    outputLine = outputLine + fileDelimiter + str(elemAsArray[pValueColNumber])
                else:
                    linePassesFilters = False

            # If this is the first line, or if the line passess our input filters, we should add it to our output file
            if(lineCounter == 0 or linePassesFilters == True):
                currentFilterOutput.write(outputLine+ fileDelimiter + "\n")

            lineCounter = lineCounter +1
        f.close()
        currentFilterOutput.close()
    print("Finished filtering data")


print("Identifying phenotype to SNP mappings...")
# Initialize SNP Mappings
snpMappingsList = {}
snpMappingsSet = {}

# If the user didn't want to filter their input data, then we just use the original directory
if(filteredDirectoryPath == ""):
    filteredDirectoryPath = directoryPath

# Create the SNP mapping from our files
for filename in os.listdir(filteredDirectoryPath):
    # Make sure to ignore any DS_Store files that are automatically generated by filtration
    if(filename.split(".")[1] != "DS_Store"):
        f = open(filteredDirectoryPath + filename, 'r')
        lines = f.readlines()
        indexforPhenotype = None
        indexforSnpName = None
        indexforPValue = None

        snpsForThisPhecodeList = []
        snpsForThisPhecodeSet = set()

        # Figure out where phenotype name, snp name, and p-value are in the filtered data that are being processed
        firstLine = lines[0]
        firstLineAsArray = firstLine.split(' ')
        for i in range(0, len(firstLineAsArray)):
            if(firstLineAsArray[i] == phecodeColName):
                indexforPhenotype = i
            elif(firstLineAsArray[i] == snpColName):
                indexforSnpName = i
            elif(firstLineAsArray[i] == pValueColName):
                indexforPValue = i

        for elem in lines:
            # Get the name of the phenotype
            phecode = elem.split(fileDelimiter)[indexforPhenotype]
            # Get the name and p-value of the SNP
            snpName = elem.split(fileDelimiter)[indexforSnpName]
            # If this isn't the first row of the file, add this tuple to our list/set
            if(snpName != snpColName):
                snpsForThisPhecodeSet.add(snpName)

            if(indexforPValue != None):
                snpPValue = elem.split(fileDelimiter)[indexforPValue]

                # If this isn't the first row of the file, add this tuple to our list/set
                if(snpName != snpColName):
                    snpPValueFloat = float(snpPValue)
                    snpsForThisPhecodeList.append((snpName, snpPValueFloat))
        
        # If SNPs have been found for this phecode, we add it to our mappings
        if(len(snpsForThisPhecodeSet) > 0):
            # Set is used for the edge map, to find overlaps of associations between phenotypes;
            snpMappingsSet[phecode] = snpsForThisPhecodeSet

        if(len(snpsForThisPhecodeList) > 0):
        # List is used for the node map to present a sorted list of nodes
            # Sort the list by p-value (i.e., by the second element of each tuple)
            snpsForThisPhecodeList.sort(key = lambda x: float(x[1]), reverse = True)
            snpsForThisPhecodeList.reverse()
            # Get a list of the now sorted SNP IDs
            snpNamesSorted = [i[0] for i in snpsForThisPhecodeList]
            snpMappingsList[phecode] = snpNamesSorted

print("Finished identifying phenotype to SNP mappings.")



# If we want to create a node map file
if(nodeMapInputPath != "" and nodeMapOutputPath != ""):
    print("Creating node map...")

    # Figure out if we need the p-value sorted list or not
    numNodesToProcess = -1
    if(len(snpMappingsList.keys()) > 0):
        numNodesToProcess = len(snpMappingsList.keys())
    else:
        numNodesToProcess = len(snpMappingsSet.keys())

    print("We will be processing " + str(numNodesToProcess) + " phenotypes in total for our node map.")

    # Read in the node map input and output files
    currentNodeMap = csv.reader(open(nodeMapInputPath, 'r'))
    nodeMapOutput = csv.writer(open(nodeMapOutputPath, 'w'))

    # Iterate through the lines of the current node map
    rowCounter = 0
    for row in currentNodeMap:
        # Get the phecode of the current row
        currentPhecode = row[0]

        # If it's the first row, append it
        if(rowCounter == 0):
            row.append('associatedSNPs')
            nodeMapOutput.writerow(row)

        # Otherwise, see if we have the phenotype in our mapping
        else:
            # If there are elements in the snpMappingsList, check there
            if(len(snpMappingsList.keys()) > 0):
                if(currentPhecode in snpMappingsList.keys()):
                    snpsForCurrPhecode = str(snpMappingsList[currentPhecode])
                    row.append(snpsForCurrPhecode)
                    nodeMapOutput.writerow(row)
            # Otherwise, refer to the unsorted snpMappingsSet
            else:
                if(currentPhecode in snpMappingsSet.keys()):
                    snpsForCurrPhecode = str(snpMappingsSet[currentPhecode])
                    row.append(snpsForCurrPhecode)
                    nodeMapOutput.writerow(row)

        print("Processed phenotype " + currentPhecode + " for node map")
        rowCounter = rowCounter + 1

    print("Finished creating node map.")
    

# If we want to create the edge map
if(edgeMapOutputPath != ""):
    print("Creating edge map...")

    # Write to an output file that will include all pairs of phecodes and their shared SNPs
    edgeMapOutput = open(edgeMapOutputPath, "a")
    edgeMapOutput.write("Source\tTarget\tWeight\tlistOfSharedSNPs\n")
    snpPairsAlreadySeen = set()

    # Do a nested for loop that will get the intersection of our lists for each pair of nodes
    i = 0
    for key1, value1 in snpMappingsSet.items():
        i = i+1
        for key2, value2 in snpMappingsSet.items():
            # If we've already seen 
            snpPairV1 = key1+"_"+key2
            snpPairV2 = key2+"_"+key1

            if((key1 != key2) and (snpPairV1 not in snpPairsAlreadySeen) and (snpPairV2 not in snpPairsAlreadySeen)):
                snpPairsAlreadySeen.add(snpPairV1)
                snpPairsAlreadySeen.add(snpPairV2)

                # Get the intersection of SNPs for the two phenotypes
                sharedSNPs = list(set.intersection(value1, value2))

                # Write this line to the output file
                if(len(list(sharedSNPs)) > 0):
                    outputLine = key1 + "\t" + key2 + "\t" + str(len(list(sharedSNPs))) + "\t" + str(list(sharedSNPs)) + "\n"
                    edgeMapOutput.write(outputLine)

        print("Processed phenotype " + key1 + ", " + str(i) + " out of " + str(len(snpMappingsSet.keys())) + " phenotypes for the edge map")

    print("Finished creating edge map.")
    edgeMapOutput.close()


print("All done!")

