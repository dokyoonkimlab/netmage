# NETMAGE: a humaN-disEase phenoType MAp GEnerator for the Visualization of PheWAS

Given genetic associations from a Phenome-Wide Association Study (PheWAS), a disease-disease network can be constructed where nodes represent phenotypes and edges represent shared genetic associations between pheno-types. To improve the accessibility of the visualization of shared genetic components across pheno-types, we developed the humaN-disEase phenoType MAp GEnerator (NETMAGE), a web-based tool that produces interactive phenotype network visualizations from summarized PheWAS results.

Our service runs at https://hdpm.biomedinfolab.com. Source code can be downloaded at https://github.com/dokyoonkimlab/netmage.

Contact: dokyoon.kim@pennmedicine.upenn.edu

## Constructing the network
We provide a script (prepDataForGephi.py) which will take in input PheWAS summary data from the user and generate corresponding "edge map" and "node map" file. Our script also provides users with the ability to filter their input PheWAS data by minor allele frequency, case count, and p-value. Instructions for how to use prepDataForGephi.py are included in the NETMAGE User Guide (https://hdpm.biomedinfolab.com/netmage/userguide/ or NETMAGE_UserGuide.docx in the GitHub repository).

#### List of arguments for prepDataForGephi.py:
```console
{--input-files} Regular expression for input data file(s)
{--edgemap-output} filepath for edge map output file (csv)
{--nodemap-input} filepath for node map input (csv)
{--nodemap-output} filepath for node map output (csv)
{--maf-threshold} threshold for minor allele frequency filtration
{--casecount-threshold} threshold for case count filtration
{--pvalue-threshold} threshold for significance of association between SNP and phenotype
{--phenotype-name} name of variable corresponding to phenotype ID
{--snp-name} name of variable corresponding to SNP ID
{--maf-name} name of variable corresponding to minor allele frequency
{--casecount-name} name of variable corresponding to case count
{--pvalue-name} name of variable corresponding to p-value of association
{--delim} delimiter used in input PheWAS datasets
```

#### Examples of how to call prepDataForGephi.py:

1) No filtration and no input node map, performed on a single input file: 
```console
python prepDataForGephi.py --input-files "./subset/allData.txt" --edgemap-output "./subsetEdgeMap.csv" --nodemap-output "./subsetNodeMap.csv" --snp-name "ID" --delim " "
```

2) Filtration and inclusion of an input node map, performed on a set of input files:
```console
python prepDataForGephi.py --input-files "./subset/*.txt" --edgemap-output "./subsetEdgeMap.csv" --nodemap-input "./phecodeDefinitions.csv" --nodemap-output "./subsetNodeMap.csv" --maf-threshold "0.01" --casecount-threshold "200" --pvalue-threshold "1e-8" --phenotype-name "phenotypeID" --snp-name "ID" --maf-name "af" --casecount-name "num_cases" --pvalue-name "pval" --delim " "
```

Following the creation of the edge and node map files, users can create their desired networks by processing data through Gephi. Instructions for creation of a Gephi network from input can also be found in the User Guide. Finally, the "data.json" file produced as output from Gephi can be plugged into NETMAGE (https://hdpm.biomedinfolab.com/netmage/) to create a corresponding disease-disease network.

## Using the network
NETMAGE is built off of the Oxford Internet Institute's InteractiveVis framework (https://github.com/oxfordinternetinstitute/InteractiveVis). Users can search maps generated through NETMAGE by a variety of attributes, and they can select nodes to view information such as related phenotypes, associated SNPs, and other network statistics. By examining the associations between phenotypes in generated map, users can gain a better understanding of the underlying genetic architecture for a variety of diseases.

#### How to interact with the network
•	Use the drop-down menu to select a variable by which to search the network map

•	If the variable selected is continuous, a 2-input slider will appear, allowing users to specify a range of values they want to be permitted. The values allowed in this slider will vary depending on the range of values for the chosen variable in the dataset

•	If the variable selected is categorical, a text-box input will appear, allowing users to search for nodes that match the search term. 

•	The Group Selector can be used to display nodes within each category. Nodes are colored by the variable chosen for category. Any node within the graph can be selected – hovering over a node will show the node’s name

•	Selection of a node will display the chosen node and all first-degree neighbors.

•	The right side of the visualization will then display a variety of network statistics, such as Modularity Class, Eccentricity, and Degree, as well as other variables included by the user (i.e associated SNPs)
