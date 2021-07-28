# NETMAGE: a Human-Disesae Phenotype Map Generator for the Visualization of PheWAS Results

Given genetic associations from a Phenome-Wide Association Study (PheWAS), a disease-disease network can be constructed where nodes represent phenotypes and edges represent shared genetic associations between phenotypes. To improve the ease of network-based analysis of shared genetic components across phenotypes, we developed the humaN-disEase phenoType MAp GEnerator (NETMAGE), a web-based tool that produces interactive phenotype network visualizations from PheWAS summary statistics.

Our service runs at https://hdpm.biomedinfolab.com/netmage. Source code can be downloaded at https://github.com/dokyoonkimlab/netmage.

Contact: dokyoon.kim@pennmedicine.upenn.edu

## Constructing the network
The NETMAGE website will automatically process uploaded PheWAS data and create a corresponding DDN. Input data can be in the format of one file per phenotype, or a single file including all phenotypes and a phenotype ID column. Each row in each file must correspond to a different SNP and have a column corresponding to SNP ID. Minor allele frequency and p-value can also be included as columns to be used for data filtration in NETMAGE. Please refer to the video at https://hdpm.biomedinfolab.com/ for a detailed demonstration of how to use our tool.

Pre-processing on input data is contained within two Python scripts: makeIntermediateNodeMap.py and createNodeAndEdgeMap.py. Both of these files are located in the "scripts" directory. They are built into the NETMAGE website and never need to be used on their own. However, if a user wishes to work directly with either script, both can be downloaded from this GitHub page and run as described later in this README. makeIntermediateNodeMap.py is intended to create a disease-to-snp mapping that can be accessed for network generation. createNodeAndEdgeMap.py will take the disease-to-snp mapping output, apply p-value and minor allele frequency filters, and generate a corresponding node map (includes attributes of each phenotype in the network) and edge map (includes a full list of all connections between phenotypes in the network).


### a) makeIntermediateNodeMap.py
This script will take input PheWAS data and create a corresponding phenotype-to-SNP dictionary. Each key of the dictionary corresponds to a unique phenotype, and each value corresponds to a list of lists. Each sub-list includes the name of the SNP, as well as its p-value and minor allele frequency if provided in the input data. The output of this script can be directly plugged into NETMAGE instead of uploading input data.

#### List of arguments for makeIntermediateNodeMap.py:
```console
{--input-files} Regular expression for input data file(s). Ex: "dir/*.csv"
{--datafile-output} filepath for data output file (tsv)
{--phenotype-name} name of column in input data corresponding to phenotype ID. Optional.
{--snp-name} name of column in input data corresponding to SNP ID. Required.
{--maf-name} name of column in input data corresponding to minor allele frequency. Optional.
{--pvalue-name} name of column in input data corresponding to p-value of association.
{--delim} delimiter used in input PheWAS datasets. Required.
```

#### Example of how to call makeIntermediateNodeMap.py:
```console
python makeIntermediateNodeMap.py --input-files "./data/*.txt" --datafile-output "./disease_snpmap.netmage" --phenotype-name "phenotype" --snp-name "ID" --maf-name "af" --pvalue-name "pval" --delim " "
```


### b) createNodeAndEdgeMap.py
Using the output of makeIntermediateNodeMap.py, this script will apply MAF and p-value filters and create a corresponding node and edge map for the network. These files can be uploaded into the Gephi application for the construction of DDNs, and they are used within NETMAGE as well for DDN generation. If a disease category input file is included, disease-category mappings will be included in the resulting node map and will lead to coloring of nodes in the final DDN according to disease category. Finally, the inclusion of an input ld-file will allow for LD-pruning of SNPs from the DDN.

#### List of arguments for createNodeAndEdgeMap.py:
```console
{--phenotypesnpmap-input} filepath for phenotype-snp map input file. Same as output of makeIntermediateNodeMap.py (tsv)
{--diseasecategory-input} filepath for disease category input file (tsv). Optional.
{--edgemap-output} filepath for edge map output (csv)
{--nodemap-output} filepath for node map output (csv)
{--maf-threshold} threshold for minor allele frequency filtration. Optional.
{--pvalue-threshold} threshold for significance of association between SNP and phenotype. Optional.
{--ld-file} Linkage disequilibrium file generated using plink --show-tags. Custom file will also work as long as the header has "SNP" and "TAGS". File should be space seperated with TAGS seperated by "|". Optional.
```

#### Example of how to call createNodeAndEdgeMap.py:
2) Filtration and inclusion of an input node map, performed on a set of input files:
```console
python createNodeAndEdgeMap.py --phenotypesnpmap-input "./intermediateNodeMap.csv" --diseasecategory-input "./diseaseCategoryMappings.csv" --edgemap-output "./edgeMap.csv" --nodemap-output "./nodeMap.csv" --maf-threshold "0.05" --pvalue-threshold "1e-5" --ld-file "./ldFile.txt"
```




## Using the network
NETMAGE is built off of the Oxford Internet Institute's InteractiveVis framework (https://github.com/oxfordinternetinstitute/InteractiveVis). Users can search maps generated through NETMAGE by a variety of attributes, and they can select nodes to view information such as related phenotypes, associated SNPs, and other network statistics. By examining the associations between phenotypes in generated map, users can gain a better understanding of the underlying genetic architecture for a variety of diseases.

#### How to interact with the network
•	Use the drop-down menu to select a variable by which to search the network map

•	If the variable selected is continuous, a 2-input slider will appear, allowing users to specify a range of values they want to be permitted. The values allowed in this slider will vary depending on the range of values for the chosen variable in the dataset

•	If the variable selected is categorical, a text-box input will appear, allowing users to search for nodes that match the search term. 

•	The Group Selector can be used to display nodes within each category. Nodes are colored by the variable chosen for category. Any node within the graph can be selected – hovering over a node will show the node’s name

•	Selection of a node will display the chosen node and all first-degree neighbors.

•	The right side of the visualization will then display a variety of network statistics, such as Modularity Class, Eccentricity, and Degree, as well as other variables included by the user (i.e associated SNPs)

#### Saving and regenerating the network
•	At the top right-hand corner of any DDN page, two options exist to download a file. First is "pdf download," which will save a screenshot of the current network being depicted. Second is "netmage download," which will save a zip file including the node and edge map (nodeMap.csv and edgeMap.csv) for the network, the Gephi-compatible two-dimensional mapping of the network (data.json), and the preprocessed NETMAGE file (disease_snpmap.netmage). 

•	nodeMap.csv and edgeMap.csv can both be uploaded in the Gephi software in order to create the matching DDN. data.json can also be uploaded to create the DDN, or it can be hosted directly by users on any web server. 

•	If the user wishes to recreate their DDN without having to re-upload all of their PheWAS data, the preprocessed NETMAGE file can be uploaded instead using the "Upload netmage file in place of summary file" button on the NETMAGE homepage. 


