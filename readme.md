# Gene Co-expression Network Construction and Analysis

This is a Python package for the generation and analysis of Gene Co-expression Networks. The source code is in the `src` folder and the executable file `GeCoNet_Tool.exe` is in: https://github.com/KSUNetSE/GeCoNet_Tool/releases or click **`Releases`** on the right side.

## Network construction

`Input data (.csv):`
select the comma delimited gene co-expression dataset, where the first row is the experimental conditions and the first column is the name of genes. For example:
|  | condition_1| condition_2|condition_3|
|---------|--------|------|-----|
|gene_1 | value | value|value|
|gene_2 | value|value|value|

default: anopheles.csv

`Remove zeros:`
check the box will remove zeros from the input dataset, this is mainly used to process RNA sequencing data.

`rescale values by log2:`
rescale each value by log2, this is mainly used to process RNA sequencing data.

`z-score columns:`
Each column will be z-score normalized.

`Drop columns (size<):`
Drop the columns with less than a certain number of values.
default: 200

`save processed data:`
The processed data will be saved as: 'input (z-scored).csv', where input is the name of input data. 

`save PCC matrix:`
The PCC between every pair of nodes will be saved as an upper triangular matrix.
default: 'input PCC matrix.csv'

`save paired element matrix:`
The number of paired elements between every pair of nodes will be saved as an upper triangular matrix.
default: 'input paired element matrix.csv'

`Optimize the sliding threshold parameters:`
Optimize the four parameters of the sliding threshold curve if the box is checked. The coefficient of determination (R squared) will be printed in the **Running status** box.
If the box is not checked, the user can input the four parameters .

`Alpha, Eta, Lambda, Eta:`
The initial parameters used to opimize the sliding  threshold curve.

`Bin size:`
 binning the PCCs into different intervals per paired element.
default: 10, for example:
|4, 5, ... , 10| 11, 12, ... , 20 | 21, 22, ... , 30 | 
|----------------|-------------------------------|-----------------------------|

In this case, the node pairs with 4, 5, ..., 10 paired elements will be categorized into the same class.

`Cutoff`
 choose the top fraction of gene pairs with the highest PCCs as the edges of the network.
default: 0.005

`save edge list:`
The edge list of the constructed network will be saved as comma delimited file: 'input edgelist.csv'. The edge list is saved in the following format:
| source | target | weight|
|---------|--------|------|
|node_1 | node_2 | edge_weight|

`save threshold curve:`
save and show the figure of the fitted sliding threshold curve.

## Network analysis

`Input edge list(.csv):`
select a comma delimited edge list in the following format:
| source | target | weight|
|---------|--------|------|
|node_1 | node_2 | edge_weight|

The third column can be neglected if the network is unweighted.
default: anopheles.csv

`Save & Analyze largest component:`
The program will save and analyze the largest connected component. The largest conneccted component will be saved as: "edgelist (LCC).csv", where edgelist represents the input name of the input network.

`weighted:`
Use weighted edges to produce network properties.

` Louvain community:`
Use the louvain algorithm to detect the communities of the network. 

` Louvain community:`
Use the Leiden algorithm to detect the communities of the network. 
 Interested Users are recommended to edit the raw code to customize the settings for communities detection.
 
`core:`
Detect the core of the network. 
The program will visualize the network with the Fruchterman Reingold algorihtm if `community` or (and) `core` is selected.

`degree:`
Calculate the degree of the nodes and produce the degree distribution, which will also be saved as a figure.

`eigenvector:`
Calculate the eigenvector centralities of the nodes.

`closeness:`
Calculate the closeness of the nodes.

`betweenness:`
Calculate the betweenness of the nodes.

After running the button, the properties table will be saved as 'edgelist property list.csv', where edgelist is the name of the input network.

## Running status
The running status of the program will be shown in the box.

# Results
## Statistics
The number of nodes, edges will be printed in this section. The number of core nodes will be printed if the core is selected in the part of network analysis.
 
## Figures
Users can visualize the figures produced in network construction and network analysis.

