# Fix 20250312

## 1. Gene name char is 20 to 50

## 2. add generator genes of seed

# Time-ordered Gene Coexpression Network (TO-GCN)

Pipeline of time-ordered gene coexpression network (TO-GCN) construction from three-dimensional (gene expression, condition, and time) data

The pipeline contains three steps: (1) Determining the cutoff values, (2) constructing eight GCNs for different coexpression types, and (3) determining time-ordered levels for nodes in a GCN of interest.

## Prepare the gene expression data

Before going to the pipeline, we need to prepare two lists of RPKM values (TF genes and all genes) at different sample points under two conditions (the three-dimensional data). In addition to the data files, you also need to prepare the information of the numbers of samples under condition 1 (n1) and condition (n2).

In the example folder (example_data), there are two data files from the study of "A Comparative Transcriptomics Method to Infer Time-ordered Gene Coexpression Networks and its Applications". The data file should be a Tab-separated values (.tsv) format that contains m rows and n+1 columns, where m is the number of genes (TF genes or all genes) and n represents the summation (n1 + n2) of sample numbers under condition 1 (n1) and condition 2 (n2).The gene IDs are listed in the first column. For each gene, the RPKM values of each sample point under condition 1 and condition 2 are listed from the second to (n1+1)-th columns and from (n1+2)-th to (n1+n2+1)-th columns, respectively. In the example data of TFs_1718.tsv, there are 1718 rows for 1718 TF genes and 27 columns for one gene ID, 13 samples of condition 1, and 13 samples of condition 2.

## Run the programs of pipeline

As mentioned above, there are three steps for the pipeline. Therefore, we provided a program for each step: (1) Cutoff, (2) GCN, and (3) TO-GCN. You can directly run the program by downloading the corresponding binary codes for different system platforms, Linux, MacOSX, or Windows. You can also download the C++ source code (.cpp) and compile to an executable one by yourself. For compiling source codes by yourself, you can use the following commands:

```sh
g++ Cutoff.cpp -o Cutoff
g++ GCN.cpp -o GCN
g++ TO-GCN.cpp -o TO-GCN
```

### (1) Cutoff: Determining the cutoff values

First of all, you need postive and negative cutoff values of Pearson’s Correlation Coefficients (PCCs) under two conditions for constructing the GCN. Our method is to calculate all the PCC values for each TF-gene pair under each condition. With all the PCC values, we generate distributions of probability density function (PDF) and cumulative density function (CDF). According to the CDF, we can suggest you the positive and negative cutoff values with p < 0.05 for each condition. To run the Cutoff program, you have to give 4 parameters: number of samples under condition 1, number of samples under condition 2, data file of TF genes, and data file of all genes. Here is the example of our study:

```sh
Cutoff 13 13 example_data/TFs_1718.tsv example_data/All_genes_25489.tsv
```

In addition to the suggested cutoff values, the program will also generate a file of PCC value distribution in the .tsv format. You can use the file to generate a histogram bar chart by Microsoft Excel or R program.

### (2) GCN: Constructing eight GCNs for different types

In the second step, we want to construct eight coexpression types of GCN under two conditions (C1 and C2): C1+C2+, C1+C20, C1+C2–, C10C2+ C1–C2+, C1–C2–, C1–C20, and C10C2–, where +, -, 0 represents the positive, negative, and no coexpression, respectively. The output file of each GCN is listed in comma-separated value (.csv) format. The five columns represent the TF gene ID, coexpression type, gene ID, PCC under condition 1, PCC under condition 2. You can import these gene pair into the network generation tool, like [Cytoscape](http://www.cytoscape.org), to get the visualization of the GCN. To run the GCN program, you have to give 4 more parameters (total 8 parameters) that indicate the positive cutoff values for conditions 1 and 2 and the negative cutoff values for condtions 1 and 2. Here is an example:

```sh
GCN 13 13 example_data/TFs_1718.tsv example_data/TFs_1718.tsv 0.84 0.84 -0.75 -0.75
```

### (3) TO-GCN: Determining time-ordered level in the interesting GCN

The final step is to determine the time-order (level) of nodes in the GCN. The time-order is assigned by the breadth-first search (BFS) algorithm, starting with a set of seed nodes you chose (listed in seeds.txt). In most case, we will select some genes as seeds that highly expressed in the first time point and lowly expressed in the following time points. In our study, we select a gene with ID, Zm00001d041056, and run the TO-GCN program to assign the time-order (level) of nodes in C1+C2+ GCN. Therefore, we only need the positive cutoffs for conditions 1 and 2 and another 2 parameters that indicate the seed node gene ID and the coexpression type (0, 1, or 2) where 0, 1, and 2 represent the C1+C2+, C1+C20, and C10C2+, respectively. The level of each node and the GCN (obtaine from 2nd step) can be both imported into the [Cytoscape](http://www.cytoscape.org).

```sh
TO-GCN 13 13 example_data/TFs_1718.tsv example_data/TFs_1718.tsv 0.84 0.84 example_data/seeds.txt 0
```

### (4) seed_generator: generate seed list from all gene

```sh
g++ seed_generator.cpp -o seed_generator

./Linux/seed_generator -h

# 基本用法（严格递减模式）
# Basic Usage (Strictly Decreasing Mode)
./Linux/seed_generator -i example_data/All_genes_25489.tsv -s 0 -e 12 -m 15.0 -o high_confidence_seeds1.txt

# examples
# 趋势递减模式（允许波动）
# Trend Decreasing Mode (allows fluctuations)
./Linux/seed_generator -i example_data/All_genes_25489.tsv -s 0 -e 12 -t -l -0.1 -o high_confidence_seeds1.txt
./Linux/seed_generator -i example_data/All_genes_25489.tsv -s 13 -e 25 -t -l -0.1 -o high_confidence_seeds2.txt

# 查找在两种条件都递减的种子基因
# Find seed genes that are decreasing under both conditions
cat high_confidence_seeds1.txt high_confidence_seeds2.txt | sort |uniq -d > target_seeds.txt
```
