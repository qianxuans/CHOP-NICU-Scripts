# Customized scripts of CHOP NICU analysis
## Hierarchical Clustering using Silhouette Scores

The customized script ```Hierarchical_Clustering.R``` uses the core-genome phylogenetic tree to partition all genomes into different groups. All possible grouping pattern, which range from 2 groups to (tip number - 1) groups.  The grouping pattern with the highest average silhouette scores is the final group. 

These groups are selected for downstream analysis based on two criteria: they contained at least 2 genomes, and the minimal SNP distance among genomes within the group was below 500. The genomes in the remaining 5 groups, which did not meet the criteria for fine-grained analysis, were considered singletons.

### Flags
```
Rscript Hierarchical_Clustering.R --tree /path/to/core-genome/phylogenetic/tree --matrix /path/to/snp/distance/matrix
```

- --tree: This flag takes the newick format of the core-genome phylogenetic tree.
  
In our analysis, our tree is namaed "NICU_core_genome.contree"

- --matrix: For now, this flag takes the RDS file that saves a single R object containing the SNP matrix.
  
In our analysis, the RDS file is "NICU_snp_matrix.RDS"
### Ouput
The output is the csv file with each row indicating the genome name and the group ID. 
- Output ```all_hierarchical_clustering_groups.csv```: This shows all groups including singletons 
- Output ```selected_hierarchical_clustering_groups.csv```: This shows only groups that meet the criteria described previsouly.
- Output ```singletons_hierarchical_clustering_groups.csv```: This shows only the singletons
### Reference
- [TreeTools](https://ms609.github.io/TreeTools/) 
- [ape](https://emmanuelparadis.github.io/) 
- [dplyr](https://github.com/tidyverse/dplyr) 
- [cluster](https://svn.r-project.org/R-packages/trunk/cluster/)

  
## Strain Determination
The script ```Strain_Determination.R``` takes the ```all_hierarchical_clustering_groups.csv``` produced by by ```Hierarchical_Clustering.R```, and only analyzed groups that meet the criteria described previsouly.
A single linkage clustering (SLC) algorithm is then used to determine closely related strains at every possible potential SNP threshold ranging from the minimal SNP distance in the group to the max threshold user defines. The output from SLC is then corrected for phylogenetic structure. The SNP threshold at which the number and composition of strains plateaued is selected to determine strains. If no plateau is found, the max threshold defined previously would be used to determine the strains for this group.

### Flags
- --subtrees: The path to directoy containing reference-based whole-genome phylogenetic trees of each selected groups.(For now the script only takes the trees built by IQtree with .contree extension, and the name of the tree should be "Group" + Group ID, like "Group1.contree")
- --matrix: The RDS file containing the SNP matrix. In our analysis, the RDS file is "NICU_snp_matrix.RDS".
- --cores: How many cores would be used for this analysis. Default is all cores available.
- --max_threshold: The max SNP threshold that the analysis used to determine strains.
- --plateau_length: The length of plateau for determining strains.
- --groups: ```all_hierarchical_clustering_groups.csv``` produced by ```Hierarchical_Clustering.R```.
- --output: The path to the output directory. 
### Output 
- ```Strains_summary.csv``` is the csv file that summarizes the determination of the strains.
  * The column "strain" indicates the strain determination before assigning new strain IDs. The number before "_" indicates the group, and the number after "_" indicate the strain ID inside the group.
    For example, ```1_2``` indicates the Strain2 in Group 1. ```1_Singleton``` indicates the singleton from Group1. ```HC_Singleton``` represents the singleton inferred using ```Hierarchical_Clustering.R```
  * The column "isolate" indicates the name of the genome.
  * Column "new_strain_id" indicates the new strain ID after the analysis
- ```All_Groups_Plateau.csv``` shows the number of genomes and the plateau threshold found in each of the groups.
  * If "No Plateau Found" is shown, this means that using the given max threshold and plateau length, no plateau is found, and the max threshold is used as the threshold to determine the strains.
- ```<GroupID>_clones_number.csv``` shows the number of clones and singletons before and after the correction by phylogenetic structure, and how many clones are broken by the public available genomes at each SNP threshold used in the analysis.
- ```<GroupID>_stats.csv``` shows the number of discrepancy genomes at each SNP threshold used in the analysis.
- ```<GroupID>_discrepancy clones.pdf``` visualizes the number of discrepancy genomes and clones at each SNP threshold used in the analysis.


### Reference
- [TreeTools](https://ms609.github.io/TreeTools/) 
- [ape](https://emmanuelparadis.github.io/) 
- [dplyr](https://github.com/tidyverse/dplyr) 
- [cluster](https://svn.r-project.org/R-packages/trunk/cluster/) 
- [ggplot2](https://ggplot2.tidyverse.org/) 
- [scale](https://scales.r-lib.org/)  
