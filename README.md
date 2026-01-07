# CHOP NICU Staphylococcus aureus Transmission Analysis Repository
## Genome assembly and quality control
Raw sequencing reads were subjected to quality control using Sunbeam (v4.3.7), which performed adapter trimming and host genome decontamination. Genome assembly was carried out using the following tools:

- SPAdes (v3.15.5) for isolate genomes.
- Anvi’o (v8) pipeline with MEGAHIT (v1.2.9) for 102 metagenomic samples.

Assembly quality was assessed using CheckM, applying the following criteria:
- ≥95% CheckM completeness
- ≤5% CheckM contamination
Lineage classification as "Staphylococcus (UID301)"

To evaluate species-level contamination, Mash was used. Assemblies were filtered to retain only those between 2.55 Mb and 3.15 Mb in size. Of the 1,670 assembled genomes, 1,446 passed the quality control criteria and were included in downstream analyses.

## Determination of clonality 
### Clonality Determination
Clonality was determined using a multi-step approach:
- Hierarchical Clustering: Genomes were stratified into distinct groups.
- Single-Linkage Clustering (SLC): Iterative SLC identified closely related genomes across SNP thresholds (minimum pairwise distance in the group to 500 SNPs).
- Phylogenetic Correction: Strain compositions were validated and corrected using reference-based maximum likelihood phylogenies. The compositions were expanded to include the smallest monophyletic clade with robust bootstrap support (≥70%).
- SNP Threshold Finilization: Final thresholds were determined where cluster composition and number plateaued.
  
This method ensured robust identification and validation of clonal relationships across the dataset.


### Multi-SNP-Threshold Plot
These plots visualize the following metrics at each tested SNP threshold:

- Number of clones and singletons.
- Number of discrepant genomes observed when SNP-only strain composition was mapped to the phylogeny.
- Bootstrap support for the strain compositions.
  
This approach ensures a comprehensive view of clonality across varying thresholds.

## Cluster analysis and visualizations
### Swimmer plot
This plot provides a temporal analysis of transmission clusters based on patient location and treatment team assignments:
- Patient Timelines: NICU section assignments are represented as colored rectangles. Admission dates are marked by triangles, and discharge dates by squares.
- Treatment Team Assignments: Shown as colored rectangles over time, corresponding to each patient’s timeline.

In both plots, sampling events are represented by dots:
Red dots: Invasive isolates.
Blue dots: Colonizing isolates.

These visualizations help illustrate the temporal and spatial dynamics within transmission clusters.

**Note: The input data is not provided as it contains sensitive patient information that cannot be shared to ensure privacy and confidentiality.**

### Cluster Floorplan
This script generates a NICU floorplan visualization(1 static and 2 animated floorplans for each cluster) to analyze patient movements in transmission dynamics. The plot highlights how proximity in time and space drives transmission events within the NICU. This visualization provides insights into spatial dynamics and helps pinpoint areas requiring targeted interventions.

**Note: The input data is not provided as it contains sensitive patient information that cannot be shared to ensure privacy and confidentiality.**

### Cluster plot
This script generates visualizations to summarize transmission clusters and their characteristics:

#### Cluster Representation: Each cluster is displayed as a box, colored based on invasive status:
-  Red: Invasive clusters.
-  Blue: Colonizing clusters.

#### Annotations:
##### Methicillin resistance status is indicated by letters: R (MRSA) or S (MSSA).
##### Patient status is represented by dots:
- Blue: Colonization only.
- Red: Infection only.
- Yellow: Both colonization and infection.
- Individual patient IDs are shown as numbers.
- Environmental isolates are marked by grey diamonds.
#### Additional Details:
Two bars below each cluster box indicate:
- The specific NICU section where the cluster was detected.
- The assigned treatment team at the time of detection.


### Spatial Distance Rescaled Phylogenetic Analysis
#### Analysis and Visualization of Spatial Distance Rescaled Phylogeny
The script analyzes clusters containing at least four genomes with maximum likelihood phylogenies and collection site coordinates from floor plans. The script calculates pairwise spatial distances between all genomes using the Pythagorean theorem based on collection site coordinates, generating a Naïve Spatial Distance Matrix measured in meters. A Spatial Distance Corrected Phylogeny is constructed using the nnls.tree() function from the phangorn R package. This function preserves the original ML phylogeny topology while rescaling branch lengths with spatial distances through ordinary least squares (OLS) regression to minimize residual sum of squares.The Phylogenetic Corrected Spatial Distance Matrix is computed using the cophenetic.phylo() function from the ape package, providing phylogenetically-informed spatial distances between samples. To track the spread distance accumulated over days post initial detection, samples collected each day are compared to all previously collected samples using phylogenetic corrected spatial distances. Single linkage clustering determines the minimum distance between each new isolate and all previously collected samples, representing the daily distance contribution.

**Note: The input data is not provided as it contains sensitive patient information that cannot be shared to ensure privacy and confidentiality.**

#### Visualization of Spatiotemporal Cluster Spread
The script generates dual-panel visualizations showing how above clusters spread over time. The upper panel displays a line plot with daily new distances (red) and accumulated distances (green) plotted against days post initial detection. The lower panel shows a floor plan with collection locations marked by colored dots, where colors represent the number of days after initial detection. The color bar positioning corresponds to the x-axis timeframe in the upper panel, providing a direct visual connection between temporal spread and spatial distribution.

**Note: The input data is not provided as it contains sensitive patient information that cannot be shared to ensure privacy and confidentiality.**

### Transmission Cluster Source Analysis
#### Molecular Clock Analysis
##### Beast2 Tree Visualization
The script generates visualizations of BEAST2 phylogenetic trees in NEXUS format, incorporating cluster metadata. Phylogenetic trees are annotated with blue bars representing the 95% highest posterior density (HPD) interval for the time to most recent common ancestor (tMRCA).

#### Transmission Tree Analysis
This script implements TransPhylo (https://pubmed.ncbi.nlm.nih.gov/28100788/) to infer transmission trees and estimate 95% confidence intervals for transmission index case dates using timed phylogenies generated using BEAST2. 

###  Transmission Cluster Source Plot
This script generates comparative visualizations plotting the time to most recent common ancestor (tMRCA) and transmission tree index cases against index patient admission dates. The analysis estimates whether transmission clusters originated from pre-colonized patients upon admission. 

**Note: The input data is not provided as it contains sensitive patient information that cannot be shared to ensure privacy and confidentiality.**


## Hidden Markov Model Analysis and Visualization
This script applies a Hidden Markov Model (HMM) to estimate the probability of failing to detect cluster presence during
the surveillance period (October 2021 – June 2024). Persistence periods were stratified into monthly time units. Two analyses were
performed: (a) probability of missing positive samples at monthly time unit within the surveillance (i.e., false negatives during
observed intervals), and (b) probability of failing to detect cluster presence outside the defined persistence period (i.e., undetected
persistence beyond first and last positive samples). Both analyses compared colonizing versus invasive clusters. Model convergence
was verified using multiple different starting matrices, all yielding identical estimates. Confidence intervals were generated via
bootstrapping (B = 100).

