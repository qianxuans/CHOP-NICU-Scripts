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

### Molecular_Clock_Analysis
#### Beast2 Tree Visualization
The script generates visualizations of BEAST2 phylogenetic trees in NEXUS format, incorporating cluster metadata. Phylogenetic trees are annotated with blue bars representing the 95% highest posterior density (HPD) interval for the time to most recent common ancestor (tMRCA).
#### Origin Analysis Plot
The script creates a comparative visualization that plots cluster tMRCA against index patient admission dates. This analysis revealed the temporal dynamics of cluster emergence, specifically demonstrating that 50% of identified clusters originated from index patient admission periods.
**Note: The input data is not provided as it contains sensitive patient information that cannot be shared to ensure privacy and confidentiality.**