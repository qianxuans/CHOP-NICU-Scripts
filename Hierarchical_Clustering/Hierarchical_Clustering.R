#!/usr/bin/env Rscript
# Function to check and install required packages
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) {
    cat("Installing missing packages:", paste(new_packages, collapse=", "), "\n")
    install.packages(new_packages, repos="https://cran.rstudio.com/")
  }
}
# List of required packages
required_packages <- c("TreeTools",
                       "ape",
                       "dplyr",
                       "cluster")

# Check and install required packages
install_if_missing(required_packages)
# Load packages
library(TreeTools)
library(ape)
library(dplyr)
library(cluster)

# Access command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Initialize variables to store flag values
tree_path <- NULL
matrix_path <- NULL

# Loop through arguments and assign values to the appropriate variables
for (i in seq_along(args)) {
  if (args[i] == "--tree") {
    tree_path <- as.character(args[i + 1])
  } else if (args[i] == "--matrix") {
    matrix_path <- as.character(args[i + 1])
  }
}

# Check if both arguments are provided
if (is.null(tree_path) || is.null(matrix_path)) {
  warning("Missing required arguments. Both --tree and --matrix must be provided.")
  quit(status = 1)
}

#load the tree
everything_tree <- read.tree(tree_path)
#mid-point root of the tree
everything_tree <- phytools::midpoint_root(everything_tree)
#pre-order the tree
everything_tree <- Preorder(everything_tree)
#Matrix of pairwise distances from tree
distance_matrix <-  cophenetic.phylo(everything_tree)
#the various methods will become options in the future
hclust_result <- hclust(as.dist(distance_matrix),
                        method = "complete")
#hierarchical clustering
#We test the range from 2 to 100 clustering patterns 
hclust_result_silhouette_scores <- numeric()

for(i in 2:(length(everything_tree$tip.label)-1)){
  clusters <- cutree(hclust_result, k = i)
  hclust_result_silhouette_scores <- c(hclust_result_silhouette_scores,
                                       mean(silhouette(clusters,distance_matrix)[, "sil_width"]))
  remove(clusters)
}
remove(i)

final_group <- cutree(hclust_result,
                      k = (which.max(hclust_result_silhouette_scores) + 1))

hc_groups_df <- data.frame(genome = names(final_group),
                           group = final_group)

hc_groups_stats_df <- as.data.frame(table(hc_groups_df$group))

colnames(hc_groups_stats_df) <- c("group",
                                  "genome")
#For the genome(s) that is(are) not in any group divided using hierarchical clustering
#They are definitely singletons
hc_singetons <- hc_groups_df$genome[hc_groups_df$group %in% hc_groups_stats_df$group[hc_groups_stats_df$genome == 1]]

#Look at each group to check the minimal pair-wise SNP distance
#If in the clade the minimal SNP distance is already more than 500(can be changed),
#The isolates are all singletons
#For those groups, although no less than 2 isolates in the clade, 
#There is no need to build the sub-tree for them.
#Those clades will be labeled "1" in the over_snp_limit column
hc_groups_stats_df <- hc_groups_stats_df %>% mutate(over_snp_limit = NA)
#Read the SNP matrix of NICU genomes comparing to NICU genomes
study_snp_df <- readRDS(matrix_path)

for(i in seq_len(nrow(hc_groups_stats_df))){
  group = hc_groups_stats_df$group[i]
  group_isolates <- hc_groups_df$genome[hc_groups_df$group == group]
  #Study-genome SNP distance matrix within this group
  if(hc_groups_stats_df$genome[i] > 1){
    group_study_snp_df <- study_snp_df %>% filter(subject %in% group_isolates &
                                                    query %in% group_isolates)
    if(unique(min(group_study_snp_df$snp)) >= 500){
      hc_groups_stats_df$over_snp_limit[i] <- 1
    }else{
      hc_groups_stats_df$over_snp_limit[i] <- 0
    }
    remove(group_study_snp_df)
  }else{
    hc_groups_stats_df$over_snp_limit[i] <- 0
  }
  remove(group)
  remove(group_isolates)
}
remove(i)
#For the isolates in the groups with over_snp_limit = 1
#They are considered singletons at this step 
hc_singetons <- c(hc_singetons,
                  hc_groups_df$genome[hc_groups_df$group %in% hc_groups_stats_df$group[hc_groups_stats_df$over_snp_limit == 1]])
#clean
remove(final_group)
remove(hclust_result_silhouette_scores)
remove(distance_matrix)
remove(hclust_result)
remove(everything_tree)
#export output
write.csv(hc_groups_df,
          "all_hierarchical_clustering_groups.csv",
          quote = FALSE,
          row.names = FALSE)

selected_hc_groups_df <- hc_groups_df %>% 
  filter(group %in% hc_groups_stats_df$group[hc_groups_stats_df$over_snp_limit == 0 &
                                               hc_groups_stats_df$genome >= 2])
  
write.csv(selected_hc_groups_df,
          "selected_hierarchical_clustering_groups.csv",
          quote = FALSE,
          row.names = FALSE)

singletons_hc_groups_df <- hc_groups_df %>% 
  filter(group %in% hc_groups_stats_df$group[hc_groups_stats_df$over_snp_limit == 1 |
                                               hc_groups_stats_df$genome < 2])
write.csv(singletons_hc_groups_df,
          "singletons_hierarchical_clustering_groups.csv",
          quote = FALSE,
          row.names = FALSE)